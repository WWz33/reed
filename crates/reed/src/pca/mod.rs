//! `reed pca`: principal components of the sample relationship matrix
//! (doc_v16 §P2).
//!
//! The data flow mirrors GCTA/plink2: standardize the genotypes per SNP
//! (binom2 by default), build the N x N relationship matrix, and extract its
//! leading eigenpairs. Missing genotypes are mean-imputed or handled
//! pairwise-complete; both are plink2 conventions and both are checked
//! against plink2 in the e2e tests.
//!
//! Output files follow plink so existing scripts can consume them:
//! `PREFIX.eigenvec` (`#FID IID PC1..PCk`), `PREFIX.eigenval`, `PREFIX.pve`.

use std::io::Write as _;
use std::path::Path;

use reed_engine::budget::ThreadBudget;

use crate::genotype::{collect_via_engine, CollectSink, GenotypeBlock, PcaJob};
use reed_reduce::pca::{pca_dense_pairwise_with, pca_dense_with, Stand};
use reed_reduce::{lanczos_top_k, Mat, StreamConfig};



/// Reads the input into a samples x sites dosage matrix (doc_v17).
///
/// Shared by `reed pca` and `reed distance`: one parser, one set of site
/// gates (biallelic, diploid GT, at least one alt allele).
///
/// # Errors
///
/// Engine and I/O failures.
pub fn read_genotype_matrix(
    mut opened: crate::OpenedInput,
    keep: Option<&[bool]>,
    budget: &ThreadBudget,
) -> Result<(Mat, Vec<String>), crate::Error> {
    let sample_names = opened.facts().sample_names.clone();
    if sample_names.is_empty() {
        return Err(crate::Error::usage(
            "the header declares no samples; nothing to compute",
        ));
    }
    let n_samples = sample_names.len();
    let keep: Vec<bool> = match keep {
        Some(mask) => {
            if mask.len() != n_samples {
                return Err(crate::Error::usage(format!(
                    "sample keep mask has {} entries but the header declares {n_samples}",
                    mask.len()
                )));
            }
            mask.to_vec()
        }
        None => vec![true; n_samples],
    };
    if !keep.iter().any(|&k| k) {
        return Err(crate::Error::usage("the sample subset selects no samples"));
    }
    let job = PcaJob { n_samples, keep: keep.clone(), allow_multiallelic: false };
    let abort = reed_engine::abort::abort_for(budget);
    let sink = CollectSink { acc: GenotypeBlock::default() };
    let block = collect_via_engine(&mut opened, &job, sink, budget, &abort)?;
    if block.sites.is_empty() {
        return Err(crate::Error::usage(
            "no biallelic sites with a diploid GT were found",
        ));
    }
    let kept_indices: Vec<usize> = (0..n_samples).filter(|&i| keep[i]).collect();
    let n_kept = kept_indices.len();
    let n_sites = block.sites.len();
    let mut x = Mat::zeros(n_kept, n_sites);
    for s in 0..n_sites {
        for col in 0..n_kept {
            x.set(col, s, block.dosages[s * n_kept + col]);
        }
    }
    let ids = kept_indices
        .iter()
        .map(|&i| sample_names[i].clone())
        .collect();
    Ok((x, ids))
}

/// Request for [`run_pca`].
#[derive(Clone, Debug)]
pub struct PcaRequest {
    /// Number of principal components.
    pub pcs: usize,
    /// `binom2` (default) or `binom`.
    pub stand: String,
    /// `pairwise` (default) or `meanimpute`.
    pub missing: String,
    /// `auto`, `dense` or `stream`.
    pub solver: String,
    /// Output prefix; writes `.eigenvec`, `.eigenval`, `.pve`.
    pub prefix: String,
    /// Optional sample keep mask (header order); `None` keeps every sample.
    pub keep_samples: Option<Vec<bool>>,
    /// Worker threads.
    pub threads: Option<usize>,
    /// Memory budget in MiB.
    pub memory: Option<u64>,
}

impl Default for PcaRequest {
    fn default() -> Self {
        Self {
            pcs: 10,
            stand: "binom2".into(),
            missing: "pairwise".into(),
            solver: "auto".into(),
            prefix: "out".into(),
            keep_samples: None,
            threads: None,
            memory: None,
        }
    }
}

/// Report of one PCA run.
#[derive(Clone, Debug)]
pub struct PcaReport {
    /// Samples used.
    pub n_samples: usize,
    /// Sites used.
    pub n_sites: usize,
    /// PCs written.
    pub n_pcs: usize,
    /// Eigenvalues (descending).
    pub eigenvalues: Vec<f64>,
    /// Solver actually used (`dense` or `stream`).
    pub solver: String,
    /// Wall time.
    pub elapsed: std::time::Duration,
}

/// Runs PCA over `path` (or stdin when `None`).
///
/// # Errors
///
/// Usage errors for a bad flag combination, and I/O or engine failures.
pub fn run_pca(
    path: Option<&std::path::Path>,
    request: &PcaRequest,
    budget: &ThreadBudget,
) -> Result<PcaReport, crate::Error> {
    let opened = crate::open_regions(path, &[])?;
    run_pca_opened(opened, request, budget)
}

/// Runs PCA over an already-opened input (the CLI opens once so stdin streams).
///
/// # Errors
///
/// As [`run_pca`].
pub fn run_pca_opened(
    mut opened: crate::OpenedInput,
    request: &PcaRequest,
    budget: &ThreadBudget,
) -> Result<PcaReport, crate::Error> {
    if request.pcs == 0 {
        return Err(crate::Error::usage("--pcs must be at least 1"));
    }
    let stand = Stand::parse(&request.stand).map_err(crate::Error::usage)?;
    let missing = request.missing.as_str();
    if !matches!(missing, "pairwise" | "meanimpute") {
        return Err(crate::Error::usage(format!(
            "--missing {missing:?} is not supported (pairwise|meanimpute)"
        )));
    }
    let solver = request.solver.as_str();
    if !matches!(solver, "auto" | "dense" | "stream") {
        return Err(crate::Error::usage(format!(
            "--solver {solver:?} is not supported (auto|dense|stream)"
        )));
    }

    let sample_names = opened.facts().sample_names.clone();
    if sample_names.is_empty() {
        return Err(crate::Error::usage(
            "the header declares no samples; nothing to compute",
        ));
    }
    let n_samples = sample_names.len();
    let keep: Vec<bool> = match &request.keep_samples {
        Some(mask) => {
            if mask.len() != n_samples {
                return Err(crate::Error::usage(format!(
                    "sample keep mask has {} entries but the header declares {n_samples}",
                    mask.len()
                )));
            }
            mask.clone()
        }
        None => vec![true; n_samples],
    };
    if !keep.iter().any(|&k| k) {
        return Err(crate::Error::usage("the sample subset selects no samples"));
    }

    let job = PcaJob {
        n_samples,
        keep: keep.clone(),
        allow_multiallelic: false,
    };
    let abort = reed_engine::abort::abort_for(budget);
    let sink = CollectSink { acc: GenotypeBlock::default() };
    let block = collect_via_engine(&mut opened, &job, sink, budget, &abort)?;

    let n_sites = block.sites.len();
    if n_sites == 0 {
        return Err(crate::Error::usage(
            "no biallelic sites with a diploid GT were found",
        ));
    }
    let kept_indices: Vec<usize> = (0..n_samples).filter(|&i| keep[i]).collect();
    let n_kept = kept_indices.len();

    // Transpose site-major -> sample-major, which is what reed_reduce wants.
    let mut x = Mat::zeros(n_kept, n_sites);
    for s in 0..n_sites {
        for col in 0..n_kept {
            x.set(col, s, block.dosages[s * n_kept + col]);
        }
    }
    let _ = kept_indices;

    let n_bytes = (n_kept as u64)
        .saturating_mul(n_kept as u64)
        .saturating_mul(8);
    let dense_fits = n_bytes <= budget.mem.bytes;
    let chosen = match solver {
        "dense" => "dense",
        "stream" => "stream",
        _ => {
            if dense_fits {
                "dense"
            } else {
                "stream"
            }
        }
    };

    let started = std::time::Instant::now();
    let (values, vectors, pve) = if chosen == "stream" {
        let st = reed_reduce::pca::standardize(&x, stand);
        let cfg = StreamConfig {
            k: request.pcs,
            tol: 1e-10,
            max_cycles: 200,
            ..Default::default()
        };
        let r = lanczos_top_k(&st.w, &cfg);
        let mut tr = 0.0;
        for i in 0..st.w.rows {
            for v in st.w.row(i) {
                tr += v * v;
            }
        }
        let total = tr / st.w.cols as f64;
        let pve = r
            .values
            .iter()
            .map(|v| if total == 0.0 { 0.0 } else { v / total })
            .collect();
        (r.values, r.vectors, pve)
    } else if missing == "pairwise" {
        let p = pca_dense_pairwise_with(&x, request.pcs, stand);
        (p.values, p.vectors, p.pve)
    } else {
        let p = pca_dense_with(&x, request.pcs, stand);
        (p.values, p.vectors, p.pve)
    };
    let elapsed = started.elapsed();

    write_outputs(
        &request.prefix,
        &sample_names,
        &kept_indices,
        &values,
        &vectors,
        &pve,
    )?;

    Ok(PcaReport {
        n_samples: n_kept,
        n_sites,
        n_pcs: values.len(),
        eigenvalues: values,
        solver: chosen.to_string(),
        elapsed,
    })
}

/// Writes `.eigenvec`, `.eigenval` and `.pve` next to `prefix`.
fn write_outputs(
    prefix: &str,
    sample_names: &[String],
    kept: &[usize],
    values: &[f64],
    vectors: &Mat,
    pve: &[f64],
) -> Result<(), crate::Error> {
    let vec_path = format!("{prefix}.eigenvec");
    let val_path = format!("{prefix}.eigenval");
    let pve_path = format!("{prefix}.pve");

    let mut f = std::fs::File::create(&vec_path)
        ?;
    // plink's header: `#FID IID PC1..PCk`; reed uses the sample name for both
    // FID and IID so single-name datasets round-trip.
    write!(f, "#FID\tIID").map_err(io_err)?;
    for c in 0..values.len() {
        write!(f, "\tPC{}", c + 1).map_err(io_err)?;
    }
    writeln!(f).map_err(io_err)?;
    for (row, &i) in kept.iter().enumerate() {
        let name = &sample_names[i];
        write!(f, "{name}\t{name}").map_err(io_err)?;
        for c in 0..values.len() {
            write!(f, "\t{:.9}", vectors.get(row, c)).map_err(io_err)?;
        }
        writeln!(f).map_err(io_err)?;
    }

    let mut f = std::fs::File::create(&val_path)
        ?;
    for v in values {
        writeln!(f, "{v:.9}").map_err(io_err)?;
    }

    let mut f = std::fs::File::create(&pve_path)
        ?;
    for v in pve {
        writeln!(f, "{v:.9}").map_err(io_err)?;
    }
    let _ = Path::new(prefix);
    Ok(())
}

fn io_err(e: std::io::Error) -> crate::Error {
    crate::Error::from(e)
}
