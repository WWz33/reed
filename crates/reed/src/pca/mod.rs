//! `reed pca`: principal components of the sample relationship matrix
//! (doc_v16 §P2).
//!
//! The data flow mirrors GCTA/plink2: standardize the genotypes per SNP
//! (binom2 by default), build the N x N relationship matrix, and extract its
//! leading eigenpairs. Missing genotypes are mean-imputed (both solvers) or
//! handled pairwise-complete (dense solver only — a plink2 convention);
//! because the stream solver has no pairwise-complete estimator, it refuses
//! `--missing pairwise` instead of silently mean-imputing (doc_v20 R10).
//! The e2e tests check the two policies and the two solvers against each
//! other on fixtures; there is no plink2 comparison.
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
    // doc_v20 review of S76: cap the collection itself, not just the final
    // matrix, so a 120 GB run refuses while scanning instead of after.
    let n_kept_hint = keep.iter().filter(|&&k| k).count();
    let exceeded = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
    let sink = CollectSink {
        acc: GenotypeBlock::default(),
        max_bytes: Some(budget.mem.bytes),
        n_samples: n_kept_hint,
        exceeded: std::sync::Arc::clone(&exceeded),
    };
    let collected = collect_via_engine(&mut opened, &job, sink, budget, &abort);
    if exceeded.load(std::sync::atomic::Ordering::Acquire) {
        // doc_v20 review of S76: the cap fires during collection, so this
        // is the actionable refusal (usage class), not the engine's abort.
        return Err(crate::Error::usage(format!(
            "the dosage matrix exceeds the memory budget of {} bytes while \
             collecting; LD-prune the input first \
             (`reed filter --indep-pairwise`) or raise --memory",
            budget.mem.bytes
        )));
    }
    let block = collected?;
    if block.n_sites == 0 {
        return Err(crate::Error::usage(
            "no biallelic sites with a diploid GT were found",
        ));
    }
    let kept_indices: Vec<usize> = (0..n_samples).filter(|&i| keep[i]).collect();
    let n_kept = kept_indices.len();
    let n_sites = block.n_sites;
    // doc_v20 S76: the dense samples x sites matrix is the command's
    // dominant allocation and must answer to `--memory` — 5,000 samples x
    // 3M unpruned WGS sites is ~120 GB, which used to be attempted anyway
    // (an OOM kill, not a usage error). LD pruning is an upstream step per
    // the README, so a refusal points there.
    let matrix_bytes = (n_kept as u64)
        .saturating_mul(n_sites as u64)
        .saturating_mul(8);
    // Peak residency is the collected block plus the matrix (and the dense
    // path adds the standardized copy), so the gate is on the peak, not on
    // the matrix alone (doc_v20 review of S76).
    let peak_bytes = matrix_bytes.saturating_mul(2);
    if peak_bytes > budget.mem.bytes {
        const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
        return Err(crate::Error::usage(format!(
            "the dosage matrix needs {} bytes ({:.1} GiB) for {n_kept} samples \
             x {n_sites} sites, over the memory budget of {} bytes \
             ({:.1} GiB); LD-prune the input first \
             (`reed filter --indep-pairwise`) or raise --memory",
            matrix_bytes,
            matrix_bytes as f64 / GIB,
            budget.mem.bytes,
            budget.mem.bytes as f64 / GIB,
        )));
    }
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
    opened: crate::OpenedInput,
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

    // doc_v20 review of R10: decide the solver refusal before reading the
    // input. `n_kept` is known from the header and the keep mask, so an
    // explicit `--solver stream --missing pairwise` (or an `auto` that will
    // pick stream because N²·8 exceeds the budget) fails in seconds instead
    // of after scanning a whole genome and building the matrix.
    {
        let n_samples = opened.facts().sample_names.len();
        let n_kept = match &request.keep_samples {
            Some(mask) => mask.iter().filter(|&&k| k).count(),
            None => n_samples,
        };
        let n_bytes = (n_kept as u64)
            .saturating_mul(n_kept as u64)
            .saturating_mul(8);
        let will_stream = match solver {
            "stream" => true,
            "dense" => false,
            _ => n_bytes > budget.mem.bytes,
        };
        if will_stream && missing == "pairwise" {
            return Err(crate::Error::usage(
                "--missing pairwise needs the dense solver; \
                 pick --solver dense, --missing meanimpute, or raise --memory",
            ));
        }
    }

    // doc_v20 S76: the collection, the sample-subset validation and the
    // dosage-matrix memory check live in [`read_genotype_matrix`] — the
    // v20 review found this copy and that reader near-duplicates, and the
    // second `Mat::zeros` is exactly where the budget check would be
    // forgotten again.
    let sample_names = opened.facts().sample_names.clone();
    let (x, _kept_names) = read_genotype_matrix(opened, request.keep_samples.as_deref(), budget)?;
    // The reader above validated the mask's length against the header, so
    // the per-sample indexing here is in-bounds by construction.
    let n_samples = sample_names.len();
    let kept_indices: Vec<usize> = match &request.keep_samples {
        Some(mask) => (0..n_samples).filter(|&i| mask[i]).collect(),
        None => (0..n_samples).collect(),
    };
    let n_kept = kept_indices.len();
    let n_sites = x.cols;

    let started = std::time::Instant::now();
    let c = solve(&x, request.pcs, stand, missing, solver, budget)?;
    let (values, vectors, pve) = (c.values, c.vectors, c.pve);
    let chosen = c.solver;
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

/// One PCA solve: the leading eigenpairs and their share of the spectrum.
#[derive(Clone, Debug)]
pub(crate) struct Components {
    /// Eigenvalues, descending (`lambda / m`).
    pub values: Vec<f64>,
    /// Eigenvectors as columns (`n_samples` x `k`), unit norm.
    pub vectors: Mat,
    /// Per-PC share of the **whole** spectrum. The entries sum to 1 when `k`
    /// covers the full rank, and to less when it does not — which is what
    /// makes `sum(pve[..n])` a meaningful "variance explained".
    pub pve: Vec<f64>,
    /// Solver actually used (`dense` or `stream`).
    pub solver: &'static str,
}

/// Runs the PCA step over a samples x sites dosage matrix.
///
/// Both `reed pca` and `reed pca --dapc` call this: the solver choice, the
/// three estimator branches and the variance normalisation must stay one
/// implementation, or the two commands drift apart.
///
/// # Errors
///
/// A usage error when the effective solver is stream and the missing policy
/// is pairwise — including the `auto` switch to stream on a tight memory
/// budget (doc_v20 R10: refuse, never silently mean-impute).
pub(crate) fn solve(
    x: &Mat,
    k: usize,
    stand: Stand,
    missing: &str,
    solver: &str,
    budget: &ThreadBudget,
) -> Result<Components, crate::Error> {
    let n = x.rows;
    let n_bytes = (n as u64).saturating_mul(n as u64).saturating_mul(8);
    let chosen = match solver {
        "dense" => "dense",
        "stream" => "stream",
        _ => {
            if n_bytes <= budget.mem.bytes {
                "dense"
            } else {
                "stream"
            }
        }
    };
    // The stream solver has no pairwise-complete estimator: it mean-imputes
    // through `standardize`. Refuse the combination instead of downgrading
    // the policy silently (doc_v20 R10).
    if chosen == "stream" && missing == "pairwise" {
        return Err(crate::Error::usage(
            "--missing pairwise needs the dense solver; \
             pick --solver dense, --missing meanimpute, or raise --memory",
        ));
    }

    if chosen == "stream" {
        let st = reed_reduce::pca::standardize(x, stand);
        let cfg = StreamConfig {
            k,
            tol: 1e-10,
            max_cycles: 200,
            ..Default::default()
        };
        let r = lanczos_top_k(&st.w, &cfg);
        // The stream solver has no full spectrum, but the trace of `W W' / m`
        // is the eigenvalue sum, so the normalisation is still the whole
        // spectrum rather than the retained part.
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
        Ok(Components {
            values: r.values,
            vectors: r.vectors,
            pve,
            solver: "stream",
        })
    } else {
        let p = if missing == "pairwise" {
            pca_dense_pairwise_with(x, k, stand)
        } else {
            pca_dense_with(x, k, stand)
        };
        Ok(Components {
            values: p.values,
            vectors: p.vectors,
            pve: p.pve,
            solver: "dense",
        })
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use reed_engine::budget::MemBudget;

    /// A plain-VCF fixture: `n_samples` diploid samples over `n_sites`
    /// biallelic sites, dosages cycling 0/1/2 so every site varies.
    fn write_fixture(dir: &std::path::Path, n_samples: usize, n_sites: usize) -> std::path::PathBuf
    {
        let path = dir.join("in.vcf");
        std::fs::write(&path, fixture_text(n_samples, n_sites)).expect("invariant: fixture");
        path
    }

    /// The fixture's text ([`write_fixture`]'s bytes).
    fn fixture_text(n_samples: usize, n_sites: usize) -> String {
        let mut text = String::from(
            "##fileformat=VCFv4.3\n##contig=<ID=1,length=1000000>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        );
        for i in 0..n_samples {
            text.push_str("\ts");
            text.push_str(&i.to_string());
        }
        text.push('\n');
        for s in 0..n_sites {
            text.push_str("1\t");
            text.push_str(&(100 + s).to_string());
            text.push_str("\t.\tA\tG\t50\tPASS\t.\tGT");
            for i in 0..n_samples {
                let gt = match (i + s) % 3 {
                    0 => "0/0",
                    1 => "0/1",
                    _ => "1/1",
                };
                text.push('\t');
                text.push_str(gt);
            }
            text.push('\n');
        }
        text
    }

    /// The same fixture text, bgzf-compressed in 4 KiB blocks: the
    /// engine's per-chunk footprint (compressed bytes x 8) stays far
    /// below the dosage matrix bytes, so a budget between the two
    /// isolates the S76 matrix check from the engine's own guard.
    fn write_fixture_bgzf(
        dir: &std::path::Path,
        n_samples: usize,
        n_sites: usize,
    ) -> std::path::PathBuf {
        let text = fixture_text(n_samples, n_sites);
        let mut deflater = reed_io::bgzf::deflate::Deflater::new(6)
            .expect("invariant: level 6 is valid");
        let mut out = Vec::new();
        for block in text.as_bytes().chunks(4_096) {
            deflater
                .deflate_block(block, &mut out)
                .expect("invariant: block fits");
        }
        out.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
        let path = dir.join("in.vcf.gz");
        std::fs::write(&path, out).expect("invariant: fixture");
        path
    }

    fn temp_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("reed-pca-r10-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        dir
    }

    /// doc_v20 R10: the stream solver mean-imputes; `--missing pairwise`
    /// asked for with it explicitly must be a usage error, not a silent
    /// downgrade to mean imputation.
    #[test]
    fn r10_stream_solver_refuses_missing_pairwise() {
        let dir = temp_dir("stream");
        let input = write_fixture(&dir, 6, 12);
        let request = PcaRequest {
            pcs: 2,
            missing: "pairwise".into(),
            solver: "stream".into(),
            prefix: dir.join("out").to_string_lossy().into_owned(),
            ..PcaRequest::default()
        };
        let budget = ThreadBudget {
            threads: 2,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        match run_pca(Some(&input), &request, &budget) {
            Ok(report) => panic!(
                "stream + pairwise must be refused, not silently mean-imputed \
                 (doc_v20 R10); got Ok with solver {}",
                report.solver
            ),
            Err(e) => {
                assert_eq!(e.kind(), crate::ErrorKind::Usage, "{e}");
                let msg = e.to_string();
                assert!(
                    msg.contains("--solver dense"),
                    "advice must offer the dense solver: {msg}"
                );
            }
        }
    }

    /// doc_v20 R10: the same refusal when `--solver auto` switches to the
    /// stream solver because the memory budget cannot hold the dense N x N
    /// relationship matrix — the large-N path that used to mean-impute
    /// silently.
    #[test]
    fn r10_auto_solver_refuses_missing_pairwise_when_memory_forces_stream() {
        let dir = temp_dir("auto");
        // 100 samples: the dense matrix needs 8 * 100^2 = 80_000 bytes,
        // above the 60_000-byte budget, while the input chunks stay well
        // under it (the engine's per-chunk footprint check must pass).
        let input = write_fixture(&dir, 100, 12);
        let request = PcaRequest {
            pcs: 2,
            missing: "pairwise".into(),
            solver: "auto".into(),
            prefix: dir.join("out").to_string_lossy().into_owned(),
            ..PcaRequest::default()
        };
        let budget = ThreadBudget {
            threads: 2,
            mem: MemBudget::from_bytes(60_000),
        };
        match run_pca(Some(&input), &request, &budget) {
            Ok(report) => panic!(
                "auto-forced stream + pairwise must be refused (doc_v20 R10); \
                 got Ok with solver {}",
                report.solver
            ),
            Err(e) => {
                assert_eq!(e.kind(), crate::ErrorKind::Usage, "{e}");
                assert!(
                    e.to_string().contains("--memory"),
                    "advice must offer raising the memory budget: {e}"
                );
            }
        }
    }

    /// doc_v20 S76: the dense samples x sites dosage matrix is the
    /// command's dominant allocation and used to be materialized with no
    /// check against `--memory` — 5,000 samples x 3M unpruned WGS sites
    /// is ~120 GB, dying by an OOM kill instead of a usage error. A
    /// budget the matrix cannot fit must be refused, naming the LD-prune
    /// escape hatch (LD pruning is an upstream step per the README).
    #[test]
    fn s76_dosage_matrix_beyond_the_budget_is_a_usage_error() {
        let dir = temp_dir("s76-refuse");
        // 8 kept samples x 1,500 sites = 96,000 bytes of matrix against a
        // 48,000-byte budget, while every bgzf chunk's engine footprint
        // (compressed x 8) stays far below 48,000.
        let input = write_fixture_bgzf(&dir, 8, 1_500);
        let request = PcaRequest {
            pcs: 2,
            prefix: dir.join("out").to_string_lossy().into_owned(),
            ..PcaRequest::default()
        };
        let budget = ThreadBudget {
            threads: 2,
            mem: MemBudget::from_bytes(48_000),
        };
        match run_pca(Some(&input), &request, &budget) {
            Ok(report) => panic!(
                "a dosage matrix over the memory budget must be a usage \
                 error, not a silent allocation (doc_v20 S76); got Ok with \
                 {} sites",
                report.n_sites
            ),
            Err(e) => {
                assert_eq!(e.kind(), crate::ErrorKind::Usage, "{e}");
                let msg = e.to_string();
                assert!(
                    msg.contains("--indep-pairwise"),
                    "advice must point at LD pruning: {msg}"
                );
                assert!(
                    msg.contains("--memory"),
                    "advice must offer raising the budget: {msg}"
                );
            }
        }
    }

    /// doc_v20 S76: the check must not fire inside the budget — the same
    /// entry point with a comfortable budget still produces a report.
    #[test]
    fn s76_comfortable_budget_still_succeeds() {
        let dir = temp_dir("s76-ok");
        let input = write_fixture_bgzf(&dir, 8, 60);
        let request = PcaRequest {
            pcs: 2,
            prefix: dir.join("out").to_string_lossy().into_owned(),
            ..PcaRequest::default()
        };
        let budget = ThreadBudget {
            threads: 2,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        let report = run_pca(Some(&input), &request, &budget)
            .expect("invariant: a comfortable budget succeeds (doc_v20 S76)");
        assert_eq!(report.n_samples, 8);
        assert_eq!(report.n_sites, 60);
        assert_eq!(report.eigenvalues.len(), 2);
    }
}
