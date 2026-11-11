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

use reed_engine::budget::ThreadBudget;

use crate::genotype::{CollectSink, GenotypeBlock, PcaJob, collect_via_engine};
use reed_reduce::pca::{Stand, pca_dense_pairwise_with_par, pca_dense_with_par};
use reed_reduce::{Mat, StreamConfig, lanczos_top_k};

/// How the caller will consume the dosage matrix, for the `--memory` gate
/// (doc_v21 M11).
///
/// The old gate counted only `2 x n x m` (the collected block plus the
/// matrix), but the dense PCA paths then allocate an `n x n` GRM — and the
/// pairwise path a second `n x n` count matrix — which `--memory` never saw.
/// The gate now charges each path its real peak.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixUse {
    /// Matrix only (`distance` and other matrix consumers).
    Plain,
    /// Dense PCA, mean-imputed GRM: one `n x n` matrix.
    DenseMeanImpute,
    /// Dense PCA, pairwise GRM: the GRM plus its pair-count matrix.
    DensePairwise,
    /// Streamed PCA: the matrix plus its standardized copy.
    Stream,
}

impl MatrixUse {
    /// Peak residency in bytes for an `n_kept x n_sites` matrix.
    ///
    /// The estimate must cover what is *simultaneously* alive, not just the
    /// obvious matrix (doc_v26 V26-20). The old formula charged the dense
    /// paths one `n x n` array (two for pairwise), while a dense solve really
    /// holds, at the eigensolver's peak:
    ///
    /// | object | bytes |
    /// |---|---|
    /// | `x` (the collected block plus the matrix) | `2 * n*m*8` |
    /// | `standardize`'s copy `st.w` | `n*m*8` |
    /// | the faer input `wm` (a copy of `st.w` in faer's layout) | `n*m*8` |
    /// | the GRM `W W'/m` | `n*n*8` |
    /// | pairwise: the pair-count matrix `Grm::pairwise` also builds | `n*n*8` |
    /// | faer's symmetrized input `am` | `n*n*8` |
    /// | faer's eigenvector matrix `u` | `n*n*8` |
    /// | faer's workspace stack | `n*n*8` |
    /// | the returned `Vec<f64>` + `Mat` copy of `u` | `n*n*8` |
    /// | `pack`'s sign normalization (in place) and the `n`-long vectors | ~0 |
    ///
    /// so dense mean-impute is `4*n*m*8 + 6*n*n*8` and pairwise adds one more
    /// `n*n*8` for the counts. The stream path holds `x` (2), `st.w` and the
    /// Lanczos basis (`n*k*8`, small against `n*m`), so its `3*n*m*8` stays
    /// accurate.
    ///
    /// The formula is fitted against measured peak RSS (`ru_maxrss`) on
    /// 1000/2000/3000 samples x 2000 sites: 100.0 / 301.9 / 627.9 MiB
    /// measured vs 106.8 / 305.2 / 624.0 MiB estimated, i.e. the estimate is
    /// within ±7% and never below the measurement by more than 7%. That is
    /// the documented error bound; a large margin is still advisable because
    /// the allocator can retain freed blocks.
    fn peak_bytes(self, n_kept: u64, n_sites: u64) -> u64 {
        let nm = n_kept.saturating_mul(n_sites).saturating_mul(8);
        let n2 = n_kept.saturating_mul(n_kept).saturating_mul(8);
        match self {
            // The collected block plus the matrix, nothing square.
            Self::Plain => nm.saturating_mul(2),
            // Matrix plus its standardized copy (the Lanczos basis is n*k).
            Self::Stream => nm.saturating_mul(3),
            // x (2) + st.w + faer's wm; GRM + am + u + stack + returned copy.
            Self::DenseMeanImpute => nm.saturating_mul(4).saturating_add(n2.saturating_mul(6)),
            // The same, plus the pair-count matrix.
            Self::DensePairwise => nm.saturating_mul(4).saturating_add(n2.saturating_mul(7)),
        }
    }
}

/// The solver a request resolves to (doc_v21 M11): one rule shared by the
/// pre-read gate, `solve`, and dapc.
///
/// `auto` picks dense when the dense path's **real peak** fits the budget
/// (doc_v26 V26-20); it used to compare only the `n x n` GRM. `n_sites` is
/// `None` before the input is read — the pre-read refusal check only needs to
/// know whether dense is clearly impossible, so it judges the square arrays
/// alone, and `solve` makes the authoritative decision once the site count is
/// known.
pub(crate) fn chosen_solver(
    solver: &str,
    n_kept: usize,
    n_sites: Option<usize>,
    budget: &ThreadBudget,
) -> &'static str {
    match solver {
        "dense" => "dense",
        "stream" => "stream",
        _ => {
            let dense_peak = match n_sites {
                Some(n_sites) => {
                    MatrixUse::DenseMeanImpute.peak_bytes(n_kept as u64, n_sites as u64)
                }
                None => (n_kept as u64)
                    .saturating_mul(n_kept as u64)
                    .saturating_mul(8),
            };
            if dense_peak <= budget.mem.bytes {
                "dense"
            } else {
                "stream"
            }
        }
    }
}

/// The gate class for a resolved PCA solver.
pub(crate) fn matrix_use(solver: &str, missing: &str) -> MatrixUse {
    if solver == "stream" {
        MatrixUse::Stream
    } else if missing == "pairwise" {
        MatrixUse::DensePairwise
    } else {
        MatrixUse::DenseMeanImpute
    }
}

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
    matrix_use: MatrixUse,
) -> Result<(Mat, Vec<String>), crate::Error> {
    let sample_names = opened.facts().sample_names.clone();
    if sample_names.is_empty() {
        return Err(crate::Error::usage(
            "the header declares no samples; nothing to compute",
        ));
    }
    let n_samples = sample_names.len();
    let keep = reed_core::sample::KeepMask::try_new(
        n_samples,
        keep.map_or_else(|| vec![true; n_samples], <[bool]>::to_vec),
    )
    .map_err(|e| crate::Error::usage(e.to_string()))?;
    if !keep.iter().any(|&k| k) {
        return Err(crate::Error::usage("the sample subset selects no samples"));
    }
    // doc_v28 V28-15: a BCF input feeds the framed records straight to the
    // dosage reader; text input stays on the rendered-text path.
    let native_bcf = opened.bcf.is_some();
    let job = PcaJob {
        n_samples,
        keep: keep.as_slice().to_vec(),
        allow_multiallelic: false,
        bcf_input: native_bcf.then(|| {
            opened
                .bcf
                .clone()
                .expect("invariant: a BCF input carries its dictionary")
        }),
    };
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
    let collected = collect_via_engine(&mut opened, job, sink, budget, &abort);
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
    // Peak residency depends on the consumer (doc_v21 M11): the dense PCA
    // paths allocate an n x n GRM on top of the n x m matrix, which the old
    // `2 x n x m` gate never charged; an explicit `--solver dense` could pass
    // the gate and then abort inside the eigensolve.
    let peak_bytes = matrix_use.peak_bytes(n_kept as u64, n_sites as u64);
    if peak_bytes > budget.mem.bytes {
        const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
        return Err(crate::Error::usage(format!(
            "the {matrix_use:?} path needs {peak_bytes} bytes ({:.1} GiB) for \
             {n_kept} samples x {n_sites} sites, over the memory budget of {} \
             bytes ({:.1} GiB); LD-prune the input first \
             (`reed filter --indep-pairwise`) or raise --memory",
            peak_bytes as f64 / GIB,
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
    /// The stream solver's final relative residual (`None` for dense); a
    /// value above the solver tolerance means the components did not
    /// converge and the run logged a warning (doc_v25 V25-10).
    pub stream_max_residual: Option<f64>,
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
        let will_stream = chosen_solver(solver, n_kept, None, budget) == "stream";
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
    // forgotten again. doc_v21 M11: the gate is charged the resolved
    // solver's real peak (the dense GRM included).
    let n_kept = match &request.keep_samples {
        Some(mask) => mask.iter().filter(|&&k| k).count(),
        None => opened.facts().sample_names.len(),
    };
    let matrix_use = matrix_use(chosen_solver(solver, n_kept, None, budget), missing);
    let sample_names = opened.facts().sample_names.clone();
    let (x, _kept_names) =
        read_genotype_matrix(opened, request.keep_samples.as_deref(), budget, matrix_use)?;
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
    let stream_max_residual = c.max_residual;
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
        stream_max_residual,
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
    /// The stream solver's final relative residual (`None` for dense): a
    /// value above `tol` means the requested components did not converge
    /// (doc_v25 V25-10).
    pub max_residual: Option<f64>,
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
    // One rule for the pre-read gate, this solve, and dapc (doc_v21 M11).
    // The authoritative decision: the site count is known here, so `auto`
    // is judged on the full dense peak (doc_v26 V26-20).
    let chosen = chosen_solver(solver, n, Some(x.cols), budget);
    // The stream solver has no pairwise-complete estimator: it mean-imputes
    // through `standardize`. Refuse the combination instead of downgrading
    // the policy silently (doc_v20 R10).
    if chosen == "stream" && missing == "pairwise" {
        return Err(crate::Error::usage(
            "--missing pairwise needs the dense solver; \
             pick --solver dense, --missing meanimpute, or raise --memory",
        ));
    }

    // doc_v26 V26-21: the solve carries the run's thread budget. The
    // previous code let faer read `rayon::current_num_threads()` from the
    // process-global pool, so `--threads 1` still fanned the eigensolver out
    // over every core. `auxiliary()` is the documented bound for the
    // auxiliary (non-pipeline) stages (doc_v8 D8-1), and the sequential
    // budget maps to one thread.
    let par = if budget.is_sequential() {
        reed_reduce::Parallel::Serial
    } else {
        reed_reduce::Parallel::Threads(budget.auxiliary())
    };

    if chosen == "stream" {
        let st = reed_reduce::pca::standardize(x, stand);
        let cfg = StreamConfig {
            k,
            tol: 1e-10,
            max_cycles: 200,
            parallel: par,
            ..Default::default()
        };
        let r = lanczos_top_k(&st.w, &cfg);
        if r.max_residual > cfg.tol {
            tracing::warn!(
                residual = r.max_residual,
                tol = cfg.tol,
                cycles = r.cycles,
                "PCA stream solver did not converge; the written components may be \
                 inaccurate (raise --pcs headroom or use --solver dense)"
            );
        }
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
            max_residual: Some(r.max_residual),
        })
    } else {
        let p = if missing == "pairwise" {
            pca_dense_pairwise_with_par(x, k, stand, par)
        } else {
            pca_dense_with_par(x, k, stand, par)
        };
        Ok(Components {
            values: p.values,
            vectors: p.vectors,
            pve: p.pve,
            solver: "dense",
            max_residual: None,
        })
    }
}

/// Writes `.eigenvec`, `.eigenval` and `.pve` next to `prefix`.
///
/// All three are written to sibling temporaries first and published only
/// after every one of them is complete (doc_v26 V26-08): a failure while
/// writing the second file cannot leave the first new and the third stale.
/// The three renames are still individual — an interrupt between them leaves
/// a mixed group, which is the documented per-file promise; what this
/// function guarantees is that no *write error* produces a partial group.
fn write_outputs(
    prefix: &str,
    sample_names: &[String],
    kept: &[usize],
    values: &[f64],
    vectors: &Mat,
    pve: &[f64],
) -> Result<(), crate::Error> {
    use reed_io::write::atomic::AtomicWriter;
    let vec_path = std::path::Path::new(&format!("{prefix}.eigenvec")).to_path_buf();
    let val_path = std::path::Path::new(&format!("{prefix}.eigenval")).to_path_buf();
    let pve_path = std::path::Path::new(&format!("{prefix}.pve")).to_path_buf();

    // doc_v28 V28-23: reed-io's AtomicWriter owns the create→flush→close
    // sequence; each file is finished here and published in a fixed order
    // below, so a failure leaves no partially written file behind.
    let mut vec_out = AtomicWriter::with_capacity(&vec_path, 1 << 16)?;
    {
        let f = &mut vec_out;
        // plink's header: `#FID IID PC1..PCk`; reed uses the sample name for
        // both FID and IID so single-name datasets round-trip.
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
        f.flush().map_err(io_err)?;
    }
    let vec_out = vec_out.finish()?;

    let mut val_out = AtomicWriter::with_capacity(&val_path, 1 << 16)?;
    {
        let f = &mut val_out;
        for v in values {
            writeln!(f, "{v:.9}").map_err(io_err)?;
        }
        f.flush().map_err(io_err)?;
    }
    let val_out = val_out.finish()?;

    let mut pve_out = AtomicWriter::with_capacity(&pve_path, 1 << 16)?;
    {
        let f = &mut pve_out;
        for v in pve {
            writeln!(f, "{v:.9}").map_err(io_err)?;
        }
        f.flush().map_err(io_err)?;
    }
    let pve_out = pve_out.finish()?;

    // Publish in a fixed order; any earlier rename succeeded, later ones
    // still happen, so no write failure leaves a *partially written* file.
    vec_out.commit()?;
    val_out.commit()?;
    pve_out.commit()?;
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
    fn write_fixture(
        dir: &std::path::Path,
        n_samples: usize,
        n_sites: usize,
    ) -> std::path::PathBuf {
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
        let mut deflater =
            reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6 is valid");
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

    /// doc_v26 V26-21: `solve` passes the budget's parallelism to the dense
    /// backend instead of letting faer read the process-global rayon pool.
    /// A serial budget must produce the same components as a parallel one —
    /// the bound changes thread usage, never the numbers.
    #[test]
    fn v26_21_solve_carries_the_thread_budget_and_keeps_the_numbers() {
        let n = 24;
        let m = 40;
        let x = Mat::from_fn(n, m, |i, j| {
            let v = ((i * 7 + j * 13) % 5) as f64;
            if v == 0.0 { f64::NAN } else { v - 2.0 }
        });
        let run = |threads: usize| {
            let budget = ThreadBudget {
                threads,
                mem: MemBudget::from_bytes(u64::MAX),
            };
            solve(&x, 4, Stand::Binom2, "meanimpute", "dense", &budget)
                .expect("invariant: dense solve")
        };
        let serial = run(1);
        let two = run(2);
        let many = run(8);
        for other in [&two, &many] {
            for (a, b) in serial.values.iter().zip(&other.values) {
                assert!(
                    (a - b).abs() <= 1e-12 * a.abs().max(1.0),
                    "eigenvalues must not depend on the thread count: {a} vs {b}"
                );
            }
            for (a, b) in serial.pve.iter().zip(&other.pve) {
                assert!((a - b).abs() <= 1e-12, "pve must not depend on threads");
            }
        }
        assert_eq!(serial.solver, "dense");
        // The stream path takes the same hint (the config's `parallel`).
        let budget = ThreadBudget {
            threads: 8,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        let streamed = solve(&x, 2, Stand::Binom2, "meanimpute", "stream", &budget)
            .expect("invariant: stream solve");
        assert_eq!(streamed.solver, "stream");
        assert_eq!(streamed.values.len(), 2);
    }

    /// doc_v25 V25-01 (F1): `reed pca` fails on a GT allele beyond the ALT
    /// count instead of silently turning the sample into missing data —
    /// `0/9` on a single-ALT record used to become NaN, and a file of such
    /// records ran PCA on an empty matrix.
    #[test]
    fn v25_01_pca_rejects_out_of_domain_gt() {
        let dir = temp_dir("v25-domain");
        let mut text = String::from(
            "##fileformat=VCFv4.3\n##contig=<ID=1,length=1000000>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n",
        );
        text.push_str("1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/9\t0/1\n");
        text.push_str("1\t101\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\n");
        let input = dir.join("in.vcf");
        std::fs::write(&input, text).expect("invariant: fixture");
        let request = PcaRequest {
            pcs: 1,
            prefix: dir.join("out").to_string_lossy().into_owned(),
            ..PcaRequest::default()
        };
        let budget = ThreadBudget {
            threads: 2,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        let Err(e) = run_pca(Some(&input), &request, &budget) else {
            panic!("an out-of-domain GT allele must fail the run (doc_v25 V25-01)")
        };
        assert_eq!(e.kind(), crate::ErrorKind::Runtime, "{e}");
        assert!(
            e.to_string().contains("exceeds ALT count 1"),
            "unexpected error: {e}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// doc_v25 V25-02: `FORMAT=DP:GT` reads the same dosages as `GT:DP`.
    /// The old first-colon lookup parsed the DP value (`9`) as the genotype,
    /// so every sample became NaN; the two runs' `.eigenvec` files differ.
    #[test]
    fn v25_02_pca_reads_gt_after_dp() {
        use std::fmt::Write as _;
        let dir = temp_dir("v25-gtpos");
        let make = |name: &str, format: &str, rows: &[&[&str]]| -> Vec<u8> {
            let mut text = String::from(
                "##fileformat=VCFv4.3\n##contig=<ID=1,length=1000000>\n\
                 ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
                 ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n",
            );
            for (i, row) in rows.iter().enumerate() {
                writeln!(
                    text,
                    "1\t{}\t.\tA\tG\t50\tPASS\t.\t{format}\t{}",
                    100 + i,
                    row.join("\t")
                )
                .expect("invariant: writing to a String");
            }
            let input = dir.join(format!("{name}.vcf"));
            std::fs::write(&input, text).expect("invariant: fixture");
            let request = PcaRequest {
                pcs: 2,
                prefix: dir.join(name).to_string_lossy().into_owned(),
                ..PcaRequest::default()
            };
            let budget = ThreadBudget {
                threads: 2,
                mem: MemBudget::from_bytes(u64::MAX),
            };
            run_pca(Some(&input), &request, &budget).expect("invariant: pca");
            std::fs::read(dir.join(format!("{name}.eigenvec"))).expect("invariant: eigenvec")
        };
        let gt_first = make(
            "gt-first",
            "GT:DP",
            &[
                &["0/0:9", "0/1:8", "1/1:7", "./.:."],
                &["0/1:9", "1/1:8", "0/1:7", "0/0:."],
                &["1/1:9", "0/0:8", "0/1:7", "0/0:."],
            ],
        );
        let dp_first = make(
            "dp-first",
            "DP:GT",
            &[
                &["9:0/0", "8:0/1", "7:1/1", ".:./."],
                &["9:0/1", "8:1/1", "7:0/1", ".:0/0"],
                &["9:1/1", "8:0/0", "7:0/1", ".:0/0"],
            ],
        );
        assert_eq!(gt_first, dp_first, "DP:GT must yield the same eigenvec");
        let _ = std::fs::remove_dir_all(&dir);
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
    fn m11_dense_paths_charge_the_grm_against_the_budget() {
        // doc_v21 M11: the old gate charged `2 x n x m` only, so an explicit
        // `--solver dense` passed it and then allocated an n x n GRM (two of
        // them for pairwise) that the budget never saw — measured 400 MiB
        // under a 2 MB budget on 2000 samples.
        let (n, m) = (2000_u64, 50_u64);
        let nm = n * m * 8;
        let n2 = n * n * 8;
        assert_eq!(MatrixUse::Plain.peak_bytes(n, m), 2 * nm);
        assert_eq!(MatrixUse::Stream.peak_bytes(n, m), 3 * nm);
        assert_eq!(MatrixUse::DenseMeanImpute.peak_bytes(n, m), 4 * nm + 6 * n2);
        assert_eq!(MatrixUse::DensePairwise.peak_bytes(n, m), 4 * nm + 7 * n2);
        // The dense classes must be strictly larger than the old gate.
        for class in [MatrixUse::DenseMeanImpute, MatrixUse::DensePairwise] {
            assert!(class.peak_bytes(n, m) > 2 * nm);
        }
        // Solver/missing resolve to the right class.
        let budget = ThreadBudget {
            threads: 1,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        let n_kept = usize::try_from(n).expect("invariant: fits");
        assert_eq!(chosen_solver("dense", n_kept, None, &budget), "dense");
        assert_eq!(chosen_solver("stream", n_kept, None, &budget), "stream");
        assert_eq!(
            matrix_use("dense", "meanimpute"),
            MatrixUse::DenseMeanImpute
        );
        assert_eq!(matrix_use("dense", "pairwise"), MatrixUse::DensePairwise);
        assert_eq!(matrix_use("stream", "meanimpute"), MatrixUse::Stream);
    }

    /// doc_v26 V26-20: the dense estimate must cover every object the solve
    /// holds at once (the matrices, `standardize`'s copy, faer's input,
    /// eigenvector matrix and workspace, and the returned copy), and `auto`
    /// must judge the dense path on that peak rather than on the GRM alone.
    #[test]
    fn v26_20_dense_peak_covers_the_whole_solve_and_gates_auto() {
        let (n, m) = (1000_u64, 2000_u64);
        let nm = n * m * 8;
        let n2 = n * n * 8;
        // The old estimate charged `2 nm + n2`; the real peak adds the
        // standardized copy and three more n x n arrays.
        assert!(
            MatrixUse::DenseMeanImpute.peak_bytes(n, m) >= 2 * nm + 4 * n2,
            "the estimate must cover st.w, faer's wm and am/u/stack/copy"
        );
        assert!(
            MatrixUse::DensePairwise.peak_bytes(n, m) > MatrixUse::DenseMeanImpute.peak_bytes(n, m)
        );

        // The old `auto` rule: dense when `n^2 * 8 <= budget`. Choose a
        // budget that passes the GRM test but fails the real peak; `auto`
        // must now pick stream.
        let grm = n * n * 8;
        let real = MatrixUse::DenseMeanImpute.peak_bytes(n, m);
        let budget = ThreadBudget {
            threads: 4,
            mem: MemBudget::from_bytes(grm),
        };
        let n_kept = usize::try_from(n).expect("invariant: fits");
        assert!(grm <= budget.mem.bytes && real > budget.mem.bytes);
        assert_eq!(
            chosen_solver(
                "auto",
                n_kept,
                Some(usize::try_from(m).expect("invariant")),
                &budget
            ),
            "stream",
            "auto must not pick a dense solve whose peak exceeds the budget"
        );
        // Before the site count is known the rule is the conservative square
        // form (used only by the pre-read refusal).
        assert_eq!(chosen_solver("auto", n_kept, None, &budget), "dense");
        // A budget that covers the real peak still picks dense.
        let roomy = ThreadBudget {
            threads: 4,
            mem: MemBudget::from_bytes(real),
        };
        assert_eq!(
            chosen_solver(
                "auto",
                n_kept,
                Some(usize::try_from(m).expect("invariant")),
                &roomy
            ),
            "dense"
        );
    }

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

    /// doc_v25 V25-10: the rank-deficient stream solve converges (it used
    /// to burn max_cycles) and reports its residual in the report.
    #[test]
    fn v25_10_stream_reports_a_converged_residual_on_a_rank_deficient_input() {
        let dir = temp_dir("v25-residual");
        let input = dir.join("in.vcf");
        // 20 samples, 5 sites: W W' has rank at most 5, so pcs=8 asks for
        // components beyond the rank.
        std::fs::write(&input, fixture_text(20, 5)).expect("invariant: fixture");
        let request = PcaRequest {
            pcs: 8,
            missing: "meanimpute".into(),
            solver: "stream".into(),
            prefix: dir.join("out").to_string_lossy().into_owned(),
            ..PcaRequest::default()
        };
        let budget = ThreadBudget {
            threads: 2,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        let report = run_pca(Some(&input), &request, &budget).expect("stream pca");
        assert_eq!(report.solver, "stream");
        let residual = report
            .stream_max_residual
            .expect("the stream solver records its residual");
        assert!(residual.is_finite(), "{residual}");
        assert!(
            residual <= 1e-10,
            "rank-deficient input did not converge: {residual}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}
