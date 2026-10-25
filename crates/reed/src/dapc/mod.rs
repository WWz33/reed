//! `reed dapc`: discriminant analysis of principal components.
//!
//! The pipeline follows `adegenet::dapc`: standardize the genotypes, take the
//! leading principal components, then run a Fisher linear discriminant
//! analysis on those components using the supplied population labels. The
//! discriminant model and its posterior probabilities come from
//! [`reed_reduce::lda`], which matches `MASS::lda` (the routine `adegenet`
//! calls).
//!
//! # Scale note
//!
//! `ade4::dudi.pca` reports scores as `U · d` while reed reports `U · sqrt(λ)`.
//! The two differ by a constant factor per axis. Discriminant *directions* and
//! posterior probabilities are invariant to that factor, so `posterior` and
//! `assign` agree with `adegenet` bit for bit; the discriminant coordinates
//! (`LD1..LDk`) differ by a per-axis scale. See
//! `docs/rounds/v17b/01-research.md` §1.5.

use std::io::Write as _;

use reed_engine::budget::ThreadBudget;
use reed_reduce::lda::{self, LdaError};
use reed_reduce::pca::Stand;
use reed_reduce::Mat;

use crate::pca::read_genotype_matrix;

/// Request for [`run_dapc`].
#[derive(Clone, Debug)]
pub struct DapcRequest {
    /// Principal components to retain; `None` keeps `min(n_samples - 1, 50)`.
    pub n_pca: Option<usize>,
    /// Discriminant axes to retain; `None` keeps every available axis.
    pub n_da: Option<usize>,
    /// `sample<TAB>population` rows (no header).
    pub populations: String,
    /// `binom2` (default) or `binom`.
    pub stand: String,
    /// `pairwise` (default) or `meanimpute`.
    pub missing: String,
    /// `auto`, `dense` or `stream`.
    pub solver: String,
    /// Output prefix: writes `.dapc.ind`, `.dapc.eig`, `.dapc.var`.
    pub prefix: String,
    /// Optional sample keep mask (header order).
    pub keep_samples: Option<Vec<bool>>,
    /// Skip population-file rows naming samples absent from the header.
    pub force_samples: bool,
    /// Worker threads.
    pub threads: Option<usize>,
    /// Memory budget in MiB.
    pub memory: Option<u64>,
    /// Run cross-validation instead of a single DAPC. `Some` writes
    /// `PREFIX.dapc.xval` and the PCA/LDA files for the selected `n_pca`.
    pub xval: Option<XvalRequest>,
}

/// Cross-validation settings (`--xval`).
#[derive(Clone, Debug)]
pub struct XvalRequest {
    /// Repetitions per candidate component count.
    pub n_rep: usize,
    /// Training fraction (reduced automatically for small groups).
    pub training_set: f64,
    /// `overall` or `group-mean`.
    pub result: String,
    /// Seed for the per-repetition sub-streams.
    pub seed: u64,
    /// Candidate component counts; `None` uses the `adegenet` rule.
    pub n_pca: Option<Vec<usize>>,
}

impl Default for XvalRequest {
    fn default() -> Self {
        Self {
            n_rep: 30,
            training_set: 0.9,
            result: "overall".into(),
            seed: 42,
            n_pca: None,
        }
    }
}

impl Default for DapcRequest {
    fn default() -> Self {
        Self {
            n_pca: None,
            n_da: None,
            populations: String::new(),
            stand: "binom2".into(),
            missing: "pairwise".into(),
            solver: "auto".into(),
            prefix: "out".into(),
            keep_samples: None,
            force_samples: false,
            threads: None,
            memory: None,
            xval: None,
        }
    }
}

/// Report of one DAPC run.
#[derive(Clone, Debug)]
pub struct DapcReport {
    /// Samples analysed.
    pub n_samples: usize,
    /// Sites used.
    pub n_sites: usize,
    /// Principal components retained.
    pub n_pca: usize,
    /// Discriminant axes written.
    pub n_da: usize,
    /// Number of populations.
    pub n_groups: usize,
    /// Variance fraction explained by the retained PCs.
    pub var: f64,
    /// Discriminant eigenvalues (`svd²`), descending.
    pub eig: Vec<f64>,
    /// Wall time.
    pub elapsed: std::time::Duration,
    /// Cross-validation outcome, when `--xval` was requested.
    pub xval: Option<XvalOutcome>,
}

/// The cross-validation result attached to a [`DapcReport`].
#[derive(Clone, Debug)]
pub struct XvalOutcome {
    /// `(n_pca, mean success)` per candidate, ascending.
    pub candidates: Vec<(usize, f64)>,
    /// The selected component count (lowest error from a perfect score).
    pub best_n_pca: usize,
    /// The component count with the highest mean success.
    pub best_mean_success_n_pca: usize,
    /// Accuracy expected from guessing by prior.
    pub random_chance: f64,
}

/// Runs DAPC over `path` (or stdin when `None`).
///
/// # Errors
///
/// Usage errors for a bad flag combination or a malformed populations file,
/// and I/O, engine or LDA failures.
pub fn run_dapc(
    path: Option<&std::path::Path>,
    request: &DapcRequest,
    budget: &ThreadBudget,
) -> Result<DapcReport, crate::Error> {
    let opened = crate::open_regions(path, &[])?;
    run_dapc_opened(opened, request, budget)
}

/// Runs DAPC over an already-opened input.
///
/// # Errors
///
/// As [`run_dapc`].
pub fn run_dapc_opened(
    opened: crate::OpenedInput,
    request: &DapcRequest,
    budget: &ThreadBudget,
) -> Result<DapcReport, crate::Error> {
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
    if request.populations.is_empty() {
        return Err(crate::Error::usage("--pop is required"));
    }

    // Group labels must be read against the *full* header before the keep
    // mask narrows the matrix, so the population file is validated against
    // every sample the header declares.
    let header_names = opened.facts().sample_names.clone();
    if header_names.is_empty() {
        return Err(crate::Error::usage(
            "the header declares no samples; nothing to compute",
        ));
    }
    let (group_of, group_names) = crate::stats::request::assign_groups(
        std::path::Path::new(&request.populations),
        &header_names,
        request.force_samples,
    )?;

    let keep = request.keep_samples.clone();
    let (x, ids) = read_genotype_matrix(opened, keep.as_deref(), budget)?;
    let n_samples = x.rows;
    let n_sites = x.cols;
    if n_samples < 3 {
        return Err(crate::Error::usage(format!(
            "DAPC needs at least 3 samples, got {n_samples}"
        )));
    }

    // Remap the header-order labels onto the kept rows.
    let kept_header: Vec<usize> = match &keep {
        Some(mask) => (0..mask.len()).filter(|&i| mask[i]).collect(),
        None => (0..header_names.len()).collect(),
    };
    let mut groups: Vec<u32> = Vec::with_capacity(n_samples);
    for &header_idx in &kept_header {
        match group_of.get(header_idx).copied().flatten() {
            Some(g) => groups.push(g),
            None => {
                return Err(crate::Error::usage(format!(
                    "sample {:?} is not in the populations file",
                    header_names[header_idx]
                )));
            }
        }
    }
    debug_assert_eq!(groups.len(), n_samples);

    // Drop groups that lost every member to the keep mask, then compact the
    // remaining indices so `lda::fit` sees a dense 0..ng range.
    let mut used = vec![false; group_names.len()];
    for &g in &groups {
        used[g as usize] = true;
    }
    let mut remap = vec![u32::MAX; group_names.len()];
    let mut levels = Vec::new();
    for (old, &u) in used.iter().enumerate() {
        if u {
            remap[old] = u32::try_from(levels.len())
                .map_err(|_| crate::Error::usage("more than u32::MAX populations"))?;
            levels.push(group_names[old].clone());
        }
    }
    for g in &mut groups {
        *g = remap[*g as usize];
    }
    if levels.len() < 2 {
        return Err(crate::Error::usage(format!(
            "DAPC needs at least 2 populations, got {}",
            levels.len()
        )));
    }

    // PCA on the standardized genotype matrix, through the same solver the
    // plain `reed pca` uses.
    let max_pcs = n_samples.saturating_sub(1).max(1);
    let want = request
        .n_pca
        .unwrap_or_else(|| max_pcs.min(50))
        .min(max_pcs);
    if want == 0 {
        return Err(crate::Error::usage("--n-pca must be at least 1"));
    }

    let started = std::time::Instant::now();
    let c = crate::pca::solve(&x, want, stand, missing, solver, budget)?;
    let (values, vectors, pve) = (c.values, c.vectors, c.pve);

    let n_pca = values.len();
    if n_pca == 0 {
        return Err(crate::Error::usage("PCA produced no components"));
    }
    // Keep only components with a positive eigenvalue, as `adegenet` does
    // (`X.rank <- sum(pcaX$eig > 1e-14)`).
    let rank = values.iter().filter(|&&v| v > 1e-14).count().max(1);
    let n_pca = n_pca.min(rank);

    // The discriminant analysis runs on the PC scores.
    let mut scores = Mat::zeros(n_samples, n_pca);
    for i in 0..n_samples {
        for (k, &lambda) in values.iter().take(n_pca).enumerate() {
            scores.set(i, k, vectors.get(i, k) * lambda.sqrt());
        }
    }

    // Cross-validation replaces the single fit: it selects n_pca, then the
    // final model is fitted with that many components.
    let mut xval_outcome = None;
    let n_pca = if let Some(xr) = &request.xval {
        let cfg = reed_reduce::xval::XvalConfig {
            training_set: xr.training_set,
            n_rep: xr.n_rep,
            result: reed_reduce::xval::XvalResult::parse(&xr.result)
                .map_err(crate::Error::usage)?,
            seed: xr.seed,
        };
        let candidates = match &xr.n_pca {
            Some(list) => list.clone(),
            None => reed_reduce::xval::candidate_counts(n_pca),
        };
        if candidates.is_empty() {
            return Err(crate::Error::usage(
                "cross-validation has no candidate component counts; \
                 pass --n-pca or use more samples",
            ));
        }
        let report = reed_reduce::xval::run(&scores, &groups, &levels, &candidates, &cfg)
            .map_err(|e| crate::Error::usage(e.to_string()))?;
        xval_outcome = Some(XvalOutcome {
            candidates: report
                .candidates
                .iter()
                .map(|c| (c.n_pca, c.mean_success()))
                .collect(),
            best_n_pca: report.best_n_pca,
            best_mean_success_n_pca: report.best_mean_success_n_pca,
            random_chance: report.random_chance,
        });
        report.best_n_pca.min(n_pca)
    } else {
        n_pca
    };

    // Re-slice the scores to the selected component count.
    let mut scores = if n_pca < scores.cols {
        let mut s = Mat::zeros(n_samples, n_pca);
        for i in 0..n_samples {
            for c in 0..n_pca {
                s.set(i, c, scores.get(i, c));
            }
        }
        s
    } else {
        scores
    };
    let _ = &mut scores;

    let model = lda::fit(&scores, &groups, &levels, 1e-30).map_err(map_lda_err)?;
    // `adegenet` keeps only the axes whose squared singular value exceeds
    // 1e-10 (`lda.dim <- sum(ldaX$svd^2 > 1e-10)`); the rest are numerical
    // null directions produced by collinear class means.
    let lda_dim = model.svd.iter().filter(|&&v| v * v > 1e-10).count().max(1);
    let n_da = request.n_da.unwrap_or(lda_dim).clamp(1, lda_dim);
    let projection = model.project(&scores, n_da).map_err(map_lda_err)?;

    // Share of the whole spectrum carried by the retained PCs — `pve` is
    // already normalised against the total, so this stays correct on the
    // stream path where only the top-k eigenvalues exist.
    let var: f64 = pve.iter().take(n_pca).sum();
    let eig: Vec<f64> = model.svd.iter().take(lda_dim).map(|v| v * v).collect();

    if let Some(xv) = &xval_outcome {
        write_xval(&request.prefix, xv)?;
    }

    write_outputs(&Outputs {
        prefix: &request.prefix,
        ids: &ids,
        levels: &levels,
        groups: &groups,
        projection: &projection,
        eig: &eig,
        n_pca,
        n_da,
        var,
    })?;

    Ok(DapcReport {
        n_samples,
        n_sites,
        n_pca,
        n_da,
        n_groups: levels.len(),
        var,
        eig,
        elapsed: started.elapsed(),
        xval: xval_outcome,
    })
}

fn map_lda_err(e: LdaError) -> crate::Error {
    crate::Error::usage(e.to_string())
}

/// Writes `.dapc.ind`, `.dapc.eig` and `.dapc.var` next to `prefix`.
struct Outputs<'a> {
    prefix: &'a str,
    ids: &'a [String],
    levels: &'a [String],
    groups: &'a [u32],
    projection: &'a reed_reduce::lda::LdaProjection,
    eig: &'a [f64],
    n_pca: usize,
    n_da: usize,
    var: f64,
}

fn write_outputs(o: &Outputs<'_>) -> Result<(), crate::Error> {
    let Outputs {
        prefix,
        ids,
        levels,
        groups,
        projection: p,
        eig,
        n_pca,
        n_da,
        var,
    } = *o;
    let ind_path = format!("{prefix}.dapc.ind");
    let mut f = std::fs::File::create(&ind_path)?;
    write!(f, "#IID\tPOP\tASSIGN").map_err(io_err)?;
    for k in 0..n_da {
        write!(f, "\tLD{}", k + 1).map_err(io_err)?;
    }
    writeln!(f).map_err(io_err)?;
    for i in 0..ids.len() {
        write!(
            f,
            "{}\t{}\t{}",
            ids[i], levels[groups[i] as usize], levels[p.assign[i] as usize]
        )
        .map_err(io_err)?;
        for k in 0..n_da {
            write!(f, "\t{:.9}", p.coords.get(i, k)).map_err(io_err)?;
        }
        writeln!(f).map_err(io_err)?;
    }

    let mut f = std::fs::File::create(format!("{prefix}.dapc.eig"))?;
    for v in eig {
        writeln!(f, "{v:.9}").map_err(io_err)?;
    }

    let mut f = std::fs::File::create(format!("{prefix}.dapc.var"))?;
    writeln!(f, "#n_pca\tn_da\tvar").map_err(io_err)?;
    writeln!(f, "{n_pca}\t{n_da}\t{var:.9}").map_err(io_err)?;
    Ok(())
}

/// Writes `PREFIX.dapc.xval`: one row per candidate plus the summary lines.
fn write_xval(prefix: &str, xv: &XvalOutcome) -> Result<(), crate::Error> {
    let mut f = std::fs::File::create(format!("{prefix}.dapc.xval"))?;
    writeln!(f, "#n_pca\tmean_success").map_err(io_err)?;
    for (n, s) in &xv.candidates {
        writeln!(f, "{n}\t{s:.9}").map_err(io_err)?;
    }
    writeln!(f, "#best_n_pca\t{}", xv.best_n_pca).map_err(io_err)?;
    writeln!(f, "#best_mean_success_n_pca\t{}", xv.best_mean_success_n_pca)
        .map_err(io_err)?;
    writeln!(f, "#random_chance\t{:.9}", xv.random_chance).map_err(io_err)?;
    Ok(())
}

fn io_err(e: std::io::Error) -> crate::Error {
    crate::Error::from(e)
}
