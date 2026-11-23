//! `reed distance`: the individual-by-individual matrix family (doc_v17).
//!
//! The CLI used to carry this pipeline itself (472 lines in
//! `reed-cli/src/translate/distance.rs`); doc_v30 P1-5 moved the matrix
//! selection, the bootstrap and the output publication behind the facade so
//! library callers get the same path and the `frontend` escape hatch can
//! stop exposing the lower crates wholesale.

use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use reed_distance::matrices::{
    AjkOptions, DiscordanceMode, DiscordanceOptions, DistanceOptions, MatrixResult, PiHatOptions,
    PiHatResult, ajk, allele_ct_distance, discordance, genome_columns, grm, ibs, king,
    one_minus_ibs, phi, pi_hat,
};
use reed_distance::output::{write_ids_to, write_matrix_to};
use reed_distance::standard::standardize;
use reed_io::write::atomic::{AtomicWriter, write_atomic};

use crate::context::RunContext;
use crate::error::Error;
use crate::options::{DistanceFormat, DistanceKind, MissingPolicy, Stand, TreeMethod};
use crate::pca::{MatrixUse, read_genotype_matrix};

/// Request for [`run_distance`].
#[derive(Clone, Debug)]
pub struct DistanceRequest {
    /// Which matrix or per-pair table to compute.
    pub kind: DistanceKind,
    /// Output layout for the matrix kinds. Per-pair tables have one layout,
    /// so they require the default (`Triangle`).
    pub format: DistanceFormat,
    /// Genotype standardization.
    pub stand: Stand,
    /// Missing-genotype policy for the GRM/covariance estimators.
    pub missing: MissingPolicy,
    /// PI_HAT: skip the `[0, 1]` clipping (plink `--genome unbounded`).
    pub unbounded: bool,
    /// PI_HAT: snap to a possible configuration when `PI_HAT² < Z2`.
    pub nudge: bool,
    /// Discordance: report the mismatch count instead of the probability.
    pub mismatch_count: bool,
    /// Discordance: phred-scaled genotyping error probability.
    pub error_probability: u32,
    /// Build a tree instead of writing the matrix; writes `PREFIX.nwk`.
    pub tree: Option<TreeMethod>,
    /// Bootstrap replicates for the tree; `None` skips the bootstrap.
    pub bootstrap: Option<usize>,
    /// Seed for the bootstrap replicates.
    pub bootstrap_seed: u64,
    /// Drop support labels below this percentage.
    pub bootstrap_cutoff: f64,
    /// Decimal digits in the newick output.
    pub digits: usize,
    /// Output prefix: writes `PREFIX.<suffix>`.
    pub prefix: String,
    /// Optional sample keep mask (header order); `None` keeps every sample.
    pub keep_samples: Option<Vec<bool>>,
}

impl Default for DistanceRequest {
    fn default() -> Self {
        Self {
            kind: DistanceKind::Grm,
            format: DistanceFormat::Triangle,
            stand: Stand::Binom2,
            missing: MissingPolicy::Pairwise,
            unbounded: false,
            nudge: false,
            mismatch_count: false,
            error_probability: 40,
            tree: None,
            bootstrap: None,
            bootstrap_seed: 42,
            bootstrap_cutoff: 0.0,
            digits: 10,
            prefix: String::new(),
            keep_samples: None,
        }
    }
}

/// Report of one distance run.
#[derive(Clone, Debug)]
pub struct DistanceReport {
    /// Samples in the matrix.
    pub n_samples: usize,
    /// Sites used.
    pub n_sites: usize,
    /// Which estimator ran.
    pub kind: DistanceKind,
    /// Output layout (matrix kinds only).
    pub format: DistanceFormat,
    /// Files written, in publication order.
    pub files: Vec<PathBuf>,
    /// Bootstrap support per internal node (tree runs only, ascending node
    /// id; the root is included, matching `ape::prop.clades`).
    pub bootstrap_support: Option<Vec<f64>>,
    /// Wall time.
    pub elapsed: Duration,
}

/// Runs `reed distance` over `path` (or stdin when `None`).
///
/// # Errors
///
/// Usage errors for a bad flag combination, and I/O, engine or graph
/// failures.
pub fn run_distance(
    path: Option<&Path>,
    request: &DistanceRequest,
    ctx: &RunContext<'_>,
) -> Result<DistanceReport, Error> {
    let opened = crate::input::open_regions(path, &[])?;
    run_distance_opened(opened, request, ctx)
}

/// Runs `reed distance` over an already-opened input.
///
/// # Errors
///
/// As [`run_distance`].
pub fn run_distance_opened(
    opened: crate::OpenedInput,
    request: &DistanceRequest,
    ctx: &RunContext<'_>,
) -> Result<DistanceReport, Error> {
    let started = Instant::now();
    validate(request)?;
    let budget = ctx.budget;
    // The estimators join whatever rayon pool is ambient; install the run's
    // budget-sized pool so `--threads` bounds the matrix and the bootstrap
    // (doc_v20 S37). The pool is installed even for one thread, because the
    // estimators are unconditionally parallel.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(budget.threads)
        .build()
        .map_err(|e| Error::usage(format!("--threads {}: {e}", budget.threads)))?;

    let (x, ids) = read_genotype_matrix(
        opened,
        request.keep_samples.as_deref(),
        ctx,
        MatrixUse::Plain,
    )?;
    if x.rows == 0 || x.cols == 0 {
        return Err(Error::usage("no biallelic sites with a diploid GT"));
    }
    let st = standardize(&x, request.stand.distance());

    let discordance_opts = DiscordanceOptions {
        mode: if request.mismatch_count {
            DiscordanceMode::Count
        } else {
            DiscordanceMode::Probability
        },
        error_probability: request.error_probability,
    };
    let (mut result, aux) = pool.install(|| match request.kind {
        DistanceKind::AlleleCountDistance => {
            (allele_ct_distance(&st, &DistanceOptions::default()), None)
        }
        DistanceKind::Ibs => (ibs(&st, &DistanceOptions::default()), None),
        DistanceKind::OneMinusIbs => (one_minus_ibs(&st, &DistanceOptions::default()), None),
        DistanceKind::Grm | DistanceKind::Covariance => (
            grm(
                &st,
                request.missing.distance(),
                request.kind == DistanceKind::Grm,
            ),
            None,
        ),
        DistanceKind::King => {
            let (m, _counts) = king(&st);
            (m, None)
        }
        DistanceKind::PiHat => {
            let opts = PiHatOptions {
                unbounded: request.unbounded,
                nudge: request.nudge,
            };
            let (m, _res) = pi_hat(&st, &opts);
            (m, None)
        }
        DistanceKind::Phi => {
            let (m, a) = phi(&st);
            (m, Some(a))
        }
        DistanceKind::Ajk => (ajk(&st, &AjkOptions::default()), None),
        DistanceKind::Discordance => {
            let (m, a) = discordance(&st, &discordance_opts);
            (m, Some(a))
        }
    });
    result.ids.clone_from(&ids);
    if aux.is_some() {
        result.aux = aux;
    }

    let mut files = Vec::new();
    let mut bootstrap_support = None;
    if let Some(method) = request.tree {
        // A tree needs a distance matrix: zero diagonal and non-negative
        // entries. Several kinds are similarities (GRM, KING, IBS, phi,
        // PI_HAT); building a tree from them would silently invent edge
        // weights, so refuse instead of guessing.
        for i in 0..result.n {
            let d = result.values[i * result.n + i];
            if d.abs() > 1e-9 {
                return Err(Error::usage(format!(
                    "--type {} is a similarity, not a distance (d({i},{i}) = {d:.6}); \
                     --tree needs a zero diagonal. Use --type 1-ibs, allele-ct, \
                     ajk or discordance",
                    request.kind.name()
                )));
            }
        }
        let mut d = crate::Mat::zeros(result.n, result.n);
        for i in 0..result.n {
            for j in 0..result.n {
                d.set(i, j, result.values[i * result.n + j]);
            }
        }
        let built = match method {
            TreeMethod::Upgma => reed_graph::upgma(&d),
            TreeMethod::Nj => reed_graph::nj(&d),
        };
        let mut tree = built.map_err(graph_error)?;
        tree.tip_labels.clone_from(&ids);

        if let Some(reps) = request.bootstrap {
            let opts = DistanceOptions::default();
            let stand = request.stand;
            let dist_fn = |m: &crate::Mat| -> crate::Mat {
                let st = standardize(m, stand.distance());
                let r = match request.kind {
                    DistanceKind::OneMinusIbs => one_minus_ibs(&st, &opts),
                    DistanceKind::Ajk => ajk(&st, &AjkOptions::default()),
                    DistanceKind::Discordance => discordance(&st, &discordance_opts).0,
                    // `allele-ct`, plus every kind the zero-diagonal gate
                    // above already rejected, falls back to allele counts.
                    _ => allele_ct_distance(&st, &opts),
                };
                let n = r.n;
                let mut out = crate::Mat::zeros(n, n);
                for i in 0..n {
                    for j in 0..n {
                        out.set(i, j, r.values[i * n + j]);
                    }
                }
                out
            };
            let cfg = reed_graph::BootstrapConfig {
                replicates: reps,
                seed: request.bootstrap_seed,
                cutoff: request.bootstrap_cutoff,
            };
            let report = pool
                .install(|| reed_graph::support(&tree, &x, method.lower(), dist_fn, &cfg))
                .map_err(graph_error)?;
            tree.node_labels = Some(
                report
                    .support
                    .iter()
                    .map(|v| {
                        if v.is_nan() {
                            String::new()
                        } else {
                            format!("{v:.0}")
                        }
                    })
                    .collect(),
            );
            bootstrap_support = Some(report.support);
        }

        let nwk = reed_graph::to_newick(&tree, request.digits).map_err(graph_error)?;
        let path = PathBuf::from(format!("{}.nwk", request.prefix));
        write_atomic(&path, nwk)?;
        files.push(path);
        return Ok(DistanceReport {
            n_samples: x.rows,
            n_sites: x.cols,
            kind: request.kind,
            format: request.format,
            files,
            bootstrap_support,
            elapsed: started.elapsed(),
        });
    }

    let suffix = request.kind.suffix();
    let out_path = PathBuf::from(format!("{}.{suffix}", request.prefix));
    let id_path = PathBuf::from(format!("{}.{suffix}.id", request.prefix));
    let mut matrix_out = AtomicWriter::with_capacity(&out_path, 1 << 16)?;
    let mut id_out = if request.kind.is_matrix() {
        Some(AtomicWriter::with_capacity(&id_path, 1 << 16)?)
    } else {
        None
    };
    if request.kind.is_matrix() {
        write_matrix_to(&mut matrix_out, &result, request.format.lower())?;
        matrix_out.flush()?;
        let matrix_out = matrix_out.finish()?;
        if let Some(mut id_out) = id_out.take() {
            write_ids_to(&mut id_out, &ids)?;
            id_out.flush()?;
            let id_out = id_out.finish()?;
            matrix_out.commit()?;
            id_out.commit()?;
            files.push(out_path);
            files.push(id_path);
        }
    } else {
        write_table_to(&mut matrix_out, &result, request.kind)?;
        matrix_out.flush()?;
        matrix_out.finish()?.commit()?;
        files.push(out_path);
    }

    Ok(DistanceReport {
        n_samples: x.rows,
        n_sites: x.cols,
        kind: request.kind,
        format: request.format,
        files,
        bootstrap_support,
        elapsed: started.elapsed(),
    })
}

fn validate(request: &DistanceRequest) -> Result<(), Error> {
    if request.prefix == "-" {
        return Err(Error::usage(
            "distance writes PREFIX files: -o - is not stdout; give a file prefix",
        ));
    }
    if request.prefix.trim().is_empty() {
        return Err(Error::usage(
            "-o needs a non-empty prefix (an empty one writes hidden files in the working directory)",
        ));
    }
    if !request.kind.is_matrix() && request.format != DistanceFormat::Triangle {
        return Err(Error::usage(format!(
            "--format applies to matrix --type values; {} is a per-pair table",
            request.kind.name()
        )));
    }
    if request.tree.is_some() && request.digits == 0 {
        return Err(Error::usage(
            "--digits must be at least 1; 0 rounds every branch length to zero",
        ));
    }
    if request.bootstrap.is_some_and(|reps| reps == 0) {
        return Err(Error::usage("--bootstrap must be at least 1"));
    }
    Ok(())
}

/// Exit-code classification for graph failures (doc_v25 V25-08): a request
/// the command line can fix is usage; a matrix the data made unusable
/// (NaN/asymmetric distances, non-finite branches, too-deep trees) is
/// runtime, so a retry loop does not spin on exit 2.
fn graph_error(e: reed_graph::GraphError) -> Error {
    use reed_graph::GraphError as G;
    match e {
        G::NotSquare { .. } | G::TooFewTaxa { .. } => Error::usage(e.to_string()),
        _ => Error::foreign(e),
    }
}

/// Writes `NA` for a missing `.genome` column value, or the formatted number.
fn write_na_or<W: std::io::Write>(f: &mut W, value: Option<String>) -> Result<(), Error> {
    match value {
        Some(text) => f.write_all(text.as_bytes()).map_err(Error::from),
        None => f.write_all(b"NA").map_err(Error::from),
    }
}

/// Writes the per-pair table formats (`.genome`, `.relatedness2`,
/// `.relatedness`, `.gtcheck`) with the reference tools' exact headers.
fn write_table_to<W: std::io::Write>(
    f: &mut W,
    m: &MatrixResult,
    kind: DistanceKind,
) -> Result<(), Error> {
    let n = m.n;
    match kind {
        DistanceKind::Phi => {
            writeln!(
                f,
                "INDV1\tINDV2\tN_AaAa\tN_AAaa\tN1_Aa\tN2_Aa\tRELATEDNESS_PHI"
            )?;
            let aux = m.aux.as_ref();
            // vcftools emits the full ordered matrix, both directions.
            for i in 0..n {
                for j in 0..n {
                    let hethet = aux.map_or(0, |a| a.n_hethet[i * n + j]);
                    let opp = aux.map_or(0, |a| a.n_opposite_hom[i * n + j]);
                    let phi_v = m.values[i * n + j];
                    let n1 = aux.map_or(0, |a| a.n_het.get(i).copied().unwrap_or(0));
                    let n2 = aux.map_or(0, |a| a.n_het.get(j).copied().unwrap_or(0));
                    writeln!(
                        f,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{:.9}",
                        m.ids[i], m.ids[j], hethet, opp, n1, n2, phi_v
                    )?;
                }
            }
        }
        DistanceKind::Ajk => {
            writeln!(f, "INDV1\tINDV2\tRELATEDNESS_AJK")?;
            for i in 0..n {
                for j in i..n {
                    writeln!(f, "{}\t{}\t{:.9}", m.ids[i], m.ids[j], m.values[i * n + j])?;
                }
            }
        }
        DistanceKind::Discordance => {
            writeln!(
                f,
                "#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t\
                 [5]Average -log P(HWE)\t[6]Number of sites compared\t\
                 [7]Number of matching genotypes"
            )?;
            let aux = m.aux.as_ref();
            for i in 0..n {
                for j in 0..i {
                    let cmp = aux.map_or(0, |a| a.n_compared[i * n + j]);
                    let mt = aux.map_or(0, |a| a.n_matching[i * n + j]);
                    writeln!(
                        f,
                        "DCv2\t{}\t{}\t{:.9e}\t0\t{}\t{}",
                        m.ids[i],
                        m.ids[j],
                        m.values[i * n + j],
                        cmp,
                        mt
                    )?;
                }
            }
        }
        DistanceKind::PiHat => {
            writeln!(
                f,
                "FID1\tIID1\tFID2\tIID2\tRT\tEZ\tZ0\tZ1\tZ2\tPI_HAT\tPHE\tDST\tPPC\tRATIO"
            )?;
            let aux = m.aux.as_ref();
            // plink writes FID and IID as separate columns; reed's ids are
            // single sample names. Writing the name in both columns keeps a
            // name containing `_` intact instead of guessing a FID/IID
            // boundary (doc_v25 V25-08, R14).
            let split = |id: &str| -> (String, String) { (id.to_string(), id.to_string()) };
            for i in 0..n {
                for j in 0..i {
                    let (f1, i1) = split(&m.ids[i]);
                    let (f2, i2) = split(&m.ids[j]);
                    let z0 = aux.map_or(f64::NAN, |a| a.z0[i * n + j]);
                    let z1 = aux.map_or(f64::NAN, |a| a.z1[i * n + j]);
                    let z2 = aux.map_or(f64::NAN, |a| a.z2[i * n + j]);
                    // DST follows plink's formula; PPC/RATIO stay NA (see
                    // `genome_columns` for why they are not reproducible).
                    let cols = aux.map(|a| {
                        genome_columns(&PiHatResult {
                            z0: a.z0[i * n + j],
                            z1: a.z1[i * n + j],
                            z2: a.z2[i * n + j],
                            pi_hat: m.values[i * n + j],
                            n_valid: a.n_valid[i * n + j],
                            ibs1: a.n_ibs1[i * n + j],
                            ibs0: a.n_ibs0[i * n + j],
                        })
                    });
                    let (dst, ppc, ratio) = cols.unwrap_or((None, None, None));
                    write!(
                        f,
                        "{f1}\t{i1}\t{f2}\t{i2}\tOT\tNA\t{z0:.4}\t{z1:.4}\t{z2:.4}\t{:.4}\t-9\t",
                        m.values[i * n + j]
                    )?;
                    write_na_or(f, dst.map(|v| format!("{v:.6}")))?;
                    write!(f, "\t")?;
                    write_na_or(f, ppc.map(|v| format!("{v:.4}")))?;
                    write!(f, "\t")?;
                    write_na_or(f, ratio.map(|v| format!("{v:.4}")))?;
                    writeln!(f)?;
                }
            }
        }
        _ => unreachable!("is_matrix covers the rest"),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorKind;

    #[test]
    fn graph_failures_split_into_usage_and_runtime() {
        use reed_graph::GraphError;
        let usage = [
            GraphError::TooFewTaxa { got: 1, need: 3 },
            GraphError::NotSquare { rows: 2, cols: 3 },
        ];
        for err in usage {
            assert_eq!(graph_error(err.clone()).kind(), ErrorKind::Usage, "{err}");
        }
        let runtime = [
            GraphError::NonFiniteDistance { row: 0, col: 1 },
            GraphError::Asymmetric { row: 0, col: 1 },
            GraphError::NonZeroDiagonal { index: 0 },
            GraphError::BadBranchLength { edge: 0 },
            GraphError::TreeTooDeep { depth: 10_001 },
        ];
        for err in runtime {
            assert_eq!(graph_error(err.clone()).kind(), ErrorKind::Runtime, "{err}");
        }
    }

    #[test]
    fn a_per_pair_table_rejects_a_matrix_format() {
        let request = DistanceRequest {
            kind: DistanceKind::PiHat,
            format: DistanceFormat::Square,
            prefix: "out".into(),
            ..DistanceRequest::default()
        };
        let err = validate(&request).expect_err("format must be rejected");
        assert_eq!(err.kind(), ErrorKind::Usage);
    }
}
