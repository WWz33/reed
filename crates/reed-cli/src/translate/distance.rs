//! `reed distance`: the individual-by-individual matrix family (doc_v17).

use std::io::Write as _;

use reed::budget::ThreadBudget;
use reed_distance::matrices::{
    ajk, allele_ct_distance, discordance, grm, ibs, king, one_minus_ibs, phi, pi_hat,
    AjkOptions, DistanceOptions, DiscordanceMode, DiscordanceOptions, MatrixKind, MatrixResult,
    MissingPolicy, PiHatOptions,
};
use reed_distance::output::{is_matrix_kind, write_ids, write_matrix};
use reed_distance::{standardize, MatrixFormat, Stand};
use reed_graph as graph;

use crate::cli::DistanceCli;

/// Which matrix the user asked for.
fn parse_kind(s: &str) -> Result<MatrixKind, String> {
    Ok(match s {
        "allele-ct" | "dist" => MatrixKind::AlleleCountDistance,
        "ibs" => MatrixKind::Ibs,
        "1-ibs" => MatrixKind::OneMinusIbs,
        "grm" | "rel" => MatrixKind::Grm,
        "cov" => MatrixKind::Covariance,
        "king" => MatrixKind::King,
        "pihat" | "genome" => MatrixKind::PiHat,
        "phi" | "relatedness2" => MatrixKind::Phi,
        "ajk" | "relatedness" => MatrixKind::Ajk,
        "discordance" | "gtcheck" => MatrixKind::Discordance,
        other => {
            return Err(format!(
                "unknown --type {other:?} \
                 (allele-ct|ibs|1-ibs|grm|cov|king|pihat|phi|ajk|discordance)"
            ));
        }
    })
}

/// Runs `reed distance`.
///
/// # Errors
///
/// Usage errors for bad flags, and I/O or engine failures.
pub(crate) fn run_distance_cli(args: DistanceCli) -> Result<(), reed::Error> {
    let kind = parse_kind(&args.matrix_type).map_err(reed::Error::usage)?;
    let format = MatrixFormat::parse(&args.format).map_err(reed::Error::usage)?;
    let stand = Stand::parse(&args.stand).map_err(reed::Error::usage)?;
    let Some(prefix) = args.output.as_deref() else {
        return Err(reed::Error::usage(
            "-o PREFIX is required: distance writes PREFIX.<suffix>",
        ));
    };

    let path = crate::units::input_path(&args.input);
    let opened = reed::open_regions(path.as_deref(), &[])?;
    let sample_names = opened.facts().sample_names.clone();
    if sample_names.is_empty() {
        return Err(reed::Error::usage(
            "the header declares no samples; nothing to compute",
        ));
    }
    let keep = super::common::parse_sample_keep(&args.samples, &sample_names)
        .map_err(reed::Error::usage)?;

    let threads = args.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .expect("invariant: parallelism")
            .get()
    });
    let budget = ThreadBudget {
        threads,
        mem: args
            .memory
            .map_or_else(reed::budget::MemBudget::default_half_memory, |m| {
                reed::budget::MemBudget::from_bytes(m.saturating_mul(1 << 20))
            }),
    };

    // doc_v20 S37: R13's reed-distance estimators join whatever rayon pool
    // is ambient, so `--threads` takes effect only through a scoped pool
    // installed around the estimator and bootstrap calls — the same
    // throwaway-pool pattern the stats/filter stages use via
    // `run_auxiliary`; unlike those, the estimators are unconditionally
    // parallel, so the pool is installed even for a single thread.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| reed::Error::usage(format!("--threads {threads}: {e}")))?;

    let (x, ids) = reed::read_genotype_matrix(opened, keep.as_deref(), &budget)?;
    if x.rows == 0 || x.cols == 0 {
        return Err(reed::Error::usage("no biallelic sites with a diploid GT"));
    }
    let st = standardize(&x, stand);

    // doc_v20 S37: built once from the user's flags so the primary matrix
    // and every bootstrap replicate below score the same metric — the
    // replicate closure used to pass `DiscordanceOptions::default()`
    // (probability mode, error 40), silently ignoring
    // `--mismatch-count`/`--error-probability` and reporting support for a
    // different metric than the tree.
    let discordance_opts = DiscordanceOptions {
        mode: if args.mismatch_count {
            DiscordanceMode::Count
        } else {
            DiscordanceMode::Probability
        },
        error_probability: args.error_probability,
    };

    let (mut result, aux) = pool.install(|| match kind {
        MatrixKind::AlleleCountDistance => {
            (allele_ct_distance(&st, &DistanceOptions::default()), None)
        }
        MatrixKind::Ibs => (ibs(&st, &DistanceOptions::default()), None),
        MatrixKind::OneMinusIbs => (one_minus_ibs(&st, &DistanceOptions::default()), None),
        MatrixKind::Grm | MatrixKind::Covariance => {
            let policy = if args.mean_impute {
                MissingPolicy::MeanImpute
            } else {
                MissingPolicy::Pairwise
            };
            (grm(&st, policy, kind == MatrixKind::Grm), None)
        }
        MatrixKind::King => {
            let (m, _counts) = king(&st);
            (m, None)
        }
        MatrixKind::PiHat => {
            let opts = PiHatOptions { unbounded: args.unbounded, nudge: args.nudge };
            let (m, _res) = pi_hat(&st, &opts);
            (m, None)
        }
        MatrixKind::Phi => {
            let (m, a) = phi(&st);
            (m, Some(a))
        }
        MatrixKind::Ajk => (ajk(&st, &AjkOptions::default()), None),
        MatrixKind::Discordance => {
            let (m, a) = discordance(&st, &discordance_opts);
            (m, Some(a))
        }
    });
    result.ids.clone_from(&ids);
    // `pi_hat` fills its own aux (the Z triplet); only overwrite when the
    // estimator returned one explicitly.
    if aux.is_some() {
        result.aux = aux;
    }

    // `--tree` builds a tree from the matrix and writes PREFIX.nwk instead
    // of the matrix file itself.
    if let Some(method) = args.tree.as_deref() {
        // A tree needs a distance matrix: d(i,i) = 0 and d(i,j) >= 0. Several
        // `--type` values are *similarities* instead (GRM, KING, IBS, phi,
        // pihat), whose diagonal is the self-relatedness rather than zero.
        // Building a tree from those without a conversion would silently
        // produce a tree over meaningless edge weights, so refuse instead of
        // guessing which transform the user meant.
        for i in 0..result.n {
            let d = result.values[i * result.n + i];
            if d.abs() > 1e-9 {
                return Err(reed::Error::usage(format!(
                    "--type {} is a similarity, not a distance (d({i},{i}) = {d:.6}); \
                     --tree needs a zero diagonal. Use --type 1-ibs, allele-ct, \
                     ajk or discordance",
                    args.matrix_type
                )));
            }
        }
        if args.digits == 0 {
            return Err(reed::Error::usage(
                "--digits must be at least 1; 0 rounds every branch length to zero",
            ));
        }
        let mut d = reed::Mat::zeros(result.n, result.n);
        for i in 0..result.n {
            for j in 0..result.n {
                d.set(i, j, result.values[i * result.n + j]);
            }
        }
        let built = match method {
            "upgma" => graph::upgma(&d),
            "nj" => graph::nj(&d),
            other => {
                return Err(reed::Error::usage(format!(
                    "unknown --tree {other:?} (upgma|nj)"
                )));
            }
        };
        let tree = built.map_err(|e| reed::Error::usage(e.to_string()))?;
        let mut tree = tree;
        tree.tip_labels.clone_from(&ids);

        if let Some(reps) = args.bootstrap {
            if reps == 0 {
                return Err(reed::Error::usage("--bootstrap must be at least 1"));
            }
            let tree_method = graph::TreeMethod::parse(method)
                .map_err(reed::Error::usage)?;
            // The bootstrap rebuilds the matrix from resampled loci, so it
            // needs the same matrix kind the tree used. Capture the kind and
            // the options; the closure re-runs the estimator on each
            // replicate. doc_v20 S37: `discordance_opts` is the very struct
            // the primary matrix above was built from, so the replicates
            // score the user's metric, not the defaults.
            let opts = DistanceOptions::default();
            let dist_fn = |m: &reed::Mat| -> reed::Mat {
                let st = standardize(m, stand);
                let r = match kind {
                    MatrixKind::OneMinusIbs => one_minus_ibs(&st, &opts),
                    MatrixKind::Ajk => ajk(&st, &AjkOptions::default()),
                    MatrixKind::Discordance => discordance(&st, &discordance_opts).0,
                    // `allele-ct`, plus every kind the zero-diagonal gate
                    // above already rejected, falls back to allele counts.
                    _ => allele_ct_distance(&st, &opts),
                };
                let n = r.n;
                let mut out = reed::Mat::zeros(n, n);
                for i in 0..n {
                    for j in 0..n {
                        out.set(i, j, r.values[i * n + j]);
                    }
                }
                out
            };
            let cfg = graph::BootstrapConfig {
                replicates: reps,
                seed: args.bootstrap_seed,
                cutoff: args.bootstrap_cutoff,
            };
            // doc_v20 S37: the replicates' estimator calls join the scoped
            // pool too, so `--threads` bounds the bootstrap as well.
            let report = pool
                .install(|| graph::support(&tree, &x, tree_method, dist_fn, &cfg))
                .map_err(|e| reed::Error::usage(e.to_string()))?;
            // Support is per internal node, ascending id; `to_newick` writes
            // node labels in that order.
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
            eprintln!(
                "reed distance --tree --bootstrap: {} replicates, {} internal nodes",
                report.replicates,
                report.support.len()
            );
        }

        let nwk = graph::to_newick(&tree, args.digits)
            .map_err(|e| reed::Error::usage(e.to_string()))?;
        std::fs::write(format!("{prefix}.nwk"), nwk).map_err(reed::Error::from)?;
        eprintln!(
            "reed distance --tree: {} samples, {} sites, method={method}",
            x.rows, x.cols
        );
        return Ok(());
    }

    let suffix = kind.suffix();
    let out_path = std::path::PathBuf::from(format!("{prefix}.{suffix}"));

    if is_matrix_kind(kind) {
        write_matrix(&out_path, &result, format).map_err(reed::Error::from)?;
        write_ids(
            &std::path::PathBuf::from(format!("{prefix}.{suffix}.id")),
            &ids,
        )
        .map_err(reed::Error::from)?;
    } else {
        write_table(&out_path, &result)?;
    }

    eprintln!(
        "reed distance: {} samples, {} sites, type={}, format={}",
        x.rows,
        x.cols,
        args.matrix_type,
        format.name()
    );
    Ok(())
}

/// Writes the per-pair table formats (`.genome`, `.relatedness2`,
/// `.relatedness`, `.gtcheck`) with the reference tools' exact headers.
fn write_table(path: &std::path::Path, m: &MatrixResult) -> Result<(), reed::Error> {
    let mut f = std::fs::File::create(path).map_err(reed::Error::from)?;
    let n = m.n;
    match m.kind {
        MatrixKind::Phi => {
            writeln!(f, "INDV1\tINDV2\tN_AaAa\tN_AAaa\tN1_Aa\tN2_Aa\tRELATEDNESS_PHI")
                .map_err(reed::Error::from)?;
            let aux = m.aux.as_ref();
            // vcftools emits the full ordered matrix, both directions.
            for i in 0..n {
                for j in 0..n {
                    let hethet = aux.map_or(0, |a| a.n_hethet[i * n + j]);
                    let opp = aux.map_or(0, |a| a.n_opposite_hom[i * n + j]);
                    // N_Aa[i] is the marginal het count; recover it from the
                    // diagonal of n_hethet is wrong, so recompute from phi's
                    // denominator is not possible here. Write the pair counts
                    // that vcftools reports and the phi value.
                    let phi_v = m.values[i * n + j];
                    let n1 = aux.map_or(0, |a| a.n_het.get(i).copied().unwrap_or(0));
                    let n2 = aux.map_or(0, |a| a.n_het.get(j).copied().unwrap_or(0));
                    writeln!(
                        f,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{:.9}",
                        m.ids[i], m.ids[j], hethet, opp, n1, n2, phi_v
                    )
                    .map_err(reed::Error::from)?;
                }
            }
        }
        MatrixKind::Ajk => {
            writeln!(f, "INDV1\tINDV2\tRELATEDNESS_AJK").map_err(reed::Error::from)?;
            for i in 0..n {
                for j in i..n {
                    writeln!(f, "{}\t{}\t{:.9}", m.ids[i], m.ids[j], m.values[i * n + j])
                        .map_err(reed::Error::from)?;
                }
            }
        }
        MatrixKind::Discordance => {
            writeln!(
                f,
                "#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t\
                 [5]Average -log P(HWE)\t[6]Number of sites compared\t\
                 [7]Number of matching genotypes"
            )
            .map_err(reed::Error::from)?;
            let aux = m.aux.as_ref();
            for i in 0..n {
                for j in 0..i {
                    let cmp = aux.map_or(0, |a| a.n_compared[i * n + j]);
                    let mt = aux.map_or(0, |a| a.n_matching[i * n + j]);
                    writeln!(
                        f,
                        "DCv2\t{}\t{}\t{:.9e}\t0\t{}\t{}",
                        m.ids[i], m.ids[j], m.values[i * n + j], cmp, mt
                    )
                    .map_err(reed::Error::from)?;
                }
            }
        }
        MatrixKind::PiHat => {
            writeln!(
                f,
                "FID1\tIID1\tFID2\tIID2\tRT\tEZ\tZ0\tZ1\tZ2\tPI_HAT\tPHE\tDST\tPPC\tRATIO"
            )
            .map_err(reed::Error::from)?;
            let aux = m.aux.as_ref();
            // plink writes FID and IID as separate columns; reed's ids are
            // the combined "FID_IID" sample names, so split on the last '_'
            // when it is present.
            let split = |id: &str| -> (String, String) {
                match id.rsplit_once('_') {
                    Some((fid, iid)) => (fid.to_string(), iid.to_string()),
                    None => (id.to_string(), id.to_string()),
                }
            };
            for i in 0..n {
                for j in 0..i {
                    let (f1, i1) = split(&m.ids[i]);
                    let (f2, i2) = split(&m.ids[j]);
                    let z0 = aux.map_or(f64::NAN, |a| a.z0[i * n + j]);
                    let z1 = aux.map_or(f64::NAN, |a| a.z1[i * n + j]);
                    let z2 = aux.map_or(f64::NAN, |a| a.z2[i * n + j]);
                    writeln!(
                        f,
                        "{f1}\t{i1}\t{f2}\t{i2}\tOT\tNA\t{z0:.4}\t{z1:.4}\t{z2:.4}\t{:.4}\t-9\t0\t0\t0",
                        m.values[i * n + j]
                    )
                    .map_err(reed::Error::from)?;
                }
            }
        }
        _ => unreachable!("is_matrix_kind covers the rest"),
    }
    Ok(())
}
