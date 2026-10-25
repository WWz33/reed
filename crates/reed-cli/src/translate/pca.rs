//! `reed pca` flag translation (doc_v16 §P2): every flag maps onto
//! [`reed::PcaRequest`]; the driver owns the semantics.

use reed::PcaRequest;
use reed::budget::{MemBudget, ThreadBudget};

use crate::cli::PcaCli;

use super::common::parse_sample_keep;

/// Runs `reed pca`.
///
/// # Errors
///
/// Usage errors from the flag parsers, and whatever the driver reports.
pub(crate) fn run_pca_cli(args: PcaCli) -> Result<(), reed::Error> {
    let Some(prefix) = args.output.as_deref() else {
        return Err(reed::Error::usage(
            "-o PREFIX is required: PCA writes PREFIX.eigenvec, PREFIX.eigenval and PREFIX.pve",
        ));
    };
    if args.pcs == 0 {
        return Err(reed::Error::usage("--pcs must be at least 1"));
    }

    // Open once so stdin streams, and to resolve the sample names for the
    // subset flags before the run starts (cli.md §0.7).
    let path = crate::units::input_path(&args.input);
    let opened = reed::open_regions(path.as_deref(), &[])?;
    let sample_names = opened.facts().sample_names.clone();
    let keep = parse_sample_keep(&args.samples, &sample_names).map_err(reed::Error::usage)?;

    let threads = args.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .expect("invariant: parallelism")
            .get()
    });
    let budget = ThreadBudget {
        threads,
        mem: args
            .memory
            .map_or_else(MemBudget::default_half_memory, |m| {
                MemBudget::from_bytes(m.saturating_mul(1 << 20))
            }),
    };

    // `--dapc` turns the run into a discriminant analysis on the components.
    if args.dapc {
        let Some(pop) = args.pop else {
            return Err(reed::Error::usage("--dapc needs --pop FILE"));
        };
        let xval = if args.xval {
            let n_pca = match &args.xval_pca {
                Some(list) => Some(
                    list.split(',')
                        .map(|v| {
                            v.trim().parse::<usize>().map_err(|_| {
                                reed::Error::usage(format!(
                                    "--xval-pca {v:?} is not a component count"
                                ))
                            })
                        })
                        .collect::<Result<Vec<_>, _>>()?,
                ),
                None => None,
            };
            Some(reed::XvalRequest {
                n_rep: args.xval_rep,
                training_set: args.xval_training,
                result: args.xval_result,
                seed: args.xval_seed,
                n_pca,
            })
        } else {
            None
        };
        let request = reed::DapcRequest {
            n_pca: Some(args.pcs),
            n_da: args.n_da,
            populations: pop,
            stand: args.stand,
            missing: args.missing,
            solver: args.solver,
            prefix: prefix.to_string(),
            keep_samples: keep,
            force_samples: args.force_samples,
            threads: Some(threads),
            memory: args.memory,
            xval,
        };
        let report = reed::run_dapc_opened(opened, &request, &budget)?;
        eprintln!(
            "reed pca --dapc: {} samples, {} sites, {} groups, \
             n_pca={}, n_da={}, var={:.4} ({:.2}s)",
            report.n_samples,
            report.n_sites,
            report.n_groups,
            report.n_pca,
            report.n_da,
            report.var,
            report.elapsed.as_secs_f64()
        );
        if let Some(xv) = &report.xval {
            eprintln!(
                "reed pca --dapc --xval: best n_pca={} (mean success {:.4}), \
                 highest mean success at n_pca={}, random chance {:.4}",
                xv.best_n_pca,
                xv.candidates
                    .iter()
                    .find(|(n, _)| *n == xv.best_n_pca)
                    .map_or(f64::NAN, |&(_, s)| s),
                xv.best_mean_success_n_pca,
                xv.random_chance
            );
            // `--xval` selects n_pca itself; say so when it overrode a
            // different `--pcs`, or the flag looks broken.
            if xv.best_n_pca != args.pcs {
                eprintln!(
                    "  note: --pcs {} was overridden by the cross-validated choice ({})",
                    args.pcs, xv.best_n_pca
                );
            }
        }
        return Ok(());
    }

    let request = PcaRequest {
        pcs: args.pcs,
        stand: args.stand,
        missing: args.missing,
        solver: args.solver,
        prefix: prefix.to_string(),
        keep_samples: keep,
        threads: Some(threads),
        memory: args.memory,
    };

    let report = reed::run_pca_opened(opened, &request, &budget)?;
    eprintln!(
        "reed pca: {} samples, {} sites, {} PCs, solver={} ({:.2}s)",
        report.n_samples,
        report.n_sites,
        report.n_pcs,
        report.solver,
        report.elapsed.as_secs_f64()
    );
    Ok(())
}
