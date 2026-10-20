//! `reed pca` flag translation (doc_v16 §P2): every flag maps onto
//! [`reed::PcaRequest`]; the driver owns the semantics.

use reed::budget::{MemBudget, ThreadBudget};
use reed::PcaRequest;

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
