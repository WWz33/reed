//! `reed distance`: the individual-by-individual matrix family (doc_v17).
//!
//! The pipeline itself lives in [`reed::run_distance_opened`] (doc_v30
//! P1-5); this front end only translates flags, opens the input once to
//! resolve the sample names, and renders the summary.

use reed::ThreadBudget;

use crate::cli::DistanceCli;

/// Runs `reed distance`.
///
/// # Errors
///
/// Usage errors for bad flags, and I/O or engine failures.
pub(crate) fn run_distance_cli(args: DistanceCli) -> Result<(), reed::Error> {
    // doc_v21 S34: distance had no --verbosity, so scripts could not silence it.
    args.verbosity.init_tracing();
    let kind = reed::DistanceKind::parse(&args.matrix_type).map_err(reed::Error::usage)?;
    // doc_v23 V23-32: `--format` only means something for the matrix kinds.
    // Silently ignoring it on `.genome`/`.relatedness2`/… let a script
    // "configure" a layout that was never used.
    if !kind.is_matrix() && args.format.is_some() {
        return Err(reed::Error::usage(format!(
            "--format applies to matrix --type values; {} is a per-pair table",
            args.matrix_type
        )));
    }
    let format = match &args.format {
        Some(text) => reed::DistanceFormat::parse(text).map_err(reed::Error::usage)?,
        None => reed::DistanceFormat::Triangle,
    };
    let stand = reed::Stand::parse(&args.stand).map_err(reed::Error::usage)?;
    let Some(prefix) = args.output.as_deref() else {
        return Err(reed::Error::usage(
            "-o PREFIX is required: distance writes PREFIX.<suffix>",
        ));
    };
    // doc_v23 V23-21: the fourth prefix command was missing both the empty
    // check (doc_v21 S19) and the `-o -` check.
    super::common::require_output_prefix(prefix, "reed distance")?;

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

    let budget = ThreadBudget::from_cli(args.threads, args.memory);
    let abort = reed::abort_for(&budget);
    let ctx = reed::RunContext::new(&budget, &abort.signal);

    let request = reed::DistanceRequest {
        kind,
        format,
        stand,
        missing: if args.mean_impute {
            reed::MissingPolicy::MeanImpute
        } else {
            reed::MissingPolicy::Pairwise
        },
        unbounded: args.unbounded,
        nudge: args.nudge,
        mismatch_count: args.mismatch_count,
        error_probability: args.error_probability,
        tree: args
            .tree
            .as_deref()
            .map(reed::TreeMethod::parse)
            .transpose()
            .map_err(reed::Error::usage)?,
        bootstrap: args.bootstrap,
        bootstrap_seed: args.bootstrap_seed,
        bootstrap_cutoff: args.bootstrap_cutoff,
        digits: args.digits,
        prefix: prefix.to_string(),
        keep_samples: keep,
    };
    let report = reed::run_distance_opened(opened, &request, &ctx)?;

    if args.verbosity.shows_summary() {
        if let Some(method) = request.tree {
            if let Some(reps) = request.bootstrap {
                eprintln!(
                    "reed distance --tree --bootstrap: {} replicates, {} internal nodes",
                    reps,
                    report.bootstrap_support.as_ref().map_or(0, Vec::len)
                );
            }
            eprintln!(
                "reed distance --tree: {} samples, {} sites, method={method}",
                report.n_samples, report.n_sites
            );
        } else {
            eprintln!(
                "reed distance: {} samples, {} sites, type={}, format={}",
                report.n_samples,
                report.n_sites,
                report.kind.name(),
                report.format.name()
            );
        }
    }
    let _ = report.files;
    Ok(())
}
