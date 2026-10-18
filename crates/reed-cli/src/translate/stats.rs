//! `reed stats` flag translation (doc_v14 02 §2.2): the per-site/per-sample
//! QC reports and the bcftools-stats summary. The windowed
//! population-genetics statistics live in `reed popg`; the flags that moved
//! are refused here with the exact rewrite (doc_v14 02 §2.4).

use std::path::PathBuf;

use reed::histogram::{HistKind, HistMetric};

use crate::cli::StatsCli;

use super::common::prelude;

/// Parses `--hist LIST` and resolves every generic tag against the header's
/// `TagTable` (cli.md §0.7: flags are checked before the first record). A tag
/// must exist and be `Integer` or `Float`; `Flag`/`String`/`Character` is a
/// usage error.
///
/// # Errors
///
/// Usage messages for an unknown metric, an unknown tag, or a non-numeric
/// tag type.
fn parse_hist(
    list: Option<&str>,
    tags: &reed::TagTable,
) -> Result<Vec<HistMetric>, reed::Error> {
    use reed::TagType;
    let Some(list) = list else {
        return Ok(Vec::new());
    };
    let mut metrics: Vec<HistMetric> = Vec::new();
    for element in list.split(',').map(str::trim).filter(|e| !e.is_empty()) {
        let mut metric = HistMetric::parse(element).map_err(reed::Error::usage)?;
        if let HistMetric::Tag { format, tag, .. } = &metric {
            let def = if *format {
                tags.format(tag)
            } else {
                tags.info(tag)
            };
            let Some(def) = def else {
                let scope = if *format { "FORMAT" } else { "INFO" };
                return Err(reed::Error::usage(format!(
                    "--hist {element:?}: the header declares no {scope}/{tag}"
                )));
            };
            let kind = match def.ty {
                TagType::Integer => HistKind::Int,
                TagType::Float => HistKind::Float,
                other => {
                    return Err(reed::Error::usage(format!(
                        "--hist {element:?}: {tag} is {other:?}; only Integer and Float tags \
                         have histograms"
                    )));
                }
            };
            metric = HistMetric::Tag {
                format: *format,
                tag: tag.clone(),
                kind,
            };
        } else {
            // The built-ins ride fields the header must declare too; without
            // this the run wrote an all-zero histogram and exited 0
            // (doc_v14 02 §3.3 promises a usage error).
            let required: Option<(&str, &str)> = match metric {
                HistMetric::Gq => Some(("FORMAT", "GQ")),
                HistMetric::Dp => Some(("FORMAT", "DP")),
                HistMetric::DpSite => Some(("INFO", "DP")),
                _ => None,
            };
            if let Some((scope, name)) = required {
                let declared = if scope == "FORMAT" {
                    tags.format(name)
                } else {
                    tags.info(name)
                };
                if declared.is_none() {
                    return Err(reed::Error::usage(format!(
                        "--hist {element:?}: the header declares no {scope}/{name}"
                    )));
                }
            }
        }
        if !metrics.contains(&metric) {
            metrics.push(metric);
        }
    }
    Ok(metrics)
}

/// The flags that moved to `reed popg`, in the order the error names them
/// (doc_v14 02 §2.4): the redirect is a usage error with the exact rewrite,
/// never a silent alias (the `--max-missing` family's pattern, cli.md §4.4).
const MOVED: &[&str] = &[
    "--pi", "--dxy", "--theta", "--tajima-d", "--fst", "--fst-global", "--pbs", "--pbs-normed",
    "--f2", "--f3", "--f4", "--dstat", "--f4-ratio", "--block-size", "--ld-decay", "--ld-pairs",
    "--sfs", "--proj", "--sfs-int", "--variant-only", "--bypass-invariant-check",
    "--contig-lengths", "--window", "--bed-file", "--step", "--levels", "--skip-empty",
];

/// The first moved flag present on the command line, if any.
fn moved_flag(args: &StatsCli) -> Option<&'static str> {
    let present: &[(&str, bool)] = &[
        ("--pi", args.pi),
        ("--dxy", args.dxy),
        ("--theta", args.theta),
        ("--tajima-d", args.tajima_d),
        ("--fst", args.fst.is_some()),
        ("--fst-global", args.fst_global),
        ("--pbs", !args.pbs.is_empty()),
        ("--pbs-normed", args.pbs_normed),
        ("--f2", !args.f2.is_empty()),
        ("--f3", !args.f3.is_empty()),
        ("--f4", !args.f4.is_empty()),
        ("--dstat", !args.dstat.is_empty()),
        ("--f4-ratio", !args.f4_ratio.is_empty()),
        ("--block-size", args.block_size.is_some()),
        ("--ld-decay", args.ld_decay.is_some()),
        ("--ld-pairs", args.ld_pairs),
        ("--sfs", args.sfs.is_some()),
        ("--proj", args.proj.is_some()),
        ("--sfs-int", args.sfs_int),
        ("--variant-only", args.variant_only),
        ("--bypass-invariant-check", args.bypass_invariant_check),
        ("--contig-lengths", args.contig_lengths.is_some()),
        ("--window", args.size.is_some()),
        ("--bed-file", args.bed_file.is_some()),
        ("--step", args.step.is_some()),
        ("--levels", args.levels.is_some()),
        ("--skip-empty", args.skip_empty),
    ];
    present
        .iter()
        .find(|(_, on)| *on)
        .map(|(name, _)| *name)
}

/// Runs the QC-report driver for the parsed flags.
///
/// # Errors
///
/// Translation (moved flags, missing `-o`, conflicting sources) and driver
/// failures.
pub(crate) fn run_stats_cli(args: StatsCli) -> Result<(), reed::Error> {
    args.verbosity.init_tracing();
    if let Some(flag) = moved_flag(&args) {
        debug_assert!(MOVED.contains(&flag));
        return Err(reed::Error::usage(format!(
            "{flag} moved to `reed popg`: the windowed population-genetics statistics live \
             there now; run `reed popg <input> --pop FILE {flag} ...`"
        )));
    }
    let _ = &args.tmp_dir;
    // `--pop` stratifies --freq/--missing/--singletons (doc_v14 02 §2.2);
    // a popfile without any of those is a silent no-op, so it is refused
    // (the D7-6 house pattern).
    if args.pop.is_some()
        && !(args.freq || args.missing || args.singletons)
    {
        return Err(reed::Error::usage(
            "--pop stratifies --freq/--missing/--singletons; pass at least one of them or drop \
             --pop",
        ));
    }
    // Both checks come before the input is opened: a usage error the
    // operator can fix must not be masked by a missing-file runtime error
    // (cli.md §0.7 validates before the first record).
    if !(args.freq
        || args.hardy
        || args.het
        || args.missing
        || args.depth
        || args.summary
        || args.singletons)
        && args.hist.is_none()
    {
        return Err(reed::Error::usage(
            "no report requested: pass --freq/--hardy/--het/--missing/--depth/--summary/--hist/\
             --singletons",
        ));
    }
    let Some(prefix) = args.out.clone() else {
        return Err(reed::Error::usage(
            "report mode writes files: pass -o PREFIX",
        ));
    };
    let p = prelude(
        "reed stats",
        &args.input,
        &args.regions_args,
        &args.sample_args,
        &args.site_args,
        args.sites_file.as_deref(),
        None,
        args.threads,
        args.memory,
    )?;
    // doc_v7 D7-6: report mode does not apply the sample subset or the
    // site/threshold filters; refuse instead of silently ignoring them.
    if p.keep_samples.is_some() || p.site.is_some() {
        return Err(reed::Error::usage(
            "report mode (--freq/--hardy/--het/--missing/--depth/--summary) does not \
             apply -s/-S/--keep/--remove or site/threshold filters yet; run \
             `reed filter` first and feed its output to reed stats",
        ));
    }
    let facts_header = p
        .opened
        .as_ref()
        .expect("invariant: the report mode consumes the first open")
        .facts()
        .header
        .clone();
    let hist = parse_hist(args.hist.as_deref(), &reed::tag_table(&facts_header))?;
    // An all-empty list (`--hist ""`, `--hist ",,"`) used to pass the guard
    // above (it tested `is_none`) and then write nothing while exiting 0.
    if args.hist.is_some() && hist.is_empty() {
        return Err(reed::Error::usage(
            "--hist: the list is empty (want qual/dp/dp-site/gq/idl/af or INFO/<TAG>/FMT/<TAG>)",
        ));
    }
    let request = reed::reports::StatsReportRequest {
        hist,
        freq: args.freq,
        hardy: args.hardy,
        het: args.het,
        missing: args.missing,
        depth: args.depth,
        summary: args.summary,
        singletons: args.singletons,
        populations: args
            .pop
            .as_deref()
            .map(|p| reed::Populations::File(std::path::PathBuf::from(p))),
        force_samples: args.force_samples,
        targets: p.targets.clone(),
        targets_exclude: p.targets_exclude.clone(),
        prefix: PathBuf::from(&prefix),
    };
    let source = p
        .opened
        .expect("invariant: the report mode consumes the first open");
    let reports_start = std::time::Instant::now();
    let outcome = reed::reports::run_stats_opened(source, &request, &p.budget)?;
    let reports_secs = reports_start.elapsed().as_secs_f64();
    let mut files = outcome.files.clone();
    if let Some(qc) = &outcome.qc {
        let mut path = std::ffi::OsString::from(&prefix);
        path.push(".stats");
        let path = std::path::PathBuf::from(path);
        let text = reed::summary::render_bcftools_stats(qc, &p.sample_names);
        reed::summary::write_bcftools_stats(&path, &text)?;
        files.push(path);
    }
    if args.verbosity.shows_summary() {
        eprintln!(
            "{} records in {:.1}s -> {}",
            outcome.n_records,
            reports_secs,
            files
                .iter()
                .map(|f| f.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    if args.verbosity.shows_timing() {
        eprintln!("# reed stats timing (reports)");
        eprintln!("phase.open_ms\t{:.3}", p.open_secs * 1e3);
        eprintln!("phase.reports_ms\t{:.3}", reports_secs * 1e3);
        eprintln!("n_records\t{}", outcome.n_records);
        super::common::print_engine_metrics(&outcome.metrics);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// doc_v14 02 §2.4: every moved flag is refused with the rewrite naming
    /// `reed popg`, and the QC flags stay accepted (guard 3/11).
    #[test]
    fn moved_flags_redirect_to_popg() {
        use clap::Parser;
        let cli = crate::cli::Cli::try_parse_from([
            "reed", "stats", "in.vcf", "--pi", "-o", "p",
        ])
        .expect("the 0.3.x spelling parses so the redirect can fire");
        let crate::cli::Command::Stats(stats) = cli.command else {
            panic!("stats subcommand");
        };
        assert_eq!(moved_flag(&stats), Some("--pi"));
        let err = run_stats_cli(stats).expect_err("moved flags must be refused");
        assert_eq!(err.kind(), reed::ErrorKind::Usage);
        assert!(err.to_string().contains("reed popg"), "{err}");
        assert!(err.to_string().contains("--pi"), "{err}");

        let cli = crate::cli::Cli::try_parse_from([
            "reed", "stats", "in.vcf", "--freq", "-o", "p",
        ])
        .expect("parse");
        let crate::cli::Command::Stats(stats) = cli.command else {
            panic!("stats subcommand");
        };
        assert_eq!(moved_flag(&stats), None, "QC flags do not redirect");
    }

    /// A bare `reed stats` (no report flag) is a usage error naming the
    /// report flags (C14-2).
    #[test]
    fn bare_stats_is_a_usage_error() {
        use clap::Parser;
        let cli = crate::cli::Cli::try_parse_from(["reed", "stats", "in.vcf", "-o", "p"])
            .expect("parse");
        let crate::cli::Command::Stats(stats) = cli.command else {
            panic!("stats subcommand");
        };
        let err = run_stats_cli(stats).expect_err("bare stats must be refused");
        assert_eq!(err.kind(), reed::ErrorKind::Usage);
        assert!(err.to_string().contains("--freq"), "{err}");
    }

    /// `--pop` without any flag it stratifies is a usage error (D7-6: no
    /// silent no-ops); with one of them the request carries the file.
    #[test]
    fn pop_requires_a_stratified_report() {
        use clap::Parser;
        let dir = std::env::temp_dir().join(format!("reed-stats-pop-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let pops = dir.join("pops.tsv");
        std::fs::write(&pops, "s0\tA\n").expect("invariant: pops");
        let parse = |argv: &[&str]| {
            let cli = crate::cli::Cli::try_parse_from(argv).expect("parse");
            let crate::cli::Command::Stats(stats) = cli.command else {
                panic!("stats subcommand");
            };
            stats
        };
        let stats = parse(&[
            "reed", "stats", "in.vcf", "--pop", pops.to_str().unwrap(), "-o", "p",
        ]);
        let err = run_stats_cli(stats).expect_err("--pop alone is a silent no-op");
        assert_eq!(err.kind(), reed::ErrorKind::Usage);
        assert!(err.to_string().contains("--freq"), "{err}");

        let stats = parse(&[
            "reed", "stats", "in.vcf", "--pop", pops.to_str().unwrap(), "--freq", "-o", "p",
        ]);
        // The input does not exist, so the run fails at open time — after
        // the --pop checks, proving the stratified route accepted it.
        let err = run_stats_cli(stats).expect_err("no such input");
        assert_ne!(err.kind(), reed::ErrorKind::Usage, "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
