//! `reed stats` flag translation (02-cli.md §5): statistics named after the
//! method, one request, one driver call.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use reed::stats::{assign_groups, default_kinds, run_stats, Levels, StatsRequest};
use reed_engine::budget::{MemBudget, ThreadBudget};

use crate::cli::StatsCli;
use crate::units::input_path;

/// Runs the statistics driver for the parsed flags.
///
/// # Errors
///
/// Translation (unknown levels/sizes) and driver failures.
pub fn run_stats_cli(args: StatsCli) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let StatsCli {

        input,
        groups,
        pi,
        dxy,
        theta,
        tajima_d,
        fst,
        pbs,
        pbs_normed,
        freq,
        hardy,
        het,
        missing,
        depth,
        summary,
        variant_only,
        bypass_invariant_check,
        force_samples,
        contig_lengths_path,
        size,
        step,
        levels,
        skip_empty,
        out,
        threads,
    } = args;

    let (_source, facts) = reed::input::open(input_path(&input).as_deref())
        .map_err(|e| -> Box<dyn std::error::Error + Send + Sync> { e.to_string().into() })?;
    let sample_names: Vec<String> = facts.sample_names.clone();

    let mut contig_lengths: HashMap<Vec<u8>, u64> = facts
        .header
        .contigs()
        .iter()
        .filter_map(|(name, contig)| {
            contig
                .length()
                .and_then(|len| u64::try_from(len).ok())
                .map(|len| (name.as_bytes().to_vec(), len))
        })
        .collect();
    if let Some(path) = &contig_lengths_path {
        for line in std::fs::read_to_string(Path::new(path))?.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut cols = line.split('\t');
            let (Some(name), Some(len)) = (cols.next(), cols.next()) else {
                continue;
            };
            let len: u64 = len
                .parse()
                .map_err(|_| format!("{path}: contig length for {name:?} is not a number"))?;
            contig_lengths.insert(name.as_bytes().to_vec(), len);
        }
    }

    let threads = threads.unwrap_or_else(|| {
        std::thread::available_parallelism().expect("invariant: parallelism").get()
    });
    let budget = ThreadBudget { threads, mem: MemBudget::default_half_memory() };

    // Per-site/per-sample report mode (02-cli.md §5.4): PREFIX.afreq,
    // .hardy, .het, .smiss/.vmiss, .stats in the reference layouts.
    if freq || hardy || het || missing || depth || summary {
        let Some(prefix) = out.clone() else {
            return Err("report mode writes files: pass -o PREFIX".into());
        };
        let mut files = Vec::new();
        if freq || hardy || het || missing || depth {
            let request = reed::reports::ReportRequest {
                freq,
                hardy,
                het,
                missing,
                depth,
                prefix: PathBuf::from(&prefix),
            };
            let outcome = reed::reports::run_reports(
                input_path(&input).as_deref(),
                &request,
                &budget,
            )?;
            files.extend(outcome.files);
            eprintln!(
                "{} records in {:.1}s -> {}",
                outcome.n_records,
                // reports are fast; report wall time via the outcome count only
                0.0,
                files.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ")
            );
        }
        if summary {
            let outcome = reed::qc_report::run_qc_report(
                input_path(&input).as_deref(),
                &reed::qc_report::QcOptions { threads, json: None },
            )?;
            let mut path = std::ffi::OsString::from(&prefix);
            path.push(".stats");
            std::fs::write(std::path::PathBuf::from(&path), &outcome.bcftools_text)?;
            files.push(std::path::PathBuf::from(path));
        }
        return Ok(());
    }

    // Windowed modes need the populations file.
    let Some(groups) = groups else {
        return Err("a windowed statistic needs --pop FILE (sample<TAB>population, no header)".into());
    };
    let (gt_groups, group_names) =
        assign_groups(Path::new(&groups), &sample_names, force_samples)?;
    eprintln!(
        "{} samples, {} groups: {}",
        sample_names.len(),
        group_names.len(),
        group_names.join(", ")
    );

    // pixy-layout mode: one command writes PREFIX_pi.txt / PREFIX_dxy.txt
    // with pixy's columns (02-cli.md §5).
    if pi || dxy || theta || tajima_d || fst.is_some() || pbs {
        let step = if step == 0 { size } else { step };
        let invariant = if variant_only {
            reed::stats::InvariantMode::VariantOnly
        } else if bypass_invariant_check {
            eprintln!(
                "reed stats: warning: --bypass-invariant-check computes over present sites \
                 only; π and dxy are overestimated"
            );
            reed::stats::InvariantMode::Bypass
        } else {
            reed::stats::InvariantMode::Require
        };
        let pixy = reed::stats::PixyOptions { invariant };
        let fst_kind = match fst.as_deref() {
            None => None,
            Some("wc84") => Some(reed::stats::FstKind::WeirCockerham84),
            Some("hudson") => Some(reed::stats::FstKind::Hudson),
            Some(other) => {
                return Err(format!("unknown --fst method {other:?} (expected wc84 or hudson)").into());
            }
        };
        let pbs_kind = pbs.then_some(if pbs_normed {
            reed::stats::PbsKind::PbsN1
        } else {
            reed::stats::PbsKind::Yi2010
        });
        let kinds = reed::stats::pixy_kinds(
            group_names.len(),
            pi,
            dxy,
            theta,
            tajima_d,
            fst_kind,
            pbs_kind,
        );
        let request = StatsRequest {
            kinds: kinds.clone(),
            gt_groups,
            window_size: size,
            step,
            levels: Levels { window: true, contig: false, genome: false },
            skip_empty: false,
            contig_lengths,
            pixy: Some(pixy),
        };
        let report = run_stats(input_path(&input).as_deref(), request, &budget)?;
        let Some(prefix) = out else {
            return Err("pixy mode writes files: pass -o PREFIX".into());
        };
        let command_line = std::env::args().collect::<Vec<_>>().join(" ");
        let written = reed::pixy::write(
            Path::new(&prefix),
            &report,
            &group_names,
            &kinds,
            pixy,
            &command_line,
            env!("CARGO_PKG_VERSION"),
        )?;
        eprintln!(
            "{} windows, {} records in {:.1}s -> {}",
            report.n_windows,
            report.n_records,
            report.elapsed.as_secs_f64(),
            written
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
        return Ok(());
    }

    let request = StatsRequest {
        kinds: default_kinds(group_names.len()),
        gt_groups,
        window_size: size,
        step: if step == 0 { size } else { step },
        levels: Levels::parse(&levels)?,
        skip_empty,
        contig_lengths,
        pixy: None,
    };
    let report = run_stats(input_path(&input).as_deref(), request, &budget)?;
    eprintln!(
        "{} windows, {} records, {} multiallelic skipped, {} chunks in {:.1}s",
        report.n_windows,
        report.n_records,
        report.n_multiallelic_skipped,
        report.n_chunks,
        report.elapsed.as_secs_f64()
    );
    match out {
        Some(path) => std::fs::write(path, report.tsv)?,
        None => print!("{}", report.tsv),
    }
    Ok(())
}
