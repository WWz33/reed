//! `reed stats` flag translation (02-cli.md §5): statistics named after the
//! method, one request, one driver call.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use reed::stats::{assign_groups, default_kinds, run_stats_opened, Levels, StatsRequest};
use reed_engine::budget::{MemBudget, ThreadBudget};

use crate::cli::StatsCli;
use crate::units::input_path;

/// Prints the engine's stage metrics (doc_v4 03 A4-5): the reader, reducer
/// and frame-chain cumulative times without which "reader-bound" is only an
/// inference.
fn print_engine_metrics(m: &reed_engine::metrics::StageMetrics) {
    eprintln!("chunks_read\t{}", m.n_chunks_read);
    eprintln!("engine.reader_ms\t{:.3}", m.reader_time.as_secs_f64() * 1e3);
    eprintln!("engine.reduce_ms\t{:.3}", m.reduce_time.as_secs_f64() * 1e3);
    eprintln!("engine.chain_wait_ms\t{:.3}", m.chain_wait.as_secs_f64() * 1e3);
    eprintln!("worker_busy_ms\t{:.3}", m.worker_busy.as_secs_f64() * 1e3);
}

/// Whether the shared site flags left the spec at its defaults (no
/// filtering requested, so the wrapper is skipped entirely).
fn site_spec_is_empty(spec: &reed_filter::spec::FilterSpec) -> bool {
    
    spec.site == reed_filter::spec::SiteRules::default()
        && spec.genotype == reed_filter::spec::GenotypeRules::default()
        && spec.site_after_mask == reed_filter::spec::TailRules::default()
}

/// The parsed region flags: indexed regions, target includes, target
/// excludes.
type ParsedRegions = (
    Vec<reed_io::chunk::region::Region>,
    Vec<reed_filter::spec::region::RegionTriple>,
    Vec<reed_filter::spec::region::RegionTriple>,
);

/// Parses the shared region flags: indexed regions (`-r/-R`, needing a
/// seekable input) and streamed targets (`-t/-T`, `^` = exclude).
///
/// # Errors
///
/// Malformed region spellings and unreadable region files.
fn parse_regions(
    args: &crate::cli::filter_flags::RegionArgs,
) -> Result<ParsedRegions, String> {
    let mut regions = Vec::new();
    for (list, from_file) in [
        (&args.regions, false),
        (&args.regions_file, true),
    ] {
        if let Some(value) = list {
            let text = if from_file {
                std::fs::read_to_string(value).map_err(|e| format!("{value}: {e}"))?
            } else {
                value.clone()
            };
            for token in text.split(',').map(str::trim).filter(|t| !t.is_empty()) {
                let parsed = reed_filter::spec::Region::parse(token)
                    .map_err(|e| format!("invalid region {token:?}: {e}"))?;
                regions.push(reed_io::chunk::region::Region {
                    contig: parsed.contig.clone(),
                    start1: parsed.start,
                    end1: parsed.end,
                });
            }
        }
    }
    let mut include = Vec::new();
    let mut exclude = Vec::new();
    for (list, from_file) in [
        (&args.targets, false),
        (&args.targets_file, true),
    ] {
        if let Some(value) = list {
            let text = if from_file {
                std::fs::read_to_string(value).map_err(|e| format!("{value}: {e}"))?
            } else {
                value.clone()
            };
            for raw in text.split(',').map(str::trim).filter(|t| !t.is_empty()) {
                let (token, include_it) = match raw.strip_prefix('^') {
                    Some(rest) => (rest, false),
                    None => (raw, true),
                };
                let parsed = reed_filter::spec::Region::parse(token)
                    .map_err(|e| format!("invalid region {token:?}: {e}"))?;
                let triple = (parsed.contig, parsed.start, parsed.end);
                if include_it {
                    include.push(triple);
                } else {
                    exclude.push(triple);
                }
            }
        }
    }
    Ok((regions, include, exclude))
}

/// Parses the shared sample flags into a keep mask, with the filter's
/// conflict rules (`-s` vs `-S/--keep`, `-s ^` vs `--remove`).
///
/// # Errors
///
/// Conflicting or unknown sample sources.
fn parse_sample_keep(
    args: &crate::cli::filter_flags::SampleArgs,
    header_samples: &[String],
) -> Result<Option<Vec<bool>>, String> {
    let index: std::collections::HashMap<&str, usize> = header_samples
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();
    let names_of = |text: &str| -> Result<Vec<String>, String> {
        Ok(text
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .map(str::to_string)
            .collect())
    };
    let mut inline_keep: Option<Vec<String>> = None;
    let mut inline_remove: Option<Vec<String>> = None;
    for list in &args.samples {
        let exclude = list.trim_start().starts_with('^');
        let names: Vec<String> = list
            .trim_start_matches('^')
            .split(',')
            .map(str::trim)
            .filter(|n| !n.is_empty())
            .map(str::to_string)
            .collect();
        if names.is_empty() {
            return Err(format!("empty sample list {list:?}"));
        }
        let slot = if exclude { &mut inline_remove } else { &mut inline_keep };
        if slot.is_some() {
            return Err(format!("repeated -s list of the same polarity ({list:?})"));
        }
        *slot = Some(names);
    }
    let file_keep = args.samples_file.as_ref().or(args.keep.as_ref());
    if inline_keep.is_some() && file_keep.is_some() {
        return Err("-s and -S/--keep are mutually exclusive: give one sample source".into());
    }
    if inline_remove.is_some() && args.remove.is_some() {
        return Err("-s ^… and --remove are mutually exclusive: give one removal source".into());
    }
    let mut keep = vec![true; header_samples.len()];
    let mut touched = false;
    if let Some(names) = inline_keep.or_else(|| {
        file_keep.and_then(|p| std::fs::read_to_string(p).ok()).and_then(|t| names_of(&t).ok())
    }) {
        keep.iter_mut().for_each(|k| *k = false);
        for name in &names {
            match index.get(name.as_str()) {
                Some(&i) => keep[i] = true,
                None => return Err(format!("sample {name:?} is not in the header")),
            }
        }
        touched = true;
    }
    let mut drop = vec![false; header_samples.len()];
    if let Some(names) = inline_remove.or_else(|| {
        args.remove
            .as_ref()
            .and_then(|p| std::fs::read_to_string(p).ok())
            .and_then(|t| names_of(&t).ok())
    }) {
        for name in &names {
            match index.get(name.as_str()) {
                Some(&i) => drop[i] = true,
                None => return Err(format!("sample {name:?} is not in the header")),
            }
        }
        touched = true;
    }
    for (k, d) in keep.iter_mut().zip(&drop) {
        if *d {
            *k = false;
        }
    }
    Ok(touched.then_some(keep))
}

/// `--pbs A,B,C`: resolves each requested triple to group indices (the
/// branch of interest first, Yi et al. 2010).
///
/// # Errors
///
/// A triple without exactly three names, or a name the group file lacks.
fn parse_pbs_triples(
    pbs: &[String],
    group_names: &[String],
) -> Result<Vec<(usize, usize, usize)>, String> {
    let mut triples = Vec::with_capacity(pbs.len());
    for triple in pbs {
        let names: Vec<&str> = triple
            .split(',')
            .map(str::trim)
            .filter(|n| !n.is_empty())
            .collect();
        if names.len() != 3 {
            return Err(format!(
                "--pbs wants exactly three population names (A,B,C), got {triple:?}"
            ));
        }
        let mut ids = [0_usize; 3];
        for (slot, name) in ids.iter_mut().zip(&names) {
            *slot = group_names
                .iter()
                .position(|g| g == name)
                .ok_or_else(|| format!("--pbs: unknown population {name:?}"))?;
        }
        triples.push((ids[0], ids[1], ids[2]));
    }
    Ok(triples)
}

/// Runs the statistics driver for the parsed flags.
///
/// # Errors
///
/// Translation (unknown levels/sizes) and driver failures.
pub(crate) fn run_stats_cli(args: StatsCli) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let StatsCli {

        input,
        pop,
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
        contig_lengths: contig_lengths_file,
        size,
        step,
        levels,
        skip_empty,
        out,
        threads,
        tmp_dir,
        memory,
        verbosity,
        regions_args,
        sample_args,
        site_args,
        // The removed / unsupported sample-level flags never reach this
        // point: their value parsers reject them at parse time.
        mind: _,
        max_missing_sample: _,
    } = args;

    // 03-architecture.md §3.3: one verbosity contract for every subcommand.
    verbosity.init_tracing();
    let started = std::time::Instant::now();
    // `-r/-R` (indexed) and `-t/-T` (streamed) share the filter's spellings
    // (04-cleanup §2.1).
    let (regions, mut targets, targets_exclude) = parse_regions(&regions_args)?;
    // Indexed ranges are coarse (whole blocks), so the same regions also act
    // as the record predicate — exactly how the filter's `-r` works.
    for region in &regions {
        targets.push((region.contig.clone(), region.start1, region.end1));
    }

    // Open once: stdin streams into the driver (03-performance.md O5). Every
    // mode reads once now — the report run accumulates the QC summary in the
    // same pass (doc_v4 O-1), so no spill is needed even on stdin.
    let _ = tmp_dir;
    let input_file = input_path(&input);
    let open_path: Option<PathBuf> = input_file.clone();
    let mut opened = Some(
        reed::input::open_regions(open_path.as_deref(), &regions)
            .map_err(|e| -> Box<dyn std::error::Error + Send + Sync> { e.to_string().into() })?,
    );
    let open_secs = started.elapsed().as_secs_f64();
    tracing::debug!(
        open_ms = open_secs * 1e3,
        input = %input,
        "reed stats opened the input"
    );
    let facts = &opened.as_ref().expect("invariant: just opened").facts;
    let sample_names: Vec<String> = facts.sample_names.clone();
    // `-s/-S/--keep/--remove`: the same conflict rules as the filter.
    let keep_samples = parse_sample_keep(&sample_args, &sample_names)?;
    // The site/threshold flags go through the filter's own compiler, so the
    // two commands cannot drift (04-cleanup §2.1).
    let site_spec = super::filter::site_spec(&site_args)?;
    let filter = if site_spec_is_empty(&site_spec) {
        None
    } else {
        // No biallelic gate here: the stats reporters decide per statistic
        // how to treat multiallelic sites (doc_v3 05-roadmap R7); gating
        // records would silently change the summary counts.
        Some(std::sync::Arc::new(reed_filter::plan::compile(
            &site_spec,
            &reed_io::header::tag_table(&facts.header),
            "",
            false,
        )?))
    };

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
    if let Some(path) = &contig_lengths_file {
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
    let budget = ThreadBudget {
        threads,
        mem: memory.map_or_else(MemBudget::default_half_memory, |m| MemBudget::from_bytes(m.saturating_mul(1 << 20))),
    };

    // Per-site/per-sample report mode (02-cli.md §5.4): PREFIX.afreq,
    // .hardy, .het, .smiss/.vmiss, .stats in the reference layouts.
    if freq || hardy || het || missing || depth || summary {
        if !targets.is_empty() || !targets_exclude.is_empty() {
            return Err(
                "stats report mode: -t/-T targets are not supported yet; use -r with an index".into(),
            );
        }
        let Some(prefix) = out.clone() else {
            return Err("report mode writes files: pass -o PREFIX".into());
        };
        let mut files = Vec::new();
        if freq || hardy || het || missing || depth {
            // O-1: `--summary` rides along in the same pass (one read for
            // both families).
            let request = reed::reports::ReportRequest {
                freq,
                hardy,
                het,
                missing,
                depth,
                summary,
                prefix: PathBuf::from(&prefix),
            };
            let source = opened
                .take()
                .expect("invariant: the first family consumes the first open");
            let reports_start = std::time::Instant::now();
            let outcome = reed::reports::run_reports_opened(source, &request, &budget)?;
            let reports_secs = reports_start.elapsed().as_secs_f64();
            files.extend(outcome.files);
            if let (Some(qc), true) = (&outcome.qc, summary) {
                let mut path = std::ffi::OsString::from(&prefix);
                path.push(".stats");
                let text = reed::qc_report::render_bcftools_stats(qc);
                std::fs::write(std::path::PathBuf::from(&path), &text)?;
                files.push(std::path::PathBuf::from(path));
            }
            if verbosity.shows_summary() {
                eprintln!(
                    "{} records in {:.1}s -> {}",
                    outcome.n_records,
                    reports_secs,
                    files.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ")
                );
            }
            if verbosity.shows_timing() {
                eprintln!("# reed stats timing (reports)");
                eprintln!("phase.open_ms\t{:.3}", open_secs * 1e3);
                eprintln!("phase.reports_ms\t{:.3}", reports_secs * 1e3);
                eprintln!("n_records\t{}", outcome.n_records);
                print_engine_metrics(&outcome.metrics);
            }
        }
        if summary && files.is_empty() {
            // Summary alone: the QC job is the whole pass.
            let source = match opened.take() {
                Some(source) => source,
                None => reed::input::open(open_path.as_deref())
                    .map_err(|e| -> Box<dyn std::error::Error + Send + Sync> { e.to_string().into() })?,
            };
            let summary_start = std::time::Instant::now();
            let outcome = reed::qc_report::run_qc_report_opened(
                source,
                &reed::qc_report::QcOptions { threads, json: None },
            )?;
            let mut path = std::ffi::OsString::from(&prefix);
            path.push(".stats");
            std::fs::write(std::path::PathBuf::from(&path), &outcome.bcftools_text)?;
            files.push(std::path::PathBuf::from(path));
            if verbosity.shows_summary() {
                eprintln!(
                    "reed stats: summary written in {:.1}s",
                    summary_start.elapsed().as_secs_f64()
                );
            }
            if verbosity.shows_timing() {
                eprintln!("# reed stats timing (summary)");
                eprintln!("phase.open_ms\t{:.3}", open_secs * 1e3);
                eprintln!("phase.summary_ms\t{:.3}", summary_start.elapsed().as_secs_f64() * 1e3);
            }
        }
        return Ok(());
    }

    // Windowed modes need the populations file.
    let Some(groups) = pop else {
        return Err("a windowed statistic needs --pop FILE (sample<TAB>population, no header)".into());
    };
    let (gt_groups, group_names) =
        assign_groups(Path::new(&groups), &sample_names, force_samples)?;
    if verbosity.shows_summary() {
        eprintln!(
            "{} samples, {} groups: {}",
            sample_names.len(),
            group_names.len(),
            group_names.join(", ")
        );
    }

    let pbs_triples = parse_pbs_triples(&pbs, &group_names)?;

    // pixy-layout mode: one command writes PREFIX_pi.txt / PREFIX_dxy.txt
    // with pixy's columns (02-cli.md §5).
    if pi || dxy || theta || tajima_d || fst.is_some() || !pbs.is_empty() {
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
        let pbs_kind = if pbs_normed {
            reed::stats::PbsKind::PbsN1
        } else {
            reed::stats::PbsKind::Yi2010
        };
        let mut kinds = reed::stats::pixy_kinds(
            group_names.len(),
            pi,
            dxy,
            theta,
            tajima_d,
            fst_kind,
            None,
        );
        for &(a, b, c) in &pbs_triples {
            kinds.push(reed::stats::StatKind::Pbs { a, b, c, kind: pbs_kind });
        }
        let request = StatsRequest {
            kinds: kinds.clone(),
            regions: regions.clone(),
            targets: targets.clone(),
            targets_exclude: targets_exclude.clone(),
            keep_samples: keep_samples.clone(),
            filter: filter.clone(),
            gt_groups,
            window_size: size,
            step,
            levels: Levels { window: true, contig: false, genome: false },
            skip_empty: false,
            contig_lengths,
            pixy: Some(pixy),
        };
        let source = opened.take().expect("invariant: the windowed mode reads once");
        let run_start = std::time::Instant::now();
        let report = run_stats_opened(source, request, &budget)?;
        let run_secs = run_start.elapsed().as_secs_f64();
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
        if verbosity.shows_summary() {
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
        }
        if verbosity.shows_timing() {
            eprintln!("# reed stats timing (windowed)");
            eprintln!("phase.open_ms\t{:.3}", open_secs * 1e3);
            eprintln!("phase.run_ms\t{:.3}", run_secs * 1e3);
            eprintln!("n_records\t{}", report.n_records);
            eprintln!("n_windows\t{}", report.n_windows);
            eprintln!("chunks\t{}", report.n_chunks);
            print_engine_metrics(&report.metrics);
        }
        return Ok(());
    }

    let request = StatsRequest {
        kinds: default_kinds(group_names.len()),
        regions,
        targets,
        targets_exclude,
        keep_samples,
        filter,
        gt_groups,
        window_size: size,
        step: if step == 0 { size } else { step },
        levels: Levels::parse(&levels)?,
        skip_empty,
        contig_lengths,
        pixy: None,
    };
    let source = opened.take().expect("invariant: the windowed mode reads once");
    let run_start = std::time::Instant::now();
    let report = run_stats_opened(source, request, &budget)?;
    let run_secs = run_start.elapsed().as_secs_f64();
    if verbosity.shows_summary() {
        eprintln!(
            "{} windows, {} records, {} multiallelic skipped, {} chunks in {:.1}s",
            report.n_windows,
            report.n_records,
            report.n_multiallelic_skipped,
            report.n_chunks,
            report.elapsed.as_secs_f64()
        );
    }
    if verbosity.shows_timing() {
        eprintln!("# reed stats timing (windowed)");
        eprintln!("phase.open_ms\t{:.3}", open_secs * 1e3);
        eprintln!("phase.run_ms\t{:.3}", run_secs * 1e3);
        eprintln!("n_records\t{}", report.n_records);
        eprintln!("n_windows\t{}", report.n_windows);
        eprintln!("chunks\t{}", report.n_chunks);
        print_engine_metrics(&report.metrics);
    }
    match out {
        Some(path) => std::fs::write(path, report.tsv)?,
        None => print!("{}", report.tsv),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::parse_pbs_triples;

    /// `--pbs A,B,C` resolves names in order and rejects malformed triples.
    #[test]
    fn pbs_triples_resolve_names() {
        let groups = ["pop1".to_string(), "pop2".to_string(), "pop3".to_string()];
        assert_eq!(
            parse_pbs_triples(&["pop3,pop2,pop1".to_string()], &groups).expect("valid"),
            vec![(2, 1, 0)]
        );
        assert!(parse_pbs_triples(&["pop1,pop2".to_string()], &groups).is_err());
        assert!(parse_pbs_triples(&["pop1,pop2,nope".to_string()], &groups).is_err());
    }
}
