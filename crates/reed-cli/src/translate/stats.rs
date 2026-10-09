//! `reed stats` flag translation (02-cli.md §5): statistics named after the
//! method, one request, one driver call.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use reed::{assign_groups, default_kinds, run_stats_opened, Levels, StatsRequest};
use reed::budget::{MemBudget, ThreadBudget};

use crate::cli::StatsCli;
use crate::units::input_path;

/// Prints the engine's stage metrics (doc_v4 03 A4-5): the reader, reducer
/// and frame-chain cumulative times without which "reader-bound" is only an
/// inference.
fn print_engine_metrics(m: &reed::metrics::StageMetrics) {
    eprintln!("chunks_read\t{}", m.n_chunks_read);
    eprintln!("engine.reader_ms\t{:.3}", m.reader_time.as_secs_f64() * 1e3);
    eprintln!("engine.reduce_ms\t{:.3}", m.reduce_time.as_secs_f64() * 1e3);
    eprintln!("engine.chain_wait_ms\t{:.3}", m.chain_wait.as_secs_f64() * 1e3);
    eprintln!("worker_busy_ms\t{:.3}", m.worker_busy.as_secs_f64() * 1e3);
}

/// Whether the shared site flags left the spec at its defaults (no
/// filtering requested, so the wrapper is skipped entirely).
fn site_spec_is_empty(spec: &reed::spec::FilterSpec) -> bool {
    
    spec.site == reed::spec::SiteRules::default()
        && spec.genotype == reed::spec::GenotypeRules::default()
        && spec.site_after_mask == reed::spec::TailRules::default()
}

/// The parsed region flags: indexed regions, target includes, target
/// excludes.
type ParsedRegions = (
    Vec<reed::chunk::region::Region>,
    Vec<reed::spec::region::RegionTriple>,
    Vec<reed::spec::region::RegionTriple>,
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
                let parsed = reed::spec::Region::parse(token)
                    .map_err(|e| format!("invalid region {token:?}: {e}"))?;
                regions.push(reed::chunk::region::Region {
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
                let parsed = reed::spec::Region::parse(token)
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
        fst_global,
        f2,
        f3,
        f4,
        dstat,
        f4_ratio,
        block_size,
        ld_decay,
        ld_pairs,
        indep_pairwise,
        pbs,
        pbs_normed,
        sfs,
        proj,
        sfs_int,
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
        reed::open_regions(open_path.as_deref(), &regions)
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
        Some(std::sync::Arc::new(reed::plan::compile(
            &site_spec,
            &reed::header::tag_table(&facts.header),
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
    // .hardy, .het, .smiss/.vmiss, .stats in the reference layouts. One pass
    // serves every requested family — the QC summary accumulates inside the
    // report job (O-1) — and `-t/-T` gate the records (F4-10).
    if freq || hardy || het || missing || depth || summary {
        let Some(prefix) = out.clone() else {
            return Err("report mode writes files: pass -o PREFIX".into());
        };
        let request = reed::reports::ReportRequest {
            freq,
            hardy,
            het,
            missing,
            depth,
            summary,
            targets: targets.clone(),
            targets_exclude: targets_exclude.clone(),
            prefix: PathBuf::from(&prefix),
        };
        let source = opened
            .take()
            .expect("invariant: the report mode consumes the first open");
        let reports_start = std::time::Instant::now();
        let outcome = reed::reports::run_reports_opened(source, &request, &budget)?;
        let reports_secs = reports_start.elapsed().as_secs_f64();
        let mut files = outcome.files.clone();
        if let Some(qc) = &outcome.qc {
            let mut path = std::ffi::OsString::from(&prefix);
            path.push(".stats");
            let text = reed::summary::render_bcftools_stats(qc);
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
        return Ok(());
    }

    // Windowed statistics need the populations file; a bare `--sfs` (no
    // windowed statistic) treats every sample as one population.
    let fstats_requested =
        !f2.is_empty() || !f3.is_empty() || !f4.is_empty() || !dstat.is_empty() || !f4_ratio.is_empty();
    let windowed = pi
        || dxy
        || theta
        || tajima_d
        || fst.is_some()
        || !pbs.is_empty()
        || fstats_requested
        || ld_decay.is_some()
        || !indep_pairwise.is_empty();
    let (gt_groups, group_names) = match &pop {
        Some(groups) => assign_groups(Path::new(groups), &sample_names, force_samples)?,
        None if sfs.is_some() && !windowed => {
            (vec![Some(0); sample_names.len()], vec!["all".to_string()])
        }
        None => {
            return Err(
                "a windowed statistic needs --pop FILE (sample<TAB>population, no header)".into()
            );
        }
    };
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
    if windowed || sfs.is_some() {
        let step = if step == 0 { size } else { step };
        // A bare `--sfs` has no window denominators, so the all-sites check
        // does not apply (05 §3.1 accepts variant-only input).
        let ld_requested = ld_decay.is_some() || !indep_pairwise.is_empty();
        let invariant = if variant_only
            || (sfs.is_some() && !windowed)
            || fstats_requested
            || ld_requested
        {
            reed::InvariantMode::VariantOnly
        } else if bypass_invariant_check {
            eprintln!(
                "reed stats: warning: --bypass-invariant-check computes over present sites \
                 only; π and dxy are overestimated"
            );
            reed::InvariantMode::Bypass
        } else {
            reed::InvariantMode::Require
        };
        let pixy = reed::PixyOptions { invariant };
        let fst_kind = match fst.as_deref() {
            None => None,
            Some("wc84") => Some(reed::FstKind::WeirCockerham84),
            Some("hudson") => Some(reed::FstKind::Hudson),
            Some(other) => {
                return Err(format!("unknown --fst method {other:?} (expected wc84 or hudson)").into());
            }
        };
        let pbs_kind = if pbs_normed {
            reed::PbsKind::PbsN1
        } else {
            reed::PbsKind::Yi2010
        };
        // `--fst-global` replaces the pairwise matrix with one r-population
        // WC84 value (docs/rounds/v5/05 §3.2); only WC84 has that form.
        if fst_global {
            if fst_kind == Some(reed::FstKind::Hudson) {
                return Err("hudson Fst has no r-population (global) form; use wc84".into());
            }
            if group_names.len() == 2 {
                eprintln!(
                    "reed stats: note: two populations: --fst-global equals the pairwise wc84 value"
                );
            }
        }
        let mut kinds = reed::pixy_kinds(
            group_names.len(),
            pi,
            dxy,
            theta,
            tajima_d,
            if fst_global { None } else { fst_kind },
            None,
        );
        if fst_global {
            kinds.push(reed::StatKind::FstGlobal);
        }
        for &(a, b, c) in &pbs_triples {
            kinds.push(reed::StatKind::Pbs { a, b, c, kind: pbs_kind });
        }
        // `--sfs [folded] --proj N[,N…]`: a genome-wide spectrum per
        // population (docs/rounds/v5/05-algorithms.md §3.1).
        let sfs_request = match sfs.as_deref() {
            None => None,
            Some(mode) => {
                let Some(raw) = proj.as_deref() else {
                    return Err("--sfs needs --proj N (the projection size in chromosomes)".into());
                };
                let mut values = Vec::new();
                for token in raw.split(',') {
                    let token = token.trim();
                    values.push(token.parse::<u32>().map_err(|_| {
                        format!("invalid --proj value {token:?} (want an integer)")
                    })?);
                }
                if values.len() == 1 {
                    values = vec![values[0]; group_names.len()];
                }
                if values.len() != group_names.len() {
                    return Err(format!(
                        "--proj gives {} value(s) but the population file declares {}",
                        values.len(),
                        group_names.len()
                    )
                    .into());
                }
                Some(reed::SfsRequest { projections: values, folded: mode == "folded" })
            }
        };
        // Patterson f-statistics (docs/rounds/v5/05-algorithms.md §3.3):
        // genome-wide, one row per test, block-jackknife SE.
        let fstats_request = build_fstats_request(
            &f2,
            &f3,
            &f4,
            &dstat,
            &f4_ratio,
            block_size,
            &group_names,
        )?;
        // LD: `--ld-decay` / `--indep-pairwise` (docs/rounds/v5/05 §3.4).
        let ld_request = build_ld_request(
            ld_decay.as_deref(),
            ld_pairs,
            &indep_pairwise.join(" "),
        )?;
        // The SFS is genome-wide: with no windowed kind, skip the window
        // grid entirely (genome level only) and give it a nominal size.
        let sfs_only = kinds.is_empty() && fstats_request.is_none();
        let _ = windowed;
        let size = if size == 0 && sfs_only { 1000 } else { size };
        let step = if step == 0 { size } else { step };
        let levels = if sfs_only {
            Levels { window: false, contig: false, genome: true }
        } else {
            // Honor the user's --levels (the default stays window-only, so
            // the pixy-layout files keep their historical shape).
            Levels::parse(&levels)?
        };
        let request = StatsRequest {
            kinds: kinds.clone(),
            regions: regions.clone(),
            targets: targets.clone(),
            targets_exclude: targets_exclude.clone(),
            keep_samples: keep_samples.clone(),
            filter: filter.clone(),
            sfs: sfs_request,
            fstats: fstats_request,
            ld: ld_request,
            gt_groups,
            window_size: size,
            step,
            levels,
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
        let mut files = Vec::new();
        if let Some(ld_report) = &report.ld {
            files.extend(reed::ld::write_ld(Path::new(&prefix), ld_report)?);
        }
        if let Some(fstats_report) = &report.fstats {
            files.extend(reed::fstats::write_fstats(
                Path::new(&prefix),
                &group_names,
                fstats_report,
            )?);
        }
        if let Some(sfs_report) = &report.sfs {
            files.extend(reed::sfs::write_projected_sfs(
                Path::new(&prefix),
                &group_names,
                sfs_report,
                sfs_int,
            )?);
        }
        // A bare `--sfs` has no window rows to write.
        if !kinds.is_empty() {
            files.extend(reed::pixy::write(
                Path::new(&prefix),
                &report,
                &group_names,
                &kinds,
                pixy,
                &command_line,
                env!("CARGO_PKG_VERSION"),
            )?);
        }
        if verbosity.shows_summary() {
            eprintln!(
                "{} windows, {} records in {:.1}s -> {}",
                report.n_windows,
                report.n_records,
                report.elapsed.as_secs_f64(),
                files
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
        sfs: None,
        fstats: None,
        ld: None,
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

/// Builds the f-statistics request from the `--f2/--f3/--f4/--dstat/
/// --f4-ratio` values (docs/rounds/v5/05-algorithms.md §3.3). Each value is
/// a semicolon-separated list of role groups; each role group is a
/// comma-separated population list.
fn build_fstats_request(
    f2: &[String],
    f3: &[String],
    f4: &[String],
    dstat: &[String],
    f4_ratio: &[String],
    block_size: Option<u64>,
    group_names: &[String],
) -> Result<Option<reed::FStatsRequest>, Box<dyn std::error::Error + Send + Sync>> {
    let mut tests = Vec::new();
    let index = |name: &str| -> Result<usize, Box<dyn std::error::Error + Send + Sync>> {
        group_names
            .iter()
            .position(|g| g == name)
            .ok_or_else(|| format!("population {name:?} is not in the populations file").into())
    };
    for spec in f2 {
        let parts = split_roles(spec, 2, "--f2")?;
        tests.push(reed::FStatTest::F2 { a: index(&parts[0])?, b: index(&parts[1])? });
    }
    for spec in f3 {
        let parts = split_roles(spec, 3, "--f3")?;
        tests.push(reed::FStatTest::F3 {
            c: index(&parts[0])?,
            a: index(&parts[1])?,
            b: index(&parts[2])?,
        });
    }
    for spec in f4 {
        let parts = split_roles(spec, 4, "--f4")?;
        tests.push(reed::FStatTest::F4 {
            a: index(&parts[0])?,
            b: index(&parts[1])?,
            c: index(&parts[2])?,
            d: index(&parts[3])?,
        });
    }
    for spec in dstat {
        let parts = split_roles(spec, 4, "--dstat")?;
        tests.push(reed::FStatTest::Dstat {
            a: index(&parts[0])?,
            b: index(&parts[1])?,
            c: index(&parts[2])?,
            d: index(&parts[3])?,
        });
    }
    for spec in f4_ratio {
        let parts = split_roles(spec, 3, "--f4-ratio")?;
        // A,B;C,D;X → f4(A, O; X, C) / f4(A, O; B, C) with O = the outgroup
        // written last in the spec's first role group.
        tests.push(reed::FStatTest::F4Ratio {
            o: index(&parts[1])?,
            a: index(&parts[0])?,
            x: index(&parts[2])?,
            b: index(&parts[0])?,
            c: index(&parts[1])?,
        });
    }
    if tests.is_empty() {
        return Ok(None);
    }
    let block_bp = block_size.unwrap_or(5_000_000);
    Ok(Some(reed::FStatsRequest { tests, block_bp }))
}

/// Builds the LD request from `--ld-decay [MAX]`, `--ld-pairs` and
/// `--indep-pairwise W S R2` (docs/rounds/v5/05-algorithms.md §3.4).
fn build_ld_request(
    ld_decay: Option<&str>,
    ld_pairs: bool,
    indep_pairwise: &str,
) -> Result<Option<reed::LdRequest>, Box<dyn std::error::Error + Send + Sync>> {
    let indep_pairwise = (!indep_pairwise.is_empty()).then_some(indep_pairwise);
    if ld_decay.is_none() && indep_pairwise.is_none() {
        return Ok(None);
    }
    let decay_max_bp = match ld_decay {
        None => None,
        Some("") => Some(100_000),
        Some(text) => Some(crate::units::parse_size(text)?),
    };
    let prune = match indep_pairwise {
        None => None,
        Some(text) => {
            let parts: Vec<&str> = text.split_whitespace().collect();
            if parts.len() != 3 {
                return Err(format!(
                    "--indep-pairwise expects W S R2 (e.g. `50 5 0.2`), got {text:?}"
                )
                .into());
            }
            let (window, kb) = if let Some(kb) = parts[0].strip_suffix("kb") {
                (crate::units::parse_size(kb)?, true)
            } else {
                (crate::units::parse_size(parts[0])?, false)
            };
            let step: usize = parts[1]
                .parse()
                .map_err(|_| format!("invalid --indep-pairwise step {:?}", parts[1]))?;
            let r2: f64 = parts[2]
                .parse()
                .map_err(|_| format!("invalid --indep-pairwise r2 {:?}", parts[2]))?;
            if kb && step != 1 {
                return Err("--indep-pairwise with a kb window requires step 1".into());
            }
            let window = if kb { reed::LdWindow::Bp(window) } else {
                reed::LdWindow::Variants(usize::try_from(window).unwrap_or(usize::MAX))
            };
            Some(reed::LdPruneRequest { window, step, r2 })
        }
    };
    Ok(Some(reed::LdRequest {
        decay_max_bp,
        decay_bin_bp: 1000,
        prune,
        pairs: ld_pairs,
    }))
}

/// Splits one f-statistic value into its role groups (semicolon-separated),
/// each naming one population; `--f3 C;A,B` and `--f3 C;A;B` both work.
fn split_roles(
    spec: &str,
    roles: usize,
    flag: &str,
) -> Result<Vec<String>, Box<dyn std::error::Error + Send + Sync>> {
    let groups: Vec<String> = spec
        .split(';')
        .map(|g| g.trim().to_string())
        .filter(|g| !g.is_empty())
        .collect();
    // Each role group is one population for f2/f3/f4/dstat (the comma form
    // names the populations in order), so flatten and count.
    let names: Vec<String> = groups.iter().flat_map(|g| g.split(',')).map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect();
    if names.len() != roles {
        return Err(format!(
            "{flag} expects {roles} population(s) (e.g. `{}`), got {} in {spec:?}",
            match roles {
                2 => "A,B",
                3 => "C;A,B",
                _ => "A,B;C,D",
            },
            names.len()
        )
        .into());
    }
    Ok(names)
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
