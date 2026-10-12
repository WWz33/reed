//! `reed stats` flag translation (02-cli.md §5): statistics named after the
//! method, one request, one driver call.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use reed::budget::{MemBudget, ThreadBudget};
use reed::{Levels, Populations, StatsRequest, default_kinds, run_stats_opened};

use crate::cli::StatsCli;
use crate::units::input_path;

/// Prints the engine's stage metrics (doc_v4 03 A4-5): the reader, reducer
/// and frame-chain cumulative times without which "reader-bound" is only an
/// inference.
fn print_engine_metrics(m: &reed::metrics::StageMetrics) {
    eprintln!("chunks_read\t{}", m.n_chunks_read);
    eprintln!("engine.reader_ms\t{:.3}", m.reader_time.as_secs_f64() * 1e3);
    eprintln!("engine.reduce_ms\t{:.3}", m.reduce_time.as_secs_f64() * 1e3);
    eprintln!(
        "engine.chain_wait_ms\t{:.3}",
        m.chain_wait.as_secs_f64() * 1e3
    );
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
fn parse_regions(args: &crate::cli::filter_flags::RegionArgs) -> Result<ParsedRegions, String> {
    let mut regions = Vec::new();
    for (list, from_file) in [(&args.regions, false), (&args.regions_file, true)] {
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
    for (list, from_file) in [(&args.targets, false), (&args.targets_file, true)] {
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
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
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
        let slot = if exclude {
            &mut inline_remove
        } else {
            &mut inline_keep
        };
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
    // A `-S`/`--keep`/`--remove` file that cannot be read is an error, not a
    // silent "no filter" (doc_v7 D7-12; the old `.ok()` chain swallowed it).
    let file_names = match file_keep {
        Some(p) => {
            let text =
                std::fs::read_to_string(p).map_err(|e| format!("cannot read sample file {p}: {e}"))?;
            Some(names_of(&text)?)
        }
        None => None,
    };
    if let Some(names) = inline_keep.or(file_names) {
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
    let remove_names = match &args.remove {
        Some(p) => {
            let text =
                std::fs::read_to_string(p).map_err(|e| format!("cannot read sample file {p}: {e}"))?;
            Some(names_of(&text)?)
        }
        None => None,
    };
    if let Some(names) = inline_remove.or(remove_names) {
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
pub(crate) fn run_stats_cli(args: StatsCli) -> Result<(), reed::Error> {
    // 03-architecture.md §3.3: one verbosity contract for every subcommand.
    args.verbosity.init_tracing();
    let started = std::time::Instant::now();
    // `-r/-R` (indexed) and `-t/-T` (streamed) share the filter's spellings
    // (04-cleanup §2.1).
    let (regions, mut targets, targets_exclude) =
        parse_regions(&args.regions_args).map_err(reed::Error::usage)?;
    // Indexed ranges are coarse (whole blocks), so the same regions also act
    // as the record predicate — exactly how the filter's `-r` works.
    for region in &regions {
        targets.push((region.contig.clone(), region.start1, region.end1));
    }

    // Open once: stdin streams into the driver (03-performance.md O5). Every
    // mode reads once now — the report run accumulates the QC summary in the
    // same pass (doc_v4 O-1), so no spill is needed even on stdin.
    let _ = &args.tmp_dir;
    let input_file = input_path(&args.input);
    let open_path: Option<PathBuf> = input_file.clone();
    let mut opened = Some(reed::open_regions(open_path.as_deref(), &regions)?);
    let open_secs = started.elapsed().as_secs_f64();
    tracing::debug!(
        open_ms = open_secs * 1e3,
        input = %args.input,
        "reed stats opened the input"
    );
    let facts = &opened.as_ref().expect("invariant: just opened").facts;
    let sample_names: Vec<String> = facts.sample_names.clone();
    // `-s/-S/--keep/--remove`: the same conflict rules as the filter.
    let keep_samples =
        parse_sample_keep(&args.sample_args, &sample_names).map_err(reed::Error::usage)?;
    // The site/threshold flags go through the filter's own compiler, so the
    // two commands cannot drift (04-cleanup §2.1).
    if !args.site_args.indep_pairwise.is_empty() {
        return Err(reed::Error::usage(
            "reed stats does not support --indep-pairwise yet (the pruning needs a second \
             pass); run `reed filter --indep-pairwise W[kb] [S] R2` first and feed its \
             output to reed stats",
        ));
    }
    let site_spec = super::filter::site_spec(&args.site_args).map_err(reed::Error::usage)?;
    // No biallelic gate here: the stats reporters decide per statistic how
    // to treat multiallelic sites (doc_v3 05-roadmap R7); gating records
    // would silently change the summary counts. Compiling the typed spec
    // against the header happens inside `run_stats_opened` (08 §6.2).
    let site = (!site_spec_is_empty(&site_spec)).then_some(site_spec);

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
    if let Some(path) = &args.contig_lengths {
        for line in std::fs::read_to_string(Path::new(path))?.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut cols = line.split('\t');
            let (Some(name), Some(len)) = (cols.next(), cols.next()) else {
                continue;
            };
            let len: u64 = len.parse().map_err(|_| {
                reed::Error::usage(format!(
                    "{path}: contig length for {name:?} is not a number"
                ))
            })?;
            contig_lengths.insert(name.as_bytes().to_vec(), len);
        }
    }

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

    // `--sites-file` narrows the counted sites for every route; parse it
    // before the report branch so report mode sees it too (doc_v7 D7-6;
    // it was parsed after the report branch's early return and so never
    // reached it).
    if let Some(path) = &args.sites_file {
        for (contig, pos1) in parse_sites(path).map_err(reed::Error::usage)? {
            targets.push((contig, Some(pos1), Some(pos1)));
        }
    }

    // Per-site/per-sample report mode (02-cli.md §5.4): PREFIX.afreq,
    // .hardy, .het, .smiss/.vmiss, .stats in the reference layouts. One pass
    // serves every requested family — the QC summary accumulates inside the
    // report job (O-1) — and `-t/-T` gate the records (F4-10).
    if args.freq || args.hardy || args.het || args.missing || args.depth || args.summary {
        let Some(prefix) = args.out.clone() else {
            return Err(reed::Error::usage(
                "report mode writes files: pass -o PREFIX",
            ));
        };
        // doc_v7 D7-6: report mode does not apply the sample subset or the
        // site/threshold filters; refuse instead of silently ignoring them.
        if keep_samples.is_some() || site.is_some() {
            return Err(reed::Error::usage(
                "report mode (--freq/--hardy/--het/--missing/--depth/--summary) does not \
                 apply -s/-S/--keep/--remove or site/threshold filters yet; run \
                 `reed filter` first and feed its output to reed stats",
            ));
        }
        let request = reed::reports::ReportRequest {
            freq: args.freq,
            hardy: args.hardy,
            het: args.het,
            missing: args.missing,
            depth: args.depth,
            summary: args.summary,
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
            let path = std::path::PathBuf::from(path);
            let text = reed::summary::render_bcftools_stats(qc, &sample_names);
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
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
        if args.verbosity.shows_timing() {
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
    // A bare `--sfs`, `--ld-decay` or `--ld-pairs` treats every sample as
    // one population (only the window rows need groups).
    let populations = match &args.pop {
        Some(groups) => Populations::File(PathBuf::from(groups)),
        None => Populations::All,
    };
    // The names for kinds/fstats/sfs construction; the assignment against
    // the header happens inside `run_stats_opened` (08 §6.2).
    let group_names = populations.group_names()?;
    if args.verbosity.shows_summary() {
        eprintln!(
            "{} samples, {} groups: {}",
            sample_names.len(),
            group_names.len(),
            group_names.join(", ")
        );
    }

    // `--bed-file`: the windows are the BED intervals, so the Bp grid is
    // replaced wholesale; `--sites-file` narrows the counted sites.
    let grid = match &args.bed_file {
        Some(path) => reed::WindowGrid::Regions(std::sync::Arc::new(
            parse_bed(path).map_err(reed::Error::usage)?,
        )),
        None => reed::WindowGrid::Bp {
            size: args.size,
            step: if args.step == 0 { args.size } else { args.step },
        },
    };

    let windowed_route = args.pi
        || args.dxy
        || args.theta
        || args.tajima_d
        || args.fst.is_some()
        || !args.pbs.is_empty()
        || !(args.f2.is_empty()
            && args.f3.is_empty()
            && args.f4.is_empty()
            && args.dstat.is_empty()
            && args.f4_ratio.is_empty());
    // Window-layout mode: one command writes PREFIX_pi.txt / PREFIX_dxy.txt
    // with the reference column layout (02-cli.md §5). LD is genome-wide (no
    // window rows), but it rides the same job as the windowed statistics; a
    // bare `--ld-decay`/`--ld-pairs` takes this path with no window kinds and
    // no `--pop` (F6-9).
    if windowed_route || args.sfs.is_some() || args.ld_decay.is_some() || args.ld_pairs {
        let ctx = WindowedContext {
            regions: regions.clone(),
            targets,
            targets_exclude,
            keep_samples,
            site,
            populations: populations.clone(),
            force_samples: args.force_samples,
            contig_lengths,
            grid,
        };
        let (request, group_names, kinds, layout) = build_windowed_request(&args, ctx)?;
        let source = opened
            .take()
            .expect("invariant: the windowed mode reads once");
        let run_start = std::time::Instant::now();
        let report = run_stats_opened(source, request, &budget)?;
        let run_secs = run_start.elapsed().as_secs_f64();
        let Some(prefix) = args.out.clone() else {
            return Err(reed::Error::usage(
                "window-layout mode writes files: pass -o PREFIX",
            ));
        };
        let command_line = std::env::args().collect::<Vec<_>>().join(" ");
        let mut files = Vec::new();
        if let Some(ld_report) = &report.ld {
            files.extend(reed::ld::write_ld(
                Path::new(&prefix),
                ld_report,
                &group_names,
            )?);
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
                args.sfs_int,
            )?);
        }
        // A bare `--sfs` has no window rows to write.
        if !kinds.is_empty() {
            files.extend(reed::windows::write(
                Path::new(&prefix),
                &report,
                &group_names,
                &kinds,
                layout,
                &command_line,
                env!("CARGO_PKG_VERSION"),
            )?);
        }
        if args.verbosity.shows_summary() {
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
        if args.verbosity.shows_timing() {
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

    let kinds = default_kinds(group_names.len());
    let request = StatsRequest {
        kinds: kinds.clone(),
        site,
        populations,
        force_samples: args.force_samples,
        regions,
        targets,
        targets_exclude,
        keep_samples,
        sfs: None,
        fstats: None,
        ld: None,
        grid,
        levels: Levels::parse(&args.levels).map_err(|e| reed::Error::usage(e.to_string()))?,
        skip_empty: args.skip_empty,
        contig_lengths,
        layout: None,
    };
    let source = opened
        .take()
        .expect("invariant: the windowed mode reads once");
    let run_start = std::time::Instant::now();
    let report = run_stats_opened(source, request, &budget)?;
    let run_secs = run_start.elapsed().as_secs_f64();
    if args.verbosity.shows_summary() {
        eprintln!(
            "{} windows, {} records, {} multiallelic skipped, {} chunks in {:.1}s",
            report.n_windows,
            report.n_records,
            report.n_multiallelic_skipped,
            report.n_chunks,
            report.elapsed.as_secs_f64()
        );
    }
    if args.verbosity.shows_timing() {
        eprintln!("# reed stats timing (windowed)");
        eprintln!("phase.open_ms\t{:.3}", open_secs * 1e3);
        eprintln!("phase.run_ms\t{:.3}", run_secs * 1e3);
        eprintln!("n_records\t{}", report.n_records);
        eprintln!("n_windows\t{}", report.n_windows);
        eprintln!("chunks\t{}", report.n_chunks);
        print_engine_metrics(&report.metrics);
    }
    let text = reed::windows::render_tsv(&report.windows, &kinds, false);
    match &args.out {
        Some(path) => std::fs::write(path, text)?,
        None => print!("{text}"),
    }
    Ok(())
}

/// The shared, already-parsed pieces a windowed request needs (the caller
/// resolved them from the region/sample/site flags before this point).
pub(crate) struct WindowedContext {
    pub regions: Vec<reed::chunk::region::Region>,
    pub targets: Vec<reed::spec::region::RegionTriple>,
    pub targets_exclude: Vec<reed::spec::region::RegionTriple>,
    pub keep_samples: Option<Vec<bool>>,
    pub site: Option<reed::spec::FilterSpec>,
    pub populations: reed::Populations,
    pub force_samples: bool,
    pub contig_lengths: HashMap<Vec<u8>, u64>,
    pub grid: reed::WindowGrid,
}

/// Builds the windowed statistics request from the parsed flags; the caller
/// routed here because a windowed statistic, `--sfs`, or LD was requested.
/// Returns the request, the population names (writers render group columns
/// from these), the statistic kinds, and the layout options.
///
/// # Errors
///
/// Usage failures: a windowed statistic without `--pop`, an unknown `--fst`
/// method, `--fst-global` with hudson, a bad `--proj`, or f-statistic names
/// not present in the population file.
pub(crate) fn build_windowed_request(
    args: &StatsCli,
    ctx: WindowedContext,
) -> Result<
    (
        StatsRequest,
        Vec<String>,
        Vec<reed::StatKind>,
        reed::WindowLayoutOptions,
    ),
    reed::Error,
> {
    let fstats_requested = !(args.f2.is_empty()
        && args.f3.is_empty()
        && args.f4.is_empty()
        && args.dstat.is_empty()
        && args.f4_ratio.is_empty());
    let windowed = args.pi
        || args.dxy
        || args.theta
        || args.tajima_d
        || args.fst.is_some()
        || !args.pbs.is_empty()
        || fstats_requested;
    if ctx.populations == reed::Populations::All && windowed {
        return Err(reed::Error::usage(
            "a windowed statistic needs --pop FILE (sample<TAB>population, no header)",
        ));
    }
    let group_names = ctx.populations.group_names()?;
    let pbs_triples = parse_pbs_triples(&args.pbs, &group_names).map_err(reed::Error::usage)?;
    // A bare `--sfs` has no window denominators, so the all-sites check
    // does not apply (05 §3.1 accepts variant-only input).
    let ld_requested = args.ld_decay.is_some();
    let invariant = if args.variant_only
        || (args.sfs.is_some() && !windowed)
        || fstats_requested
        || ld_requested
    {
        reed::InvariantMode::VariantOnly
    } else if args.bypass_invariant_check {
        eprintln!(
            "reed stats: warning: --bypass-invariant-check computes over present sites \
                 only; π and dxy are overestimated"
        );
        reed::InvariantMode::Bypass
    } else {
        reed::InvariantMode::Require
    };
    let layout = reed::WindowLayoutOptions { invariant };
    let fst_kind = match args.fst.as_deref() {
        None => None,
        Some("wc84") => Some(reed::FstKind::WeirCockerham84),
        Some("hudson") => Some(reed::FstKind::Hudson),
        Some(other) => {
            return Err(reed::Error::usage(format!(
                "unknown --fst method {other:?} (expected wc84 or hudson)"
            )));
        }
    };
    let pbs_kind = if args.pbs_normed {
        reed::PbsKind::PbsN1
    } else {
        reed::PbsKind::Yi2010
    };
    // `--fst-global` replaces the pairwise matrix with one r-population
    // WC84 value (docs/rounds/v5/05 §3.2); only WC84 has that form.
    if args.fst_global {
        if fst_kind == Some(reed::FstKind::Hudson) {
            return Err(reed::Error::usage(
                "hudson Fst has no r-population (global) form; use wc84",
            ));
        }
        if group_names.len() == 2 {
            eprintln!(
                "reed stats: note: two populations: --fst-global equals the pairwise wc84 value"
            );
        }
    }
    let mut kinds = reed::window_kinds(
        group_names.len(),
        args.pi,
        args.dxy,
        args.theta,
        args.tajima_d,
        if args.fst_global { None } else { fst_kind },
        None,
    );
    if args.fst_global {
        kinds.push(reed::StatKind::FstGlobal);
    }
    for &(a, b, c) in &pbs_triples {
        kinds.push(reed::StatKind::Pbs {
            a,
            b,
            c,
            kind: pbs_kind,
        });
    }
    // `--sfs [folded] --proj N[,N…]`: a genome-wide spectrum per
    // population (docs/rounds/v5/05-algorithms.md §3.1).
    let sfs_request = match args.sfs.as_deref() {
        None => None,
        Some(mode) => {
            let Some(raw) = args.proj.as_deref() else {
                return Err(reed::Error::usage(
                    "--sfs needs --proj N (the projection size in chromosomes)",
                ));
            };
            let mut values = Vec::new();
            for token in raw.split(',') {
                let token = token.trim();
                values.push(token.parse::<u32>().map_err(|_| {
                    reed::Error::usage(format!("invalid --proj value {token:?} (want an integer)"))
                })?);
            }
            if values.len() == 1 {
                values = vec![values[0]; group_names.len()];
            }
            if values.len() != group_names.len() {
                return Err(reed::Error::usage(format!(
                    "--proj gives {} value(s) but the population file declares {}",
                    values.len(),
                    group_names.len()
                )));
            }
            Some(reed::SfsRequest {
                projections: values,
                folded: mode == "folded",
            })
        }
    };
    // Patterson f-statistics (docs/rounds/v5/05-algorithms.md §3.3):
    // genome-wide, one row per test, block-jackknife SE.
    let fstats_request = build_fstats_request(
        &args.f2,
        &args.f3,
        &args.f4,
        &args.dstat,
        &args.f4_ratio,
        args.block_size,
        &group_names,
    )?;
    // LD: `--ld-decay` / `--indep-pairwise` (docs/rounds/v5/05 §3.4).
    let ld_request = build_ld_request(args.ld_decay.as_deref(), args.ld_pairs)?;
    // The SFS is genome-wide: with no windowed kind, skip the window
    // grid entirely (genome level only) and give it a nominal size.
    let sfs_only = kinds.is_empty() && fstats_request.is_none();
    let _ = windowed;
    let levels = if sfs_only {
        Levels {
            window: false,
            contig: false,
            genome: true,
        }
    } else {
        // Honor the user's --levels (the default stays window-only, so
        // the window-layout files keep their historical shape).
        Levels::parse(&args.levels).map_err(|e| reed::Error::usage(e.to_string()))?
    };
    let request = StatsRequest {
        kinds: kinds.clone(),
        site: ctx.site,
        populations: ctx.populations,
        force_samples: ctx.force_samples,
        regions: ctx.regions,
        targets: ctx.targets,
        targets_exclude: ctx.targets_exclude,
        keep_samples: ctx.keep_samples,
        sfs: sfs_request,
        fstats: fstats_request,
        ld: ld_request,
        grid: ctx.grid.clone(),
        levels,
        skip_empty: args.skip_empty,
        contig_lengths: ctx.contig_lengths,
        layout: Some(layout),
    };
    Ok((request, group_names, kinds, layout))
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
) -> Result<Option<reed::FStatsRequest>, reed::Error> {
    let mut tests = Vec::new();
    let index = |name: &str| -> Result<usize, reed::Error> {
        group_names.iter().position(|g| g == name).ok_or_else(|| {
            reed::Error::usage(format!(
                "population {name:?} is not in the populations file"
            ))
        })
    };
    for spec in f2 {
        let parts = split_roles(spec, 2, "--f2")?;
        tests.push(reed::FStatTest::F2 {
            a: index(&parts[0])?,
            b: index(&parts[1])?,
        });
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
        let parts = split_roles(spec, 5, "--f4-ratio")?;
        // A,B;C,D;X → f4(A, O; X, C) / f4(A, O; B, C): A = parts[0],
        // O = parts[1] (the outgroup, last of the first role group),
        // C = parts[2], B = parts[3], X = parts[4] (doc_v7 D7-8; the CLI
        // used to take only three names and hardcode b = a, c = o).
        tests.push(reed::FStatTest::F4Ratio {
            a: index(&parts[0])?,
            o: index(&parts[1])?,
            c: index(&parts[2])?,
            b: index(&parts[3])?,
            x: index(&parts[4])?,
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
) -> Result<Option<reed::LdRequest>, reed::Error> {
    if ld_decay.is_none() {
        return Ok(None);
    }
    let decay_max_bp = match ld_decay {
        None => None,
        Some("") => Some(100_000),
        Some(text) => Some(crate::units::parse_size(text).map_err(reed::Error::usage)?),
    };
    Ok(Some(reed::LdRequest {
        decay_max_bp,
        decay_bin_bp: 1000,
        pairs: ld_pairs,
    }))
}

/// Splits one f-statistic value into its role groups (semicolon-separated),
/// each naming one population; `--f3 C;A,B` and `--f3 C;A;B` both work.
fn split_roles(spec: &str, roles: usize, flag: &str) -> Result<Vec<String>, reed::Error> {
    let groups: Vec<String> = spec
        .split(';')
        .map(|g| g.trim().to_string())
        .filter(|g| !g.is_empty())
        .collect();
    // Each role group is one population for f2/f3/f4/dstat (the comma form
    // names the populations in order), so flatten and count.
    let names: Vec<String> = groups
        .iter()
        .flat_map(|g| g.split(','))
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if names.len() != roles {
        return Err(reed::Error::usage(format!(
            "{flag} expects {roles} population(s) (e.g. `{}`), got {} in {spec:?}",
            match roles {
                2 => "A,B",
                3 => "C;A,B",
                5 => "A,B;C,D;X",
                _ => "A,B;C,D",
            },
            names.len()
        )));
    }
    Ok(names)
}

/// Parses a BED3 file (`chrom<TAB>start0<TAB>end0`, 0-based half-open) into
/// the per-contig window regions of a `--bed-file` grid. `track`/`browser`/
/// `#` lines are skipped; a fourth column (name) is accepted and ignored.
///
/// # Errors
///
/// A string message for unreadable files, rows with fewer than three
/// fields, non-integer or reversed coordinates (exit 2 upstream).
fn parse_bed(path: &str) -> Result<reed::ContigRegions, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|err| format!("cannot read --bed-file {path}: {err}"))?;
    let mut entries: Vec<(Vec<u8>, Vec<reed::Interval0>)> = Vec::new();
    for (line_no, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let mut cols = line.split('\t');
        let (Some(chrom), Some(start), Some(end)) = (cols.next(), cols.next(), cols.next()) else {
            return Err(format!(
                "{path}:{}: expected BED3 (chrom, start, end), got {line:?}",
                line_no + 1
            ));
        };
        let start: u64 = start
            .trim()
            .parse()
            .map_err(|_| format!("{path}:{}: start {start:?} is not an integer", line_no + 1))?;
        let end: u64 = end
            .trim()
            .parse()
            .map_err(|_| format!("{path}:{}: end {end:?} is not an integer", line_no + 1))?;
        if end <= start {
            return Err(format!(
                "{path}:{}: end ({end}) must be greater than start ({start})",
                line_no + 1
            ));
        }
        let contig = chrom.trim().as_bytes().to_vec();
        let interval = reed::Interval0::new(reed::Pos0(start), reed::Pos0(end));
        match entries.iter_mut().find(|(name, _)| *name == contig) {
            Some((_, intervals)) => intervals.push(interval),
            None => entries.push((contig, vec![interval])),
        }
    }
    reed::ContigRegions::new(entries).map_err(|err| format!("{path}: {err}"))
}

/// Parses a sites file (`CHROM<TAB>POS`, 1-based) into per-site targets.
///
/// # Errors
///
/// A string message for unreadable files or malformed rows.
fn parse_sites(path: &str) -> Result<Vec<(String, u64)>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|err| format!("cannot read --sites-file {path}: {err}"))?;
    let mut sites = Vec::new();
    for (line_no, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut cols = line.split('\t');
        let (Some(chrom), Some(pos)) = (cols.next(), cols.next()) else {
            return Err(format!(
                "{path}:{}: expected `chrom<TAB>pos`, got {line:?}",
                line_no + 1
            ));
        };
        let pos: u64 = pos
            .trim()
            .parse()
            .map_err(|_| format!("{path}:{}: POS {pos:?} is not an integer", line_no + 1))?;
        if pos == 0 {
            return Err(format!("{path}:{}: POS is 0 (1-based)", line_no + 1));
        }
        sites.push((chrom.trim().to_string(), pos));
    }
    Ok(sites)
}

#[cfg(test)]
mod tests {
    use super::{WindowedContext, build_windowed_request, parse_pbs_triples};

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

    /// A representative windowed request renders to flags clap accepts, and
    /// rebuilding from the parsed flags reproduces the request field by
    /// field (doc_v5 08 §6.2 text contract, the stats half).
    #[test]
    fn stats_request_round_trips_through_the_flags() {
        use clap::Parser;
        let dir = std::env::temp_dir().join(format!("reed-stats-rt-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let pops = dir.join("pops.tsv");
        std::fs::write(&pops, "s0\tA\ns1\tB\ns2\tC\n").expect("invariant: pops");
        let names = || vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let mut kinds = reed::window_kinds(
            3,
            true,
            true,
            true,
            true,
            Some(reed::FstKind::WeirCockerham84),
            None,
        );
        kinds.push(reed::StatKind::Pbs {
            a: 0,
            b: 1,
            c: 2,
            kind: reed::PbsKind::Yi2010,
        });
        let request = reed::StatsRequest {
            kinds: kinds.clone(),
            site: None,
            populations: reed::Populations::File(pops.clone()),
            force_samples: false,
            grid: reed::WindowGrid::Bp {
                size: 100_000,
                step: 50_000,
            },
            levels: reed::Levels {
                window: true,
                contig: true,
                genome: false,
            },
            skip_empty: false,
            contig_lengths: std::collections::HashMap::new(),
            layout: Some(reed::WindowLayoutOptions {
                // fstats/ld force VariantOnly in the builder; the request
                // must carry the CLI-shaped layout for the round trip.
                invariant: reed::InvariantMode::VariantOnly,
            }),
            regions: Vec::new(),
            targets: Vec::new(),
            targets_exclude: Vec::new(),
            keep_samples: None,
            sfs: Some(reed::SfsRequest {
                projections: vec![8, 8, 8],
                folded: false,
            }),
            fstats: Some(reed::FStatsRequest {
                tests: vec![reed::FStatTest::F3 { c: 0, a: 1, b: 2 }],
                block_bp: 5_000_000,
            }),
            ld: Some(reed::LdRequest {
                decay_max_bp: Some(100_000),
                decay_bin_bp: 1000,
                pairs: true,
            }),
        };
        let mut argv: Vec<String> = vec!["reed".into(), "stats".into(), "in.vcf".into()];
        argv.extend(request.to_args(&names()));
        let cli = crate::cli::Cli::try_parse_from(argv).expect("to_args must parse");
        let crate::cli::Command::Stats(stats) = cli.command else {
            panic!("stats subcommand");
        };
        let ctx = WindowedContext {
            regions: Vec::new(),
            targets: Vec::new(),
            targets_exclude: Vec::new(),
            keep_samples: None,
            site: None,
            populations: reed::Populations::File(pops),
            force_samples: false,
            contig_lengths: std::collections::HashMap::new(),
            grid: reed::WindowGrid::Bp {
                size: 100_000,
                step: 50_000,
            },
        };
        let (rebuilt, rebuilt_names, rebuilt_kinds, layout) =
            build_windowed_request(&stats, ctx).expect("rebuild");
        assert_eq!(rebuilt.kinds, kinds);
        assert_eq!(rebuilt_kinds, kinds);
        assert_eq!(rebuilt_names, names());
        assert_eq!(rebuilt.sfs, request.sfs);
        assert_eq!(rebuilt.fstats, request.fstats);
        assert_eq!(rebuilt.ld, request.ld);
        assert_eq!(rebuilt.levels, request.levels);
        assert!(!rebuilt.skip_empty);
        assert_eq!(rebuilt.populations, request.populations);
        assert!(!rebuilt.force_samples);
        assert_eq!(
            layout.invariant,
            reed::InvariantMode::VariantOnly,
            "the layout rides the builder return for the writers"
        );
        match (&rebuilt.grid, &request.grid) {
            (
                reed::WindowGrid::Bp { size: a, step: b },
                reed::WindowGrid::Bp { size: c, step: d },
            ) => assert_eq!((a, b), (c, d)),
            _ => panic!("grid must round-trip as Bp"),
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}
