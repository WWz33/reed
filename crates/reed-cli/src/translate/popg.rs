//! `reed popg` flag translation (doc_v14 02 §2.3): the windowed
//! population-genetics statistics. Rows are per population (pair/triple) per
//! window or genome-wide; the QC reports live in `reed stats`.

use std::collections::HashMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use reed::{Levels, Populations, PopgenRequest, default_kinds, run_popgen_opened};

use crate::cli::PopgCli;

use super::common::{Prelude, prelude, print_engine_metrics};

/// Runs the windowed-statistics driver for the parsed flags.
///
/// # Errors
///
/// Translation (unknown levels/sizes) and driver failures.
pub(crate) fn run_popg_cli(args: PopgCli) -> Result<(), reed::Error> {
    args.verbosity.init_tracing();
    // Usage checks before anything is read (cli.md §0.7): the window-layout
    // route always writes PREFIX files, so a missing `-o` is refused here,
    // not after the whole run (doc_v15 D15-13).
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
    let writes_files =
        windowed_route || args.sfs.is_some() || args.ld_decay.is_some() || args.ld_pairs;
    // doc_v21 S17: four flags were silently dropped when their companion was
    // missing; `convert` already rejects the same shapes. Check before any
    // input is read.
    if args.proj.is_some() && args.sfs.is_none() {
        return Err(reed::Error::usage("--proj needs --sfs"));
    }
    if args.sfs_int && args.sfs.is_none() {
        return Err(reed::Error::usage("--sfs-int needs --sfs"));
    }
    if args.block_size.is_some()
        && args.f2.is_empty()
        && args.f3.is_empty()
        && args.f4.is_empty()
        && args.dstat.is_empty()
        && args.f4_ratio.is_empty()
    {
        return Err(reed::Error::usage(
            "--block-size needs an f-statistic (--f2/--f3/--f4/--dstat/--f4-ratio)",
        ));
    }
    if args.pbs_normed && args.pbs.is_empty() {
        return Err(reed::Error::usage("--pbs-normed needs --pbs"));
    }
    if writes_files && args.out.is_none() {
        return Err(reed::Error::usage(
            "window-layout mode writes files: pass -o PREFIX",
        ));
    }
    let Prelude {
        regions,
        targets,
        targets_exclude,
        mut opened,
        sample_names,
        keep_samples,
        site,
        contig_lengths,
        budget,
        open_secs,
    } = prelude(
        "reed popg",
        &args.input,
        &args.regions_args,
        &args.sample_args,
        &args.site_args,
        args.sites_file.as_deref(),
        args.contig_lengths.as_deref(),
        args.threads,
        args.memory,
    )?;

    // Windowed statistics need the populations file; a bare `--sfs` (no
    // windowed statistic) treats every sample as one population.
    // A bare `--sfs`, `--ld-decay` or `--ld-pairs` treats every sample as
    // one population (only the window rows need groups).
    let populations = match &args.pop {
        Some(groups) => Populations::File(PathBuf::from(groups)),
        None => Populations::All,
    };
    // The names for kinds/fstats/sfs construction; resolved with the header
    // by the run's own parser so `--force-samples` cannot shift the list
    // (doc_v21 M10).
    let group_names = populations.group_names(&sample_names, args.force_samples)?;
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

    // Window-layout mode: one command writes PREFIX_pi.txt / PREFIX_dxy.txt
    // with the reference column layout (02-cli.md §5). LD is genome-wide (no
    // window rows), but it rides the same job as the windowed statistics; a
    // bare `--ld-decay` takes this path with no window kinds and no `--pop`
    // (F6-9). `--ld-pairs` needs `--ld-decay` (the builder refuses it alone).
    if writes_files {
        let ctx = WindowedContext {
            regions: regions.clone(),
            targets,
            targets_exclude,
            keep_samples,
            site,
            populations: populations.clone(),
            sample_names: sample_names.clone(),
            force_samples: args.force_samples,
            contig_lengths,
            grid,
        };
        let (request, group_names, kinds, layout) = build_windowed_request(&args, ctx)?;
        let source = opened
            .take()
            .expect("invariant: the windowed mode reads once");
        let run_start = std::time::Instant::now();
        let report = run_popgen_opened(source, &request, &budget)?;
        let run_secs = run_start.elapsed().as_secs_f64();
        let prefix = args
            .out
            .clone()
            .expect("invariant: -o checked before the run");
        if prefix == "-" {
            return Err(reed::Error::usage(
                "window-layout mode writes PREFIX files: -o - is not stdout",
            ));
        }
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
            eprintln!("# reed popg timing (windowed)");
            eprintln!("phase.open_ms\t{:.3}", open_secs * 1e3);
            eprintln!("phase.run_ms\t{:.3}", run_secs * 1e3);
            eprintln!("n_records\t{}", report.n_records);
            eprintln!("n_windows\t{}", report.n_windows);
            eprintln!("chunks\t{}", report.n_chunks);
            print_engine_metrics(&report.metrics);
        }
        return Ok(());
    }

    // The bare route counts pi/dxy too: honor the same invariant-site
    // contract as the windowed route (doc_v15 D15-15).
    let invariant = if args.variant_only {
        reed::InvariantMode::VariantOnly
    } else if args.bypass_invariant_check {
        eprintln!(
            "reed popg: warning: --bypass-invariant-check computes over present sites \
                 only; π and dxy are overestimated"
        );
        reed::InvariantMode::Bypass
    } else {
        reed::InvariantMode::Require
    };
    let kinds = default_kinds(group_names.len());
    let request = PopgenRequest {
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
        levels: Levels::parse(&args.levels)?,
        skip_empty: args.skip_empty,
        contig_lengths,
        layout: Some(reed::WindowLayoutOptions { invariant }),
    };
    let source = opened
        .take()
        .expect("invariant: the windowed mode reads once");
    let run_start = std::time::Instant::now();
    let report = run_popgen_opened(source, &request, &budget)?;
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
        eprintln!("# reed popg timing (windowed)");
        eprintln!("phase.open_ms\t{:.3}", open_secs * 1e3);
        eprintln!("phase.run_ms\t{:.3}", run_secs * 1e3);
        eprintln!("n_records\t{}", report.n_records);
        eprintln!("n_windows\t{}", report.n_windows);
        eprintln!("chunks\t{}", report.n_chunks);
        print_engine_metrics(&report.metrics);
    }
    let text = reed::windows::render_tsv(&report.windows, &kinds, false);
    match args.out.as_deref() {
        Some(path) if path != "-" => std::fs::write(path, text)?,
        // `-o -` and no `-o` both stream to stdout; write_all turns a closed
        // pipe into a normal error instead of a panic (doc_v15 D15-14/D15-39).
        _ => std::io::stdout().write_all(text.as_bytes())?,
    }
    Ok(())
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
        if ids[0] == ids[1] || ids[0] == ids[2] {
            return Err(format!(
                "--pbs {triple:?}: population {:?} repeats (want three distinct populations)",
                names[0]
            ));
        }
        if ids[1] == ids[2] {
            return Err(format!(
                "--pbs {triple:?}: population {:?} repeats (want three distinct populations)",
                names[1]
            ));
        }
        triples.push((ids[0], ids[1], ids[2]));
    }
    Ok(triples)
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
    /// Header sample names (doc_v21 M10: group names are resolved by the same
    /// parser the run uses, which needs the header).
    pub sample_names: Vec<String>,
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
    args: &PopgCli,
    ctx: WindowedContext,
) -> Result<
    (
        PopgenRequest,
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
    let group_names = ctx
        .populations
        .group_names(&ctx.sample_names, ctx.force_samples)?;
    let pbs_triples = parse_pbs_triples(&args.pbs, &group_names).map_err(reed::Error::usage)?;
    // A bare `--sfs` has no window denominators, so the all-sites check
    // does not apply (05 §3.1 accepts variant-only input).
    let ld_requested = args.ld_decay.is_some();
    // π/dxy/θ count the all-sites denominator; the f-statistic / LD
    // variant-only convention must not silently switch it when both are
    // requested in one run (doc_v15 D15-18). A bare `--sfs` keeps its
    // documented variant-only acceptance (05 §3.1).
    let all_sites_kinds = args.pi || args.dxy || args.theta || args.tajima_d;
    let invariant = if args.variant_only
        || (args.sfs.is_some() && !windowed)
        || ((fstats_requested || ld_requested) && !all_sites_kinds)
    {
        reed::InvariantMode::VariantOnly
    } else if args.bypass_invariant_check {
        eprintln!(
            "reed popg: warning: --bypass-invariant-check computes over present sites \
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
                "reed popg: note: two populations: --fst-global equals the pairwise wc84 value"
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
    // population (docs/rounds/v5/05-algorithms.md §3.1). The request
    // construction is shared with `reed convert` (doc_v13 D13-6).
    let sfs_request =
        super::build_sfs_request(args.sfs.as_deref(), args.proj.as_deref(), &group_names)?;
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
        Levels::parse(&args.levels)?
    };
    let request = PopgenRequest {
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
///
/// # Errors
///
/// `--ld-pairs` without `--ld-decay`: the pair table is that table's detail
/// level, and alone the flag used to produce no output at all.
fn build_ld_request(
    ld_decay: Option<&str>,
    ld_pairs: bool,
) -> Result<Option<reed::LdRequest>, reed::Error> {
    if ld_decay.is_none() {
        if ld_pairs {
            return Err(reed::Error::usage(
                "--ld-pairs writes the per-pair table of --ld-decay; pass --ld-decay [MAX_BP] \
                 as well",
            ));
        }
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
        let request = reed::PopgenRequest {
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
        let mut argv: Vec<String> = vec!["reed".into(), "popg".into(), "in.vcf".into()];
        argv.extend(request.to_args(&names()).expect("to_args"));
        let cli = crate::cli::Cli::try_parse_from(argv).expect("to_args must parse");
        let crate::cli::Command::Popg(popg) = cli.command else {
            panic!("popg subcommand");
        };
        let ctx = WindowedContext {
            regions: Vec::new(),
            targets: Vec::new(),
            targets_exclude: Vec::new(),
            keep_samples: None,
            site: None,
            populations: reed::Populations::File(pops),
            sample_names: vec!["s0".to_string(), "s1".to_string(), "s2".to_string()],
            force_samples: false,
            contig_lengths: std::collections::HashMap::new(),
            grid: reed::WindowGrid::Bp {
                size: 100_000,
                step: 50_000,
            },
        };
        let (rebuilt, rebuilt_names, rebuilt_kinds, layout) =
            build_windowed_request(&popg, ctx).expect("rebuild");
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
