//! Shared translation pieces for `reed stats` and `reed popg` (doc_v14
//! 02 §2.5): the region/sample/site flag parsing both commands share, so the
//! two front ends cannot drift (04-cleanup §2.1).

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use reed::budget::{MemBudget, ThreadBudget};

use crate::units::input_path;

/// Prints the engine's stage metrics (doc_v4 03 A4-5): the reader, reducer
/// and frame-chain cumulative times without which "reader-bound" is only an
/// inference.
pub(crate) fn print_engine_metrics(m: &reed::metrics::StageMetrics) {
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
pub(crate) fn parse_sample_keep(
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
            let text = std::fs::read_to_string(p)
                .map_err(|e| format!("cannot read sample file {p}: {e}"))?;
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
            let text = std::fs::read_to_string(p)
                .map_err(|e| format!("cannot read sample file {p}: {e}"))?;
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


/// The parsed-and-opened state both statistics commands need before they
/// diverge (regions resolved, the input opened once, the shared site/sample
/// flags compiled, the thread budget built).
pub(crate) struct Prelude {
    pub regions: Vec<reed::chunk::region::Region>,
    pub targets: Vec<reed::spec::region::RegionTriple>,
    pub targets_exclude: Vec<reed::spec::region::RegionTriple>,
    pub opened: Option<reed::OpenedInput>,
    pub sample_names: Vec<String>,
    pub keep_samples: Option<Vec<bool>>,
    pub site: Option<reed::spec::FilterSpec>,
    pub contig_lengths: HashMap<Vec<u8>, u64>,
    pub budget: ThreadBudget,
    pub open_secs: f64,
}

/// The shared prologue (doc_v14 02 §2.5): parse `-r/-R/-t/-T`, open the
/// input once (stdin streams, 03-performance.md O5), resolve `-s/-S/
/// --keep/--remove`, compile the site/threshold flags through the filter's
/// own compiler, load the contig lengths and build the thread budget.
///
/// # Errors
///
/// Usage failures from any of the flag parsers, and I/O failures opening
/// the input.
#[allow(clippy::too_many_arguments)]
pub(crate) fn prelude(
    verb: &str,
    input: &str,
    regions_args: &crate::cli::filter_flags::RegionArgs,
    sample_args: &crate::cli::filter_flags::SampleArgs,
    site_args: &crate::cli::filter_flags::SiteFilterArgs,
    sites_file: Option<&str>,
    contig_lengths_arg: Option<&str>,
    threads_arg: Option<usize>,
    memory_arg: Option<u64>,
) -> Result<Prelude, reed::Error> {
    let started = std::time::Instant::now();
    let (regions, mut targets, targets_exclude) =
        parse_regions(regions_args).map_err(reed::Error::usage)?;
    // Indexed ranges are coarse (whole blocks), so the same regions also act
    // as the record predicate — exactly how the filter's `-r` works.
    for region in &regions {
        targets.push((region.contig.clone(), region.start1, region.end1));
    }
    let input_file = input_path(input);
    let open_path: Option<PathBuf> = input_file.clone();
    let opened = Some(reed::open_regions(open_path.as_deref(), &regions)?);
    let open_secs = started.elapsed().as_secs_f64();
    tracing::debug!(open_ms = open_secs * 1e3, input, "{verb} opened the input");
    let facts = opened.as_ref().expect("invariant: just opened").facts();
    let sample_names: Vec<String> = facts.sample_names.clone();
    let keep_samples =
        parse_sample_keep(sample_args, &sample_names).map_err(reed::Error::usage)?;
    if !site_args.indep_pairwise.is_empty() {
        return Err(reed::Error::usage(format!(
            "{verb} does not support --indep-pairwise yet (the pruning needs a second pass);              run `reed filter --indep-pairwise W[kb] [S] R2` first and feed its output to {verb}"
        )));
    }
    let site_spec = super::filter::site_spec(site_args).map_err(reed::Error::usage)?;
    // No biallelic gate here: the reporters decide per statistic how to treat
    // multiallelic sites (doc_v3 05-roadmap R7); gating records would
    // silently change the summary counts. Compiling the typed spec against
    // the header happens inside the driver (08 §6.2).
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
    if let Some(path) = contig_lengths_arg {
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
                reed::Error::usage(format!("{path}: contig length for {name:?} is not a number"))
            })?;
            contig_lengths.insert(name.as_bytes().to_vec(), len);
        }
    }

    let threads = threads_arg.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .expect("invariant: parallelism")
            .get()
    });
    let budget = ThreadBudget {
        threads,
        mem: memory_arg.map_or_else(MemBudget::default_half_memory, |m| {
            MemBudget::from_bytes(m.saturating_mul(1 << 20))
        }),
    };

    // `--sites-file` narrows the counted sites for every route (doc_v7 D7-6).
    if let Some(path) = sites_file {
        for (contig, pos1) in parse_sites(path).map_err(reed::Error::usage)? {
            targets.push((contig, Some(pos1), Some(pos1)));
        }
    }

    Ok(Prelude {
        regions,
        targets,
        targets_exclude,
        opened,
        sample_names,
        keep_samples,
        site,
        contig_lengths,
        budget,
        open_secs,
    })
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

