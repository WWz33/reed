//! `reed` CLI (`docs/architecture-spec.md` M0: clap; M4: `reed stats`).
//! The `filter` and `qc` subcommands follow with M5/M6.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};
use reed::stats::{assign_groups, default_kinds, run_stats, Levels, StatsRequest};
use reed_engine::budget::{MemBudget, ThreadBudget};

#[derive(Parser)]
#[command(
    name = "reed",
    version,
    about = "Population-genetics statistics and hard filtering, in pure Rust"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Quality-control report (§10): bcftools-stats-compatible output plus
    /// the structured JSON report.
    Qc {
        /// Input BGZF VCF (.vcf.gz)
        input: String,
        /// Write the structured JSON report here
        #[arg(long)]
        json: Option<String>,
        /// Output path for the bcftools-stats text (omit for stdout)
        #[arg(short, long)]
        out: Option<String>,
        /// Worker threads (default: all cores)
        #[arg(long)]
        threads: Option<usize>,
    },
    /// Streaming hard filter (v2 engine, §11): gold-parity funnel with
    /// worker-local BGZF compression and synchronized indexing.
    Filter {
        /// Input BGZF VCF (.vcf.gz)
        input: String,
        /// Filter spec as TOML (§11.1)
        #[arg(long, conflicts_with = "conf")]
        spec: Option<String>,
        /// Filter spec as the gold KEY=VALUE conf (compat path)
        #[arg(long)]
        conf: Option<String>,
        /// Output VCF (.vcf.gz) — omit for plain VCF on the path/stdout
        #[arg(short, long)]
        out: Option<String>,
        /// Write an index next to the output: csi (default) or tbi
        #[arg(short = 'W', long = "write-index",
              value_name = "csi|tbi", num_args = 0..=1, default_missing_value = "csi")]
        write_index: Option<String>,
        /// BGZF compression level (default 6)
        #[arg(long, default_value_t = 6)]
        compress_level: u8,
        /// JSON report path
        #[arg(long)]
        report: Option<String>,
        /// Test-only strategy override (single_pass|two_pass_reparse)
        #[arg(long, hide = true)]
        strategy: Option<String>,
        /// Worker threads (default: all cores)
        #[arg(long)]
        threads: Option<usize>,
    },

    /// Windowed statistics (π, dxy, Fst, PBS, SFS, θ_W, Tajima's D) over a
    /// BGZF VCF, as additive accumulators with byte-deterministic output.
    Stats {
        /// Input BGZF VCF (.vcf.gz); `-` for stdin
        input: String,
        /// Groups file: `sample<TAB>population` rows, no header (pixy
        /// populations-file format)
        #[arg(short, long, alias = "pop", alias = "populations")]
        groups: String,
        /// π per population (Nei & Li 1979; missing data: Korunes & Samuk
        /// 2021 = pixy)
        #[arg(long)]
        pi: bool,
        /// dxy per population pair (Nei 1987; missing data: pixy)
        #[arg(long)]
        dxy: bool,
        /// Watterson's θ per population (Watterson 1975; per-allele-class
        /// mutation counts, pixy)
        #[arg(long)]
        theta: bool,
        /// Tajima's D per population (Tajima 1989; denominator once at the
        /// window-average allele count, pixy)
        #[arg(long = "tajima-d")]
        tajima_d: bool,
        /// Fst for all population pairs; `wc84` (Weir & Cockerham 1984,
        /// default) or `hudson` (Hudson 1992 / Bhatia 2013)
        #[arg(long, value_name = "wc84|hudson", num_args = 0..=1,
              default_missing_value = "wc84")]
        fst: Option<String>,
        /// PBS for the first three populations (Yi et al. 2010); needs
        /// three groups
        #[arg(long)]
        pbs: bool,
        /// Report PBSn1 (Malaspinas et al. 2016) instead of the raw PBS
        #[arg(long)]
        pbs_normed: bool,
        /// Declare an all-variant input: window positions the file does not
        /// contain count as monomorphic with all samples called
        #[arg(long)]
        variant_only: bool,
        /// Compute over the sites the file contains (overestimates π/dxy)
        #[arg(long)]
        bypass_invariant_check: bool,
        /// Skip group-file rows naming samples absent from the header
        /// (default: error, pixy behavior)
        #[arg(long)]
        force_samples: bool,
        /// Contig lengths as a .fai-style file (name, length), for inputs
        /// whose header lacks ##contig lines
        #[arg(long)]
        contig_lengths_path: Option<String>,
        /// Window size in bp (accepts 100k, 1M, 100kb)
        #[arg(long = "window", alias = "size", default_value = "100000", value_parser = parse_size)]
        size: u64,
        /// Window step in bp (default: the window size)
        #[arg(long, default_value = "0", value_parser = parse_size)]
        step: u64,
        /// Output levels, comma-separated: window,contig,genome
        #[arg(long, default_value = "window")]
        levels: String,
        /// Skip windows that no site fell into
        #[arg(long)]
        skip_empty: bool,
        /// Write the TSV here instead of stdout
        #[arg(short, long)]
        out: Option<String>,
        /// Worker threads (default: all cores)
        #[arg(long)]
        threads: Option<usize>,
    },
}

/// Length parser for window/step: plain bp, or k/kb/M/Mb suffixes
/// (02-cli.md §2).
fn parse_size(text: &str) -> Result<u64, String> {
    let lower = text.trim().to_ascii_lowercase();
    let (digits, scale) = if let Some(rest) = lower.strip_suffix("kb") {
        (rest, 1_000_u64)
    } else if let Some(rest) = lower.strip_suffix('k') {
        (rest, 1_000)
    } else if let Some(rest) = lower.strip_suffix("mb") {
        (rest, 1_000_000)
    } else if let Some(rest) = lower.strip_suffix('m') {
        (rest, 1_000_000)
    } else if let Some(rest) = lower.strip_suffix("bp") {
        (rest, 1)
    } else {
        (lower.as_str(), 1)
    };
    let value: u64 = digits
        .trim()
        .parse()
        .map_err(|_| format!("{text:?} is not a length (use bp, k/kb, M/Mb)"))?;
    value
        .checked_mul(scale)
        .ok_or_else(|| format!("{text:?} is too large"))
}

/// `-` means stdin; everything else is a file path.
fn input_path(arg: &str) -> Option<PathBuf> {
    (arg != "-").then(|| PathBuf::from(arg))
}

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let command = Cli::parse().command;
    if let Command::Qc { input, json, out, threads } = &command {
        let threads = threads
            .or_else(|| std::env::var("REED_THREADS").ok().and_then(|t| t.parse().ok()))
            .unwrap_or_else(|| {
                std::thread::available_parallelism().expect("invariant: parallelism").get()
            });
        let started = std::time::Instant::now();
        let outcome = reed::qc_report::run_qc_report(
            input_path(input).as_deref(),
            &reed::qc_report::QcOptions {
                threads,
                json: json.clone().map(PathBuf::from),
            },
        )?;
        match out {
            Some(path) => std::fs::write(path, &outcome.bcftools_text)?,
            None => print!("{}", outcome.bcftools_text),
        }
        eprintln!(
            "reed qc: {} records, {} samples in {:.1}s",
            outcome.report.summary.n_records,
            outcome.report.n_samples,
            started.elapsed().as_secs_f64()
        );
        return Ok(());
    }
    if let Command::Filter {
        input, spec, conf, out, write_index, compress_level, report, strategy, threads,
    } = command
    {
        return run_filter_cli(FilterCli {
            input,
            spec,
            conf,
            out,
            write_index,
            compress_level,
            report,
            strategy,
            threads,
        });
    }
    let Command::Stats {
        input,
        groups,
        pi,
        dxy,
        theta,
        tajima_d,
        fst,
        pbs,
        pbs_normed,
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
    } = command
    else {
        unreachable!("the filter subcommand returned above");
    };

    let (_source, facts) = reed::input::open(input_path(&input).as_deref())
        .map_err(|e| -> Box<dyn std::error::Error + Send + Sync> { e.to_string().into() })?;
    let sample_names: Vec<String> = facts.sample_names.clone();
    let (gt_groups, group_names) =
        assign_groups(Path::new(&groups), &sample_names, force_samples)?;
    eprintln!(
        "{} samples, {} groups: {}",
        sample_names.len(),
        group_names.len(),
        group_names.join(", ")
    );

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

    let threads = threads
        .or_else(|| std::env::var("REED_THREADS").ok().and_then(|t| t.parse().ok()))
        .unwrap_or_else(|| {
            std::thread::available_parallelism().expect("invariant: parallelism").get()
        });
    let budget = ThreadBudget { threads, mem: MemBudget::default_half_memory() };

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

struct FilterCli {
    input: String,
    spec: Option<String>,
    conf: Option<String>,
    out: Option<String>,
    write_index: Option<String>,
    compress_level: u8,
    report: Option<String>,
    strategy: Option<String>,
    threads: Option<usize>,
}

fn run_filter_cli(
    FilterCli {
        input,
        spec,
        conf,
        out,
        write_index,
        compress_level,
        report,
        strategy,
        threads,
    }: FilterCli,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let spec_v2 = match (spec, conf) {
        (Some(path), _) => reed::filter::read_spec_toml(Path::new(&path))?,
        (None, Some(path)) => reed_filter::spec::compat::from_gold_conf(Path::new(&path))?,
        (None, None) => {
            return Err("one of --spec <filter.toml> or --conf <filters.conf> is required".into());
        }
    };
    let strategy = strategy
        .map(|s| match s.as_str() {
            "single_pass" => Ok(reed_filter::plan::ExecStrategy::SinglePass),
            "two_pass_reparse" => Ok(reed_filter::plan::ExecStrategy::TwoPassReparse),
            "sketch_then_emit" => Ok(reed_filter::plan::ExecStrategy::SketchThenEmit),
            other => Err(format!("unknown strategy {other:?}")),
        })
        .transpose()?;
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let threads = threads
        .or_else(|| std::env::var("REED_THREADS").ok().and_then(|t| t.parse().ok()))
        .unwrap_or_else(|| {
            std::thread::available_parallelism().expect("invariant: parallelism").get()
        });
    let opts = reed::filter::FilterOptions {
        threads,
        output: out.map(PathBuf::from),
        write_index: reed::filter::WriteIndex::parse(write_index.as_deref())?,
        compress_level,
        strategy,
        command_line,
        quiet: false,
        report_json: report.map(PathBuf::from),
    };
    let report = reed::filter::run_filter(input_path(&input).as_deref(), &spec_v2, &opts)?;
    eprintln!(
        "reed filter: {} -> {} records in {:.1}s",
        report.n_input_records, report.n_output_records, report.elapsed_secs
    );
    Ok(())
}
