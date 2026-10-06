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
        /// Input BGZF VCF (.vcf.gz)
        input: String,
        /// Groups file: `sample<TAB>population` rows, no header (pixy
        /// populations-file format)
        #[arg(short, long)]
        groups: String,
        /// Skip group-file rows naming samples absent from the header
        /// (default: error, pixy behavior)
        #[arg(long)]
        force_samples: bool,
        /// Contig lengths as a .fai-style file (name, length), for inputs
        /// whose header lacks ##contig lines
        #[arg(long)]
        contig_lengths_path: Option<String>,
        /// Window size in bp
        #[arg(long, default_value_t = 100_000)]
        size: u64,
        /// Window step in bp
        #[arg(long, default_value_t = 50_000)]
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

    let request = StatsRequest {
        kinds: default_kinds(group_names.len()),
        gt_groups,
        window_size: size,
        step,
        levels: Levels::parse(&levels)?,
        skip_empty,
        contig_lengths,
    };
    let threads = threads
        .or_else(|| std::env::var("REED_THREADS").ok().and_then(|t| t.parse().ok()))
        .unwrap_or_else(|| {
            std::thread::available_parallelism().expect("invariant: parallelism").get()
        });
    let budget = ThreadBudget { threads, mem: MemBudget::default_half_memory() };

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
