//! `reed` CLI (`docs/architecture-spec.md` M0: clap; M4: `reed stats`).
//! The `filter` and `qc` subcommands follow with M5/M6.

use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};
use reed::stats::{assign_groups, default_kinds, run_stats, Levels, StatsRequest};
use reed_engine::budget::{MemBudget, ThreadBudget};
use reed_io::record::record_header;

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
        /// Write a tabix index next to the output
        #[arg(long)]
        write_index: bool,
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
        /// Groups file: `sample<TAB>population` rows, first line skipped
        #[arg(short, long)]
        groups: String,
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

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let command = Cli::parse().command;
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

    let header = record_header(Path::new(&input))?;
    let sample_names: Vec<String> =
        header.sample_names().iter().map(ToString::to_string).collect();
    let (gt_groups, group_names) = assign_groups(Path::new(&groups), &sample_names)?;
    eprintln!(
        "{} samples, {} groups: {}",
        sample_names.len(),
        group_names.len(),
        group_names.join(", ")
    );

    let request = StatsRequest {
        kinds: default_kinds(group_names.len()),
        gt_groups,
        window_size: size,
        step,
        levels: Levels::parse(&levels)?,
        skip_empty,
    };
    let threads = threads
        .or_else(|| std::env::var("REED_THREADS").ok().and_then(|t| t.parse().ok()))
        .unwrap_or_else(|| {
            std::thread::available_parallelism().expect("invariant: parallelism").get()
        });
    let budget = ThreadBudget { threads, mem: MemBudget::from_bytes(u64::MAX) };

    let report = run_stats(Path::new(&input), request, &budget)?;
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
    write_index: bool,
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
        write_index,
        compress_level,
        strategy,
        command_line,
        quiet: false,
        report_json: report.map(PathBuf::from),
    };
    let report = reed::filter::run_filter(Path::new(&input), &spec_v2, &opts)?;
    eprintln!(
        "reed filter: {} -> {} records in {:.1}s",
        report.n_input_records, report.n_output_records, report.elapsed_secs
    );
    Ok(())
}
