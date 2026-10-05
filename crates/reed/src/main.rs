//! `reed` CLI (`docs/architecture-spec.md` M0: clap; M4: `reed stats`).
//! The `filter` and `qc` subcommands follow with M5/M6.

use std::path::Path;

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
    let Command::Stats {
        input,
        groups,
        size,
        step,
        levels,
        skip_empty,
        out,
        threads,
    } = Cli::parse().command;

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
