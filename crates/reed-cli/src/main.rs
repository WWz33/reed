//! The `reed` binary: clap parsing plus dispatch (03 §7). Everything the
//! command line knows lives here; the driver library `reed` has no clap
//! dependency.

mod cli;
mod translate;
mod units;

use clap::Parser;

use cli::{Cli, Command};

fn main() {
    if let Err(err) = run() {
        // Display, not Debug: a `String` error would otherwise print its
        // Debug form with surrounding quotes.
        eprintln!("reed: {err}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    match Cli::parse().command {
        Command::Filter(args) => translate::filter::run_filter_cli(args),
        Command::Stats(args) => translate::stats::run_stats_cli(args),
        Command::Index(args) => translate::index::run_index_cli(args),
    }
}
