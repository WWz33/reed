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
        // cli.md §9: usage errors exit 2 (the operator can fix the command
        // line), runtime failures exit 1. clap's own usage errors exit 2
        // inside `Cli::parse` and never reach this point.
        let code = match err.kind() {
            reed::ErrorKind::Usage => 2,
            reed::ErrorKind::Runtime => 1,
        };
        std::process::exit(code);
    }
}

fn run() -> Result<(), reed::Error> {
    match Cli::parse().command {
        Command::Filter(args) => translate::filter::run_filter_cli(args),
        Command::Stats(args) => translate::stats::run_stats_cli(args),
        Command::Index(args) => translate::index::run_index_cli(args),
        Command::Convert(args) => translate::convert::run_convert_cli(args),
    }
}
