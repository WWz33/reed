//! The `reed` binary: clap parsing plus dispatch (03 §7). Everything the
//! command line knows lives here; the driver library `reed` has no clap
//! dependency.

mod cli;
mod translate;
mod units;

use clap::Parser;

use cli::{Cli, Command};


fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    match Cli::parse().command {
        Command::Qc(args) => translate::qc::run_qc(args),
        Command::Filter(args) => translate::filter::run_filter_cli(args),
        Command::Stats(args) => translate::stats::run_stats_cli(args),
        Command::Index(args) => translate::index::run_index_cli(args),
    }
}
