//! The `reed` binary: clap parsing plus dispatch (03 §7). Everything the
//! command line knows lives here; the driver library `reed` has no clap
//! dependency.

mod cli;
mod translate;
mod units;

use clap::{CommandFactory, FromArgMatches};

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

/// Parses the command line with the grouped `--help` (doc_v14 02 §1.3): the
/// group list is dynamic (`GROUPS` + clap's own `about` texts), so it rides
/// `after_help` on the command clap parses with.
fn parse_cli() -> Cli {
    let command = Cli::command().after_help(cli::render_grouped_commands());
    let matches = command.get_matches();
    Cli::from_arg_matches(&matches).expect("invariant: matches come from Cli::command()")
}

fn run() -> Result<(), reed::Error> {
    match parse_cli().command {
        Command::Filter(args) => translate::filter::run_filter_cli(args),
        Command::Stats(args) => translate::stats::run_stats_cli(args),
        Command::Popg(args) => translate::popg::run_popg_cli(args),
        Command::Index(args) => translate::index::run_index_cli(args),
        Command::Convert(args) => translate::convert::run_convert_cli(args),
        Command::Rename(args) => translate::rename::run_rename_cli(args),
        Command::Pca(args) => translate::pca::run_pca_cli(args),
        Command::Distance(args) => translate::distance::run_distance_cli(args),
    }
}
