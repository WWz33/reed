//! reed-filter CLI thin shell.
//!
//! Same philosophy as bcftools for the user-facing surface: zero strategy
//! flags, only resource and IO knobs. The execution strategy is chosen by
//! the engine (the ExecOptions/strategy parity hook is not exposed on the
//! CLI).
//!
//! Exit codes: 0 success; 1 runtime failure (FilterError); 2 usage error.
//! The RunReport (funnel) is written to stderr; stdout is the data stream.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{error::ErrorKind, CommandFactory, Parser};

use reed_filter::{emitter, exec, FilterError, FilterSpec};

#[derive(Debug, Parser)]
#[command(
    name = "reed-filter",
    version,
    about = "reed streaming hard-filter engine"
)]
struct Cli {
    /// Filter spec (gold filters.conf KEY=VALUE format; default = built-in
    /// default spec).
    #[arg(short = 'c', long = "conf")]
    conf: Option<PathBuf>,
    /// Output path; default is plain text VCF on stdout (.vcf.gz and .bcf
    /// are dispatched by file extension).
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,
    /// MT decompression/compression pool worker count (default = available
    /// parallelism).
    #[arg(long)]
    threads: Option<usize>,
    /// Input VCF/BCF path.
    input: PathBuf,
}

fn main() -> ExitCode {
    // Diagnostics go to stderr; stdout is the data stream. Verbosity via
    // RUST_LOG (e.g. RUST_LOG=info for stage timing, debug for per-stage).
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .with_writer(std::io::stderr)
        .init();
    match run() {
        Ok(code) => code,
        Err(e) => {
            eprintln!("reed-filter: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<ExitCode, FilterError> {
    run_from(std::env::args().skip(1))
}

/// Testable core of [`run`]: consumes an explicit argument iterator.
///
/// A usage error is an exit code 2; runtime failures propagate as
/// [`FilterError`] and become exit code 1 in main.
fn run_from<I>(args: I) -> Result<ExitCode, FilterError>
where
    I: IntoIterator<Item = String>,
{
    let argv: Vec<String> = args.into_iter().collect();
    if argv.is_empty() {
        // Legacy behaviour: no arguments means help, not a usage error.
        let _ = Cli::command().print_help();
        return Ok(ExitCode::SUCCESS);
    }
    let cli = match Cli::try_parse_from(std::iter::once("reed-filter".to_owned()).chain(argv)) {
        Ok(cli) => cli,
        Err(err) => {
            if matches!(err.kind(), ErrorKind::DisplayHelp | ErrorKind::DisplayVersion) {
                let _ = err.print();
                return Ok(ExitCode::SUCCESS);
            }
            eprintln!("reed-filter: {err}");
            return Ok(ExitCode::from(2));
        }
    };

    let spec = match cli.conf {
        Some(p) => FilterSpec::from_conf(&p)?,
        None => FilterSpec::default(),
    };

    let mut opts = exec::ExecOptions::default();
    if let Some(t) = cli.threads {
        opts.threads = t;
    }

    let mut em = emitter::open_sink(cli.output.as_deref(), opts.threads)?;
    let report = exec::execute(&cli.input, &spec, em.as_mut(), &opts)?;
    em.finish()?;

    // Funnel report to stderr (stdout is the data stream, never mixed in).
    eprintln!(
        "reed-filter: input={} sites_out={} samples_out={}/{} masked_gt={} elapsed={:.1}s",
        report.input.display(),
        report.n_output_records,
        report.n_samples_out,
        report.n_samples_in,
        report.n_masked_gt,
        report.elapsed.as_secs_f64(),
    );
    for (rule, kept) in &report.funnel {
        eprintln!("  funnel {rule:<16} {kept}");
    }
    Ok(ExitCode::SUCCESS)
}

#[cfg(test)]
mod tests {
    //! Parsing and pre-engine error paths. The end-to-end smoke (synthetic
    //! VCF.gz through run()) needs the record decoder and is covered by the
    //! gold parity gate instead.

    use super::*;

    fn parse(v: &[&str]) -> Result<Cli, clap::Error> {
        Cli::try_parse_from(std::iter::once("reed-filter".to_owned()).chain(v.iter().map(|s| (*s).to_owned())))
    }

    #[test]
    fn parse_positional_only() {
        let cli = parse(&["in.vcf.gz"]).ok().unwrap();
        assert_eq!(cli.input, PathBuf::from("in.vcf.gz"));
        assert!(cli.conf.is_none() && cli.output.is_none() && cli.threads.is_none());
    }

    #[test]
    fn parse_all_flags() {
        let cli = parse(&["-c", "conf.txt", "-o", "out.bcf", "--threads", "8", "in.bcf"])
            .ok()
            .unwrap();
        assert_eq!(cli.conf, Some(PathBuf::from("conf.txt")));
        assert_eq!(cli.output, Some(PathBuf::from("out.bcf")));
        assert_eq!(cli.threads, Some(8));
        assert_eq!(cli.input, PathBuf::from("in.bcf"));
    }

    #[test]
    fn parse_long_flags_equivalent() {
        let cli = parse(&[
            "--conf", "c.txt", "--output", "o.vcf.gz", "--threads", "4", "in.vcf.gz",
        ])
        .ok()
        .unwrap();
        assert_eq!(cli.conf, Some(PathBuf::from("c.txt")));
        assert_eq!(cli.output, Some(PathBuf::from("o.vcf.gz")));
        assert_eq!(cli.threads, Some(4));
    }

    #[test]
    fn parse_help_wins_over_everything() {
        assert_eq!(parse(&["-h"]).err().unwrap().kind(), ErrorKind::DisplayHelp);
        assert_eq!(parse(&["--help"]).err().unwrap().kind(), ErrorKind::DisplayHelp);
        assert_eq!(
            parse(&["in.vcf.gz", "-h"]).err().unwrap().kind(),
            ErrorKind::DisplayHelp
        );
    }

    #[test]
    fn parse_missing_flag_value() {
        assert!(parse(&["-c"]).is_err());
        assert!(parse(&["--conf"]).is_err());
        assert!(parse(&["-o"]).is_err());
        assert!(parse(&["--output"]).is_err());
        assert!(parse(&["--threads"]).is_err());
    }

    #[test]
    fn parse_unknown_flag() {
        assert!(parse(&["--bogus", "in.vcf.gz"]).is_err());
        assert!(parse(&["-x", "in.vcf.gz"]).is_err());
    }

    #[test]
    fn parse_extra_positional() {
        assert!(parse(&["a.vcf.gz", "b.vcf.gz"]).is_err());
    }

    #[test]
    fn parse_threads_not_a_number() {
        assert!(parse(&["--threads", "abc", "in.vcf.gz"]).is_err());
    }

    #[test]
    fn run_usage_error_exits_two() {
        let code = run_from(["--bogus".to_owned()]).ok().unwrap();
        assert_eq!(code, ExitCode::from(2));
        let code = run_from(["a.vcf.gz".to_owned(), "b.vcf.gz".to_owned()])
            .ok()
            .unwrap();
        assert_eq!(code, ExitCode::from(2));
    }

    #[test]
    fn run_no_args_prints_help_with_success() {
        let code = run_from([]).ok().unwrap();
        assert_eq!(code, ExitCode::SUCCESS);
    }

    #[test]
    fn run_missing_conf_is_runtime_failure() {
        // Stops at FilterSpec::from_conf, well before the engine runs.
        let e = run_from([
            "-c".to_owned(),
            "/nonexistent/reed-filter-test.conf".to_owned(),
            "in.vcf.gz".to_owned(),
        ])
        .err()
        .unwrap();
        assert!(matches!(e, FilterError::Io(_)));
    }
}
