//! `reed qc` dispatch (02-cli.md §1; the report itself is `reed-qc`).

use std::path::PathBuf;

use reed::qc_report::{run_qc_report_opened, QcOptions};

use crate::cli::QcCli;
use crate::units::input_path;

/// Runs the QC report.
///
/// # Errors
///
/// Driver and output failures.
pub(crate) fn run_qc(args: QcCli) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    eprintln!(
        "reed qc: warning: `reed qc` is deprecated and will be removed; use \
         `reed stats --summary` instead"
    );
    let QcCli { input, json, out, threads } = args;
    let threads = threads.unwrap_or_else(|| {
        std::thread::available_parallelism().expect("invariant: parallelism").get()
    });
    let started = std::time::Instant::now();
    let opened = reed::input::open(input_path(&input).as_deref())
        .map_err(|e| -> Box<dyn std::error::Error + Send + Sync> { e.to_string().into() })?;
    let outcome = run_qc_report_opened(
        opened,
        &QcOptions { threads, json: json.map(PathBuf::from) },
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
    Ok(())
}
