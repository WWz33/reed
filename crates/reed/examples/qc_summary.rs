//! Runs `run_qc` (the QC `summary` section) over a BGZF VCF; the M3
//! end-to-end consumer used for the full-file acceptance run.
//!
//! Usage: `qc_summary <input.vcf.gz> [--threads N]`.

use reed::qc::run_qc;
use reed_engine::budget::ThreadBudget;

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut args = std::env::args().skip(1);
    let Some(input) = args.next() else {
        eprintln!("usage: qc_summary <input.vcf.gz> [--threads N]");
        std::process::exit(2);
    };
    let threads = match (args.next().as_deref(), args.next().as_deref()) {
        (Some("--threads"), Some(n)) => n.parse::<usize>()?,
        _ => std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get),
    };
    let budget = ThreadBudget {
        threads,
        mem: reed_engine::budget::MemBudget::default_half_memory(),
    };
    let started = std::time::Instant::now();
    let report = run_qc(std::path::Path::new(&input), &budget)?;
    eprintln!(
        "qc: n_records={} n_samples={} n_chunks={} elapsed={:.2?} (total {:.2?})",
        report.summary.n_records,
        report.n_samples,
        report.n_chunks,
        report.elapsed,
        started.elapsed()
    );
    println!(
        "{}",
        serde_json::to_string_pretty(&report.summary)
            .expect("invariant: summary serializes")
    );
    Ok(())
}
