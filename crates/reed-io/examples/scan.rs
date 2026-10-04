//! Smoke scanner: stream a VCF/BCF, print per-window pi and Fst between the
//! first and second halves of the samples.
//!
//! Usage: `scan <vcf|bcf> [window_size]`

use std::path::PathBuf;

use reed_core::estimators::fst::{fst, FstEstimator};
use reed_core::estimators::pi::{pi, PiEstimator};
use reed_core::{Ascertainment, PopulationView, SiteCompleteness};
use reed_io::{for_each_window, header_info, ReadError, ReadOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(
        args.next()
            .ok_or_else(|| "usage: scan <vcf|bcf> [window_size]".to_string())?,
    );
    let window_size: u32 = match args.next() {
        Some(s) => s.parse()?,
        None => 100_000,
    };

    let info = header_info(&path)?;
    let n = info.n_samples;
    let missing_lengths = info
        .contigs
        .iter()
        .filter(|(_, l)| l.is_none())
        .count();
    println!(
        "samples={n} contigs={} contigs_without_length={missing_lengths}",
        info.contigs.len()
    );
    if n < 2 {
        return Err("need at least 2 samples".into());
    }

    let half = n / 2;
    let assignment: Vec<u32> = (0..n).map(|i| if i < half { 0 } else { 1 }).collect();
    let pops = PopulationView::from_assignment(n, &assignment, 2)?;

    let opts = ReadOptions {
        window_size,
        step: window_size,
        manifest: SiteCompleteness::VariantOnly(Ascertainment::Complete),
    };

    let mut printed = 0usize;
    let summary = for_each_window(&path, &opts, |w| {
        let pi_a = pi(&w, &pops, 0, PiEstimator::PerPairComparable)
            .map_err(|e| ReadError::User(e.to_string()))?;
        let f = fst(&w, &pops, 0, 1, FstEstimator::Hudson)
            .map_err(|e| ReadError::User(e.to_string()))?;
        if printed < 10 {
            println!(
                "window contig={} [{}, {}) sites={} piA={:?} fst={:?}",
                w.contig.0,
                w.start,
                w.end,
                w.n_sites(),
                pi_a.value(),
                f.value()
            );
            printed += 1;
        }
        Ok(())
    })?;

    println!("windows={} sites={}", summary.windows, summary.sites);
    Ok(())
}
