//! LD writers (docs/rounds/v5/05-algorithms.md §3.4): `PREFIX.ld.decay`
//! (`bin_start bin_end n_pairs mean_r2`) and the plink2-style per-pair
//! `PREFIX.ld` (`CHROM_A POS_A ID_A CHROM_B POS_B ID_B UNPHASED_R2`).
//! Pruning moved to `reed filter --indep-pairwise` (docs/rounds/v5/05
//! §3.4.1).

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use reed_engine::error::EngineError;

use crate::stats::LdReport;

/// Writes the LD files the report carries; returns the paths written.
///
/// # Errors
///
/// [`EngineError::Io`] for write failures.
pub fn write_ld(prefix: &Path, report: &LdReport) -> Result<Vec<PathBuf>, EngineError> {
    let mut files = Vec::new();
    if !report.decay.is_empty() {
        let mut out = String::from("bin_start\tbin_end\tn_pairs\tmean_r2\n");
        for (start, end, n, mean) in &report.decay {
            let cell = match mean {
                Some(v) => format!("{v}"),
                None => "NA".to_string(),
            };
            writeln!(out, "{start}\t{end}\t{n}\t{cell}").expect("invariant: writing to String");
        }
        files.push(write(prefix, ".ld.decay", out)?);
    }
    if !report.pairs.is_empty() {
        let mut out = String::from("#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n");
        for (contig, a, b, r2) in &report.pairs {
            writeln!(out, "{contig}\t{a}\t.\t{contig}\t{b}\t.\t{r2}")
                .expect("invariant: writing to String");
        }
        files.push(write(prefix, ".ld", out)?);
    }
    Ok(files)
}

fn write(prefix: &Path, suffix: &str, text: String) -> Result<PathBuf, EngineError> {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(suffix);
    let path = PathBuf::from(path);
    std::fs::write(&path, text).map_err(EngineError::Io)?;
    Ok(path)
}
