//! LD writers (docs/rounds/v5/05-algorithms.md §3.4): `PREFIX.ld.decay`
//! (`bin_start bin_end n_pairs mean_r2`), the plink2-style per-pair
//! `PREFIX.ld` (`CHROM_A POS_A ID_A CHROM_B POS_B ID_B UNPHASED_R2`) and the
//! pruning lists `PREFIX.prune.in` / `PREFIX.prune.out` (one `CHROM_POS`
//! variant id per line).

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
            out.push_str(&format!("{start}\t{end}\t{n}\t{cell}\n"));
        }
        files.push(write(prefix, ".ld.decay", out)?);
    }
    if !report.pairs.is_empty() {
        let mut out = String::from(
            "#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n",
        );
        for (contig, a, b, r2) in &report.pairs {
            out.push_str(&format!("{contig}\t{a}\t.\t{contig}\t{b}\t.\t{r2}\n"));
        }
        files.push(write(prefix, ".ld", out)?);
    }
    if !report.prune_in.is_empty() || !report.prune_out.is_empty() {
        let key = |(contig, pos): &(String, u64)| format!("{contig}_{pos}");
        let in_text: String = report.prune_in.iter().map(|k| format!("{}\n", key(k))).collect();
        let out_text: String =
            report.prune_out.iter().map(|k| format!("{}\n", key(k))).collect();
        files.push(write(prefix, ".prune.in", in_text)?);
        files.push(write(prefix, ".prune.out", out_text)?);
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
