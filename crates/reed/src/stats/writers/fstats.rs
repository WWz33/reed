//! Patterson f-statistics writer (docs/rounds/v5/05-algorithms.md §3.3):
//! `PREFIX.fstats` with `stat pops estimate se z n_blocks n_snps`, one row
//! per test. `pops` names the populations in the test's role order
//! (`f3` prints `C;A,B`-style roles joined by commas for a flat column).

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use reed_engine::error::EngineError;

use crate::stats::FStatsReport;

/// Writes `PREFIX.fstats`; returns the path written.
///
/// # Errors
///
/// [`EngineError::Io`] for write failures.
pub fn write_fstats(
    prefix: &Path,
    group_names: &[String],
    report: &FStatsReport,
) -> Result<Vec<PathBuf>, EngineError> {
    let mut out = String::from("stat\tpops\testimate\tse\tz\tn_blocks\tn_snps\n");
    for row in &report.rows {
        let pops = row
            .test
            .pops()
            .iter()
            .map(|&g| group_names.get(g).map_or("?", String::as_str))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            out,
            "{}\t{pops}\t{}\t{}\t{}\t{}\t{}",
            row.test.id(),
            cell(row.estimate),
            cell(row.se),
            cell(row.z),
            row.n_blocks,
            row.n_sites,
        )
        .expect("invariant: writing to String");
    }
    let mut path = prefix.as_os_str().to_os_string();
    path.push(".fstats");
    let path = PathBuf::from(path);
    std::fs::write(&path, out).map_err(EngineError::Io)?;
    Ok(vec![path])
}

/// `NA` for non-finite values (fewer than two blocks, zero denominator).
fn cell(v: f64) -> String {
    if v.is_finite() {
        format!("{v}")
    } else {
        "NA".to_string()
    }
}
