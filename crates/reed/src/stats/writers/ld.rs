//! LD writers (docs/rounds/v5/05-algorithms.md §3.4): `PREFIX.ld.decay`
//! (`pop bin_start bin_end n_pairs mean_r2` — the cohort under `all`, then
//! one block per population when `--pop` names more than one) and the
//! plink2-style per-pair `PREFIX.ld` (`CHROM_A POS_A ID_A CHROM_B POS_B
//! ID_B UNPHASED_R2`). Pruning moved to `reed filter --indep-pairwise`
//! (docs/rounds/v5/05 §3.4.1).

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use reed_engine::error::EngineError;

use crate::stats::LdReport;

/// The per-pair table's header line. One definition for both writers — the
/// in-memory render below and the streaming sink
/// (`stats::sink::write_pairs`) — so a `--ld-pairs` run is byte-identical
/// whether its rows were accumulated or streamed (doc_v23 V23-14).
pub const PAIRS_HEADER: &[u8] = b"#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n";

/// Appends one per-pair row in the table's format. The streaming writer and
/// the in-memory renderer both call this, so the bytes cannot drift.
pub fn render_pair_row(out: &mut String, contig: &str, a: u64, b: u64, r2: f64) {
    writeln!(out, "{contig}\t{a}\t.\t{contig}\t{b}\t.\t{r2}")
        .expect("invariant: writing to String");
}

/// Writes the LD files the report carries; returns the paths written.
///
/// # Errors
///
/// [`EngineError::Io`] for write failures.
pub fn write_ld(
    prefix: &Path,
    report: &LdReport,
    group_names: &[String],
    pairs: bool,
) -> Result<Vec<PathBuf>, EngineError> {
    let mut files = Vec::new();
    if !report.decay.is_empty() {
        let mut out = String::from("pop\tbin_start\tbin_end\tn_pairs\tmean_r2\n");
        let mut row = |pop: &str, bins: &[(u64, u64, u64, Option<f64>)]| {
            for (start, end, n, mean) in bins {
                let cell = match mean {
                    Some(v) => format!("{v}"),
                    None => "NA".to_string(),
                };
                writeln!(out, "{pop}\t{start}\t{end}\t{n}\t{cell}")
                    .expect("invariant: writing to String");
            }
        };
        row("all", &report.decay);
        // With more than one population the per-population curves follow
        // (a single population's curve is the cohort's).
        if report.group_decay.len() > 1 {
            for (g, bins) in report.group_decay.iter().enumerate() {
                let pop = group_names.get(g).map_or("group", String::as_str);
                row(pop, bins);
            }
        }
        files.push(write(prefix, ".ld.decay", out)?);
    }
    // `pairs` is a request, not a nonempty-content check. An accumulated
    // zero-pair run writes a header-only table; a decay-only request writes
    // no `.ld` at all. The CLI's streaming form owns this path itself
    // (doc_v27 V27-05).
    if !pairs {
        return Ok(files);
    }
    if report.pairs.is_empty() {
        files.push(write(
            prefix,
            ".ld",
            String::from_utf8_lossy(PAIRS_HEADER).into_owned(),
        )?);
    } else {
        let mut out = String::from_utf8_lossy(PAIRS_HEADER).into_owned();
        for (contig, a, b, r2) in &report.pairs {
            render_pair_row(&mut out, contig, *a, *b, *r2);
        }
        files.push(write(prefix, ".ld", out)?);
    }
    Ok(files)
}

fn write(prefix: &Path, suffix: &str, text: String) -> Result<PathBuf, EngineError> {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(suffix);
    let path = PathBuf::from(path);
    reed_io::write::atomic::write_atomic(&path, text).map_err(EngineError::Io)?;
    Ok(path)
}
