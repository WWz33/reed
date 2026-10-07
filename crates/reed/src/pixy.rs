//! pixy-layout window writers (02-cli.md §5.4): `PREFIX_pi.txt`,
//! `PREFIX_dxy.txt`, `PREFIX.log`. Columns follow pixy so downstream pixy
//! plotting scripts work unchanged; missing values print `NA`.

use std::fmt::Write as _;
use std::path::Path;

use reed_stats::request::StatKind;

use crate::stats::{InvariantMode, PixyOptions, StatsReport};

/// Writes the requested pixy-layout files plus the run log.
///
/// # Errors
///
/// Any output-file write failure.
pub fn write(
    prefix: &Path,
    report: &StatsReport,
    group_names: &[String],
    kinds: &[StatKind],
    pixy: PixyOptions,
    command: &str,
    version: &str,
) -> Result<Vec<std::path::PathBuf>, std::io::Error> {
    let mut pi_rows = String::from(
        "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing\n",
    );
    let mut dxy_rows = String::from(
        "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_dxy\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing\n",
    );
    let variant_only = pixy.invariant == InvariantMode::VariantOnly;
    for row in &report.windows {
        for (kind, acc) in kinds.iter().zip(&row.accs) {
            let (Some((diffs, comparable, missing)), Some(n_sites)) =
                (acc.pair_counts(), acc.n_present())
            else {
                continue;
            };
            // Empty windows (no site of the file falls in them) print NA,
            // like pixy (04-algorithms.md §4).
            let value = if n_sites == 0 {
                "NA".to_string()
            } else {
                match value_of(acc, kind, row.bp, variant_only) {
                    Some(v) => format!("{v}"),
                    None => "NA".to_string(),
                }
            };
            match kind {
                StatKind::Pi { group, .. } => {
                    let name = group_names.get(*group).map_or("?", String::as_str);
                    let _ = writeln!(
                        pi_rows,
                        "{name}\t{}\t{}\t{}\t{value}\t{n_sites}\t{diffs}\t{comparable}\t{missing}",
                        row.contig, row.start1, row.end0
                    );
                }
                StatKind::Dxy { a, b, .. } => {
                    let name_a = group_names.get(*a).map_or("?", String::as_str);
                    let name_b = group_names.get(*b).map_or("?", String::as_str);
                    let _ = writeln!(
                        dxy_rows,
                        "{name_a}\t{name_b}\t{}\t{}\t{}\t{value}\t{n_sites}\t{diffs}\t{comparable}\t{missing}",
                        row.contig, row.start1, row.end0
                    );
                }
                _ => {}
            }
        }
    }
    let mut written = Vec::new();
    if kinds.iter().any(|k| matches!(k, StatKind::Pi { .. })) {
        let path = with_suffix(prefix, "_pi.txt");
        std::fs::write(&path, pi_rows)?;
        written.push(path);
    }
    if kinds.iter().any(|k| matches!(k, StatKind::Dxy { .. })) {
        let path = with_suffix(prefix, "_dxy.txt");
        std::fs::write(&path, dxy_rows)?;
        written.push(path);
    }
    let log = format!(
        "reed {version}\ncommand: {command}\ninvariant-sites: {:?}\nwindows: {}\nrecords: {}\n\
         invariant sites seen: {}\nmethods:\n  pi   Nei & Li 1979; missing data: Korunes & Samuk 2021 (pixy)\n  \
         dxy  Nei 1987; missing data: Korunes & Samuk 2021 (pixy)\n",
        pixy.invariant, report.n_windows, report.n_records, report.n_invariant_sites
    );
    let path = with_suffix(prefix, ".log");
    std::fs::write(&path, log)?;
    written.push(path);
    Ok(written)
}

/// `PREFIX` + suffix appended to the full path (no extension juggling).
fn with_suffix(prefix: &Path, suffix: &str) -> std::path::PathBuf {
    let mut s = prefix.as_os_str().to_os_string();
    s.push(suffix);
    std::path::PathBuf::from(s)
}

/// The window value for `kind` under the declared invariant policy.
fn value_of(
    acc: &reed_stats::request::StatAcc,
    kind: &StatKind,
    bp: u64,
    variant_only: bool,
) -> Option<f64> {
    if variant_only {
        acc.value_variant_only(kind, bp)
    } else {
        acc.value(kind, bp)
    }
}
