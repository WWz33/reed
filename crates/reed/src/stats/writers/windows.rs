//! Window-layout writers (02-cli.md §5.4): `PREFIX_pi.txt`,
//! `PREFIX_dxy.txt`, `PREFIX.log`. The column layout is the one the
//! reference implementation emits — verified column by column in
//! `gold/popgen/compare.py` — so downstream scripts written for that
//! layout work unchanged; missing values print `NA`.

use std::fmt::Write as _;
use std::path::Path;

use reed_stats::request::{StatKey, StatKind};

use crate::stats::{InvariantMode, PopgenReport, WindowLayoutOptions, WindowRow};

/// Renders the flat wide table (the `layout: None` debug output): a header
/// (`contig`, `start`, `end`, then one column per statistic key) and one row
/// per window. Cell rules match the window-layout writers: SFS bins are
/// `;`-joined, a missing value is `.`, and `variant_only` picks the
/// variant-only denominator.
#[must_use]
pub fn render_tsv(rows: &[WindowRow], kinds: &[StatKind], variant_only: bool) -> String {
    let keys: Vec<StatKey> = kinds.iter().map(StatKind::key).collect();
    let mut out = String::from("contig\tstart\tend");
    for key in &keys {
        out.push('\t');
        out.push_str(key.as_str());
    }
    out.push('\n');
    for row in rows {
        out.push_str(&row.contig);
        out.push('\t');
        out.push_str(&row.start1.to_string());
        out.push('\t');
        out.push_str(&row.end.to_string());
        for (kind, acc) in kinds.iter().zip(&row.accs) {
            out.push('\t');
            match acc.sfs_bins() {
                Some(bins) => {
                    for (i, bin) in bins.iter().enumerate() {
                        if i > 0 {
                            out.push(';');
                        }
                        out.push_str(&bin.to_string());
                    }
                }
                None => {
                    let value = if variant_only {
                        acc.value_variant_only(kind, row.bp)
                    } else {
                        acc.value(kind, row.bp)
                    };
                    match value {
                        Some(value) => {
                            write!(out, "{value}").expect("invariant: writing to String");
                        }
                        None => out.push('.'),
                    }
                }
            }
        }
        out.push('\n');
    }
    out
}

/// Writes the requested window-layout files plus the run log.
///
/// # Errors
///
/// Any output-file write failure.
pub fn write(
    prefix: &Path,
    report: &PopgenReport,
    group_names: &[String],
    kinds: &[StatKind],
    layout: WindowLayoutOptions,
    command: &str,
    version: &str,
) -> Result<Vec<std::path::PathBuf>, std::io::Error> {
    let mut pi_rows = String::from(
        "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing\n",
    );
    let mut dxy_rows = String::from(
        "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_dxy\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing\n",
    );
    let mut theta_rows = String::from(
        "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_watterson_theta\tno_sites\traw_watterson_theta\tno_var_sites\tweighted_no_sites\n",
    );
    let mut tajima_rows = String::from(
        "pop\tchromosome\twindow_pos_1\twindow_pos_2\ttajima_d\tno_sites\traw_pi\traw_watterson_theta\ttajima_d_stdev\n",
    );
    // The reference layout prints `avg_wc_fst` / `avg_hudson_fst`
    // (02-cli.md §5.4 spells the WC column avg_wc84_fst; the layout's own
    // name wins so downstream scripts work unchanged).
    let fst_method = kinds.iter().find_map(|kind| match kind {
        StatKind::Fst { kind, .. } => Some(match kind {
            reed_stats::fst::FstKind::Hudson => "hudson",
            reed_stats::fst::FstKind::WeirCockerham84 => "wc",
        }),
        _ => None,
    });
    // One buffer (and later one file) per PBS triple (doc_v15 D15-8): a
    // single triple keeps the documented single-file name, while multiple
    // triples each get a header naming their own groups — a shared header
    // used to name only the first triple (`--pbs A,B,C --pbs C,D,A`).
    let pbs_slots: Vec<Option<usize>> = {
        let mut slots = vec![None; kinds.len()];
        let mut next = 0;
        for (i, kind) in kinds.iter().enumerate() {
            if matches!(kind, StatKind::Pbs { .. }) {
                slots[i] = Some(next);
                next += 1;
            }
        }
        slots
    };
    let mut pbs_buffers: Vec<String> = Vec::new();
    let mut pbs_stems: Vec<String> = Vec::new();
    for kind in kinds {
        let StatKind::Pbs { a, b, c, .. } = kind else {
            continue;
        };
        let name = |i: usize| group_names.get(i).map_or("?", String::as_str);
        let (na, nb, nc) = (name(*a), name(*b), name(*c));
        pbs_buffers.push(format!(
            "{na}\t{nb}\t{nc}\tchromosome\twindow_pos_1\twindow_pos_2\tpbs\tno_snps\tfst_12\tfst_13\tfst_23\n"
        ));
        let base = pbs_stem(na, nb, nc);
        let mut stem = base.clone();
        let mut n = 2;
        while pbs_stems.contains(&stem) {
            stem = format!("{base}-{n}");
            n += 1;
        }
        pbs_stems.push(stem);
    }
    let mut fst_rows = match fst_method {
        Some(method) => format!(
            "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_{method}_fst\tno_snps\n"
        ),
        None => String::new(),
    };
    // Global (r-population) WC84 Fst: one row per window, `pops` = the
    // comma-joined population names (docs/rounds/v5/05 §3.2).
    let pops_cell = group_names.join(",");
    let mut fst_global_rows = if kinds.iter().any(|k| matches!(k, StatKind::FstGlobal)) {
        String::from("pops\tchromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst\tno_snps\n")
    } else {
        String::new()
    };
    let variant_only = layout.invariant == InvariantMode::VariantOnly;
    for row in &report.windows {
        for (kind_idx, (kind, acc)) in kinds.iter().zip(&row.accs).enumerate() {
            match kind {
                StatKind::Pi { .. } | StatKind::Dxy { .. } => {
                    let (Some((diffs, comparable, missing)), Some(n_sites)) =
                        (acc.pair_counts(), acc.n_present())
                    else {
                        continue;
                    };
                    // Empty windows (no site of the file falls in them) print
                    // NA, like the reference layout (04-algorithms.md §4).
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
                                row.contig, row.start1, row.end
                            );
                        }
                        StatKind::Dxy { a, b, .. } => {
                            let name_a = group_names.get(*a).map_or("?", String::as_str);
                            let name_b = group_names.get(*b).map_or("?", String::as_str);
                            let _ = writeln!(
                                dxy_rows,
                                "{name_a}\t{name_b}\t{}\t{}\t{}\t{value}\t{n_sites}\t{diffs}\t{comparable}\t{missing}",
                                row.contig, row.start1, row.end
                            );
                        }
                        _ => unreachable!("guarded above"),
                    }
                }
                StatKind::ThetaW { group, .. } => {
                    let Some((raw, no_sites, no_var, weighted)) = acc.theta_columns() else {
                        continue;
                    };
                    let name = group_names.get(*group).map_or("?", String::as_str);
                    let avg = match value_of(acc, kind, row.bp, variant_only) {
                        Some(v) => format!("{v}"),
                        None => "NA".to_string(),
                    };
                    let _ = writeln!(
                        theta_rows,
                        "{name}\t{}\t{}\t{}\t{avg}\t{no_sites}\t{raw}\t{no_var}\t{weighted}",
                        row.contig, row.start1, row.end
                    );
                }
                StatKind::TajimaD { group, .. } => {
                    let Some((raw_pi, raw_theta, stdev, no_sites)) = acc.tajima_columns() else {
                        continue;
                    };
                    let name = group_names.get(*group).map_or("?", String::as_str);
                    let cell = |v: Option<f64>| match v {
                        Some(v) => format!("{v}"),
                        None => "NA".to_string(),
                    };
                    let d = match acc.value(kind, row.bp) {
                        Some(v) => format!("{v}"),
                        None => "NA".to_string(),
                    };
                    let _ = writeln!(
                        tajima_rows,
                        "{name}\t{}\t{}\t{}\t{d}\t{no_sites}\t{}\t{}\t{}",
                        row.contig,
                        row.start1,
                        row.end,
                        cell(raw_pi),
                        cell(raw_theta),
                        cell(stdev)
                    );
                }
                StatKind::Fst { a, b, .. } => {
                    let Some((value, no_snps)) = acc.fst_columns(kind) else {
                        continue;
                    };
                    let name_a = group_names.get(*a).map_or("?", String::as_str);
                    let name_b = group_names.get(*b).map_or("?", String::as_str);
                    let cell = match value {
                        Some(v) => format!("{v}"),
                        None => "NA".to_string(),
                    };
                    let _ = writeln!(
                        fst_rows,
                        "{name_a}\t{name_b}\t{}\t{}\t{}\t{cell}\t{no_snps}",
                        row.contig, row.start1, row.end
                    );
                }
                StatKind::FstGlobal => {
                    let Some((value, no_snps)) = acc.fst_global_columns(kind) else {
                        continue;
                    };
                    let cell = match value {
                        Some(v) => format!("{v}"),
                        None => "NA".to_string(),
                    };
                    let _ = writeln!(
                        fst_global_rows,
                        "{pops_cell}\t{}\t{}\t{}\t{cell}\t{no_snps}",
                        row.contig, row.start1, row.end
                    );
                }
                StatKind::Pbs { a, b, c, .. } => {
                    let Some((value, no_snps, (f_ab, f_ac, f_bc))) = acc.pbs_columns(kind) else {
                        continue;
                    };
                    let cell = |v: Option<f64>| match v {
                        Some(v) => format!("{v}"),
                        None => "NA".to_string(),
                    };
                    let name = |i: usize| group_names.get(i).map_or("?", String::as_str);
                    let slot = pbs_slots[kind_idx].expect("invariant: PBS kinds were indexed");
                    let _ = writeln!(
                        pbs_buffers[slot],
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{no_snps}\t{}\t{}\t{}",
                        name(*a),
                        name(*b),
                        name(*c),
                        row.contig,
                        row.start1,
                        row.end,
                        cell(value),
                        cell(f_ab),
                        cell(f_ac),
                        cell(f_bc)
                    );
                }
                // The window layout carries pi/dxy plus the extended columns
                // handled above; the SFS accumulator has no window row.
                StatKind::Sfs { .. } => {}
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
    if kinds.iter().any(|k| matches!(k, StatKind::ThetaW { .. })) {
        let path = with_suffix(prefix, "_watterson_theta.txt");
        std::fs::write(&path, theta_rows)?;
        written.push(path);
    }
    if pbs_buffers.len() == 1 {
        let path = with_suffix(prefix, "_pbs.txt");
        std::fs::write(&path, &pbs_buffers[0])?;
        written.push(path);
    } else if !pbs_buffers.is_empty() {
        for (stem, text) in pbs_stems.iter().zip(&pbs_buffers) {
            let path = with_suffix(prefix, &format!("_pbs.{stem}.txt"));
            std::fs::write(&path, text)?;
            written.push(path);
        }
    }
    if fst_method.is_some() {
        let path = with_suffix(prefix, "_fst.txt");
        std::fs::write(&path, fst_rows)?;
        written.push(path);
    }
    if kinds.iter().any(|k| matches!(k, StatKind::TajimaD { .. })) {
        let path = with_suffix(prefix, "_tajima_d.txt");
        std::fs::write(&path, tajima_rows)?;
        written.push(path);
    }
    if kinds.iter().any(|k| matches!(k, StatKind::FstGlobal)) {
        let path = with_suffix(prefix, "_fst_global.txt");
        std::fs::write(&path, fst_global_rows)?;
        written.push(path);
    }
    let log = format!(
        "reed {version}\ncommand: {command}\ninvariant-sites: {:?}\nwindows: {}\nrecords: {}\n\
         invariant sites seen: {}\nmethods:\n  pi   Nei & Li 1979; missing data: Korunes & Samuk 2021\n  \
         dxy  Nei 1987; missing data: Korunes & Samuk 2021\n",
        layout.invariant, report.n_windows, report.n_records, report.n_invariant_sites
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

/// A file-name stem for one PBS triple (`A-B-C`, sanitized): characters
/// outside `[A-Za-z0-9_.-]` become `_` (the naming is documented in
/// cli.md §7; doc_v15 D15-8).
fn pbs_stem(a: &str, b: &str, c: &str) -> String {
    let mut stem = String::with_capacity(a.len() + b.len() + c.len() + 2);
    for (i, name) in [a, b, c].iter().enumerate() {
        if i > 0 {
            stem.push('-');
        }
        for ch in name.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '.' | '-') {
                stem.push(ch);
            } else {
                stem.push('_');
            }
        }
    }
    stem
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
