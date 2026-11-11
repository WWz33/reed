//! Projected-SFS writers (docs/rounds/v5/05-algorithms.md §3.1): the dadi
//! and fastsimcoal2 layouts, matching easySFS's actual output.
//!
//! dadi: `<n+1> <unfolded|folded>` / the spectrum / a mask line. A folded
//! spectrum keeps the `n+1` length with the upper half zeroed, and the mask
//! marks k = 0 and the structurally-zero top (easySFS convention).
//! fastsimcoal2: `1 observation` / `d0_0 … d0_n` / the spectrum, under
//! `<pop>_DAFpop0.obs` (unfolded) or `<pop>_MAFpop0.obs` (folded).
//!
//! Values print at full round-trip precision; `--sfs-int` rounds them.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use reed_engine::error::EngineError;

use crate::stats::{JointSfsReport, SfsReport};

/// Turns a group name into a filesystem-safe stem: every byte outside
/// `[A-Za-z0-9_.-]` becomes `_` (the PBS convention, doc_v15 D15-8), and the
/// `.`/`..` names become `_`/`__` so a group cannot escape its directory or
/// address a parent.
///
/// Group names come from the operator's populations file and were appended
/// verbatim, so `../../ESCAPED` wrote outside `-o` and `A/B` created nested
/// directories (doc_v27 V27-08); `convert --sfs` shares this writer.
fn stem(name: &str) -> String {
    if name.is_empty() {
        return "_".to_string();
    }
    // A leading dot would make the file hidden (and `.`/`..` would name a
    // directory), so dots at either end — and any `..` run anywhere — map to
    // `_`. Interior single dots (`.`, as in `pop.a`) stay, matching the PBS
    // stem's whitelist.
    let mapped: String = name
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '.' | '-') {
                ch
            } else {
                '_'
            }
        })
        .collect();
    let mut out = String::with_capacity(mapped.len());
    for ch in mapped.chars() {
        if ch == '.' && (out.is_empty() || out.ends_with('.')) {
            out.push('_');
        } else {
            out.push(ch);
        }
    }
    while out.len() > 1 && out.ends_with('.') {
        out.pop();
        out.push('_');
    }
    if out.is_empty() {
        out.push('_');
    }
    out
}

/// Allocates the stems for a whole run: sanitized, de-duplicated in input
/// order with `-2`, `-3`, … suffixes. Two groups that sanitize to the same
/// stem (or just share a name) would otherwise share a file and one
/// spectrum would silently overwrite the other (doc_v27 V27-09).
fn unique_stems(names: &[String]) -> Vec<String> {
    let mut out: Vec<String> = Vec::with_capacity(names.len());
    for name in names {
        let base = stem(name);
        let mut candidate = base.clone();
        let mut n = 2_u32;
        while out.contains(&candidate) {
            candidate = format!("{base}-{n}");
            n += 1;
        }
        out.push(candidate);
    }
    out
}

/// Claims `path` for this run, refusing a second writer the same file:
/// one-dimensional `<name>-<n>.sfs` and pairwise `<A>-<B>.sfs` live in the
/// same `_sfs/dadi` namespace, so group `A` projected to 4 and the pair
/// `A`+`4` used to address the same path and one spectrum silently
/// overwrote the other (doc_v27 V27-09). Refusing beats picking a winner.
///
/// # Errors
///
/// [`EngineError::Job`] on a duplicate path.
fn claim(path: PathBuf, seen: &mut Vec<PathBuf>) -> Result<PathBuf, EngineError> {
    if seen.contains(&path) {
        return Err(EngineError::job(format!(
            "two SFS outputs would write {}; rename the populations (or change --proj) so the \
             file names differ",
            path.display()
        )));
    }
    seen.push(path.clone());
    Ok(path)
}

/// Writes the dadi and fastsimcoal2 files for every projected spectrum;
/// returns the paths written.
///
/// # Errors
///
/// [`EngineError::Job`] when the group names and spectra disagree,
/// [`EngineError::Io`] for write failures.
pub fn write_projected_sfs(
    prefix: &Path,
    group_names: &[String],
    report: &SfsReport,
    int_dtype: bool,
) -> Result<Vec<PathBuf>, EngineError> {
    if group_names.len() != report.spectra.len() {
        return Err(EngineError::job(format!(
            "{} group names but {} spectra",
            group_names.len(),
            report.spectra.len()
        )));
    }
    let stems = unique_stems(group_names);
    let mut files = Vec::new();
    let mut seen: Vec<PathBuf> = Vec::new();
    for (stem, (spectrum, &n)) in stems
        .iter()
        .zip(report.spectra.iter().zip(&report.projections))
    {
        // easySFS keeps the full n+1 length for folded spectra.
        let mut values = spectrum.clone();
        if report.folded {
            values.resize(n as usize + 1, 0.0);
        }
        files.push(write_dadi(
            prefix,
            stem,
            n,
            report.folded,
            &values,
            int_dtype,
            &mut seen,
        )?);
        files.push(write_fsc2(
            prefix,
            stem,
            n,
            report.folded,
            &values,
            int_dtype,
            &mut seen,
        )?);
    }
    if let Some(joint) = &report.joint {
        files.extend(write_joint_sfs(
            prefix,
            group_names,
            joint,
            int_dtype,
            &mut seen,
        )?);
    }
    Ok(files)
}

/// Writes the pairwise joint (2D) spectra in the dadi and fastsimcoal2
/// layouts, matching easySFS's `dadi_twoD_sfs_combinations`: the dadi file
/// is `<A>-<B>.sfs` (no projection suffix) and the fsc2 file is
/// `<prefix>_joint{M,D}AFpop{idxB}_{idxA}.obs` with group indices from the
/// population file (0-based).
///
/// # Errors
///
/// [`EngineError::Job`] when a pair names a group without a name,
/// [`EngineError::Io`] for write failures.
pub fn write_joint_sfs(
    prefix: &Path,
    group_names: &[String],
    joint: &JointSfsReport,
    int_dtype: bool,
    seen: &mut Vec<PathBuf>,
) -> Result<Vec<PathBuf>, EngineError> {
    let stems = unique_stems(group_names);
    let mut files = Vec::with_capacity(joint.pairs.len() * 2);
    for (idx, (&(g1, g2), (values, &(n1, n2)))) in joint
        .pairs
        .iter()
        .zip(joint.spectra.iter().zip(&joint.projections))
        .enumerate()
    {
        let (Some(name1), Some(name2)) = (stems.get(g1), stems.get(g2)) else {
            return Err(EngineError::job(format!(
                "joint spectrum {idx} names groups ({g1}, {g2}) but only {} groups have names",
                group_names.len()
            )));
        };
        files.push(write_joint_dadi(
            prefix,
            name1,
            name2,
            n1,
            n2,
            joint.folded,
            values,
            int_dtype,
            seen,
        )?);
        files.push(write_joint_fsc2(
            prefix,
            g1,
            g2,
            n1,
            n2,
            joint.folded,
            values,
            int_dtype,
            seen,
        )?);
    }
    Ok(files)
}

#[allow(clippy::too_many_arguments)]
fn write_joint_dadi(
    prefix: &Path,
    name1: &str,
    name2: &str,
    n1: u32,
    n2: u32,
    folded: bool,
    values: &[f64],
    int_dtype: bool,
    seen: &mut Vec<PathBuf>,
) -> Result<PathBuf, EngineError> {
    let path = claim(joint_dadi_path(prefix, name1, name2), seen)?;
    let mut out = String::new();
    writeln!(
        out,
        "{} {} {}",
        n1 + 1,
        n2 + 1,
        if folded { "folded" } else { "unfolded" }
    )
    .expect("invariant: writing to String");
    out.push_str(&line(values, int_dtype));
    out.push('\n');
    // Corners are always masked; folded spectra additionally mask the
    // mirror region `i + j > (n1 + n2) / 2` (dadi `Spectrum.fold`).
    let mask: Vec<&str> = (0..=n1)
        .flat_map(|i| {
            (0..=n2).map(move |j| {
                let masked = if folded {
                    (i == 0 && j == 0) || i + j > u32::midpoint(n1, n2)
                } else {
                    (i == 0 && j == 0) || (i == n1 && j == n2)
                };
                if masked { "1" } else { "0" }
            })
        })
        .collect();
    out.push_str(&mask.join(" "));
    out.push('\n');
    write_file(&path, out)?;
    Ok(path)
}

/// The fastsimcoal2 joint layout: columns = group `g1`'s bins (its
/// projection `n1`), rows = group `g2`'s bins, both labelled by the
/// population-file order (`d{idx}_{k}`), matching easySFS's transposed
/// rendering of the dadi matrix.
#[allow(clippy::too_many_arguments)]
fn write_joint_fsc2(
    prefix: &Path,
    g1: usize,
    g2: usize,
    n1: u32,
    n2: u32,
    folded: bool,
    values: &[f64],
    int_dtype: bool,
    seen: &mut Vec<PathBuf>,
) -> Result<PathBuf, EngineError> {
    // easySFS names the joint fsc2 file `<prefix>_joint{M,D}AFpop{idx2}_{idx1}.obs`
    // where its prefix is the input basename; reed's analog is the `-o` prefix.
    let path = claim(joint_fsc2_path(prefix, g1, g2, folded), seen)?;
    let mut out = String::new();
    out.push_str("1 observation\n");
    let header: Vec<String> = (0..=n1).map(|k| format!("d{g1}_{k}")).collect();
    out.push('\t');
    out.push_str(&header.join("\t"));
    out.push('\n');
    // Row-major source `i * (n2 + 1) + j` (i over group 1) is transposed:
    // row = group 2's bin, columns = group 1's bins; values space-separated
    // (easySFS writes `row_head + "\t" + " ".join(values)`).
    for j in 0..=n2 {
        let _ = write!(out, "d{g2}_{j}\t");
        for (k, i) in (0..=n1).enumerate() {
            if k > 0 {
                out.push(' ');
            }
            let v = values[i as usize * (n2 as usize + 1) + j as usize];
            if int_dtype {
                #[allow(clippy::cast_possible_truncation)]
                let n = v.round() as i64;
                let _ = write!(out, "{n}");
            } else {
                let _ = write!(out, "{v}");
            }
        }
        out.push('\n');
    }
    write_file(&path, out)?;
    Ok(path)
}

fn write_dadi(
    prefix: &Path,
    name: &str,
    n: u32,
    folded: bool,
    values: &[f64],
    int_dtype: bool,
    seen: &mut Vec<PathBuf>,
) -> Result<PathBuf, EngineError> {
    let path = claim(dadi_path(prefix, name, n), seen)?;
    let mut out = String::new();
    writeln!(
        out,
        "{} {}",
        n + 1,
        if folded { "folded" } else { "unfolded" }
    )
    .expect("invariant: writing to String");
    out.push_str(&line(values, int_dtype));
    out.push('\n');
    let mask: Vec<&str> = (0..=n as usize)
        .map(|k| {
            let non_seg = k == 0
                || if folded {
                    k > n as usize / 2
                } else {
                    k == n as usize
                };
            if non_seg { "1" } else { "0" }
        })
        .collect();
    out.push_str(&mask.join(" "));
    out.push('\n');
    write_file(&path, out)?;
    Ok(path)
}

fn write_fsc2(
    prefix: &Path,
    name: &str,
    n: u32,
    folded: bool,
    values: &[f64],
    int_dtype: bool,
    seen: &mut Vec<PathBuf>,
) -> Result<PathBuf, EngineError> {
    let path = claim(fsc2_path(prefix, name, folded), seen)?;
    let mut out = String::new();
    out.push_str("1 observation\n");
    let header: Vec<String> = (0..=n).map(|k| format!("d0_{k}")).collect();
    out.push_str(&header.join("\t"));
    out.push('\n');
    out.push_str(&line(values, int_dtype));
    out.push('\n');
    write_file(&path, out)?;
    Ok(path)
}

/// One value line: space-separated, rounded when `int_dtype`.
fn line(values: &[f64], int_dtype: bool) -> String {
    values
        .iter()
        .map(|v| {
            if int_dtype {
                // Projected spectrum values are non-negative counts bounded
                // by the projection sample size; the cast is exact.
                #[allow(clippy::cast_possible_truncation)]
                let n = v.round() as i64;
                format!("{n}")
            } else {
                format!("{v}")
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn prefixed(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(suffix);
    PathBuf::from(path)
}

/// The one-dimensional dadi path (`<stem>-<n>.sfs`).
fn dadi_path(prefix: &Path, name: &str, n: u32) -> PathBuf {
    prefixed(prefix, &format!("_sfs/dadi/{name}-{n}.sfs"))
}

/// The one-dimensional fastsimcoal2 path (`<stem>_{D,M}AFpop0.obs`).
fn fsc2_path(prefix: &Path, name: &str, folded: bool) -> PathBuf {
    let kind = if folded { "MAF" } else { "DAF" };
    prefixed(prefix, &format!("_sfs/fastsimcoal2/{name}_{kind}pop0.obs"))
}

/// The joint dadi path (`<name1>-<name2>.sfs`).
fn joint_dadi_path(prefix: &Path, name1: &str, name2: &str) -> PathBuf {
    prefixed(prefix, &format!("_sfs/dadi/{name1}-{name2}.sfs"))
}

/// The joint fastsimcoal2 path, named from the `-o` prefix's basename
/// (`<base>_joint{D,M}AFpop<idx2>_<idx1>.obs`).
fn joint_fsc2_path(prefix: &Path, g1: usize, g2: usize, folded: bool) -> PathBuf {
    let kind = if folded { "MAF" } else { "DAF" };
    let base = prefix
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("reed");
    prefixed(
        prefix,
        &format!("_sfs/fastsimcoal2/{base}_joint{kind}pop{g2}_{g1}.obs"),
    )
}

/// Every path [`write_projected_sfs`] would write for `plan`, in writer
/// order. The convert driver validates this set before any writer opens
/// (v29 V29-03); the helpers above are shared with the writers so planning
/// and writing cannot drift.
pub(crate) fn planned_paths(
    prefix: &Path,
    group_names: &[String],
    plan: &crate::stats::request::SfsPlan,
) -> Vec<PathBuf> {
    let stems = unique_stems(group_names);
    let mut paths = Vec::with_capacity(2 * plan.projections.len() + 2 * plan.joint_pairs.len());
    for (stem, &n) in stems.iter().zip(plan.projections.iter()) {
        paths.push(dadi_path(prefix, stem, n));
        paths.push(fsc2_path(prefix, stem, plan.folded));
    }
    for &(g1, g2) in plan.joint_pairs.iter() {
        if let (Some(name1), Some(name2)) = (stems.get(g1), stems.get(g2)) {
            paths.push(joint_dadi_path(prefix, name1, name2));
            paths.push(joint_fsc2_path(prefix, g1, g2, plan.folded));
        }
    }
    paths
}

/// Writes one SFS file, creating its `_sfs/{dadi,fastsimcoal2}` parent.
///
/// This is the only output family whose parent directories are part of the
/// layout rather than the operator's `-o`, so it is the one caller that asks
/// [`reed_io::write::atomic::AtomicWriter::create_with_parents`] to create
/// them (doc_v26 V26-08, doc_v28 V28-23).
fn write_file(path: &Path, text: String) -> Result<(), EngineError> {
    use std::io::Write as _;
    let mut writer = reed_io::write::atomic::AtomicWriter::create_with_parents(
        path,
        reed_io::write::atomic::DEFAULT_BUFFER,
    )
    .map_err(EngineError::Io)?;
    writer.write_all(text.as_bytes()).map_err(EngineError::Io)?;
    writer.commit().map_err(EngineError::Io)
}
