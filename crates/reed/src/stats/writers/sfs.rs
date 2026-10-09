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

use std::path::{Path, PathBuf};

use reed_engine::error::EngineError;

use crate::stats::SfsReport;

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
    let mut files = Vec::new();
    for (name, (spectrum, &n)) in group_names.iter().zip(report.spectra.iter().zip(&report.projections))
    {
        // easySFS keeps the full n+1 length for folded spectra.
        let mut values = spectrum.clone();
        if report.folded {
            values.resize(n as usize + 1, 0.0);
        }
        files.push(write_dadi(prefix, name, n, report.folded, &values, int_dtype)?);
        files.push(write_fsc2(prefix, name, n, report.folded, &values, int_dtype)?);
    }
    Ok(files)
}

fn write_dadi(
    prefix: &Path,
    name: &str,
    n: u32,
    folded: bool,
    values: &[f64],
    int_dtype: bool,
) -> Result<PathBuf, EngineError> {
    let path = prefixed(prefix, &format!("_sfs/dadi/{name}-{n}.sfs"));
    let mut out = String::new();
    out.push_str(&format!("{} {}\n", n + 1, if folded { "folded" } else { "unfolded" }));
    out.push_str(&line(values, int_dtype));
    out.push('\n');
    let mask: Vec<&str> = (0..=n as usize)
        .map(|k| {
            let non_seg = k == 0 || if folded { k > n as usize / 2 } else { k == n as usize };
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
) -> Result<PathBuf, EngineError> {
    let kind = if folded { "MAF" } else { "DAF" };
    let path = prefixed(prefix, &format!("_sfs/fastsimcoal2/{name}_{kind}pop0.obs"));
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
                format!("{}", v.round() as i64)
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

fn write_file(path: &Path, text: String) -> Result<(), EngineError> {
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(EngineError::Io)?;
    }
    std::fs::write(path, text).map_err(EngineError::Io)
}
