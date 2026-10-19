//! Length/`-` argument parsing (02-cli.md §2).

use std::path::PathBuf;

/// Length parser for window/step: plain bp, or k/kb/M/Mb suffixes
/// (02-cli.md §2).
///
/// # Errors
///
/// A non-numeric length or one whose scaled value overflows `u64`.
pub(crate) fn parse_size(text: &str) -> Result<u64, String> {
    let lower = text.trim().to_ascii_lowercase();
    let (digits, scale) = if let Some(rest) = lower.strip_suffix("kb") {
        (rest, 1_000_u64)
    } else if let Some(rest) = lower.strip_suffix('k') {
        (rest, 1_000)
    } else if let Some(rest) = lower.strip_suffix("mb") {
        (rest, 1_000_000)
    } else if let Some(rest) = lower.strip_suffix('m') {
        (rest, 1_000_000)
    } else if let Some(rest) = lower.strip_suffix("bp") {
        (rest, 1)
    } else {
        (lower.as_str(), 1)
    };
    let value: u64 = digits
        .trim()
        .parse()
        .map_err(|_| format!("{text:?} is not a length (use bp, k/kb, M/Mb)"))?;
    value
        .checked_mul(scale)
        .ok_or_else(|| format!("{text:?} is too large"))
}

/// `-` means stdin; everything else is a file path.
#[must_use]
pub(crate) fn input_path(arg: &str) -> Option<PathBuf> {
    (arg != "-").then(|| PathBuf::from(arg))
}

/// `-` means stdout; everything else is a file path (doc_v15 D15-39).
#[must_use]
pub(crate) fn output_path(arg: &str) -> Option<PathBuf> {
    (arg != "-").then(|| PathBuf::from(arg))
}
