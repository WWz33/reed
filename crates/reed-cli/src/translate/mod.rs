//! Flag → typed request translation and subcommand dispatch.

pub(crate) mod common;
pub(crate) mod convert;
pub(crate) mod filter;
pub(crate) mod index;
pub(crate) mod popg;
pub(crate) mod rename;
pub(crate) mod stats;

/// Builds the projected-SFS request from `--sfs [mode]` and `--proj`
/// (docs/rounds/v5/05-algorithms.md §3.1). This is the single
/// implementation: both `reed stats` and `reed convert` call it, so the two
/// entries cannot drift (doc_v13 D13-6).
///
/// # Errors
///
/// `--sfs` without `--proj`, a non-integer projection value, or a projection
/// count that contradicts the population file.
pub(crate) fn build_sfs_request(
    sfs: Option<&str>,
    proj: Option<&str>,
    group_names: &[String],
) -> Result<Option<reed::SfsRequest>, reed::Error> {
    let Some(mode) = sfs else {
        return Ok(None);
    };
    let Some(raw) = proj else {
        return Err(reed::Error::usage(
            "--sfs needs --proj N (the projection size in chromosomes)",
        ));
    };
    let mut values = Vec::new();
    for token in raw.split(',') {
        let token = token.trim();
        values.push(token.parse::<u32>().map_err(|_| {
            reed::Error::usage(format!("invalid --proj value {token:?} (want an integer)"))
        })?);
    }
    if values.len() == 1 {
        values = vec![values[0]; group_names.len()];
    }
    if values.len() != group_names.len() {
        return Err(reed::Error::usage(format!(
            "--proj gives {} value(s) but the population file declares {}",
            values.len(),
            group_names.len()
        )));
    }
    Ok(Some(reed::SfsRequest {
        projections: values,
        folded: mode == "folded",
    }))
}
