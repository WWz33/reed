//! `reed index` dispatch (02-cli.md §1): the scan lives in `reed::index`.

use std::path::Path;

use reed::index::{run_index, WriteIndex};

use crate::cli::IndexCli;

/// Builds the index.
///
/// # Errors
///
/// Unknown index types, transport failures and unsorted input.
pub fn run_index_cli(args: IndexCli) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let IndexCli { input, index_type, out } = args;
    let flavor = match index_type.as_str() {
        "csi" => WriteIndex::Csi,
        "tbi" => WriteIndex::Tbi,
        other => return Err(format!("unknown index type {other:?} (expected csi or tbi)").into()),
    };
    let started = std::time::Instant::now();
    let outcome = run_index(Path::new(&input), flavor, out.as_deref().map(Path::new))?;
    eprintln!(
        "reed index: {} records -> {} in {:.1}s",
        outcome.n_records,
        outcome.path.display(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}
