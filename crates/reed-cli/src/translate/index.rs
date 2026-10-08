//! `reed index` dispatch (02-cli.md §1): the scan lives in `reed::index`.

use std::path::Path;


use reed::index::{run_index, WriteIndex};

use crate::cli::IndexCli;

/// Builds the index.
///
/// # Errors
///
/// Unknown index types, transport failures and unsorted input.
pub(crate) fn run_index_cli(args: IndexCli) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let IndexCli { input, csi: _, tbi, output, force, threads } = args;
    let flavor = if tbi { WriteIndex::Tbi } else { WriteIndex::Csi };
    let out = output.as_deref().map(Path::new);
    let index_path = match out {
        Some(p) => p.to_path_buf(),
        None => WriteIndex::path_for(flavor, Path::new(&input)),
    };
    if index_path.exists() && !force {
        return Err(format!(
            "{}: index exists (use -f/--force to overwrite)",
            index_path.display()
        )
        .into());
    }
    let threads = threads.unwrap_or_else(|| {
        std::thread::available_parallelism().expect("invariant: parallelism").get()
    });
    let started = std::time::Instant::now();
    let outcome = run_index(Path::new(&input), flavor, out, threads)?;
    eprintln!(
        "reed index: {} records -> {} in {:.1}s",
        outcome.n_records,
        outcome.path.display(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}
