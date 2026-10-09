//! `reed index` dispatch (02-cli.md §1): the scan lives in `reed::index`.

use std::path::Path;

use reed::index::{run_index, WriteIndex};

use crate::cli::IndexCli;

/// Builds the index.
///
/// # Errors
///
/// Unknown index types, stdin input (the index needs a seekable file),
/// transport failures and unsorted input.
pub(crate) fn run_index_cli(args: IndexCli) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let IndexCli { input, csi: _, tbi, output, force, threads, verbosity } = args;
    verbosity.init_tracing();
    // The index stores virtual offsets against the file's end and re-reads
    // the input; stdin can do neither (doc_v4 01 F4-4). Fail with the reason
    // instead of a `No such file` from the transport.
    if input == "-" {
        return Err(
            "reed index needs a seekable file: stdin is not supported (the index stores \
             the file's end offset and re-reads it)"
                .into(),
        );
    }
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
    tracing::debug!(flavor = flavor.suffix(), threads, input, "reed index starting");
    let started = std::time::Instant::now();
    let outcome = run_index(Path::new(&input), flavor, out, threads)?;
    if verbosity.shows_summary() {
        eprintln!(
            "reed index: {} records -> {} in {:.1}s",
            outcome.n_records,
            outcome.path.display(),
            started.elapsed().as_secs_f64()
        );
    }
    if verbosity.shows_timing() {
        eprintln!("# reed index timing");
        eprintln!("flavor\t{}", flavor.suffix());
        eprintln!("threads\t{threads}");
        eprintln!("n_records\t{}", outcome.n_records);
        eprintln!("total_ms\t{:.3}", started.elapsed().as_secs_f64() * 1e3);
    }
    Ok(())
}
