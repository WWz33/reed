//! `reed rename` flag translation (doc_v14 02 §4): the identifier edits —
//! sample names, contig names, variant IDs and header lines. The flag
//! validation lives in the driver (`run_rename`), which checks every plane
//! before the first record (02 §4.3).

use std::path::PathBuf;

use reed::budget::{MemBudget, ThreadBudget};
use reed::{ChrStyle, RenameOutput, RenameSpec};

use crate::cli::RenameCli;
use crate::units::input_path;

/// Runs `reed rename` for the parsed flags.
///
/// # Errors
///
/// Usage failures from the flag validation and driver failures.
pub(crate) fn run_rename_cli(args: RenameCli) -> Result<(), reed::Error> {
    args.verbosity.init_tracing();
    let started = std::time::Instant::now();
    let spec = RenameSpec {
        header: args.header.as_ref().map(PathBuf::from),
        fai: args.fai.as_ref().map(PathBuf::from),
        add_meta: args.add_meta.clone(),
        remove_meta: args.remove_meta.clone(),
        rename_chrs: args.rename_chrs.as_ref().map(PathBuf::from),
        chr_style: args
            .chr_style
            .as_deref()
            .map(ChrStyle::parse)
            .transpose()
            .map_err(reed::Error::usage)?,
        set_id: args.set_id.clone(),
        set_missing_var_ids: args.set_missing_var_ids.clone(),
        update_name: args.update_name.as_ref().map(PathBuf::from),
        idspace_to: args.idspace_to.clone(),
        rm_dup_ids: args.rm_dup_ids,
        clear_ids: args.clear_ids,
        samples: args.samples.clone(),
        samples_file: args.samples_file.as_ref().map(PathBuf::from),
        replace_sample_name: args.replace_sample_name.clone(),
    };
    // `-O`/`-o`: the same rules as filter's `-O v|z[level]` (convert's
    // hapmap→VCF target shares the spelling).
    let (compressed, level) = match args.output_type.as_deref() {
        None => {
            // The extension decides when `-O` is absent.
            let compressed = args
                .out
                .as_deref()
                .is_some_and(|o| o.ends_with(".gz") || o.ends_with(".bgz"));
            (compressed, 6_u8)
        }
        Some(text) => {
            let (kind, level) = text.split_at(1);
            let level = if level.is_empty() {
                6
            } else {
                level.parse::<u8>().map_err(|_| {
                    reed::Error::usage(format!("invalid -O level in {text:?} (want 1-9)"))
                })?
            };
            match kind {
                "v" => (false, level),
                "z" => (true, level),
                "b" | "u" => {
                    return Err(reed::Error::usage(
                        "reed rename writes VCF text: -O accepts v or z (BCF input/output is \
                         not supported yet)",
                    ));
                }
                other => {
                    return Err(reed::Error::usage(format!(
                        "unknown -O type {other:?} (expected v or z)"
                    )));
                }
            }
        }
    };
    let threads = args.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .expect("invariant: parallelism")
            .get()
    });
    let budget = ThreadBudget {
        threads,
        mem: args
            .memory
            .map_or_else(MemBudget::default_half_memory, |m| {
                MemBudget::from_bytes(m.saturating_mul(1 << 20))
            }),
    };
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let write_index = reed::WriteIndex::parse(args.write_index.as_deref())?;
    let output = RenameOutput {
        path: args.out.as_ref().map(PathBuf::from),
        compressed,
        level,
        command_line,
        no_version: args.no_version,
        write_index: (write_index != reed::WriteIndex::None).then_some(write_index),
    };
    let input_file = input_path(&args.input);
    let report = reed::run_rename(input_file.as_deref(), &spec, &output, &budget)?;
    if let Some(stale) = &report.stale_index {
        eprintln!(
            "reed rename: warning: {} is stale after renaming (the CHROM/contig names \
             changed); run `reed index -f {}`",
            stale.display(),
            report
                .files
                .first()
                .map_or_else(|| "-".to_string(), |p| p.display().to_string())
        );
    }
    if args.verbosity.shows_summary() {
        let secs = started.elapsed().as_secs_f64();
        eprintln!(
            "{} records, {} contig(s), {} sample(s), {} ID(s){} in {:.1}s -> {}",
            report.n_records,
            report.n_contigs,
            report.n_samples,
            report.n_ids,
            if report.n_dup_ids > 0 {
                format!(" ({} deduplicated)", report.n_dup_ids)
            } else {
                String::new()
            },
            secs,
            if report.files.is_empty() {
                "stdout".to_string()
            } else {
                report
                    .files
                    .iter()
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            }
        );
    }
    Ok(())
}
