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
            let (kind, level) = super::common::split_output_type(
                text,
                "v or z with an optional 1-9 level, e.g. z6",
            )?;
            let level = if level.is_empty() {
                6
            } else {
                let parsed = level.parse::<u8>().map_err(|_| {
                    reed::Error::usage(format!("invalid -O level in {text:?} (want 1-9)"))
                })?;
                // doc_v21 S15: an out-of-range level is a usage error, not a
                // writer-time "invalid compression level" (exit 1).
                if !(1..=9).contains(&parsed) {
                    return Err(reed::Error::usage(format!(
                        "invalid -O level in {text:?} (want 1-9)"
                    )));
                }
                parsed
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
    // doc_v21 S16: `-W` with no possible index target used to be a silent
    // no-op. Checked here, not in the library, so library callers are
    // unaffected.
    if write_index != reed::WriteIndex::None {
        if args.out.is_none() {
            return Err(reed::Error::usage(
                "-W needs -o FILE: an index is written next to the output file, not stdout",
            ));
        }
        if !compressed {
            return Err(reed::Error::usage(
                "-W needs a BGZF output (-Oz): a plain VCF cannot be indexed",
            ));
        }
    }
    let output = RenameOutput {
        path: args.out.as_deref().and_then(crate::units::output_path),
        compressed,
        level,
        command_line,
        no_version: args.no_version,
        no_command: args.no_command,
        write_index: (write_index != reed::WriteIndex::None).then_some(write_index),
    };
    let input_file = input_path(&args.input);
    // (doc_v20 R5) refuse `-o` == input before the writer opens.
    super::common::refuse_output_eq_input(input_file.as_deref(), output.path.as_deref())?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    /// Parses a flag vector through the real clap definitions (the filter
    /// module's `cli_from` shape).
    fn cli_from(args: &[String]) -> RenameCli {
        let mut argv = vec![
            "reed".to_string(),
            "rename".to_string(),
            "in.vcf".to_string(),
        ];
        argv.extend_from_slice(args);
        match crate::cli::Cli::try_parse_from(argv)
            .expect("invariant: flags parse")
            .command
        {
            crate::cli::Command::Rename(cli) => cli,
            _ => unreachable!("invariant: the rename subcommand was given"),
        }
    }

    /// doc_v20 R4: `-O ''` (clap accepts the empty value) is a usage error
    /// naming the accepted spellings — the parser used to `split_at(1)` the
    /// empty string and abort with exit 101.
    #[test]
    fn r4_empty_output_type_is_a_usage_error() {
        let err = run_rename_cli(cli_from(&["-O".to_string(), String::new()]))
            .expect_err("the empty -O must be refused");
        assert_eq!(err.kind(), reed::ErrorKind::Usage, "{err}");
        assert!(err.to_string().contains("invalid -O ''"), "{err}");
    }
}
