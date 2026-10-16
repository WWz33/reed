//! `reed convert` flag translation (docs/rounds/v13/02 §1): the gate rules,
//! the request construction, the run and the output log.

use std::path::{Path, PathBuf};

use reed::budget::{MemBudget, ThreadBudget};
use reed::{
    ConvertRequest, ConvertTarget, MsaRegion, PhyloFormats, PhyloOptions, Populations, run_convert,
    run_convert_opened,
};

use crate::cli::ConvertCli;

/// The default output prefix: the input file name with its known extensions
/// stripped (`in.vcf.gz` → `in`, v13 §3.6). `None` for stdin.
fn default_prefix(path: &Path) -> Option<PathBuf> {
    let name = path.file_name()?.to_str()?;
    let mut stem = name;
    loop {
        let mut changed = false;
        for ext in [".gz", ".bgz", ".vcf", ".bcf"] {
            if let Some(rest) = stem.strip_suffix(ext) {
                stem = rest;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    (!stem.is_empty()).then(|| PathBuf::from(stem))
}

/// Parses one `--msa` region spelling (`CHR`, `CHR:B-E`, `CHR:B-`, `CHR:-E`).
///
/// # Errors
///
/// A malformed region spelling.
fn parse_msa_region(text: &str) -> Result<MsaRegion, reed::Error> {
    let parsed = reed::spec::Region::parse(text)
        .map_err(|e| reed::Error::usage(format!("invalid --region {text:?}: {e}")))?;
    Ok(MsaRegion {
        contig: parsed.contig,
        start1: parsed.start,
        end1: parsed.end,
    })
}

/// Parses one `-O` spelling for `--vcf`: `v` (plain text) or `z[LEVEL]`
/// (BGZF, default level 6). The `b`/`u` (BCF) spellings are rejected.
///
/// # Errors
///
/// An unsupported type letter or an out-of-range level.
fn parse_vcf_output_type(spec: &str) -> Result<(bool, u8), reed::Error> {
    let (kind, rest) = spec.split_at(1);
    match kind {
        "v" if rest.is_empty() => Ok((false, 6)),
        "z" => {
            let level = if rest.is_empty() {
                6
            } else {
                rest.parse::<u8>().map_err(|_| {
                    reed::Error::usage(format!("invalid -O {spec:?}: the level must be 1-9"))
                })?
            };
            if !(1..=9).contains(&level) {
                return Err(reed::Error::usage(format!(
                    "invalid -O {spec:?}: the level must be 1-9"
                )));
            }
            Ok((true, level))
        }
        _ => Err(reed::Error::usage(
            "hapmap→VCF supports -O v or -O z only (no BCF output)",
        )),
    }
}

/// Runs the hapmap→VCF reverse conversion (`reed convert --vcf`, doc_v13
/// §3.12). The input is a hapmap (`.hmp.txt`), so every VCF-input target and
/// its knobs are rejected up front.
fn run_vcf_cli(args: ConvertCli) -> Result<(), reed::Error> {
    // --vcf reads a hapmap: exclusive with every VCF-input target (rule 6).
    if args.phylip
        || args.fasta
        || args.nexus
        || args.nexus_binary
        || args.sfs.is_some()
        || args.count_012
        || args.hapmap
        || args.msa
    {
        return Err(reed::Error::usage(
            "cannot mix --vcf (a hapmap input) with VCF-input targets",
        ));
    }
    // Sample subsetting is not supported: the hapmap header is the sample set.
    if !args.sample_args.samples.is_empty()
        || args.sample_args.samples_file.is_some()
        || args.sample_args.keep.is_some()
        || args.sample_args.remove.is_some()
    {
        return Err(reed::Error::usage(
            "--vcf converts every sample in the hapmap header; \
             -s/-S/--keep/--remove are not supported",
        ));
    }
    // Every other target's knob is inert here; reject rather than ignore.
    if args.proj.is_some()
        || args.sfs_int
        || args.pop.is_some()
        || args.force_samples
        || args.outgroup.is_some()
        || args.min_samples_locus.is_some()
        || args.resolve_iupac
        || args.write_used_sites
        || args.reference.is_some()
        || !args.region.is_empty()
        || args.flank > 0
        || args.collapse_suffix.is_some()
    {
        return Err(reed::Error::usage(
            "--vcf (a hapmap input) takes no VCF-input options; only -o, -O, \
             --no-version and --no-command apply",
        ));
    }

    // -O wins; otherwise infer from the -o extension (.gz/.bgz → BGZF z6).
    let (compressed, level) = match args.output_type.as_deref() {
        Some(spec) => parse_vcf_output_type(spec)?,
        None => {
            let gz = args.out.as_deref().is_some_and(|p| {
                let lower = p.to_ascii_lowercase();
                lower.ends_with(".gz") || lower.ends_with(".bgz")
            });
            (gz, 6)
        }
    };

    let command_line = if args.no_command {
        String::new()
    } else {
        std::env::args().collect::<Vec<_>>().join(" ")
    };
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Vcf {
            compressed,
            level,
            output: args.out.as_deref().map(PathBuf::from),
            command_line,
            no_version: args.no_version,
        }],
        ..ConvertRequest::default()
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

    let input_file = crate::units::input_path(&args.input);
    let run_start = std::time::Instant::now();
    let report = run_convert(input_file.as_deref().map(Path::new), &request, &budget)?;
    let run_secs = run_start.elapsed().as_secs_f64();

    if args.verbosity.shows_summary() {
        let dest = args.out.as_deref().unwrap_or("stdout");
        eprintln!(
            "{} records in {:.1}s -> {}",
            report.n_records,
            report.elapsed.as_secs_f64(),
            dest
        );
    }
    if args.verbosity.shows_timing() {
        eprintln!("# reed convert timing");
        eprintln!("phase.run_ms\t{:.3}", run_secs * 1e3);
        eprintln!("n_records\t{}", report.n_records);
        eprintln!("chunks\t{}", report.metrics.n_chunks_read);
        eprintln!(
            "engine.reader_ms\t{:.3}",
            report.metrics.reader_time.as_secs_f64() * 1e3
        );
        eprintln!(
            "engine.reduce_ms\t{:.3}",
            report.metrics.reduce_time.as_secs_f64() * 1e3
        );
    }
    if let Some(counters) = &report.counters
        && counters.malformed > 0
    {
        eprintln!(
            "warning: skipped {} malformed hapmap row(s)",
            counters.malformed
        );
    }
    Ok(())
}

/// Runs the conversion for the parsed flags.
///
/// # Errors
///
/// Translation (gate rules) and driver failures.
pub(crate) fn run_convert_cli(args: ConvertCli) -> Result<(), reed::Error> {
    args.verbosity.init_tracing();
    let started = std::time::Instant::now();

    // ---- The hapmap→VCF reverse target (doc_v13 §3.12) -------------------
    if args.to_vcf {
        return run_vcf_cli(args);
    }
    if args.output_type.is_some() {
        return Err(reed::Error::usage("-O/--output-type is only used by --vcf"));
    }
    if args.no_version || args.no_command {
        return Err(reed::Error::usage(
            "--no-version/--no-command are only used by --vcf",
        ));
    }

    // ---- Gate rules (v13 §1.1) -------------------------------------------
    let phylo_target = args.phylip || args.fasta || args.nexus;
    let phylo_adjusted = phylo_target || args.nexus_binary;
    if !(phylo_target
        || args.nexus_binary
        || args.sfs.is_some()
        || args.count_012
        || args.hapmap
        || args.msa)
    {
        return Err(reed::Error::usage(
            "nothing to convert: give --sfs, a phylo format (--phylip/--fasta/--nexus), \
             --012, --hapmap, --msa or --vcf",
        ));
    }
    if args.proj.is_some() && args.sfs.is_none() {
        return Err(reed::Error::usage("--proj needs --sfs"));
    }
    if (args.pop.is_some() || args.force_samples) && args.sfs.is_none() {
        return Err(reed::Error::usage(
            "--pop/--force-samples are only used by --sfs",
        ));
    }
    if args.outgroup.is_some() && !phylo_adjusted {
        return Err(reed::Error::usage(
            "--outgroup needs a phylo output (--phylip/--fasta/--nexus/--nexus-binary)",
        ));
    }
    if args.resolve_iupac && !phylo_target {
        return Err(reed::Error::usage(
            "--resolve-iupac needs a nucleotide output (--phylip/--fasta/--nexus)",
        ));
    }
    if args.write_used_sites && !phylo_target {
        return Err(reed::Error::usage(
            "--write-used-sites needs a nucleotide output (--phylip/--fasta/--nexus)",
        ));
    }
    if args.min_samples_locus.is_some() && !(phylo_adjusted || args.count_012 || args.hapmap) {
        return Err(reed::Error::usage(
            "--min-samples-locus needs a locus-stream output \
             (--phylip/--fasta/--nexus/--nexus-binary/--012/--hapmap)",
        ));
    }
    if args.msa && args.reference.is_none() {
        return Err(reed::Error::usage("--msa needs --reference FILE (a FASTA)"));
    }
    if args.reference.is_some() && !args.msa {
        return Err(reed::Error::usage("--reference is only used by --msa"));
    }
    if !args.msa && (!args.region.is_empty() || args.flank > 0 || args.collapse_suffix.is_some()) {
        return Err(reed::Error::usage(
            "--region/--flank/--collapse-suffix are only used by --msa",
        ));
    }

    // ---- Open once (stdin streams) ---------------------------------------
    let input_file = crate::units::input_path(&args.input);
    let prefix = match &args.out {
        Some(p) => PathBuf::from(p),
        None => input_file
            .as_deref()
            .and_then(default_prefix)
            .ok_or_else(|| {
                reed::Error::usage(
                    "conversion writes files: pass -o PREFIX (stdin inputs have no name)",
                )
            })?,
    };
    let mut opened = Some(reed::open_regions(input_file.as_deref(), &[])?);
    let facts = opened.as_ref().expect("invariant: just opened").facts();
    let sample_names: Vec<String> = facts.sample_names.clone();
    let keep_samples = super::stats::parse_sample_keep(&args.sample_args, &sample_names)
        .map_err(reed::Error::usage)?;

    // ---- Targets ---------------------------------------------------------
    let populations = match &args.pop {
        Some(groups) => Populations::File(PathBuf::from(groups)),
        None => Populations::All,
    };
    let group_names = populations.group_names()?;
    let mut targets = Vec::new();
    if let Some(sfs_request) =
        super::build_sfs_request(args.sfs.as_deref(), args.proj.as_deref(), &group_names)?
    {
        targets.push(ConvertTarget::Sfs {
            request: sfs_request,
            int_dtype: args.sfs_int,
        });
    }
    if phylo_target {
        targets.push(ConvertTarget::Phylo {
            formats: PhyloFormats {
                phylip: args.phylip,
                fasta: args.fasta,
                nexus: args.nexus,
            },
            options: PhyloOptions {
                outgroup: args.outgroup.clone(),
                resolve_iupac: args.resolve_iupac,
                write_used_sites: args.write_used_sites,
            },
        });
    }
    if args.nexus_binary {
        targets.push(ConvertTarget::NexusBin);
    }
    if args.count_012 {
        targets.push(ConvertTarget::CountMatrix);
    }
    if args.hapmap {
        targets.push(ConvertTarget::Hapmap);
    }
    if args.msa {
        let mut regions = Vec::with_capacity(args.region.len());
        for text in &args.region {
            regions.push(parse_msa_region(text)?);
        }
        targets.push(ConvertTarget::Msa {
            reference: PathBuf::from(args.reference.as_deref().expect("invariant: checked above")),
            regions,
            flank: args.flank,
            collapse_suffix: args.collapse_suffix.clone(),
        });
    }
    let request = ConvertRequest {
        targets,
        prefix: prefix.clone(),
        populations,
        force_samples: args.force_samples,
        keep_samples,
        min_samples_locus: args.min_samples_locus,
        seed: 0,
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

    let source = opened.take().expect("invariant: convert reads once");
    let run_start = std::time::Instant::now();
    let report = run_convert_opened(source, &request, &budget)?;
    let run_secs = run_start.elapsed().as_secs_f64();

    // ---- The run log (v13 §5.2) ------------------------------------------
    let mut files = report.files.clone();
    {
        use std::fmt::Write as _;
        let mut log = String::new();
        let _ = writeln!(log, "reed {}", env!("CARGO_PKG_VERSION"));
        let _ = writeln!(
            log,
            "command: {}",
            std::env::args().collect::<Vec<_>>().join(" ")
        );
        let _ = writeln!(
            log,
            "targets: {}",
            request.to_args().unwrap_or_default().join(" ")
        );
        if let Some(counters) = &report.counters {
            let _ = writeln!(
                log,
                "counters: total_records={} excess_missing={} mnps={} accepted={} \
                 biallelic_for_bin={} malformed={}",
                counters.total_records,
                counters.excess_missing,
                counters.mnps,
                counters.accepted,
                counters.biallelic_for_bin,
                counters.malformed
            );
        }
        let _ = writeln!(
            log,
            "outgroup: {}",
            args.outgroup.as_deref().unwrap_or("none")
        );
        let _ = writeln!(log, "records: {}", report.n_records);
        log.push_str(
            "methods:\n  phylo  PHYLIP/FASTA/NEXUS per Ortiz 2019 (vcf2phylip semantics)\n  \
             012    vcftools 0/1/2 matrix\n  hapmap TASSEL .hmp.txt\n  \
             msa    reference-driven alignment (vcf2msa semantics, no external aligner)\n  \
             sfs    Marth 2004; Gutenkunst 2009\n",
        );
        let mut path = prefix.clone().into_os_string();
        path.push(".log");
        let path = PathBuf::from(path);
        std::fs::write(&path, log)?;
        files.push(path);
    }

    if args.verbosity.shows_summary() {
        eprintln!(
            "{} records in {:.1}s -> {}",
            report.n_records,
            report.elapsed.as_secs_f64(),
            files
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }
    if args.verbosity.shows_timing() {
        eprintln!("# reed convert timing");
        eprintln!("phase.open_ms\t{:.3}", 0.0);
        eprintln!("phase.run_ms\t{:.3}", run_secs * 1e3);
        eprintln!("n_records\t{}", report.n_records);
        eprintln!("chunks\t{}", report.metrics.n_chunks_read);
        eprintln!(
            "engine.reader_ms\t{:.3}",
            report.metrics.reader_time.as_secs_f64() * 1e3
        );
        eprintln!(
            "engine.reduce_ms\t{:.3}",
            report.metrics.reduce_time.as_secs_f64() * 1e3
        );
        eprintln!(
            "worker_busy_ms\t{:.3}",
            report.metrics.worker_busy.as_secs_f64() * 1e3
        );
    }
    let _ = started;
    Ok(())
}
