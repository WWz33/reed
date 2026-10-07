//! `reed` CLI (`docs/architecture-spec.md` M0: clap; M4: `reed stats`).
//! The `filter` and `qc` subcommands follow with M5/M6.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};
use reed::stats::{assign_groups, default_kinds, run_stats, Levels, StatsRequest};
use reed_engine::budget::{MemBudget, ThreadBudget};

#[derive(Parser)]
#[command(
    name = "reed",
    version,
    about = "Population-genetics statistics and hard filtering, in pure Rust"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Quality-control report (§10): bcftools-stats-compatible output plus
    /// the structured JSON report.
    Qc {
        /// Input BGZF VCF (.vcf.gz)
        input: String,
        /// Write the structured JSON report here
        #[arg(long)]
        json: Option<String>,
        /// Output path for the bcftools-stats text (omit for stdout)
        #[arg(short, long)]
        out: Option<String>,
        /// Worker threads (default: all cores)
        #[arg(long)]
        threads: Option<usize>,
    },
    /// Streaming hard filter (v2 engine, §11): gold-parity funnel with
    /// worker-local BGZF compression and synchronized indexing.
    Filter {
        /// Input BGZF VCF (.vcf.gz)
        input: String,
        /// Samples to keep: comma-separated names (bcftools `-s`)
        #[arg(short = 's', long = "samples", value_name = "LIST")]
        samples: Option<String>,
        /// Samples to keep: one name per line (bcftools `-S`)
        #[arg(short = 'S', long = "samples-file", value_name = "FILE")]
        samples_file: Option<String>,
        /// Samples to keep (plink2/vcftools spelling)
        #[arg(long, value_name = "FILE")]
        keep: Option<String>,
        /// Samples to drop
        #[arg(long, value_name = "FILE")]
        remove: Option<String>,
        /// Sample lists naming samples absent from the header only warn
        #[arg(long)]
        force_samples: bool,
        /// Variant types to keep: snps,mnps,indels,other,ref (bcftools `-v`)
        #[arg(short = 'v', long = "types", value_name = "LIST")]
        types: Option<String>,
        /// Keep only SNP sites (`=just-acgt` restricts ALT to A/C/G/T)
        #[arg(long, value_name = "just-acgt", num_args = 0..=1,
              default_missing_value = "just-acgt")]
        snps_only: Option<String>,
        /// Minimum allele count, REF+ALT (bcftools `-m`)
        #[arg(short = 'm', long = "min-alleles", value_name = "INT")]
        min_alleles: Option<u32>,
        /// Maximum allele count, REF+ALT (bcftools `-M`)
        #[arg(short = 'M', long = "max-alleles", value_name = "INT")]
        max_alleles: Option<u32>,
        /// Keep records whose FILTER is in LIST (bcftools `-f`)
        #[arg(short = 'f', long = "apply-filters", value_name = "LIST")]
        apply_filters: Option<String>,
        /// bcftools expression: keep records that satisfy it (repeatable;
        /// 02-cli.md §4 subset)
        #[arg(short = 'i', long = "include", value_name = "EXPR")]
        include: Vec<String>,
        /// bcftools expression: drop records that satisfy it (repeatable)
        #[arg(short = 'e', long = "exclude", value_name = "EXPR")]
        exclude: Vec<String>,
        /// Mask genotypes with FMT/DP below INT (missing DP never masks)
        #[arg(long, value_name = "INT")]
        min_dp: Option<u32>,
        /// Mask genotypes with FMT/DP above INT
        #[arg(long, value_name = "INT")]
        max_dp: Option<u32>,
        /// Mask genotypes with FMT/GQ below INT
        #[arg(long, value_name = "INT")]
        min_gq: Option<u32>,
        /// Drop sites whose missing rate (after masking) exceeds FLOAT
        #[arg(long, value_name = "FLOAT")]
        geno: Option<f64>,
        /// Minimum minor allele frequency (after masking)
        #[arg(long, value_name = "FLOAT")]
        maf: Option<f64>,
        /// Minimum minor allele count (after masking)
        #[arg(long, value_name = "INT")]
        mac: Option<u32>,
        /// HWE exact-test p lower bound (`P`, optional `midp`)
        #[arg(long, value_name = "P [midp]", num_args = 1..=2)]
        hwe: Vec<String>,
        /// Remove these INFO tags (bcftools `-x`)
        #[arg(short = 'x', long = "remove-info", value_name = "LIST")]
        remove_info: Option<String>,
        /// Recompute and write these INFO tags (bcftools `+fill-tags`)
        #[arg(long, value_name = "LIST")]
        fill_tags: Option<String>,
        /// Output type: v (plain VCF) or z (BGZF VCF), with an optional
        /// compression level (`-Oz6`); b/u are refused until the BCF encoder
        /// lands
        #[arg(short = 'O', long = "output-type", value_name = "v|z[0-9]")]
        output_type: Option<String>,
        /// Output VCF (.vcf.gz) — omit for plain VCF on the path/stdout
        #[arg(short, long)]
        out: Option<String>,
        /// Write an index next to the output: csi (default) or tbi
        #[arg(short = 'W', long = "write-index",
              value_name = "csi|tbi", num_args = 0..=1, default_missing_value = "csi")]
        write_index: Option<String>,
        /// BGZF compression level (default 6)
        #[arg(long, default_value_t = 6)]
        compress_level: u8,
        /// JSON report path
        #[arg(long)]
        report: Option<String>,
        /// Test-only strategy override (single_pass|two_pass_reparse)
        #[arg(long, hide = true)]
        strategy: Option<String>,
        /// Worker threads (default: all cores)
        #[arg(long)]
        threads: Option<usize>,
    },

    /// Build a CSI (default) or TBI index for an existing BGZF VCF, like
    /// `bcftools index`.
    Index {
        /// Input BGZF VCF (.vcf.gz)
        input: String,
        /// Index flavor (default csi)
        #[arg(short = 't', long = "type", value_name = "csi|tbi", default_value = "csi")]
        index_type: String,
        /// Write the index here instead of next to the input
        #[arg(short, long)]
        out: Option<String>,
    },

    /// Windowed statistics (π, dxy, Fst, PBS, SFS, θ_W, Tajima's D) over a
    /// BGZF VCF, as additive accumulators with byte-deterministic output.
    Stats {
        /// Input BGZF VCF (.vcf.gz); `-` for stdin
        input: String,
        /// Groups file: `sample<TAB>population` rows, no header (pixy
        /// populations-file format); required by the windowed/pixy modes
        #[arg(short, long, alias = "pop", alias = "populations")]
        groups: Option<String>,
        /// π per population (Nei & Li 1979; missing data: Korunes & Samuk
        /// 2021 = pixy)
        #[arg(long)]
        pi: bool,
        /// dxy per population pair (Nei 1987; missing data: pixy)
        #[arg(long)]
        dxy: bool,
        /// Watterson's θ per population (Watterson 1975; per-allele-class
        /// mutation counts, pixy)
        #[arg(long)]
        theta: bool,
        /// Tajima's D per population (Tajima 1989; denominator once at the
        /// window-average allele count, pixy)
        #[arg(long = "tajima-d")]
        tajima_d: bool,
        /// Fst for all population pairs; `wc84` (Weir & Cockerham 1984,
        /// default) or `hudson` (Hudson 1992 / Bhatia 2013)
        #[arg(long, value_name = "wc84|hudson", num_args = 0..=1,
              default_missing_value = "wc84")]
        fst: Option<String>,
        /// PBS for the first three populations (Yi et al. 2010); needs
        /// three groups
        #[arg(long)]
        pbs: bool,
        /// Per-site allele frequencies (PREFIX.afreq, plink2 layout)
        #[arg(long)]
        freq: bool,
        /// Per-site HWE exact test and observed/expected heterozygosity
        /// (PREFIX.hardy, plink2 layout)
        #[arg(long)]
        hardy: bool,
        /// Per-sample inbreeding coefficient F (PREFIX.het, plink2 layout)
        #[arg(long)]
        het: bool,
        /// Per-sample and per-site missingness (PREFIX.smiss, PREFIX.vmiss)
        #[arg(long)]
        missing: bool,
        /// Per-sample and per-site mean depth (PREFIX.idepth,
        /// PREFIX.ldepth.mean, vcftools layouts; reads FMT/DP)
        #[arg(long)]
        depth: bool,
        /// QC summary in bcftools-stats format (PREFIX.stats)
        #[arg(long)]
        summary: bool,
        /// Report PBSn1 (Malaspinas et al. 2016) instead of the raw PBS
        #[arg(long)]
        pbs_normed: bool,
        /// Declare an all-variant input: window positions the file does not
        /// contain count as monomorphic with all samples called
        #[arg(long)]
        variant_only: bool,
        /// Compute over the sites the file contains (overestimates π/dxy)
        #[arg(long)]
        bypass_invariant_check: bool,
        /// Skip group-file rows naming samples absent from the header
        /// (default: error, pixy behavior)
        #[arg(long)]
        force_samples: bool,
        /// Contig lengths as a .fai-style file (name, length), for inputs
        /// whose header lacks ##contig lines
        #[arg(long)]
        contig_lengths_path: Option<String>,
        /// Window size in bp (accepts 100k, 1M, 100kb)
        #[arg(long = "window", alias = "size", default_value = "100000", value_parser = parse_size)]
        size: u64,
        /// Window step in bp (default: the window size)
        #[arg(long, default_value = "0", value_parser = parse_size)]
        step: u64,
        /// Output levels, comma-separated: window,contig,genome
        #[arg(long, default_value = "window")]
        levels: String,
        /// Skip windows that no site fell into
        #[arg(long)]
        skip_empty: bool,
        /// Write the TSV here instead of stdout
        #[arg(short, long)]
        out: Option<String>,
        /// Worker threads (default: all cores)
        #[arg(long)]
        threads: Option<usize>,
    },
}

/// Length parser for window/step: plain bp, or k/kb/M/Mb suffixes
/// (02-cli.md §2).
fn parse_size(text: &str) -> Result<u64, String> {
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
fn input_path(arg: &str) -> Option<PathBuf> {
    (arg != "-").then(|| PathBuf::from(arg))
}

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let command = Cli::parse().command;
    if let Command::Qc { input, json, out, threads } = &command {
        let threads = threads.unwrap_or_else(|| {
            std::thread::available_parallelism().expect("invariant: parallelism").get()
        });
        let started = std::time::Instant::now();
        let outcome = reed::qc_report::run_qc_report(
            input_path(input).as_deref(),
            &reed::qc_report::QcOptions {
                threads,
                json: json.clone().map(PathBuf::from),
            },
        )?;
        match out {
            Some(path) => std::fs::write(path, &outcome.bcftools_text)?,
            None => print!("{}", outcome.bcftools_text),
        }
        eprintln!(
            "reed qc: {} records, {} samples in {:.1}s",
            outcome.report.summary.n_records,
            outcome.report.n_samples,
            started.elapsed().as_secs_f64()
        );
        return Ok(());
    }
    if let Command::Index { input, index_type, out } = &command {
        let flavor = match index_type.as_str() {
            "csi" => reed::index::WriteIndex::Csi,
            "tbi" => reed::index::WriteIndex::Tbi,
            other => return Err(format!("unknown index type {other:?} (expected csi or tbi)").into()),
        };
        let started = std::time::Instant::now();
        let outcome = reed::index::run_index(
            Path::new(input),
            flavor,
            out.as_deref().map(Path::new),
        )?;
        eprintln!(
            "reed index: {} records -> {} in {:.1}s",
            outcome.n_records,
            outcome.path.display(),
            started.elapsed().as_secs_f64()
        );
        return Ok(());
    }
    if let Command::Filter {
        input,
        samples,
        samples_file,
        keep,
        remove,
        force_samples,
        types,
        snps_only,
        min_alleles,
        max_alleles,
        apply_filters,
        include,
        exclude,
        min_dp,
        max_dp,
        min_gq,
        geno,
        maf,
        mac,
        hwe,
        remove_info,
        fill_tags,
        output_type,
        out,
        write_index,
        compress_level,
        report,
        strategy,
        threads,
    } = command
    {
        return run_filter_cli(FilterCli {
            input,
            samples,
            samples_file,
            keep,
            remove,
            force_samples,
            types,
            snps_only,
            min_alleles,
            max_alleles,
            apply_filters,
            include,
            exclude,
            min_dp,
            max_dp,
            min_gq,
            geno,
            maf,
            mac,
            hwe,
            remove_info,
            fill_tags,
            output_type,
            out,
            write_index,
            compress_level,
            report,
            strategy,
            threads,
        });
    }
    let Command::Stats {
        input,
        groups,
        pi,
        dxy,
        theta,
        tajima_d,
        fst,
        pbs,
        pbs_normed,
        freq,
        hardy,
        het,
        missing,
        depth,
        summary,
        variant_only,
        bypass_invariant_check,
        force_samples,
        contig_lengths_path,
        size,
        step,
        levels,
        skip_empty,
        out,
        threads,
    } = command
    else {
        unreachable!("the filter subcommand returned above");
    };

    let (_source, facts) = reed::input::open(input_path(&input).as_deref())
        .map_err(|e| -> Box<dyn std::error::Error + Send + Sync> { e.to_string().into() })?;
    let sample_names: Vec<String> = facts.sample_names.clone();

    let mut contig_lengths: HashMap<Vec<u8>, u64> = facts
        .header
        .contigs()
        .iter()
        .filter_map(|(name, contig)| {
            contig
                .length()
                .and_then(|len| u64::try_from(len).ok())
                .map(|len| (name.as_bytes().to_vec(), len))
        })
        .collect();
    if let Some(path) = &contig_lengths_path {
        for line in std::fs::read_to_string(Path::new(path))?.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut cols = line.split('\t');
            let (Some(name), Some(len)) = (cols.next(), cols.next()) else {
                continue;
            };
            let len: u64 = len
                .parse()
                .map_err(|_| format!("{path}: contig length for {name:?} is not a number"))?;
            contig_lengths.insert(name.as_bytes().to_vec(), len);
        }
    }

    let threads = threads.unwrap_or_else(|| {
        std::thread::available_parallelism().expect("invariant: parallelism").get()
    });
    let budget = ThreadBudget { threads, mem: MemBudget::default_half_memory() };

    // Per-site/per-sample report mode (02-cli.md §5.4): PREFIX.afreq,
    // .hardy, .het, .smiss/.vmiss, .stats in the reference layouts.
    if freq || hardy || het || missing || depth || summary {
        let Some(prefix) = out.clone() else {
            return Err("report mode writes files: pass -o PREFIX".into());
        };
        let mut files = Vec::new();
        if freq || hardy || het || missing || depth {
            let request = reed::reports::ReportRequest {
                freq,
                hardy,
                het,
                missing,
                depth,
                prefix: PathBuf::from(&prefix),
            };
            let outcome = reed::reports::run_reports(
                input_path(&input).as_deref(),
                &request,
                &budget,
            )?;
            files.extend(outcome.files);
            eprintln!(
                "{} records in {:.1}s -> {}",
                outcome.n_records,
                // reports are fast; report wall time via the outcome count only
                0.0,
                files.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ")
            );
        }
        if summary {
            let outcome = reed::qc_report::run_qc_report(
                input_path(&input).as_deref(),
                &reed::qc_report::QcOptions { threads, json: None },
            )?;
            let mut path = std::ffi::OsString::from(&prefix);
            path.push(".stats");
            std::fs::write(std::path::PathBuf::from(&path), &outcome.bcftools_text)?;
            files.push(std::path::PathBuf::from(path));
        }
        return Ok(());
    }

    // Windowed modes need the populations file.
    let Some(groups) = groups else {
        return Err("a windowed statistic needs --pop FILE (sample<TAB>population, no header)".into());
    };
    let (gt_groups, group_names) =
        assign_groups(Path::new(&groups), &sample_names, force_samples)?;
    eprintln!(
        "{} samples, {} groups: {}",
        sample_names.len(),
        group_names.len(),
        group_names.join(", ")
    );

    // pixy-layout mode: one command writes PREFIX_pi.txt / PREFIX_dxy.txt
    // with pixy's columns (02-cli.md §5).
    if pi || dxy || theta || tajima_d || fst.is_some() || pbs {
        let step = if step == 0 { size } else { step };
        let invariant = if variant_only {
            reed::stats::InvariantMode::VariantOnly
        } else if bypass_invariant_check {
            eprintln!(
                "reed stats: warning: --bypass-invariant-check computes over present sites \
                 only; π and dxy are overestimated"
            );
            reed::stats::InvariantMode::Bypass
        } else {
            reed::stats::InvariantMode::Require
        };
        let pixy = reed::stats::PixyOptions { invariant };
        let fst_kind = match fst.as_deref() {
            None => None,
            Some("wc84") => Some(reed::stats::FstKind::WeirCockerham84),
            Some("hudson") => Some(reed::stats::FstKind::Hudson),
            Some(other) => {
                return Err(format!("unknown --fst method {other:?} (expected wc84 or hudson)").into());
            }
        };
        let pbs_kind = pbs.then_some(if pbs_normed {
            reed::stats::PbsKind::PbsN1
        } else {
            reed::stats::PbsKind::Yi2010
        });
        let kinds = reed::stats::pixy_kinds(
            group_names.len(),
            pi,
            dxy,
            theta,
            tajima_d,
            fst_kind,
            pbs_kind,
        );
        let request = StatsRequest {
            kinds: kinds.clone(),
            gt_groups,
            window_size: size,
            step,
            levels: Levels { window: true, contig: false, genome: false },
            skip_empty: false,
            contig_lengths,
            pixy: Some(pixy),
        };
        let report = run_stats(input_path(&input).as_deref(), request, &budget)?;
        let Some(prefix) = out else {
            return Err("pixy mode writes files: pass -o PREFIX".into());
        };
        let command_line = std::env::args().collect::<Vec<_>>().join(" ");
        let written = reed::pixy::write(
            Path::new(&prefix),
            &report,
            &group_names,
            &kinds,
            pixy,
            &command_line,
            env!("CARGO_PKG_VERSION"),
        )?;
        eprintln!(
            "{} windows, {} records in {:.1}s -> {}",
            report.n_windows,
            report.n_records,
            report.elapsed.as_secs_f64(),
            written
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
        return Ok(());
    }

    let request = StatsRequest {
        kinds: default_kinds(group_names.len()),
        gt_groups,
        window_size: size,
        step: if step == 0 { size } else { step },
        levels: Levels::parse(&levels)?,
        skip_empty,
        contig_lengths,
        pixy: None,
    };
    let report = run_stats(input_path(&input).as_deref(), request, &budget)?;
    eprintln!(
        "{} windows, {} records, {} multiallelic skipped, {} chunks in {:.1}s",
        report.n_windows,
        report.n_records,
        report.n_multiallelic_skipped,
        report.n_chunks,
        report.elapsed.as_secs_f64()
    );
    match out {
        Some(path) => std::fs::write(path, report.tsv)?,
        None => print!("{}", report.tsv),
    }
    Ok(())
}

struct FilterCli {
    input: String,
    samples: Option<String>,
    samples_file: Option<String>,
    keep: Option<String>,
    remove: Option<String>,
    force_samples: bool,
    types: Option<String>,
    snps_only: Option<String>,
    min_alleles: Option<u32>,
    max_alleles: Option<u32>,
    apply_filters: Option<String>,
    include: Vec<String>,
    exclude: Vec<String>,
    min_dp: Option<u32>,
    max_dp: Option<u32>,
    min_gq: Option<u32>,
    geno: Option<f64>,
    maf: Option<f64>,
    mac: Option<u32>,
    hwe: Vec<String>,
    remove_info: Option<String>,
    fill_tags: Option<String>,
    output_type: Option<String>,
    out: Option<String>,
    write_index: Option<String>,
    compress_level: u8,
    report: Option<String>,
    strategy: Option<String>,
    threads: Option<usize>,
}

/// Translates the CLI flags into a typed [`FilterSpecV2`] (02-cli.md §3, §6).
///
/// # Errors
///
/// Unknown tokens (variant types, output type) and malformed expressions.
fn spec_from_flags(cli: &FilterCli) -> Result<reed_filter::spec::FilterSpecV2, String> {
    use reed_filter::spec::{
        CohortRule, FilterSpecV2, GenotypeMaskRule, OnMissing,
    };
    let mut spec = FilterSpecV2::default();

    // ---- samples (cohort) ----
    if let Some(list) = &cli.samples {
        let names: Vec<String> = list
            .split(',')
            .map(str::trim)
            .filter(|n| !n.is_empty())
            .map(str::to_string)
            .collect();
        spec.cohort = Some(CohortRule::KeepNames(names));
    } else if let Some(path) = cli.samples_file.as_ref().or(cli.keep.as_ref()) {
        spec.cohort = Some(CohortRule::KeepList(path.clone()));
    } else if let Some(path) = &cli.remove {
        spec.cohort = Some(CohortRule::RemoveList(path.clone()));
    }

    // ---- site: variant class / alleles / FILTER ----
    if let Some(list) = &cli.types {
        for token in list.split(',') {
            let class = match token.trim() {
                "snps" | "snp" | "snv" | "snvs" => "snv",
                "mnps" | "mnp" | "mnv" => "mnv",
                "indels" | "indel" => "indel",
                "other" | "symbolic" => "symbolic",
                "ref" => "ref",
                other => return Err(format!("unknown variant type {other:?}")),
            };
            spec.site.variant_class.push(class.to_string());
        }
    }
    if let Some(modifier) = &cli.snps_only {
        if !matches!(modifier.as_str(), "just-acgt") {
            return Err(format!("unknown --snps-only modifier {modifier:?}"));
        }
        if !spec.site.variant_class.iter().any(|c| c == "snv") {
            spec.site.variant_class.push("snv".to_string());
        }
    }
    // bcftools counts REF+ALT, the spec counts ALTs.
    if let Some(m) = cli.min_alleles {
        spec.site.min_alt_alleles = m.saturating_sub(1);
    }
    if let Some(m) = cli.max_alleles {
        spec.site.max_alt_alleles = m.saturating_sub(1);
    }
    if let Some(list) = &cli.apply_filters {
        spec.site.require_filter = list
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .collect();
    }

    // ---- expressions (02-cli.md §4): route by genotype dependence ----
    for (text, include) in cli
        .include
        .iter()
        .map(|t| (t, true))
        .chain(cli.exclude.iter().map(|t| (t, false)))
    {
        let expr = reed_filter::expr::parse(text)
            .map_err(|e| format!("{} {text:?}: {e}", if include { "-i" } else { "-e" }))?;
        let genotype_derived = expr.references_genotypes();
        match (include, genotype_derived) {
            (true, false) => spec.site.include.push(text.clone()),
            (true, true) => spec.site_after_mask.include.push(text.clone()),
            (false, false) => spec.site.exclude.push(text.clone()),
            (false, true) => spec.site_after_mask.exclude.push(text.clone()),
        }
    }

    // ---- genotype masking (plink2 semantics: missing never masks) ----
    let mut mask = |field: &str, when: String| {
        spec.genotype.mask_if.push(GenotypeMaskRule {
            field: field.to_string(),
            when,
            on_missing: OnMissing::Keep,
        });
    };
    if let Some(v) = cli.min_dp {
        mask("DP", format!("< {v}"));
    }
    if let Some(v) = cli.max_dp {
        mask("DP", format!("> {v}"));
    }
    if let Some(v) = cli.min_gq {
        mask("GQ", format!("< {v}"));
    }

    // ---- after-mask site rules ----
    if let Some(g) = cli.geno {
        spec.site_after_mask.max_missing = Some(g);
    }
    if let Some(m) = cli.maf {
        spec.site_after_mask.min_maf = Some(m);
    }
    if let Some(m) = cli.mac {
        spec.site_after_mask.min_mac = m;
    }
    if !cli.hwe.is_empty() {
        let p: f64 = cli.hwe[0]
            .parse()
            .map_err(|_| format!("invalid --hwe p-value {:?}", cli.hwe[0]))?;
        spec.site_after_mask.min_hwe_p = Some(p);
        if cli.hwe.get(1).is_some_and(|m| m == "midp") {
            spec.site_after_mask.hwe_mode = reed_filter::spec::HweMode::MidP;
        }
    }

    // ---- INFO ----
    if let Some(list) = &cli.remove_info {
        spec.recompute.drop_stale = list
            .split(',')
            .map(str::trim)
            .map(|t| t.strip_prefix("INFO/").unwrap_or(t).to_string())
            .filter(|t| !t.is_empty())
            .collect();
    }
    if let Some(list) = &cli.fill_tags {
        spec.recompute.tags = list
            .split(',')
            .map(str::trim)
            .filter(|t| !t.is_empty())
            .map(str::to_string)
            .collect();
    }
    Ok(spec)
}

fn run_filter_cli(cli: FilterCli) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let FilterCli {
        input,
        out,
        write_index,
        compress_level,
        report,
        strategy,
        threads,
        output_type,
        force_samples,
        ..
    } = &cli;
    let spec_v2 = spec_from_flags(&cli)?;
    let strategy = strategy
        .as_deref()
        .map(|s| match s {
            "single_pass" => Ok(reed_filter::plan::ExecStrategy::SinglePass),
            "two_pass_reparse" => Ok(reed_filter::plan::ExecStrategy::TwoPassReparse),
            "sketch_then_emit" => Ok(reed_filter::plan::ExecStrategy::SketchThenEmit),
            other => Err(format!("unknown strategy {other:?}")),
        })
        .transpose()?;
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let threads = threads.unwrap_or_else(|| {
        std::thread::available_parallelism().expect("invariant: parallelism").get()
    });
    // -O v|z[level]: explicit output type; the extension is inferred when
    // it is absent. b/u are refused until the BCF encoder lands (R6).
    let (output, compress_level) = match output_type.as_deref() {
        None => (out.clone().map(PathBuf::from), *compress_level),
        Some(spec) => {
            let (kind, level) = spec.split_at(1);
            let level = if level.is_empty() {
                *compress_level
            } else {
                level.parse::<u8>().map_err(|_| format!("invalid -O level in {spec:?}"))?
            };
            match kind {
                "v" => (None, level),
                "z" => {
                    let path = out.clone().ok_or("the compressed output type needs -o FILE")?;
                    let path = PathBuf::from(path);
                    let path = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
                        path
                    } else {
                        PathBuf::from(format!("{}.gz", path.display()))
                    };
                    (Some(path), level)
                }
                other => {
                    return Err(format!(
                        "output type {other:?} is not supported yet (use v or z; BCF is R6)"
                    )
                    .into())
                }
            }
        }
    };
    let opts = reed::filter::FilterOptions {
        threads,
        output,
        write_index: reed::filter::WriteIndex::parse(write_index.as_deref())?,
        compress_level,
        strategy,
        command_line,
        quiet: false,
        report_json: report.as_ref().map(PathBuf::from),
        force_samples: *force_samples,
    };
    let report = reed::filter::run_filter(input_path(&input).as_deref(), &spec_v2, &opts)?;
    eprintln!(
        "reed filter: {} -> {} records in {:.1}s",
        report.n_input_records, report.n_output_records, report.elapsed_secs
    );
    Ok(())
}
