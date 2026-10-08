//! The clap argument groups (02-cli.md §1–§5): one `Args` struct per
//! subcommand, no translation logic.

use clap::{Args, Parser, Subcommand};

use crate::units::parse_size;

pub mod filter_flags;

/// The `reed` command line.
#[derive(Parser)]
#[command(
    name = "reed",
    version,
    about = "Population-genetics statistics and hard filtering, in pure Rust"
)]
pub struct Cli {
    /// The subcommand to run.
    #[command(subcommand)]
    pub command: Command,
}

/// Subcommands.
// Clap's derive builds one owned variant per subcommand; the enum is parsed
// exactly once at startup, so the padding between variants is irrelevant.
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
pub enum Command {
    /// Deprecated alias of `stats --summary` (hidden; removed after one
    /// version per 04-cleanup §C4)
    #[command(hide = true)]
    Qc(QcCli),
    /// Hard filtering
    Filter(FilterCli),
    /// Window statistics
    Stats(StatsCli),
    /// Build a CSI/TBI index
    Index(IndexCli),
}

/// Deprecated `reed qc` (§10): the bcftools-stats-compatible report now
/// lives behind `reed stats --summary` (04-cleanup §C4).
#[derive(Args)]
pub struct QcCli {
    /// Input VCF/VCF.gz/BCF, or `-` for stdin
    pub input: String,
    /// Write the structured JSON report here
    #[arg(long)]
    pub json: Option<String>,
    /// Output path for the bcftools-stats text (omit for stdout)
    #[arg(short = 'o', long = "output", visible_alias = "out", value_name = "FILE")]
    pub out: Option<String>,
    /// Worker threads (default: all cores)
    #[arg(long)]
    pub threads: Option<usize>,
}


/// Hard filtering (02-cli.md §3): the bcftools `view` + `filter` + `annotate -x`
/// flags in one command.
#[derive(Args)]
pub struct FilterCli {
    /// Input BGZF VCF (.vcf.gz)
    pub input: String,
    #[command(flatten)]
    pub regions_args: filter_flags::RegionArgs,
    #[command(flatten)]
    pub sample_args: filter_flags::SampleArgs,
    /// Drop samples whose missing rate exceeds this (plink2 semantics:
    /// computed after the genotype mask, over the sites that pass the
    /// pre-mask rules; applied before --geno/--maf/--hwe)
    #[arg(long, value_name = "FLOAT")]
    pub mind: Option<f64>,
    /// Sample lists naming samples absent from the header only warn
    #[arg(long)]
    pub force_samples: bool,
    #[command(flatten)]
    pub site_args: filter_flags::SiteFilterArgs,
    /// Do not update INFO/AC,INFO/AN when samples are dropped
    /// (bcftools `view -I`; the update is the default)
    #[arg(short = 'I', long = "no-update")]
    pub no_update: bool,
    /// Remove these INFO tags (bcftools `-x`)
    #[arg(short = 'x', long = "remove-info", value_name = "LIST")]
    pub remove_info: Option<String>,
    /// Recompute and write these INFO tags (bcftools `+fill-tags`)
    #[arg(long, value_name = "LIST")]
    pub fill_tags: Option<String>,
    /// Output type: v (plain VCF), z (BGZF VCF), b (BCF) or u
    /// (uncompressed BCF); optional compression level (`-Oz6`, `-Ob1`)
    #[arg(short = 'O', long = "output-type", value_name = "v|z|b|u[0-9]")]
    pub output_type: Option<String>,
    /// Output file (default: stdout); the format comes from -O or the
    /// extension
    #[arg(short = 'o', long = "output", visible_alias = "out", value_name = "FILE")]
    pub out: Option<String>,
    /// Write an index next to the output: csi (default) or tbi
    #[arg(short = 'W', long = "write-index",
          value_name = "csi|tbi", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,
    /// BGZF compression level (default 6); prefer `-Oz<level>`/`-Ob<level>`
    #[arg(long, default_value_t = 6, hide = true)]
    pub compress_level: u8,
    /// JSON report path
    #[arg(long)]
    pub report: Option<String>,
    /// Test-only strategy override (single_pass|two_pass_reparse)
    #[arg(long, hide = true)]
    pub strategy: Option<String>,
    /// Worker threads (default: all cores)
    #[arg(long)]
    pub threads: Option<usize>,
    /// In-flight memory budget in MiB (default: half of physical memory)
    #[arg(long, value_name = "MB")]
    pub memory: Option<u64>,
    /// Directory for temporary files (the stdin spill of multi-pass runs)
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<String>,
    /// 0 silent, 1 summary (default), 2 per-stage, 3 debug
    #[arg(long, value_name = "INT", default_value_t = 1)]
    pub verbosity: u8,
    /// Do not write the ##reed_* provenance header lines
    #[arg(long)]
    pub no_version: bool,
}


/// Build a CSI (default) or TBI index for an existing BGZF VCF.
#[derive(Args)]
pub struct IndexCli {
    /// Input VCF/VCF.gz/BCF, or `-` for stdin
    pub input: String,
    /// Write a CSI index (default; bcftools `-c`)
    #[arg(short = 'c', long = "csi", conflicts_with = "tbi")]
    pub csi: bool,
    /// Write a TBI index instead (bcftools `-t`)
    #[arg(short = 't', long = "tbi")]
    pub tbi: bool,
    /// Write the index here instead of next to the input (bcftools `-o`)
    #[arg(short, long)]
    pub output: Option<String>,
    /// Overwrite an existing index (bcftools `-f`)
    #[arg(short, long)]
    pub force: bool,
    /// Worker threads (default: all cores)
    #[arg(long)]
    pub threads: Option<usize>,
}


/// Window statistics (02-cli.md §5): the statistics are named after the method,
/// never after a tool.
#[derive(Args)]
pub struct StatsCli {
    /// Input VCF/VCF.gz/BCF, or `-` for stdin
    pub input: String,
    #[command(flatten)]
    pub regions_args: filter_flags::RegionArgs,
    #[command(flatten)]
    pub sample_args: filter_flags::SampleArgs,
    #[command(flatten)]
    pub site_args: filter_flags::SiteFilterArgs,
    /// Groups file: `sample<TAB>population` rows, no header (pixy
    /// populations-file format); required by the windowed/pixy modes
    #[arg(long = "pop", visible_alias = "populations", short = 'g',
          hide_short_help = true, alias = "groups", value_name = "FILE")]
    pub pop: Option<String>,
    /// π per population (Nei & Li 1979; missing data: Korunes & Samuk
    /// 2021 = pixy)
    #[arg(long)]
    pub pi: bool,
    /// dxy per population pair (Nei 1987; missing data: pixy)
    #[arg(long)]
    pub dxy: bool,
    /// Watterson's θ per population (Watterson 1975; per-allele-class
    /// mutation counts, pixy)
    #[arg(long)]
    pub theta: bool,
    /// Tajima's D per population (Tajima 1989; denominator once at the
    /// window-average allele count, pixy)
    #[arg(long = "tajima-d")]
    pub tajima_d: bool,
    /// Fst for all population pairs; `wc84` (Weir & Cockerham 1984,
    /// default) or `hudson` (Hudson 1992 / Bhatia 2013)
    #[arg(long, value_name = "wc84|hudson", num_args = 0..=1,
          default_missing_value = "wc84")]
    pub fst: Option<String>,
    /// PBS for the named triple `A,B,C` (A is the branch of interest,
    /// Yi et al. 2010); repeatable
    #[arg(long, value_name = "A,B,C")]
    pub pbs: Vec<String>,
    /// Per-site allele frequencies (PREFIX.afreq, plink2 layout)
    #[arg(long)]
    pub freq: bool,
    /// Per-site HWE exact test and observed/expected heterozygosity
    /// (PREFIX.hardy, plink2 layout)
    #[arg(long)]
    pub hardy: bool,
    /// Per-sample inbreeding coefficient F (PREFIX.het, plink2 layout)
    #[arg(long)]
    pub het: bool,
    /// Per-sample and per-site missingness (PREFIX.smiss, PREFIX.vmiss)
    #[arg(long)]
    pub missing: bool,
    /// Per-sample and per-site mean depth (PREFIX.idepth,
    /// PREFIX.ldepth.mean, vcftools layouts; reads FMT/DP)
    #[arg(long)]
    pub depth: bool,
    /// QC summary in bcftools-stats format (PREFIX.stats)
    #[arg(long)]
    pub summary: bool,
    /// Report PBSn1 (Malaspinas et al. 2016) instead of the raw PBS
    #[arg(long)]
    pub pbs_normed: bool,
    /// Declare an all-variant input: window positions the file does not
    /// contain count as monomorphic with all samples called
    #[arg(long)]
    pub variant_only: bool,
    /// Compute over the sites the file contains (overestimates π/dxy)
    #[arg(long)]
    pub bypass_invariant_check: bool,
    /// Skip group-file rows naming samples absent from the header
    /// (default: error, pixy behavior)
    #[arg(long)]
    pub force_samples: bool,
    /// Contig lengths as a .fai-style file (name, length), for inputs
    /// whose header lacks ##contig lines
    #[arg(long = "contig-lengths", alias = "contig-lengths-path", value_name = "FILE")]
    pub contig_lengths: Option<String>,
    /// Window size in bp (accepts 100k, 1M, 100kb)
    #[arg(long = "window", alias = "size", default_value = "100000", value_parser = parse_size)]
    pub size: u64,
    /// Window step in bp (default: the window size)
    #[arg(long, default_value = "0", value_parser = parse_size)]
    pub step: u64,
    /// Output levels, comma-separated: window,contig,genome
    #[arg(long, default_value = "window")]
    pub levels: String,
    /// Skip windows that no site fell into
    #[arg(long)]
    pub skip_empty: bool,
    /// Output prefix: PREFIX_pi.txt, PREFIX.afreq, ... (or the TSV on stdout)
    #[arg(short = 'o', long = "output", visible_alias = "out", value_name = "PREFIX")]
    pub out: Option<String>,
    /// Worker threads (default: all cores)
    #[arg(long)]
    pub threads: Option<usize>,
    /// Directory for temporary files (the stdin spill of report runs that
    /// read twice)
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<String>,
    /// In-flight memory budget in MiB (default: half of physical memory)
    #[arg(long, value_name = "MB")]
    pub memory: Option<u64>,
    /// 0 silent, 1 summary (default), 2 per-stage, 3 debug
    #[arg(long, value_name = "INT", default_value_t = 1)]
    pub verbosity: u8,
}

