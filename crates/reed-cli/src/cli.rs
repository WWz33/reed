//! The clap argument groups (02-cli.md §1–§5): one `Args` struct per
//! subcommand, no translation logic.

use clap::{Args, Parser, Subcommand};

use crate::units::parse_size;

/// The shared `<INPUT>` help text (doc_v4 03-architecture.md §3.6): one
/// constant for `filter` and `stats` so the accepted formats cannot drift
/// (01 F4-13). `index` states its file-only restriction separately.
const INPUT_HELP: &str = "Input VCF/VCF.gz/BCF, or `-` for stdin";

pub(crate) mod filter_flags;
pub(crate) mod verbosity;

use verbosity::Verbosity;

/// The `reed` command line.
#[derive(Parser)]
#[command(
    name = "reed",
    version,
    about = "Population-genetics statistics and hard filtering, in pure Rust"
)]
pub(crate) struct Cli {
    /// The subcommand to run.
    #[command(subcommand)]
    pub command: Command,
}

/// Subcommands.
// Clap's derive builds one owned variant per subcommand; the enum is parsed
// exactly once at startup, so the padding between variants is irrelevant.
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
pub(crate) enum Command {
    /// Hard filtering
    Filter(FilterCli),
    /// Window statistics
    Stats(StatsCli),
    /// Build a CSI/TBI index
    Index(IndexCli),
}

/// Hard filtering (02-cli.md §3): the bcftools `view` + `filter` + `annotate -x`
/// flags in one command.
#[derive(Args)]
pub(crate) struct FilterCli {
    #[arg(help = INPUT_HELP)]
    pub input: String,
    #[command(flatten)]
    pub regions_args: filter_flags::RegionArgs,
    #[command(flatten)]
    pub sample_args: filter_flags::SampleArgs,
    /// Keep samples whose missing rate (missing alleles / 2·sites, after the
    /// genotype mask, over the sites that pass the pre-mask rules) is at
    /// most FLOAT; applied before --geno/--maf/--hwe
    #[arg(long, value_name = "FLOAT", value_parser = reed::spec::MissingRate::parse)]
    pub mind: Option<reed::spec::MissingRate>,
    /// Withdrawn flag: it stated the fraction that must be called; the
    /// replacement names the fraction that may be missing
    #[arg(long, value_name = "FLOAT", hide = true,
          value_parser = filter_flags::removed_max_missing_sample)]
    pub max_missing_sample: Option<String>,
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
    /// Write an index next to the output: `-W` or `-W=csi` (default), or
    /// `-W=tbi`; omit for no index (bcftools semantics)
    #[arg(short = 'W', long = "write-index",
          value_name = "csi|tbi", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,
    /// BGZF compression level (default 6); prefer `-Oz<level>`/`-Ob<level>`
    #[arg(long, default_value_t = 6, hide = true)]
    pub compress_level: u8,
    /// JSON report path
    #[arg(long)]
    pub report: Option<String>,
    /// Write the per-phase timing report here (03-performance.md O10)
    #[arg(long, value_name = "FILE")]
    pub timing: Option<String>,
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
    /// Output level: 0 silent, 1 summary (default), 2 stage timing, 3 debug
    #[arg(long, value_name = "INT", default_value_t = Verbosity::Summary,
          value_parser = Verbosity::parse)]
    pub verbosity: Verbosity,
    /// Do not write the ##reed_* provenance header lines
    #[arg(long)]
    pub no_version: bool,
    /// Do not write the ##reed_command line (the body stays deterministic;
    /// the recorded argv carries machine-specific paths)
    #[arg(long)]
    pub no_command: bool,
}


/// Build a CSI (default) or TBI index for an existing BGZF VCF.
#[derive(Args)]
pub(crate) struct IndexCli {
    /// Input VCF/VCF.gz/BCF (file only: the index needs a seekable input)
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
    /// Output level: 0 silent, 1 summary (default), 2 stage timing, 3 debug
    #[arg(long, value_name = "INT", default_value_t = Verbosity::Summary,
          value_parser = Verbosity::parse)]
    pub verbosity: Verbosity,
}


/// Window statistics (02-cli.md §5): the statistics are named after the method,
/// never after a tool.
#[derive(Args)]
pub(crate) struct StatsCli {
    #[arg(help = INPUT_HELP)]
    pub input: String,
    #[command(flatten)]
    pub regions_args: filter_flags::RegionArgs,
    #[command(flatten)]
    pub sample_args: filter_flags::SampleArgs,
    #[command(flatten)]
    pub site_args: filter_flags::SiteFilterArgs,
    /// Not available here (doc_v4 02 回退版 §3.2): filter the samples with
    /// `reed filter --mind` first
    #[arg(long, value_name = "FLOAT", hide = true,
          value_parser = filter_flags::unsupported_stats_mind)]
    pub mind: Option<String>,
    /// Withdrawn flag: it stated the fraction that must be called
    #[arg(long, value_name = "FLOAT", hide = true,
          value_parser = filter_flags::removed_max_missing_sample)]
    pub max_missing_sample: Option<String>,
    /// Groups file: `sample<TAB>population` rows, no header (pixy
    /// populations-file format); required by the windowed/pixy modes
    #[arg(long = "pop", visible_alias = "populations", alias = "groups",
          value_name = "FILE")]
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
    /// Report the global (r-population) WC84 Fst instead of the pairwise
    /// values: one value over all --pop populations (docs/rounds/v5/05 §3.2)
    #[arg(long, requires = "fst")]
    pub fst_global: bool,
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
    /// Patterson f2(A, B); repeatable
    #[arg(long, value_name = "A,B")]
    pub f2: Vec<String>,
    /// Patterson f3(C; A, B); repeatable
    #[arg(long, value_name = "C;A,B")]
    pub f3: Vec<String>,
    /// Patterson f4(A, B; C, D); repeatable
    #[arg(long, value_name = "A,B;C,D")]
    pub f4: Vec<String>,
    /// Patterson's D(A, B; C, D); repeatable
    #[arg(long, value_name = "A,B;C,D")]
    pub dstat: Vec<String>,
    /// f4-ratio A,B;C,D;X (admixtools alpha estimate); repeatable
    #[arg(long = "f4-ratio", value_name = "A,B;C,D;X")]
    pub f4_ratio: Vec<String>,
    /// Block-jackknife block size in bp for the f-statistics (default 5M)
    #[arg(long = "block-size", value_name = "SIZE", value_parser = crate::units::parse_size)]
    pub block_size: Option<u64>,
    /// Projected site-frequency spectrum, genome-wide (docs/rounds/v5/05
    /// §3.1): writes PREFIX_sfs/dadi and PREFIX_sfs/fastsimcoal2; needs
    /// --proj. Default unfolded (REF = ancestral); `--sfs folded` folds it
    #[arg(long, value_name = "folded", num_args = 0..=1, default_missing_value = "unfolded",
          value_parser = ["unfolded", "folded"])]
    pub sfs: Option<String>,
    /// Projection sizes in chromosomes (--sfs): one value broadcasts to
    /// every population, or one per population in first-appearance order
    #[arg(long, value_name = "N[,N...]")]
    pub proj: Option<String>,
    /// Round the projected spectrum to integers (easySFS `--dtype int`)
    #[arg(long)]
    pub sfs_int: bool,
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
    /// Directory for temporary files (the stdin spill for report runs)
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<String>,
    /// In-flight memory budget in MiB (default: half of physical memory)
    #[arg(long, value_name = "MB")]
    pub memory: Option<u64>,
    /// Output level: 0 silent, 1 summary (default), 2 stage timing, 3 debug
    #[arg(long, value_name = "INT", default_value_t = Verbosity::Summary,
          value_parser = Verbosity::parse)]
    pub verbosity: Verbosity,
}

