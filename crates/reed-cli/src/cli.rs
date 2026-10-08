//! The clap argument groups (02-cli.md §1–§5): one `Args` struct per
//! subcommand, no translation logic.

use clap::{Args, Parser, Subcommand};

use crate::units::parse_size;

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
    /// Quality-control report
    Qc(QcCli),
    /// Hard filtering
    Filter(FilterCli),
    /// Window statistics
    Stats(StatsCli),
    /// Build a CSI/TBI index
    Index(IndexCli),
}

/// Quality-control report (§10): bcftools-stats-compatible output plus
/// the structured JSON report.
#[derive(Args)]
pub struct QcCli {
    /// Input BGZF VCF (.vcf.gz)
    pub input: String,
    /// Write the structured JSON report here
    #[arg(long)]
    pub json: Option<String>,
    /// Output path for the bcftools-stats text (omit for stdout)
    #[arg(short, long)]
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
    /// Regions to keep, comma-separated (`Chr01`, `Chr01:100-200`);
    /// needs an index next to the input (bcftools `-r`)
    #[arg(short = 'r', long = "regions", value_name = "LIST")]
    pub regions: Option<String>,
    /// Regions to keep, one per line; needs an index (bcftools `-R`)
    #[arg(short = 'R', long = "regions-file", value_name = "FILE")]
    pub regions_file: Option<String>,
    /// Regions to keep, streamed without an index (bcftools `-t`)
    #[arg(short = 't', long = "targets", value_name = "LIST")]
    pub targets: Option<String>,
    /// Regions to keep, one per line (bcftools `-T`); `^` prefix
    /// excludes
    #[arg(short = 'T', long = "targets-file", value_name = "FILE")]
    pub targets_file: Option<String>,
    /// Samples to keep: comma-separated names; `^` prefix excludes
    /// (`^A,B` = keep all but A and B; bcftools `-s`)
    #[arg(short = 's', long = "samples", value_name = "[^]LIST")]
    pub samples: Vec<String>,
    /// Samples to keep: one name per line (bcftools `-S`)
    #[arg(short = 'S', long = "samples-file", value_name = "FILE")]
    pub samples_file: Option<String>,
    /// Samples to keep, one name per line (plink2/vcftools spelling)
    #[arg(long, value_name = "FILE")]
    pub keep: Option<String>,
    /// Samples to drop, one name per line
    #[arg(long, value_name = "FILE")]
    pub remove: Option<String>,
    /// Drop samples whose missing rate exceeds this (plink2 semantics:
    /// computed after the genotype mask, over the sites that pass the
    /// pre-mask rules; applied before --geno/--maf/--hwe)
    #[arg(long, value_name = "FLOAT")]
    pub mind: Option<f64>,
    /// Sample lists naming samples absent from the header only warn
    #[arg(long)]
    pub force_samples: bool,
    /// Variant types to keep: snps,mnps,indels,other,ref (bcftools `-v`)
    #[arg(short = 'v', long = "types", value_name = "LIST")]
    pub types: Option<String>,
    /// Variant types to drop (bcftools `-V/--exclude-types`)
    #[arg(short = 'V', long = "exclude-types", value_name = "LIST")]
    pub exclude_types: Option<String>,
    /// Keep only SNP sites (`=just-acgt` restricts ALT to A/C/G/T)
    #[arg(long, value_name = "just-acgt", num_args = 0..=1,
          default_missing_value = "just-acgt")]
    pub snps_only: Option<String>,
    /// Minimum allele count, REF+ALT (bcftools `-m`)
    #[arg(short = 'm', long = "min-alleles", value_name = "INT")]
    pub min_alleles: Option<u32>,
    /// Maximum allele count, REF+ALT (bcftools `-M`)
    #[arg(short = 'M', long = "max-alleles", value_name = "INT")]
    pub max_alleles: Option<u32>,
    /// Keep records whose FILTER is in LIST (bcftools `-f`)
    #[arg(short = 'f', long = "apply-filters", value_name = "LIST")]
    pub apply_filters: Option<String>,
    /// bcftools expression: keep records that satisfy it (repeatable;
    /// 02-cli.md §4 subset)
    #[arg(short = 'i', long = "include", value_name = "EXPR")]
    pub include: Vec<String>,
    /// bcftools expression: drop records that satisfy it (repeatable)
    #[arg(short = 'e', long = "exclude", value_name = "EXPR")]
    pub exclude: Vec<String>,
    /// Mask genotypes with FMT/DP below INT (missing DP never masks)
    #[arg(long, value_name = "INT")]
    pub min_dp: Option<u32>,
    /// Mask genotypes with FMT/DP above INT
    #[arg(long, value_name = "INT")]
    pub max_dp: Option<u32>,
    /// Mask genotypes with FMT/GQ below INT
    #[arg(long, value_name = "INT")]
    pub min_gq: Option<u32>,
    /// Mask from the -i/-e FORMAT expressions: set the genotype to this
    /// value (only `.`, bcftools filter -S .)
    #[arg(long, value_name = ".")]
    pub set_gts: Option<String>,
    /// Drop sites whose missing rate (after masking) exceeds FLOAT
    #[arg(long, value_name = "FLOAT")]
    pub geno: Option<f64>,
    /// Minimum minor allele frequency (after masking)
    #[arg(long, value_name = "FLOAT")]
    pub maf: Option<f64>,
    /// Maximum minor allele frequency (after masking)
    #[arg(long, value_name = "FLOAT")]
    pub max_maf: Option<f64>,
    /// Minimum minor allele count (after masking)
    #[arg(long, value_name = "INT")]
    pub mac: Option<u32>,
    /// Maximum minor allele count (after masking)
    #[arg(long, value_name = "INT")]
    pub max_mac: Option<u32>,
    /// Do not update INFO/AC,INFO/AN when samples are dropped
    /// (bcftools `view -I`; the update is the default)
    #[arg(short = 'I', long = "no-update")]
    pub no_update: bool,
    /// HWE exact-test p lower bound (`P`, optional `midp`)
    #[arg(long, value_name = "P [midp]", num_args = 1..=2)]
    pub hwe: Vec<String>,
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
    /// Output VCF (.vcf.gz) — omit for plain VCF on the path/stdout
    #[arg(short, long)]
    pub out: Option<String>,
    /// Write an index next to the output: csi (default) or tbi
    #[arg(short = 'W', long = "write-index",
          value_name = "csi|tbi", num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,
    /// BGZF compression level (default 6)
    #[arg(long, default_value_t = 6)]
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
    /// Directory for temporary files (the stdin spill of multi-pass runs)
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<String>,
}


/// Build a CSI (default) or TBI index for an existing BGZF VCF.
#[derive(Args)]
pub struct IndexCli {
    /// Input VCF/VCF.gz
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
    /// Input BGZF VCF (.vcf.gz); `-` for stdin
    pub input: String,
    /// Groups file: `sample<TAB>population` rows, no header (pixy
    /// populations-file format); required by the windowed/pixy modes
    #[arg(short, long, alias = "pop", alias = "populations")]
    pub groups: Option<String>,
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
    /// PBS for the first three populations (Yi et al. 2010); needs
    /// three groups
    #[arg(long)]
    pub pbs: bool,
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
    #[arg(long)]
    pub contig_lengths_path: Option<String>,
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
    /// Write the TSV here instead of stdout
    #[arg(short, long)]
    pub out: Option<String>,
    /// Worker threads (default: all cores)
    #[arg(long)]
    pub threads: Option<usize>,
    /// Directory for temporary files (the stdin spill of report runs that
    /// read twice)
    #[arg(long, value_name = "DIR")]
    pub tmp_dir: Option<String>,
}

