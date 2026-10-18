//! The clap argument groups (02-cli.md §1–§5): one `Args` struct per
//! subcommand, no translation logic.

use clap::{Args, Parser, Subcommand};

use crate::units::parse_size;

/// `--ld-decay` accepts a bare flag (the default span) or a size.
fn parse_optional_size(text: &str) -> Result<String, String> {
    if text.is_empty() {
        return Ok(String::new());
    }
    parse_size(text).map(|_| text.to_string())
}

/// The shared `<INPUT>` help text (doc_v4 03-architecture.md §3.6): one
/// constant for `filter` and `stats` so the accepted formats cannot drift
/// (01 F4-13). `index` states its file-only restriction separately.
const INPUT_HELP: &str = "Input VCF/VCF.gz/BCF, or `-` for stdin";

pub(crate) mod filter_flags;
pub(crate) mod verbosity;

use verbosity::Verbosity;

/// The `reed` command line.
///
/// `--version` prints the crate version plus the exact source state the
/// binary was built from (E6-11): `reed 0.3.0 (v0.3.0-12-g8709f1f 2026-10-10,
/// rustc 1.99.0)`. The pieces come from `build.rs`; a source tarball without
/// git shows `unknown` for the commit/date.
const LONG_VERSION: &str = concat!(
    env!("CARGO_PKG_VERSION"),
    " (",
    env!("REED_BUILD_COMMIT"),
    " ",
    env!("REED_BUILD_DATE"),
    ", ",
    env!("REED_BUILD_RUSTC"),
    ")"
);

/// The `reed` command line.
#[derive(Parser)]
#[command(
    name = "reed",
    version,
    long_version = LONG_VERSION,
    about = "Population-genetics statistics and hard filtering, in pure Rust",
    // doc_v14 02 §1.3: the flat `{subcommands}` section is replaced by the
    // grouped list (`render_grouped_commands`, set as `after_help` by the
    // parser below); `{options}` still renders clap's own options.
    help_template = "{before-help}{about-with-newline}\n{usage-heading} {usage}\n\n{after-help}{options}"
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
    /// Hard filtering: selection, genotype masking, QC thresholds, LD pruning
    Filter(FilterCli),
    /// Per-site and per-sample QC: freq, hardy, het, missing, depth, summary
    Stats(StatsCli),
    /// Windowed population genetics: pi, dxy, Fst, PBS, SFS, f-stats, LD decay
    #[command(visible_alias = "popgen")]
    Popg(PopgCli),
    /// Build a CSI/TBI index
    Index(IndexCli),
    /// Convert to phylogenetic matrices, genotype matrices, hapmap, MSA, SFS
    Convert(ConvertCli),
    /// Edit identifiers: sample names, contig names, variant IDs, header lines
    #[command(visible_alias = "reheader")]
    Rename(RenameCli),
}

/// The help groups (doc_v14 02 §1.2): the single source of the grouping.
/// The order here is the order `reed --help` prints; a test asserts every
/// `Command` variant appears in exactly one group (plus clap's automatic
/// `help`, which is not an enum variant).
pub(crate) const GROUPS: &[(&str, &[&str])] = &[
    ("Filtering", &["filter"]),
    ("Statistics", &["stats", "popg"]),
    ("Transform", &["convert", "rename"]),
    ("File Operations", &["index"]),
    ("Miscellaneous", &["help"]),
];

/// Renders the grouped command list for `reed --help` (doc_v14 02 §1.3):
/// each group's title, then `name` padded to the width clap uses, then the
/// subcommand's own `about` text (taken from clap, never copied here, so the
/// two cannot drift).
pub(crate) fn render_grouped_commands() -> String {
    use clap::CommandFactory as _;
    use std::fmt::Write as _;
    let root = Cli::command();
    let mut out = String::new();
    for (title, names) in GROUPS {
        let _ = writeln!(out, "Commands for {title}:");
        for name in *names {
            let about = if *name == "help" {
                // clap's auto-generated help subcommand's own about is long;
                // the group list keeps one line (doc_v14 02 §1.1).
                "Print this message or a subcommand's help".to_string()
            } else {
                let cmd = root.get_subcommands().find(|c| c.get_name() == *name);
                let about = cmd
                    .and_then(|c| c.get_about())
                    .map_or_else(String::new, |a| a.to_string());
                // The command list is the only place the aliases can be seen:
                // the group template drops clap's `{subcommands}` block, which
                // is what used to append them.
                let aliases = cmd
                    .map(|c| {
                        c.get_visible_aliases()
                            .map(str::to_string)
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default();
                if aliases.is_empty() {
                    about
                } else {
                    format!("{about} [aliases: {}]", aliases.join(", "))
                }
            };
            let _ = writeln!(out, "  {name:<12} {about}");
        }
    }
    out
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
    /// Remove these INFO tags (bcftools `-x`). Defaults to the stale list
    /// AC,AN,AF,MLEAC,MLEAF,ExcessHet,InbreedingCoeff,DP (the gold
    /// pipeline's `annotate -x`); pass an empty value (`-x ''`) to keep
    /// every tag
    #[arg(short = 'x', long = "remove-info", value_name = "LIST")]
    pub remove_info: Option<String>,
    /// Recompute and write these INFO tags (bcftools `+fill-tags`):
    /// AC, AN, NS, MAF or HWE (case-insensitive). Unsupported tags are an
    /// error, like bcftools
    #[arg(long, value_name = "LIST")]
    pub fill_tags: Option<String>,
    /// Output type: v (plain VCF), z (BGZF VCF), b (BCF) or u
    /// (uncompressed BCF); optional compression level (`-Oz6`, `-Ob1`)
    #[arg(short = 'O', long = "output-type", value_name = "v|z|b|u[0-9]")]
    pub output_type: Option<String>,
    /// Output file (default: stdout); the format comes from -O or the
    /// extension
    #[arg(
        short = 'o',
        long = "output",
        visible_alias = "out",
        value_name = "FILE"
    )]
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

/// Quality-control reports (doc_v14 02 §2.2): per-site and per-sample rows
/// (`--freq/--hardy/--het/--missing/--depth`), the bcftools-stats summary and
/// the histogram family. The windowed population-genetics statistics moved to
/// `reed popg` (doc_v14); the moved flags stay here as hidden parse-time
/// redirects so a 0.3.x command line fails with the exact rewrite.
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
    /// Groups file: `sample<TAB>population` rows, no header. Stratifies
    /// `--freq`/`--missing` (`PREFIX.afreq.pop`, `PREFIX.vmiss.pop`) and adds
    /// the POP column to `PREFIX.singletons`; the cohort files are unchanged
    /// (doc_v14 02 §2.2)
    #[arg(
        long = "pop",
        visible_alias = "populations",
        alias = "groups",
        value_name = "FILE"
    )]
    pub pop: Option<String>,
    /// Skip --pop rows naming samples absent from the header (they warn)
    #[arg(long)]
    pub force_samples: bool,
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
    /// Histograms: comma-separated built-ins (qual, dp, dp-site, gq, idl,
    /// af) and INFO/<TAG> or FMT/<TAG> requests; each writes PREFIX.<m>.hist
    #[arg(long, value_name = "LIST")]
    pub hist: Option<String>,
    /// Per-allele singletons and private doubletons (PREFIX.singletons, the
    /// vcftools layout)
    #[arg(long)]
    pub singletons: bool,
    /// Sites file (`CHROM<TAB>POS`, 1-based): only those sites are reported
    /// (the QC reports see the same gated records, doc_v7 D7-6)
    #[arg(long = "sites-file", value_name = "FILE")]
    pub sites_file: Option<String>,
    /// Output prefix: PREFIX.afreq, PREFIX.stats, ... (the reports write files)
    #[arg(
        short = 'o',
        long = "output",
        visible_alias = "out",
        value_name = "PREFIX"
    )]
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

    // ---- moved to `reed popg` (doc_v14 02 §2.4): hidden redirects ----
    // The fields keep their 0.3.x types so a 0.3.x command line parses and
    // `translate::stats` can refuse it with the exact rewrite before any I/O.
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub pi: bool,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub dxy: bool,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub theta: bool,
    /// Moved to `reed popg`.
    #[arg(long = "tajima-d", hide = true)]
    pub tajima_d: bool,
    /// Moved to `reed popg`.
    #[arg(long, hide = true, num_args = 0..=1, default_missing_value = "wc84")]
    pub fst: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub fst_global: bool,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub pbs: Vec<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub pbs_normed: bool,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub f2: Vec<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub f3: Vec<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub f4: Vec<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub dstat: Vec<String>,
    /// Moved to `reed popg`.
    #[arg(long = "f4-ratio", hide = true)]
    pub f4_ratio: Vec<String>,
    /// Moved to `reed popg`.
    #[arg(long = "block-size", hide = true)]
    pub block_size: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long = "ld-decay", hide = true, num_args = 0..=1, default_missing_value = "")]
    pub ld_decay: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long = "ld-pairs", hide = true)]
    pub ld_pairs: bool,
    /// Moved to `reed popg`.
    #[arg(long, hide = true, num_args = 0..=1, default_missing_value = "unfolded")]
    pub sfs: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub proj: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub sfs_int: bool,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub variant_only: bool,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub bypass_invariant_check: bool,
    /// Moved to `reed popg`.
    #[arg(long = "contig-lengths", alias = "contig-lengths-path", hide = true)]
    pub contig_lengths: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long = "window", alias = "size", hide = true)]
    pub size: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long = "bed-file", hide = true)]
    pub bed_file: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub step: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub levels: Option<String>,
    /// Moved to `reed popg`.
    #[arg(long, hide = true)]
    pub skip_empty: bool,
}

/// `reed popg`: the windowed population-genetics statistics (doc_v14 02 §2.1).
/// Report rows are per population (pair/triple) per window or genome-wide.
#[derive(Args)]
pub(crate) struct PopgCli {
    #[arg(help = INPUT_HELP)]
    pub input: String,
    #[command(flatten)]
    pub regions_args: filter_flags::RegionArgs,
    #[command(flatten)]
    pub sample_args: filter_flags::SampleArgs,
    #[command(flatten)]
    pub site_args: filter_flags::SiteFilterArgs,
    /// Withdrawn flag: it stated the fraction that must be called
    #[arg(long, value_name = "FLOAT", hide = true,
          value_parser = filter_flags::removed_max_missing_sample)]
    pub max_missing_sample: Option<String>,
    /// Groups file: `sample<TAB>population` rows, no header (the
    /// populations-file convention); required by the windowed modes
    #[arg(
        long = "pop",
        visible_alias = "populations",
        alias = "groups",
        value_name = "FILE"
    )]
    pub pop: Option<String>,
    /// π per population (Nei & Li 1979; missing data: Korunes & Samuk
    /// 2021)
    #[arg(long)]
    pub pi: bool,
    /// dxy per population pair (Nei 1987; missing data: Korunes & Samuk 2021)
    #[arg(long)]
    pub dxy: bool,
    /// Watterson's θ per population (Watterson 1975; per-allele-class
    /// mutation counts, Korunes & Samuk 2021)
    #[arg(long)]
    pub theta: bool,
    /// Tajima's D per population (Tajima 1989; denominator once at the
    /// window-average allele count, Korunes & Samuk 2021)
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
    /// LD decay: pairwise r² within MAX_BP (default 100k), 1 kb bins
    #[arg(long = "ld-decay", value_name = "MAX_BP", num_args = 0..=1,
          default_missing_value = "", value_parser = parse_optional_size)]
    pub ld_decay: Option<String>,
    /// Also write the per-pair LD table (PREFIX.ld)
    #[arg(long = "ld-pairs")]
    pub ld_pairs: bool,
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
    /// (default: error, the reference behavior)
    #[arg(long)]
    pub force_samples: bool,
    /// Contig lengths as a .fai-style file (name, length), for inputs
    /// whose header lacks ##contig lines
    #[arg(
        long = "contig-lengths",
        alias = "contig-lengths-path",
        value_name = "FILE"
    )]
    pub contig_lengths: Option<String>,
    /// Window size in bp (accepts 100k, 1M, 100kb)
    #[arg(long = "window", alias = "size", default_value = "100000", value_parser = parse_size,
          conflicts_with = "bed_file")]
    pub size: u64,
    /// BED3 file giving the windows (0-based half-open; mutually exclusive
    /// with --window/--step)
    #[arg(long = "bed-file", value_name = "BED")]
    pub bed_file: Option<String>,
    /// Sites file (`CHROM<TAB>POS`, 1-based): only those sites are counted
    #[arg(long = "sites-file", value_name = "FILE")]
    pub sites_file: Option<String>,
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
    #[arg(
        short = 'o',
        long = "output",
        visible_alias = "out",
        value_name = "PREFIX"
    )]
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

/// Conversion targets (doc_v13 §1): phylogenetic matrices (vcf2phylip
/// semantics), the vcftools 0/1/2 matrix, TASSEL hapmap, a reference-driven
/// MSA (vcf2msa semantics) and the projected-SFS export.
#[derive(Args)]
pub(crate) struct ConvertCli {
    #[arg(help = INPUT_HELP)]
    pub input: String,
    /// Projected SFS (mirrors `reed popg --sfs`): writes PREFIX_sfs/dadi
    /// and PREFIX_sfs/fastsimcoal2; needs --proj; `--sfs folded` folds by
    /// the minor allele
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
    /// Groups file for --sfs: `sample<TAB>population` rows, no header
    #[arg(
        long = "pop",
        visible_alias = "populations",
        alias = "groups",
        value_name = "FILE"
    )]
    pub pop: Option<String>,
    /// Skip group-file rows naming samples absent from the header
    #[arg(long)]
    pub force_samples: bool,
    /// PHYLIP matrix (relaxed layout; Ortiz 2019, vcf2phylip semantics)
    #[arg(long)]
    pub phylip: bool,
    /// FASTA matrix
    #[arg(long)]
    pub fasta: bool,
    /// NEXUS (DNA) matrix
    #[arg(long)]
    pub nexus: bool,
    /// Binary NEXUS for SNAPP (diploid biallelic sites only)
    #[arg(long = "nexus-binary")]
    pub nexus_binary: bool,
    /// Name of the sample written first in the alignments
    #[arg(long, value_name = "NAME")]
    pub outgroup: Option<String>,
    /// Minimum number of samples present at a locus (default 4 for the
    /// phylogenetic formats; --012/--hapmap apply the gate only when given)
    #[arg(long = "min-samples-locus", value_name = "N")]
    pub min_samples_locus: Option<u32>,
    /// Resolve ambiguous (heterozygous) genotypes to one allele; the draw is
    /// deterministic (seed 0), unlike vcf2phylip's
    #[arg(long = "resolve-iupac")]
    pub resolve_iupac: bool,
    /// Write the coordinates of the sites used in the alignments
    /// (PREFIX.used_sites.tsv)
    #[arg(long = "write-used-sites")]
    pub write_used_sites: bool,
    /// 0/1/2/-1 genotype matrix (vcftools `--012` layout):
    /// PREFIX.012, PREFIX.012.pos, PREFIX.012.indv
    #[arg(long = "012")]
    pub count_012: bool,
    /// TASSEL hapmap matrix (PREFIX.hmp.txt)
    #[arg(long)]
    pub hapmap: bool,
    /// Reference-driven whole-site MSA (vcf2msa semantics):
    /// PREFIX_<contig>.fasta per contig or region
    #[arg(long)]
    pub msa: bool,
    /// Reference genome FASTA (required by --msa; must match the VCF)
    #[arg(long, value_name = "FASTA")]
    pub reference: Option<String>,
    /// MSA region `CHR:B-E` (1-based inclusive; repeatable; default: whole
    /// contigs)
    #[arg(long, value_name = "CHR:B-E")]
    pub region: Vec<String>,
    /// Bases added on both sides of each --region
    #[arg(long, default_value_t = 0, value_parser = crate::units::parse_size)]
    pub flank: u64,
    /// Collapse sample names at the first occurrence of STR (--msa)
    #[arg(long = "collapse-suffix", value_name = "STR")]
    pub collapse_suffix: Option<String>,
    /// Convert a hapmap (`.hmp.txt`) input back to VCF (exclusive with every
    /// VCF-input target; doc_v13 §3.12)
    #[arg(long = "vcf")]
    pub to_vcf: bool,
    /// Output type for --vcf: v (plain VCF text) or z (BGZF), optional level
    /// (`-Oz6`)
    #[arg(short = 'O', long = "output-type", value_name = "v|z[0-9]")]
    pub output_type: Option<String>,
    /// Do not write the ##reed_* provenance header lines (--vcf)
    #[arg(long)]
    pub no_version: bool,
    /// Do not write the ##reed_command header line (--vcf)
    #[arg(long)]
    pub no_command: bool,
    #[command(flatten)]
    pub sample_args: filter_flags::SampleArgs,
    /// Output prefix: PREFIX.phy, PREFIX_sfs/, ... (default: the input's
    /// stem)
    #[arg(
        short = 'o',
        long = "output",
        visible_alias = "out",
        value_name = "PREFIX"
    )]
    pub out: Option<String>,
    /// Worker threads (default: all cores)
    #[arg(long)]
    pub threads: Option<usize>,
    /// In-flight memory budget in MiB (default: half of physical memory)
    #[arg(long, value_name = "MB")]
    pub memory: Option<u64>,
    /// Output level: 0 silent, 1 summary (default), 2 stage timing, 3 debug
    #[arg(long, value_name = "INT", default_value_t = Verbosity::Summary,
          value_parser = Verbosity::parse)]
    pub verbosity: Verbosity,
}

/// Identifier edits (doc_v14 02 §4.2): sample names, contig names, variant
/// IDs and header lines. The command never decodes genotypes (D14-4).
///
/// `-h` is bcftools `reheader -h`'s header replacement, so the automatic
/// `-h` help short is disabled here; `--help` stays.
#[derive(Args)]
#[command(disable_help_flag = true)]
pub(crate) struct RenameCli {
    #[arg(help = INPUT_HELP)]
    pub input: String,
    /// Replace the whole header (must keep the same samples)
    #[arg(short = 'h', long, value_name = "FILE")]
    pub header: Option<String>,
    /// Set `##contig` lengths from a `.fai` (name<TAB>length)
    #[arg(short = 'f', long, value_name = "FILE")]
    pub fai: Option<String>,
    /// Append a `##key=value` header line; repeatable
    #[arg(long = "add-meta", value_name = "LINE")]
    pub add_meta: Vec<String>,
    /// Remove every `##KEY=` header line; repeatable
    #[arg(long = "remove-meta", value_name = "KEY")]
    pub remove_meta: Vec<String>,
    /// Contig renames: `old<TAB>new` lines (the `##contig` lines and the
    /// CHROM column change together)
    #[arg(long = "rename-chrs", value_name = "FILE")]
    pub rename_chrs: Option<String>,
    /// Normalize contig names: plain (1, MT), ucsc (chr1, chrMT) or
    /// mt-is-m (chr1, chrM)
    #[arg(long = "chr-style", value_name = "STYLE")]
    pub chr_style: Option<String>,
    /// Rewrite every variant ID with a template: `%CHROM %POS %REF %ALT`
    #[arg(long = "set-id", value_name = "PATTERN")]
    pub set_id: Option<String>,
    /// Fill in IDs that are `.` with a template
    #[arg(long = "set-missing-var-ids", value_name = "PATTERN")]
    pub set_missing_var_ids: Option<String>,
    /// ID renames: `old<TAB>new` lines
    #[arg(long = "update-name", value_name = "FILE")]
    pub update_name: Option<String>,
    /// Replace `:` in IDs with SEP (GCTA/treemix-safe IDs)
    #[arg(long = "idspace-to", value_name = "SEP")]
    pub idspace_to: Option<String>,
    /// Suffix repeated IDs with `_1`, `_2`, … (never drop a record)
    #[arg(long = "rm-dup-ids")]
    pub rm_dup_ids: bool,
    /// Set every ID to `.`
    #[arg(long = "clear-ids")]
    pub clear_ids: bool,
    /// Sample renames as `old new,old2 new2`
    #[arg(short = 's', long, value_name = "PAIRS")]
    pub samples: Option<String>,
    /// Sample renames: `old<TAB>new` lines
    #[arg(short = 'N', long = "samples-file", value_name = "FILE")]
    pub samples_file: Option<String>,
    /// Rewrite matching sample names: a `*` prefix/suffix strips, else the
    /// text is removed wherever it occurs
    #[arg(long = "replace-sample-name", value_name = "PATTERN")]
    pub replace_sample_name: Option<String>,
    /// Output file (default: stdout); the format comes from -O or the
    /// extension
    #[arg(short = 'o', long = "output", visible_alias = "out", value_name = "FILE")]
    pub out: Option<String>,
    /// Output type: v (plain VCF) or z (BGZF VCF); optional compression
    /// level (`-Oz6`)
    #[arg(short = 'O', long = "output-type", value_name = "v|z[0-9]")]
    pub output_type: Option<String>,
    /// Do not write the ##reed_* provenance header lines
    #[arg(long)]
    pub no_version: bool,
    /// Do not write the ##reed_command line
    #[arg(long)]
    pub no_command: bool,
    /// Write an index next to the output: -W or -W=csi (default), or
    /// -W=tbi; omit for no index (bcftools semantics)
    #[arg(short = 'W', long = "write-index", value_name = "csi|tbi",
          num_args = 0..=1, default_missing_value = "csi")]
    pub write_index: Option<String>,
    /// Print help (the `-h` short belongs to --header here)
    #[arg(long = "help", action = clap::ArgAction::Help)]
    pub help: Option<bool>,
    /// Worker threads (default: all cores)
    #[arg(long)]
    pub threads: Option<usize>,
    /// In-flight memory budget in MiB (default: half of physical memory)
    #[arg(long, value_name = "MB")]
    pub memory: Option<u64>,
    /// Output level: 0 silent, 1 summary (default), 2 stage timing, 3 debug
    #[arg(long, value_name = "INT", default_value_t = Verbosity::Summary,
          value_parser = Verbosity::parse)]
    pub verbosity: Verbosity,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// doc_v14 02 §1.4 (guard 1/2): every `Command` variant appears in
    /// exactly one `GROUPS` entry, and the groups name nothing else. The
    /// automatic `help` subcommand is not a variant but is listed, so the
    /// comparison is against the variants plus `help`.
    #[test]
    fn groups_cover_every_command_exactly_once() {
        use clap::CommandFactory as _;
        // `build()` injects clap's automatic `help` subcommand; without it
        // the enumeration misses it and the guard would be vacuous for it.
        let mut root = Cli::command();
        root.build();
        let names: Vec<String> = root
            .get_subcommands()
            .map(|c| c.get_name().to_string())
            .collect();
        let mut listed: Vec<&str> = GROUPS
            .iter()
            .flat_map(|(_, members)| members.iter().copied())
            .collect();
        let total = listed.len();
        listed.sort_unstable();
        listed.dedup();
        assert_eq!(total, listed.len(), "a command is listed in two groups");
        for name in &names {
            assert!(
                listed.contains(&name.as_str()),
                "{name} is in no group"
            );
        }
        for name in &listed {
            assert!(
                names.iter().any(|n| n == name),
                "group lists unknown command {name}"
            );
        }
    }

    /// doc_v14 02 §1.4 (guard 13): `reed --help` renders the five group
    /// titles and every command's about text comes from clap.
    #[test]
    fn grouped_help_lists_the_five_groups() {
        let text = render_grouped_commands();
        for title in ["Filtering", "Statistics", "Transform", "File Operations", "Miscellaneous"] {
            assert!(text.contains(&format!("Commands for {title}:")), "{title} missing:\n{text}");
        }
        for name in ["filter", "stats", "popg", "convert", "index"] {
            assert!(text.contains(name), "{name} missing:\n{text}");
        }
    }


    /// A representative convert request renders to flags clap accepts, and
    /// the parsed flags rebuild the same request field by field (doc_v13
    /// §4.3: the `stats --sfs` / `convert --sfs` text contract).
    #[test]
    fn convert_args_round_trip() {
        use reed::{ConvertRequest, ConvertTarget, PhyloFormats, PhyloOptions, SfsRequest};
        let request = ConvertRequest {
            targets: vec![
                ConvertTarget::Sfs {
                    request: SfsRequest {
                        projections: vec![6],
                        folded: true,
                    },
                    int_dtype: true,
                },
                ConvertTarget::Phylo {
                    formats: PhyloFormats {
                        phylip: true,
                        fasta: false,
                        nexus: true,
                    },
                    options: PhyloOptions {
                        outgroup: Some("s3".into()),
                        resolve_iupac: true,
                        write_used_sites: true,
                    },
                },
                ConvertTarget::NexusBin,
                ConvertTarget::CountMatrix,
                ConvertTarget::Hapmap,
            ],
            prefix: std::path::PathBuf::from("out"),
            min_samples_locus: Some(7),
            ..ConvertRequest::default()
        };
        let mut argv: Vec<String> = vec!["reed".into(), "convert".into(), "in.vcf".into()];
        argv.extend(request.to_args().expect("to_args"));
        let cli = Cli::try_parse_from(argv).expect("to_args must parse");
        let Command::Convert(convert) = cli.command else {
            panic!("convert subcommand");
        };
        assert_eq!(convert.sfs.as_deref(), Some("folded"));
        assert_eq!(convert.proj.as_deref(), Some("6"));
        assert!(convert.sfs_int);
        assert!(convert.phylip && !convert.fasta && convert.nexus);
        assert_eq!(convert.outgroup.as_deref(), Some("s3"));
        assert!(convert.resolve_iupac);
        assert!(convert.write_used_sites);
        assert!(convert.nexus_binary);
        assert!(convert.count_012);
        assert!(convert.hapmap);
        assert!(!convert.msa);
        assert_eq!(convert.min_samples_locus, Some(7));
    }

    /// The `--vcf` reverse target renders to flags clap accepts, and the
    /// parsed flags rebuild it (doc_v13 §3.12: the hapmap→VCF text contract).
    #[test]
    fn convert_vcf_args_round_trip() {
        use reed::{ConvertRequest, ConvertTarget};
        let request = ConvertRequest {
            targets: vec![ConvertTarget::Vcf {
                compressed: true,
                level: 4,
                output: Some(std::path::PathBuf::from("back.vcf.gz")),
                command_line: "reed convert in.hmp.txt --vcf".into(),
                no_version: true,
            }],
            ..ConvertRequest::default()
        };
        let mut argv: Vec<String> = vec!["reed".into(), "convert".into(), "in.hmp.txt".into()];
        argv.extend(request.to_args().expect("to_args"));
        let cli = Cli::try_parse_from(argv).expect("to_args must parse");
        let Command::Convert(convert) = cli.command else {
            panic!("convert subcommand");
        };
        assert!(convert.to_vcf);
        assert_eq!(convert.output_type.as_deref(), Some("z4"));
        assert_eq!(convert.out.as_deref(), Some("back.vcf.gz"));
        assert!(convert.no_version);
        assert!(!convert.phylip && !convert.hapmap && !convert.msa);

        // Plain text renders `-O v` and no provenance switch.
        let plain = ConvertRequest {
            targets: vec![ConvertTarget::Vcf {
                compressed: false,
                level: 6,
                output: None,
                command_line: String::new(),
                no_version: false,
            }],
            ..ConvertRequest::default()
        };
        assert_eq!(plain.to_args().expect("to_args"), vec!["--vcf", "-O", "v"]);
    }

    /// `reed convert --sfs` and `reed popg --sfs` render the SFS portion
    /// identically (D13-6: one request, two entries).
    #[test]
    fn convert_and_stats_sfs_spellings_agree() {
        use reed::{ConvertRequest, ConvertTarget, SfsRequest};
        let sfs = SfsRequest {
            projections: vec![20, 20],
            folded: false,
        };
        let convert = ConvertRequest {
            targets: vec![ConvertTarget::Sfs {
                request: sfs.clone(),
                int_dtype: false,
            }],
            ..ConvertRequest::default()
        };
        let convert_args = convert.to_args().expect("to_args");
        assert_eq!(convert_args, vec!["--sfs", "unfolded", "--proj", "20,20"]);

        let mut argv: Vec<String> = vec![
            "reed".into(),
            "stats".into(),
            "in.vcf".into(),
            "--pop".into(),
            "pops.tsv".into(),
        ];
        argv.extend(convert_args);
        let cli = Cli::try_parse_from(argv).expect("the convert spelling must parse for stats");
        let Command::Stats(stats) = cli.command else {
            panic!("stats subcommand");
        };
        assert_eq!(stats.sfs.as_deref(), Some("unfolded"));
        assert_eq!(stats.proj.as_deref(), Some("20,20"));
    }

    /// `WriteIndex::to_args` renders flags `reed index` accepts, and the
    /// parsed `IndexCli` maps back to the same flavor (doc_v5 08 §6.2 text
    /// contract; the flavor mapping is `translate/index.rs`'s rule).
    #[test]
    fn index_args_round_trip() {
        for (flavor, args) in [
            (reed::WriteIndex::Csi, vec!["--csi"]),
            (reed::WriteIndex::Tbi, vec!["--tbi"]),
        ] {
            let mut argv = vec!["reed", "index", "in.vcf.gz"];
            argv.extend(args);
            let cli = Cli::try_parse_from(argv).expect("to_args must parse");
            let Command::Index(index) = cli.command else {
                panic!("index subcommand");
            };
            let rebuilt = if index.tbi {
                reed::WriteIndex::Tbi
            } else {
                reed::WriteIndex::Csi
            };
            assert_eq!(rebuilt, flavor);
        }
    }
}
