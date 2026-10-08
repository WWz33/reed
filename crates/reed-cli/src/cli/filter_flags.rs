//! Shared filtering flag groups (04-cleanup §2.1): `reed filter` and
//! `reed stats` accept the same region and sample selection spellings, so
//! the definitions live here once.

use clap::Args;

/// Region selection (`doc_v2/02-cli.md` §2): indexed and streamed regions.
#[derive(Args, Clone, Debug, Default)]
pub(crate) struct RegionArgs {
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
    /// Regions to keep, one per line (bcftools `-T`); `^` prefix excludes
    #[arg(short = 'T', long = "targets-file", value_name = "FILE")]
    pub targets_file: Option<String>,
}

/// Sample selection (`doc_v2/02-cli.md` §3): keep and remove lists.
#[derive(Args, Clone, Debug, Default)]
pub(crate) struct SampleArgs {
    /// Samples to keep: comma-separated names; `^` prefix excludes
    /// (`^A,B` = keep all but A and B; bcftools `-s`)
    #[arg(short = 's', long = "samples", value_name = "[^]LIST")]
    pub samples: Vec<String>,
    /// Samples to keep, one name per line (bcftools `-S`)
    #[arg(short = 'S', long = "samples-file", value_name = "FILE")]
    pub samples_file: Option<String>,
    /// Samples to keep, one name per line (plink2/vcftools spelling)
    #[arg(long, value_name = "FILE")]
    pub keep: Option<String>,
    /// Samples to drop, one name per line
    #[arg(long, value_name = "FILE")]
    pub remove: Option<String>,
}

/// Site-level selection and genotype thresholds (04-cleanup §2.1): the
/// bcftools `view`/`filter` spellings both `reed filter` and `reed stats`
/// accept. INFO *updates* (`-x`, `--fill-tags`, `-I`) stay filter-only.
#[derive(Args, Clone, Debug, Default)]
pub(crate) struct SiteFilterArgs {
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
    /// Mask genotypes with FMT/DP below INT (missing DP never masks);
    /// plink2 spells this --vcf-min-dp
    #[arg(long, visible_alias = "vcf-min-dp", value_name = "INT")]
    pub min_dp: Option<u32>,
    /// Mask genotypes with FMT/DP above INT (plink2 --vcf-max-dp)
    #[arg(long, visible_alias = "vcf-max-dp", value_name = "INT")]
    pub max_dp: Option<u32>,
    /// Mask genotypes with FMT/GQ below INT (plink2 --vcf-min-gq)
    #[arg(long, visible_alias = "vcf-min-gq", value_name = "INT")]
    pub min_gq: Option<u32>,
    /// Mask from the -i/-e FORMAT expressions: set the genotype to this
    /// value (only `.`, bcftools filter -S .)
    #[arg(long, value_name = ".")]
    pub set_gts: Option<String>,
    /// Minimum non-reference allele frequency (bcftools -q): the INFO AC/AN
    /// values are used when both are present, else the genotype counts
    #[arg(short = 'q', long = "min-af", value_name = "FLOAT")]
    pub min_af: Option<f64>,
    /// Maximum non-reference allele frequency (bcftools -Q)
    #[arg(short = 'Q', long = "max-af", value_name = "FLOAT")]
    pub max_af: Option<f64>,
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
    /// HWE exact-test p lower bound (`P`, optional `midp`)
    #[arg(long, value_name = "P [midp]", num_args = 1..=2)]
    pub hwe: Vec<String>,
}
