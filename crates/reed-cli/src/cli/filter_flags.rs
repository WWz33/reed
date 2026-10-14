//! Shared filtering flag groups (04-cleanup §2.1): `reed filter` and
//! `reed stats` accept the same region and sample selection spellings, so
//! the definitions live here once.

use clap::Args;
use reed::spec::MissingRate;

/// The error text for a withdrawn keep-direction flag (doc_v4 02 回退版
/// §5.2): the replacement names the fraction that **may be missing**, `1 − X`
/// in exact decimal (`--max-missing 0.9` → `--geno 0.1`, not
/// `0.09999999999999998`). There is no silent alias — the message is the
/// only place the withdrawn spelling appears.
fn removed_keep_direction(flag: &str, replacement: &str, text: &str) -> Result<String, String> {
    let rate = MissingRate::parse(text).map_err(|e| e.to_string())?;
    Err(format!(
        "{flag} was removed: it stated the fraction that must be called; use {replacement} {} \
         (the fraction that may be missing)",
        rate.complement().render()
    ))
}

/// `--max-missing X` → `--geno 1−X`.
pub(crate) fn removed_max_missing(text: &str) -> Result<String, String> {
    removed_keep_direction("--max-missing", "--geno", text)
}

/// `--max-missing-sample X` → `--mind 1−X`.
pub(crate) fn removed_max_missing_sample(text: &str) -> Result<String, String> {
    removed_keep_direction("--max-missing-sample", "--mind", text)
}

/// `--max-missing-count C` has no counterpart: the family is two parameters
/// (doc_v4 02 回退版 §3.1), and a count cap is expressible as a fraction.
pub(crate) fn removed_max_missing_count(_text: &str) -> Result<String, String> {
    Err(
        "--max-missing-count was removed: the missing-threshold family is --geno/--mind only; \
         express a count cap C as the fraction C/2N for --geno (N is the retained sample count)"
            .to_string(),
    )
}

/// `reed stats` does not take the sample-level rule (`--mind`): it needs a
/// full pass before the statistics (doc_v4 02 回退版 §3.2).
/// Parses `FLOAT[:type]` for `-q/-Q` (bcftools `:nref/:alt1/:minor/:major/
/// `:nonmajor`; docs/rounds/v5/02-cli.md §3.2). The type defaults to `nref`.
pub(crate) fn parse_af_with_type(text: &str) -> Result<(f64, reed::spec::AfKind), String> {
    let (value, kind) = match text.split_once(':') {
        None => (text, reed::spec::AfKind::Nref),
        Some((v, t)) => (v, reed::spec::AfKind::parse(t).map_err(|e| e.to_string())?),
    };
    let value: f64 = value
        .trim()
        .parse()
        .map_err(|_| format!("invalid allele-frequency value {value:?}"))?;
    if !(0.0..=1.0).contains(&value) {
        return Err(format!("allele frequency out of range [0,1]: {value}"));
    }
    Ok((value, kind))
}

pub(crate) fn unsupported_stats_mind(_text: &str) -> Result<String, String> {
    Err(
        "reed stats does not support --mind (the sample rule needs a full pass before the \
         statistics); run `reed filter --mind X` first and feed its output to reed stats"
            .to_string(),
    )
}

/// Region selection (`docs/rounds/v2/02-cli.md` §2): indexed and streamed regions.
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

/// Sample selection (`docs/rounds/v2/02-cli.md` §3): keep and remove lists.
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
    /// Minimum allele frequency (bcftools -q): the INFO AC/AN values are
    /// used when both are present, else the genotype counts; an optional
    /// `:type` suffix selects the allele (`nref` default, `alt1`, `minor`,
    /// `major`, `nonmajor`)
    #[arg(short = 'q', long = "min-af", value_name = "FLOAT[:type]")]
    pub min_af: Option<String>,
    /// Maximum allele frequency (bcftools -Q), with the same `:type`
    /// suffix as `-q`
    #[arg(short = 'Q', long = "max-af", value_name = "FLOAT[:type]")]
    pub max_af: Option<String>,
    /// Keep sites whose missing rate (missing alleles / 2·samples, after
    /// masking) is at most FLOAT (plink2 semantics; vcftools spells the
    /// complement)
    #[arg(long, value_name = "FLOAT", value_parser = MissingRate::parse)]
    pub geno: Option<MissingRate>,
    /// Withdrawn flag: it stated the fraction that must be called; the
    /// replacement names the fraction that may be missing
    #[arg(long, value_name = "FLOAT", hide = true, value_parser = removed_max_missing)]
    pub max_missing: Option<String>,
    /// Withdrawn flag: the missing-threshold family is `--geno`/`--mind`
    #[arg(long, value_name = "INT", hide = true, value_parser = removed_max_missing_count)]
    pub max_missing_count: Option<String>,
    /// Minimum minor allele frequency (after masking)
    #[arg(long, value_name = "FLOAT")]
    pub maf: Option<f64>,
    /// Maximum minor allele frequency (after masking)
    #[arg(long, value_name = "FLOAT")]
    pub max_maf: Option<f64>,
    /// Minimum minor allele count (after masking)
    #[arg(long, value_name = "INT")]
    pub mac: Option<u32>,
    /// LD pruning (plink2 algorithm): `W[kb] [S] R2`, e.g. `50 5 0.2` or
    /// `50kb 0.2`; the kept sites are the ones the command emits
    #[arg(long = "indep-pairwise", value_name = "W[kb] [S] R2", num_args = 2..=3)]
    pub indep_pairwise: Vec<String>,
    /// Maximum minor allele count (after masking)
    #[arg(long, value_name = "INT")]
    pub max_mac: Option<u32>,
    /// HWE exact-test p lower bound (`P`, optional `midp`)
    #[arg(long, value_name = "P [midp]", num_args = 1..=2)]
    pub hwe: Vec<String>,
    /// Keep sites with/without a genotype class: [^](hom|het|miss); a
    /// half-call counts as missing (bcftools -g; repeated flags AND —
    /// bcftools keeps only the last)
    #[arg(short = 'g', long = "genotype", value_name = "[^](hom|het|miss)")]
    pub genotype: Vec<String>,
}
