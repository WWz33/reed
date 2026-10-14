//! The filter specification (§11.1): typed, serde-serializable, validated
//! at plan-compile time against the input header.
//!
//! One source feeds the spec: the CLI flags (02-cli.md §3), translated
//! into these types by `reed-cli`; tests construct it directly. There is no
//! TOML or conf file format (02-cli.md §0).

pub mod cohort;
pub mod mask;
pub mod missing_rate;
pub mod predicate;
pub mod recompute;
pub mod region;
pub mod render;

use serde::{Deserialize, Serialize};

use crate::FilterError;

pub use cohort::{SampleList, SampleSelect};
pub use mask::GenotypeMaskRule;
pub use missing_rate::MissingRate;
pub use predicate::{Cmp, OnMissing, Quantifier, Typed};
pub use recompute::RecomputeSpec;
pub use region::{Region, RegionRule};

/// The filter specification (§11.1). One typed value; the flag vector is
/// its other spelling (`to_args`, `spec/render.rs`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct FilterSpec {
    /// Static sample subset (`-s/-S/--keep/--remove`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_select: Option<SampleSelect>,
    /// Dynamic sample rule: keep samples whose **missing rate** is at most
    /// this (`--mind`, doc_v4 02-cli-semantics.md 回退版 §3.2). Computed
    /// after masking, over the alleles of the sites that pass the head
    /// rules; applied before the tail rules.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mind: Option<MissingRate>,
    /// Streaming region predicates (bcftools `-t` semantics).
    #[serde(default)]
    pub targets: RegionRule,
    /// Site-level rules evaluated on the fixed columns + INFO.
    #[serde(default)]
    pub site: SiteRules,
    /// Genotype-level mask rules (any hit → GT becomes `./.`).
    #[serde(default)]
    pub genotype: GenotypeRules,
    /// Site-level rules recomputed after masking.
    #[serde(default)]
    pub site_after_mask: TailRules,
    /// INFO tags recomputed and stale tags dropped.
    #[serde(default)]
    pub recompute: RecomputeSpec,
    /// Ploidy / multiallelic policy.
    #[serde(default)]
    pub policy: PolicyRules,
    /// `--indep-pairwise W S R2`: LD pruning as a site rule
    /// (docs/rounds/v5/05-algorithms.md §3.4.1). `None` = no pruning.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ld_prune: Option<LdPruneSpec>,
}

/// Which allele frequency a `-q/-Q` threshold applies to (bcftools
/// `:type` suffix; docs/rounds/v5/02-cli.md §3.2). The default is `nref`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum AfKind {
    /// Frequency of non-reference alleles (bcftools' default).
    Nref,
    /// Frequency of the first ALT allele (bcftools `:alt1`).
    Alt1,
    /// Frequency of the minor allele (bcftools `:minor`).
    Minor,
    /// Frequency of the major allele (bcftools `:major`).
    Major,
    /// `1 − minor` (bcftools `:nonmajor`).
    Nonmajor,
}

impl AfKind {
    /// Parses the `:type` suffix (without the colon).
    ///
    /// # Errors
    ///
    /// [`FilterError::Spec`] listing the supported tokens for an unknown one.
    pub fn parse(text: &str) -> Result<Self, FilterError> {
        Ok(match text {
            "nref" => Self::Nref,
            "alt1" => Self::Alt1,
            "minor" => Self::Minor,
            "major" => Self::Major,
            "nonmajor" => Self::Nonmajor,
            other => {
                return Err(FilterError::Spec(format!(
                    "unknown allele-frequency type {other:?} \
                     (expected nref, alt1, minor, major or nonmajor)"
                )));
            }
        })
    }

    /// The stable suffix spelling.
    #[must_use]
    pub const fn id(self) -> &'static str {
        match self {
            Self::Nref => "nref",
            Self::Alt1 => "alt1",
            Self::Minor => "minor",
            Self::Major => "major",
            Self::Nonmajor => "nonmajor",
        }
    }
}

/// The `--indep-pairwise` rule (docs/rounds/v5/05-algorithms.md §3.4.1).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct LdPruneSpec {
    /// Window size: `W` variants, or `W` bp when `kb`.
    pub window: u64,
    /// Whether the window is in kb (plink2 `W kb`).
    pub kb: bool,
    /// Step in variants (1 when the window is in kb).
    pub step: usize,
    /// r² threshold.
    pub r2: f64,
}

impl FilterSpec {
    /// Whether multiallelic inputs may survive into the output
    /// (`docs/spec/architecture.md` §11.6: `multiallelic != "reject"` or
    /// `max_alt_alleles > 1`). When false, the filter's plan gates records
    /// with more than one ALT ([`crate::plan::HeadRule::BiallelicGate`]),
    /// and the recomputed `AC`/`MAF`/`HWE` values stay per-record valid.
    #[must_use]
    pub fn multiallelics_allowed(&self) -> bool {
        self.policy.multiallelic == MultiallelicPolicy::Allow || self.site.max_alt_alleles > 1
    }
}

/// Site-level rules (§11.1 `[site]`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SiteRules {
    /// Keep only these variant classes (`snv | mnv | indel | symbolic | ref`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub variant_class: Vec<String>,
    /// Drop these variant classes (bcftools `-V/--exclude-types`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exclude_variant_class: Vec<String>,
    /// `-i` expressions that do not read genotype-derived quantities
    /// (02-cli.md §4): the record must satisfy every one of them.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub include: Vec<String>,
    /// `-e` expressions that do not read genotype-derived quantities: a
    /// record satisfying any of them is dropped.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exclude: Vec<String>,
    /// Minimum ALT allele count (`0` = off) — bcftools' `-m` counts all
    /// alleles, so `-m 2` becomes `min_alt_alleles = 1`.
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub min_alt_alleles: u32,
    /// Maximum ALT allele count (`0`/absent = unlimited).
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub max_alt_alleles: u32,
    /// FILTER membership (`.` never passes).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub require_filter: Vec<String>,
    /// QUAL lower bound.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_qual: Option<f32>,
    /// Arbitrary numeric INFO predicates, evaluated in declaration order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub info: Vec<InfoPredicate>,
}

/// One numeric INFO predicate (`[[site.info]]`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct InfoPredicate {
    /// The INFO tag (validated against the header at compile time).
    pub tag: String,
    /// The comparison, e.g. `">= 4.0"`.
    pub keep_if: String,
    /// Missing-tag behaviour (default `pass`, bcftools semantics).
    #[serde(default)]
    pub on_missing: OnMissing,
    /// Required when the header declares `Number=A` and multiallelic
    /// records are possible: `any` | `all`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub quantifier: Option<Quantifier>,
}

/// Genotype mask rules (`[genotype]`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct GenotypeRules {
    /// Any hit masks the sample's GT to `./.` (evaluated in order).
    #[serde(default, rename = "mask_if", skip_serializing_if = "Vec::is_empty")]
    pub mask_if: Vec<GenotypeMaskRule>,
    /// `--set-gts .` (bcftools `filter -S .`): mask from the `FMT/`
    /// expressions in `site_after_mask` instead of dropping the site.
    /// `Some(".")` is the only accepted value.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub set_gts: Option<String>,
}

/// Post-mask site rules (`[site_after_mask]`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct TailRules {
    /// `-i` expressions that read genotype-derived quantities (evaluated
    /// after the mask stage, 02-cli.md §4).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub include: Vec<String>,
    /// `-e` expressions that read genotype-derived quantities.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exclude: Vec<String>,
    /// Site missing-rate cap (`--geno`: keep sites whose missing-allele
    /// fraction is at most this; denominator: kept samples' allele slots,
    /// numerator: missing alleles — vcftools/plink2 semantics, doc_v4 02
    /// 回退版 §3.1).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_missing: Option<MissingRate>,
    /// MAF lower bound (recomputed from GT).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_maf: Option<f64>,
    /// MAF upper bound (plink2 `--max-maf`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_maf: Option<f64>,
    /// Minor-allele-count lower bound (`0` = off).
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub min_mac: u32,
    /// Minor-allele-count upper bound (`0` = off).
    #[serde(default, skip_serializing_if = "is_zero_u32")]
    pub max_mac: u32,
    /// Non-reference allele frequency lower bound (bcftools `-q/--min-af`);
    /// the INFO AC/AN values are used when both are present.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_af: Option<f64>,
    /// `-q`'s `:type` suffix (bcftools `:nref/:alt1/:minor/:major/:nonmajor`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_af_kind: Option<AfKind>,
    /// Non-reference allele frequency upper bound (bcftools `-Q/--max-af`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_af: Option<f64>,
    /// `-Q`'s `:type` suffix.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_af_kind: Option<AfKind>,
    /// HWE exact-test p lower bound (`0`/absent = off).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_hwe_p: Option<f64>,
    /// HWE test flavour (`exact` only for M5; `midp` is M7).
    #[serde(default)]
    pub hwe_mode: HweMode,
    /// Genotype-class rules (`-g/--genotype [^](hom|het|miss)`): the record
    /// must satisfy **every** listed rule (repeated flags AND, a deliberate
    /// superset of bcftools, which keeps the last `-g` only). Compiled after
    /// masking, over the kept samples; a half-call (`0/.`) reads as missing
    /// (bcftools `bcf_gt_type` UNKN). Haploid calls follow the ploidy
    /// policy: rejected by default, read as missing (half-call) under
    /// `allow_haploid` — bcftools instead reads them as hom.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub gt_class: Vec<GtClassRule>,
}

/// One genotype-class site rule (bcftools `-g/--genotype`; `TailRules::gt_class`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum GtClassRule {
    /// `hom`: at least one homozygous call.
    #[serde(rename = "hom")]
    Hom,
    /// `het`: at least one heterozygous call.
    #[serde(rename = "het")]
    Het,
    /// `miss`: at least one missing genotype — `./.`, `.`, **or a
    /// half-call** (`0/.`; bcftools `bcf_gt_type` UNKN).
    #[serde(rename = "miss")]
    Miss,
    /// `^hom`: no homozygous call.
    #[serde(rename = "^hom")]
    NoHom,
    /// `^het`: no heterozygous call.
    #[serde(rename = "^het")]
    NoHet,
    /// `^miss`: no missing genotype (half-calls included).
    #[serde(rename = "^miss")]
    NoMiss,
}

impl GtClassRule {
    /// Parses the CLI spelling: `[^](hom|het|miss)`, case-insensitive.
    ///
    /// # Errors
    ///
    /// [`FilterError::Spec`] listing the supported tokens for an unknown one.
    pub fn parse(text: &str) -> Result<Self, FilterError> {
        let (exclude, name) = match text.strip_prefix('^') {
            Some(rest) => (true, rest),
            None => (false, text),
        };
        // Exact spelling, case-insensitive only (bcftools `strcasecmp`;
        // no whitespace is tolerated).
        let class = name.to_ascii_lowercase();
        Ok(match (exclude, class.as_str()) {
            (false, "hom") => Self::Hom,
            (false, "het") => Self::Het,
            (false, "miss") => Self::Miss,
            (true, "hom") => Self::NoHom,
            (true, "het") => Self::NoHet,
            (true, "miss") => Self::NoMiss,
            _ => {
                return Err(FilterError::Spec(format!(
                    "unknown genotype class {text:?} \
                     (expected hom, het, miss, ^hom, ^het or ^miss)"
                )));
            }
        })
    }

    /// The canonical CLI spelling.
    #[must_use]
    pub const fn id(self) -> &'static str {
        match self {
            Self::Hom => "hom",
            Self::Het => "het",
            Self::Miss => "miss",
            Self::NoHom => "^hom",
            Self::NoHet => "^het",
            Self::NoMiss => "^miss",
        }
    }

    /// The class this rule tests (`hom`, `het` or `miss`).
    #[must_use]
    pub const fn class(self) -> &'static str {
        match self {
            Self::Hom | Self::NoHom => "hom",
            Self::Het | Self::NoHet => "het",
            Self::Miss | Self::NoMiss => "miss",
        }
    }
}

/// HWE test flavour.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum HweMode {
    /// Wigginton exact test (the default).
    #[default]
    #[serde(rename = "exact")]
    Exact,
    /// Mid-p adjustment (M7).
    #[serde(rename = "midp")]
    MidP,
}

/// Ploidy / multiallelic policy (`[policy]`).
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolicyRules {
    /// `require_diploid` (default) rejects haploid/polyploid calls.
    #[serde(default)]
    pub ploidy: PloidyPolicy,
    /// `reject` (default) rejects multiallelic records.
    #[serde(default)]
    pub multiallelic: MultiallelicPolicy,
}

/// Ploidy policy.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum PloidyPolicy {
    /// Haploid or polyploid calls are a per-record error.
    #[default]
    #[serde(rename = "require_diploid")]
    RequireDiploid,
    /// Haploid calls are tolerated (counted as half-calls).
    #[serde(rename = "allow_haploid")]
    AllowHaploid,
}

/// Multiallelic policy.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum MultiallelicPolicy {
    /// Multiallelic records are a per-record error.
    #[default]
    #[serde(rename = "reject")]
    Reject,
    /// Multiallelic records pass the biallelic gate (recompute handles
    /// `Number=A` semantics, §11.6).
    #[serde(rename = "allow")]
    Allow,
}

fn is_zero_u32(v: &u32) -> bool {
    *v == 0
}

#[cfg(test)]
mod tests {
    use super::GtClassRule as R;

    /// `-g/--genotype` spellings: case-insensitive class names, `^`
    /// excludes; `id()` is the canonical spelling the CLI renders back.
    #[test]
    fn gt_class_rule_parse_and_render() {
        for (text, rule) in [
            ("hom", R::Hom),
            ("Het", R::Het),
            ("MISS", R::Miss),
            ("^hom", R::NoHom),
            ("^Het", R::NoHet),
            ("^MISS", R::NoMiss),
        ] {
            assert_eq!(R::parse(text).expect("valid"), rule, "{text}");
            // `id()` round-trips and its class part agrees with `class()`.
            assert_eq!(R::parse(rule.id()).expect("canonical"), rule);
            assert_eq!(rule.id().trim_start_matches('^'), rule.class());
        }
        // No whitespace tolerance (bcftools `strcasecmp` is exact).
        for bad in ["", "^", "homo", "hom,het", "^hom^", " HOM ", "hom "] {
            assert!(R::parse(bad).is_err(), "{bad:?}");
        }
    }
}
