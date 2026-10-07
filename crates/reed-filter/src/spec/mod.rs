//! The filter specification (§11.1): typed, serde-serializable, validated
//! at plan-compile time against the input header.
//!
//! One source feeds the v2 spec: the CLI flags (02-cli.md §3), translated
//! into these types by `reed-cli`; tests construct it directly. There is no
//! TOML or conf file format (02-cli.md §0).

pub mod cohort;
pub mod mask;
pub mod predicate;
pub mod recompute;
pub mod region;

use serde::{Deserialize, Serialize};

pub use cohort::CohortRule;
pub use mask::GenotypeMaskRule;
pub use predicate::{Cmp, OnMissing, Quantifier, Typed};
pub use region::{Region, RegionRule};
pub use recompute::RecomputeSpec;

/// The v2 filter specification (§11.1). Field names mirror the TOML keys.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct FilterSpecV2 {
    /// Sample-level rule (the pass0 global reduction).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cohort: Option<CohortRule>,
    /// Streaming region predicates (bcftools `-t` semantics).
    #[serde(default)]
    pub targets: RegionRule,
    /// Site-level rules evaluated on the fixed columns + INFO.
    #[serde(default)]
    pub site: SiteRulesV2,
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
}

/// Site-level rules (§11.1 `[site]`).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SiteRulesV2 {
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
    /// Site missing-rate cap (denominator: kept samples).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_missing: Option<f64>,
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
    /// HWE exact-test p lower bound (`0`/absent = off).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_hwe_p: Option<f64>,
    /// HWE test flavour (`exact` only for M5; `midp` is M7).
    #[serde(default)]
    pub hwe_mode: HweMode,
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


