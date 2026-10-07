//! Ascertainment v2 derivation and write-back (§11.7): the filter declares
//! what it did to the data, in the output header, so the statistics layer
//! and the filter agree by construction.

use reed_core::manifest::v2::{Ascertainment, CohortState, GenotypeState, InvariantSet, SiteSet};

/// Derives the ascertainment from the compiled spec (§11.7 table). The
/// cohort dimension is `Full` here and completed at runtime once the keep
/// decision exists.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn derive_ascertainment(
    has_site_rules: bool,
    variant_class: &[String],
    exclude_variant_class: &[String],
    min_maf: Option<f64>,
    max_maf: Option<f64>,
    min_mac: u32,
    max_mac: u32,
    mask_rules: &[crate::plan::MaskRule],
    plan_hash: &str,
) -> Ascertainment {
    let sites = if has_site_rules {
        SiteSet::Filtered { plan_hash: plan_hash.to_string() }
    } else {
        SiteSet::Complete
    };
    let invariants = if variant_class.iter().any(|c| c != "ref")
        || exclude_variant_class.iter().any(|c| c != "ref")
        || min_maf.is_some_and(|m| m > 0.0)
        || max_maf.is_some_and(|m| m < 0.5)
        || min_mac > 0
        || max_mac > 0
    {
        InvariantSet::Absent
    } else {
        InvariantSet::Present
    };
    let genotypes = if mask_rules.is_empty() {
        GenotypeState::Original
    } else {
        GenotypeState::Masked { plan_hash: plan_hash.to_string() }
    };
    Ascertainment { sites, invariants, genotypes, cohort: CohortState::Full }
}

/// Completes the cohort dimension once the keep decision exists (§11.7:
/// "cohort 规则实际删除了样本").
#[must_use]
pub fn with_cohort_subset(
    mut ascertainment: Ascertainment,
    keep_hash: &str,
    n_samples: u32,
) -> Ascertainment {
    ascertainment.cohort = CohortState::Subset { keep_hash: keep_hash.to_string(), n_samples };
    ascertainment
}
