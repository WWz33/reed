//! Input manifest: what sites the file actually contains.
//!
//! The manifest is a one-time declaration about the data, not a per-call knob.
//! Estimators derive their availability from it; the denominator convention of
//! each estimator is sealed into the estimator definition itself.

/// How the input came to contain the sites it contains.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ascertainment {
    /// Every variant call was retained; no site was dropped by filtering.
    Complete,
    /// True variants were removed by filtering (e.g. GATK `SNP_FILTER`). In the
    /// file, removed true variants are indistinguishable from invariant sites.
    Filtered,
    /// The caller has not established whether filtering removed sites.
    Unknown,
}

/// Whether the input explicitly represents monomorphic sites.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SiteCompleteness {
    /// All-sites input: monomorphic records are present (all-sites VCF or an
    /// expanded gVCF). Per-pair denominators are exact.
    WithInvariants,
    /// Only variant sites are present in the file.
    VariantOnly(Ascertainment),
}

/// What an estimator can do under a given manifest.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Availability {
    /// Computable, no known bias from the manifest.
    Available,
    /// Computable, but the estimate carries a documented, unresolvable bias.
    Biased(&'static str),
    /// Refused: the manifest cannot support this estimator's denominator.
    Unavailable(&'static str),
}
