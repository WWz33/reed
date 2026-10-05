//! Input manifest: what sites the file actually contains.
//!
//! The manifest is a one-time declaration about the data, not a per-call knob.
//! Estimators derive their availability from it; the denominator convention of
//! each estimator is sealed into the estimator definition itself.
//!
//! The v1 types below are what the legacy estimators consume; [`v2`] is the
//! target model (`docs/architecture-spec.md` §5.8) with a `From` compat
//! mapping. When the legacy estimators retire in reed-stats, the v1 types
//! move under `legacy` and v2 promotes to this level.

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


/// Ascertainment v2 (`docs/architecture-spec.md` §5.8): the four-dimension
/// manifest the target architecture consumes. Lives in a submodule while the
/// legacy estimators still own the v1 names at this level; promotes to the
/// crate root namespace when they retire (M4).
pub mod v2 {
/// Ascertainment v2: the four dimensions of what the input actually is
/// (`docs/architecture-spec.md` §5.8). Machine-readable (JSON in output
/// headers and sidecars).
#[derive(Clone, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Ascertainment {
    /// Which sites the input retains.
    pub sites: SiteSet,
    /// Whether invariant (monomorphic) sites are present.
    pub invariants: InvariantSet,
    /// Whether genotypes were masked after calling.
    pub genotypes: GenotypeState,
    /// Whether the cohort is the full sample set.
    pub cohort: CohortState,
}

/// Site dimension of [`Ascertainment`].
#[derive(Clone, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum SiteSet {
    /// Every variant call retained.
    #[default]
    Complete,
    /// Sites removed by filtering, identified by the filter plan hash.
    Filtered {
        /// sha256 of the canonical filter plan.
        plan_hash: String,
    },
    /// Nothing established.
    Unknown,
}

/// Invariant-site dimension of [`Ascertainment`].
#[derive(Clone, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum InvariantSet {
    /// Monomorphic records present.
    #[default]
    Present,
    /// Some invariant sites removed, identified by the filter plan hash.
    PartiallyRemoved {
        /// sha256 of the canonical filter plan.
        plan_hash: String,
    },
    /// No invariant sites in the input.
    Absent,
}

/// Genotype dimension of [`Ascertainment`].
#[derive(Clone, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum GenotypeState {
    /// Calls as the caller produced them.
    #[default]
    Original,
    /// Genotypes masked to missing, identified by the filter plan hash.
    Masked {
        /// sha256 of the canonical filter plan.
        plan_hash: String,
    },
}

/// Cohort dimension of [`Ascertainment`].
#[derive(Clone, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum CohortState {
    /// All header samples present.
    #[default]
    Full,
    /// A subset retained, identified by the keep mask hash.
    Subset {
        /// sha256 identifying the keep decision.
        keep_hash: String,
        /// Retained sample count.
        n_samples: u32,
    },
}

/// Availability under the v2 manifest: machine-readable bias and refusal
/// codes replace the v1 free-text payloads (fixing legacy issue A23).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Availability {
    /// Computable without manifest-induced bias.
    Available,
    /// Computable with a documented bias.
    Biased(BiasCode),
    /// Refused.
    Unavailable(RefusalCode),
}

/// Machine-readable bias codes (why an estimate is biased, not broken).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum BiasCode {
    /// Denominator counts only variant sites while the input also carries
    /// invariant sites the estimator would need.
    VariantOnlyDenominator,
    /// The numerator lost sites to filtering.
    FilteredNumerator,
    /// Genotypes were masked after calling.
    MaskedGenotypes,
}

impl BiasCode {
    /// Human-readable explanation (also rendered in TSV `bias` columns).
    #[must_use]
    pub const fn message(self) -> &'static str {
        match self {
            Self::VariantOnlyDenominator => {
                "denominator counts variant sites only; input also carries invariant sites"
            }
            Self::FilteredNumerator => "some sites were removed by filtering",
            Self::MaskedGenotypes => "genotypes were masked to missing after calling",
        }
    }

    /// Stable token for JSON sidecars.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::VariantOnlyDenominator => "variant_only_denominator",
            Self::FilteredNumerator => "filtered_numerator",
            Self::MaskedGenotypes => "masked_genotypes",
        }
    }
}

/// Machine-readable refusal codes (why an estimator cannot run at all).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum RefusalCode {
    /// The input lacks the invariant sites the estimator's denominator needs.
    RequiresInvariants,
    /// The manifest does not establish whether filtering removed sites and
    /// the estimator cannot tolerate the ambiguity.
    UnknownSiteCompleteness,
}

impl RefusalCode {
    /// Human-readable explanation.
    #[must_use]
    pub const fn message(self) -> &'static str {
        match self {
            Self::RequiresInvariants => {
                "estimator requires invariant sites, which the input does not carry"
            }
            Self::UnknownSiteCompleteness => {
                "cannot establish whether filtering removed sites"
            }
        }
    }

    /// Stable token for JSON sidecars.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::RequiresInvariants => "requires_invariants",
            Self::UnknownSiteCompleteness => "unknown_site_completeness",
        }
    }
}

impl From<super::SiteCompleteness> for Ascertainment {
    fn from(v1: super::SiteCompleteness) -> Self {
        match v1 {
            super::SiteCompleteness::WithInvariants => Self {
                sites: SiteSet::Complete,
                invariants: InvariantSet::Present,
                genotypes: GenotypeState::Original,
                cohort: CohortState::Full,
            },
            super::SiteCompleteness::VariantOnly(v1_asc) => Self {
                sites: match v1_asc {
                    super::Ascertainment::Complete => SiteSet::Complete,
                    super::Ascertainment::Filtered | super::Ascertainment::Unknown => {
                        SiteSet::Unknown
                    }
                },
                invariants: InvariantSet::Absent,
                genotypes: GenotypeState::Original,
                cohort: CohortState::Full,
            },
        }
    }
}
}
