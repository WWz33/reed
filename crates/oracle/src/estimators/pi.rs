//! Nucleotide diversity (π), per group.

use crate::manifest::{Ascertainment, Availability, SiteCompleteness};
use crate::population::PopulationView;
use crate::result::{PairCounts, PiResult};
use crate::window::WindowBuffer;

use super::{samples_of, Unavailable};

/// Sealed π estimator definitions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PiEstimator {
    /// Pairwise differences over comparable pairs at the sites present
    /// (Korunes & Samuk 2021 semantics). Requires all-sites input for an
    /// exact denominator.
    PerPairComparable,
    /// Mean pairwise difference per site, summed over sites and divided by the
    /// window length in base pairs (scikit-allel `sequence_diversity`
    /// semantics). Sites with no comparable pair contribute zero to the
    /// numerator.
    AccessibleBp,
}

impl PiEstimator {
    /// Stable identifier for reports and golden tests.
    pub fn id(self) -> &'static str {
        match self {
            PiEstimator::PerPairComparable => "pi.per_pair_comparable",
            PiEstimator::AccessibleBp => "pi.accessible_bp",
        }
    }

    /// Literature the definition follows.
    pub fn citation(self) -> &'static str {
        match self {
            PiEstimator::PerPairComparable => {
                "Korunes & Samuk 2023 eLife; per-site counts as in scikit-allel \
                 mean_pairwise_difference"
            }
            PiEstimator::AccessibleBp => {
                "scikit-allel sequence_diversity (mean pairwise difference per accessible bp)"
            }
        }
    }

    /// Whether this estimator is computable under a manifest.
    pub fn availability(self, manifest: SiteCompleteness) -> Availability {
        use Ascertainment::{Complete, Filtered, Unknown};
        use Availability::{Available, Biased, Unavailable};
        match (self, manifest) {
            (_, SiteCompleteness::WithInvariants) => Available,
            (PiEstimator::AccessibleBp, SiteCompleteness::VariantOnly(Complete)) => Available,
            (PiEstimator::AccessibleBp, SiteCompleteness::VariantOnly(Filtered)) => Biased(
                "filtered variant-only input: removed true variants are indistinguishable \
                 from invariant sites; the numerator is undercounted",
            ),
            (PiEstimator::AccessibleBp, SiteCompleteness::VariantOnly(Unknown)) => Unavailable(
                "variant-only input with unknown ascertainment: declare Complete/Filtered or \
                 provide an all-sites VCF",
            ),
            (PiEstimator::PerPairComparable, SiteCompleteness::VariantOnly(Complete)) => Biased(
                "variant-only input: the denominator counts only variant sites; invariant \
                 sites are absent and pi is overestimated",
            ),
            (PiEstimator::PerPairComparable, SiteCompleteness::VariantOnly(Filtered)) => {
                Unavailable(
                    "filtered variant-only input: the per-pair denominator cannot count \
                     invariant sites and removed variants are not recoverable",
                )
            }
            (PiEstimator::PerPairComparable, SiteCompleteness::VariantOnly(Unknown)) => {
                Unavailable(
                    "variant-only input with unknown ascertainment: per-pair denominators \
                     require monomorphic sites, or a declared Complete file",
                )
            }
        }
    }
}

/// Computes π for one group.
pub fn pi(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
    estimator: PiEstimator,
) -> Result<PiResult, Unavailable> {
    let samples = samples_of(pops, window, group)?;

    let bias = match estimator.availability(window.manifest) {
        Availability::Available => None,
        Availability::Biased(reason) => Some(reason),
        Availability::Unavailable(reason) => return Err(Unavailable { reason }),
    };

    let n_chrom = 2 * samples.len() as u64;
    let total_pairs = n_chrom * n_chrom.saturating_sub(1) / 2;

    let mut counts = PairCounts::default();
    let mut sum_site_means = 0.0f64;
    let mut n_variant_sites = 0u64;
    let mut n_invariant_sites = 0u64;

    for site in 0..window.n_sites() {
        if window.is_variant.get(site) {
            n_variant_sites += 1;
        } else {
            n_invariant_sites += 1;
        }
        let ac = window.allele_counts(site, samples);
        let n = ac.ref_count + ac.alt_count;
        let comps = n * n.saturating_sub(1) / 2;
        let diffs = ac.ref_count * ac.alt_count;
        counts.diffs += diffs;
        counts.comparable += comps;
        counts.missing += total_pairs - comps;
        if comps > 0 {
            sum_site_means += diffs as f64 / comps as f64;
        }
    }

    Ok(PiResult {
        estimator,
        contig: window.contig,
        start: window.start,
        end: window.end,
        group,
        counts,
        sum_site_means,
        n_variant_sites,
        n_invariant_sites,
        bias,
    })
}
