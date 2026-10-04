//! Absolute between-group divergence (dxy).

use crate::manifest::{Ascertainment, Availability, SiteCompleteness};
use crate::population::PopulationView;
use crate::result::{DxyResult, PairCounts};
use crate::window::WindowBuffer;

use super::Unavailable;

/// Sealed dxy estimator definitions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DxyEstimator {
    /// Between-group pairwise differences over comparable cross-pairs at the
    /// sites present (pixy semantics).
    PerPairComparable,
    /// Mean cross-pairwise difference per site, summed over sites and divided
    /// by the window length in base pairs (scikit-allel `sequence_divergence`
    /// semantics).
    AccessibleBp,
}

impl DxyEstimator {
    pub fn id(self) -> &'static str {
        match self {
            DxyEstimator::PerPairComparable => "dxy.per_pair_comparable",
            DxyEstimator::AccessibleBp => "dxy.accessible_bp",
        }
    }

    pub fn citation(self) -> &'static str {
        match self {
            DxyEstimator::PerPairComparable => {
                "pixy, Korunes & Samuk 2023 eLife (between-group formulation)"
            }
            DxyEstimator::AccessibleBp => {
                "scikit-allel sequence_divergence (mean cross-pairwise difference per \
                 accessible bp)"
            }
        }
    }

    pub fn availability(self, manifest: SiteCompleteness) -> Availability {
        use Ascertainment::{Complete, Filtered, Unknown};
        use Availability::{Available, Biased, Unavailable};
        match (self, manifest) {
            (_, SiteCompleteness::WithInvariants) => Available,
            (DxyEstimator::AccessibleBp, SiteCompleteness::VariantOnly(Complete)) => Available,
            (DxyEstimator::AccessibleBp, SiteCompleteness::VariantOnly(Filtered)) => Biased(
                "filtered variant-only input: removed true variants are indistinguishable \
                 from invariant sites; the numerator is undercounted",
            ),
            (DxyEstimator::AccessibleBp, SiteCompleteness::VariantOnly(Unknown)) => Unavailable(
                "variant-only input with unknown ascertainment: declare Complete/Filtered or \
                 provide an all-sites VCF",
            ),
            (DxyEstimator::PerPairComparable, SiteCompleteness::VariantOnly(Complete)) => Biased(
                "variant-only input: the denominator counts only variant sites; invariant \
                 sites are absent and dxy is overestimated",
            ),
            (DxyEstimator::PerPairComparable, SiteCompleteness::VariantOnly(Filtered)) => {
                Unavailable(
                    "filtered variant-only input: the per-pair denominator cannot count \
                     invariant sites and removed variants are not recoverable",
                )
            }
            (DxyEstimator::PerPairComparable, SiteCompleteness::VariantOnly(Unknown)) => {
                Unavailable(
                    "variant-only input with unknown ascertainment: per-pair denominators \
                     require monomorphic sites, or a declared Complete file",
                )
            }
        }
    }
}

/// Computes dxy between two groups.
pub fn dxy(
    window: &WindowBuffer,
    pops: &PopulationView,
    group_a: u32,
    group_b: u32,
    estimator: DxyEstimator,
) -> Result<DxyResult, Unavailable> {
    if group_a == group_b {
        return Err(Unavailable {
            reason: "dxy requires two distinct groups",
        });
    }
    let bias = match estimator.availability(window.manifest) {
        Availability::Available => None,
        Availability::Biased(reason) => Some(reason),
        Availability::Unavailable(reason) => return Err(Unavailable { reason }),
    };

    let samples_a = pops.samples(group_a);
    let samples_b = pops.samples(group_b);
    let n_chrom_a = 2 * samples_a.len() as u64;
    let n_chrom_b = 2 * samples_b.len() as u64;
    let total_pairs = n_chrom_a * n_chrom_b;

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
        let a = window.allele_counts(site, samples_a);
        let b = window.allele_counts(site, samples_b);
        let cross_diffs = a.ref_count * b.alt_count + a.alt_count * b.ref_count;
        let comparable = (a.ref_count + a.alt_count) * (b.ref_count + b.alt_count);
        counts.diffs += cross_diffs;
        counts.comparable += comparable;
        counts.missing += total_pairs - comparable;
        if comparable > 0 {
            sum_site_means += cross_diffs as f64 / comparable as f64;
        }
    }

    Ok(DxyResult {
        estimator,
        contig: window.contig,
        start: window.start,
        end: window.end,
        group_a,
        group_b,
        counts,
        sum_site_means,
        n_variant_sites,
        n_invariant_sites,
        bias,
    })
}
