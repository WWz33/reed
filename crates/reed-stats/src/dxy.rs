//! Between-group absolute divergence dxy (`docs/spec/architecture.md` §9.1).

use reed_core::count::SiteGroupCounts;
use reed_core::math::sum::NeumaierSum;

use crate::pi::PairAcc;

/// dxy accumulator for one group pair (same two estimator conventions as
/// [`crate::pi::PiAcc`]).
#[derive(Clone, Debug)]
pub struct DxyAcc {
    group_a: usize,
    group_b: usize,
    total_pairs: u64,
    counts: PairAcc,
    sum_site_means: NeumaierSum,
    n_variant_sites: u64,
    n_invariant_sites: u64,
}

impl DxyAcc {
    /// Accumulator for the `group_a`/`group_b` pair.
    #[must_use]
    pub fn new(group_a: usize, group_b: usize, n_group_a: u32, n_group_b: u32) -> Self {
        Self {
            group_a,
            group_b,
            total_pairs: u64::from(n_group_a) * 2 * (u64::from(n_group_b) * 2),
            counts: PairAcc::default(),
            sum_site_means: NeumaierSum::new(),
            n_variant_sites: 0,
            n_invariant_sites: 0,
        }
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts], all_variant: bool) {
        let (a, b) = (&row[self.group_a], &row[self.group_b]);
        let cross_diffs = u64::from(a.n_ref_alleles) * u64::from(b.n_alt_alleles)
            + u64::from(a.n_alt_alleles) * u64::from(b.n_ref_alleles);
        let comparable = u64::from(a.n_ref_alleles + a.n_alt_alleles)
            * u64::from(b.n_ref_alleles + b.n_alt_alleles);
        self.counts.add(cross_diffs, comparable, self.total_pairs);
        if comparable > 0 {
            self.sum_site_means
                .add(cross_diffs as f64 / comparable as f64);
        }
        if all_variant {
            self.n_variant_sites += 1;
        } else {
            self.n_invariant_sites += 1;
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.counts.merge(&other.counts);
        self.sum_site_means.merge(&other.sum_site_means);
        self.n_variant_sites += other.n_variant_sites;
        self.n_invariant_sites += other.n_invariant_sites;
    }

    /// Integer audit trail.
    #[must_use]
    pub const fn counts(&self) -> &PairAcc {
        &self.counts
    }

    /// `PerPairComparable` value.
    #[must_use]
    pub fn value_per_pair(&self) -> Option<f64> {
        (self.counts.comparable > 0)
            .then(|| self.counts.diffs as f64 / self.counts.comparable as f64)
    }

    /// `AccessibleBp` value.
    #[must_use]
    pub fn value_accessible_bp(&self, window_bp: u64) -> Option<f64> {
        (self.counts.comparable > 0 && window_bp > 0)
            .then(|| self.sum_site_means.value() / window_bp as f64)
    }

    /// Sites the window's file actually contained (`|S|`).
    #[must_use]
    pub const fn n_present(&self) -> u64 {
        self.n_variant_sites + self.n_invariant_sites
    }

    /// dxy with the `--variant-only` denominator completion: absent window
    /// positions count as monomorphic with all called, adding
    /// `(L − |S|) · 2N_A · 2N_B` comparable cross pairs (04-algorithms.md §3).
    #[must_use]
    pub fn value_variant_only(&self, window_bp: u64) -> Option<f64> {
        let absent = window_bp.saturating_sub(self.n_present());
        let comparable = self.counts.comparable + absent * self.total_pairs;
        (comparable > 0).then(|| self.counts.diffs as f64 / comparable as f64)
    }
}
