//! Nucleotide diversity π (`docs/spec/architecture.md` §9.1), the exact
//! accumulators behind the oracle's `estimators::pi`.

use reed_core::count::SiteGroupCounts;
use reed_core::math::sum::NeumaierSum;

/// Additive pair counts (the oracle's `PairCounts` audit trail): all three
/// components are exact integers.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PairAcc {
    /// Called REF×ALT product summed over sites.
    pub diffs: u64,
    /// Comparable pairs summed over sites.
    pub comparable: u64,
    /// Pairs lost to missing alleles.
    pub missing: u64,
}

impl PairAcc {
    /// Absorbs one site's contribution.
    pub fn add(&mut self, diffs: u64, comparable: u64, total_pairs: u64) {
        self.diffs += diffs;
        self.comparable += comparable;
        self.missing += total_pairs - comparable;
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.diffs += other.diffs;
        self.comparable += other.comparable;
        self.missing += other.missing;
    }
}

/// π accumulator for one group. `value_per_pair` is the reference-style
/// PerPairComparable estimate; `value_accessible_bp` the scikit-allel
/// sequence-diversity convention (per window bp).
#[derive(Clone, Debug)]
pub struct PiAcc {
    group: usize,
    n_chrom: u64,
    total_pairs: u64,
    counts: PairAcc,
    sum_site_means: NeumaierSum,
    n_variant_sites: u64,
    n_invariant_sites: u64,
}

impl PiAcc {
    /// Accumulator for `group` (index into the count row).
    #[must_use]
    pub fn new(group: usize, n_group: u32) -> Self {
        let n_chrom = u64::from(n_group) * 2;
        Self {
            group,
            n_chrom,
            total_pairs: n_chrom * n_chrom.saturating_sub(1) / 2,
            counts: PairAcc::default(),
            sum_site_means: NeumaierSum::new(),
            n_variant_sites: 0,
            n_invariant_sites: 0,
        }
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts], all_variant: bool) {
        let c = &row[self.group];
        let n = u64::from(c.n_ref_alleles + c.n_alt_alleles);
        let comps = n * n.saturating_sub(1) / 2;
        let diffs = u64::from(c.n_ref_alleles) * u64::from(c.n_alt_alleles);
        self.counts.add(diffs, comps, self.total_pairs);
        if comps > 0 {
            self.sum_site_means.add(diffs as f64 / comps as f64);
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

    /// Called-allele count bookkeeping (for tests): the group's chromosome
    /// count.
    #[must_use]
    pub const fn n_chrom(&self) -> u64 {
        self.n_chrom
    }

    /// Integer audit trail.
    #[must_use]
    pub const fn counts(&self) -> &PairAcc {
        &self.counts
    }

    /// PerPairComparable: `Σ diffs / Σ comparable` over the window.
    #[must_use]
    pub fn value_per_pair(&self) -> Option<f64> {
        (self.counts.comparable > 0)
            .then(|| self.counts.diffs as f64 / self.counts.comparable as f64)
    }

    /// AccessibleBp: `Σ per-site mean pairwise difference / window bp`.
    #[must_use]
    pub fn value_accessible_bp(&self, window_bp: u64) -> Option<f64> {
        (self.counts.comparable > 0 && window_bp > 0)
            .then(|| self.sum_site_means.value() / window_bp as f64)
    }

    /// Site tally (variant / invariant under the site-wide flag).
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_variant_sites, self.n_invariant_sites)
    }

    /// Sites the window's file actually contained (`|S|` in 04-algorithms §3).
    #[must_use]
    pub const fn n_present(&self) -> u64 {
        self.n_variant_sites + self.n_invariant_sites
    }

    /// π with the `--variant-only` denominator completion: every position of
    /// the window the file does not contain counts as monomorphic with all
    /// `2N` chromosomes called, i.e. `(L − |S|) · C(2N, 2)` extra comparable
    /// pairs (04-algorithms.md §3; same convention as vcftools --window-pi).
    #[must_use]
    pub fn value_variant_only(&self, window_bp: u64) -> Option<f64> {
        let absent = window_bp.saturating_sub(self.n_present());
        let comparable = self.counts.comparable + absent * self.total_pairs;
        (comparable > 0).then(|| self.counts.diffs as f64 / comparable as f64)
    }
}
