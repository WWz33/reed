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
    /// Pairs lost to missing alleles. Only maintained by accumulators that
    /// know their cohort totals ([`PairAcc::add`]); fst/PBS audits use
    /// [`PairAcc::add_partial`] and leave it 0 (doc_v15 D15-38).
    pub missing: u64,
}

impl PairAcc {
    /// Absorbs one site's contribution.
    pub fn add(&mut self, diffs: u64, comparable: u64, total_pairs: u64) {
        self.diffs += diffs;
        self.comparable += comparable;
        self.missing += total_pairs - comparable;
    }

    /// Absorbs one site's contribution when the total pair count is not
    /// known: fst/PBS audits carry no cohort totals, so `missing` is not
    /// pretended to be tracked (doc_v15 D15-38).
    pub fn add_partial(&mut self, diffs: u64, comparable: u64) {
        self.diffs += diffs;
        self.comparable += comparable;
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
    /// Sites the file contains but the statistic cannot count (a skipped
    /// multiallelic record): present for the variant-only denominator,
    /// contributing no comparable pairs (doc_v7 D7-11).
    n_skipped_present: u64,
}

impl PiAcc {
    /// Accumulator for `group` (index into the count row).
    #[must_use]
    pub fn new(group: usize, n_group: u32) -> Self {
        let n_chrom = u64::from(n_group) * 2;
        Self {
            group,
            n_chrom,
            total_pairs: n_chrom.saturating_mul(n_chrom.saturating_sub(1)) / 2,
            counts: PairAcc::default(),
            sum_site_means: NeumaierSum::new(),
            n_variant_sites: 0,
            n_invariant_sites: 0,
            n_skipped_present: 0,
        }
    }

    /// Registers a site the file contains but the statistic cannot count
    /// (skipped multiallelic): it is present for the variant-only
    /// denominator but adds no pairs (doc_v7 D7-11).
    pub fn note_skipped_site(&mut self) {
        self.n_skipped_present += 1;
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts], all_variant: bool) {
        let c = &row[self.group];
        let n = u64::from(c.n_ref_alleles + c.n_alt_alleles);
        let comps = n.saturating_mul(n.saturating_sub(1)) / 2;
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
        self.n_skipped_present += other.n_skipped_present;
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

    /// Sites the file contained, including the ones this statistic skipped
    /// (multiallelic): the variant-only `absent` completion is over the
    /// file's own positions, not only the counted ones (doc_v7 D7-11).
    #[must_use]
    pub const fn n_present_file(&self) -> u64 {
        self.n_present() + self.n_skipped_present
    }

    /// π with the `--variant-only` denominator completion: every position of
    /// the window the file does not contain counts as monomorphic with all
    /// `2N` chromosomes called, i.e. `(L − |S|) · C(2N, 2)` extra comparable
    /// pairs (04-algorithms.md §3; same convention as vcftools --window-pi).
    #[must_use]
    pub fn value_variant_only(&self, window_bp: u64) -> Option<f64> {
        let absent = window_bp.saturating_sub(self.n_present_file());
        // `absent · total_pairs` can exceed u64 well before the sample count
        // is physically absurd (≈10⁵ samples × a whole-chromosome window);
        // accumulate in u128 so the denominator never wraps (doc_v7 D7-21).
        let comparable =
            u128::from(self.counts.comparable) + u128::from(absent) * u128::from(self.total_pairs);
        (comparable > 0).then(|| self.counts.diffs as f64 / comparable as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// doc_v7 D7-21: a large sample count times a whole-chromosome window
    /// overflows the old `absent * total_pairs` u64 multiplication (a debug
    /// build panicked); the u128 accumulation keeps it finite.
    #[test]
    fn variant_only_denominator_does_not_overflow() {
        let acc = PiAcc::new(0, 100_000); // 200 000 chromosomes
        let value = acc
            .value_variant_only(1_000_000_000)
            .expect("invariant: comparable > 0");
        assert_eq!(value, 0.0); // no sites pushed → no diffs
    }

    /// doc_v7 D7-11: a skipped multiallelic site is present in the file, so
    /// the variant-only denominator must not treat its position as absent.
    #[test]
    fn skipped_multiallelic_is_present_for_variant_only() {
        let row = [SiteGroupCounts {
            n_ref_alleles: 2,
            n_alt_alleles: 2,
            ..SiteGroupCounts::default()
        }];
        let mut counted = PiAcc::new(0, 2);
        counted.push(&row, true);
        let baseline = counted.value_variant_only(100).expect("invariant: defined");
        assert_eq!(counted.n_present(), 1);
        assert_eq!(counted.n_present_file(), 1);

        counted.note_skipped_site();
        assert_eq!(counted.n_present(), 1, "the used-site tally is unchanged");
        assert_eq!(counted.n_present_file(), 2, "the file carries both sites");
        let with_skip = counted.value_variant_only(100).expect("invariant: defined");
        // One fewer absent position → a smaller comparable denominator →
        // a larger π (same diffs).
        assert!(
            with_skip > baseline,
            "the skipped site must shrink the absent count: {with_skip} vs {baseline}"
        );
    }
}
