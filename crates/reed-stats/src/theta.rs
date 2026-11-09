//! Watterson's θ in the reference missing-data form (Watterson 1975;
//! mutation-count per observed-allele class, Korunes & Samuk 2021,
//! `calc_watterson_theta`).
//!
//! Sites are grouped by their observed allele count `n_i`; each variant
//! site contributes `k_i − 1` mutations (k = distinct observed alleles, so
//! 1 for a biallelic segregating site) and each class contributes
//! `mutations(n) / a(n)` with `a(n) = Σ_{j=1}^{n−1} 1/j`.

use reed_core::count::SiteGroupCounts;
use reed_core::math::harmonic::a1;

/// θ accumulator for one group.
#[derive(Clone, Debug)]
pub struct ThetaAcc {
    group: usize,
    n_chrom: u64,
    /// `mutations[n]` = Σ (k_i − 1) over variant sites with `n` observed alleles.
    mutations: Vec<u64>,
    /// `sites[n]` = number of sites with `n` observed alleles.
    sites: Vec<u64>,
    /// Σ n_i over sites with n_i > 0 (for the window-average n).
    total_allele_count: u64,
    /// Sites with at least one observed allele.
    n_sites: u64,
    /// Variant sites (any ALT observed).
    n_var_sites: u64,
    /// Terminal scalar form used for CLI window rows. The distribution
    /// vectors are released once no later merge is possible (doc_v27 V27-18).
    compact: Option<ThetaCompact>,
}

#[derive(Clone, Debug)]
struct ThetaCompact {
    raw_theta: f64,
    weighted_no_sites: f64,
    n_mutations: u64,
}

impl ThetaAcc {
    /// Accumulator for `group`; `None` when the group has < 2 chromosomes.
    #[must_use]
    pub fn new(group: usize, n_group: u32) -> Option<Self> {
        let n_chrom = u64::from(n_group) * 2;
        if n_chrom < 2 {
            return None;
        }
        let width = usize::try_from(n_chrom)
            .unwrap_or(usize::MAX)
            .saturating_add(1);
        Some(Self {
            group,
            n_chrom,
            mutations: vec![0; width],
            sites: vec![0; width],
            total_allele_count: 0,
            n_sites: 0,
            n_var_sites: 0,
            compact: None,
        })
    }

    /// Absorbs one site (all observed-allele counts; missing data included).
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let c = &row[self.group];
        let n = u64::from(c.n_ref_alleles) + u64::from(c.n_alt_alleles);
        if n == 0 {
            return;
        }
        let n_idx = usize::try_from(n).unwrap_or(usize::MAX);
        if n_idx >= self.sites.len() {
            self.sites.resize(n_idx + 1, 0);
            self.mutations.resize(n_idx + 1, 0);
        }
        self.sites[n_idx] += 1;
        self.total_allele_count += n;
        self.n_sites += 1;
        let alt = u64::from(c.n_alt_alleles);
        if alt > 0 {
            self.n_var_sites += 1;
        }
        // the reference counts mutations as k − 1 (k = distinct observed alleles), so a
        // site fixed for ALT has k = 1 and contributes 0 — it is not
        // segregating. Biallelic sites with both alleles observed contribute 1.
        // Multi-allelic sites need per-allele counts (04-algorithms.md §3),
        // which the count layer does not carry yet; the driver excludes such
        // records for now.
        if alt > 0 && alt < n {
            self.mutations[n_idx] += 1;
        }
    }

    /// Component-wise merge. Histogram lengths are storage shapes, not the
    /// request identity: omitted high-allele classes are zero. A caller must
    /// never merge a different group/chromosome request or a compacted row.
    pub fn merge(&mut self, other: &Self) {
        assert!(
            self.group == other.group && self.n_chrom == other.n_chrom,
            "ThetaAcc: merge operands have different group/chromosome semantics ({}:{} vs {}:{})",
            self.group,
            self.n_chrom,
            other.group,
            other.n_chrom
        );
        assert!(
            self.compact.is_none() && other.compact.is_none(),
            "a compacted ThetaAcc is a terminal output, not a merge input"
        );
        assert_eq!(
            self.sites.len(),
            self.mutations.len(),
            "ThetaAcc: sites/mutations storage shapes differ"
        );
        assert_eq!(
            other.sites.len(),
            other.mutations.len(),
            "ThetaAcc: sites/mutations storage shapes differ"
        );
        let len = self.sites.len().max(other.sites.len());
        self.sites.resize(len, 0);
        self.mutations.resize(len, 0);
        for (dst, src) in self.sites.iter_mut().zip(&other.sites) {
            *dst += src;
        }
        for (dst, src) in self.mutations.iter_mut().zip(&other.mutations) {
            *dst += src;
        }
        self.total_allele_count += other.total_allele_count;
        self.n_sites += other.n_sites;
        self.n_var_sites += other.n_var_sites;
    }

    /// `Σ_n mutations(n) / a(n)` — the window's raw Watterson θ.
    #[must_use]
    pub fn raw_theta(&self) -> f64 {
        if let Some(compact) = &self.compact {
            return compact.raw_theta;
        }
        let mut raw = 0.0;
        for (n, &mutations) in self.mutations.iter().enumerate() {
            if mutations == 0 || n < 2 {
                continue;
            }
            let a = a1(n);
            if a > 0.0 {
                raw += mutations as f64 / a;
            }
        }
        raw
    }

    /// Window-**average** θ = `raw / n_sites` (the CLI's
    /// `avg_watterson_theta` column). `raw_theta()` is the window total
    /// (the `raw_watterson_theta` column; the oracle's `total_theta`); the
    /// old name `value_total` was misleading — this is a per-site mean
    /// (doc_v7 D7-17).
    #[must_use]
    pub fn value_avg(&self) -> Option<f64> {
        (self.n_sites > 0).then(|| self.raw_theta() / self.n_sites as f64)
    }

    /// `--variant-only` form: `raw / L` (04-algorithms.md §3).
    #[must_use]
    pub fn value_variant_only(&self, window_bp: u64) -> Option<f64> {
        (window_bp > 0).then(|| self.raw_theta() / window_bp as f64)
    }

    /// Watterson θ per physical bp: `θ_raw / L`. This is the wide-TSV `bp`
    /// quantity, not the window-average θ divided by L a second time (V27-38).
    #[must_use]
    pub fn value_per_bp(&self, window_bp: u64) -> Option<f64> {
        (self.n_sites > 0 && window_bp > 0).then(|| self.raw_theta() / window_bp as f64)
    }

    /// The reference `weighted_no_sites`: `Σ_n sites(n) · n / max(n)`.
    #[must_use]
    pub fn weighted_no_sites(&self) -> f64 {
        if let Some(compact) = &self.compact {
            return compact.weighted_no_sites;
        }
        let Some(max_n) = self.sites.iter().rposition(|&count| count > 0) else {
            return 0.0;
        };
        if max_n == 0 {
            return 0.0;
        }
        let mut weighted = 0.0;
        for (n, &count) in self.sites.iter().enumerate() {
            if count > 0 {
                weighted += count as f64 * n as f64 / max_n as f64;
            }
        }
        weighted
    }

    /// Sites with ≥ 1 observed allele.
    #[must_use]
    pub const fn n_sites(&self) -> u64 {
        self.n_sites
    }

    /// Variant sites.
    #[must_use]
    pub const fn n_var_sites(&self) -> u64 {
        self.n_var_sites
    }

    /// Σ observed alleles over sites with calls (Tajima's mean-n input).
    #[must_use]
    pub const fn total_allele_count(&self) -> u64 {
        self.total_allele_count
    }

    /// `Σ (k_i − 1)` over variant sites — Tajima's `S`.
    #[must_use]
    pub fn n_mutations(&self) -> u64 {
        self.compact
            .as_ref()
            .map_or_else(|| self.mutations.iter().sum(), |c| c.n_mutations)
    }

    /// Converts a finished window's accumulator to its scalar output form.
    /// The caller must not merge it afterwards; contig/genome aggregation
    /// still uses the original full accumulators. This is what lets the CLI
    /// discard O(sample-count) θ vectors per emitted window (V27-18).
    pub fn compact(&mut self) {
        if self.compact.is_some() {
            return;
        }
        self.compact = Some(ThetaCompact {
            raw_theta: self.raw_theta(),
            weighted_no_sites: self.weighted_no_sites(),
            n_mutations: self.n_mutations(),
        });
        // `clear()` keeps the allocation; replacing the vectors actually
        // releases the per-window distribution storage (V27-18).
        self.mutations = Vec::new();
        self.sites = Vec::new();
    }

    /// Contributing / skipped site tally (skipped is always 0 under the reference
    /// semantics: missing data is incorporated, not dropped).
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, 0)
    }

    /// Segregating-site count.
    #[must_use]
    pub const fn segregating(&self) -> u64 {
        self.n_var_sites
    }

    /// Chromosome count.
    #[must_use]
    pub const fn n_chrom(&self) -> u64 {
        self.n_chrom
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(n_ref_alleles: u32, n_alt_alleles: u32) -> SiteGroupCounts {
        SiteGroupCounts {
            n_ref_alleles,
            n_alt_alleles,
            ..Default::default()
        }
    }

    fn trim_zero_tail(acc: &mut ThetaAcc) {
        assert_eq!(acc.sites.len(), acc.mutations.len());
        let len = (0..acc.sites.len())
            .rev()
            .find(|&i| acc.sites[i] != 0 || acc.mutations[i] != 0)
            .map_or(0, |i| i + 1);
        acc.sites.truncate(len);
        acc.mutations.truncate(len);
    }

    fn assert_same_distribution(actual: &ThetaAcc, expected: &ThetaAcc) {
        let len = actual.sites.len().max(expected.sites.len());
        for i in 0..len {
            assert_eq!(
                actual.sites.get(i).copied().unwrap_or(0),
                expected.sites.get(i).copied().unwrap_or(0),
                "sites class {i}"
            );
            assert_eq!(
                actual.mutations.get(i).copied().unwrap_or(0),
                expected.mutations.get(i).copied().unwrap_or(0),
                "mutations class {i}"
            );
        }
        assert_eq!(actual.raw_theta(), expected.raw_theta());
        assert_eq!(actual.weighted_no_sites(), expected.weighted_no_sites());
        assert_eq!(actual.n_sites(), expected.n_sites());
        assert_eq!(actual.n_var_sites(), expected.n_var_sites());
        assert_eq!(actual.total_allele_count(), expected.total_allele_count());
        assert_eq!(actual.n_mutations(), expected.n_mutations());
        assert_eq!(actual.value_avg(), expected.value_avg());
    }

    #[test]
    fn merge_accepts_equivalent_zero_trimmed_histogram_tails() {
        let low = site(1, 1);
        let high = site(2, 4);
        let mut whole = ThetaAcc::new(0, 4).expect("invariant: nonempty group");
        whole.push(std::slice::from_ref(&low));
        whole.push(std::slice::from_ref(&high));

        let mut a = ThetaAcc::new(0, 4).expect("invariant: nonempty group");
        a.push(std::slice::from_ref(&low));
        let mut b = ThetaAcc::new(0, 4).expect("invariant: nonempty group");
        b.push(std::slice::from_ref(&high));
        trim_zero_tail(&mut a);
        trim_zero_tail(&mut b);
        assert_ne!(a.sites.len(), b.sites.len(), "failing-first fixture");

        let mut forward = a.clone();
        forward.merge(&b);
        assert_same_distribution(&forward, &whole);
        let mut reverse = b.clone();
        reverse.merge(&a);
        assert_same_distribution(&reverse, &whole);
    }

    #[test]
    #[should_panic(expected = "different group/chromosome semantics")]
    fn merge_rejects_a_different_group() {
        let mut a = ThetaAcc::new(0, 4).expect("invariant: nonempty group");
        let b = ThetaAcc::new(1, 4).expect("invariant: nonempty group");
        a.merge(&b);
    }

    #[test]
    #[should_panic(expected = "different group/chromosome semantics")]
    fn merge_rejects_a_different_chromosome_count() {
        let mut a = ThetaAcc::new(0, 4).expect("invariant: nonempty group");
        let b = ThetaAcc::new(0, 5).expect("invariant: nonempty group");
        a.merge(&b);
    }

    #[test]
    #[should_panic(expected = "terminal output, not a merge input")]
    fn merge_rejects_a_compacted_accumulator() {
        let mut a = ThetaAcc::new(0, 4).expect("invariant: nonempty group");
        a.compact();
        let b = ThetaAcc::new(0, 4).expect("invariant: nonempty group");
        a.merge(&b);
    }

    #[test]
    #[should_panic(expected = "sites/mutations storage shapes differ")]
    fn merge_rejects_an_internally_invalid_histogram() {
        let mut a = ThetaAcc::new(0, 4).expect("invariant: nonempty group");
        a.sites.pop();
        let b = ThetaAcc::new(0, 4).expect("invariant: nonempty group");
        a.merge(&b);
    }
}
