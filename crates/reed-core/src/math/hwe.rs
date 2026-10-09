//! Hardy-Weinberg exact test (Wigginton et al. 2005) and the mid-p variant
//! (Graffelman & Moreno 2013), with a bounded memo shared by the filter and
//! QC consumers (`docs/spec/architecture.md` §9.6).
//!
//! Given `n` called diploid genotypes and `mac` minor-allele copies, the
//! conditional probability of observing `nh` heterozygotes under random
//! mating is proportional to `n!/(hm!·nh!·hM!) · 2^nh` with
//! `hm = (mac − nh)/2` and `hM = n − nh − hm`, summed over all valid `nh`.
//! The two-sided p-value adds every probability at or below the observed
//! one; mid-p subtracts half the observed probability.
//!
//! Degenerate inputs (monomorphic, single valid `nh`) return `p = 1`.
//! Impossible combinations (`mac > 2n`, `nh` of wrong parity or outside the
//! support) return `NaN` — callers never produce them from real counts.

use std::collections::HashMap;

/// Memo capacity (`docs/spec/architecture.md` §13.8).
pub const HWE_MEMO_CAPACITY: usize = 8192;

/// Normalized conditional probabilities of every valid het count, ascending
/// from the smallest (`mac % 2`), with that smallest count.
#[derive(Clone, Debug)]
pub struct HweDistribution {
    /// Het count of the first entry.
    pub nh_min: u32,
    /// Normalized probabilities, one per valid `nh`, ascending in steps of 2.
    pub probs: Vec<f64>,
}

impl HweDistribution {
    /// Probability of observing `nh` heterozygotes (0 for invalid `nh`).
    #[must_use]
    pub fn prob(&self, nh: u32) -> f64 {
        if nh < self.nh_min || !(nh - self.nh_min).is_multiple_of(2) {
            return 0.0;
        }
        let idx = ((nh - self.nh_min) / 2) as usize;
        if idx >= self.probs.len() {
            0.0
        } else {
            self.probs[idx]
        }
    }

    /// Two-sided exact p-value for the observed het count.
    #[must_use]
    pub fn p_exact(&self, nh: u32) -> f64 {
        let obs = self.prob(nh);
        self.probs.iter().filter(|&&p| p <= obs).sum()
    }

    /// Mid-p two-sided value: exact p minus half the observed probability.
    #[must_use]
    pub fn p_mid(&self, nh: u32) -> f64 {
        self.p_exact(nh) - 0.5 * self.prob(nh)
    }
}

/// Conditional genotype-count distribution; impossible inputs give `None`.
//
// The recurrence runs in u64 counts; every value is bounded by
// `2 * n_called` (≤ 2^33 for u32 inputs), far below f64's exact-integer
// range (2^53) and `usize` on every supported target.
#[allow(clippy::cast_possible_truncation, clippy::cast_precision_loss)]
#[must_use]
pub fn hwe_distribution(n_called: u32, n_minor: u32) -> Option<HweDistribution> {
    let (n, mac) = (u64::from(n_called), u64::from(n_minor));
    if n == 0 || mac == 0 || mac == 2 * n {
        // Degenerate: no het variation to test against.
        return Some(HweDistribution {
            nh_min: 0,
            probs: vec![1.0],
        });
    }
    if mac > 2 * n {
        return None;
    }
    let nh_min = mac % 2;
    // nh ≤ mac (each het carries one minor allele) and
    // nh ≤ 2n − mac (enough major alleles to fill the rest).
    let nh_max = mac.min(2 * n - mac);
    if nh_min > nh_max {
        return None;
    }

    // Relative log-probabilities via the recurrence
    // P(nh+2)/P(nh) = 4·hm·hM / ((nh+1)(nh+2)), where hm = (mac − nh)/2 and
    // hM = n − nh − hm (both taken *before* the increment: going nh → nh+2
    // turns one hom-minor into hets and one hom-major into a het).
    let mut log_weights = Vec::with_capacity(((nh_max - nh_min) / 2 + 1) as usize);
    let mut log_w = 0.0_f64;
    let mut nh = nh_min;
    while nh <= nh_max {
        log_weights.push(log_w);
        let hm = (mac - nh) / 2;
        let hmaj = n - nh - hm;
        if hm == 0 || hmaj == 0 {
            break;
        }
        log_w += 4.0_f64.ln() + (hm as f64).ln() + (hmaj as f64).ln()
            - (nh as f64 + 1.0).ln()
            - (nh as f64 + 2.0).ln();
        nh += 2;
    }
    let max_log = log_weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = log_weights.iter().map(|&lw| (lw - max_log).exp()).sum();
    let probs = log_weights.iter().map(|&lw| (lw - max_log).exp() / total).collect();
    Some(HweDistribution {
        nh_min: u32::try_from(nh_min).ok()?,
        probs,
    })
}

/// Exact-test p-value without memoization.
#[must_use]
pub fn hwe_exact_p(n_called: u32, n_minor: u32, n_het: u32) -> f64 {
    hwe_distribution(n_called, n_minor)
        .map_or(f64::NAN, |d| d.p_exact(n_het))
}

/// Mid-p p-value without memoization.
#[must_use]
pub fn hwe_midp_p(n_called: u32, n_minor: u32, n_het: u32) -> f64 {
    hwe_distribution(n_called, n_minor)
        .map_or(f64::NAN, |d| d.p_mid(n_het))
}

/// Bounded memo for `(n_called, n_minor) → distribution`, shared within one
/// worker. Capacity is [`HWE_MEMO_CAPACITY`]; on overflow the whole cache is
/// dropped (simple and bounded; a true LRU replaces it if hit-rate profiling
/// ever demands).
///
/// # ponytail: clear-on-full memo, swap for a real LRU if hit rate matters.
#[derive(Debug, Default)]
pub struct HweCache {
    map: HashMap<(u32, u32), HweDistribution>,
}

impl HweCache {
    /// Creates an empty cache.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Distribution for `(n_called, n_minor)`, cached.
    #[must_use]
    pub fn distribution(&mut self, n_called: u32, n_minor: u32) -> Option<HweDistribution> {
        if let Some(d) = self.map.get(&(n_called, n_minor)) {
            return Some(d.clone());
        }
        let d = hwe_distribution(n_called, n_minor)?;
        if self.map.len() >= HWE_MEMO_CAPACITY {
            self.map.clear();
        }
        self.map.insert((n_called, n_minor), d.clone());
        Some(d)
    }

    /// Exact p-value for `(n_called, n_minor, n_het)`.
    #[must_use]
    pub fn p_exact(&mut self, n_called: u32, n_minor: u32, n_het: u32) -> f64 {
        self.distribution(n_called, n_minor).map_or(f64::NAN, |d| d.p_exact(n_het))
    }

    /// Mid-p value for `(n_called, n_minor, n_het)`.
    #[must_use]
    pub fn p_mid(&mut self, n_called: u32, n_minor: u32, n_het: u32) -> f64 {
        self.distribution(n_called, n_minor).map_or(f64::NAN, |d| d.p_mid(n_het))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn degenerate_sites_get_p_one() {
        assert_eq!(hwe_exact_p(10, 0, 0), 1.0);
        assert_eq!(hwe_exact_p(10, 20, 0), 1.0);
        // n=1: the only valid het count is 1 → nothing more extreme exists.
        assert_eq!(hwe_exact_p(1, 1, 1), 1.0);
    }

    #[test]
    fn balanced_hets_are_not_significant() {
        // n=100, mac=100, nh=50 (hm=25, hM=25): the maximum-probability
        // configuration, p must be close to 1.
        let p = hwe_exact_p(100, 100, 50);
        assert!(p > 0.99, "p={p}");
    }

    #[test]
    fn het_excess_and_deficit_are_significant() {
        // All 100 minor alleles sit in heterozygotes of n=100 individuals
        // (expected het count is 50): the most extreme configuration.
        let p = hwe_exact_p(100, 100, 100);
        assert!(p < 1e-6, "p={p}");
        // Complete het deficit: zero heterozygotes carrying 100 minor alleles.
        let p0 = hwe_exact_p(1000, 100, 0);
        assert!(p0 < 1e-6, "p={p0}");
    }

    #[test]
    fn mid_p_below_exact() {
        let d = hwe_distribution(100, 100).unwrap();
        assert!(d.p_mid(50) <= d.p_exact(50));
    }

    #[test]
    fn probabilities_normalized() {
        let d = hwe_distribution(50, 30).unwrap();
        let total: f64 = d.probs.iter().sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cache_matches_uncached() {
        let mut cache = HweCache::new();
        for (n, mac, nh) in [(100_u32, 100_u32, 50_u32), (1000, 100, 100), (50, 30, 15)] {
            assert_eq!(cache.p_exact(n, mac, nh), hwe_exact_p(n, mac, nh));
        }
    }
}
