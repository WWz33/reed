//! Mergeable distribution structures (§10.2): dense exact histograms for
//! integer fields (DP, GQ, indel lengths), fixed-boundary histograms for
//! floats, and per-sample vectors. All satisfy the accumulator laws of
//! §8.4, so QC shares the chunk-parallel driver with the stats layer.

// Casts in this file are bounded by construction: Histogram bins are bounded by the sample count.
#![allow(clippy::cast_possible_truncation)]

use serde::Serialize;

/// Dense exact integer histogram over `0..=max` plus one overflow bucket
/// (§10.2: quantiles are exact for integer fields).
#[derive(Clone, Debug, Serialize)]
pub struct IntHistogram {
    max: u32,
    /// `bins[0..=max]` are the exact values; `bins[max+1]` is overflow.
    bins: Vec<u64>,
}

impl IntHistogram {
    /// A histogram of `0..=max` with an overflow bucket.
    #[must_use]
    pub fn new(max: u32) -> Self {
        Self {
            max,
            bins: vec![0; max as usize + 2],
        }
    }

    /// Absorbs one observation (values above `max` land in the overflow
    /// bucket `max + 1`).
    pub fn push(&mut self, value: u32) {
        let idx = if value > self.max {
            self.bins.len() - 1
        } else {
            value as usize
        };
        self.bins[idx] += 1;
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        debug_assert_eq!(self.max, other.max, "invariant: same histogram range");
        for (dst, src) in self.bins.iter_mut().zip(&other.bins) {
            *dst += src;
        }
    }

    /// Total observations.
    #[must_use]
    pub fn total(&self) -> u64 {
        self.bins.iter().sum()
    }

    /// Sum of observations × value (the mean numerator; overflow values
    /// count as `max + 1`).
    #[must_use]
    pub fn weighted_sum(&self) -> u128 {
        self.bins
            .iter()
            .enumerate()
            .map(|(i, &n)| {
                let value = if i + 1 == self.bins.len() {
                    u128::from(self.max) + 1
                } else {
                    u128::try_from(i).expect("invariant: fits")
                };
                u128::from(n) * value
            })
            .sum()
    }

    /// The `q`-quantile (0..=1); `None` when empty. Exact; the overflow
    /// bucket reads as `max + 1`.
    #[must_use]
    pub fn quantile(&self, q: f64) -> Option<u32> {
        let total = self.total();
        if total == 0 {
            return None;
        }
        let target = ((q * total as f64).ceil() as u64).clamp(1, total);
        let mut acc = 0_u64;
        for (i, &n) in self.bins.iter().enumerate() {
            acc += n;
            if acc >= target {
                let value = u32::try_from(i).ok().unwrap_or(self.max + 1);
                return Some(value.min(if i + 1 == self.bins.len() {
                    self.max + 1
                } else {
                    value
                }));
            }
        }
        None
    }

    /// Iterates `(value, count)` pairs including the overflow bucket as
    /// `max + 1`.
    pub fn iter(&self) -> impl Iterator<Item = (u32, u64)> + '_ {
        self.bins.iter().enumerate().map(|(i, &n)| {
            let value = if i + 1 == self.bins.len() {
                self.max + 1
            } else {
                u32::try_from(i).expect("invariant: fits")
            };
            (value, n)
        })
    }
}

/// bcftools QUAL distribution — an exact replica of `dist_t` from
/// bcftools's `dist.c` as used for QUAL (`dist_init(5)`): values with
/// `iqual ≤ 10^5` are exact (0.1 QUAL resolution); above that, bins widen
/// logarithmically (level `lv−4`, step `10^level`). `iqual` is bcftools's
/// `1 + (int)(qual*10)` computed in f32; 0 = missing/NaN/negative.
#[derive(Clone, Debug, Default, Serialize)]
pub struct QualDist {
    /// dist bin index → count.
    bins: std::collections::BTreeMap<usize, u64>,
}

/// The exact/aggregate boundary (`nexact = 10^npow`, npow = 5).
const NEXACT: usize = 100_000;
/// `nlevel = nexact − 10^(npow−1)` — the number of level-1 bins.
const NLEVEL: usize = 90_000;

impl QualDist {
    /// An empty distribution.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// bcftools's `dist_insert` bin for an iqual value.
    fn ibin(iqual: u32) -> usize {
        let value = u64::from(iqual);
        if value <= NEXACT as u64 {
            return value as usize;
        }
        let lv = (value as f64).log10().floor() as u32;
        let level = lv - 4;
        let step = 10_u64.pow(level);
        NEXACT + NLEVEL * (level as usize - 1) + ((value - 10_u64.pow(lv)) / step) as usize
    }

    /// bcftools's `dist_get` range for a bin index: `[beg, end)` in iqual
    /// units (exact bins span one 0.1 QUAL step).
    #[must_use]
    pub fn bin_range(index: usize) -> (u32, u32) {
        if index < NEXACT {
            return (
                u32::try_from(index).unwrap_or(u32::MAX),
                u32::try_from(index + 1).unwrap_or(u32::MAX),
            );
        }
        let level = (index - NEXACT) / NLEVEL + 1;
        let bin = index - NEXACT - NLEVEL * (level - 1);
        let step = 10_u64.pow(u32::try_from(level).unwrap_or(0));
        let value = 10_u64.pow(u32::try_from(level + 4).unwrap_or(63)) + step * bin as u64;
        (
            u32::try_from(value).unwrap_or(u32::MAX),
            u32::try_from(value + step).unwrap_or(u32::MAX),
        )
    }

    /// Absorbs one QUAL's iqual (`0` = missing/NaN/negative).
    pub fn push_iqual(&mut self, iqual: u32) {
        *self.bins.entry(Self::ibin(iqual)).or_insert(0) += 1;
    }

    /// Absorbs one QUAL (`None` = missing/NaN/negative).
    pub fn push(&mut self, qual: Option<f32>) {
        let iqual = match qual {
            Some(q) if q >= 0.0 => 1 + (q * 10.0) as u32,
            _ => 0,
        };
        self.push_iqual(iqual);
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        for (bin, n) in &other.bins {
            *self.bins.entry(*bin).or_insert(0) += n;
        }
    }

    /// Count at a bin index (`0` beyond the inserted range).
    #[must_use]
    pub fn count_at(&self, index: usize) -> u64 {
        self.bins.get(&index).copied().unwrap_or(0)
    }

    /// Number of bins (`nbins` = highest inserted index + 1, like
    /// bcftools's realloc-grown array).
    #[must_use]
    pub fn nbins(&self) -> usize {
        self.bins.keys().next_back().map_or(0, |k| k + 1)
    }

    /// Total observations.
    #[must_use]
    pub fn total(&self) -> u64 {
        self.bins.values().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn int_histogram_counts_and_quantiles() {
        let mut h = IntHistogram::new(10);
        for v in [0_u32, 1, 2, 3, 10, 11, 100] {
            h.push(v);
        }
        assert_eq!(h.total(), 7);
        // bins[v] holds v; overflow (11, 100) → the last bucket.
        assert_eq!(h.bins[10], 1); // value 10
        assert_eq!(h.bins[11], 2); // overflow
        assert_eq!(h.quantile(1.0), Some(11));
        assert_eq!(h.quantile(0.5), Some(3));
        assert_eq!(h.quantile(0.0), Some(0));
        // Merge is additive.
        let mut a = IntHistogram::new(10);
        a.push(5);
        let mut b = IntHistogram::new(10);
        b.push(5);
        b.push(6);
        a.merge(&b);
        assert_eq!(a.bins[5], 2);
        assert_eq!(a.bins[6], 1);
    }

    #[test]
    fn qual_dist_matches_bcftools_semantics() {
        let mut h = QualDist::new();
        h.push(None);
        h.push(Some(-1.0));
        h.push(Some(0.05_f32));
        h.push(Some(50.0_f32));
        h.push(Some(123.45_f32));
        h.push(Some(9999.95_f32));
        h.push(Some(10001.01_f32));
        // Exact range: iqual 1 (0.0, incl. missing), 2 (0.1), 501 (50.0),
        // 1235 (123.4).
        assert_eq!(h.count_at(0), 2);
        assert_eq!(h.count_at(1), 1);
        assert_eq!(h.count_at(501), 1);
        assert_eq!(h.count_at(1235), 1);
        // Aggregate range: 9999.95 → iqual 100000 → bin 100000 (spanning
        // [10000.0, 10001.0) QUAL); 10001.01 → bin 100001.
        assert_eq!(h.count_at(100_000), 1);
        assert_eq!(h.count_at(100_001), 1);
        assert_eq!(QualDist::bin_range(501), (501, 502));
        assert_eq!(QualDist::bin_range(100_000), (100_000, 100_010));
        // The printed QUAL is 0.1×(beg−1): 50.0 and 9999.9.
        let (beg, _) = QualDist::bin_range(100_000);
        assert_eq!(format!("{:.1}", 0.1 * f64::from(beg - 1)), "9999.9");
        assert_eq!(h.total(), 7);
    }
}
