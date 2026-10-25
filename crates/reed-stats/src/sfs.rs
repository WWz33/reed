//! Site frequency spectrum over complete-case sites
//! (`docs/spec/architecture.md` §9.4, the oracle's `CompleteCase` semantics;
//! the projected estimator is M7).

// Casts in this file are bounded by construction: Sample counts and SFS bins fit the accumulator dimensions.
#![allow(clippy::cast_possible_truncation)]

use reed_core::count::SiteGroupCounts;
use reed_core::math::harmonic::{a1, a2};
use reed_core::math::sum::NeumaierSum;

/// Unfolded (or folded) SFS accumulator for one group. Complete-case: any
/// missing allele in the group skips the site (counted separately).
#[derive(Clone, Debug)]
pub struct SfsAcc {
    group: usize,
    n_chrom: u64,
    folded: bool,
    bins: Vec<u64>,
    n_sites: u64,
    n_skipped_missing: u64,
}

impl SfsAcc {
    /// Accumulator for `group`. Errors when `n_chrom < 2` (no spectrum).
    ///
    /// # Errors
    ///
    /// [`reed_core::CoreError::InvalidWindow`] when the group has fewer than
    /// two chromosomes.
    pub fn new(group: usize, n_group: u32, folded: bool) -> Result<Self, reed_core::CoreError> {
        let n_chrom = u64::from(n_group) * 2;
        if n_chrom < 2 {
            return Err(reed_core::CoreError::InvalidWindow);
        }
        let bins = if folded {
            (n_chrom / 2) as usize
        } else {
            n_chrom.saturating_sub(1) as usize
        };
        Ok(Self {
            group,
            n_chrom,
            folded,
            bins: vec![0; bins],
            n_sites: 0,
            n_skipped_missing: 0,
        })
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let c = &row[self.group];
        let n_group = self.n_chrom / 2;
        if c.n_missing_alleles(n_group as u32) != 0 {
            self.n_skipped_missing += 1;
            return;
        }
        self.n_sites += 1;
        let alt = u64::from(c.n_alt_alleles);
        if alt > 0 && alt < self.n_chrom {
            let index = if self.folded {
                alt.min(self.n_chrom - alt)
            } else {
                alt
            };
            self.bins[(index - 1) as usize] += 1;
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        for (dst, src) in self.bins.iter_mut().zip(&other.bins) {
            *dst += src;
        }
        self.n_sites += other.n_sites;
        self.n_skipped_missing += other.n_skipped_missing;
    }

    /// Bin counts (index `i` = allele count `i+1`, or its minor count when
    /// folded).
    #[must_use]
    pub fn bins(&self) -> &[u64] {
        &self.bins
    }

    /// Contributing / skipped site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped_missing)
    }

    /// Chromosome count the spectrum is defined over.
    #[must_use]
    pub const fn n_chrom(&self) -> u64 {
        self.n_chrom
    }
}

// ---------------------------------------------------------------- projected

/// Hypergeometric projection weights for one site: `P(k | n, m, alt) =
/// C(alt,k)·C(m−alt, n−k) / C(m, n)` (Marth 2004; dadi `project`), built
/// with the ratio recursion `P(k+1)/P(k) = (alt−k)(n−k) /
/// ((k+1)(m−alt−n+k+1))` anchored at the distribution's mode and run
/// outward in both directions, normalized, and truncated at a relative
/// weight of 1e−15 (§9.4). The mode anchor keeps every running weight ≤ 1
/// (doc_v20 R7: the old left-edge anchor multiplied up-hill and hit inf
/// for n ≳ 1024 at mid p, silently dropping the site's mass). Returns
/// `(k, weight)` pairs over the **full** support `k ∈ 0..=n`: the
/// non-segregating mass at k = 0/n is part of the distribution and is what
/// the dadi/fastsimcoal2 writers report
/// (docs/rounds/v5/05-algorithms.md §3.1).
#[must_use]
pub fn projection_weights(m: u32, alt: u32, n: u32) -> Vec<(u32, f64)> {
    let (m, alt, n) = (u64::from(m), u64::from(alt), u64::from(n));
    // Full support, including k = 0 and k = n: the projection is a true
    // probability distribution, and its non-segregating mass (k = 0/n)
    // is simply not stored in the spectrum — never renormalized away.
    let k_lo = n.saturating_sub(m - alt);
    let k_hi = n.min(alt);
    if k_lo > k_hi {
        return Vec::new();
    }
    let mut out = Vec::with_capacity((k_hi - k_lo + 1) as usize);
    // Anchor the recursion at the mode k* = ⌊(alt+1)(n+1)/(m+2)⌋, clamped
    // into the support (doc_v20 R7): every ratio stepping away from the
    // mode is ≤ 1, so the running weights never exceed 1, where the old
    // left-edge anchor climbed up-hill through a peak that grows like
    // e^{n·ln 2} and hit inf for --proj ≳ 1024 at mid p. The quotient is
    // evaluated in u128 — (alt+1)(n+1) overflows u64 at large m.
    let k_star = u64::try_from(
        ((alt + 1) as u128 * (n + 1) as u128 / (m as u128 + 2)).clamp(k_lo as u128, k_hi as u128),
    )
    .unwrap_or(k_hi);
    // Downward from the mode: P(k−1)/P(k) = k(m−alt−n+k) / ((alt−k+1)(n−k+1))
    // — the factored inverse of the upward ratio below, so it is one small
    // quotient per step, never a quotient of two accumulated weights.
    let mut down = Vec::new();
    let mut w = 1.0_f64;
    let mut k = k_star;
    while k > k_lo {
        // The denominator's inner difference goes through i64 — the u64
        // intermediate (m−alt)−n can be negative before k is added.
        let num = k as f64 * (m as i64 - alt as i64 - n as i64 + k as i64) as f64;
        let den = (alt - k + 1) as f64 * (n - k + 1) as f64;
        w *= num / den;
        k -= 1;
        down.push((u32::try_from(k).unwrap_or(u32::MAX), w));
    }
    down.reverse();
    out.append(&mut down);
    out.push((u32::try_from(k_star).unwrap_or(u32::MAX), 1.0));
    // Upward from the mode: P(k+1)/P(k) = (alt−k)(n−k) / ((k+1)(m−alt−n+k+1)).
    let mut w = 1.0_f64;
    let mut k = k_star;
    while k < k_hi {
        // Same i64 note as the downward loop: the denominator's inner
        // difference can be negative in u64 before k is added.
        let num = (alt - k) as f64 * (n - k) as f64;
        let den = (k + 1) as f64 * (m as i64 - alt as i64 - n as i64 + k as i64 + 1) as f64;
        w *= num / den;
        k += 1;
        out.push((u32::try_from(k).unwrap_or(u32::MAX), w));
    }
    let sum: f64 = out.iter().map(|(_, w)| *w).sum();
    if sum <= 0.0 || !sum.is_finite() {
        return Vec::new();
    }
    for (_, w) in &mut out {
        *w /= sum;
    }
    // Numerical guard (§9.4): drop weights below 1e−15 of the PEAK
    // without renormalizing — the distribution is unimodal, so the dropped
    // tails sit at the recursion's noise floor and their mass (≤ 1e−15 ×
    // count) is negligible.
    let peak = out.iter().map(|(_, w)| *w).fold(0.0_f64, f64::max);
    out.retain(|(_, w)| *w >= 1e-15 * peak);
    out
}

/// Worker-side memo for projection weights, keyed by `(m, alt)` for a
/// fixed `n` (§9.4: LRU; TODO: clear-on-full like `HweCache`).
/// A shared weight slice (the cache hands out `Arc`s so sites reuse one
/// allocation).
pub type WeightsArc = std::sync::Arc<Vec<(u32, f64)>>;

/// Per-`n` projection-weight memo (§9.4): cleared wholesale whenever it
/// reaches [`PROJECTION_CACHE_CAPACITY`] entries.
#[derive(Clone, Debug, Default)]
pub struct ProjectionCache {
    n: u32,
    map: std::collections::HashMap<(u32, u32), WeightsArc>,
}

/// Capacity of [`ProjectionCache`] before a wholesale clear.
pub const PROJECTION_CACHE_CAPACITY: usize = 8192;

impl ProjectionCache {
    /// A cache for a fixed projection size `n`.
    #[must_use]
    pub fn new(n: u32) -> Self {
        Self {
            n,
            map: std::collections::HashMap::new(),
        }
    }

    /// The projection size.
    #[must_use]
    pub const fn n(&self) -> u32 {
        self.n
    }

    /// Weights for `(m, alt)`, cached.
    #[must_use]
    pub fn weights(&mut self, m: u32, alt: u32) -> WeightsArc {
        if let Some(w) = self.map.get(&(m, alt)) {
            return std::sync::Arc::clone(w);
        }
        let w: WeightsArc = std::sync::Arc::new(projection_weights(m, alt, self.n));
        if self.map.len() >= PROJECTION_CACHE_CAPACITY {
            self.map.clear();
        }
        self.map.insert((m, alt), std::sync::Arc::clone(&w));
        w
    }
}

/// Projected site-frequency spectrum for one group (§9.4 `Projected { n }`):
/// sites with `m ≥ n` called chromosomes distribute their alt count over
/// bins `ξ_k` by hypergeometric weights; `m < n` sites are skipped and
/// counted. Float bins are [`NeumaierSum`]s so chunk merges stay within the
/// §8.4 tolerance.
///
/// `bins` spans the **full** support `k ∈ 0..=n` (folded: `0..=n/2`), so
/// the dadi/fastsimcoal2 writers can report the non-segregating mass
/// (docs/rounds/v5/05-algorithms.md §3.1); the ratio statistics
/// ([`Self::pi`], [`Self::segregating`]) still use only `1..n`.
#[derive(Clone, Debug)]
pub struct ProjectedSfsAcc {
    n: u32,
    folded: bool,
    /// `bins[k]` = ξ_k for `k` in `0..=n` (folded: `0..=n/2`).
    bins: Vec<NeumaierSum>,
    n_sites: u64,
    n_skipped_m_lt_n: u64,
}

impl ProjectedSfsAcc {
    /// Accumulator projected to `n` chromosomes; `None` when `n < 2`.
    #[must_use]
    pub fn new(n: u32, folded: bool) -> Option<Self> {
        (n >= 2).then(|| Self {
            n,
            folded,
            bins: vec![NeumaierSum::new(); n as usize + 1],
            n_sites: 0,
            n_skipped_m_lt_n: 0,
        })
    }

    /// The projection size.
    #[must_use]
    pub const fn n(&self) -> u32 {
        self.n
    }

    /// Absorbs one site with pre-fetched weights (`None` = the site was
    /// skipped for `m < n`; the caller decides via its cache).
    pub fn push_weighted(&mut self, weights: Option<&[(u32, f64)]>) {
        if weights.is_none() {
            self.n_skipped_m_lt_n += 1;
            return;
        }
        self.n_sites += 1;
        for &(k, w) in weights.unwrap_or(&[]) {
            // Folded spectra put k and n−k in one bin; the whole support
            // (including k = 0/n) is stored.
            let k = if self.folded { k.min(self.n - k) } else { k };
            if k <= self.n {
                self.bins[k as usize].add(w);
            }
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        debug_assert_eq!(self.n, other.n, "invariant: same projection size");
        debug_assert_eq!(self.folded, other.folded, "invariant: same folding");
        for (dst, src) in self.bins.iter_mut().zip(&other.bins) {
            dst.merge(src);
        }
        self.n_sites += other.n_sites;
        self.n_skipped_m_lt_n += other.n_skipped_m_lt_n;
    }

    /// Bin ξ_k (k in 0..=n; folded: 0..=n/2).
    #[must_use]
    pub fn xi(&self, k: u32) -> f64 {
        self.bins[k as usize].value()
    }

    /// All bins as a slice view (index = k, length n+1).
    #[must_use]
    pub fn bins(&self) -> &[NeumaierSum] {
        &self.bins
    }

    /// Whether the spectrum is folded.
    #[must_use]
    pub const fn folded(&self) -> bool {
        self.folded
    }

    /// The spectrum as the dadi/fastsimcoal2 array: `k = 0..=n` (unfolded)
    /// or `k = 0..=n/2` (folded), docs/rounds/v5/05-algorithms.md §3.1.
    #[must_use]
    pub fn spectrum(&self) -> Vec<f64> {
        let top = if self.folded { self.n / 2 } else { self.n };
        (0..=top).map(|k| self.bins[k as usize].value()).collect()
    }

    /// θ_π = Σ_k k(n−k)·ξ_k / C(n,2).
    #[must_use]
    pub fn pi(&self) -> f64 {
        let n = u64::from(self.n);
        let num: f64 = (1..self.n)
            .map(|k| {
                let k = u64::from(k);
                (k * (n - k)) as f64 * self.bins[k as usize].value()
            })
            .sum();
        let comps = n * (n - 1) / 2;
        num / comps as f64
    }

    /// S = Σ_{k=1}^{n−1} ξ_k (segregating sites on the projected spectrum).
    #[must_use]
    pub fn segregating(&self) -> f64 {
        (1..self.n as usize).map(|k| self.bins[k].value()).sum()
    }

    /// θ_W = S / a₁(n) on the projected spectrum.
    #[must_use]
    pub fn theta_w(&self) -> Option<f64> {
        let n = usize::try_from(self.n).ok()?;
        (n >= 2).then(|| self.segregating() / a1(n))
    }

    /// Tajima's D on the projected spectrum (fixed n: the standard moment
    /// variance applies; `min_sites` guards S).
    #[must_use]
    pub fn tajima_d(&self, min_sites: u64) -> Option<f64> {
        let n = u64::from(self.n);
        if n < 2 {
            return None;
        }
        let n_us = usize::try_from(n).ok()?;
        let a1n = a1(n_us);
        let a2n = a2(n_us);
        let s = self.segregating();
        if s < min_sites as f64 {
            return None;
        }
        let theta_w = if a1n > 0.0 { s / a1n } else { 0.0 };
        let nf = n as f64;
        let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
        let b2 = 2.0 * (nf * nf + nf + 3.0) / (9.0 * nf * (nf - 1.0));
        let c1 = b1 - 1.0 / a1n;
        let c2 = b2 - (nf + 2.0) / (a1n * nf) + a2n / (a1n * a1n);
        let e1 = c1 / a1n;
        let e2 = c2 / (a1n * a1n + a2n);
        let variance = e1 * s + e2 * s * (s - 1.0);
        (variance > 0.0).then(|| (self.pi() - theta_w) / variance.sqrt())
    }

    /// Site tallies: projected / skipped (m < n).
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped_m_lt_n)
    }
}

// ------------------------------------------------------------------- joint

// ------------------------------------------------------------- joint projected

/// The projected joint (2D) spectrum: each site's per-group 1D Marth weights
/// combine as an outer product, and folding follows dadi's **global** rule
/// (`Spectrum.fold`: `i + j` against `total/2`, with the diagonal split
/// 0.5/0.5) — not a per-axis minor projection. Layout matches the dadi 2D
/// data line (row-major `i * (n2 + 1) + j`).
#[derive(Clone, Debug)]
pub struct JointProjectedAcc {
    n1: u32,
    n2: u32,
    folded: bool,
    bins: Vec<NeumaierSum>,
    n_sites: u64,
    n_skipped: u64,
}

impl JointProjectedAcc {
    /// Accumulator projected to `(n1, n2)`; `None` when either is below 2.
    #[must_use]
    pub fn new(n1: u32, n2: u32, folded: bool) -> Option<Self> {
        (n1 >= 2 && n2 >= 2).then(|| Self {
            n1,
            n2,
            folded,
            bins: vec![NeumaierSum::new(); (n1 as usize + 1) * (n2 as usize + 1)],
            n_sites: 0,
            n_skipped: 0,
        })
    }

    /// The projection sizes.
    #[must_use]
    pub const fn shape(&self) -> (u32, u32) {
        (self.n1, self.n2)
    }

    /// Whether the spectrum is folded.
    #[must_use]
    pub const fn folded(&self) -> bool {
        self.folded
    }

    /// Absorbs one site from the two groups' pre-fetched 1D weight vectors;
    /// `None` on either side = that group's `m < n`, so the site is skipped
    /// (easySFS/dadi semantics).
    pub fn push_weighted(&mut self, w1: Option<&[(u32, f64)]>, w2: Option<&[(u32, f64)]>) {
        let (Some(w1), Some(w2)) = (w1, w2) else {
            self.n_skipped += 1;
            return;
        };
        self.n_sites += 1;
        let stride = self.n2 as usize + 1;
        for &(i, a) in w1 {
            for &(j, b) in w2 {
                let w = a * b;
                if self.folded {
                    self.add_folded(i, j, w, stride);
                } else {
                    self.bins[i as usize * stride + j as usize].add(w);
                }
            }
        }
    }

    /// dadi's `Spectrum.fold`: mirror entries above the anti-diagonal, split
    /// the diagonal (only when `n1 + n2` is even) 0.5/0.5 between mirrors.
    fn add_folded(&mut self, i: u32, j: u32, w: f64, stride: usize) {
        let total = self.n1 + self.n2;
        let here = i as usize * stride + j as usize;
        let mirror = (self.n1 - i) as usize * stride + (self.n2 - j) as usize;
        match (i + j) * 2 {
            s if s < total => self.bins[here].add(w),
            s if s == total => {
                self.bins[here].add(0.5 * w);
                self.bins[mirror].add(0.5 * w);
            }
            _ => self.bins[mirror].add(w),
        }
    }

    fn index(&self, i: u32, j: u32) -> usize {
        i as usize * (self.n2 as usize + 1) + j as usize
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        debug_assert_eq!(
            (self.n1, self.n2, self.folded),
            (other.n1, other.n2, other.folded),
            "invariant: same shape and folding"
        );
        for (dst, src) in self.bins.iter_mut().zip(&other.bins) {
            dst.merge(src);
        }
        self.n_sites += other.n_sites;
        self.n_skipped += other.n_skipped;
    }

    /// Bin (i, j).
    #[must_use]
    pub fn at(&self, i: u32, j: u32) -> f64 {
        self.bins[self.index(i, j)].value()
    }

    /// Site tallies: absorbed / skipped (either group `m < n`).
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped)
    }

    /// The spectrum in the dadi 2D data-line order (row-major).
    #[must_use]
    pub fn values(&self) -> Vec<f64> {
        self.bins.iter().map(NeumaierSum::value).collect()
    }
}
