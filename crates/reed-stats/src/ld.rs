//! Linkage disequilibrium on genotype planes (§9.7): integer Rogers-Huff
//! r² between unphased genotype dosage pairs, computed with popcounts over
//! the shared-valid mask — i64 sums throughout, the single final step in
//! f64. The two-stage topology (plane-store a contig, then pair blocks that
//! look back `max_bp`/`max_sites`) is the driver's job; this module is the
//! pure math over contiguous plane slices.

/// Integer sums for one (x, y) site pair: enough to form r² exactly.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PairSums {
    /// Valid samples: full_x ∧ full_y.
    pub n: i64,
    /// Σx over valid samples (dosage 0/1/2 per sample).
    pub sum_x: i64,
    /// Σx² = Σx + 2·(slot0_alt ∧ slot1_alt) — a² = a for bits.
    pub sum_x2: i64,
    /// Σy over valid samples.
    pub sum_y: i64,
    /// Σy².
    pub sum_y2: i64,
    /// Σxy over valid samples.
    pub sum_xy: i64,
}

/// Computes the pair sums from two sites' planes (§9.7's twelve
/// popcounts per pair, shared mask first).
///
/// The `x`/`y` arrays are each site's four planes
/// `[slot0_alt, slot0_miss, slot1_alt, slot1_miss]` as `W`-word slices;
/// `n_samples` bounds the real samples (the word padding past it is not
/// part of the population and never enters any mask).
#[must_use]
pub fn pair_sums(x: [&[u64]; 4], y: [&[u64]; 4], n_samples: usize) -> PairSums {
    let mut s = PairSums::default();
    for w in 0..x[0].len().min(n_samples.div_ceil(64)) {
        let (x0a, x0m, x1a, x1m) = (x[0][w], x[1][w], x[2][w], x[3][w]);
        let (y0a, y0m, y1a, y1m) = (y[0][w], y[1][w], y[2][w], y[3][w]);
        // Real-sample mask for this word (padding bits excluded).
        let real = if (w + 1) * 64 <= n_samples {
            !0_u64
        } else {
            (1_u64 << (n_samples % 64)) - 1
        };
        // Valid: both slots of both sites called.
        let vx = real & !(x0m | x1m);
        let vy = real & !(y0m | y1m);
        let v = vx & vy;
        let n = v.count_ones() as i64;
        if n == 0 {
            continue;
        }
        // Dosage sums: alt bits of each slot, masked valid.
        let sx = (v & x0a).count_ones() as i64 + (v & x1a).count_ones() as i64;
        let sy = (v & y0a).count_ones() as i64 + (v & y1a).count_ones() as i64;
        // Σx² = Σx + 2·(both slots alt): a bit squared is itself, doubled
        // for the 1+1 case.
        let both_x = (v & x0a & x1a).count_ones() as i64;
        let both_y = (v & y0a & y1a).count_ones() as i64;
        // Σxy: dosage products over the four slot combinations.
        let sxy = (v & x0a & y0a).count_ones() as i64
            + (v & x0a & y1a).count_ones() as i64
            + (v & x1a & y0a).count_ones() as i64
            + (v & x1a & y1a).count_ones() as i64;
        s.n += n;
        s.sum_x += sx;
        s.sum_x2 += sx + 2 * both_x;
        s.sum_y += sy;
        s.sum_y2 += sy + 2 * both_y;
        s.sum_xy += sxy;
    }
    s
}

impl PairSums {
    /// Rogers-Huff r² over the valid samples: the squared Pearson
    /// correlation of the dosage vectors. The cross-products run through
    /// i128 (exact for any real sample count), one f64 step at the very
    /// end (§9.7). `None` when the variance is zero on either side or
    /// fewer than 2 valid samples.
    #[must_use]
    pub fn r2(&self) -> Option<f64> {
        if self.n < 2 {
            return None;
        }
        let n = i128::from(self.n);
        let cov = i128::from(self.sum_xy) * n - i128::from(self.sum_x) * i128::from(self.sum_y);
        let var_x = i128::from(self.sum_x2) * n - i128::from(self.sum_x) * i128::from(self.sum_x);
        let var_y = i128::from(self.sum_y2) * n - i128::from(self.sum_y) * i128::from(self.sum_y);
        if var_x <= 0 || var_y <= 0 {
            return None;
        }
        let num = cov as f64 * cov as f64;
        let den = var_x as f64 * var_y as f64;
        Some(num / den)
    }

    /// The covariance E[xy] − E[x]·E[y] over the valid samples (the
    /// sample-size-n denominator applies — Rogers-Huff D).
    #[must_use]
    pub fn covariance(&self) -> Option<f64> {
        if self.n == 0 {
            return None;
        }
        let n = i128::from(self.n);
        let cov = i128::from(self.sum_xy) * n - i128::from(self.sum_x) * i128::from(self.sum_y);
        Some(cov as f64 / (n * n) as f64)
    }
}

/// LD-decay bins: Σr² and counts by distance (additive; §9.7's decay
/// output). Distances are site-index distances or bp — the caller decides
/// and supplies the bin edges.
#[derive(Clone, Debug)]
pub struct LdDecayAcc {
    /// Inclusive lower edges (bins: `[edges[i], edges[i+1])`, the last
    /// closed).
    edges: Vec<u64>,
    sum_r2: Vec<reed_core::math::sum::NeumaierSum>,
    counts: Vec<u64>,
}

impl LdDecayAcc {
    /// Bins over the given ascending edges.
    #[must_use]
    pub fn new(edges: Vec<u64>) -> Self {
        let n = edges.len().max(1);
        Self {
            edges,
            sum_r2: (0..n).map(|_| reed_core::math::sum::NeumaierSum::new()).collect(),
            counts: vec![0; n],
        }
    }

    fn bin_of(&self, distance: u64) -> usize {
        // Last bin with edge <= distance.
        match self.edges.binary_search(&distance) {
            Ok(i) => i,
            Err(i) => i.saturating_sub(1),
        }
    }

    /// Absorbs one pair's (distance, r²).
    pub fn push(&mut self, distance: u64, r2: f64) {
        let b = self.bin_of(distance);
        if b < self.counts.len() {
            self.sum_r2[b].add(r2);
            self.counts[b] += 1;
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        for (dst, src) in self.sum_r2.iter_mut().zip(&other.sum_r2) {
            dst.merge(src);
        }
        for (dst, src) in self.counts.iter_mut().zip(&other.counts) {
            *dst += src;
        }
    }

    /// Mean r² per bin (`None` for empty bins).
    #[must_use]
    pub fn mean_r2(&self) -> Vec<Option<f64>> {
        self.counts
            .iter()
            .zip(&self.sum_r2)
            .map(|(&c, s)| (c > 0).then(|| s.value() / c as f64))
            .collect()
    }

    /// Counts per bin.
    #[must_use]
    pub fn counts(&self) -> &[u64] {
        &self.counts
    }
}

/// Ordered-carry LD pruning (plink2 semantics): sites are consumed in
/// order; a site is KEPT when its r² with every already-kept site within
/// `max_sites` is strictly below `threshold`; otherwise it is pruned.
#[derive(Clone, Debug)]
pub struct LdPruner {
    threshold: f64,
    max_sites: usize,
    /// The ring of kept sites' (site, planes) — bounded by max_sites.
    kept: std::collections::VecDeque<(usize, Vec<u64>)>,
    /// Per kept site: the four planes (W words each, flattened).
    w: usize,
    n_samples: usize,
}

impl LdPruner {
    /// A pruner at `threshold` (r²) looking back `max_sites` kept sites.
    #[must_use]
    pub fn new(threshold: f64, max_sites: usize, w: usize, n_samples: usize) -> Self {
        Self {
            threshold,
            max_sites,
            kept: std::collections::VecDeque::new(),
            w,
            n_samples,
        }
    }

    /// Whether this site is kept given its four planes; absorbs it into
    /// the kept set when it is.
    pub fn push(&mut self, site: usize, planes: [&[u64]; 4]) -> bool {
        for (kept_site, kept_planes) in &self.kept {
            let kp: [&[u64]; 4] = [
                &kept_planes[0..self.w],
                &kept_planes[self.w..2 * self.w],
                &kept_planes[2 * self.w..3 * self.w],
                &kept_planes[3 * self.w..4 * self.w],
            ];
            if let Some(r2) = pair_sums(planes, kp, self.n_samples).r2()
                && r2 >= self.threshold
            {
                return false;
            }
            let _ = kept_site;
        }
        let mut flat = Vec::with_capacity(4 * self.w);
        for p in planes {
            flat.extend_from_slice(p);
        }
        self.kept.push_back((site, flat));
        if self.kept.len() > self.max_sites {
            self.kept.pop_front();
        }
        true
    }

    /// The kept site indices in consumption order.
    #[must_use]
    pub fn kept_sites(&self) -> Vec<usize> {
        self.kept.iter().map(|(s, _)| *s).collect()
    }
}
