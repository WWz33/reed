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

/// The pruning window (plink2 `--indep-pairwise`): a variant count, or a
/// span in bp (`W kb`; plink2 then requires step 1).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LdWindow {
    /// Look back this many kept variants.
    Variants(usize),
    /// Look back this many bp.
    Bp(u64),
}

/// plink2 `--indep-pairwise` pruning (docs/rounds/v5/05-algorithms.md
/// §3.4.1; `plink2_ld.cc` `IndepPairwiseThread`, read 2026-10-10).
///
/// Sites arrive in contig/position order. The window spans `window` (a
/// variant count, or bp) **in arrival space**: every arriving variant
/// occupies a slot, including the monomorphic ones (pruned on arrival, like
/// plink2's `FillVaggs` check) and the ones pruned earlier. After each
/// window the start advances by `step` variants and the pruned slots leave.
///
/// Within a window every surviving pair is tested with the exact integer
/// judgment `cov² > R2 · var1 · var2` (strictly greater, so a pair with no
/// common observations never prunes); when it fires, the variant with the
/// **higher major-allele frequency** is pruned (a tie, within `1 + 1e-6`,
/// prunes the later one). The scan runs nearest-pair-first (plink2's
/// default order).
#[derive(Clone, Debug)]
pub struct LdPruner {
    threshold: f64,
    window: LdWindow,
    step: usize,
    w: usize,
    n_samples: usize,
    /// The window's live slots, in arrival order (pruned ones removed).
    sites: Vec<PruneSite>,
    /// Slots seen since the last contig reset (the arrival space).
    seen: usize,
    /// The arrival index the current window starts at.
    window_start: usize,
    /// Index into `sites` where the current window's NEW arrivals begin:
    /// only those take the `second` role in the scan (plink2's
    /// `winpos_split`), so an old pair is never tested twice.
    new_from: usize,
    /// Global ordinals kept so far.
    kept: Vec<usize>,
    /// Global ordinals pruned.
    pruned: Vec<usize>,
    /// Position of the last arrival (the bp window's right edge).
    last_pos0: u64,
    /// The bp window's start position (0 until the first evaluation).
    window_start_bp: u64,
}

/// The per-site aggregates plink2's `FillVaggs` produces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SitePruneStats {
    /// Major-allele frequency over all samples (`--freq` convention).
    pub maj_freq: f64,
    /// Called homozygotes carrying the reference allele.
    pub n_hom_ref: u64,
    /// Called homozygotes carrying the alternate allele.
    pub n_hom_alt: u64,
    /// Called heterozygotes.
    pub n_het: u64,
}

impl SitePruneStats {
    /// plink2's monomorphic test (`FillVaggs`' three cases): no homozygotes
    /// at all (all-het), all-reference, or all-alternate. Such a site is
    /// pruned on arrival.
    #[must_use]
    pub const fn monomorphic(&self) -> bool {
        let n_hom = self.n_hom_ref + self.n_hom_alt;
        let n_called = n_hom + self.n_het;
        n_hom == 0
            || (self.n_het == 0 && self.n_hom_alt == 0)
            || (self.n_het == 0 && self.n_hom_ref == 0)
            || n_called == 0
    }
}

/// One live window slot.
#[derive(Clone, Debug)]
struct PruneSite {
    ordinal: usize,
    /// The arrival index (window slots are counted in arrival space).
    arrival: usize,
    pos0: u64,
    planes: Vec<u64>,
    maj_freq: f64,
}

impl LdPruner {
    /// A pruner at `threshold` (r²) over `window`, stepping `step` variants
    /// (`step` is read as 1 when 0).
    #[must_use]
    pub fn new(threshold: f64, window: LdWindow, step: usize, w: usize, n_samples: usize) -> Self {
        Self {
            threshold,
            window,
            step: step.max(1),
            w,
            n_samples,
            sites: Vec::new(),
            seen: 0,
            window_start: 0,
            new_from: 0,
            kept: Vec::new(),
            pruned: Vec::new(),
            last_pos0: 0,
            window_start_bp: 0,
        }
    }

    /// Feeds one site: `ordinal` is its global index among the sites this
    /// pruner sees, `stats` its plink2 aggregates. Returns whether it is
    /// kept.
    pub fn push(
        &mut self,
        ordinal: usize,
        pos0: u64,
        planes: [&[u64]; 4],
        stats: SitePruneStats,
    ) -> bool {
        self.seen += 1;
        self.last_pos0 = pos0;
        // plink2 prunes monomorphic variants on arrival: they occupy a slot
        // and never enter the pair scan.
        if stats.monomorphic() {
            self.pruned.push(ordinal);
            self.maybe_evaluate();
            return false;
        }
        let mut flat = Vec::with_capacity(4 * self.w);
        for p in planes {
            flat.extend_from_slice(p);
        }
        self.sites.push(PruneSite {
            ordinal,
            arrival: self.seen - 1,
            pos0,
            planes: flat,
            maj_freq: stats.maj_freq,
        });
        self.maybe_evaluate();
        // The site survives unless this window's scan removed it.
        !self.pruned.contains(&ordinal)
    }

    /// Evaluates the window once its span is complete, then advances it.
    fn maybe_evaluate(&mut self) {
        let full = match self.window {
            LdWindow::Variants(max) => self.seen >= self.window_start + max,
            // plink2's bp window: the end is the first variant whose position
            // exceeds start_bp + W; the window is complete when such a
            // variant arrives.
            LdWindow::Bp(bp) => {
                let start_bp = self.window_start_bp;
                self.last_pos0 > start_bp.saturating_add(bp)
            }
        };
        if full {
            self.evaluate_window();
        }
    }

    /// Runs the window's pairwise scan, then advances the window by `step`.
    fn evaluate_window(&mut self) {
        let n = self.sites.len();
        let mut removed = vec![false; n];
        // Reverse scan over the NEW slots, nearest pairs first (plink2's
        // `winpos_split`: only the sites that entered since the last
        // evaluation take the `second` role).
        let stop = if self.new_from == 0 { 1 } else { self.new_from };
        let mut second = n;
        while second > stop {
            second -= 1;
            // plink2 does NOT skip a removed `second` (only `first` carries
            // the `IsSet(cur_window_removed, first_winpos)` check): a site
            // pruned earlier in this window still pairs as the second
            // member, which is what lets a later pair prune an older site.
            let mut first = second;
            while first > 0 {
                first -= 1;
                if removed[first] {
                    continue;
                }
                if self.prunes(first, second) {
                    if self.sites[first].maj_freq <= self.sites[second].maj_freq * (1.0 + 1e-6) {
                        removed[second] = true;
                        break;
                    }
                    removed[first] = true;
                }
            }
        }
        // Record the decisions and drop every decided slot (kept ones leave
        // the window too: they are behind the next window's start). The
        // slot's arrival index — not its position in the compacted vec —
        // decides which window it belongs to.
        let next_start = self.window_start + self.step;
        let mut live = Vec::with_capacity(n);
        for (i, site) in self.sites.drain(..).enumerate() {
            if removed[i] {
                self.pruned.push(site.ordinal);
                continue;
            }
            if site.arrival < next_start {
                self.kept.push(site.ordinal);
                continue;
            }
            live.push(site);
        }
        self.sites = live;
        self.new_from = self.sites.len();
        self.window_start = next_start;
        // A bp window's start advances until the end can move by at least
        // one variant (plink2's `LdPruneNextWindow`: the new start is the
        // first variant whose position is >= end_bp − W).
        if let LdWindow::Bp(bp) = self.window {
            let floor = self.last_pos0.saturating_sub(bp);
            let first_live = self
                .sites
                .iter()
                .find(|s| s.pos0 >= floor)
                .map_or(self.last_pos0, |s| s.pos0);
            self.window_start_bp = first_live;
        }
    }

    /// plink2's exact judgment: `cov² > R2 · var1 · var2` on integer
    /// cross-products (a pair with zero common observations never prunes).
    fn prunes(&self, first: usize, second: usize) -> bool {
        fn planes_of(s: &PruneSite, w: usize) -> [&[u64]; 4] {
            [
                &s.planes[0..w],
                &s.planes[w..2 * w],
                &s.planes[2 * w..3 * w],
                &s.planes[3 * w..4 * w],
            ]
        }
        let sums = pair_sums(
            planes_of(&self.sites[first], self.w),
            planes_of(&self.sites[second], self.w),
            self.n_samples,
        );
        let n = i128::from(sums.n);
        let cov = n * i128::from(sums.sum_xy) - i128::from(sums.sum_x) * i128::from(sums.sum_y);
        let var1 = n * i128::from(sums.sum_x2) - i128::from(sums.sum_x) * i128::from(sums.sum_x);
        let var2 = n * i128::from(sums.sum_y2) - i128::from(sums.sum_y) * i128::from(sums.sum_y);
        let lhs = (cov as f64) * (cov as f64);
        let rhs = self.threshold * (var1 as f64) * (var2 as f64);
        lhs > rhs
    }

    /// Ends the current contig: the trailing partial window is scanned
    /// (plink2 evaluates it when the subcontig's last variant arrives), then
    /// every remaining live slot is kept.
    pub fn finish_contig(&mut self) {
        if !self.sites.is_empty() {
            self.evaluate_window();
        }
        for site in self.sites.drain(..) {
            self.kept.push(site.ordinal);
        }
        self.seen = 0;
        self.window_start = 0;
        self.new_from = 0;
        self.window_start_bp = 0;
    }

    /// The kept global ordinals, ascending.
    #[must_use]
    pub fn kept_sites(&self) -> Vec<usize> {
        let mut out = self.kept.clone();
        out.sort_unstable();
        out
    }

    /// The pruned global ordinals, ascending.
    #[must_use]
    pub fn pruned_sites(&self) -> Vec<usize> {
        let mut out = self.pruned.clone();
        out.sort_unstable();
        out
    }
}
