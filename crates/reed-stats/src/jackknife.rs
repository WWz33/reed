//! Block jackknife (§8.6): delete-one-block standard errors for any
//! statistic that reduces to per-site additive (numerator, denominator)
//! sums. Blocks live on the *global site ordinal* grid — the sink-side
//! accumulator numbers sites as chunk partials arrive in order, so the
//! block partition (and therefore the SE) is independent of the chunk
//! partition and thread count (§7.6). Blocks are fixed site counts
//! (equal-sized), so the plain delete-1 jackknife applies. The physical-bp
//! implementation below instead uses Busing/ADMIXTOOLS weighting (V27-35).

use reed_core::math::sum::NeumaierSum;

/// A chunk partial: per-site `(numerator, denominator)` values in site
/// order. Small (two f64 per site) and trivially additive at the sink.
#[derive(Clone, Debug, Default)]
pub struct BlockPartial {
    /// `(numerator, denominator)` per site, in site order.
    pub values: Vec<(f64, f64)>,
}

/// The sink-side jackknife accumulator: assigns each arriving site its
/// global ordinal, aggregates into blocks of `block_size`, and keeps the
/// delete-one-block machinery.
#[derive(Clone, Debug)]
pub struct BlockJackknife {
    block_size: u64,
    next_ordinal: u64,
    /// The open (not yet full) block: (global block index, sums, count).
    open: Option<(u64, NeumaierSum, NeumaierSum, u64)>,
    blocks: std::collections::BTreeMap<u64, (NeumaierSum, NeumaierSum, u64)>,
    total: (NeumaierSum, NeumaierSum),
    n_sites: u64,
}

impl BlockJackknife {
    /// A jackknife over blocks of `block_size` sites.
    #[must_use]
    pub fn new(block_size: u64) -> Self {
        Self {
            block_size: block_size.max(1),
            next_ordinal: 0,
            open: None,
            blocks: std::collections::BTreeMap::new(),
            total: (NeumaierSum::new(), NeumaierSum::new()),
            n_sites: 0,
        }
    }

    /// Absorbs one site's (num, den) at the next global ordinal.
    pub fn push(&mut self, num: f64, den: f64) {
        let idx = self.next_ordinal / self.block_size;
        self.next_ordinal += 1;
        self.total.0.add(num);
        self.total.1.add(den);
        self.n_sites += 1;
        match &mut self.open {
            Some((oidx, num_sum, den_sum, n)) if *oidx == idx => {
                num_sum.add(num);
                den_sum.add(den);
                *n += 1;
                if *n >= self.block_size {
                    let (_, num_sum, den_sum, n) = self.open.take().expect("invariant: open");
                    self.blocks.insert(idx, (num_sum, den_sum, n));
                }
            }
            _ => {
                let mut num_sum = NeumaierSum::new();
                num_sum.add(num);
                let mut den_sum = NeumaierSum::new();
                den_sum.add(den);
                // A brand-new block can already be full (block_size == 1).
                if self.block_size <= 1 {
                    self.blocks.insert(idx, (num_sum, den_sum, 1));
                } else {
                    self.open = Some((idx, num_sum, den_sum, 1));
                }
            }
        }
    }

    /// Absorbs a chunk partial site by site (ordered sink semantics).
    pub fn accept(&mut self, partial: &BlockPartial) {
        for (num, den) in &partial.values {
            self.push(*num, *den);
        }
    }

    /// The point estimate Σnum / Σden.
    #[must_use]
    pub fn estimate(&self) -> Option<f64> {
        let den = self.total.1.value();
        (self.n_sites > 0 && den != 0.0).then(|| self.total.0.value() / den)
    }

    /// (Σnum, Σden) — the inputs of a ratio statistic.
    #[must_use]
    pub const fn totals(&self) -> (&NeumaierSum, &NeumaierSum) {
        (&self.total.0, &self.total.1)
    }

    /// The delete-one-block jackknife: (estimate, se, per-block delete-one
    /// values). `None` with fewer than 2 blocks.
    #[must_use]
    pub fn jackknife(&self) -> Option<(f64, f64, Vec<f64>)> {
        let n_blocks = self.blocks.len() + usize::from(self.open.is_some());
        if n_blocks < 2 {
            return None;
        }
        let total_num = self.total.0.value();
        let total_den = self.total.1.value();
        if total_den == 0.0 {
            return None;
        }
        let mut vj = Vec::with_capacity(n_blocks);
        let iter_blocks = || -> Vec<(NeumaierSum, NeumaierSum, u64)> {
            let mut out: Vec<(NeumaierSum, NeumaierSum, u64)> =
                self.blocks.values().copied().collect();
            if let Some((_, n, d, c)) = &self.open {
                out.push((*n, *d, *c));
            }
            out
        };
        for (num_b, den_b, _) in &iter_blocks() {
            let num = total_num - num_b.value();
            let den = total_den - den_b.value();
            if den != 0.0 {
                vj.push(num / den);
            }
        }
        if vj.len() < 2 {
            return None;
        }
        let n = vj.len() as f64;
        let mean = vj.iter().sum::<f64>() / n;
        let sv = ((n - 1.0) / n) * vj.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>();
        Some((total_num / total_den, sv.sqrt(), vj))
    }

    /// Number of blocks.
    #[must_use]
    pub fn n_blocks(&self) -> usize {
        self.blocks.len() + usize::from(self.open.is_some())
    }

    /// The delete-one-block jackknife for an additive statistic
    /// (`den ≡ 1`): `(Σnum, se, n_blocks)` (doc_v7 D7-9).
    #[must_use]
    pub fn jackknife_sum(&self) -> Option<(f64, f64, usize)> {
        let n_blocks = self.n_blocks();
        if n_blocks < 2 {
            return None;
        }
        let total = self.total.0.value();
        let mut vj = Vec::with_capacity(n_blocks);
        let mut each = |num_b: f64| vj.push(total - num_b);
        for (num_b, _, _) in self.blocks.values() {
            each(num_b.value());
        }
        if let Some((_, num_b, _, _)) = &self.open {
            each(num_b.value());
        }
        let n = vj.len() as f64;
        let mean = vj.iter().sum::<f64>() / n;
        let sv = ((n - 1.0) / n) * vj.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>();
        Some((total, sv.sqrt(), vj.len()))
    }
}

#[derive(Clone, Debug)]
struct BpObservation {
    contig: std::sync::Arc<[u8]>,
    pos0: u64,
    num: f64,
    den: f64,
    usable: bool,
}

#[derive(Clone, Copy)]
struct BpSite<'a> {
    contig: &'a [u8],
    pos0: u64,
    num: f64,
    den: f64,
    usable: bool,
}

#[derive(Clone, Debug)]
struct BpBlock {
    num: NeumaierSum,
    den: NeumaierSum,
    n_sites: u64,
}

/// ADMIXTOOLS-style observed-site block jackknife. A new block starts at a
/// contig change or when the position is at least `block_bp` from the current
/// block's first position, not its previous position. Physical positions use
/// 1 Morgan = 100 Mb, so the CLI's 5 Mb matches ADMIXTOOLS' 0.05-Morgan
/// physical-position mode (doc_v27 V27-35).
///
/// The input positions form one grid shared by every f-statistic. Each test
/// stores its small per-site `(num, den)` contributions and counts its own
/// informative sites as weights. The ordered sink joins the grid once and
/// evaluates blocks from the joined positions, so a local chunk can neither
/// move a block boundary nor misattribute a statistic's contributions.
#[derive(Clone, Debug)]
pub struct BpBlockJackknife {
    block_bp: u64,
    observations: Vec<BpObservation>,
    /// Shared-grid statistics store only their per-site values. The caller
    /// supplies the joined position grid to `jackknife_shared` at finish.
    shared: Option<Vec<(f64, f64, bool)>>,
}

impl BpBlockJackknife {
    /// A jackknife with a maximum observed span of `block_bp` bp per block.
    #[must_use]
    pub fn new(block_bp: u64) -> Self {
        Self {
            block_bp: block_bp.max(1),
            observations: Vec::new(),
            shared: None,
        }
    }

    /// Uses the caller's one input position grid instead of retaining another
    /// copy for this statistic. Call before the first `push`.
    pub fn share_grid(&mut self) {
        assert!(
            self.observations.is_empty(),
            "invariant: grid shared before input"
        );
        self.shared = Some(Vec::new());
    }

    /// Absorbs one usable grid record. `(0, 0)` is a valid zero contribution,
    /// not a skip: it still increments the block's informative-site weight.
    pub fn push(&mut self, contig: &[u8], pos0: u64, num: f64, den: f64) {
        self.push_usable(contig, pos0, num, den, true);
    }

    /// Marks whether this statistic can use the record. A skipped record
    /// still delimits the common grid but never contributes a weight.
    pub fn push_usable(&mut self, contig: &[u8], pos0: u64, num: f64, den: f64, usable: bool) {
        if let Some(values) = &mut self.shared {
            values.push((num, den, usable));
            return;
        }
        if let Some(previous) = self.observations.last() {
            assert!(
                previous.contig.as_ref() != contig || previous.pos0 <= pos0,
                "invariant: jackknife grid arrives in contig/position order"
            );
        }
        self.observations.push(BpObservation {
            contig: std::sync::Arc::from(contig),
            pos0,
            num,
            den,
            usable,
        });
    }

    /// The natural point estimate Σnum / Σden.
    #[must_use]
    pub fn estimate(&self) -> Option<f64> {
        let (num, den) = self.totals();
        (self.n_sites() > 0 && den != 0.0).then(|| num / den)
    }

    /// Number of records in the grid, including records a test cannot use.
    #[must_use]
    pub fn n_sites(&self) -> u64 {
        self.shared
            .as_ref()
            .map_or(self.observations.len(), Vec::len) as u64
    }

    /// Number of blocks in the caller's joined grid for a shared statistic.
    /// This is independent of whether a particular test can use every site.
    #[must_use]
    pub fn n_blocks_shared(&self, grid: &Self) -> usize {
        assert!(grid.shared.is_none(), "invariant: owned grid");
        assert_eq!(self.block_bp, grid.block_bp, "invariant: same block size");
        assert_eq!(self.n_sites(), grid.n_sites(), "invariant: complete grid");
        grid.n_blocks()
    }

    /// Number of blocks for an owned grid. A shared-grid statistic has no
    /// block grid of its own, so this returns 0 until the caller evaluates it
    /// with `jackknife_shared` or `jackknife_sum_shared`.
    #[must_use]
    pub fn n_blocks(&self) -> usize {
        if self.shared.is_none() {
            self.blocks(None).len()
        } else {
            0
        }
    }

    fn totals(&self) -> (f64, f64) {
        let mut num = NeumaierSum::new();
        let mut den = NeumaierSum::new();
        if let Some(values) = &self.shared {
            for &(n, d, _) in values {
                num.add(n);
                den.add(d);
            }
        } else {
            for observation in &self.observations {
                num.add(observation.num);
                den.add(observation.den);
            }
        }
        (num.value(), den.value())
    }

    /// Materializes the one global block grid and this statistic's totals.
    fn blocks(&self, grid: Option<&Self>) -> Vec<BpBlock> {
        let sites: Vec<BpSite<'_>> = match &self.shared {
            Some(values) => {
                let grid = grid.expect("invariant: shared statistic receives its joined grid");
                assert!(grid.shared.is_none(), "invariant: owned position grid");
                assert_eq!(
                    values.len(),
                    grid.observations.len(),
                    "invariant: complete shared grid"
                );
                grid.observations
                    .iter()
                    .zip(values)
                    .map(|(site, &(num, den, usable))| BpSite {
                        contig: &site.contig,
                        pos0: site.pos0,
                        num,
                        den,
                        usable,
                    })
                    .collect()
            }
            None => self
                .observations
                .iter()
                .map(|site| BpSite {
                    contig: &site.contig,
                    pos0: site.pos0,
                    num: site.num,
                    den: site.den,
                    usable: site.usable,
                })
                .collect(),
        };
        let mut out: Vec<BpBlock> = Vec::new();
        let mut origin: Option<(&[u8], u64)> = None;
        for site in sites {
            let new_block = origin.is_none_or(|(start_contig, start_pos)| {
                start_contig != site.contig || site.pos0 - start_pos >= self.block_bp
            });
            if new_block {
                out.push(BpBlock {
                    num: NeumaierSum::new(),
                    den: NeumaierSum::new(),
                    n_sites: 0,
                });
                origin = Some((site.contig, site.pos0));
            }
            let block = out.last_mut().expect("invariant: opened above");
            block.num.add(site.num);
            block.den.add(site.den);
            block.n_sites += u64::from(site.usable);
        }
        out
    }

    /// Busing/ADMIXTOOLS weighted jackknife (Busing et al. 1999). The
    /// natural estimate is retained; only the SE uses weighted pseudovalues.
    /// Ratio-form delete-one denominators have ADMIXTOOLS' `1e-10` floor;
    /// a zero full-data denominator still has no valid point estimate.
    fn weighted(mean: f64, deleted: Vec<(f64, u64)>) -> Option<(f64, f64, usize)> {
        let usable: Vec<(f64, f64)> = deleted
            .into_iter()
            .filter(|(value, weight)| *weight > 0 && value.is_finite())
            .map(|(value, weight)| (value, weight as f64))
            .collect();
        if usable.len() < 2 || !mean.is_finite() {
            return None;
        }
        let total_weight: f64 = usable.iter().map(|(_, weight)| weight).sum();
        let pseudo_mean = usable
            .iter()
            .map(|(deleted, weight)| mean - deleted + weight * deleted / total_weight)
            .sum::<f64>();
        let variance = usable
            .iter()
            .map(|(deleted, weight)| {
                let h = total_weight / weight;
                let pseudo = h * mean - (h - 1.0) * deleted - pseudo_mean;
                pseudo * pseudo / (h - 1.0)
            })
            .sum::<f64>()
            / usable.len() as f64;
        Some((mean, variance.sqrt(), usable.len()))
    }

    /// Ratio-of-sums weighted delete-one-block jackknife for an owned grid.
    ///
    /// # Panics
    /// Panics for a shared-grid statistic until `jackknife_shared` supplies
    /// the caller's joined positions.
    #[must_use]
    pub fn jackknife(&self) -> Option<(f64, f64, usize)> {
        self.jackknife_impl(None)
    }

    /// Ratio-of-sums jackknife over the caller's joined shared grid.
    #[must_use]
    pub fn jackknife_shared(&self, grid: &Self) -> Option<(f64, f64, usize)> {
        self.jackknife_impl(Some(grid))
    }

    fn jackknife_impl(&self, grid: Option<&Self>) -> Option<(f64, f64, usize)> {
        assert!(
            self.shared.is_none() || grid.is_some(),
            "shared statistic requires its joined grid"
        );
        let (total_num, total_den) = self.totals();
        if total_den == 0.0 {
            return None;
        }
        let deleted = self
            .blocks(grid)
            .iter()
            .map(|block| {
                (
                    (total_num - block.num.value()) / (total_den - block.den.value()).max(1e-10),
                    block.n_sites,
                )
            })
            .collect();
        Self::weighted(total_num / total_den, deleted).or_else(|| {
            let n_blocks = self
                .blocks(grid)
                .into_iter()
                .filter(|block| block.n_sites > 0)
                .count();
            Some((total_num / total_den, f64::NAN, n_blocks))
        })
    }

    /// Sum-scale weighted jackknife for reed's established f2/f3/f4 = Σnum
    /// contract. Estimate and SE remain on the same scale.
    ///
    /// # Panics
    /// Panics for a shared-grid statistic until `jackknife_sum_shared`
    /// receives the caller's joined grid.
    #[must_use]
    pub fn jackknife_sum(&self) -> Option<(f64, f64, usize)> {
        self.jackknife_sum_impl(None)
    }

    /// Sum-scale jackknife over the caller's joined shared grid.
    #[must_use]
    pub fn jackknife_sum_shared(&self, grid: &Self) -> Option<(f64, f64, usize)> {
        self.jackknife_sum_impl(Some(grid))
    }

    fn jackknife_sum_impl(&self, grid: Option<&Self>) -> Option<(f64, f64, usize)> {
        assert!(
            self.shared.is_none() || grid.is_some(),
            "shared statistic requires its joined grid"
        );
        let (total, _) = self.totals();
        let deleted = self
            .blocks(grid)
            .iter()
            .map(|block| (total - block.num.value(), block.n_sites))
            .collect();
        Self::weighted(total, deleted).or_else(|| {
            let n_blocks = self
                .blocks(grid)
                .into_iter()
                .filter(|block| block.n_sites > 0)
                .count();
            Some((total, f64::NAN, n_blocks))
        })
    }

    /// Joins a later owned grid. The actual block partition is computed at
    /// finish, after all positions are known; no worker-local boundaries leak.
    pub fn merge_grid(&mut self, other: &Self) {
        assert!(
            self.shared.is_none() && other.shared.is_none(),
            "invariant: owned grid"
        );
        assert_eq!(self.block_bp, other.block_bp, "invariant: same block size");
        if let (Some(left), Some(right)) = (self.observations.last(), other.observations.first()) {
            assert!(
                left.contig != right.contig || left.pos0 <= right.pos0,
                "invariant: ordered jackknife partials do not overlap"
            );
        }
        self.observations.extend_from_slice(&other.observations);
    }

    /// Joins per-site values once the caller has merged their shared grid.
    pub fn merge_shared(&mut self, other: &Self, grid: &Self) {
        assert!(
            self.shared.is_some() && other.shared.is_some(),
            "invariant: shared statistic"
        );
        assert!(grid.shared.is_none(), "invariant: owned grid");
        assert_eq!(self.block_bp, other.block_bp, "invariant: same block size");
        assert_eq!(
            self.block_bp, grid.block_bp,
            "invariant: same grid block size"
        );
        assert_eq!(
            self.n_sites() + other.n_sites(),
            grid.n_sites(),
            "invariant: one complete grid"
        );
        self.shared
            .as_mut()
            .expect("invariant: shared values")
            .extend_from_slice(other.shared.as_deref().expect("invariant: shared values"));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_is_global_not_chunk_local() {
        // 100 sites of deterministic pseudo-data, fed in two partials split
        // at an arbitrary point (55) — the block grid is the global ordinal
        // grid, so the result matches the single-partial run exactly.
        let values: Vec<(f64, f64)> = (0..100)
            .map(|i| {
                (
                    (f64::from(i) * 7.0 % 13.0) - 3.0,
                    5.0 + (f64::from(i) * 11.0 % 17.0),
                )
            })
            .collect();
        let mut whole = BlockJackknife::new(10);
        whole.accept(&BlockPartial {
            values: values.clone(),
        });

        let mut split = BlockJackknife::new(10);
        split.accept(&BlockPartial {
            values: values[..55].to_vec(),
        });
        split.accept(&BlockPartial {
            values: values[55..].to_vec(),
        });

        assert_eq!(split.n_blocks(), whole.n_blocks());
        let (est_w, se_w, vj_w) = whole.jackknife().expect("invariant: blocks");
        let (est_s, se_s, vj_s) = split.jackknife().expect("invariant: blocks");
        assert!((est_w - est_s).abs() < 1e-12);
        assert!((se_w - se_s).abs() < 1e-12);
        for (x, y) in vj_w.iter().zip(&vj_s) {
            assert!((x - y).abs() < 1e-12);
        }
    }

    #[test]
    fn se_matches_brute_force_delete_one() {
        let values: Vec<(f64, f64)> = (0..200)
            .map(|i| {
                (
                    (f64::from(i) * 7.0 % 13.0) - 3.0,
                    5.0 + (f64::from(i) * 11.0 % 17.0),
                )
            })
            .collect();
        let block_size = 20_u64;
        let mut jk = BlockJackknife::new(block_size);
        for (num, den) in &values {
            jk.push(*num, *den);
        }
        let (est, se, vj) = jk.jackknife().expect("invariant: ≥ 2 blocks");
        let n_blocks = values.len().div_ceil(usize::try_from(block_size).unwrap());
        let mut sums = vec![(0.0, 0.0); n_blocks];
        for (i, (num, den)) in values.iter().enumerate() {
            let (n, d) = &mut sums[i / usize::try_from(block_size).unwrap()];
            *n += num;
            *d += den;
        }
        let total_num: f64 = sums.iter().map(|(n, _)| n).sum();
        let total_den: f64 = sums.iter().map(|(_, d)| d).sum();
        assert!((est - total_num / total_den).abs() < 1e-9);
        let mut vj_brute = Vec::new();
        for (bn, bd) in &sums {
            let num = total_num - bn;
            let den = total_den - bd;
            if den != 0.0 {
                vj_brute.push(num / den);
            }
        }
        let n = vj_brute.len() as f64;
        let mean = vj_brute.iter().sum::<f64>() / n;
        let se_brute = (((n - 1.0) / n)
            * vj_brute
                .iter()
                .map(|v| (v - mean) * (v - mean))
                .sum::<f64>())
        .sqrt();
        assert_eq!(vj.len(), vj_brute.len());
        assert!((se - se_brute).abs() < 1e-9, "{se} vs {se_brute}");
    }

    #[test]
    fn needs_two_blocks() {
        let mut jk = BlockJackknife::new(10);
        jk.push(1.0, 2.0);
        jk.push(3.0, 4.0);
        assert_eq!(jk.n_blocks(), 1, "one open block");
        assert!(jk.jackknife().is_none());
        // block_size = 1: every site is its own block.
        let mut jk = BlockJackknife::new(1);
        jk.push(1.0, 2.0);
        assert!(jk.jackknife().is_none());
        jk.push(3.0, 4.0);
        assert_eq!(jk.n_blocks(), 2, "two sealed singleton blocks");
        assert!(jk.jackknife().is_some());
    }

    #[test]
    fn bp_grid_uses_the_first_site_of_the_current_block() {
        let mut jk = BpBlockJackknife::new(6);
        for pos in [0, 5, 10, 15] {
            jk.push(b"chr1", pos, 0.0, 1.0);
        }
        let boundaries: Vec<_> = jk
            .observations
            .iter()
            .map(|observation| (observation.contig.as_ref(), observation.pos0))
            .collect();
        assert_eq!(
            boundaries,
            [
                (b"chr1".as_slice(), 0),
                (b"chr1", 5),
                (b"chr1", 10),
                (b"chr1", 15)
            ]
        );
        assert_eq!(jk.n_blocks(), 2);
    }

    #[test]
    fn bp_grid_is_independent_of_the_ordered_chunk_partition() {
        let sites = [
            (b"chr1".as_slice(), 0),
            (b"chr1", 5),
            (b"chr1", 10),
            (b"chr1", 15),
            (b"chr2", 0),
            (b"chr2", 5),
        ];
        let mut whole = BpBlockJackknife::new(6);
        for &(contig, pos) in &sites {
            whole.push(contig, pos, pos as f64 + 1.0, 2.0);
        }
        for split in 1..sites.len() {
            let mut a = BpBlockJackknife::new(6);
            let mut b = BpBlockJackknife::new(6);
            for &(contig, pos) in &sites[..split] {
                a.push(contig, pos, pos as f64 + 1.0, 2.0);
            }
            for &(contig, pos) in &sites[split..] {
                b.push(contig, pos, pos as f64 + 1.0, 2.0);
            }
            a.merge_grid(&b);
            assert_eq!(a.n_blocks(), whole.n_blocks(), "split {split}");
            assert_eq!(a.jackknife(), whole.jackknife(), "split {split}");
            fn observations(jk: &BpBlockJackknife) -> Vec<(&[u8], u64, f64, f64)> {
                jk.observations
                    .iter()
                    .map(|site| (site.contig.as_ref(), site.pos0, site.num, site.den))
                    .collect()
            }
            assert_eq!(observations(&a), observations(&whole), "split {split}");
        }
    }

    #[test]
    fn unequal_blocks_use_informative_site_weights() {
        let mut jk = BpBlockJackknife::new(10);
        // Eight sites in block 0, two in block 1, with per-site f2 values
        // 0 and 1 respectively. Deleted sums are 2 and 0, h=1.25 and 5;
        // pseudo values are 2 and 2, their ADMIXTOOLS center is 3.6, and
        // the variance is ((2-3.6)^2/0.25 + (2-3.6)^2/4) / 2 = 3.2^2.
        for pos in 0..8 {
            jk.push(b"chr1", pos, 0.0, 1.0);
        }
        for pos in 10..12 {
            jk.push(b"chr1", pos, 1.0, 1.0);
        }
        let (estimate, se, n_blocks) = jk.jackknife_sum().expect("invariant: two blocks");
        assert_eq!(estimate, 2.0);
        assert_eq!(n_blocks, 2);
        assert!((se - 3.2).abs() < 1e-12, "{se}");
    }

    #[test]
    fn one_informative_block_has_an_undefined_se() {
        let mut jk = BpBlockJackknife::new(6);
        jk.push_usable(b"chr1", 0, 0.0, 0.0, false);
        jk.push_usable(b"chr1", 1, 0.0, 0.0, false);
        jk.push(b"chr1", 10, 2.0, 1.0);
        let (_, se, n_blocks) = jk.jackknife().expect("invariant: valid natural estimate");
        assert!(!se.is_finite());
        assert_eq!(n_blocks, 1);
        assert_eq!(jk.n_blocks(), 2, "skips still delimit the shared grid");
    }

    #[test]
    fn shared_grid_matches_the_whole_at_every_split() {
        let sites = [
            (b"chr1".as_slice(), 0, 0.0, 0.0),
            (b"chr1", 5, 1.0, 2.0),
            (b"chr1", 10, 0.0, 0.0),
            (b"chr1", 15, 3.0, 4.0),
            (b"chr2", 0, 2.0, 1.0),
            (b"chr2", 5, 4.0, 2.0),
        ];
        let mut whole = BpBlockJackknife::new(6);
        for &(contig, pos, num, den) in &sites {
            whole.push(contig, pos, num, den);
        }
        let mut whole_grid = BpBlockJackknife::new(6);
        whole_grid.share_grid();
        for &(contig, pos, num, den) in &sites {
            whole_grid.push(contig, pos, num, den);
        }
        let expected = whole_grid.jackknife_shared(&whole);
        let blocks = whole.blocks(None);
        assert_eq!(
            blocks[0].n_sites, 2,
            "a zero contribution is still informative"
        );
        assert_eq!(blocks[1].n_sites, 2, "each block retains its valid count");
        assert_eq!(whole_grid.n_sites(), 6, "every record is retained");
        assert_eq!(
            whole_grid.n_blocks_shared(&whole),
            3,
            "the shared grid owns block boundaries"
        );
        for split in 1..sites.len() {
            let mut grid_a = BpBlockJackknife::new(6);
            let mut grid_b = BpBlockJackknife::new(6);
            let mut stat_a = BpBlockJackknife::new(6);
            let mut stat_b = BpBlockJackknife::new(6);
            stat_a.share_grid();
            stat_b.share_grid();
            for &(contig, pos, num, den) in &sites[..split] {
                grid_a.push(contig, pos, 0.0, 0.0);
                stat_a.push(contig, pos, num, den);
            }
            for &(contig, pos, num, den) in &sites[split..] {
                grid_b.push(contig, pos, 0.0, 0.0);
                stat_b.push(contig, pos, num, den);
            }
            grid_a.merge_grid(&grid_b);
            stat_a.merge_shared(&stat_b, &grid_a);
            assert_eq!(stat_a.jackknife_shared(&grid_a), expected, "split {split}");
        }
    }

    #[test]
    fn a_sparse_block_span_is_not_a_fixed_coordinate_bin() {
        let mut jk = BpBlockJackknife::new(5_000_000);
        // 1-based 100, 4.9M and 5,000,001: the final span is 4,999,901 bp,
        // so ADMIXTOOLS' current-block-first rule keeps one block.
        for pos in [99, 4_899_999, 5_000_000] {
            jk.push(b"chr1", pos, 1.0, 1.0);
        }
        assert_eq!(jk.n_blocks(), 1, "{:?}", jk.observations);
        let (_, se, n_blocks) = jk.jackknife().expect("invariant: natural ratio");
        assert!(!se.is_finite());
        assert_eq!(n_blocks, 1);
        // Exactly block_bp from the current block's first site opens a new
        // block (the reference uses `>=`, not `>`).
        jk.push(b"chr1", 5_000_099, 1.0, 1.0);
        assert_eq!(jk.n_blocks(), 2);
        assert!(jk.jackknife().is_some());
    }

    #[test]
    #[should_panic(expected = "shared statistic requires its joined grid")]
    fn a_shared_statistic_cannot_evaluate_without_its_grid() {
        let mut stat = BpBlockJackknife::new(6);
        stat.share_grid();
        stat.push(b"chr1", 0, 0.0, 0.0);
        let _ = stat.jackknife();
    }

    /// doc_v23 V23-28: merging jackknives with different block sizes panics
    /// in release, not just debug.
    #[test]
    #[should_panic(expected = "same block size")]
    fn bp_merge_rejects_a_different_block_size() {
        let mut a = BpBlockJackknife::new(100);
        let b = BpBlockJackknife::new(200);
        a.merge_grid(&b);
    }
}
