//! Block jackknife (§8.6): delete-one-block standard errors for any
//! statistic that reduces to per-site additive (numerator, denominator)
//! sums. Blocks live on the *global site ordinal* grid — the sink-side
//! accumulator numbers sites as chunk partials arrive in order, so the
//! block partition (and therefore the SE) is independent of the chunk
//! partition and thread count (§7.6). Blocks are fixed site counts
//! (equal-sized), so the plain delete-1 jackknife applies; unequal blocks
//! (Busing-weighted) are future work (# ponytail).

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
                self.blocks.values().cloned().collect();
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
}

/// A bp-blocked jackknife (docs/rounds/v5/05-algorithms.md §3.3): blocks are
/// fixed genomic spans (`--block-size`, default 5M bp) rather than fixed
/// site counts, which is what admixtools/Patterson block on for real data.
/// Sites arrive in coordinate order; the accumulator keys blocks by
/// `(contig, pos0 / block_bp)` so the partition is independent of the chunk
/// partition and thread count.
#[derive(Clone, Debug)]
pub struct BpBlockJackknife {
    block_bp: u64,
    /// Sealed blocks in key order: `((contig, block), Σnum, Σden, n_sites)`.
    blocks: std::collections::BTreeMap<(Vec<u8>, u64), (NeumaierSum, NeumaierSum, u64)>,
    total: (NeumaierSum, NeumaierSum),
    n_sites: u64,
}

impl BpBlockJackknife {
    /// A jackknife over blocks of `block_bp` base pairs (`0` is read as 1).
    #[must_use]
    pub fn new(block_bp: u64) -> Self {
        Self {
            block_bp: block_bp.max(1),
            blocks: std::collections::BTreeMap::new(),
            total: (NeumaierSum::new(), NeumaierSum::new()),
            n_sites: 0,
        }
    }

    /// Absorbs one site at `(contig, pos0)`.
    pub fn push(&mut self, contig: &[u8], pos0: u64, num: f64, den: f64) {
        let key = (contig.to_vec(), pos0 / self.block_bp);
        self.total.0.add(num);
        self.total.1.add(den);
        self.n_sites += 1;
        let entry = self.blocks.entry(key).or_insert_with(|| {
            (NeumaierSum::new(), NeumaierSum::new(), 0)
        });
        entry.0.add(num);
        entry.1.add(den);
        entry.2 += 1;
    }

    /// The point estimate Σnum / Σden.
    #[must_use]
    pub fn estimate(&self) -> Option<f64> {
        let den = self.total.1.value();
        (self.n_sites > 0 && den != 0.0).then(|| self.total.0.value() / den)
    }

    /// The delete-one-block jackknife: `(estimate, se, n_blocks)`.
    #[must_use]
    pub fn jackknife(&self) -> Option<(f64, f64, usize)> {
        if self.blocks.len() < 2 {
            return None;
        }
        let total_num = self.total.0.value();
        let total_den = self.total.1.value();
        if total_den == 0.0 {
            return None;
        }
        let mut vj = Vec::with_capacity(self.blocks.len());
        for (num_b, den_b, _) in self.blocks.values() {
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
        Some((total_num / total_den, sv.sqrt(), vj.len()))
    }

    /// Number of blocks seen.
    #[must_use]
    pub fn n_blocks(&self) -> usize {
        self.blocks.len()
    }

    /// Contributing site count.
    #[must_use]
    pub const fn n_sites(&self) -> u64 {
        self.n_sites
    }

    /// Merges another partial (blocks keyed by `(contig, block)` add).
    pub fn merge(&mut self, other: &Self) {
        debug_assert_eq!(self.block_bp, other.block_bp, "invariant: same block size");
        for (key, (num, den, n)) in &other.blocks {
            let entry = self.blocks.entry(key.clone()).or_insert_with(|| {
                (NeumaierSum::new(), NeumaierSum::new(), 0)
            });
            entry.0.merge(num);
            entry.1.merge(den);
            entry.2 += n;
        }
        self.total.0.merge(&other.total.0);
        self.total.1.merge(&other.total.1);
        self.n_sites += other.n_sites;
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
            .map(|i| ((f64::from(i) * 7.0 % 13.0) - 3.0, 5.0 + (f64::from(i) * 11.0 % 17.0)))
            .collect();
        let mut whole = BlockJackknife::new(10);
        whole.accept(&BlockPartial { values: values.clone() });

        let mut split = BlockJackknife::new(10);
        split.accept(&BlockPartial { values: values[..55].to_vec() });
        split.accept(&BlockPartial { values: values[55..].to_vec() });

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
            .map(|i| ((f64::from(i) * 7.0 % 13.0) - 3.0, 5.0 + (f64::from(i) * 11.0 % 17.0)))
            .collect();
        let block_size = 20_u64;
        let mut jk = BlockJackknife::new(block_size);
        for (num, den) in &values {
            jk.push(*num, *den);
        }
        let (est, se, vj) = jk.jackknife().expect("invariant: ≥ 2 blocks");
        let n_blocks = values.len().div_ceil(block_size as usize);
        let mut sums = vec![(0.0, 0.0); n_blocks];
        for (i, (num, den)) in values.iter().enumerate() {
            let (n, d) = &mut sums[i / block_size as usize];
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
            * vj_brute.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>())
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
}
