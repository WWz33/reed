//! Binned window accumulation (`docs/spec/architecture.md` §8.5): a `Bp`
//! grid is cut into bins of width `w = gcd(size, step)`, keyed by the
//! *absolute* bin index `pos0 / w` so chunks accumulate independently
//! without anchor negotiation. Each site is pushed into exactly one bin's
//! accumulators; a window is the merge of its bins.

use std::collections::BTreeMap;

use reed_core::count::SiteGroupCounts;
use reed_core::grid::{BpWindow, WindowGrid};

use crate::request::StatAcc;

/// Bin layout for one contig: absolute bins of a fixed bp width.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BpBinning {
    width: u64,
}

impl BpBinning {
    /// A bin layout with an explicit width (must be the grid's
    /// `gcd(size, step)`; prefer [`BpBinning::from_grid`]).
    #[must_use]
    pub const fn new(width: u64) -> Self {
        Self { width }
    }

    /// The bin layout of a `Bp` grid; `None` for non-Bp grids.
    #[must_use]
    pub fn from_grid(grid: &WindowGrid) -> Option<Self> {
        grid.bin_width_bp().map(|width| Self { width })
    }

    /// Bin width in bp (`gcd(size, step)` of the grid).
    #[must_use]
    pub const fn width(&self) -> u64 {
        self.width
    }

    /// Absolute bin index of a 0-based position.
    #[must_use]
    pub const fn bin_of(&self, pos0: u64) -> u64 {
        pos0 / self.width
    }
}

/// Sparse per-contig bin accumulators, mergeable across chunks in any
/// partition of the site sequence (the merge is component-wise addition).
///
/// Untouched bins inside a finished window behave as zero accumulators, so
/// empty windows evaluate to `None` values with zero counts (§8.5).
#[derive(Clone, Debug, Default)]
pub struct ContigBins {
    width: u64,
    seed: Vec<StatAcc>,
    bins: BTreeMap<u64, Vec<StatAcc>>,
    n_sites: u64,
}

impl ContigBins {
    /// Empty bin set for `binning`, seeding new bins with fresh copies of
    /// `seed` (one accumulator per request, in request order).
    #[must_use]
    pub fn new(binning: BpBinning, seed: &[StatAcc]) -> Self {
        Self {
            width: binning.width(),
            seed: seed.to_vec(),
            bins: BTreeMap::new(),
            n_sites: 0,
        }
    }

    /// Absorbs one site into its bin (created on first touch).
    pub fn push_site(&mut self, pos0: u64, row: &[SiteGroupCounts], all_variant: bool) {
        let bin = pos0 / self.width;
        let seed = &self.seed;
        let accs = self.bins.entry(bin).or_insert_with(|| seed.clone());
        for acc in accs {
            acc.push(row, all_variant);
        }
        self.n_sites += 1;
    }

    /// Registers a present-but-uncountable site in its bin (doc_v7 D7-11).
    pub fn push_skipped(&mut self, pos0: u64) {
        let bin = pos0 / self.width;
        let seed = &self.seed;
        let accs = self.bins.entry(bin).or_insert_with(|| seed.clone());
        for acc in accs {
            acc.note_skipped();
        }
    }

    /// Merges another partial's bins into this one, consuming `other`.
    pub fn merge_bins(&mut self, other: &mut Self) {
        debug_assert_eq!(self.width, other.width, "invariant: same bin width");
        for (bin, accs) in std::mem::take(&mut other.bins) {
            match self.bins.get_mut(&bin) {
                Some(dst) => {
                    for (dst, src) in dst.iter_mut().zip(&accs) {
                        dst.merge(src);
                    }
                }
                None => {
                    self.bins.insert(bin, accs);
                }
            }
        }
        self.n_sites += other.n_sites;
    }

    /// Sites absorbed so far.
    #[must_use]
    pub const fn n_sites(&self) -> u64 {
        self.n_sites
    }

    /// Whether any touched bin falls in `range` (a window with no touched
    /// bin is an empty window).
    #[must_use]
    pub fn any_bin_in(&self, range: std::ops::Range<u64>) -> bool {
        self.bins.range(range).next().is_some()
    }

    /// Merges every touched bin's accumulators into `totals`, in bin order
    /// (the contig/genome aggregation of §8.7).
    pub fn fold(&self, totals: &mut [StatAcc]) {
        for accs in self.bins.values() {
            for (dst, src) in totals.iter_mut().zip(accs) {
                dst.merge(src);
            }
        }
    }

    /// Frees every bin strictly below `bin` — windows that start at or
    /// after `bin` can never read them again (the watermark of §8.3).
    pub fn drop_bins_below(&mut self, bin: u64) {
        let tail = self.bins.split_off(&bin);
        self.bins = tail;
    }

    /// Finishes one window: clones and merges its bins' accumulators in
    /// bin order (missing bins contribute zeros). The result lines up with
    /// the request order the seed was built in.
    #[must_use]
    pub fn finish_window(&self, window: &BpWindow) -> Vec<StatAcc> {
        self.finish_bins(window.bins.start.0..window.bins.end.0)
    }

    /// Finishes an absolute bin range (the driver's watermark path builds
    /// window ranges directly, without [`BpWindow`] values).
    #[must_use]
    pub fn finish_bins(&self, bins: std::ops::Range<u64>) -> Vec<StatAcc> {
        let mut out = self.seed.clone();
        for (_, accs) in self.bins.range(bins) {
            for (dst, src) in out.iter_mut().zip(accs) {
                dst.merge(src);
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fst::FstKind;
    use crate::request::{StatKind, build_accs};

    fn seed() -> Vec<StatAcc> {
        build_accs(
            &[
                StatKind::Pi {
                    group: 0,
                    per_bp: false,
                },
                StatKind::Fst {
                    a: 0,
                    b: 1,
                    kind: FstKind::Hudson,
                },
            ],
            &[2, 2],
        )
        .unwrap()
    }

    fn row(refs: [u32; 2], alts: [u32; 2]) -> Vec<reed_core::count::SiteGroupCounts> {
        [
            reed_core::count::SiteGroupCounts {
                n_ref_alleles: refs[0],
                n_alt_alleles: alts[0],
                ..Default::default()
            },
            reed_core::count::SiteGroupCounts {
                n_ref_alleles: refs[1],
                n_alt_alleles: alts[1],
                ..Default::default()
            },
        ]
        .into()
    }

    #[test]
    fn bins_partition_sites_and_windows_merge_their_range() {
        let grid = WindowGrid::Bp {
            size: 100,
            step: 100,
        };
        let binning = BpBinning::from_grid(&grid).unwrap();
        let mut bins = ContigBins::new(binning, &seed());
        // Two sites in bin 0 (pos 0..50), one in bin 1 (pos 150).
        bins.push_site(10, &row([3, 1], [1, 3]), true);
        bins.push_site(40, &row([4, 0], [0, 0]), false);
        bins.push_site(150, &row([2, 2], [2, 2]), true);
        assert_eq!(bins.n_sites(), 3);

        let windows = grid.bp_windows(10, 150, 1000).unwrap();
        assert_eq!(windows.len(), 2);
        let w0 = bins.finish_window(&windows[0]);
        let w1 = bins.finish_window(&windows[1]);
        // Window 0 saw both groups' counts at two sites; window 1 at one.
        let StatAcc::Pi(pi0) = &w0[0] else {
            panic!("pi first")
        };
        let StatAcc::Pi(pi1) = &w1[0] else {
            panic!("pi first")
        };
        assert_eq!(pi0.counts().comparable, 6 + 6);
        assert_eq!(pi1.counts().comparable, 6);
        // Window 2 (beyond the last site) is all zeros → value None.
        let empty_window = grid.bp_windows(900, 950, 1000).unwrap();
        let w2 = bins.finish_window(&empty_window[0]);
        let StatAcc::Pi(pi2) = &w2[0] else {
            panic!("pi first")
        };
        assert_eq!(pi2.counts().comparable, 0);
        assert!(pi2.value_per_pair().is_none());
    }

    #[test]
    fn merge_partitions_equal_whole() {
        let binning = BpBinning { width: 50 };
        let mut whole = ContigBins::new(binning, &seed());
        for (pos, refs, alts) in [
            (10_u64, [3_u32, 1], [1_u32, 3]),
            (60, [4, 0], [0, 0]),
            (70, [2, 2], [2, 2]),
        ] {
            whole.push_site(pos, &row(refs, alts), true);
        }
        // Same sites split into two partials at an arbitrary point.
        let mut part_a = ContigBins::new(binning, &seed());
        let mut part_b = ContigBins::new(binning, &seed());
        part_a.push_site(10, &row([3, 1], [1, 3]), true);
        part_b.push_site(60, &row([4, 0], [0, 0]), true);
        part_b.push_site(70, &row([2, 2], [2, 2]), true);
        part_a.merge_bins(&mut part_b);
        assert_eq!(part_a.n_sites(), whole.n_sites());

        let grid = WindowGrid::Bp {
            size: 100,
            step: 100,
        };
        let windows = grid.bp_windows(10, 70, 1000).unwrap();
        let merged = part_a.finish_window(&windows[0]);
        let direct = whole.finish_window(&windows[0]);
        let (StatAcc::Pi(a), StatAcc::Pi(b)) = (&merged[0], &direct[0]) else {
            panic!("pi first");
        };
        assert_eq!(a.counts(), b.counts());
        assert_eq!(a.value_per_pair(), b.value_per_pair());
    }
}

/// Per-region accumulators for a `Regions` grid (`--bed-file`): one
/// accumulator set per interval, in input order; a site stabs every
/// covering region (overlaps each accumulate it — BED-window semantics).
/// Merge is component-wise per region index; workers on the same contig
/// share one `RegionSet` (Arc), so the indices line up by construction.
///
/// `totals` receives each site at most once (a site inside several
/// overlapping regions still folds once) — the contig-level aggregation.
#[derive(Clone, Debug)]
pub struct RegionAccs {
    regions: std::sync::Arc<reed_core::grid::RegionSet>,
    accs: Vec<Vec<StatAcc>>,
    totals: Vec<StatAcc>,
    touched: Vec<bool>,
    n_sites: u64,
}

impl RegionAccs {
    /// Fresh accumulators for every region of `regions`.
    #[must_use]
    pub fn new(regions: std::sync::Arc<reed_core::grid::RegionSet>, seed: &[StatAcc]) -> Self {
        let n = regions.windows().len();
        Self {
            regions,
            accs: vec![seed.to_vec(); n],
            totals: seed.to_vec(),
            touched: vec![false; n],
            n_sites: 0,
        }
    }

    /// Absorbs one site into every covering region (and once into the
    /// contig totals when it landed in at least one).
    pub fn push_site(&mut self, pos0: u64, row: &[SiteGroupCounts], all_variant: bool) {
        let mut any = false;
        for idx in self.regions.stab(pos0) {
            for acc in &mut self.accs[idx] {
                acc.push(row, all_variant);
            }
            self.touched[idx] = true;
            any = true;
        }
        if any {
            for acc in &mut self.totals {
                acc.push(row, all_variant);
            }
        }
        self.n_sites += 1;
    }

    /// Registers a present-but-uncountable site in every covering region
    /// (doc_v7 D7-11).
    pub fn push_skipped(&mut self, pos0: u64) {
        let mut any = false;
        for idx in self.regions.stab(pos0) {
            for acc in &mut self.accs[idx] {
                acc.note_skipped();
            }
            self.touched[idx] = true;
            any = true;
        }
        if any {
            for acc in &mut self.totals {
                acc.note_skipped();
            }
        }
        self.n_sites += 1;
    }

    /// Merges another partial's regions (same `RegionSet` by construction).
    pub fn merge_regions(&mut self, other: &Self) {
        debug_assert!(
            std::sync::Arc::ptr_eq(&self.regions, &other.regions),
            "invariant: one RegionSet per contig"
        );
        for ((dst, src), (dt, st)) in self
            .accs
            .iter_mut()
            .zip(&other.accs)
            .zip(self.touched.iter_mut().zip(&other.touched))
        {
            for (d, s) in dst.iter_mut().zip(src) {
                d.merge(s);
            }
            *dt |= *st;
        }
        for (d, s) in self.totals.iter_mut().zip(&other.totals) {
            d.merge(s);
        }
        self.n_sites += other.n_sites;
    }

    /// Whether any site landed in region `idx`.
    #[must_use]
    pub fn touched(&self, idx: usize) -> bool {
        self.touched[idx]
    }

    /// Region `idx`'s accumulators (caller decides empty-window policy).
    #[must_use]
    pub fn accs(&self, idx: usize) -> &[StatAcc] {
        &self.accs[idx]
    }

    /// Merges the contig totals into `totals` (each site counted once).
    pub fn fold(&self, totals: &mut [StatAcc]) {
        for (dst, src) in totals.iter_mut().zip(&self.totals) {
            dst.merge(src);
        }
    }

    /// The region set.
    #[must_use]
    pub fn regions(&self) -> &reed_core::grid::RegionSet {
        &self.regions
    }

    /// Sites absorbed.
    #[must_use]
    pub const fn n_sites(&self) -> u64 {
        self.n_sites
    }
}
