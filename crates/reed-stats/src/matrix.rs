//! The site-major bit-plane matrix and its pairwise LD scan
//! (docs/rounds/v6/02-models.md §2.2).
//!
//! A [`ContigMatrix`] holds one contig's analysed sites, in position order,
//! as flat site-major plane blocks — one block per contributing pass0 run
//! (`[slot0_alt, slot0_miss, slot1_alt, slot1_miss]`, `w` words each). The
//! matrix is the read-only input of the global analyses: `--indep-pairwise`
//! (contig-parallel, `reed`), `--ld-decay` (positional-block parallel) and
//! `--mind` (column popcounts). This crate stays pure; the parallelism
//! lives with the caller.
//!
//! Samples are laid out **group-major**: each group is padded to a `u64`
//! word boundary, so a group's pair sums read one contiguous word slice and
//! the cohort (all stored samples) sums are the exactly-additive sum of the
//! group sums. Per-population LD decay therefore costs about the same as
//! the whole-cohort scan instead of a multiple of it.

use std::ops::Range;
use std::sync::Arc;

use crate::ld::{LdDecayAcc, PairSums, pair_sums};

/// The plane bit layout of the samples: one slice per group (word-aligned)
/// plus a trailing slice for kept-but-ungrouped samples (they still count
/// in the cohort).
#[derive(Clone, Debug)]
pub struct SampleLayout {
    /// Plane bit of each original sample (`u32::MAX` = not stored).
    dest: Vec<u32>,
    slices: Vec<Slice>,
    /// How many leading `slices` are named groups.
    n_named: usize,
    w: usize,
}

/// One contiguous, word-aligned sample slice in the plane layout.
#[derive(Clone, Copy, Debug)]
struct Slice {
    /// First word of the slice in each plane.
    word: usize,
    /// Samples in the slice (the slice's tail word masks the rest).
    n: usize,
}

impl SampleLayout {
    /// Builds the layout from the per-sample group assignment: group `g`
    /// owns exactly the samples with `Some(g)`; `None` samples land in a
    /// trailing slice. Groups are packed in order, each padded to a word.
    ///
    /// `keep` (when given) drops the samples `--keep`/`--remove` selected
    /// out before anything else: a dropped sample gets no plane bit, so it
    /// never enters the cohort or a group sum (doc_v7 D7-5).
    #[must_use]
    pub fn from_assignments(
        assignments: &[Option<u32>],
        n_groups: usize,
        keep: Option<&[bool]>,
    ) -> Self {
        let kept = |sample: usize| keep.is_none_or(|k| k.get(sample).copied().unwrap_or(true));
        let mut dest = vec![u32::MAX; assignments.len()];
        let mut slices = Vec::with_capacity(n_groups + 1);
        let mut cursor = 0_usize;
        for group in 0..n_groups {
            let start = cursor;
            let want = u32::try_from(group).expect("invariant: group count fits u32");
            for (sample, slot) in assignments.iter().enumerate() {
                if *slot == Some(want) && kept(sample) {
                    dest[sample] = u32::try_from(cursor).expect("invariant: sample count fits u32");
                    cursor += 1;
                }
            }
            slices.push(Slice {
                word: start / 64,
                n: cursor - start,
            });
            cursor = cursor.div_ceil(64) * 64;
        }
        let start = cursor;
        for (sample, slot) in assignments.iter().enumerate() {
            if slot.is_none() && kept(sample) {
                dest[sample] = u32::try_from(cursor).expect("invariant: sample count fits u32");
                cursor += 1;
            }
        }
        if cursor > start {
            slices.push(Slice {
                word: start / 64,
                n: cursor - start,
            });
        }
        Self {
            dest,
            slices,
            n_named: n_groups,
            w: cursor.div_ceil(64).max(1),
        }
    }

    /// Words per plane in the laid-out matrix.
    #[must_use]
    pub const fn w(&self) -> usize {
        self.w
    }

    /// Named groups (the trailing ungrouped slice is not counted).
    #[must_use]
    pub const fn n_groups(&self) -> usize {
        self.n_named
    }

    /// The plane bit of each original sample (`u32::MAX` = not stored).
    #[must_use]
    pub fn dest(&self) -> &[u32] {
        &self.dest
    }

    /// The plane bit of one original sample (`None` = not stored).
    #[must_use]
    pub fn bit_of(&self, sample: usize) -> Option<usize> {
        let d = *self.dest.get(sample)?;
        (d != u32::MAX).then(|| usize::try_from(d).ok()).flatten()
    }
}

/// A contig's analysed sites in position order, as flat site-major plane
/// blocks.
#[derive(Clone, Debug)]
pub struct ContigMatrix {
    contig: Vec<u8>,
    layout: Arc<SampleLayout>,
    blocks: Vec<Vec<u64>>,
    /// First site index of each block; `base.len() == blocks.len() + 1`.
    base: Vec<usize>,
    pos0: Vec<u64>,
}

impl ContigMatrix {
    /// An empty matrix for one contig over a shared sample layout.
    #[must_use]
    pub fn new(contig: Vec<u8>, layout: Arc<SampleLayout>) -> Self {
        Self {
            contig,
            layout,
            blocks: Vec::new(),
            base: vec![0],
            pos0: Vec::new(),
        }
    }

    /// Appends one pass0 run: `flat` is `positions.len() · 4 · w` words.
    pub fn push_run(&mut self, flat: Vec<u64>, positions: Vec<u64>) {
        if positions.is_empty() {
            return;
        }
        debug_assert_eq!(flat.len(), positions.len() * 4 * self.layout.w());
        let n = self.pos0.len() + positions.len();
        self.blocks.push(flat);
        self.pos0.extend(positions);
        self.base.push(n);
    }

    /// Sites in the matrix.
    #[must_use]
    pub fn n_sites(&self) -> usize {
        self.pos0.len()
    }

    /// The contig name bytes.
    #[must_use]
    pub fn contig(&self) -> &[u8] {
        &self.contig
    }

    /// Every site's 0-based position (ascending).
    #[must_use]
    pub fn pos0(&self) -> &[u64] {
        &self.pos0
    }

    /// The shared sample layout.
    #[must_use]
    pub fn layout(&self) -> &SampleLayout {
        &self.layout
    }

    fn block_of(&self, i: usize) -> usize {
        self.base.partition_point(|&b| b <= i) - 1
    }

    /// The four full-layout plane slices of site `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= n_sites()`.
    #[must_use]
    pub fn planes(&self, i: usize) -> [&[u64]; 4] {
        let b = self.block_of(i);
        let w = self.layout.w();
        let k = (i - self.base[b]) * 4 * w;
        let blk = &self.blocks[b];
        [
            &blk[k..k + w],
            &blk[k + w..k + 2 * w],
            &blk[k + 2 * w..k + 3 * w],
            &blk[k + 3 * w..k + 4 * w],
        ]
    }
}

/// One per-pair row: 1-based positions and the cohort r².
pub type LdPairRow = (u64, u64, f64);

/// The accumulators one positional block produced.
#[derive(Clone, Debug)]
pub struct LdDecaySet {
    /// The whole-cohort curve.
    pub cohort: LdDecayAcc,
    /// One curve per named group.
    pub groups: Vec<LdDecayAcc>,
    /// Per-pair rows (cohort r²) when requested.
    pub pairs: Vec<LdPairRow>,
}

impl LdDecaySet {
    /// Empty accumulators over `edges` for `n_groups` named groups.
    #[must_use]
    pub fn new(edges: &[u64], n_groups: usize) -> Self {
        Self {
            cohort: LdDecayAcc::new(edges.to_vec()),
            groups: (0..n_groups)
                .map(|_| LdDecayAcc::new(edges.to_vec()))
                .collect(),
            pairs: Vec::new(),
        }
    }

    /// Component-wise, order-preserving merge (`other` must come later in
    /// the site order for a deterministic result).
    pub fn merge(&mut self, other: &Self) {
        self.cohort.merge(&other.cohort);
        for (dst, src) in self.groups.iter_mut().zip(&other.groups) {
            dst.merge(src);
        }
        self.pairs.extend_from_slice(&other.pairs);
    }
}

/// The slice of one group's words out of a full-layout plane set.
fn group_planes(x: [&[u64]; 4], slice: Slice) -> [&[u64]; 4] {
    let wg = slice.n.div_ceil(64);
    [
        &x[0][slice.word..slice.word + wg],
        &x[1][slice.word..slice.word + wg],
        &x[2][slice.word..slice.word + wg],
        &x[3][slice.word..slice.word + wg],
    ]
}

/// Scans the pairs whose later site lies in `range`, each earlier site
/// within `max_bp` — exactly the pair set the serial carry produced, in
/// `(later site, earlier site)` order. The cohort curve is derived from the
/// group sums (integer-exact, so its r² values are the ones a full-sample
/// pass would compute).
///
/// # Panics
///
/// Panics only on a corrupt matrix (block boundaries out of step with
/// `pos0`), which `push_run` maintains by construction.
#[must_use]
pub fn decay_range(
    m: &ContigMatrix,
    range: Range<usize>,
    max_bp: u64,
    edges: &[u64],
    collect_pairs: bool,
) -> LdDecaySet {
    let layout = m.layout();
    let mut set = LdDecaySet::new(edges, layout.n_groups());
    let n = m.n_sites();
    let end = range.end.min(n);
    let start = range.start.min(end);
    if end.saturating_sub(start) < 2 {
        return set;
    }
    let pos = m.pos0();
    // The lookback span: the first site within `max_bp` of the block start.
    let floor = pos[start].saturating_sub(max_bp);
    let span = pos[..end].partition_point(|&p| p < floor);
    // Per-slice site planes, hoisted out of the pair loop: the slices are
    // disjoint word ranges, so the per-pair scans share nothing but the
    // plane data itself.
    let per_slice: Vec<Vec<[&[u64]; 4]>> = layout
        .slices
        .iter()
        .map(|&sl| {
            if sl.word == 0 && sl.n.div_ceil(64) == layout.w() {
                // The slice spans the whole plane: the planes as stored.
                (span..end).map(|i| m.planes(i)).collect()
            } else {
                (span..end).map(|i| group_planes(m.planes(i), sl)).collect()
            }
        })
        .collect();
    let slice_n: Vec<usize> = layout.slices.iter().map(|sl| sl.n).collect();
    let named = layout.n_groups();
    let mut lo = 0_usize;
    for i in start..end {
        let limit = pos[i].saturating_sub(max_bp);
        while span + lo < i && pos[span + lo] < limit {
            lo += 1;
        }
        for j in (span + lo)..i {
            let dist = pos[i] - pos[j];
            // One bin index serves every accumulator (they share edges).
            let b = set.cohort.bin_index(dist);
            let mut total = PairSums::default();
            for (g, sets) in per_slice.iter().enumerate() {
                let sums = pair_sums(sets[j - span], sets[i - span], slice_n[g]);
                total.add_assign(&sums);
                if g < named
                    && let Some(r2) = sums.r2()
                {
                    set.groups[g].push_bin(b, r2);
                }
            }
            if let Some(r2) = total.r2() {
                set.cohort.push_bin(b, r2);
                if collect_pairs {
                    set.pairs.push((pos[j] + 1, pos[i] + 1, r2));
                }
            }
        }
    }
    set
}
