//! Window grids and binning (`docs/spec/architecture.md` §5.9, §8.5).
//!
//! `Bp` windows keep the legacy anchoring exactly: the first window starts at
//! `(first_pos / step) · step`, windows advance by `step`, the last window
//! ends at or before the contig length, and empty windows are emitted. Bp
//! windows are stored as *bins* of width `gcd(size, step)`, so each site is
//! accumulated once and a window is the merge of its bins.

use std::sync::Arc;

use crate::coord::{Interval0, Pos0};
use crate::error::CoreError;

/// A bin index: the unit of accumulation underneath windows.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BinId(pub u64);

/// How windows are laid out over a contig.
#[derive(Clone, Debug)]
pub enum WindowGrid {
    /// Fixed base-pair windows with a stride (may overlap when `step < size`).
    Bp {
        /// Window width in bp.
        size: u64,
        /// Stride in bp.
        step: u64,
    },
    /// Fixed site-count windows (each contig restarts the count).
    Sites {
        /// Sites per window.
        n_sites: u32,
        /// Stride in sites.
        step: u32,
    },
    /// Windows given explicitly as regions (BED-style, may overlap),
    /// per contig.
    Regions(Arc<ContigRegions>),
}

/// The per-contig window regions of a `Regions` grid (`--bed-file`): each
/// interval is one window, kept in input order; contigs carry the names the
/// BED file used.
#[derive(Clone, Debug, Default)]
pub struct ContigRegions {
    by_contig: std::collections::HashMap<Vec<u8>, Arc<RegionSet>>,
}

impl ContigRegions {
    /// Builds the map; empty or reversed intervals are rejected.
    pub fn new(
        entries: impl IntoIterator<Item = (Vec<u8>, Vec<Interval0>)>,
    ) -> Result<Self, CoreError> {
        let mut by_contig = std::collections::HashMap::new();
        for (contig, intervals) in entries {
            let set = RegionSet::new(intervals)?;
            by_contig.insert(contig, Arc::new(set));
        }
        Ok(Self { by_contig })
    }

    /// The windows of one contig, in input order (`None` = no regions).
    #[must_use]
    pub fn get(&self, contig: &[u8]) -> Option<Arc<RegionSet>> {
        self.by_contig.get(contig).cloned()
    }

    /// Number of contigs with regions.
    #[must_use]
    pub fn len(&self) -> usize {
        self.by_contig.len()
    }

    /// Whether there are no regions at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.by_contig.is_empty()
    }
}

impl Default for WindowGrid {
    /// The historical default grid: 1 bp windows (callers override before
    /// running; a request built field-by-field needs *some* value).
    fn default() -> Self {
        Self::Bp { size: 1, step: 1 }
    }
}

impl WindowGrid {
    /// Checks size/step positivity and region validity.
    pub fn validate(&self) -> Result<(), CoreError> {
        match self {
            Self::Bp { size, step } if *size == 0 || *step == 0 => Err(CoreError::InvalidWindow),
            Self::Sites { n_sites, step } if *n_sites == 0 || *step == 0 => {
                Err(CoreError::InvalidWindow)
            }
            Self::Regions(regions) => {
                for set in regions.by_contig.values() {
                    set.validate()?;
                }
                Ok(())
            }
            _ => Ok(()),
        }
    }

    /// Bin width in bp (Bp grids only).
    #[must_use]
    pub fn bin_width_bp(&self) -> Option<u64> {
        match *self {
            Self::Bp { size, step } => Some(gcd(size, step)),
            _ => None,
        }
    }


    /// Enumerates the Bp windows touched by sites in
    /// `[first_pos0, last_pos0]`, in order, with their bin ranges.
    ///
    /// Returns `None` for non-Bp grids.
    #[must_use]
    pub fn bp_windows(
        &self,
        first_pos0: u64,
        last_pos0: u64,
        contig_len: u64,
    ) -> Option<Vec<BpWindow>> {
        let Self::Bp { size, step } = *self else {
            return None;
        };
        let width = gcd(size, step);
        let anchor = (first_pos0 / step) * step;
        let mut windows = Vec::new();
        let mut index = 0_u64;
        let mut start = anchor;
        while start <= last_pos0 && start < contig_len {
            let end = (start + size).min(contig_len);
            windows.push(BpWindow {
                index,
                interval: Interval0::new(Pos0(start), Pos0(end)),
                // The anchor is a multiple of `step`, hence of `width`.
                bins: BinId(start / width)..BinId(end.div_ceil(width)),
            });
            index += 1;
            start += step;
        }
        Some(windows)
    }

    /// Enumerates the site-count windows touched by site ranks
    /// `[first_rank, last_rank]`. Site ranks start at 0 per contig.
    /// Returns `None` for non-Sites grids.
    #[must_use]
    pub fn sites_windows(&self, first_rank: u64, last_rank: u64) -> Option<Vec<SitesWindow>> {
        let Self::Sites { n_sites, step } = *self else {
            return None;
        };
        let (n_sites, step) = (u64::from(n_sites), u64::from(step));
        let mut windows = Vec::new();
        let mut index = 0_u64;
        let mut start = (first_rank / step) * step;
        while start <= last_rank {
            windows.push(SitesWindow {
                index,
                first_rank: start,
                last_rank_exclusive: start + n_sites,
            });
            index += 1;
            start += step;
        }
        Some(windows)
    }
}

/// One base-pair window: its 0-based half-open span and the bins whose merge
/// produces it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BpWindow {
    /// Window number within the contig (0-based, stride order).
    pub index: u64,
    /// 0-based half-open span, clipped to the contig length.
    pub interval: Interval0,
    /// Bins to merge, in order.
    pub bins: std::ops::Range<BinId>,
}

/// One site-count window over site ranks within a contig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SitesWindow {
    /// Window number within the contig.
    pub index: u64,
    /// First site rank in the window.
    pub first_rank: u64,
    /// One past the last site rank.
    pub last_rank_exclusive: u64,
}

/// A set of 0-based half-open intervals, kept in input order for window
/// identity plus a sorted merged copy — each merged span carrying the
/// input indices of its members — for probes.
#[derive(Clone, Debug, Default)]
pub struct RegionSet {
    intervals: Vec<Interval0>,
    merged: Vec<Interval0>,
    members: Vec<Vec<usize>>,
}

impl RegionSet {
    /// Builds a region set; empty or reversed intervals are rejected.
    pub fn new(intervals: Vec<Interval0>) -> Result<Self, CoreError> {
        for interval in &intervals {
            if interval.end <= interval.start {
                return Err(CoreError::InvalidInterval {
                    start0: interval.start.0,
                    end0: interval.end.0,
                });
            }
        }
        let mut sorted: Vec<(usize, Interval0)> =
            intervals.iter().copied().enumerate().collect();
        sorted.sort_by_key(|(_, interval)| (interval.start, interval.end));
        let mut merged: Vec<Interval0> = Vec::with_capacity(sorted.len());
        let mut members: Vec<Vec<usize>> = Vec::with_capacity(sorted.len());
        for (idx, interval) in sorted {
            match (merged.last_mut(), members.last_mut()) {
                (Some(last), Some(member)) if interval.start <= last.end => {
                    last.end = Pos0(last.end.0.max(interval.end.0));
                    member.push(idx);
                }
                _ => {
                    merged.push(interval);
                    members.push(vec![idx]);
                }
            }
        }
        // Members arrive in (start, end) order; stab promises input order,
        // so sort each span's list once here instead of on every probe.
        for member in &mut members {
            member.sort_unstable();
        }
        Ok(Self {
            intervals,
            merged,
            members,
        })
    }

    // Placeholder for future region invariants (kept symmetric with
    // `WindowGrid::validate`, which dispatches here); nothing to check yet,
    // so `self` is unused and `Ok` is the only outcome.
    #[allow(clippy::unused_self, clippy::unnecessary_wraps)]
    fn validate(&self) -> Result<(), CoreError> {
        Ok(())
    }

    /// Regions as given (window order), each region being one window.
    #[must_use]
    pub fn windows(&self) -> &[Interval0] {
        &self.intervals
    }

    /// Indices of the regions covering `pos0`, in input order. The merged
    /// spans are disjoint and sorted, so a binary search finds the single
    /// span that can contain `pos0`; only its members are then checked
    /// with the exact coverage predicate (a member of a span need not
    /// cover the whole span — nested regions only cover their own
    /// subrange). O(log spans + members of the span) per probe instead of
    /// a per-site linear scan over every region (doc_v20 R14). Overlapping
    /// regions each match (the BED-window convention).
    #[must_use]
    pub fn stab(&self, pos0: u64) -> Vec<usize> {
        // Number of spans starting at or before pos0; the last of them is
        // the only candidate, earlier spans end before it starts.
        let candidate = self.merged.partition_point(|span| span.start.0 <= pos0);
        let Some(span) = candidate
            .checked_sub(1)
            .filter(|&span| pos0 < self.merged[span].end.0)
        else {
            return Vec::new();
        };
        // Regions live in exactly one span and each is contained in it, so
        // regions of every other span are already excluded; members arrive
        // input-sorted, hence the filtered answer is too.
        self.members[span]
            .iter()
            .filter(|&&idx| {
                let interval = &self.intervals[idx];
                interval.start.0 <= pos0 && pos0 < interval.end.0
            })
            .copied()
            .collect()
    }

    /// Total covered base pairs inside `span`.
    #[must_use]
    pub fn accessible_bp(&self, span: Interval0) -> u64 {
        self.merged
            .iter()
            .map(|region| {
                let start = region.start.0.max(span.start.0);
                let end = region.end.0.min(span.end.0);
                end.saturating_sub(start)
            })
            .sum()
    }
}

/// Accessibility mask: the callable/accessible regions of a genome
/// (`docs/rounds/v1/design.md` §12.1); makes `AccessibleBp` denominators use covered bp.
#[derive(Clone, Debug, Default)]
pub struct AccessibilityMask(pub Arc<RegionSet>);

impl AccessibilityMask {
    /// Accessible bp inside a window.
    #[must_use]
    pub fn accessible_bp(&self, span: Interval0) -> u64 {
        self.0.accessible_bp(span)
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bp_windows_match_legacy_anchoring() {
        let grid = WindowGrid::Bp {
            size: 100,
            step: 100,
        };
        // First site at 150 anchors the grid at 100.
        let windows = grid.bp_windows(150, 450, 1000).unwrap();
        let starts: Vec<u64> = windows.iter().map(|w| w.interval.start.0).collect();
        assert_eq!(starts, [100, 200, 300, 400]);
        assert!(windows.iter().all(|w| w.interval.end.0 <= 1000));
        // Every position falls in exactly one window.
        assert!(windows.iter().any(|w| w.interval.contains(Pos0(150))));
        assert!(windows.iter().any(|w| w.interval.contains(Pos0(450))));
        // Bins partition positions exactly.
        let width = grid.bin_width_bp().unwrap();
        assert_eq!(width, 100);
    }

    #[test]
    fn overlapping_bp_windows_share_bins() {
        let grid = WindowGrid::Bp {
            size: 100,
            step: 50,
        };
        assert_eq!(grid.bin_width_bp(), Some(50));
        let windows = grid.bp_windows(10, 200, 1000).unwrap();
        // Anchor 0: windows [0,100), [50,150), [100,200) — site 10 in window 0.
        assert_eq!(windows[0].interval, Interval0::new(Pos0(0), Pos0(100)));
        assert_eq!(windows[0].bins, BinId(0)..BinId(2));
        assert_eq!(windows[1].bins, BinId(1)..BinId(3));
    }

    #[test]
    fn last_window_clips_to_contig_length() {
        let grid = WindowGrid::Bp {
            size: 100,
            step: 100,
        };
        let windows = grid.bp_windows(0, 950, 960).unwrap();
        let last = windows.last().unwrap();
        assert_eq!(last.interval, Interval0::new(Pos0(900), Pos0(960)));
        assert_eq!(last.bins, BinId(9)..BinId(10));
    }

    #[test]
    fn sites_windows_enumerate_ranks() {
        let grid = WindowGrid::Sites {
            n_sites: 10,
            step: 10,
        };
        let windows = grid.sites_windows(0, 25).unwrap();
        assert_eq!(windows.len(), 3);
        assert_eq!(windows[2].first_rank, 20);
        assert_eq!(windows[2].last_rank_exclusive, 30);
        assert!(grid.bp_windows(0, 10, 100).is_none());
    }

    #[test]
    fn regions_validate_and_count_accessible_bp() {
        let set = RegionSet::new(vec![
            Interval0::new(Pos0(100), Pos0(200)),
            Interval0::new(Pos0(150), Pos0(250)),
            Interval0::new(Pos0(300), Pos0(400)),
        ])
        .unwrap();
        assert_eq!(set.windows().len(), 3);
        // Overlap merges: covered = [100,250) + [300,400).
        assert_eq!(set.accessible_bp(Interval0::new(Pos0(0), Pos0(1000))), 250);
        assert_eq!(set.accessible_bp(Interval0::new(Pos0(0), Pos0(150))), 50);
        assert!(RegionSet::new(vec![Interval0::new(Pos0(5), Pos0(5))]).is_err());
    }

    /// Linear-scan `stab` reference (doc_v20 R14): the pre-fix
    /// implementation, kept verbatim as the oracle for the binary-search
    /// rewrite. The bug was complexity, not wrong answers — O(intervals)
    /// per site made an exome BED over WGS ~10^12 comparisons — and
    /// complexity is not observable from a result, so the executable proof
    /// is this exact-equivalence sweep over every boundary (plus
    /// `r14_stab_large_bed_fast`, which is intractable for the scan) and
    /// the removed `TODO(#1)`.
    fn stab_linear_reference(set: &RegionSet, pos0: u64) -> Vec<usize> {
        set.windows()
            .iter()
            .enumerate()
            .filter(|(_, iv)| iv.start.0 <= pos0 && pos0 < iv.end.0)
            .map(|(idx, _)| idx)
            .collect()
    }

    #[test]
    fn r14_stab_matches_linear_reference() {
        let sets: Vec<Vec<Interval0>> = vec![
            // Empty set: every probe is empty, no panic.
            vec![],
            // Single interval.
            vec![Interval0::new(Pos0(100), Pos0(200))],
            // Overlapping pair, unsorted input.
            vec![
                Interval0::new(Pos0(50), Pos0(150)),
                Interval0::new(Pos0(100), Pos0(120)),
            ],
            // Nested + duplicated spans + shared boundaries, unsorted.
            vec![
                Interval0::new(Pos0(10), Pos0(30)),
                Interval0::new(Pos0(15), Pos0(25)),
                Interval0::new(Pos0(10), Pos0(30)),
                Interval0::new(Pos0(0), Pos0(10)),
            ],
            // Adjacent (touching) intervals coalesce into one merged span;
            // a probe at the seam must hit only the right-hand region.
            vec![
                Interval0::new(Pos0(10), Pos0(20)),
                Interval0::new(Pos0(20), Pos0(30)),
                Interval0::new(Pos0(30), Pos0(40)),
            ],
            // One long span overlapping many short ones: a single merged
            // span with several members.
            (0..8)
                .map(|i| Interval0::new(Pos0(100 + i * 10), Pos0(105 + i * 10)))
                .chain([Interval0::new(Pos0(90), Pos0(200))])
                .collect(),
        ];
        for intervals in sets {
            let set = RegionSet::new(intervals).unwrap();
            // Sweep before-first, after-last and around every boundary
            // (start-1, start, end-1, end, end+1).
            let mut probes: Vec<u64> = vec![0, u64::MAX];
            for iv in set.windows() {
                probes.extend_from_slice(&[
                    iv.start.0.saturating_sub(1),
                    iv.start.0,
                    iv.end.0 - 1,
                    iv.end.0,
                    iv.end.0.saturating_add(1),
                ]);
            }
            probes.sort_unstable();
            probes.dedup();
            for pos0 in probes {
                assert_eq!(
                    set.stab(pos0),
                    stab_linear_reference(&set, pos0),
                    "intervals {:?}, pos0 {pos0}",
                    set.windows()
                );
            }
        }
    }

    /// doc_v20 R14 smoke: an exome-scale BED (200k spans) probed 2k times
    /// must finish in milliseconds; the linear scan cost ~4×10^8
    /// comparisons (minutes in a debug build). The 10 s ceiling is
    /// generous so slow CI cannot flake — the scan overshoots it by an
    /// order of magnitude, the binary search sits far below it.
    #[test]
    fn r14_stab_large_bed_fast() {
        let intervals: Vec<Interval0> = (0..200_000_u64)
            .map(|i| Interval0::new(Pos0(i * 10), Pos0(i * 10 + 5)))
            .collect();
        let set = RegionSet::new(intervals).unwrap();
        // One exact probe pins the binary search mid-array: span i covers
        // [10i, 10i+5), so 1_000_003 lands in input interval 100_000.
        assert_eq!(set.stab(1_000_003), vec![100_000]);
        let start = std::time::Instant::now();
        let mut hits = 0_usize;
        for q in 0..2_000_u64 {
            let pos0 = q * 1_000;
            if !set.stab(pos0).is_empty() {
                hits += 1;
            }
        }
        let elapsed = start.elapsed();
        assert!(hits > 0, "the probe sweep must land inside some spans");
        assert!(
            elapsed < std::time::Duration::from_secs(10),
            "stab sweep took {elapsed:?} (linear scan: ~4×10^8 comparisons)"
        );
    }

    #[test]
    fn validate_rejects_zero_size_or_step() {
        assert!(WindowGrid::Bp { size: 0, step: 10 }.validate().is_err());
        assert!(
            WindowGrid::Sites {
                n_sites: 10,
                step: 0
            }
            .validate()
            .is_err()
        );
        assert!(WindowGrid::Bp { size: 10, step: 10 }.validate().is_ok());
    }
}
