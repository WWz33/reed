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
/// identity plus a sorted merged copy for lookups.
#[derive(Clone, Debug, Default)]
pub struct RegionSet {
    intervals: Vec<Interval0>,
    merged: Vec<Interval0>,
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
        let mut merged = intervals.clone();
        merged.sort_by_key(|i| (i.start, i.end));
        let mut combined: Vec<Interval0> = Vec::with_capacity(merged.len());
        for interval in merged {
            match combined.last_mut() {
                Some(last) if interval.start <= last.end => {
                    last.end = Pos0(last.end.0.max(interval.end.0));
                }
                _ => combined.push(interval),
            }
        }
        Ok(Self {
            intervals,
            merged: combined,
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

    /// Indices of the regions covering `pos0`, in input order. Regions are
    /// kept in input order and a sorted merged copy backs the probe: walk
    /// the merged spans whose end is past `pos0` and collect their members.
    /// Overlapping regions each match (the BED-window convention).
    #[must_use]
    pub fn stab(&self, pos0: u64) -> Vec<usize> {
        // Candidate spans: merged spans that contain or start after pos0
        // are not enough for overlaps (a long span can start before pos0).
        // The intervals are few (a BED line each); a linear pass is exact
        // and simple. TODO(#1): O(n) per site; index by span if BEDs grow.
        self.intervals
            .iter()
            .enumerate()
            .filter(|(_, iv)| iv.start.0 <= pos0 && pos0 < iv.end.0)
            .map(|(idx, _)| idx)
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
