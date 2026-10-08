//! Grid integration properties: windows and bins agree with each other and
//! with the legacy anchoring (`docs/architecture-spec.md` §19 M1, grid half).

use proptest::prelude::*;

use reed_core::grid::{BinId, WindowGrid};
use reed_core::coord::Pos0;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    #[test]
    fn bp_windows_and_bins_partition_the_touched_range(
        size in 1_u64..500,
        step in 1_u64..500,
        first_pos0 in 0_u64..2_000,
        last_delta in 0_u64..2_000,
        contig_len in 2_000_u64..8_000,
    ) {
        let last_pos0 = (first_pos0 + last_delta).min(contig_len - 1);
        let grid = WindowGrid::Bp { size, step };
        let width = grid.bin_width_bp().unwrap();
        let windows = grid.bp_windows(first_pos0, last_pos0, contig_len).unwrap();
        prop_assert!(!windows.is_empty());

        // Windows are ordered, within the contig, and contain both endpoints.
        for pair in windows.windows(2) {
            prop_assert!(pair[0].interval.start < pair[1].interval.start);
        }
        // With step <= size the windows cover the whole touched span; with
        // step > size a position may legitimately fall in a gap (legacy
        // behaviour: gap sites belong to no window).
        if step <= size {
            prop_assert!(windows.first().unwrap().interval.contains(Pos0(first_pos0)));
            prop_assert!(windows.last().unwrap().interval.contains(Pos0(last_pos0)));
        } else {
            prop_assert!(windows.first().unwrap().interval.start.0 <= first_pos0);
        }

        // With a non-overlapping grid (step >= size) every position maps to
        // exactly one window and the bin ranges partition the touched span.
        if step >= size {
            for (i, w) in windows.iter().enumerate() {
                prop_assert_eq!(w.index as usize, i);
                prop_assert!(w.bins.start <= w.bins.end);
                let covered = w.bins.end.0 - w.bins.start.0;
                prop_assert!(covered <= size / width + 1);
                if i + 1 < windows.len() {
                    prop_assert!(w.bins.end <= windows[i + 1].bins.start,
                        "bins must not overlap across windows");
                }
            }
            // Positions inside one bin share the bin index (probe the first
            // bin, which is aligned by construction).
            let anchor = windows[0].interval.start.0;
            let b1 = grid.bin_of_pos0(anchor, anchor).unwrap();
            let b2 = grid.bin_of_pos0(anchor + width - 1, anchor).unwrap();
            prop_assert_eq!(b1, b2);
        }
    }
}

#[test]
fn legacy_anchoring_is_preserved_end_to_end() {
    // The legacy reader anchors the first window at floor(pos/step)*step and
    // clips the final window to the contig length (reed-io lib.rs semantics).
    let grid = WindowGrid::Bp { size: 1_000, step: 1_000 };
    let windows = grid.bp_windows(3_700, 9_900, 10_000).unwrap();
    let spans: Vec<(u64, u64)> = windows
        .iter()
        .map(|w| (w.interval.start.0, w.interval.end.0))
        .collect();
    assert_eq!(
        spans,
        vec![(3_000, 4_000), (4_000, 5_000), (5_000, 6_000), (6_000, 7_000), (7_000, 8_000), (8_000, 9_000), (9_000, 10_000)]
    );
    assert_eq!(windows[0].bins, BinId(3)..BinId(4));
}

#[test]
fn overlapping_windows_share_bins_without_reaccumulation() {
    // step < size: a bin belongs to several windows; the ranges must reflect
    // that sharing so the reducer can merge shared bins once.
    let grid = WindowGrid::Bp { size: 100, step: 50 };
    let windows = grid.bp_windows(0, 100, 10_000).unwrap();
    assert_eq!(windows[0].bins, BinId(0)..BinId(2));
    assert_eq!(windows[1].bins, BinId(1)..BinId(3));
    // Shared middle bin: window 0's last bin is window 1's first bin.
    assert_eq!(windows[0].bins.end.0 - 1, windows[1].bins.start.0);
}
