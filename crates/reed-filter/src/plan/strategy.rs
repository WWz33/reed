//! Execution strategies (§11.4). Strategy choice affects speed only: all
//! strategies produce byte-identical output and funnels for the same input
//! (locked by `plan_strategy` tests and gold).

use serde::{Deserialize, Serialize};

/// The three strategies. `SketchThenEmit` selection arrives with M5b; the
/// variant already exists so M5b changes no public shape.
///
/// Test-only override (doc_v11 N10): production code never sets it — the
/// driver picks the strategy from the plan. Hidden from the docs so the
/// facade's published surface does not carry an internal hook.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExecStrategy {
    /// One full pass: head → mask → tail → edit → compress, all in the
    /// worker. Triggered when no cohort rule is active, or the cohort is an
    /// explicit keep/remove list.
    #[default]
    #[serde(rename = "single_pass")]
    SinglePass,
    /// A dynamic sample-cohort rule with the sketch fitting memory (M5b).
    #[serde(rename = "sketch_then_emit")]
    SketchThenEmit,
    /// Two full passes: pass0 reduces per-sample missing rates; pass1 is a
    /// SinglePass over the surviving sample columns.
    #[serde(rename = "two_pass_reparse")]
    TwoPassReparse,
}

impl ExecStrategy {
    /// The stable name (provenance + report).
    #[must_use]
    pub const fn id(self) -> &'static str {
        match self {
            Self::SinglePass => "single_pass",
            Self::SketchThenEmit => "sketch_then_emit",
            Self::TwoPassReparse => "two_pass_reparse",
        }
    }
}

/// Sketch size upper bound for `n_records` records, `n_samples` samples
/// and `n_rules` mask rules: `(5 + R) × W` words per head-surviving site,
/// plus the head bitmap — assuming every record survives the head stage
/// (an upper bound, §11.4). The fifth word row is the per-sample mask
/// plane the pruner adds (doc_v15 D15-16).
#[must_use]
pub fn sketch_size_bytes(n_records: u64, n_samples: usize, n_rules: usize) -> u64 {
    let w = u64::try_from(n_samples.div_ceil(64)).unwrap_or(u64::MAX);
    n_records
        .saturating_mul((5 + u64::try_from(n_rules).unwrap_or(u64::MAX)) * w * 8)
        .saturating_add(n_records.div_ceil(64) * 8)
}

/// Picks the strategy for a `--mind` cohort: the sketch path
/// when the upper-bounded sketch fits the memory budget, TwoPassReparse
/// otherwise (§11.4; the record count comes from the first-chunk
/// records/byte extrapolation, spec-sanctioned when no index is available).
#[must_use]
pub fn choose_for_mind(
    estimated_records: u64,
    n_samples: usize,
    n_rules: usize,
    mem_budget_bytes: u64,
) -> ExecStrategy {
    if sketch_size_bytes(estimated_records, n_samples, n_rules) <= mem_budget_bytes {
        ExecStrategy::SketchThenEmit
    } else {
        ExecStrategy::TwoPassReparse
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// doc_v15 D15-16: the upper bound counts the pruner's per-sample mask
    /// plane, i.e. `(5 + R)` words per site, not `(4 + R)`.
    #[test]
    fn sketch_model_includes_the_mask_plane() {
        // One 64-bit word of samples: (5 + 0) words for the site, one word
        // for the head bitmap.
        assert_eq!(sketch_size_bytes(1, 64, 0), 5 * 8 + 8);
        // Two sample words per row with 65 samples.
        assert_eq!(sketch_size_bytes(1, 65, 0), 5 * 2 * 8 + 8);
        // Two mask rules add two word rows.
        assert_eq!(
            sketch_size_bytes(1000, 64, 2),
            1000 * 7 * 8 + 1000_u64.div_ceil(64) * 8
        );
    }

    /// The strategy decision switches exactly at the estimate (the `<=`
    /// boundary), which is what `--mind`'s fallback rides on.
    #[test]
    fn choose_for_mind_switches_at_the_estimate() {
        let estimate = sketch_size_bytes(10_000, 128, 1);
        assert_eq!(
            choose_for_mind(10_000, 128, 1, estimate),
            ExecStrategy::SketchThenEmit
        );
        assert_eq!(
            choose_for_mind(10_000, 128, 1, estimate - 1),
            ExecStrategy::TwoPassReparse
        );
    }
}
