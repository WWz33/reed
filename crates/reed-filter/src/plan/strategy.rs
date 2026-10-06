//! Execution strategies (§11.4). Strategy choice affects speed only: all
//! strategies produce byte-identical output and funnels for the same input
//! (locked by `plan_strategy` tests and gold).

use serde::{Deserialize, Serialize};

/// The three strategies. `SketchThenEmit` selection arrives with M5b; the
/// variant already exists so M5b changes no public shape.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExecStrategy {
    /// One full pass: head → mask → tail → edit → compress, all in the
    /// worker. Triggered when no cohort rule is active, or the cohort is an
    /// explicit keep/remove list.
    #[default]
    #[serde(rename = "single_pass")]
    SinglePass,
    /// `max_missing` cohort with the sketch fitting memory (M5b).
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
/// and `n_rules` mask rules: `(4 + R) × W` words per head-surviving site,
/// plus the head bitmap — assuming every record survives the head stage
/// (an upper bound, §11.4).
#[must_use]
pub fn sketch_size_bytes(n_records: u64, n_samples: usize, n_rules: usize) -> u64 {
    let w = u64::try_from(n_samples.div_ceil(64)).unwrap_or(u64::MAX);
    n_records.saturating_mul((4 + u64::try_from(n_rules).unwrap_or(u64::MAX)) * w * 8)
        .saturating_add(n_records.div_ceil(64) * 8)
}

/// Picks the strategy for a `max_missing` cohort: the sketch path when the
/// upper-bounded sketch fits the memory budget, TwoPassReparse otherwise
/// (§11.4; the record count comes from the first-chunk records/byte
/// extrapolation, spec-sanctioned when no index is available).
#[must_use]
pub fn choose_for_max_missing(
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
