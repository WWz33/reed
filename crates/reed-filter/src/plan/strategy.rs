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
