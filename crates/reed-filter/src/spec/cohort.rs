//! Sample-level rules (§11.1 `[cohort]`): the pass0 global reduction.

use serde::{Deserialize, Serialize};

/// Which samples survive the cohort stage.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum CohortRule {
    /// Drop samples whose genome-wide raw missing rate exceeds the value.
    #[serde(rename = "max_missing")]
    MaxMissing(f64),
    /// Keep exactly the listed samples (one per line; `#` comments allowed).
    #[serde(rename = "keep_list")]
    KeepList(String),
    /// Keep exactly these sample names (bcftools' inline `-s` list).
    #[serde(rename = "keep_names")]
    KeepNames(Vec<String>),
    /// Drop the listed samples.
    #[serde(rename = "remove_list")]
    RemoveList(String),
}
