//! Sample-level rules: the static subset (`-s/-S/--keep/--remove`) and the
//! dynamic `--mind` threshold. They are modelled separately (04 §C0b) so
//! both can be given at once: the static subset is applied first, `--mind`
//! second (plink2's order).

use serde::{Deserialize, Serialize};

/// One sample-list source.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum SampleList {
    /// Inline names (`-s A,B`, or the exclusion form `-s ^A,B`).
    #[serde(rename = "names")]
    Names(Vec<String>),
    /// One name per line (`-S FILE`, `--keep FILE`, `--remove FILE`).
    #[serde(rename = "file")]
    File(String),
}

/// The static sample subset: keep ∩ (not remove). `None` fields mean "all"
/// and "none" respectively.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SampleSelect {
    /// Samples to keep (`-s`, `-S`, `--keep`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub keep: Option<SampleList>,
    /// Samples to drop (`-s ^…`, `--remove`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remove: Option<SampleList>,
}

impl SampleSelect {
    /// Whether the selection is the identity (no samples touched).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.keep.is_none() && self.remove.is_none()
    }
}
