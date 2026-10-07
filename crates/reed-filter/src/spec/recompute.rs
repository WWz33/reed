//! Recompute specification (§11.1 `[recompute]`, semantics §11.6).

use serde::{Deserialize, Serialize};

/// Which INFO tags are recomputed after masking, and which stale tags are
/// dropped from every emitted record.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RecomputeSpec {
    /// Subset of `AC | AN | NS | MAF | HWE` (§11.6 semantics).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tags: Vec<String>,
    /// Tags removed from every record regardless of `tags` (the gold
    /// `annotate -x` list by default).
    #[serde(default = "default_drop_stale", skip_serializing_if = "Vec::is_empty")]
    pub drop_stale: Vec<String>,
    /// bcftools `view -I`: do not update `INFO/AC,INFO/AN` when samples are
    /// dropped (the update is the default, 02-cli.md §3.5).
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub no_update: bool,
}

impl Default for RecomputeSpec {
    fn default() -> Self {
        Self { tags: Vec::new(), drop_stale: default_drop_stale(), no_update: false }
    }
}

fn default_drop_stale() -> Vec<String> {
    ["AC", "AN", "AF", "MLEAC", "MLEAF", "ExcessHet", "InbreedingCoeff", "DP"]
        .iter()
        .map(ToString::to_string)
        .collect()
}

impl RecomputeSpec {
    /// The gold conf default: every supported tag.
    #[must_use]
    pub fn all() -> Self {
        Self {
            tags: ["AC", "AN", "NS", "MAF", "HWE"].iter().map(ToString::to_string).collect(),
            drop_stale: default_drop_stale(),
            no_update: false,
        }
    }

    /// Whether a tag is requested.
    #[must_use]
    pub fn wants(&self, tag: &str) -> bool {
        self.tags.iter().any(|t| t == tag)
    }

    /// The stale list as byte slices for the editor.
    #[must_use]
    pub fn stale_bytes(&self) -> Vec<Vec<u8>> {
        self.drop_stale.iter().map(String::as_bytes).map(<[u8]>::to_vec).collect()
    }
}
