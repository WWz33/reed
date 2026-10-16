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
    /// Whether `drop_stale` came from an explicit `-x` (the CLI sets it)
    /// rather than from the default list. The automatic AC/AN update treats
    /// its own entries as exempt on GT-less records (there is no value to
    /// put back, doc_v12 D12-2); an explicit `-x INFO/AC` is the user's
    /// request and applies unconditionally.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub drop_stale_explicit: bool,
    /// bcftools `view -I`: do not update `INFO/AC,INFO/AN` when samples are
    /// dropped (the update is the default, 02-cli.md §3.5).
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub no_update: bool,
}

impl Default for RecomputeSpec {
    fn default() -> Self {
        Self {
            tags: Vec::new(),
            drop_stale: default_drop_stale(),
            drop_stale_explicit: false,
            no_update: false,
        }
    }
}

fn default_drop_stale() -> Vec<String> {
    [
        "AC",
        "AN",
        "AF",
        "MLEAC",
        "MLEAF",
        "ExcessHet",
        "InbreedingCoeff",
        "DP",
    ]
    .iter()
    .map(ToString::to_string)
    .collect()
}

/// The `--fill-tags` tags this build implements, in canonical spelling
/// (bcftools `+fill-tags` matches case-insensitively, `strcasecmp`,
/// plugins/fill-tags.c:505-537). AF and F_MISSING are in bcftools' set and
/// in `cli.md`'s prose but are not implemented here (doc_v12 N10).
pub const SUPPORTED_FILL_TAGS: [&str; 5] = ["AC", "AN", "NS", "MAF", "HWE"];

/// Normalises one `--fill-tags` entry to its canonical spelling, or reports
/// it as unsupported. bcftools refuses the unknown ones with
/// `the tag "X" is not supported` and exits (plugins/fill-tags.c:530-535);
/// silently accepting them used to turn `--fill-tags AF` into a no-op that
/// then *deleted* the input's AF via the default `-x` (doc_v12 D12-11).
///
/// # Errors
///
/// A message naming the tag and listing the supported set.
pub fn canonical_fill_tag(tag: &str) -> Result<&'static str, String> {
    let tag = tag.trim();
    SUPPORTED_FILL_TAGS
        .iter()
        .find(|t| t.eq_ignore_ascii_case(tag))
        .copied()
        .ok_or_else(|| {
            format!(
                "unsupported --fill-tags tag {tag:?}; supported: {}",
                SUPPORTED_FILL_TAGS.join(", ")
            )
        })
}

impl RecomputeSpec {
    /// The gold conf default: every supported tag.
    #[must_use]
    pub fn all() -> Self {
        Self {
            tags: ["AC", "AN", "NS", "MAF", "HWE"]
                .iter()
                .map(ToString::to_string)
                .collect(),
            drop_stale: default_drop_stale(),
            drop_stale_explicit: false,
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
        self.drop_stale
            .iter()
            .map(String::as_bytes)
            .map(<[u8]>::to_vec)
            .collect()
    }
}
