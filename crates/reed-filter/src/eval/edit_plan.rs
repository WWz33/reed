//! Survivor → EditPlan assembly (§4 eval layer): the recomputed tag values
//! rendered from the post-mask counts, in fill-tags order. The caller
//! (worker) owns the scratch buffers and constructs the borrowed
//! edit plan inline — no leaked allocations, no unsafe.

// Casts in this file are bounded by construction: INFO Float values are f32 by the VCF header type.
#![allow(clippy::cast_possible_truncation)]

use reed_core::Bitset;
use reed_core::count::SiteCounts;
use reed_core::math::hwe::HweCache;
use reed_core::sample::SampleMask;
use reed_core::vcf::line_edit::{EditPlan, InfoValueOut};

use crate::plan::RecomputeFlags;

/// One survivor's edits as a pure value (doc_v30 R14): the INFO tags to
/// drop, the recomputed values to write, and the sample keep/mask planes.
///
/// The three emit paths (single pass, sketch emit, filtered passthrough)
/// used to assemble `EditPlan { .. }` inline with the same four fields;
/// this value is the single constructor, so a field added to the edit
/// cannot be forgotten in one path.
#[derive(Debug)]
pub struct SurvivorEdit<'a> {
    /// `None` = keep every sample column.
    pub keep_samples: Option<&'a SampleMask>,
    /// Samples whose GT subfield becomes `./.`.
    pub mask_gt: Option<&'a Bitset>,
    /// Stale INFO tags to drop from this record.
    pub info_remove: &'a [&'a [u8]],
    /// Recomputed tags (fill order).
    pub info_set: &'a [(&'a [u8], InfoValueOut)],
}

impl<'a> SurvivorEdit<'a> {
    /// The byte-level edit to apply.
    #[must_use]
    pub fn edit_plan(&self) -> EditPlan<'a> {
        EditPlan {
            keep_samples: self.keep_samples,
            mask_gt: self.mask_gt,
            info_remove: self.info_remove,
            info_set: self.info_set,
        }
    }
}

/// Fills the worker's reusable `info_set` with this site's recomputed
/// values (fill-tags order: AC, AN, NS, MAF, HWE). `L_tags` compares
/// values, not byte positions, but a stable order keeps outputs
/// reproducible. `hwe` memoizes the exact-test distribution per
/// `(n_called, mac)` (A6-16); the p-values are identical to the uncached
/// `SiteCounts::hwe_p`.
///
/// A record whose FORMAT carries no GT key has no genotype-derived counts
/// to recompute from: `counts` is all zeros there, and writing them would
/// replace the input's AC/AN with `AC=0;AN=0`. htslib takes the same way
/// out in both of its update sites — `vcfview.c:343`
/// (`!bcf_get_fmt(…,"GT") → update_ac = 0`) and `plugins/fill-tags.c:763`
/// (`if (!fmt_gt) return;`) — so a GT-less record keeps its input values
/// (doc_v12 D12-2).
pub fn fill_info_set(
    counts: &SiteCounts,
    recompute: &RecomputeFlags,
    hwe: &mut HweCache,
    info_set: &mut Vec<(&'static [u8], InfoValueOut)>,
) {
    info_set.clear();
    if !counts.has_gt {
        return;
    }
    for tag in &recompute.tags {
        info_set.push((tag.name.as_bytes(), (tag.compute)(counts, hwe)));
    }
}

/// The stale tags for one record: borrowed from the plan on the common path
/// (GT present), owned only when the automatic AC/AN entries must be kept
/// (no GT; doc_v12 D12-2). This used to allocate a `Vec` per record even in
/// the borrowed case (doc_v24 V24-07).
#[derive(Debug)]
pub enum StaleTags<'s, 'a> {
    /// The plan's full list.
    Borrowed(&'s [&'a [u8]]),
    /// The list minus the automatic update entries.
    Owned(Vec<&'a [u8]>),
}

impl<'a> StaleTags<'_, 'a> {
    /// The tags to hand to the edit plan.
    #[must_use]
    pub fn as_slice(&self) -> &[&'a [u8]] {
        match self {
            Self::Borrowed(tags) => tags,
            Self::Owned(tags) => tags,
        }
    }
}

/// The stale tags to drop from one record: the plan's list, minus the
/// automatic AC/AN update entries when the record has no GT (there is no
/// value to put back, so dropping them would destroy the input's INFO —
/// doc_v12 D12-2). A user's explicit `-x INFO/AC` still applies: the caller
/// passes the tags that came from `--fill-tags`/`-x` separately.
#[must_use]
pub fn stale_for_record<'s, 'a>(
    stale: &'s [&'a [u8]],
    auto_update: &[&'a [u8]],
    has_gt: bool,
) -> StaleTags<'s, 'a> {
    if has_gt {
        return StaleTags::Borrowed(stale);
    }
    StaleTags::Owned(
        stale
            .iter()
            .copied()
            .filter(|t| !auto_update.iter().any(|a| a == t))
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// doc_v30 R14: the value carries exactly the four edit planes through
    /// to `EditPlan` — no field is dropped or reordered on the way.
    #[test]
    fn survivor_edit_maps_every_plane_to_the_edit_plan() {
        let set = [(b"AC".as_slice(), InfoValueOut::Int(2))];
        let remove: [&[u8]; 1] = [b"DP".as_slice()];
        let edit = SurvivorEdit {
            keep_samples: None,
            mask_gt: None,
            info_remove: &remove,
            info_set: &set,
        };
        let plan = edit.edit_plan();
        assert_eq!(plan.info_set, &set);
        assert_eq!(plan.info_remove, &remove);
        assert!(plan.keep_samples.is_none());
        assert!(plan.mask_gt.is_none());
    }
}
