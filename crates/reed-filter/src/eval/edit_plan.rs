//! Survivor → EditPlan assembly (§4 eval layer): the recomputed tag values
//! rendered from the post-mask counts, in fill-tags order. The caller
//! (worker) owns the scratch buffers and constructs the borrowed
//! edit plan inline — no leaked allocations, no unsafe.

// Casts in this file are bounded by construction: INFO Float values are f32 by the VCF header type.
#![allow(clippy::cast_possible_truncation)]

use reed_core::count::SiteCounts;
use reed_core::math::hwe::HweCache;
use reed_io::vcf::line_edit::InfoValueOut;

use crate::plan::RecomputeFlags;

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
    if recompute.ac {
        info_set.push((b"AC".as_slice(), InfoValueOut::Int(counts.ac as i64)));
    }
    if recompute.an {
        info_set.push((b"AN".as_slice(), InfoValueOut::Int(counts.an as i64)));
    }
    if recompute.ns {
        info_set.push((b"NS".as_slice(), InfoValueOut::Int(counts.ns as i64)));
    }
    if recompute.maf {
        info_set.push((b"MAF".as_slice(), InfoValueOut::Float(counts.maf() as f32)));
    }
    if recompute.hwe {
        let p = counts
            .hwe_test_input()
            .map_or(1.0, |(n, mac, het)| hwe.p_exact(n, mac, het));
        info_set.push((b"HWE".as_slice(), InfoValueOut::Float(p as f32)));
    }
}

/// The stale tags to drop from one record: the plan's list, minus the
/// automatic AC/AN update entries when the record has no GT (there is no
/// value to put back, so dropping them would destroy the input's INFO —
/// doc_v12 D12-2). A user's explicit `-x INFO/AC` still applies: the caller
/// passes the tags that came from `--fill-tags`/`-x` separately.
#[must_use]
pub fn stale_for_record<'a>(
    stale: &[&'a [u8]],
    auto_update: &[&'a [u8]],
    has_gt: bool,
) -> Vec<&'a [u8]> {
    if has_gt {
        return stale.to_vec();
    }
    stale
        .iter()
        .copied()
        .filter(|t| !auto_update.iter().any(|a| a == t))
        .collect()
}
