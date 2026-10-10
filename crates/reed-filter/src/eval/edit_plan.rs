//! Survivor → EditPlan assembly (§4 eval layer): the recomputed tag values
//! rendered from the post-mask counts, in fill-tags order. The caller
//! (worker) owns the scratch buffers and constructs the borrowed
//! edit plan inline — no leaked allocations, no unsafe.

// Casts in this file are bounded by construction: INFO Float values are f32 by the VCF header type.
#![allow(clippy::cast_possible_truncation)]

use reed_core::count::SiteCounts;
use reed_io::vcf::line_edit::InfoValueOut;

use crate::plan::RecomputeFlags;

/// Fills the worker's reusable `info_set` with this site's recomputed
/// values (fill-tags order: AC, AN, NS, MAF, HWE). `L_tags` compares
/// values, not byte positions, but a stable order keeps outputs
/// reproducible.
pub fn fill_info_set(
    counts: &SiteCounts,
    recompute: &RecomputeFlags,
    info_set: &mut Vec<(&'static [u8], InfoValueOut)>,
) {
    info_set.clear();
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
        info_set.push((
            b"HWE".as_slice(),
            InfoValueOut::Float(counts.hwe_p() as f32),
        ));
    }
}
