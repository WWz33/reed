//! Aggregate allele counts per group.
//!
//! Census semantics: these tallies count exactly what the window contains —
//! there is no denominator convention, hence no manifest gate (§4.2 位点内
//! 统计豁免); sites absent from the input are simply not counted.

use crate::estimators::Unavailable;
use crate::population::PopulationView;
use crate::window::{ContigId, WindowBuffer};

use super::ERR_POP_MISMATCH;

/// Aggregate allele tallies for one group over one window.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AlleleCountsResult {
    pub contig: ContigId,
    pub start: u64,
    pub end: u64,
    pub group: u32,
    pub ref_count: u64,
    pub alt_count: u64,
    /// Missing alleles summed over sites (per-allele deletion).
    pub missing: u64,
    pub n_variant_sites: u64,
    pub n_invariant_sites: u64,
}

/// Computes aggregate allele counts for every group.
pub fn allele_counts(
    window: &WindowBuffer,
    pops: &PopulationView,
) -> Result<Vec<AlleleCountsResult>, Unavailable> {
    // The view must cover exactly the window's samples (§12.2-8).
    pops.validate_against(window)
        .map_err(|_| Unavailable { reason: ERR_POP_MISMATCH })?;
    Ok((0..pops.n_groups() as u32)
        .map(|group| {
            let samples = pops.samples(group);
            let mut result = AlleleCountsResult {
                contig: window.contig,
                start: window.start,
                end: window.end,
                group,
                ref_count: 0,
                alt_count: 0,
                missing: 0,
                n_variant_sites: 0,
                n_invariant_sites: 0,
            };
            for site in 0..window.n_sites() {
                if window.is_variant.get(site) {
                    result.n_variant_sites += 1;
                } else {
                    result.n_invariant_sites += 1;
                }
                let ac = window.allele_counts(site, samples);
                result.ref_count += ac.ref_count;
                result.alt_count += ac.alt_count;
                result.missing += ac.missing;
            }
            result
        })
        .collect())
}
