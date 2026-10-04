//! Aggregate allele counts per group.

use crate::population::PopulationView;
use crate::window::WindowBuffer;

/// Aggregate allele tallies for one group over one window.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AlleleCountsResult {
    pub group: u32,
    pub ref_count: u64,
    pub alt_count: u64,
    /// Missing alleles summed over sites (per-allele deletion).
    pub missing: u64,
    pub n_variant_sites: u64,
    pub n_invariant_sites: u64,
}

/// Computes aggregate allele counts for every group.
pub fn allele_counts(window: &WindowBuffer, pops: &PopulationView) -> Vec<AlleleCountsResult> {
    (0..pops.n_groups() as u32)
        .map(|group| {
            let samples = pops.samples(group);
            let mut result = AlleleCountsResult {
                group,
                ..Default::default()
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
        .collect()
}
