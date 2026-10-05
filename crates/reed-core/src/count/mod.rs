//! Per-site, per-group integer counts (`docs/architecture-spec.md` §5.6).
//!
//! [`SiteGroupCounts`] is the single shared input of every statistic: π,
//! dxy, Fst, PBS, SFS, θ, Tajima's D, f-statistics and HWE all derive from
//! these six integers, so each site is decoded exactly once (fixing legacy
//! issue A4). The legacy VCF-tag view (`SiteCounts` in the old
//! `estimators::site_counts`) becomes a derived view in reed-filter.

pub mod kernel;

use crate::genotype::PlaneBlock;
use crate::sample::SampleMask;

pub use kernel::{count_site_group, count_site_group_scalar};

/// Genotype and allele tallies for one group at one site (24 bytes).
///
/// Derived quantities (not stored, to keep the struct at 6 × `u32`):
/// `n_called_alleles = 2·n_group − missing`, `n_full_called = hom_ref + het +
/// hom_alt`, `missing = n_half_called + 2·n_both_missing`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SiteGroupCounts {
    /// Called REF allele count.
    pub n_ref_alleles: u32,
    /// Called ALT allele count (under the site's multi-allelic policy).
    pub n_alt_alleles: u32,
    /// Full calls with two REF alleles.
    pub n_hom_ref: u32,
    /// Full calls with one REF and one ALT allele.
    pub n_het: u32,
    /// Full calls with two ALT alleles.
    pub n_hom_alt: u32,
    /// Genotypes with exactly one called allele (half calls).
    pub n_half_called: u32,
}

impl SiteGroupCounts {
    /// Missing allele count, given the group's selected-sample count.
    #[must_use]
    pub const fn n_missing_alleles(&self, n_group: u32) -> u32 {
        2 * n_group - self.n_ref_alleles - self.n_alt_alleles
    }

    /// Called allele count, given the group's selected-sample count.
    #[must_use]
    pub const fn n_called_alleles(&self, n_group: u32) -> u32 {
        2 * n_group - self.n_missing_alleles(n_group)
    }

    /// Genotypes with both alleles called.
    #[must_use]
    pub const fn n_full_called(&self) -> u32 {
        self.n_hom_ref + self.n_het + self.n_hom_alt
    }
}

/// Column-major store of [`SiteGroupCounts`]: `counts[site * n_groups +
/// group]`.
#[derive(Clone, Debug, Default)]
pub struct CountTable {
    n_groups: u16,
    counts: Vec<SiteGroupCounts>,
}

impl CountTable {
    /// Allocates a zeroed table for `n_sites` sites and `n_groups` groups.
    #[must_use]
    pub fn new(n_sites: usize, n_groups: usize) -> Self {
        Self {
            n_groups: u16::try_from(n_groups).unwrap_or(u16::MAX),
            counts: vec![SiteGroupCounts::default(); n_sites * n_groups],
        }
    }

    /// Number of sites.
    #[must_use]
    pub fn n_sites(&self) -> usize {
        if self.n_groups == 0 {
            0
        } else {
            self.counts.len() / self.n_groups as usize
        }
    }

    /// Number of groups.
    #[must_use]
    pub const fn n_groups(&self) -> usize {
        self.n_groups as usize
    }

    /// The counts of one site across all groups.
    ///
    /// # Panics
    ///
    /// Panics if `site >= n_sites` (caller invariant).
    #[must_use]
    pub fn row(&self, site: usize) -> &[SiteGroupCounts] {
        let n_groups = self.n_groups as usize;
        assert!(site * n_groups + n_groups <= self.counts.len(), "invariant: site index out of range");
        &self.counts[site * n_groups..site * n_groups + n_groups]
    }

    /// Mutable row of one site.
    ///
    /// # Panics
    ///
    /// Panics if `site >= n_sites`.
    pub fn row_mut(&mut self, site: usize) -> &mut [SiteGroupCounts] {
        let n_groups = self.n_groups as usize;
        assert!(site * n_groups + n_groups <= self.counts.len(), "invariant: site index out of range");
        &mut self.counts[site * n_groups..site * n_groups + n_groups]
    }

    /// The count of one (site, group) pair.
    #[must_use]
    pub fn get(&self, site: usize, group: usize) -> Option<SiteGroupCounts> {
        self.row(site).get(group).copied()
    }

    /// Fills `sites` (0-based, ascending, within the block) from the block's
    /// planes over the given group masks.
    ///
    /// # Panics
    ///
    /// Panics if the site range exceeds either the table or the block.
    pub fn fill(&mut self, sites: std::ops::Range<usize>, block: &PlaneBlock, masks: &[SampleMask]) {
        assert!(
            sites.end <= self.n_sites() && sites.end <= block.n_sites(),
            "invariant: fill range exceeds table or block"
        );
        assert_eq!(self.n_groups(), masks.len(), "invariant: mask count matches table groups");
        for site_idx in sites {
            let row = self.row_mut(site_idx);
            let planes = block.site_planes(site_idx);
            for (mask, counts) in masks.iter().zip(row.iter_mut()) {
                count_site_group(&planes, mask.words(), mask.n_selected(), counts);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype::write_genotype;
    use crate::id::GroupId;
    use crate::sample::{SampleMask, SampleTable};

    /// Builds a 4-sample block from explicit genotypes, one site per row of
    /// `genotypes` (each entry = `(slot0, slot1)` with `None` = missing).
    fn block_from(genotypes: &[[(Option<u16>, Option<u16>); 4]]) -> PlaneBlock {
        let mut block = PlaneBlock::new(genotypes.len(), 4, false);
        let policies = crate::genotype::GenotypePolicies::default();
        let mut counters = crate::genotype::PolicyCounters::default();
        for (site_idx, site) in genotypes.iter().enumerate() {
            for (sample_idx, &(a, b)) in site.iter().enumerate() {
                write_genotype(
                    &mut block,
                    site_idx,
                    sample_idx,
                    [a, b],
                    crate::genotype::Ploidy::Diploid,
                    false,
                    &policies,
                    &mut counters,
                )
                .expect("invariant: test genotypes are biallelic diploid");
            }
        }
        block
    }

    #[test]
    fn counts_match_hand_computation() {
        // site 0: s0=0/1 s1=1/1 s2=./. s3=0/0
        // site 1: s0=1/. s1=0/0 s2=0/1 s3=1/1
        let block = block_from(&[
            [(Some(0), Some(1)), (Some(1), Some(1)), (None, None), (Some(0), Some(0))],
            [(Some(1), None), (Some(0), Some(0)), (Some(0), Some(1)), (Some(1), Some(1))],
        ]);
        let samples = SampleTable::from_names(["s0", "s1", "s2", "s3"]).unwrap();
        let groups = crate::group::GroupTable::from_assignment(
            &samples,
            &[("s0", "A"), ("s1", "A"), ("s2", "B"), ("s3", "B")],
        )
        .unwrap();
        let masks = groups.masks();

        let mut table = CountTable::new(2, 2);
        table.fill(0..2, &block, &masks);

        // Group A at site 0: s0=0/1, s1=1/1 → ref=1, alt=3, het=1, hom_alt=1.
        let a0 = table.get(0, 0).unwrap();
        assert_eq!((a0.n_ref_alleles, a0.n_alt_alleles), (1, 3));
        assert_eq!((a0.n_het, a0.n_hom_alt, a0.n_hom_ref), (1, 1, 0));
        assert_eq!(a0.n_full_called(), 2);
        assert_eq!(a0.n_missing_alleles(2), 0);
        // Group B at site 0: s2=./., s3=0/0 → ref=2, alt=0, both-missing=1.
        let b0 = table.get(0, 1).unwrap();
        assert_eq!((b0.n_ref_alleles, b0.n_alt_alleles), (2, 0));
        assert_eq!(b0.n_full_called(), 1);
        assert_eq!(b0.n_missing_alleles(2), 2);
        // Group A at site 1: s0=1/. (half call), s1=0/0 → alt=1, ref=2, half=1.
        let a1 = table.get(1, 0).unwrap();
        assert_eq!((a1.n_ref_alleles, a1.n_alt_alleles), (2, 1));
        assert_eq!(a1.n_half_called, 1);
        assert_eq!(a1.n_full_called(), 1);
    }

    #[test]
    fn group_id_lookup_via_table() {
        let samples = SampleTable::from_names(["a", "b"]).unwrap();
        let groups = crate::group::GroupTable::from_assignment(&samples, &[("a", "P")]).unwrap();
        assert_eq!(groups.by_name("P"), Some(GroupId(0)));
    }

    #[test]
    fn all_ref_site_counts() {
        let block = block_from(&[
            [(Some(0), Some(0)), (Some(0), Some(0)), (Some(0), Some(0)), (Some(0), Some(0))],
        ]);
        let mask = SampleMask::all(4);
        let mut table = CountTable::new(1, 1);
        table.fill(0..1, &block, &[mask]);
        let c = table.get(0, 0).unwrap();
        assert_eq!((c.n_ref_alleles, c.n_alt_alleles), (8, 0));
        assert_eq!((c.n_hom_ref, c.n_het, c.n_hom_alt), (4, 0, 0));
        assert_eq!(c.n_half_called, 0);
    }
}
