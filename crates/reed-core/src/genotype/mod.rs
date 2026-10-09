//! Genotype representation: one-byte codes and explicit decoding policies.
//!
//! All consumers (statistics, filtering, QC) share this one decoder so that
//! half calls, haploid genotypes and multi-allelic sites mean the same thing
//! everywhere (`docs/spec/architecture.md` §5.4, fixing legacy issue A6).
//! [`GtCode::from_alleles`] is the only place the policies execute;
//! [`write_genotype`] is the plane-writing adapter on top of it.

pub mod gt_code;
pub mod planes;
pub mod policy;

use crate::error::CoreError;

pub use gt_code::{GtCode, STATE_ALT, STATE_MISSING, STATE_REF};
pub use planes::{
    PlaneBlock, SitePlanes, N_BASE_PLANES, PLANE_PHASED, PLANE_SLOT0_ALT, PLANE_SLOT0_MISS,
    PLANE_SLOT1_ALT, PLANE_SLOT1_MISS,
};
pub use policy::{
    GenotypePolicies, HalfCallPolicy, MultiallelicPolicy, PloidyPolicy, PolicyCounters,
};

/// Ploidy as observed in the input text (the decoder knows whether a genotype
/// had a separator, which allele indices alone cannot express).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ploidy {
    /// One allele (`1` in VCF text).
    Haploid,
    /// Two allele slots (`0/1`, `./1`, `.`).
    Diploid,
}

/// Decodes one genotype and writes it into the planes of `site_idx` for
/// `sample_idx`.
///
/// `slots` holds the allele index per slot with `None` for a missing slot
/// (allele `0` = REF, `1..` = ALT). `ploidy` is what the input text showed:
/// haploid genotypes go through [`PloidyPolicy`], diploid genotypes with one
/// missing slot go through [`HalfCallPolicy`], and any ALT index above 1 goes
/// through [`MultiallelicPolicy`]. Every non-default policy trigger is counted
/// into `counters`.
//
// The signature mirrors the fields of one VCF genotype (block, site, sample,
// two slots, ploidy, phasing, policies, counters); bundling them into a
// struct would only move the count elsewhere without clarifying anything.
#[allow(clippy::too_many_arguments)]
pub fn write_genotype(
    block: &mut PlaneBlock,
    site_idx: usize,
    sample_idx: usize,
    slots: [Option<u16>; 2],
    ploidy: Ploidy,
    phased: bool,
    policies: &GenotypePolicies,
    counters: &mut PolicyCounters,
) -> Result<(), CoreError> {
    let code = GtCode::from_alleles(slots, ploidy, phased, policies, counters)?;
    block.set_code(site_idx, sample_idx, code);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn block() -> PlaneBlock {
        PlaneBlock::new(2, 8, true)
    }

    #[test]
    fn diploid_full_calls_write_planes() {
        let mut b = block();
        let mut c = PolicyCounters::default();
        let p = GenotypePolicies::default();
        write_genotype(&mut b, 0, 0, [Some(0), Some(1)], Ploidy::Diploid, true, &p, &mut c).unwrap();
        let planes = b.site_planes(0);
        assert!(planes.slot0_alt[0] & 1 == 0);
        assert!(planes.slot1_alt[0] & 1 == 1);
        assert!(planes.phased.expect("invariant: phased plane")[0] & 1 == 1);
        assert_eq!(c, PolicyCounters::default());
    }

    #[test]
    fn fully_missing_write_both_miss_planes() {
        let mut b = block();
        let mut c = PolicyCounters::default();
        let p = GenotypePolicies::default();
        write_genotype(&mut b, 0, 3, [None, None], Ploidy::Diploid, false, &p, &mut c).unwrap();
        let planes = b.site_planes(0);
        assert!(planes.slot0_miss[0] >> 3 & 1 == 1);
        assert!(planes.slot1_miss[0] >> 3 & 1 == 1);
    }

    #[test]
    fn half_call_semantics_per_policy() {
        let p_allele = GenotypePolicies::default();
        let p_missing = GenotypePolicies {
            half_call: HalfCallPolicy::TreatAsMissing,
            ..GenotypePolicies::default()
        };
        let mut b = block();
        let mut c = PolicyCounters::default();
        // `./1` keeps the called slot.
        write_genotype(&mut b, 0, 0, [None, Some(1)], Ploidy::Diploid, false, &p_allele, &mut c).unwrap();
        let planes = b.site_planes(0);
        assert!(planes.slot1_alt[0] & 1 == 1);
        assert!(planes.slot0_miss[0] & 1 == 1);
        assert!(planes.slot1_miss[0] & 1 == 0);
        assert_eq!(c.n_half_call_as_missing, 0);
        // Same genotype under TreatAsMissing zeroes both slots.
        let mut b2 = block();
        write_genotype(&mut b2, 0, 0, [None, Some(1)], Ploidy::Diploid, false, &p_missing, &mut c).unwrap();
        let planes = b2.site_planes(0);
        assert!(planes.slot0_miss[0] & 1 == 1);
        assert!(planes.slot1_miss[0] & 1 == 1);
        assert!(planes.slot1_alt[0] & 1 == 0);
        assert_eq!(c.n_half_call_as_missing, 1);
    }

    #[test]
    fn haploid_policies() {
        let mut c = PolicyCounters::default();
        // Default rejects.
        let mut b = block();
        let err = write_genotype(&mut b, 0, 0, [Some(1), None], Ploidy::Haploid, false, &GenotypePolicies::default(), &mut c).unwrap_err();
        assert!(matches!(err, CoreError::PolicyRejected(_)));
        // Half-call conversion.
        let p_half = GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            ..GenotypePolicies::default()
        };
        let mut b = block();
        write_genotype(&mut b, 0, 1, [Some(1), None], Ploidy::Haploid, false, &p_half, &mut c).unwrap();
        let planes = b.site_planes(0);
        assert!(planes.slot0_alt[0] >> 1 & 1 == 1);
        assert!(planes.slot1_miss[0] >> 1 & 1 == 1);
        assert_eq!(c.n_haploid_as_half_call, 1);
        // Homozygous duplication.
        let p_hom = GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHomozygous,
            ..GenotypePolicies::default()
        };
        let mut b = block();
        write_genotype(&mut b, 0, 2, [Some(1), None], Ploidy::Haploid, false, &p_hom, &mut c).unwrap();
        let planes = b.site_planes(0);
        assert!(planes.slot0_alt[0] >> 2 & 1 == 1);
        assert!(planes.slot1_alt[0] >> 2 & 1 == 1);
        assert!(planes.slot1_miss[0] >> 2 & 1 == 0);
        assert_eq!(c.n_haploid_as_homozygous, 1);
    }

    #[test]
    fn multiallelic_policies() {
        let mut c = PolicyCounters::default();
        // Default rejects allele index 2.
        let mut b = block();
        let err = write_genotype(&mut b, 0, 0, [Some(0), Some(2)], Ploidy::Diploid, false, &GenotypePolicies::default(), &mut c).unwrap_err();
        assert!(matches!(err, CoreError::PolicyRejected(_)));
        // CollapseAlt folds it onto the ALT bit and counts the trigger.
        let p_collapse = GenotypePolicies {
            multiallelic: MultiallelicPolicy::CollapseAlt,
            ..GenotypePolicies::default()
        };
        let mut b = block();
        write_genotype(&mut b, 0, 0, [Some(0), Some(2)], Ploidy::Diploid, false, &p_collapse, &mut c).unwrap();
        let planes = b.site_planes(0);
        assert!(planes.slot1_alt[0] & 1 == 1);
        assert_eq!(c.n_multiallelic_collapsed, 1);
    }
}
