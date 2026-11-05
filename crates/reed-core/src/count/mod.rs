//! Per-site, per-group integer counts (`docs/spec/architecture.md` §5.6).
//!
//! [`SiteGroupCounts`] is the single shared input of every statistic: π,
//! dxy, Fst, PBS, SFS, θ, Tajima's D, f-statistics and HWE all derive from
//! these six integers, so each site is decoded exactly once (fixing legacy
//! issue A4). The legacy VCF-tag view (`SiteCounts` in the old
//! `estimators::site_counts`) becomes a derived view in reed-filter.

pub mod site_counts;

pub use site_counts::SiteCounts;

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

    /// Adds one decoded genotype: the [`GtCode`](crate::genotype::GtCode)
    /// twin of the text-path counting for the consumers that decode codes
    /// directly (docs/rounds/v2/03 §3.3).
    pub fn add_code(&mut self, code: crate::genotype::GtCode) {
        // The class and allele counts are shared with `SiteCounts::from_codes`
        // (doc_v25 V25-07a); only the destination fields differ.
        use crate::genotype::GtClass;
        self.n_alt_alleles += code.n_alt();
        self.n_ref_alleles += code.n_ref();
        match code.class() {
            GtClass::HomRef => self.n_hom_ref += 1,
            GtClass::Het => self.n_het += 1,
            GtClass::HomAlt => self.n_hom_alt += 1,
            GtClass::HalfCall => self.n_half_called += 1,
            GtClass::Missing => {}
        }
    }
}

/// Whether a decoded genotype carries an ALT allele in either slot.
#[must_use]
pub fn code_has_alt(code: crate::genotype::GtCode) -> bool {
    code.slot_alt(0) || code.slot_alt(1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype::{GtCode, STATE_ALT, STATE_MISSING, STATE_REF};

    fn code(slot0: u8, slot1: u8) -> GtCode {
        GtCode::from_states(slot0, slot1, false, false)
    }

    #[test]
    fn add_code_matches_hand_computation() {
        // s0=0/1 s1=1/1 s2=./. s3=0/0 → ref=3, alt=3, het=1, hom_alt=1,
        // hom_ref=1, both-missing=1.
        let mut c = SiteGroupCounts::default();
        for gt in [
            code(STATE_REF, STATE_ALT),
            code(STATE_ALT, STATE_ALT),
            code(STATE_MISSING, STATE_MISSING),
            code(STATE_REF, STATE_REF),
        ] {
            c.add_code(gt);
        }
        assert_eq!((c.n_ref_alleles, c.n_alt_alleles), (3, 3));
        assert_eq!((c.n_het, c.n_hom_alt, c.n_hom_ref), (1, 1, 1));
        assert_eq!(c.n_full_called(), 3);
        assert_eq!(c.n_missing_alleles(4), 2);
        assert_eq!(c.n_called_alleles(4), 6);
        assert_eq!(c.n_half_called, 0);
    }

    #[test]
    fn add_code_half_call_counts_the_present_allele() {
        // s0=1/. (half call) → alt=1, ref=0, half=1.
        let mut c = SiteGroupCounts::default();
        c.add_code(code(STATE_ALT, STATE_MISSING));
        assert_eq!((c.n_ref_alleles, c.n_alt_alleles), (0, 1));
        assert_eq!(c.n_half_called, 1);
        assert_eq!(c.n_full_called(), 0);
        assert_eq!(c.n_missing_alleles(1), 1);
    }

    #[test]
    fn code_has_alt_sees_either_slot() {
        assert!(code_has_alt(code(STATE_ALT, STATE_MISSING)));
        assert!(code_has_alt(code(STATE_REF, STATE_ALT)));
        assert!(!code_has_alt(code(STATE_REF, STATE_MISSING)));
        assert!(!code_has_alt(code(STATE_MISSING, STATE_MISSING)));
    }

    #[test]
    fn all_ref_site_counts() {
        let mut c = SiteGroupCounts::default();
        for _ in 0..4 {
            c.add_code(code(STATE_REF, STATE_REF));
        }
        assert_eq!((c.n_ref_alleles, c.n_alt_alleles), (8, 0));
        assert_eq!((c.n_hom_ref, c.n_het, c.n_hom_alt), (4, 0, 0));
        assert_eq!(c.n_half_called, 0);
    }
}
