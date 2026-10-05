//! The `summary` section (§10.1): record and variant-class counts, the
//! bcftools-stats `SN` analogue.

use reed_core::site::VariantClass;
use serde::{Deserialize, Serialize};

/// Additive counts of the summary section.
///
/// Merge is component-wise addition: exact integers, order-independent, so
/// per-chunk partials reduce deterministically (§8.4).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct SummaryAcc {
    /// Total records seen.
    pub n_records: u64,
    /// Records whose primary ALT is a single-nucleotide substitution.
    pub n_snv: u64,
    /// Records whose primary ALT is a multi-nucleotide substitution.
    pub n_mnv: u64,
    /// Records whose primary ALT is an insertion or deletion.
    pub n_indel: u64,
    /// Records whose primary ALT is symbolic (`<...>`, breakends).
    pub n_symbolic: u64,
    /// Records with no ALT (`.` or REF-equal): monomorphic sites.
    pub n_ref: u64,
    /// Records with any other primary ALT class (e.g. spanning `*`).
    pub n_other: u64,
    /// Records carrying more than one ALT allele.
    pub n_multiallelic: u64,
}

impl SummaryAcc {
    /// Absorbs one record.
    ///
    /// `n_alt` is the number of ALT alleles before any policy (`.` → 0).
    pub fn push(&mut self, class: VariantClass, n_alt: u8) {
        self.n_records += 1;
        match class {
            VariantClass::Snv => self.n_snv += 1,
            VariantClass::Mnv => self.n_mnv += 1,
            VariantClass::Indel => self.n_indel += 1,
            VariantClass::Symbolic => self.n_symbolic += 1,
            VariantClass::Ref => self.n_ref += 1,
            VariantClass::Other => self.n_other += 1,
            _ => self.n_other += 1,
        }
        if n_alt > 1 {
            self.n_multiallelic += 1;
        }
    }

    /// Merges another partial (component-wise addition).
    pub fn merge(&mut self, other: &Self) {
        self.n_records += other.n_records;
        self.n_snv += other.n_snv;
        self.n_mnv += other.n_mnv;
        self.n_indel += other.n_indel;
        self.n_symbolic += other.n_symbolic;
        self.n_ref += other.n_ref;
        self.n_other += other.n_other;
        self.n_multiallelic += other.n_multiallelic;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn class_routing_and_multiallelic() {
        let mut acc = SummaryAcc::default();
        acc.push(VariantClass::Snv, 1);
        acc.push(VariantClass::Indel, 2);
        acc.push(VariantClass::Ref, 0);
        assert_eq!(
            acc,
            SummaryAcc {
                n_records: 3,
                n_snv: 1,
                n_indel: 1,
                n_ref: 1,
                n_multiallelic: 1,
                ..SummaryAcc::default()
            }
        );
    }

    #[test]
    fn merge_is_addition() {
        let mut a = SummaryAcc { n_records: 2, n_snv: 2, ..SummaryAcc::default() };
        let b = SummaryAcc { n_records: 3, n_indel: 3, ..SummaryAcc::default() };
        a.merge(&b);
        assert_eq!(a.n_records, 5);
        assert_eq!(a.n_snv, 2);
        assert_eq!(a.n_indel, 3);
    }
}
