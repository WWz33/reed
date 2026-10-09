//! Site identity and classification (`docs/spec/architecture.md` §5.7).

use bitflags::bitflags;

use crate::coord::Pos0;
use crate::id::ContigId;

/// What kind of variation a record carries, per ALT allele.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum VariantClass {
    /// Single-nucleotide substitution (equal length 1).
    Snv,
    /// Multi-nucleotide substitution (equal length > 1; no per-position
    /// analysis is performed, so compound SNVs are classified here too).
    Mnv,
    /// Insertion or deletion (unequal lengths).
    Indel,
    /// Symbolic allele (`<...>`) or breakend notation.
    Symbolic,
    /// No ALT: the record is monomorphic (ALT `.` or equal to REF).
    #[default]
    Ref,
    /// Anything else (e.g. the spanning-deletion `*` allele).
    Other,
}

impl VariantClass {
    /// Stable lower-case token used in CLI lists and reports.
    #[must_use]
    pub const fn as_str(&self) -> &'static str {
        match self {
            Self::Snv => "snv",
            Self::Mnv => "mnv",
            Self::Indel => "indel",
            Self::Symbolic => "symbolic",
            Self::Ref => "ref",
            Self::Other => "other",
        }
    }
}

/// Classifies one REF/ALT allele pair.
///
/// `alt_allele` must be a single allele (commas already split by the caller).
#[must_use]
pub fn classify_variant(reference: &[u8], alt_allele: &[u8]) -> VariantClass {
    if alt_allele == b"." || alt_allele == reference {
        return VariantClass::Ref;
    }
    if alt_allele.starts_with(b"<")
        || alt_allele.ends_with(b">")
        || alt_allele.contains(&b'[')
        || alt_allele.contains(&b']')
    {
        return VariantClass::Symbolic;
    }
    if alt_allele == b"*" {
        return VariantClass::Other;
    }
    if reference.len() == alt_allele.len() {
        if reference.len() == 1 {
            VariantClass::Snv
        } else {
            VariantClass::Mnv
        }
    } else {
        VariantClass::Indel
    }
}

bitflags! {
    /// Per-site boolean attributes.
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
    pub struct SiteFlags: u8 {
        /// FILTER column is `PASS` (or absent, matching VCF semantics).
        const IS_PASS = 1;
        /// Segregating within the analysed sample set (any called ALT allele).
        const IS_SEGREGATING = 1 << 1;
    }
}

/// Identity of one site: where it is, what it is, and sample-set flags.
///
/// Position is 0-based and stored as `u64` (large-genome safe).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SiteKey {
    /// Contig in header order.
    pub contig: ContigId,
    /// 0-based position.
    pub pos: Pos0,
    /// Classification of the (primary) ALT allele.
    pub class: VariantClass,
    /// Number of ALT alleles before any policy collapses or splits them.
    pub n_alt: u8,
    /// FILTER / segregating flags.
    pub flags: SiteFlags,
}

impl SiteKey {
    /// Whether the record passes FILTER.
    #[must_use]
    pub const fn is_pass(&self) -> bool {
        self.flags.contains(SiteFlags::IS_PASS)
    }

    /// Whether the site segregates within the analysed samples.
    #[must_use]
    pub const fn is_segregating(&self) -> bool {
        self.flags.contains(SiteFlags::IS_SEGREGATING)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classification_matches_vcf_conventions() {
        assert_eq!(classify_variant(b"A", b"G"), VariantClass::Snv);
        assert_eq!(classify_variant(b"AG", b"CT"), VariantClass::Mnv);
        assert_eq!(classify_variant(b"A", b"AGG"), VariantClass::Indel);
        assert_eq!(classify_variant(b"AGG", b"A"), VariantClass::Indel);
        assert_eq!(classify_variant(b"A", b"."), VariantClass::Ref);
        assert_eq!(classify_variant(b"A", b"A"), VariantClass::Ref);
        assert_eq!(classify_variant(b"A", b"<DEL>"), VariantClass::Symbolic);
        assert_eq!(classify_variant(b"A", b"*"), VariantClass::Other);
        assert_eq!(classify_variant(b"A", b"A]chr:1]"), VariantClass::Symbolic);
    }

    #[test]
    fn class_tokens() {
        assert_eq!(VariantClass::Snv.as_str(), "snv");
        assert_eq!(VariantClass::default(), VariantClass::Ref);
    }
}
