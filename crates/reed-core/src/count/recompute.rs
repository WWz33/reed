//! The `--fill-tags` registry (doc_v30 R13): one static table owns each
//! tag's name, header type, Number, description and compute function.
//!
//! Before this, the same five tags were enumerated in six places (the
//! accepted set, the plan's flags, the fill loop, the `##INFO` renderer, the
//! CLI help and the docs). Adding a tag now means adding one row here; the
//! filter plan, the header writer and the CLI help all derive from it.

// Float tags are f32 by the VCF header type (`InfoTagType::Float`), so the
// f64→f32 casts below are the declared precision, not an accident.
#![allow(clippy::cast_possible_truncation)]

use crate::count::SiteCounts;
use crate::math::hwe::HweCache;
use crate::vcf::line_edit::InfoValueOut;

/// The VCF type of a recomputed INFO tag.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InfoTagType {
    /// Integer (`AC`/`AN`/`NS`).
    Integer,
    /// Float (`MAF`/`HWE`).
    Float,
}

/// The VCF Number of a recomputed INFO tag.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InfoTagNumber {
    /// One value per ALT allele (`A`).
    A,
    /// Exactly one value.
    One,
}

/// One `--fill-tags` tag (doc_v30 R13).
#[derive(Clone, Copy, Debug)]
pub struct RecomputeTag {
    /// Canonical spelling (bcftools matches case-insensitively).
    pub name: &'static str,
    /// Header `Type`.
    pub ty: InfoTagType,
    /// Header `Number`.
    pub number: InfoTagNumber,
    /// Header `Description`.
    pub description: &'static str,
    /// The value from the post-mask counts.
    pub compute: fn(&SiteCounts, &mut HweCache) -> InfoValueOut,
}

/// A tag's identity is its canonical name: two rows with the same name are
/// the same tag (the registry guarantees names are unique). Comparing the
/// compute function pointers instead would tie plan equality to codegen.
impl PartialEq for RecomputeTag {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

impl Eq for RecomputeTag {}

impl RecomputeTag {
    /// Finds a tag by name, case-insensitively.
    #[must_use]
    pub fn find(name: &str) -> Option<&'static Self> {
        RECOMPUTE_TAGS
            .iter()
            .find(|tag| tag.name.eq_ignore_ascii_case(name.trim()))
    }

    /// The canonical name, when `name` is a supported tag.
    #[must_use]
    pub fn canonical(name: &str) -> Option<&'static str> {
        Self::find(name).map(|tag| tag.name)
    }
}

/// The implemented `--fill-tags` set, in fill order.
///
/// AF and F_MISSING are in bcftools' set and in `docs/spec/cli.md`'s prose
/// but are not implemented here (doc_v12 N10); adding them is one row.
pub const RECOMPUTE_TAGS: [RecomputeTag; 5] = [
    RecomputeTag {
        name: "AC",
        ty: InfoTagType::Integer,
        number: InfoTagNumber::A,
        description: "Allele count in genotypes",
        compute: |counts, _| InfoValueOut::Int(counts.ac as i64),
    },
    RecomputeTag {
        name: "AN",
        ty: InfoTagType::Integer,
        number: InfoTagNumber::One,
        description: "Total number of alleles in called genotypes",
        compute: |counts, _| InfoValueOut::Int(counts.an as i64),
    },
    RecomputeTag {
        name: "NS",
        ty: InfoTagType::Integer,
        number: InfoTagNumber::One,
        description: "Number of samples with data",
        compute: |counts, _| InfoValueOut::Int(counts.ns as i64),
    },
    RecomputeTag {
        name: "MAF",
        ty: InfoTagType::Float,
        number: InfoTagNumber::One,
        description: "Frequency of the second most common allele",
        compute: |counts, _| InfoValueOut::Float(counts.maf() as f32),
    },
    RecomputeTag {
        name: "HWE",
        ty: InfoTagType::Float,
        number: InfoTagNumber::A,
        description: "HWE test (PMID:15789306); 1=good, 0=bad",
        compute: |counts, hwe| {
            let p = counts
                .hwe_test_input()
                .map_or(1.0, |(n, mac, het)| hwe.p_exact(n, mac, het));
            InfoValueOut::Float(p as f32)
        },
    },
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_registry_names_and_types_are_unique_and_canonical() {
        for (i, tag) in RECOMPUTE_TAGS.iter().enumerate() {
            assert_ne!(tag.name, "");
            assert_eq!(tag.name, tag.name.to_ascii_uppercase());
            assert!(
                RECOMPUTE_TAGS[..i].iter().all(|t| t.name != tag.name),
                "duplicate tag {}",
                tag.name
            );
        }
        assert_eq!(RecomputeTag::canonical("ac"), Some("AC"));
        assert_eq!(RecomputeTag::canonical("Hwe"), Some("HWE"));
        assert_eq!(RecomputeTag::canonical("af"), None);
    }
}
