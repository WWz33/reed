//! The IUPAC nucleotide ambiguity table and the binary-NEXUS genotype table
//! (doc_v13 §3.2–§3.4): the single implementation of "genotype allele set →
//! matrix letter" shared by every convert target (D13-7).
//!
//! The tables are factual data from the IUPAC nucleotide ambiguity codes
//! (`Nomenclature for incompletely specified bases in nucleic acid
//! sequences`, 1985) plus the `*` (GATK deletion) participations the
//! reference implementation's table carries. Lowercase entries mark a set
//! that includes `N` or `*` (soft-masked consensus), matching vcf2phylip's
//! table exactly.

/// One allele-set → letter entry.
pub type AmbigEntry = (&'static [u8], u8);

/// The 62-entry ambiguity table, sorted by key. Keys are the *sorted,
/// de-duplicated* concatenation of the allele letters the genotype uses;
/// a key that is absent means the allele set is not a valid nucleotide
/// combination (the reference implementation raises here).
pub const AMBIG: &[AmbigEntry] = &[
    (b"*", b'-'),
    (b"*A", b'a'),
    (b"*AC", b'm'),
    (b"*ACG", b'v'),
    (b"*ACGN", b'v'),
    (b"*ACGNT", b'N'),
    (b"*ACGT", b'N'),
    (b"*ACN", b'm'),
    (b"*ACNT", b'h'),
    (b"*ACT", b'h'),
    (b"*AG", b'r'),
    (b"*AGN", b'r'),
    (b"*AGNT", b'd'),
    (b"*AGT", b'd'),
    (b"*AN", b'a'),
    (b"*ANT", b'w'),
    (b"*AT", b'w'),
    (b"*C", b'c'),
    (b"*CG", b's'),
    (b"*CGN", b's'),
    (b"*CGNT", b'b'),
    (b"*CGT", b'b'),
    (b"*CN", b'c'),
    (b"*CNT", b'y'),
    (b"*CT", b'y'),
    (b"*G", b'g'),
    (b"*GN", b'g'),
    (b"*GNT", b'k'),
    (b"*GT", b'k'),
    (b"*N", b'n'),
    (b"*NT", b't'),
    (b"*T", b't'),
    (b"A", b'A'),
    (b"AC", b'M'),
    (b"ACG", b'V'),
    (b"ACGN", b'v'),
    (b"ACGT", b'N'),
    (b"ACN", b'm'),
    (b"ACNT", b'h'),
    (b"ACT", b'H'),
    (b"AG", b'R'),
    (b"AGN", b'r'),
    (b"AGNT", b'd'),
    (b"AGT", b'D'),
    (b"AN", b'a'),
    (b"ANT", b'w'),
    (b"AT", b'W'),
    (b"C", b'C'),
    (b"CG", b'S'),
    (b"CGN", b's'),
    (b"CGNT", b'b'),
    (b"CGT", b'B'),
    (b"CN", b'c'),
    (b"CNT", b'y'),
    (b"CT", b'Y'),
    (b"G", b'G'),
    (b"GN", b'g'),
    (b"GNT", b'k'),
    (b"GT", b'K'),
    (b"N", b'N'),
    (b"NT", b't'),
    (b"T", b'T'),
];

/// The diploid biallelic genotype → SNAPP state table (vcf2phylip `GEN_BIN`):
/// `0` = homozygous REF, `1` = heterozygous, `2` = homozygous ALT, `?` =
/// anything else (the caller's fallback).
pub const GEN_BIN: &[(&[u8], u8)] = &[
    (b"./.", b'?'),
    (b".|.", b'?'),
    (b"0/0", b'0'),
    (b"0/1", b'1'),
    (b"0|0", b'0'),
    (b"0|1", b'1'),
    (b"1/0", b'1'),
    (b"1/1", b'2'),
    (b"1|0", b'1'),
    (b"1|1", b'2'),
];

/// Looks up an allele-set key (already sorted and de-duplicated) in
/// [`AMBIG`].
///
/// `None` means the key is absent — the reference implementation raises a
/// `KeyError` here (it crashes on a REF/ALT letter outside the table); reed
/// treats the row as malformed and skips it with a counter instead
/// (doc_v13 B13-9).
#[must_use]
pub fn ambig_lookup(key: &[u8]) -> Option<u8> {
    AMBIG
        .binary_search_by(|(k, _)| k.cmp(&key))
        .ok()
        .map(|i| AMBIG[i].1)
}

/// Looks up a raw genotype token (e.g. `0/1`, `1|1`) in [`GEN_BIN`];
/// `None` falls back to `?` (vcf2phylip's `else` branch).
#[must_use]
pub fn gen_bin_lookup(gt: &[u8]) -> Option<u8> {
    GEN_BIN
        .binary_search_by(|(k, _)| k.cmp(&gt))
        .ok()
        .map(|i| GEN_BIN[i].1)
}

/// One letter → its allele set (the inverse of [`AMBIG`], for hapmap → VCF,
/// doc_v13 §3.12). `None` marks a missing/unknown cell (`N`, `?`); every
/// other letter expands to the base set it stands for — `-` and `+` are
/// TASSEL's gap/insert alleles (kept as their own single-letter sets), and
/// `0` is TASSEL's `-+` diploid (two alleles).
#[must_use]
pub fn ambig_expand(letter: u8) -> Option<&'static [u8]> {
    Some(match letter.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'-' => b"-",
        b'+' => b"+",
        b'0' => b"-+",
        // N and ? are the missing codes (checked first by
        // `is_hapmap_missing`); they expand to no allele. Every other letter
        // outside the IUPAC set is unknown to the caller.
        _ => return None,
    })
}

/// Whether a hapmap cell reads as missing data (doc_v13 §3.12: only `N` and
/// `?` are the missing codes; every other letter is an allele).
#[must_use]
pub fn is_hapmap_missing(cell: u8) -> bool {
    matches!(cell, b'?') || cell == b'N'
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The table is sorted (the binary searches depend on it) and carries
    /// the exact 62 reference entries (doc_v13 01 §6 guard 4).
    #[test]
    fn tables_are_sorted_and_complete() {
        assert_eq!(AMBIG.len(), 62);
        assert_eq!(GEN_BIN.len(), 10);
        assert!(AMBIG.windows(2).all(|w| w[0].0 < w[1].0));
        assert!(GEN_BIN.windows(2).all(|w| w[0].0 < w[1].0));
    }

    /// Spot checks of the IUPAC letters and the `*` participations.
    #[test]
    fn ambiguity_letters_match_iupac() {
        for (key, want) in [
            (&b"A"[..], b'A'),
            (&b"AG"[..], b'R'),
            (&b"CT"[..], b'Y'),
            (&b"ACGT"[..], b'N'),
            (&b"*"[..], b'-'),
            (&b"*ACGT"[..], b'N'),
            (&b"*A"[..], b'a'),
            (&b"AN"[..], b'a'),
        ] {
            assert_eq!(ambig_lookup(key), Some(want), "{key:?}");
        }
        // A letter outside the table is absent, not a fallback.
        assert_eq!(ambig_lookup(b"AR"), None);
    }

    /// `ambig_expand` inverts the IUPAC letters and marks the missing codes.
    #[test]
    fn expand_inverts_iupac() {
        assert_eq!(ambig_expand(b'A'), Some(&b"A"[..]));
        assert_eq!(ambig_expand(b'R'), Some(&b"AG"[..]));
        assert_eq!(ambig_expand(b'a'), Some(&b"A"[..]), "case-insensitive");
        assert_eq!(ambig_expand(b'-'), Some(&b"-"[..]));
        assert_eq!(ambig_expand(b'0'), Some(&b"-+"[..]));
        // The expansion agrees with AMBIG on every pure two-base code.
        for (set, letter) in [
            (&b"A"[..], b'A'),
            (&b"AG"[..], b'R'),
            (&b"CT"[..], b'Y'),
            (&b"CG"[..], b'S'),
            (&b"AT"[..], b'W'),
            (&b"GT"[..], b'K'),
            (&b"AC"[..], b'M'),
        ] {
            assert_eq!(ambig_lookup(set), Some(letter));
            assert_eq!(ambig_expand(letter), Some(set));
        }
        assert_eq!(ambig_expand(b'?'), None);
        assert!(is_hapmap_missing(b'?'));
        assert!(is_hapmap_missing(b'N'));
        assert!(!is_hapmap_missing(b'A'));
        assert!(!is_hapmap_missing(b'-'));
    }

    /// The binary-NEXUS states.
    #[test]
    fn gen_bin_states() {
        assert_eq!(gen_bin_lookup(b"0/0"), Some(b'0'));
        assert_eq!(gen_bin_lookup(b"0|1"), Some(b'1'));
        assert_eq!(gen_bin_lookup(b"1|1"), Some(b'2'));
        assert_eq!(gen_bin_lookup(b"./."), Some(b'?'));
        assert_eq!(gen_bin_lookup(b"2/2"), None);
        assert_eq!(gen_bin_lookup(b"0"), None);
    }
}
