//! The BCF dictionaries (hts-specs §6.2.1): FILTER/INFO/FORMAT share one
//! string dictionary in header order, contigs have their own in `##contig`
//! order. An explicit `IDX=` pins an entry's index; otherwise the next free
//! index is used — the same rule htslib applies when it reads a BCF header.
//!
//! The string dictionary is a sparse map keyed by index: a header `IDX=`
//! is attacker-controlled input, and a dense `Vec` sized by it turns a few
//! header bytes into an allocation of the attacker's choosing. htslib
//! bounds the same input with `idx >= INT_MAX - 1 → skip` and a fallible
//! resize (`vcf.c:800-830`, `vcf.c:877-891`); the map removes the
//! allocation entirely (doc_v10 D10-3).

use std::collections::BTreeMap;

use crate::error::IoError;

/// Index → name lookups for one BCF header.
#[derive(Clone, Debug, Default)]
pub struct BcfDict {
    /// String dictionary: index → tag name (indices are sparse).
    strings: BTreeMap<u32, Vec<u8>>,
    /// Contig names in `##contig` order.
    contigs: Vec<Vec<u8>>,
    /// Contig dictionary index → position in `contigs` (an explicit `IDX=`
    /// may disagree with the line order; doc_v15 D15-24).
    contig_idx: BTreeMap<u32, usize>,
    /// The header declares VCF ≥ 4.4 (`##fileformat`): GT rendering applies
    /// the prefixed-phasing rules (A6-14, htslib `bcf_format_gt_v2`).
    vcf44: bool,
}

/// htslib's IDX ceiling (`idx >= INT_MAX - 1` is skipped, vcf.c:877-891).
const IDX_MAX: u64 = (i32::MAX - 1) as u64;

/// One header entry's dictionary index: its explicit `IDX=` when present
/// and in range, else the next free slot (`next` advances). The same rule
/// serves the tag dictionary and the contig dictionary (doc_v15 D15-24).
fn entry_idx(rest: &[u8], next: &mut u32, what: &str) -> u32 {
    field(rest, b"IDX")
        .and_then(|v| std::str::from_utf8(v).ok())
        .and_then(|v| v.parse::<u64>().ok())
        .filter(|&idx| {
            if idx >= IDX_MAX {
                tracing::warn!(
                    "BCF header {what}IDX={idx} is out of range; assigning the next free index"
                );
                false
            } else {
                true
            }
        })
        .map_or_else(
            || {
                let i = *next;
                *next += 1;
                i
            },
            |idx| u32::try_from(idx).expect("invariant: bounded by IDX_MAX"),
        )
}

impl BcfDict {
    /// Builds the dictionaries from the header text embedded in the BCF.
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when a header line is malformed (no `ID=`).
    pub fn from_header(text: &[u8]) -> Result<Self, IoError> {
        let mut dict = Self::default();
        let mut next = 0_u32;
        let mut next_contig = 0_u32;
        for line in text.split(|&b| b == b'\n') {
            if let Some(version) = line.strip_prefix(b"##fileformat=VCFv") {
                dict.vcf44 = parse_version(version).is_some_and(|v| v >= (4, 4));
                continue;
            }
            if let Some(rest) = line.strip_prefix(b"##contig=<") {
                let name = field(rest, b"ID");
                if let Some(name) = name {
                    // Contig lines honor an explicit `IDX=` exactly like tag
                    // lines (htslib's rule; doc_v15 D15-24).
                    let idx = entry_idx(rest, &mut next_contig, "contig ");
                    dict.contig_idx.insert(idx, dict.contigs.len());
                    dict.contigs.push(name.to_vec());
                    next_contig = next_contig.max(idx + 1);
                }
                continue;
            }
            let kind = if line.starts_with(b"##FILTER=<") {
                Some(b"##FILTER=<".len())
            } else if line.starts_with(b"##INFO=<") {
                Some(b"##INFO=<".len())
            } else if line.starts_with(b"##FORMAT=<") {
                Some(b"##FORMAT=<".len())
            } else {
                None
            };
            let Some(prefix) = kind else { continue };
            let rest = &line[prefix..];
            let Some(name) = field(rest, b"ID") else {
                return Err(IoError::Header(format!(
                    "BCF header line has no ID=: {}",
                    String::from_utf8_lossy(&line[..line.len().min(60)])
                )));
            };
            let idx = entry_idx(rest, &mut next, "");
            dict.strings.insert(idx, name.to_vec());
            next = next.max(idx + 1);
        }
        Ok(dict)
    }

    /// The tag name at `index` (`None` outside the dictionary).
    #[must_use]
    pub fn tag(&self, index: u32) -> Option<&[u8]> {
        self.strings.get(&index).map(Vec::as_slice)
    }

    /// The contig name at `index` (`None` outside the dictionary). The
    /// dictionary index, not the line order: an explicit `IDX=` on a
    /// `##contig` line wins (doc_v15 D15-24).
    #[must_use]
    pub fn contig(&self, index: i32) -> Option<&[u8]> {
        u32::try_from(index)
            .ok()
            .and_then(|i| self.contig_idx.get(&i))
            .and_then(|&pos| self.contigs.get(pos))
            .map(Vec::as_slice)
    }

    /// Number of contigs.
    #[must_use]
    pub fn n_contigs(&self) -> usize {
        self.contigs.len()
    }

    /// Whether the header declares VCF ≥ 4.4 (prefixed GT phasing, A6-14).
    #[must_use]
    pub const fn vcf44(&self) -> bool {
        self.vcf44
    }

    /// The dictionary index of `tag`, when present. The fast path resolves
    /// `GT` once per run and then compares integers (doc_v3 02-bcf.md §2.1).
    #[must_use]
    pub fn key_of(&self, tag: &[u8]) -> Option<u32> {
        self.strings
            .iter()
            .find(|(_, name)| name.as_slice() == tag)
            .map(|(&i, _)| i)
    }
}

/// The `(major, minor)` of a `##fileformat=VCFvX.Y` payload.
fn parse_version(text: &[u8]) -> Option<(u32, u32)> {
    let text = std::str::from_utf8(text).ok()?.trim();
    let (major, minor) = text.split_once('.')?;
    Some((major.parse().ok()?, minor.parse().ok()?))
}

/// The top-level `KEY=value` pairs inside a `##TAG=<...>` payload.
///
/// The payload may contain quoted commas and escaped quotes
/// (`Description="a, b \"c\""`); a `>` inside quotes is data, not the end.
/// One tokenizer serves every header line so the hand-written `field()`
/// lookup can no longer match `ID=` inside a quoted description
/// (doc_v23 V23-35, doc_v24 V24-03). Values are returned verbatim, quotes
/// included; callers only need `ID`/`IDX`, which are never quoted.
pub(crate) struct AnglePairs<'a> {
    rest: &'a [u8],
}

impl<'a> Iterator for AnglePairs<'a> {
    type Item = (&'a [u8], &'a [u8]);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Skip separators between pairs.
            while let Some((&b, tail)) = self.rest.split_first() {
                if b == b',' || b == b' ' {
                    self.rest = tail;
                } else {
                    break;
                }
            }
            if self.rest.is_empty() || self.rest[0] == b'>' {
                return None;
            }
            // Key: up to '='; a bare token (no '=') is skipped.
            let key_end = self
                .rest
                .iter()
                .position(|&b| b == b'=' || b == b',' || b == b'>')?;
            if self.rest[key_end] != b'=' {
                self.rest = &self.rest[key_end..];
                continue;
            }
            let key = &self.rest[..key_end];
            let mut i = key_end + 1;
            let mut quoted = false;
            let mut escaped = false;
            while i < self.rest.len() {
                let b = self.rest[i];
                if escaped {
                    escaped = false;
                    i += 1;
                    continue;
                }
                match b {
                    b'\\' if quoted => {
                        escaped = true;
                        i += 1;
                    }
                    b'"' => {
                        quoted = !quoted;
                        i += 1;
                    }
                    b',' | b'>' if !quoted => break,
                    _ => i += 1,
                }
            }
            let value = &self.rest[key_end + 1..i];
            // Continue after the separator (a '>' leaves the iterator done).
            self.rest = if i < self.rest.len() && self.rest[i] == b',' {
                &self.rest[i + 1..]
            } else {
                &self.rest[i..]
            };
            return Some((key, value));
        }
    }
}

/// [`AnglePairs`] over `rest`.
pub(crate) fn angle_pairs(rest: &[u8]) -> AnglePairs<'_> {
    AnglePairs { rest }
}

/// The `KEY=value` inside a `##TAG=<...>` payload (value up to a top-level
/// `,` or `>`), using the quote-aware tokenizer above.
fn field<'a>(rest: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    angle_pairs(rest)
        .find(|(k, _)| *k == key)
        .map(|(_, value)| value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_the_dictionary_in_header_order() {
        let text = b"##fileformat=VCFv4.2\n\
                     ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
                     ##contig=<ID=chr1,length=1000>\n\
                     ##FILTER=<ID=LowQual,Description=\"Low quality\">\n\
                     ##INFO=<ID=AC,Number=A,Type=Integer,Description=\"a\">\n\
                     ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
                     ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
                     ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q, with comma\">\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";
        let dict = BcfDict::from_header(text).expect("invariant: parses");
        assert_eq!(dict.tag(0).expect("0"), b"PASS");
        assert_eq!(dict.tag(1).expect("1"), b"LowQual");
        assert_eq!(dict.tag(2).expect("2"), b"AC");
        assert_eq!(dict.tag(3).expect("3"), b"DP");
        assert_eq!(dict.tag(4).expect("4"), b"GT");
        assert_eq!(dict.tag(5).expect("5"), b"GQ");
        assert_eq!(dict.contig(0).expect("c0"), b"chr1");
        assert!(dict.tag(6).is_none());
    }

    #[test]
    fn honors_explicit_idx() {
        let text = b"##FILTER=<ID=PASS,Description=\"p\">\n\
                     ##INFO=<ID=X,Number=1,Type=Integer,IDX=7,Description=\"x\">\n\
                     ##INFO=<ID=Y,Number=1,Type=Integer,Description=\"y\">\n";
        let dict = BcfDict::from_header(text).expect("invariant: parses");
        assert_eq!(dict.tag(7).expect("7"), b"X");
        assert_eq!(dict.tag(8).expect("8"), b"Y");
        assert!(dict.tag(1).is_none());
    }

    /// doc_v10 D10-3: a huge `IDX=` is attacker-controlled input. The old
    /// dense `Vec::resize(idx + 1)` turned `IDX=4000000000` into an ~96 GB
    /// allocation (abort on failure). The sparse map plus the htslib
    /// ceiling (`idx >= INT_MAX - 1` skipped, vcf.c:877-891) must keep the
    /// header parseable and tiny.
    #[test]
    fn huge_idx_does_not_allocate_and_falls_back() {
        let text = b"##FILTER=<ID=PASS,Description=\"p\">\n\
                     ##INFO=<ID=X,Number=1,Type=Integer,IDX=4000000000,Description=\"x\">\n\
                     ##INFO=<ID=Y,Number=1,Type=Integer,Description=\"y\">\n";
        let dict = BcfDict::from_header(text).expect("invariant: parses");
        // The out-of-range IDX is skipped; X takes the next free index (1),
        // and Y follows (2).
        assert_eq!(dict.tag(1).expect("1"), b"X");
        assert_eq!(dict.tag(2).expect("2"), b"Y");
        assert!(dict.tag(4_000_000_000).is_none());

        // A large-but-legal IDX stays sparse: no allocation proportional to
        // the index, and the lookup still resolves.
        let text = b"##INFO=<ID=Z,Number=1,Type=Integer,IDX=1000000,Description=\"z\">\n";
        let dict = BcfDict::from_header(text).expect("invariant: parses");
        assert_eq!(dict.tag(1_000_000).expect("1e6"), b"Z");
        assert_eq!(dict.key_of(b"Z"), Some(1_000_000));
    }

    /// doc_v24 V24-03: the tokenizer honors quoted commas, escaped quotes and
    /// a `>` inside a quoted value, and `ID=` inside a description is not a
    /// key.
    #[test]
    fn angle_pairs_respect_quotes_and_escapes() {
        let rest = br#"Description="see ID=FAKE, a \"quote\" and > here",ID=REAL,Number=1"#;
        let pairs: Vec<(&[u8], &[u8])> = angle_pairs(rest).collect();
        assert_eq!(pairs[0].0, b"Description");
        assert!(
            pairs[0].1.starts_with(br#""see ID=FAKE"#),
            "value = {:?}",
            String::from_utf8_lossy(pairs[0].1)
        );
        assert_eq!(pairs[1], (&b"ID"[..], &b"REAL"[..]));
        assert_eq!(pairs[2], (&b"Number"[..], &b"1"[..]));
        assert!(
            field(rest, b"ID") == Some(&b"REAL"[..]),
            "ID= inside a quoted description must not win"
        );
    }

    /// The dictionary itself uses the same tokenizer: a quoted `ID=FAKE`
    /// before the real `ID` must not register a bogus tag.
    #[test]
    fn dict_ignores_id_inside_a_quoted_description() {
        let text = b"##fileformat=VCFv4.4\n\
                     ##INFO=<Description=\"see ID=FAKE\",ID=REAL,Number=1,Type=Integer>\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let dict = BcfDict::from_header(text).expect("invariant: parses");
        assert_eq!(dict.tag(0).expect("0"), b"REAL");
        assert!(dict.key_of(b"FAKE").is_none());
    }
}
