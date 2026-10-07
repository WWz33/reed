//! Site-level record view (`doc_v2/03-architecture.md` §3.2): CHROM resolved
//! through a contig table with a last-hit cache, POS, REF/ALT, QUAL, FILTER
//! and zero-allocation INFO lookups.
//!
//! Consumers that address the index (the filter's tid per record) resolve the
//! contig id once here instead of scanning the contig table per record.

use std::collections::HashMap;

use crate::error::IoError;
use crate::vcf::fixed::parse_pos1;
use crate::vcf::info::info_fragments;
use crate::vcf::line::LineRef;
use crate::Locus;

/// Contig-name to `tid` table with a last-hit cache (records from one contig
/// arrive in runs).
#[derive(Debug, Default)]
pub struct ContigTable {
    ids: HashMap<Box<[u8]>, u32>,
    last: (Vec<u8>, u32),
}

impl ContigTable {
    /// Builds the table from the header's `##contig` order (the tid space).
    #[must_use]
    pub fn new(names: impl IntoIterator<Item = Vec<u8>>) -> Self {
        let mut ids = HashMap::new();
        for (i, name) in names.into_iter().enumerate() {
            if let Ok(id) = u32::try_from(i) {
                ids.insert(name.into_boxed_slice(), id);
            }
        }
        Self { ids, last: (Vec::new(), u32::MAX) }
    }

    /// The contig's id, or `None` when the header does not declare it.
    pub fn id_of(&mut self, chrom: &[u8]) -> Option<u32> {
        if self.last.0 == chrom {
            return Some(self.last.1);
        }
        let id = self.ids.get(chrom).copied()?;
        self.last = (chrom.to_vec(), id);
        Some(id)
    }
}

/// One record's site-level columns.
#[derive(Clone, Copy, Debug)]
pub struct SiteView<'a> {
    /// The contig's tid, resolved when the record arrived through
    /// [`SiteView::parse_declared`] (or left unresolved).
    pub contig: Option<u32>,
    /// CHROM bytes.
    pub contig_name: &'a [u8],
    /// POS (1-based).
    pub pos1: u64,
    /// REF column.
    pub reference: &'a [u8],
    /// ALT column (comma-separated, unparsed).
    pub alts: &'a [u8],
    /// QUAL, `None` when `.`.
    pub qual: Option<f32>,
    /// FILTER column.
    pub filter: &'a [u8],
    /// INFO column.
    pub info: &'a [u8],
}

impl<'a> SiteView<'a> {
    /// Parses the site-level columns without resolving the contig id.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for an unparsable POS or QUAL.
    pub fn parse(line: &LineRef<'a>) -> Result<Self, IoError> {
        let locus = Locus {
            contig_name: String::from_utf8_lossy(line.chrom()).into_owned(),
            pos1: 0,
        };
        let pos1 = parse_pos1(line.pos_text(), &locus)?;
        let locus = Locus { pos1, ..locus };
        let qual = crate::vcf::fixed::parse_qual(line.qual_text(), &locus)?;
        Ok(Self {
            contig: None,
            contig_name: line.chrom(),
            pos1,
            reference: line.reference(),
            alts: line.alts(),
            qual,
            filter: line.filter_text(),
            info: line.info_text(),
        })
    }

    /// Parses the site-level columns and resolves the contig id through
    /// `contigs`.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for an unparsable POS/QUAL and for a
    /// CHROM the header does not declare.
    pub fn parse_declared(
        line: &LineRef<'a>,
        contigs: &mut ContigTable,
    ) -> Result<Self, IoError> {
        let mut site = Self::parse(line)?;
        let Some(contig) = contigs.id_of(site.contig_name) else {
            return Err(IoError::MalformedRecord {
                locus: site.locus(),
                reason: "contig is not declared in the header",
            });
        };
        site.contig = Some(contig);
        Ok(site)
    }

    /// The record's locus, for error messages.
    #[must_use]
    pub fn locus(&self) -> Locus {
        Locus {
            contig_name: String::from_utf8_lossy(self.contig_name).into_owned(),
            pos1: self.pos1,
        }
    }

    /// Number of ALT alleles (`ALT=.` = 0).
    #[must_use]
    pub fn n_alt(&self) -> u32 {
        if self.alts == b"." {
            return 0;
        }
        let commas = self.alts.iter().filter(|&&b| b == b',').count();
        u32::try_from(commas + 1).unwrap_or(u32::MAX)
    }

    /// Whether ALT is `.` (pixy's all-sites marker).
    #[must_use]
    pub fn is_invariant(&self) -> bool {
        self.alts == b"."
    }

    /// The `i`-th ALT allele text.
    #[must_use]
    pub fn alt_at(&self, i: usize) -> Option<&'a [u8]> {
        if self.is_invariant() {
            return None;
        }
        self.alts.split(|&b| b == b',').nth(i)
    }

    /// Zero-allocation INFO lookup: the fragment's value bytes (`None` for a
    /// flag or an absent tag).
    #[must_use]
    pub fn info_value(&self, tag: &[u8]) -> Option<&'a [u8]> {
        info_fragments(self.info)
            .find(|(key, _)| *key == tag)
            .and_then(|(_, value)| value)
    }

    /// INFO tag as an integer (`None` when absent or unparsable).
    #[must_use]
    pub fn info_int(&self, tag: &[u8]) -> Option<i64> {
        let text = self.info_value(tag)?;
        std::str::from_utf8(text).ok()?.trim().parse().ok()
    }

    /// INFO tag as a float (`None` when absent or unparsable).
    #[must_use]
    pub fn info_float(&self, tag: &[u8]) -> Option<f64> {
        let text = self.info_value(tag)?;
        std::str::from_utf8(text).ok()?.trim().parse().ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(text: &[u8]) -> LineRef<'_> {
        LineRef::parse(text).expect("invariant: 10 columns")
    }

    #[test]
    fn parses_site_columns_and_alts() {
        let lr = line(b"chr1\t10\t.\tA\tG,T\t50\tPASS\tDP=5;AF=0.25\tGT\t0/1\n");
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        assert_eq!(site.contig_name, b"chr1");
        assert_eq!(site.pos1, 10);
        assert_eq!(site.reference, b"A");
        assert_eq!(site.n_alt(), 2);
        assert_eq!(site.alt_at(1), Some(b"T".as_slice()));
        assert_eq!(site.qual, Some(50.0));
        assert_eq!(site.filter, b"PASS");
        assert!(!site.is_invariant());
        assert_eq!(site.info_int(b"DP"), Some(5));
        assert_eq!(site.info_float(b"AF"), Some(0.25));
        assert_eq!(site.info_value(b"AF"), Some(b"0.25".as_slice()));
        assert_eq!(site.info_value(b"ABSENT"), None);
    }

    #[test]
    fn invariant_alt_and_missing_qual() {
        let lr = line(b"chr2\t7\t.\tC\t.\t.\t.\t.\tGT\t./.\n");
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        assert_eq!(site.n_alt(), 0);
        assert!(site.is_invariant());
        assert_eq!(site.alt_at(0), None);
        assert_eq!(site.qual, None);
    }

    #[test]
    fn contig_ids_follow_header_order_and_errors() {
        let mut contigs = ContigTable::new(vec![b"chr1".to_vec(), b"chr2".to_vec()]);
        assert_eq!(contigs.id_of(b"chr1"), Some(0));
        assert_eq!(contigs.id_of(b"chr2"), Some(1));
        assert_eq!(contigs.id_of(b"chrX"), None);

        let lr = line(b"chr2\t1\t.\tC\tT\t.\t.\t.\tGT\t0/0\n");
        let site = SiteView::parse_declared(&lr, &mut contigs).expect("invariant: declared");
        assert_eq!(site.contig, Some(1));
        let lr = line(b"chrX\t1\t.\tC\tT\t.\t.\t.\tGT\t0/0\n");
        let err = SiteView::parse_declared(&lr, &mut contigs).unwrap_err();
        assert!(matches!(err, IoError::MalformedRecord { .. }));
    }

    #[test]
    fn bad_pos_is_an_error() {
        let lr = line(b"chr1\tX\t.\tC\tT\t.\t.\t.\tGT\t0/0\n");
        assert!(SiteView::parse(&lr).is_err());
    }
}
