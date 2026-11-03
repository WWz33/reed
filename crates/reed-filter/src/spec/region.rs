//! Region rules (§11.1 `[targets]`): bcftools `--targets` streaming
//! semantics — a site is kept when it falls in any include region and in
//! no exclude region. Never needs an index.
//!
//! `-r/-R` (indexed regions) land in [`RegionRule::indexed`], a **separate**
//! list that is ANDed with the `-t/-T` target rule: bcftools intersects the
//! two families (`-r chr2 -t chr1:10` yields nothing, `-r chr1,chr2 -t
//! chr2:10` yields chr2 only), and so does reed (doc_v23 V23-61). Within one
//! family the CLI refuses two sources instead of copying bcftools'
//! order-dependent overwrite (`vcfview.c:651-652` writes both flags into one
//! slot, so `-r`+`-R` silently keeps the last one *and* leaves `is_file`
//! set from either).

use serde::{Deserialize, Serialize};

use crate::FilterError;

/// Include/exclude region lists.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct RegionRule {
    /// `-r/-R`: the indexed region family. ANDed with [`Self::include`] when
    /// both are present (bcftools' `-r`+`-t` intersection).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub indexed: Vec<String>,
    /// `-t/-T`: empty = keep everything.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub include: Vec<String>,
    /// A site in any exclude region is dropped (`^` entries).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exclude: Vec<String>,
}

impl RegionRule {
    /// Whether the rule keeps everything (all lists empty).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.indexed.is_empty() && self.include.is_empty() && self.exclude.is_empty()
    }
}

/// A parsed genomic interval (1-based closed; open ends allowed).
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    /// Contig name (exact match).
    pub contig: String,
    /// Inclusive 1-based lower bound (`None` = open).
    pub start: Option<u64>,
    /// Inclusive 1-based upper bound (`None` = open).
    pub end: Option<u64>,
}

/// A parsed region triple: `(contig, start1, end1)` with open ends as
/// `None` (the shape the filter's head rule stores).
pub type RegionTriple = (String, Option<u64>, Option<u64>);

/// Parses a list of region spellings into triples, using the header's
/// contig vocabulary to resolve names that contain `:` (doc_v24 V24-04).
///
/// # Errors
///
/// [`FilterError::Spec`] for a malformed region.
pub fn parse_regions_with(
    list: &[String],
    contigs: &[String],
) -> Result<Vec<RegionTriple>, FilterError> {
    list.iter()
        .map(|s| {
            let r = Region::parse_in(s, contigs)?;
            Ok((r.contig, r.start, r.end))
        })
        .collect()
}

/// Header-less variant of [`parse_regions_with`] (tests and preview paths).
/// Production callers must pass the header's contigs.
///
/// # Errors
///
/// [`FilterError::Spec`] for a malformed region.
pub fn parse_regions(list: &[String]) -> Result<Vec<RegionTriple>, FilterError> {
    parse_regions_with(list, &[])
}

/// Like [`parse_regions_with`], but an indexed region whose contig is not in
/// the header is an error. `-r/-R` must name a contig the index can seek to;
/// `-t/-T` stays lenient (a target list may name contigs absent from one
/// file, and bcftools simply selects nothing there).
///
/// # Errors
///
/// [`FilterError::Spec`] for a malformed region or an undeclared contig.
pub fn parse_regions_strict(
    list: &[String],
    contigs: &[String],
) -> Result<Vec<RegionTriple>, FilterError> {
    list.iter()
        .map(|s| {
            let r = Region::parse_in(s, contigs)?;
            if !contigs.is_empty() && !contigs.iter().any(|c| c == &r.contig) {
                return Err(FilterError::Spec(format!(
                    "region {s:?}: contig {:?} is not declared in the header",
                    r.contig
                )));
            }
            Ok((r.contig, r.start, r.end))
        })
        .collect()
}

/// bcftools `-t/-T` semantics: kept when inside any include region (or the
/// include list is empty) and in no exclude region.
#[must_use]
pub fn region_allows(
    include: &[RegionTriple],
    exclude: &[RegionTriple],
    chrom: &[u8],
    pos1: u64,
) -> bool {
    let inside = |list: &[RegionTriple]| {
        list.iter().any(|(c, s, e)| {
            c.as_bytes() == chrom && s.is_none_or(|s| pos1 >= s) && e.is_none_or(|e| pos1 <= e)
        })
    };
    (include.is_empty() || inside(include)) && !inside(exclude)
}

impl Region {
    /// Header-less heuristic parse: `Chr01` / `Chr01:100` /
    /// `Chr01:100-200` / `Chr01:100-` / `Chr01:-200`.
    ///
    /// This cannot disambiguate a contig whose name ends in `:<number>`; the
    /// production entry point is [`Self::parse_in`] with the header's contig
    /// list (doc_v24 V24-04). Kept for tests, `convert --region` and other
    /// preview paths that have no header at hand.
    ///
    /// # Errors
    ///
    /// [`FilterError::Spec`] for a malformed region.
    pub fn parse(s: &str) -> Result<Self, FilterError> {
        Self::parse_in(s, &[])
    }

    /// Header-aware parse. A whole string that names a declared contig wins
    /// over splitting (`HLA-A*01:01` is one contig); otherwise the **last**
    /// `:` with a valid interval splits name from range, so
    /// `HLA-A*01:01:100-200` works. A `:` with no valid interval is an
    /// error rather than a silently wrong contig.
    ///
    /// # Errors
    ///
    /// [`FilterError::Spec`] for an empty name or a malformed interval.
    pub fn parse_in(s: &str, contigs: &[String]) -> Result<Self, FilterError> {
        let s = s.trim();
        if s.is_empty() {
            return Err(FilterError::Spec("empty region".into()));
        }
        if contigs.iter().any(|c| c == s) {
            return Ok(Self {
                contig: s.to_string(),
                start: None,
                end: None,
            });
        }
        let Some(pos) = s.rfind(':') else {
            return Ok(Self {
                contig: s.to_string(),
                start: None,
                end: None,
            });
        };
        let (contig, interval) = (&s[..pos], &s[pos + 1..]);
        if contig.is_empty() {
            return Err(FilterError::Spec(format!("empty contig in region {s:?}")));
        }
        if !valid_interval(interval) {
            return Err(FilterError::Spec(format!(
                "region {s:?} is neither a declared contig nor <contig>:<interval>"
            )));
        }
        let (start, end) = parse_interval(interval, s)?;
        Ok(Self {
            contig: contig.to_string(),
            start,
            end,
        })
    }

    /// 1-based closed-interval containment test.
    #[must_use]
    pub fn contains(&self, contig: &str, pos: u64) -> bool {
        self.contig == contig
            && self.start.is_none_or(|s| pos >= s)
            && self.end.is_none_or(|e| pos <= e)
    }
}

/// `100`, `100-200`, `100-`, `-200` (at most one dash, at least one digit).
fn valid_interval(text: &str) -> bool {
    let mut dashes = 0usize;
    let mut digits = 0usize;
    for c in text.chars() {
        match c {
            '-' => {
                dashes += 1;
                if dashes > 1 {
                    return false;
                }
            }
            c if c.is_ascii_digit() => digits += 1,
            _ => return false,
        }
    }
    digits > 0
}

fn parse_interval(interval: &str, whole: &str) -> Result<(Option<u64>, Option<u64>), FilterError> {
    match interval.split_once('-') {
        Some((a, b)) => {
            let p = |x: &str, what: &str| -> Result<Option<u64>, FilterError> {
                if x.is_empty() {
                    return Ok(None);
                }
                x.parse::<u64>()
                    .map(Some)
                    .map_err(|_| FilterError::Spec(format!("region {what} is not a number: {whole:?}")))
            };
            Ok((p(a, "start")?, p(b, "end")?))
        }
        None => {
            let pos = interval.parse::<u64>().map_err(|_| {
                FilterError::Spec(format!("region position is not a number: {whole:?}"))
            })?;
            Ok((Some(pos), Some(pos)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn m5_bare_position_is_a_single_site() {
        // doc_v21 M5: bcftools accepts `-t CHROM:POS`; reed rejected it, so
        // every two-column target file and every single-position script failed.
        let r = Region::parse("chr1:200").expect("bare position");
        assert_eq!((r.contig.as_str(), r.start, r.end), ("chr1", Some(200), Some(200)));
        assert!(r.contains("chr1", 200));
        assert!(!r.contains("chr1", 199));
        assert!(!r.contains("chr1", 201));
        // Ranges and open ends keep their meaning.
        let r = Region::parse("chr1:100-200").expect("range");
        assert_eq!((r.start, r.end), (Some(100), Some(200)));
        let r = Region::parse("chr1:100-").expect("open end");
        assert_eq!((r.start, r.end), (Some(100), None));
        let r = Region::parse("chr1:-200").expect("open start");
        assert_eq!((r.start, r.end), (None, Some(200)));
        let r = Region::parse("chr1").expect("contig");
        assert_eq!((r.start, r.end), (None, None));
    }

    #[test]
    fn m5_non_numeric_position_is_rejected() {
        assert!(Region::parse("chr1:abc").is_err());
        assert!(Region::parse("chr1:10-x").is_err());
        assert!(Region::parse(":100").is_err());
    }

    /// doc_v24 V24-04: a declared contig wins over the `:` split, and a
    /// contig that contains `:` still accepts an interval.
    #[test]
    fn declared_contigs_disambiguate_colons() {
        let contigs = vec!["HLA-A*01:01".to_string(), "chr1:100".to_string()];
        let whole = Region::parse_in("HLA-A*01:01", &contigs).expect("whole contig");
        assert_eq!((whole.contig.as_str(), whole.start, whole.end), ("HLA-A*01:01", None, None));
        let ranged = Region::parse_in("HLA-A*01:01:100-200", &contigs).expect("interval");
        assert_eq!(
            (ranged.contig.as_str(), ranged.start, ranged.end),
            ("HLA-A*01:01", Some(100), Some(200))
        );
        // A colon-containing contig also wins as a whole name.
        let whole = Region::parse_in("chr1:100", &contigs).expect("whole contig");
        assert_eq!((whole.contig.as_str(), whole.start, whole.end), ("chr1:100", None, None));
        // Without the vocabulary the same text is an interval.
        let heur = Region::parse("chr1:100").expect("heuristic");
        assert_eq!((heur.contig.as_str(), heur.start, heur.end), ("chr1", Some(100), Some(100)));
    }

    /// A colon with no valid interval is loud, not a silent whole-contig
    /// match: `HLA-A*01:abc` is a typo, not a contig the header declared.
    #[test]
    fn an_undeclared_colon_name_is_an_error() {
        let contigs = vec!["HLA-A*01:01".to_string()];
        let err = Region::parse_in("HLA-A*01:bc", &contigs).expect_err("malformed");
        assert!(err.to_string().contains("neither a declared contig"), "{err}");
    }
}
