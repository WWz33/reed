//! Region rules (§11.1 `[targets]`): bcftools `--targets` streaming
//! semantics — a site is kept when it falls in any include region and in
//! no exclude region. Never needs an index.

use serde::{Deserialize, Serialize};

use crate::FilterError;

/// Include/exclude region lists.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct RegionRule {
    /// Empty = keep everything.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub include: Vec<String>,
    /// A site in any exclude region is dropped.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exclude: Vec<String>,
}

impl RegionRule {
    /// Whether the rule keeps everything (both lists empty).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.include.is_empty() && self.exclude.is_empty()
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

/// Parses a list of region spellings into triples.
///
/// # Errors
///
/// [`FilterError::Spec`] for a malformed region.
pub fn parse_regions(list: &[String]) -> Result<Vec<RegionTriple>, FilterError> {
    list.iter()
        .map(|s| {
            let r = Region::parse(s)?;
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
    /// Parses `Chr01` / `Chr01:100` / `Chr01:100-200` / `Chr01:100-` /
    /// `Chr01:-200`.
    ///
    /// A bare `Chr01:100` is a single position (`start == end`), matching
    /// bcftools `-t/-r` (doc_v21 M5); such positions are what a two-column
    /// `-t`/`-T` file carries.
    pub fn parse(s: &str) -> Result<Self, FilterError> {
        let s = s.trim();
        let (contig, interval) = match s.split_once(':') {
            Some((c, iv)) => (c, Some(iv)),
            None => (s, None),
        };
        if contig.is_empty() {
            return Err(FilterError::Spec(format!("empty contig in region {s:?}")));
        }
        let (start, end) = match interval {
            None => (None, None),
            Some(iv) => match iv.split_once('-') {
                Some((a, b)) => {
                    let p = |x: &str, what: &str| -> Result<Option<u64>, FilterError> {
                        if x.is_empty() {
                            return Ok(None);
                        }
                        x.parse::<u64>().map(Some).map_err(|_| {
                            FilterError::Spec(format!("region {what} is not a number: {s:?}"))
                        })
                    };
                    (p(a, "start")?, p(b, "end")?)
                }
                // A single 1-based position: `start == end`.
                None => {
                    let pos = iv.parse::<u64>().map_err(|_| {
                        FilterError::Spec(format!("region position is not a number: {s:?}"))
                    })?;
                    (Some(pos), Some(pos))
                }
            },
        };
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
}
