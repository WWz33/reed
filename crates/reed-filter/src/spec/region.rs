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
    /// Parses `Chr01` / `Chr01:100-200` / `Chr01:100-` / `Chr01:-200`.
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
            Some(iv) => {
                let (a, b) = iv
                    .split_once('-')
                    .ok_or_else(|| FilterError::Spec(format!("region is missing '-': {s:?}")))?;
                let p = |x: &str, what: &str| -> Result<Option<u64>, FilterError> {
                    if x.is_empty() {
                        return Ok(None);
                    }
                    x.parse::<u64>()
                        .map(Some)
                        .map_err(|_| FilterError::Spec(format!("region {what} is not a number: {s:?}")))
                };
                (p(a, "start")?, p(b, "end")?)
            }
        };
        Ok(Self { contig: contig.to_string(), start, end })
    }

    /// 1-based closed-interval containment test.
    #[must_use]
    pub fn contains(&self, contig: &str, pos: u64) -> bool {
        self.contig == contig
            && self.start.is_none_or(|s| pos >= s)
            && self.end.is_none_or(|e| pos <= e)
    }
}
