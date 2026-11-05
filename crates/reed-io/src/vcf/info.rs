//! INFO fragment iteration (`docs/spec/architecture.md` §6.4).
//!
//! Typed access by `TagId` lands with the header tag dictionary (M5 filter
//! head predicates); for now this is the zero-copy fragment view.

/// Iterates INFO fragments as `(key, value)`; a Flag yields `value = None`.
/// Empty INFO (`.`) yields nothing.
pub fn info_fragments(info: &[u8]) -> impl Iterator<Item = (&[u8], Option<&[u8]>)> {
    info.split(|&b| b == b';')
        .filter(|f| !f.is_empty() && f != b".")
        .map(|fragment| match memchr::memchr(b'=', fragment) {
            Some(eq) => (&fragment[..eq], Some(&fragment[eq + 1..])),
            None => (fragment, None),
        })
}

/// One INFO value element parsed by the header's declared type.
///
/// The single implementation of "INFO text -> number" (doc_v12 D12-14):
/// the BCF encoder, the expression evaluator and the head rules used three
/// different rules (strict decimal, `str::parse::<f64>`, and a saturated
/// u32), so the same record yielded three different values — `XI=1e3`
/// filtered as 1000, wrote `.` to BCF and kept `1e3` in text.
///
/// The grammar is the strict one, deliberately **not** htslib's prefix
/// parser (`hts_str2int` reads `5abc` as 5, `1e3` as 1, and `0x10` as 0 in
/// Integer but 16 in Float — see 01-review §2.5-B3/B4, a permanent
/// divergence). No trimming: VCFv4.5's numeric grammar has no whitespace
/// (§1.3; values may contain spaces, but that applies to strings, §1.6.1).
///
/// - Integer: optional `-` or `+` (the spec's Float regex allows `[-+]?`,
///   tex:97) then digits; out of `i32` range is the caller's concern (the
///   encoder maps it to missing, doc_v12 N3);
/// - Float: Rust's `f64` grammar, which matches the spec regex (accepts
///   `+5`, `5.0`, `1e3`, `inf`, `nan`; rejects `0x10` and whitespace).
///
/// `None` for `.`, empty or unparsable text.
#[must_use]
pub fn parse_info_number(text: &[u8], float: bool) -> Option<f64> {
    if text.is_empty() || text == b"." {
        return None;
    }
    if float {
        return std::str::from_utf8(text).ok()?.parse::<f64>().ok();
    }
    parse_info_int(text).map(|v| v as f64)
}

/// The integer half of [`parse_info_number`]: strict decimal with an
/// optional sign, `None` on anything else.
#[must_use]
pub fn parse_info_int(text: &[u8]) -> Option<i64> {
    if text.is_empty() || text == b"." {
        return None;
    }
    let (neg, digits) = match text[0] {
        b'-' => (true, &text[1..]),
        b'+' => (false, &text[1..]),
        _ => (false, text),
    };
    if digits.is_empty() {
        return None;
    }
    let mut value: i64 = 0;
    for &b in digits {
        if !b.is_ascii_digit() {
            return None;
        }
        value = value.checked_mul(10)?.checked_add(i64::from(b - b'0'))?;
    }
    Some(if neg { -value } else { value })
}

/// Why an INFO column cannot be indexed (doc_v26 V26-04).
///
/// A repeated tag is not a warning: the index reader, the AC/AN reader and
/// the expression evaluator each picked a different element (first or last),
/// so `AC=1;AC=3;AN=4` filtered differently depending on which rule ran. The
/// column is rejected once, at the shared index, instead of every consumer
/// choosing a convention. htslib's own readers take the first match
/// (`bcf_get_info_id`, `vcf.c:6207-6216`), which is what the retained single
/// entry then delivers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DuplicateInfoKey {
    /// Offset and length of the repeated key within the INFO column.
    pub key_start: usize,
    /// Key length in bytes.
    pub key_len: usize,
}

/// The most INFO tags one record may carry (review item 2).
///
/// The duplicate check is O(k²) in the tag count (see [`validate_info`]), so
/// an unbounded k left the worst case open: a record with tens of thousands
/// of unique tags made one `classify_line` call scan quadratically. Real VCFs
/// carry a few dozen (`DP`, `AC`, `AN`, `MQ`, `QD`, …); htslib's own ceiling is
/// `UINT16_MAX` tags (`vcf.c:3950`), which still allows a pathological scan.
///
/// 4096 is far above any real record and makes the worst case finite —
/// ~8.4M short `memcmp`s per record, milliseconds, not seconds. A record past
/// it is [`IoError::MalformedRecord`]-class malformed input, not a slow one.
pub const MAX_INFO_TAGS: usize = 4096;

/// Validates an INFO column for the record contract shared by every consumer
/// (doc_v26 V26-04): the no-duplicate-tag rule plus the tag-count ceiling.
///
/// # Cost
///
/// **Allocation-free.** The first version built a `Vec<&[u8]>` and used
/// `contains`, so every record with INFO paid a heap allocation and an
/// O(k²) scan (review item 5; its doc comment claimed "allocation-free and
/// linear", which the implementation contradicted). This version walks the
/// fragments and compares each key against the keys *before* it by re-slicing
/// `info`, so it allocates nothing.
///
/// It is O(k²) in the number of tags, which is the honest description: k is a
/// handful (`DP`, `AC`, `AN`, `MQ`, …) in every real call, and the constant is
/// a `memcmp` of short keys. A linear check would need per-worker scratch and
/// a hash of the keys, which costs more than it saves at this k. A caller
/// that already builds an [`InfoIndex`] gets the check for free from
/// [`InfoIndex::build`] (same rule, one walk); [`crate::vcf::line::classify_line`]
/// uses this function because it has no scratch of its own.
///
/// # Errors
///
/// [`DuplicateInfoKey`] when a tag repeats, and [`InfoTagOverflow`] when the
/// record carries more than [`MAX_INFO_TAGS`] tags.
pub fn validate_info(info: &[u8]) -> Result<(), InfoColumnError> {
    let mut tags = 0_usize;
    for (tag, _) in info_fragments(info) {
        tags += 1;
        // Checked before the quadratic scan below, so the bound really does
        // cap the work (review item 2).
        if tags > MAX_INFO_TAGS {
            return Err(InfoColumnError::TooManyTags {
                limit: MAX_INFO_TAGS,
            });
        }
        let key_start = tag.as_ptr() as usize - info.as_ptr() as usize;
        // Compare against every earlier key: the earlier region is
        // `info[..key_start]`, so re-walking that prefix keeps the check
        // allocation-free and identical to the index's rule.
        if prefix_has_key(&info[..key_start], tag) {
            return Err(InfoColumnError::Duplicate(DuplicateInfoKey {
                key_start,
                key_len: tag.len(),
            }));
        }
    }
    Ok(())
}

/// Why an INFO column fails validation (review item 2).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InfoColumnError {
    /// A tag repeats.
    Duplicate(DuplicateInfoKey),
    /// The record declares more tags than [`MAX_INFO_TAGS`].
    TooManyTags {
        /// The ceiling that was exceeded.
        limit: usize,
    },
}

impl From<DuplicateInfoKey> for InfoColumnError {
    fn from(err: DuplicateInfoKey) -> Self {
        Self::Duplicate(err)
    }
}

/// Whether `prefix` (the INFO text before the current fragment) already
/// carries the fragment `tag`.
///
/// `prefix` ends with the `;` that starts `tag`, so it is a run of
/// `key[=value];` entries; each entry's key is its text up to `=` or `;`.
fn prefix_has_key(prefix: &[u8], tag: &[u8]) -> bool {
    for (key, _) in info_fragments(prefix) {
        if key == tag {
            return true;
        }
    }
    false
}

/// One record's INFO fragments as byte ranges — built once per record and
/// shared by every INFO reader (typed head rules and expressions, A6-18).
#[derive(Default, Debug)]
pub struct InfoIndex {
    frags: Vec<(std::ops::Range<usize>, Option<std::ops::Range<usize>>)>,
}

impl InfoIndex {
    /// Rebuilds the ranges over `info` (empty or `.` INFO yields no
    /// entries).
    ///
    /// # Errors
    ///
    /// [`InfoColumnError`] when a tag repeats (doc_v26 V26-04) or the record
    /// carries more than [`MAX_INFO_TAGS`] tags (review item 2: the same
    /// ceiling as [`validate_info`], so the index cannot be built for a
    /// record the classifier would have rejected). The index is left empty in
    /// that case: a caller that keeps going must not filter on a partially
    /// built table.
    ///
    /// # Panics
    ///
    /// Panics only on a slice/subslice relationship the caller cannot break
    /// (`info_fragments` yields borrows of `info`).
    pub fn build(&mut self, info: &[u8]) -> Result<(), InfoColumnError> {
        self.frags.clear();
        for (tag, value) in info_fragments(info) {
            let tag_start = tag.as_ptr() as usize - info.as_ptr() as usize;
            if self.frags.len() >= MAX_INFO_TAGS {
                self.frags.clear();
                return Err(InfoColumnError::TooManyTags {
                    limit: MAX_INFO_TAGS,
                });
            }
            // Duplicate check: the fragment list is short (a handful of
            // tags), so a linear scan beats any map here.
            if self
                .frags
                .iter()
                .any(|(range, _)| &info[range.clone()] == tag)
            {
                let duplicate = DuplicateInfoKey {
                    key_start: tag_start,
                    key_len: tag.len(),
                };
                self.frags.clear();
                return Err(duplicate.into());
            }
            let value = value.map(|v| {
                let start = v.as_ptr() as usize - info.as_ptr() as usize;
                start..start + v.len()
            });
            self.frags.push((tag_start..tag_start + tag.len(), value));
        }
        Ok(())
    }

    /// Looks one tag up: `Some(None)` = flag present, `Some(Some(v))` = the
    /// value bytes, `None` = tag absent.
    pub fn get<'a>(&self, info: &'a [u8], tag: &[u8]) -> Option<Option<&'a [u8]>> {
        for (t, v) in &self.frags {
            if info[t.clone()] == *tag {
                return Some(v.as_ref().map(|r| &info[r.clone()]));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fragments_split_and_flags() {
        let fragments: Vec<_> = info_fragments(b"DP=5;SOMETHING;MQ=58.2").collect();
        assert_eq!(fragments[0], (b"DP".as_slice(), Some(b"5".as_slice())));
        assert_eq!(fragments[1], (b"SOMETHING".as_slice(), None));
        assert_eq!(fragments[2], (b"MQ".as_slice(), Some(b"58.2".as_slice())));
    }

    #[test]
    fn missing_info_is_empty() {
        assert_eq!(info_fragments(b".").count(), 0);
        assert_eq!(info_fragments(b"").count(), 0);
    }

    #[test]
    fn index_get_matches_the_scan() {
        let info = b"DP=5;SOMETHING;MQ=58.2;AC=1,2";
        let mut idx = InfoIndex::default();
        idx.build(info).expect("invariant: unique tags");
        assert_eq!(idx.get(info, b"DP"), Some(Some(b"5".as_slice())));
        assert_eq!(idx.get(info, b"SOMETHING"), Some(None));
        assert_eq!(idx.get(info, b"MQ"), Some(Some(b"58.2".as_slice())));
        assert_eq!(idx.get(info, b"AC"), Some(Some(b"1,2".as_slice())));
        assert_eq!(idx.get(info, b"NOPE"), None);
        // A rebuilt index serves the next record's bytes.
        let info2 = b"AF=0.5";
        idx.build(info2).expect("invariant: unique tags");
        assert_eq!(idx.get(info2, b"AF"), Some(Some(b"0.5".as_slice())));
        assert_eq!(idx.get(info2, b"DP"), None);
    }

    /// doc_v26 V26-04: a repeated tag is rejected by the shared validator and
    /// the index; no reader may pick the first or last element.
    #[test]
    fn duplicate_tags_are_rejected() {
        let info = b"AC=1;AN=4;AC=3";
        let err = validate_info(info).expect_err("duplicate AC must be rejected");
        let InfoColumnError::Duplicate(key) = err else {
            panic!("expected a duplicate, got {err:?}");
        };
        assert_eq!(&info[key.key_start..key.key_start + key.key_len], b"AC");
        let mut idx = InfoIndex::default();
        assert!(idx.build(info).is_err(), "the index must reject duplicates");
        assert_eq!(idx.get(info, b"AN"), None, "a failed build leaves no table");
        // The first occurrence is the one reported.
        let again = b"DP=1;DP=2";
        let err = validate_info(again).expect_err("duplicate DP must be rejected");
        let InfoColumnError::Duplicate(key) = err else {
            panic!("expected a duplicate, got {err:?}");
        };
        assert_eq!(key.key_start, 5, "DP starts right after the first `;`");
    }

    /// Review item 2: an unbounded tag count made the O(k²) duplicate scan
    /// quadratic in an attacker-controlled k. A record past `MAX_INFO_TAGS`
    /// is rejected before the scan can be exploited.
    #[test]
    fn an_excessive_tag_count_is_rejected() {
        // Exactly at the ceiling: accepted (all unique).
        use std::fmt::Write as _;
        let mut at = String::new();
        for i in 0..MAX_INFO_TAGS {
            if i > 0 {
                at.push(';');
            }
            write!(at, "T{i}=1").expect("invariant: writing to a String");
        }
        assert!(validate_info(at.as_bytes()).is_ok(), "the ceiling itself is legal");

        // One past it: rejected with the histogram-free error.
        let mut over = at.clone();
        over.push_str(";EXTRA=1");
        assert_eq!(
            validate_info(over.as_bytes()),
            Err(InfoColumnError::TooManyTags {
                limit: MAX_INFO_TAGS
            })
        );

        // The index enforces the same rule, so a caller cannot build a table
        // the classifier would have rejected.
        let mut idx = InfoIndex::default();
        assert_eq!(
            idx.build(over.as_bytes()),
            Err(InfoColumnError::TooManyTags {
                limit: MAX_INFO_TAGS
            }),
            "the index must not accept a record the classifier rejects"
        );
        assert_eq!(idx.get(over.as_bytes(), b"T0"), None, "a failed build leaves no table");
        // And the legal case still builds.
        let mut idx = InfoIndex::default();
        idx.build(at.as_bytes()).expect("invariant: at the ceiling");
        assert!(idx.get(at.as_bytes(), b"T0").is_some());
    }

    /// Flags and a repeated flag are covered as well (a Flag has no `=`).
    #[test]
    fn duplicate_flags_are_rejected() {
        assert!(validate_info(b"DB;DB").is_err());
        assert!(validate_info(b"DB;SOMETHING").is_ok());
    }
}
