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

/// One record's INFO fragments as byte ranges — built once per record and
/// shared by every INFO reader (typed head rules and expressions, A6-18).
#[derive(Default, Debug)]
pub struct InfoIndex {
    frags: Vec<(std::ops::Range<usize>, Option<std::ops::Range<usize>>)>,
}

impl InfoIndex {
    /// Rebuilds the ranges over `info` (empty or `.` INFO yields no
    /// entries).
    pub fn build(&mut self, info: &[u8]) {
        self.frags.clear();
        for (tag, value) in info_fragments(info) {
            let tag_start = tag.as_ptr() as usize - info.as_ptr() as usize;
            let value = value.map(|v| {
                let start = v.as_ptr() as usize - info.as_ptr() as usize;
                start..start + v.len()
            });
            self.frags.push((tag_start..tag_start + tag.len(), value));
        }
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
        idx.build(info);
        assert_eq!(idx.get(info, b"DP"), Some(Some(b"5".as_slice())));
        assert_eq!(idx.get(info, b"SOMETHING"), Some(None));
        assert_eq!(idx.get(info, b"MQ"), Some(Some(b"58.2".as_slice())));
        assert_eq!(idx.get(info, b"AC"), Some(Some(b"1,2".as_slice())));
        assert_eq!(idx.get(info, b"NOPE"), None);
        // A rebuilt index serves the next record's bytes.
        let info2 = b"AF=0.5";
        idx.build(info2);
        assert_eq!(idx.get(info2, b"AF"), Some(Some(b"0.5".as_slice())));
        assert_eq!(idx.get(info2, b"DP"), None);
    }
}
