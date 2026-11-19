//! Fixed-column value parsing that the policy layer also needs (doc_v30
//! R12).
//!
//! The transport keeps the located parsers (`parse_pos1`, the error-returning
//! `parse_qual` wrapper) in `reed-io::vcf::fixed`; this module owns the value
//! grammar so the filter and the transport cannot disagree about what a QUAL
//! string means.

/// What the QUAL column held.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum QualValue {
    /// `.` or an empty column: no quality score.
    Missing,
    /// A float as written (the VCF grammar is Rust's `f32` parser).
    Value(f32),
    /// Not ASCII, or not a float.
    Malformed,
}

/// Parses the POS column as a 1-based unsigned decimal (`None` for empty,
/// over-long or non-digit text).
///
/// The grammar lives here so the record classifier and the transport parser
/// agree on what POS text is valid; `reed_io::vcf::fixed::parse_pos1` adds the
/// located error the transport reports.
#[must_use]
pub fn parse_pos1_number(text: &[u8]) -> Option<u64> {
    if text.is_empty() || text.len() > 19 {
        return None;
    }
    let mut pos1: u64 = 0;
    for &b in text {
        if !b.is_ascii_digit() {
            return None;
        }
        // At most 19 digits: the product cannot overflow `u64`.
        pos1 = pos1 * 10 + u64::from(b - b'0');
    }
    Some(pos1)
}

/// Parses the QUAL column (`.`/empty = [`QualValue::Missing`]).
#[must_use]
pub fn parse_qual(text: &[u8]) -> QualValue {
    if text == b"." || text.is_empty() {
        return QualValue::Missing;
    }
    let Ok(text) = std::str::from_utf8(text) else {
        return QualValue::Malformed;
    };
    match text.parse::<f32>() {
        Ok(value) => QualValue::Value(value),
        Err(_) => QualValue::Malformed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pos_grammar_is_shared() {
        assert_eq!(parse_pos1_number(b"12345"), Some(12_345));
        assert_eq!(parse_pos1_number(b"0"), Some(0));
        assert_eq!(parse_pos1_number(b""), None);
        assert_eq!(parse_pos1_number(b"12x"), None);
        assert_eq!(parse_pos1_number(b"12345678901234567890"), None);
    }

    #[test]
    fn qual_grammar_is_shared() {
        assert_eq!(parse_qual(b"."), QualValue::Missing);
        assert_eq!(parse_qual(b""), QualValue::Missing);
        assert_eq!(parse_qual(b"99.5"), QualValue::Value(99.5));
        assert_eq!(parse_qual(b"0"), QualValue::Value(0.0));
        assert_eq!(parse_qual(b"abc"), QualValue::Malformed);
        assert_eq!(parse_qual(&[0xFF]), QualValue::Malformed);
    }
}
