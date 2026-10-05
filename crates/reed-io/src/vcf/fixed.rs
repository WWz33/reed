//! Fixed-column fast parsers (`docs/architecture-spec.md` §6.4): POS is a
//! hand-written decimal; QUAL uses std float parsing; FILTER is a byte
//! comparison.

use crate::error::IoError;

/// Parses the POS column (1-based unsigned decimal).
///
/// # Errors
///
/// [`IoError::MalformedRecord`] on empty or non-digit text.
pub fn parse_pos1(text: &[u8], locus: &crate::Locus) -> Result<u64, IoError> {
    if text.is_empty() || text.len() > 19 {
        return Err(IoError::MalformedRecord {
            locus: locus.clone(),
            reason: "POS is empty or too long",
        });
    }
    let mut pos1: u64 = 0;
    for &b in text {
        if !b.is_ascii_digit() {
            return Err(IoError::MalformedRecord {
                locus: locus.clone(),
                reason: "POS contains a non-digit",
            });
        }
        pos1 = pos1 * 10 + u64::from(b - b'0');
    }
    Ok(pos1)
}

/// Parses the QUAL column; `.` is no quality score.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] when the text is neither `.` nor a float.
pub fn parse_qual(text: &[u8], locus: &crate::Locus) -> Result<Option<f32>, IoError> {
    if text == b"." || text.is_empty() {
        return Ok(None);
    }
    let s = std::str::from_utf8(text)
        .map_err(|_| IoError::MalformedRecord { locus: locus.clone(), reason: "QUAL is not ASCII" })?;
    s.parse::<f32>()
        .map(Some)
        .map_err(|_| IoError::MalformedRecord { locus: locus.clone(), reason: "QUAL is not a number" })
}

/// Whether the FILTER column passes (`PASS`, byte-exact).
#[must_use]
pub fn is_pass(text: &[u8]) -> bool {
    text == b"PASS"
}

#[cfg(test)]
mod tests {
    use super::*;

    fn locus() -> crate::Locus {
        crate::Locus { contig_name: "chr1".into(), pos1: 1 }
    }

    #[test]
    fn pos_parses_and_rejects() {
        let l = locus();
        assert_eq!(parse_pos1(b"12345", &l).unwrap(), 12_345);
        assert_eq!(parse_pos1(b"0", &l).unwrap(), 0);
        assert!(parse_pos1(b"", &l).is_err());
        assert!(parse_pos1(b"12x", &l).is_err());
        assert!(parse_pos1(b"12345678901234567890", &l).is_err());
    }

    #[test]
    fn qual_dot_is_none() {
        let l = locus();
        assert_eq!(parse_qual(b".", &l).unwrap(), None);
        assert_eq!(parse_qual(b"99.5", &l).unwrap(), Some(99.5));
        assert_eq!(parse_qual(b"0", &l).unwrap(), Some(0.0));
        assert!(parse_qual(b"abc", &l).is_err());
    }

    #[test]
    fn pass_is_byte_exact() {
        assert!(is_pass(b"PASS"));
        assert!(!is_pass(b"pass"));
        assert!(!is_pass(b"PASS;LowQual"));
        assert!(!is_pass(b"."));
    }
}
