//! Fixed-column fast parsers (`docs/spec/architecture.md` §6.4): POS is a
//! hand-written decimal; QUAL uses std float parsing; FILTER is a byte
//! comparison.

use crate::error::IoError;

/// Parses the POS column (1-based unsigned decimal). `contig_name` only
/// feeds the error locus, which is built when the error actually happens
/// (A6-11: no per-record `String` on the happy path).
///
/// # Errors
///
/// [`IoError::MalformedRecord`] on empty or non-digit text.
pub fn parse_pos1(text: &[u8], contig_name: &[u8]) -> Result<u64, IoError> {
    let locus = || crate::Locus {
        contig_name: String::from_utf8_lossy(contig_name).into_owned(),
        pos1: 0,
    };
    if text.is_empty() || text.len() > 19 {
        return Err(IoError::MalformedRecord {
            locus: locus(),
            reason: "POS is empty or too long",
        });
    }
    let mut pos1: u64 = 0;
    for &b in text {
        if !b.is_ascii_digit() {
            return Err(IoError::MalformedRecord {
                locus: locus(),
                reason: "POS contains a non-digit",
            });
        }
        pos1 = pos1 * 10 + u64::from(b - b'0');
    }
    Ok(pos1)
}

/// Parses the QUAL column; `.` is no quality score. `locus` is the record's
/// already-parsed `(contig, POS)` pair, borrowed for the error path only.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] when the text is neither `.` nor a float.
pub fn parse_qual(text: &[u8], locus: (&[u8], u64)) -> Result<Option<f32>, IoError> {
    if text == b"." || text.is_empty() {
        return Ok(None);
    }
    let err = |reason: &'static str| IoError::MalformedRecord {
        locus: crate::Locus {
            contig_name: String::from_utf8_lossy(locus.0).into_owned(),
            pos1: locus.1,
        },
        reason,
    };
    let s = std::str::from_utf8(text).map_err(|_| err("QUAL is not ASCII"))?;
    s.parse::<f32>()
        .map(Some)
        .map_err(|_| err("QUAL is not a number"))
}

/// Whether the FILTER column passes (`PASS`, byte-exact).
#[must_use]
pub fn is_pass(text: &[u8]) -> bool {
    text == b"PASS"
}

#[cfg(test)]
mod tests {
    use super::*;

    const CHROM: &[u8] = b"chr1";

    #[test]
    fn pos_parses_and_rejects() {
        assert_eq!(parse_pos1(b"12345", CHROM).unwrap(), 12_345);
        assert_eq!(parse_pos1(b"0", CHROM).unwrap(), 0);
        assert!(parse_pos1(b"", CHROM).is_err());
        assert!(parse_pos1(b"12x", CHROM).is_err());
        assert!(parse_pos1(b"12345678901234567890", CHROM).is_err());
    }

    #[test]
    fn qual_dot_is_none() {
        assert_eq!(parse_qual(b".", (CHROM, 1)).unwrap(), None);
        assert_eq!(parse_qual(b"99.5", (CHROM, 1)).unwrap(), Some(99.5));
        assert_eq!(parse_qual(b"0", (CHROM, 1)).unwrap(), Some(0.0));
        assert!(parse_qual(b"abc", (CHROM, 1)).is_err());
    }

    #[test]
    fn pass_is_byte_exact() {
        assert!(is_pass(b"PASS"));
        assert!(!is_pass(b"pass"));
        assert!(!is_pass(b"PASS;LowQual"));
        assert!(!is_pass(b"."));
    }
}
