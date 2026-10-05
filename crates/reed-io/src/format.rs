//! Format detection by magic bytes (`docs/architecture-spec.md` §6.1).
//!
//! Extension names are only for naming outputs; the bytes decide how to
//! read. The caller hands over the first BGZF/gzip block's raw bytes and,
//! when it is BGZF, the inflated first block (needed to tell BCF from
//! BGZF-compressed VCF).

use crate::error::IoError;

/// The supported input formats.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    /// Uncompressed VCF text (mmap-able, parallel by byte ranges).
    Vcf,
    /// BGZF-compressed VCF (block-aligned parallel parsing).
    VcfBgzf,
    /// Regular gzip VCF (serial inflate only; warn and suggest `bgzip`).
    VcfGzip,
    /// BCF (binary; needs an index or serial delimiting).
    Bcf,
}

impl Format {
    /// Detects the format from the first block's raw bytes plus, for BGZF,
    /// its inflated content. `raw` must start at the file offset where
    /// detection began (offset 0 in practice).
    ///
    /// # Errors
    ///
    /// [`IoError::UnsupportedInput`] when no known magic matches.
    pub fn detect(raw_first_block: &[u8], inflated_first_block: Option<&[u8]>) -> Result<Self, IoError> {
        if raw_first_block.starts_with(b"##fileformat=VCF") {
            return Ok(Self::Vcf);
        }
        if raw_first_block.len() >= 2 && raw_first_block[0] == 0x1f && raw_first_block[1] == 0x8b {
            if has_bc_extra_field(raw_first_block) {
                // BGZF: the inflated content tells VCF from BCF.
                let Some(text) = inflated_first_block else {
                    return Err(IoError::UnsupportedInput(
                        "BGZF member detected but no inflated content supplied for format sniffing",
                    ));
                };
                if text.starts_with(b"BCF\x02") {
                    return Ok(Self::Bcf);
                }
                if text.starts_with(b"##fileformat=VCF") {
                    return Ok(Self::VcfBgzf);
                }
                return Err(IoError::UnsupportedInput(
                    "BGZF member is neither BCF nor VCF",
                ));
            }
            // Plain gzip, not BGZF.
            return Ok(Self::VcfGzip);
        }
        Err(IoError::UnsupportedInput(
            "input is neither VCF text nor gzip/BGZF",
        ))
    }

    /// Extension-name fallback for naming outputs (never for reading).
    #[must_use]
    pub fn from_extension(name: &str) -> Option<Self> {
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".bcf") {
            Some(Self::Bcf)
        } else if lower.ends_with(".vcf.gz") || lower.ends_with(".vcf.bgz") {
            Some(Self::VcfBgzf)
        } else if lower.ends_with(".gz") {
            Some(Self::VcfGzip)
        } else if lower.ends_with(".vcf") {
            Some(Self::Vcf)
        } else {
            None
        }
    }
}

/// Whether the gzip member carries the BGZF `BC` extra subfield
/// (FEXTRA flag set, subfield SI1='B' SI2='C' anywhere in the extra field).
fn has_bc_extra_field(raw: &[u8]) -> bool {
    if raw.len() < 18 || raw[3] & 0x04 == 0 {
        return false;
    }
    let xlen = u16::from_le_bytes([raw[10], raw[11]]) as usize;
    let extra_end = 12 + xlen.min(raw.len().saturating_sub(12));
    let mut i = 12;
    while i + 4 <= extra_end {
        let slen = u16::from_le_bytes([raw[i + 2], raw[i + 3]]) as usize;
        if raw[i] == b'B' && raw[i + 1] == b'C' && slen == 2 {
            return true;
        }
        i += 4 + slen;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal 18-byte BGZF header template (BSIZE placeholder).
    fn bgzf_header() -> [u8; 18] {
        let mut h = [0; 18];
        h[0] = 0x1f;
        h[1] = 0x8b;
        h[2] = 8;
        h[3] = 4; // FEXTRA
        h[10] = 6; // XLEN
        h[12] = b'B';
        h[13] = b'C';
        h[14] = 2;
        h[15] = 0;
        h
    }

    #[test]
    fn detects_vcf_bgzf() {
        let mut raw = bgzf_header().to_vec();
        raw.resize(60, 0);
        assert_eq!(
            Format::detect(&raw, Some(b"##fileformat=VCFv4.3\n")).unwrap(),
            Format::VcfBgzf
        );
    }

    #[test]
    fn detects_bcf() {
        let mut raw = bgzf_header().to_vec();
        raw.resize(60, 0);
        assert_eq!(
            Format::detect(&raw, Some(b"BCF\x02\x02")).unwrap(),
            Format::Bcf
        );
    }

    #[test]
    fn detects_plain_gzip_and_text() {
        let mut raw = [0_u8; 60];
        raw[0] = 0x1f;
        raw[1] = 0x8b;
        raw[2] = 8;
        raw[3] = 0; // no FEXTRA
        assert_eq!(Format::detect(&raw, None).unwrap(), Format::VcfGzip);
        assert_eq!(
            Format::detect(b"##fileformat=VCFv4.3", None).unwrap(),
            Format::Vcf
        );
        assert!(Format::detect(b"not a vcf at all........", None).is_err());
    }

    #[test]
    fn gzip_with_non_bc_extra_is_plain_gzip() {
        let mut raw = [0_u8; 60];
        raw[0] = 0x1f;
        raw[1] = 0x8b;
        raw[2] = 8;
        raw[3] = 4; // FEXTRA present
        raw[10] = 4; // XLEN=4: subfield XX len 0
        raw[12] = b'X';
        raw[13] = b'X';
        raw[14] = 0;
        raw[15] = 0;
        assert_eq!(Format::detect(&raw, None).unwrap(), Format::VcfGzip);
    }

    #[test]
    fn extension_fallback() {
        assert_eq!(Format::from_extension("x.vcf.gz"), Some(Format::VcfBgzf));
        assert_eq!(Format::from_extension("x.BCF"), Some(Format::Bcf));
        assert_eq!(Format::from_extension("x.vcf"), Some(Format::Vcf));
        assert_eq!(Format::from_extension("x.gz"), Some(Format::VcfGzip));
        assert_eq!(Format::from_extension("x.txt"), None);
    }
}
