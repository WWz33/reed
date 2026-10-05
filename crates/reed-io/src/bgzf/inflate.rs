//! libdeflate-based BGZF block decompression with worker-local reuse
//! (`docs/architecture-spec.md` §6.2).

use crate::bgzf::block::{parse_block_size, BGZF_INFLATE_MAX};
use crate::error::IoError;

/// Worker-local decompressor: allocate once per worker, reuse per block.
pub struct Inflater {
    decompressor: libdeflater::Decompressor,
}

impl Default for Inflater {
    fn default() -> Self {
        Self::new()
    }
}

impl Inflater {
    /// Creates a decompressor (libdeflate verifies the gzip CRC).
    #[must_use]
    pub fn new() -> Self {
        Self {
            decompressor: libdeflater::Decompressor::new(),
        }
    }

    /// Inflates one BGZF block, appending the decompressed bytes to `dst`.
    ///
    /// Returns the number of bytes appended.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedBgzf`] on a corrupt member or a block whose
    /// inflated size exceeds the BGZF limit.
    pub fn inflate_block(&mut self, block: &[u8], dst: &mut Vec<u8>) -> Result<usize, IoError> {
        let Some(total) = parse_block_size(block) else {
            return Err(IoError::MalformedBgzf {
                locus: None,
                reason: "not a BGZF block header",
            });
        };
        if total != block.len() {
            return Err(IoError::MalformedBgzf {
                locus: None,
                reason: "block length disagrees with BSIZE",
            });
        }
        let start = dst.len();
        dst.resize(start + BGZF_INFLATE_MAX, 0);
        // libdeflate parses the gzip member itself: pass the whole block
        // (header included), not just the deflate payload.
        let written = self.decompressor.gzip_decompress(&block[..total], &mut dst[start..])?;
        dst.truncate(start + written);
        Ok(written)
    }
}

#[cfg(test)]
mod tests {
    use super::super::deflate::Deflater;
    use super::*;

    #[test]
    fn inflate_roundtrips_deflated_block() {
        let mut deflater = Deflater::new(6).expect("invariant: level 6 is valid");
        let mut compressed = Vec::new();
        deflater.deflate_block(b"hello bgzf world", &mut compressed).expect("invariant: deflate of 16 bytes fits");

        let mut inflater = Inflater::new();
        let mut out = Vec::new();
        let n = inflater.inflate_block(&compressed, &mut out).expect("invariant: roundtrip");
        assert_eq!(n, 16);
        assert_eq!(out, b"hello bgzf world");
    }

    #[test]
    fn rejects_non_bgzf_block() {
        let mut inflater = Inflater::new();
        let mut out = Vec::new();
        assert!(inflater.inflate_block(b"garbage bytes here!!", &mut out).is_err());
    }
}
