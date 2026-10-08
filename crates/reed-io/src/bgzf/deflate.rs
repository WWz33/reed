//! libdeflate-based BGZF block compression (`docs/architecture-spec.md`
//! §6.6): one call produces one complete BGZF member.

// Casts in this file are bounded by construction: BGZF blocks are at most 65,536 bytes by construction.
#![allow(clippy::cast_possible_truncation)]

use crate::bgzf::block::{BGZF_MAX_BLOCK_DATA, EOF_BLOCK};
use crate::error::IoError;

/// Worker-local compressor: allocate once per worker, reuse per block.
pub struct Deflater {
    compressor: libdeflater::Compressor,
}

impl Deflater {
    /// Creates a compressor at the given deflate level (6 = bgzip default).
    ///
    /// # Errors
    ///
    /// Fails for levels outside libdeflate's 1..=12 range.
    pub fn new(level: i32) -> Result<Self, IoError> {
        let lvl = libdeflater::CompressionLvl::new(level)
            .map_err(|_| IoError::UnsupportedInput("invalid compression level"))?;
        Ok(Self {
            compressor: libdeflater::Compressor::new(lvl),
        })
    }

    /// Compresses `src` into one BGZF block appended to `dst`.
    ///
    /// Returns the total block length written. Blocks never exceed the
    /// 64 KiB gzip-member ceiling: callers keep `src` at or below
    /// [`crate::bgzf::block::BGZF_MAX_BLOCK_DATA`] uncompressed bytes.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedBgzf`] when `src` exceeds
    /// [`BGZF_MAX_BLOCK_DATA`] — writing it would corrupt the member's
    /// 16-bit BSIZE (01-review D11) — or if libdeflate reports an output
    /// overflow.
    pub fn deflate_block(&mut self, src: &[u8], dst: &mut Vec<u8>) -> Result<usize, IoError> {
        if src.len() > BGZF_MAX_BLOCK_DATA {
            return Err(IoError::MalformedBgzf {
                locus: None,
                reason: "block data exceeds the BGZF 65,280-byte member limit",
            });
        }
        let start = dst.len();
        // Raw deflate payload with libdeflate's own worst-case bound, plus
        // the 18-byte header and the 8-byte CRC/ISIZE trailer.
        let payload_start = start + 18;
        let bound = self.compressor.deflate_compress_bound(src.len());
        dst.resize(payload_start + bound + 8, 0);
        let n_deflated = self
            .compressor
            .deflate_compress(src, &mut dst[payload_start..])
            .map_err(|_| IoError::MalformedBgzf {
                locus: None,
                reason: "deflate output overflow (block data too large)",
            })?;
        let total = 18 + n_deflated + 8;
        // Header: gzip with FEXTRA carrying the BC subfield (XLEN = 6).
        dst[start..start + 18].copy_from_slice(&[
            0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0, 0,
        ]);
        let bsize = u16::try_from(total - 1).map_err(|_| IoError::MalformedBgzf {
            locus: None,
            reason: "BGZF block size overflows the 16-bit BSIZE field",
        })?;
        dst[start + 16..start + 18].copy_from_slice(&bsize.to_le_bytes());
        // Trailer: CRC32 of the uncompressed data + ISIZE mod 2^32.
        let crc = libdeflater::crc32(src);
        let trailer = start + 18 + n_deflated;
        dst[trailer..trailer + 4].copy_from_slice(&crc.to_le_bytes());
        dst[trailer + 4..trailer + 8].copy_from_slice(&(src.len() as u32).to_le_bytes());
        dst.truncate(start + total);
        Ok(total)
    }
}

/// The 28-byte EOF block every BGZF file must end with.
#[must_use]
pub const fn eof_block() -> [u8; 28] {
    EOF_BLOCK
}

#[cfg(test)]
mod tests {
    use super::super::inflate::Inflater;
    use super::*;

    #[test]
    fn oversize_input_is_rejected_not_truncated() {
        // Regression (01-review D11): 65,281 uncompressed bytes used to be
        // written with a truncated BSIZE, corrupting the block.
        let mut deflater = Deflater::new(6).expect("invariant: level 6 is valid");
        let mut dst = Vec::new();
        let src = vec![b'x'; BGZF_MAX_BLOCK_DATA + 1];
        let err = deflater.deflate_block(&src, &mut dst).unwrap_err();
        assert!(matches!(err, IoError::MalformedBgzf { .. }), "got {err:?}");
        // Exactly at the limit still works.
        dst.clear();
        let src = vec![b'x'; BGZF_MAX_BLOCK_DATA];
        deflater.deflate_block(&src, &mut dst).expect("invariant: max-size block fits");
    }

    #[test]
    fn block_structure_is_self_consistent() {
        let mut deflater = Deflater::new(6).expect("invariant: level 6 is valid");
        let mut dst = Vec::new();
        let total = deflater.deflate_block(b"0123456789", &mut dst).expect("invariant: small block");
        assert_eq!(total, dst.len());
        assert_eq!(total, usize::from(u16::from_le_bytes([dst[16], dst[17]])) + 1);
        assert_eq!(&dst[..4], &[0x1f, 0x8b, 0x08, 0x04]);
        assert_eq!(&dst[12..14], b"BC");
    }

    #[test]
    fn empty_block_deflates_and_reinflates() {
        let mut deflater = Deflater::new(6).expect("invariant: level 6 is valid");
        let mut dst = Vec::new();
        deflater.deflate_block(b"", &mut dst).expect("invariant: empty block");
        let mut inflater = Inflater::new();
        let mut out = Vec::new();
        let n = inflater.inflate_block(&dst, &mut out).expect("invariant: empty roundtrip");
        assert_eq!(n, 0);
    }
}

#[cfg(test)]
mod noodles_crosscheck {
    use std::io::Read;

    #[test]
    fn noodles_reads_our_block() {
        let mut deflater = super::Deflater::new(6).expect("invariant: level 6 is valid");
        let mut dst = Vec::new();
        deflater.deflate_block(b"noodles crosscheck payload", &mut dst).expect("invariant: small block");
        let mut reader = noodles::bgzf::io::Reader::new(&dst[..]);
        let mut out = Vec::new();
        reader.read_to_end(&mut out).expect("invariant: noodles can read our block");
        assert_eq!(out, b"noodles crosscheck payload");
    }
}
