//! BGZF block header parsing and the EOF marker.
//!
//! A BGZF block is one gzip member whose FEXTRA field carries the `BC`
//! subfield with the total block size. Only the 18-byte fixed header is
//! parsed on the reader's hot path — never inflated.
//!
//! [`block_entry`] is the one place that turns a BGZF member into the
//! `(compressed length, inflated length)` pair every chunk source records
//! (doc_v26 V26-01). It validates the declared ISIZE and, for a member that
//! claims to be empty, that the payload really inflates to nothing: a
//! damaged block with ISIZE zeroed used to be skipped without any CRC check,
//! silently dropping its records.

use crate::error::IoError;

/// Byte length of the fixed BGZF block header (through BSIZE).
pub const BGZF_HEADER_LEN: usize = 18;

/// Maximum uncompressed data per BGZF block, matching bgzip
/// (`docs/spec/architecture.md` §6.6).
pub const BGZF_MAX_BLOCK_DATA: usize = 65_280;

/// Safe upper bound for inflating one block (htslib uses 64 KiB buffers).
pub const BGZF_INFLATE_MAX: usize = 65_536;

/// The canonical 28-byte BGZF EOF marker (empty last block).
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Whether `block` is exactly the BGZF EOF marker.
#[must_use]
pub fn is_eof_block(block: &[u8]) -> bool {
    block == EOF_BLOCK
}

/// Whether `block` is an empty gzip member (ISIZE = 0). Empty members carry
/// no data; BGZF files concatenated with `cat` (or `bcftools concat
/// --naive`) place the 28-byte EOF marker at every join, so readers must
/// skip them instead of ending the stream (01-review D2).
#[must_use]
pub fn is_empty_member(block: &[u8]) -> bool {
    block.len() >= 8 && block[block.len() - 4..] == [0, 0, 0, 0]
}

/// The inflated size (ISIZE, the last 4 footer bytes) of a full BGZF block.
///
/// # Panics
///
/// Panics when `block` is shorter than the 8-byte footer (a caller
/// invariant: `block` always comes from [`parse_block_size`]).
#[must_use]
pub fn block_inflated_size(block: &[u8]) -> u32 {
    let n = block.len();
    u32::from_le_bytes([block[n - 4], block[n - 3], block[n - 2], block[n - 1]])
}

/// Parses the fixed header; returns the **total** block size (BSIZE + 1).
///
/// Returns `None` when the bytes are not a BGZF header (callers treat this
/// as end-of-stream or corrupt data, per context), and when BSIZE declares
/// fewer than [`BGZF_HEADER_LEN`] bytes — htslib rejects that block length
/// before touching the payload (`bgzf.c:1219-1225`); without the check the
/// callers slice or underflow on the short block (doc_v10 D10-4).
#[must_use]
pub fn parse_block_size(header: &[u8]) -> Option<usize> {
    if header.len() < BGZF_HEADER_LEN {
        return None;
    }
    let gzip_ok = header[0] == 0x1f && header[1] == 0x8b && header[2] == 0x08;
    let fextra_ok = header[3] & 0x04 != 0;
    let bc_ok = header[12] == b'B' && header[13] == b'C';
    if !(gzip_ok && fextra_ok && bc_ok) {
        return None;
    }
    let bsize = u16::from_le_bytes([header[16], header[17]]) as usize;
    if bsize + 1 < BGZF_HEADER_LEN {
        return None;
    }
    Some(bsize + 1)
}

/// One BGZF member's `(compressed length, inflated length)`, validated.
///
/// The chunk sources (streaming, memory-mapped, region) all record this pair
/// in [`ChunkBytes::blocks`](crate::chunk::ChunkBytes) and every later
/// consumer trusts it: `inflate_chunk` sizes its output window from the
/// inflated length and skips members whose length is zero. That trust used to
/// be misplaced (doc_v26 V26-01):
///
/// - a member whose ISIZE was damaged to `0` was treated as empty and skipped
///   **without inflating it**, so its records vanished while the run
///   succeeded;
/// - an ISIZE larger than [`BGZF_INFLATE_MAX`] made `inflate_chunk` grow its
///   scratch buffer to the declared size *before* the decompressor could
///   reject the block, so a 64 KiB member could request an arbitrary
///   allocation.
///
/// This entry therefore checks the declared size against the BGZF limit and,
/// for a member that claims to be empty, inflates it into a zero-length
/// window: libdeflate verifies the CRC and ISIZE of the member, an empty
/// payload yields `Ok(0)`, and any real payload fails with
/// `InsufficientSpace`. Nothing is skipped on the strength of the footer
/// alone.
///
/// # Errors
///
/// [`IoError::MalformedBgzf`] when `block` is not one whole member, declares
/// an inflated size above [`BGZF_INFLATE_MAX`], or claims to be empty while
/// its payload inflates to bytes.
pub fn block_entry(block: &[u8], inflater: &mut crate::bgzf::inflate::Inflater) -> Result<(u32, u32), IoError> {
    let Some(total) = parse_block_size(block) else {
        return Err(IoError::MalformedBgzf {
            locus: None,
            reason: "not a BGZF block header".to_string(),
        });
    };
    if total != block.len() {
        return Err(IoError::MalformedBgzf {
            locus: None,
            reason: "block length disagrees with BSIZE".to_string(),
        });
    }
    let inflated = block_inflated_size(block);
    if inflated as usize > BGZF_INFLATE_MAX {
        return Err(IoError::MalformedBgzf {
            locus: None,
            reason: format!(
                "block declares an inflated size of {inflated} bytes, over the BGZF limit of \
                 {BGZF_INFLATE_MAX}"
            ),
        });
    }
    if inflated == 0 {
        // A member claiming to be empty must really be one; the zero-length
        // window makes libdeflate prove it (CRC and ISIZE checked).
        inflater.verify_empty_member(block)?;
    }
    Ok((
        u32::try_from(total).unwrap_or(u32::MAX),
        inflated,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eof_block_is_recognized() {
        assert!(is_eof_block(&EOF_BLOCK));
        assert!(!is_eof_block(&EOF_BLOCK[..27]));
    }

    #[test]
    fn header_with_bsize_roundtrip() {
        // Header of a block whose total size is 100 bytes: BSIZE = 99.
        let mut h = [0_u8; 18];
        h[0] = 0x1f;
        h[1] = 0x8b;
        h[2] = 0x08;
        h[3] = 0x04;
        h[10] = 6;
        h[12] = b'B';
        h[13] = b'C';
        h[14] = 2;
        h[15] = 0;
        h[16] = 99;
        h[17] = 0;
        assert_eq!(parse_block_size(&h), Some(100));
        // Truncated / non-BGZF inputs yield None.
        assert_eq!(parse_block_size(&h[..17]), None);
        let mut plain = h;
        plain[3] = 0; // no FEXTRA
        assert_eq!(parse_block_size(&plain), None);
    }

    // ---- doc_v26 V26-01: the shared validated block entry ----

    use crate::bgzf::deflate::Deflater;
    use crate::bgzf::inflate::Inflater;

    /// A real BGZF member of `text`.
    fn member(text: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        Deflater::new(6)
            .expect("invariant: level 6 is valid")
            .deflate_block(text, &mut out)
            .expect("invariant: block fits");
        out
    }

    /// Rewrites a member's ISIZE footer (its last four bytes).
    fn set_isize(block: &mut [u8], isize: u32) {
        let n = block.len();
        block[n - 4..].copy_from_slice(&isize.to_le_bytes());
    }

    /// A valid member whose payload is non-empty is accepted and reports its
    /// real inflated length.
    #[test]
    fn valid_member_entry_reports_both_lengths() {
        let block = member(b"records here");
        let (clen, ilen) =
            block_entry(&block, &mut Inflater::new()).expect("invariant: valid member");
        assert_eq!(clen as usize, block.len());
        assert_eq!(ilen, 12);
    }

    /// A genuinely empty member (the EOF marker and every concatenation
    /// join) is accepted with an inflated length of zero.
    #[test]
    fn valid_empty_member_entry_is_zero() {
        let (clen, ilen) =
            block_entry(&EOF_BLOCK, &mut Inflater::new()).expect("invariant: EOF marker");
        assert_eq!(clen as usize, EOF_BLOCK.len());
        assert_eq!(ilen, 0);
        let empty = member(b"");
        assert_eq!(
            block_entry(&empty, &mut Inflater::new())
                .expect("invariant: empty member")
                .1,
            0
        );
    }

    /// The core V26-01 defect: zeroing a non-empty member's ISIZE used to
    /// make every reader skip it and drop its records. The footer is not
    /// trusted — the payload must prove it is empty.
    #[test]
    fn isize_zeroed_on_a_nonempty_member_is_malformed() {
        let mut block = member(b"chr1\t10\t.\tA\tG\t1\tPASS\t.\tGT\t0/0\n");
        set_isize(&mut block, 0);
        let err = block_entry(&block, &mut Inflater::new())
            .expect_err("a non-empty payload behind ISIZE=0 must be MalformedBgzf");
        assert!(matches!(err, IoError::MalformedBgzf { .. }), "{err:?}");
    }

    /// An ISIZE above the BGZF limit must fail at the entry, before any
    /// consumer sizes a buffer from it.
    #[test]
    fn oversized_isize_is_malformed_before_any_allocation() {
        let mut block = member(b"0123456789");
        set_isize(&mut block, 16 * 1024 * 1024);
        let err = block_entry(&block, &mut Inflater::new())
            .expect_err("ISIZE over the BGZF limit must be MalformedBgzf");
        assert!(matches!(err, IoError::MalformedBgzf { .. }), "{err:?}");
    }

    /// A non-empty member's payload bytes are not decompressed a second time
    /// just to cut chunks (that would double the inflate cost on the reader
    /// thread). Its CRC is verified by the inflater when the chunk is
    /// decoded, which every record path runs; this test pins that division of
    /// labour so it cannot be mistaken for "validation never happens".
    #[test]
    fn nonempty_payload_crc_is_checked_at_inflate() {
        let text = b"a readable record line\n";
        let mut block = member(text);
        // Damage the CRC footer (4 bytes before ISIZE), leaving the sizes.
        let n = block.len();
        block[n - 8] ^= 0xff;
        assert!(
            block_entry(&block, &mut Inflater::new()).is_ok(),
            "the chunk-cutting entry must not inflate a non-empty block"
        );
        let mut out = vec![0_u8; text.len()];
        assert!(
            Inflater::new().inflate_window(&block, &mut out).is_err(),
            "a damaged payload must fail when the chunk is decoded"
        );
    }

    /// A member that is not one whole block (trailing or missing bytes)
    /// fails at the entry.
    #[test]
    fn truncated_member_is_malformed() {
        let block = member(b"whole member");
        let err = block_entry(&block[..block.len() - 1], &mut Inflater::new())
            .expect_err("a short member must be MalformedBgzf");
        assert!(matches!(err, IoError::MalformedBgzf { .. }), "{err:?}");
    }
}
