//! BGZF block header parsing and the EOF marker.
//!
//! A BGZF block is one gzip member whose FEXTRA field carries the `BC`
//! subfield with the total block size. Only the 18-byte fixed header is
//! parsed on the reader's hot path — never inflated.

/// Byte length of the fixed BGZF block header (through BSIZE).
pub const BGZF_HEADER_LEN: usize = 18;

/// Maximum uncompressed data per BGZF block, matching bgzip
/// (`docs/architecture-spec.md` §6.6).
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
/// as end-of-stream or corrupt data, per context).
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
    Some(bsize + 1)
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
}
