//! BGZF virtual offsets and the per-chunk block table that maps inflated text
//! back to them (doc_v28 V28-27).
//!
//! A virtual offset packs a block's compressed file offset and an offset
//! inside that block's inflated payload into one `u64`:
//!
//! ```text
//! (compressed_offset << 16) | in_block_offset
//! ```
//!
//! The low 16 bits hold the in-block offset, which is why a BGZF payload may
//! not exceed 65,536 bytes. Two exact copies of the text-offset → virtual
//! offset walk existed (the BCF record mapper and the index job's); this module
//! is their single home. It knows BGZF addressing only — not chunks, records,
//! indexes, or regions (those stay with the callers).

/// A block's compressed length and its inflated length, in file order.
///
/// Empty members (the EOF marker, a concatenation join) appear as
/// `(clen, 0)` so compressed offsets stay contiguous.
pub type BlockEntry = (u32, u32);

/// The virtual offset of a position in the input file.
///
/// Ordering is `(compressed, in-block)`, which for well-formed BGZF is the
/// same as raw `u64` ordering because the in-block half is a `u16`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VirtualOffset(u64);

impl VirtualOffset {
    /// Packs a compressed offset and an in-block offset.
    ///
    /// # Panics
    ///
    /// Panics when `in_block` does not fit the low 16 bits — a caller
    /// invariant: a BGZF block's inflated length is at most 65,536 bytes.
    #[must_use]
    pub const fn new(compressed: u64, in_block: u16) -> Self {
        Self((compressed << 16) | in_block as u64)
    }

    /// The offset at the start of a block (or of the file, for an EOF
    /// position): the compressed half shifted up, in-block zero.
    #[must_use]
    pub const fn at_block_start(compressed: u64) -> Self {
        Self(compressed << 16)
    }

    /// The offset just past a block whose inflated length is `inflated_len`.
    ///
    /// # Panics
    ///
    /// Panics when `inflated_len` does not fit the low 16 bits — a caller
    /// invariant: a BGZF block inflates to at most 65,536 bytes.
    #[must_use]
    pub const fn at_block_end(compressed: u64, inflated_len: u16) -> Self {
        Self((compressed << 16) | inflated_len as u64)
    }

    /// The raw packed value (the representation indexes serialise).
    #[must_use]
    pub const fn from_raw(raw: u64) -> Self {
        Self(raw)
    }

    /// The packed value.
    #[must_use]
    pub const fn raw(self) -> u64 {
        self.0
    }

    /// The block's compressed file offset (the high 48 bits).
    #[must_use]
    pub const fn compressed(self) -> u64 {
        self.0 >> 16
    }

    /// The offset inside the block's inflated payload (the low 16 bits).
    #[must_use]
    // The low 16 bits are exactly the in-block field; the cast is the mask.
    #[allow(clippy::cast_possible_truncation)]
    pub const fn in_block(self) -> u16 {
        self.0 as u16
    }
}

/// A read-only view over one chunk's block table.
///
/// The table lists every block in the chunk in file order, including empty
/// members, so walking it reproduces the writer's compressed offsets exactly.
#[derive(Clone, Copy, Debug)]
pub struct BlockTable<'a> {
    blocks: &'a [BlockEntry],
}

impl<'a> BlockTable<'a> {
    /// A view over `blocks`, which must be in file order.
    #[must_use]
    pub const fn new(blocks: &'a [BlockEntry]) -> Self {
        Self { blocks }
    }

    /// The virtual offset of `text_off` in the chunk's inflated text, where
    /// `chunk_start` is the compressed offset of the chunk's first block.
    ///
    /// An offset past the table clamps to the block after the last one — the
    /// chunk's end — which is what both previous copies did for an empty
    /// trailing chunk.
    #[must_use]
    pub fn voff(&self, chunk_start: u64, text_off: usize) -> VirtualOffset {
        let mut remaining = text_off as u64;
        let mut compressed = chunk_start;
        for &(clen, ilen) in self.blocks {
            if remaining < u64::from(ilen) {
                // SAFETY: a BGZF block's inflated size is at most 65,536, so
                // any in-block offset below it fits the low 16 bits.
                #[allow(clippy::cast_possible_truncation)]
                let in_block = remaining as u16;
                return VirtualOffset::new(compressed, in_block);
            }
            remaining -= u64::from(ilen);
            compressed += u64::from(clen);
        }
        VirtualOffset::from_raw(compressed << 16)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packing_round_trips() {
        let voff = VirtualOffset::new(0x0123_4567_89ab, 0xcdef);
        assert_eq!(voff.compressed(), 0x0123_4567_89ab);
        assert_eq!(voff.in_block(), 0xcdef);
        assert_eq!(VirtualOffset::from_raw(voff.raw()), voff);
    }

    /// A text offset inside a block resolves to that block's compressed
    /// start plus the offset; past the last block it clamps to the end.
    #[test]
    fn voff_walks_the_table_and_clamps_past_it() {
        // Two blocks: 100 compressed -> 10 inflated, 200 -> 20.
        let table = BlockTable::new(&[(100, 10), (200, 20)]);
        assert_eq!(table.voff(1000, 0).raw(), 1000 << 16);
        assert_eq!(table.voff(1000, 9).raw(), (1000 << 16) | 0x09);
        assert_eq!(table.voff(1000, 10).raw(), 1100 << 16);
        assert_eq!(table.voff(1000, 29).raw(), (1100 << 16) | 0x13);
        // Past the table: the chunk's end (after both blocks).
        assert_eq!(table.voff(1000, 30).raw(), 1300 << 16);
        assert_eq!(table.voff(1000, 1 << 20).raw(), 1300 << 16);
    }

    /// An empty member `(clen, 0)` advances only the compressed offset, so
    /// the virtual offsets of later records stay correct.
    #[test]
    fn empty_members_keep_the_compressed_offsets_contiguous() {
        let table = BlockTable::new(&[(50, 5), (28, 0), (70, 7)]);
        // Text offset 5 skips the empty member and lands at the start of the
        // third block, whose compressed offset is 50 + 28.
        assert_eq!(table.voff(0, 5).raw(), (50 + 28) << 16);
        // Offset 11 is 6 bytes into that block's 7-byte payload.
        assert_eq!(table.voff(0, 11).raw(), ((50 + 28) << 16) | 0x06);
        // Offset 12 is exactly the end, which has no in-block position: the
        // walk falls through to the block after the last one.
        assert_eq!(table.voff(0, 12).raw(), (50 + 28 + 70) << 16);
        // Well past the table: the same chunk-end clamp.
        assert_eq!(table.voff(0, 12 + (1 << 20)).raw(), (50 + 28 + 70) << 16);
    }
}
