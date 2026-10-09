//! Worker-local BGZF block writer (`docs/rounds/v2/03-architecture.md` §5).
//!
//! Block boundaries are independent of record boundaries, exactly like
//! htslib: records are appended to a pending buffer and cut whenever it
//! reaches the 65,280-byte member limit, so a record wider than one member
//! spans blocks instead of corrupting one (01-review D11).
//!
//! Every record's index entry anchors at the block and offset where its
//! bytes start; a record whose start falls after a cut is re-based into the
//! continuation block.

use reed_core::ContigId;

use crate::bgzf::deflate::Deflater;
use crate::error::IoError;

/// BGZF maximum uncompressed payload per block (htslib `BGZF_BLOCK_SIZE`).
pub const BGZF_BLOCK_SIZE: usize = 65_280;

/// Where one record starts in the index key space.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IndexKey {
    /// Contig id (header `##contig` order).
    pub contig: ContigId,
    /// 1-based record start.
    pub start1: u64,
    /// Reference allele length: the index bin span is
    /// `[start1, start1 + ref_len − 1]`, the htslib convention.
    pub ref_len: u32,
}

/// One record start anchored in an output block.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IndexRecord {
    /// The record's index key.
    pub key: IndexKey,
    /// Byte offset of the record's first byte within the block.
    pub offset_in_block: u32,
}

/// One finished output block.
#[derive(Clone, Debug, Default)]
pub struct OutBlock {
    /// Compressed BGZF member (empty on the plain path).
    pub compressed: Vec<u8>,
    /// Uncompressed payload length.
    pub uncompressed_len: u32,
    /// Record starts that fall in this block, in record order.
    pub records: Vec<IndexRecord>,
}

/// One chunk's output: blocks (compressed path) or raw bytes (plain path).
#[derive(Clone, Debug, Default)]
pub struct ChunkOutput {
    /// Compressed blocks in cut order.
    pub blocks: Vec<OutBlock>,
    /// Plain bytes (uncompressed output).
    pub raw: Vec<u8>,
}

/// What the writer emits: BGZF-compressed VCF text (`-Oz`, `-Ob`),
/// BGZF level-0 (`-Ou`) or plain text (`-Ov`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    /// Plain VCF text.
    Plain,
    /// BGZF blocks at the configured deflate level.
    Bgzf,
    /// BGZF level-0 blocks (uncompressed BCF's container, htslib's `-Ou`).
    BgzfLevel0,
}

/// Append-only block builder.
pub struct BlockWriter {
    compressed: bool,
    deflater: Option<Deflater>,
    pending: Vec<u8>,
    records: Vec<IndexRecord>,
    blocks: Vec<OutBlock>,
    raw: Vec<u8>,
}

impl BlockWriter {
    /// Creates a writer; `compressed = false` degrades it to a plain byte
    /// buffer with the same interface (`-O v/u`).
    ///
    /// # Errors
    ///
    /// [`IoError::Io`] when the deflate level is invalid.
    pub fn new(level: u8, compressed: bool) -> Result<Self, IoError> {
        Ok(Self {
            compressed,
            deflater: if compressed {
                Some(Deflater::new(i32::from(level)).map_err(|e| IoError::Io(e.to_string()))?)
            } else {
                None
            },
            pending: Vec::with_capacity(BGZF_BLOCK_SIZE + (1 << 14)),
            records: Vec::new(),
            blocks: Vec::new(),
            raw: Vec::new(),
        })
    }

    /// A writer for an explicit output format (BCF's `-Ou` wants level-0
    /// BGZF blocks, which is still "compressed" framing).
    ///
    /// # Errors
    ///
    /// [`IoError::Io`] from the deflater.
    pub fn for_format(format: OutputFormat, level: u8) -> Result<Self, IoError> {
        match format {
            OutputFormat::Plain => Self::new(level, false),
            OutputFormat::Bgzf => Self::new(level, true),
            OutputFormat::BgzfLevel0 => Self::new(0, true),
        }
    }

    /// Appends one record and cuts a block whenever the pending buffer
    /// reached [`BGZF_BLOCK_SIZE`] bytes.
    ///
    /// # Errors
    ///
    /// [`IoError`] from the deflater.
    pub fn push_record(&mut self, bytes: &[u8], key: IndexKey) -> Result<(), IoError> {
        if !self.compressed {
            self.raw.extend_from_slice(bytes);
            return Ok(());
        }
        let offset = u32::try_from(self.pending.len())
            .map_err(|_| IoError::Io("pending block exceeds u32::MAX bytes".into()))?;
        self.records.push(IndexRecord { key, offset_in_block: offset });
        self.pending.extend_from_slice(bytes);
        while self.pending.len() >= BGZF_BLOCK_SIZE {
            self.cut()?;
        }
        Ok(())
    }

    /// Closes the chunk: cuts the remaining pending bytes into a final block
    /// and hands the output over.
    ///
    /// # Errors
    ///
    /// [`IoError`] from the deflater.
    pub fn finish_chunk(&mut self) -> Result<ChunkOutput, IoError> {
        self.cut()?;
        Ok(ChunkOutput {
            blocks: std::mem::take(&mut self.blocks),
            raw: std::mem::take(&mut self.raw),
        })
    }

    /// Whether nothing has been appended since the last `finish_chunk`.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pending.is_empty() && self.blocks.is_empty() && self.raw.is_empty()
    }

    /// Deflates one block of at most [`BGZF_BLOCK_SIZE`] bytes. Records that
    /// start before the cut stay with this block; later starts are re-based
    /// into the continuation block.
    fn cut(&mut self) -> Result<(), IoError> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let Some(deflater) = &mut self.deflater else {
            return Ok(());
        };
        let cut = self.pending.len().min(BGZF_BLOCK_SIZE);
        let cut32 = u32::try_from(cut).map_err(|_| IoError::Io("block cut exceeds u32".into()))?;
        let mut compressed = Vec::with_capacity(cut / 2 + 64);
        deflater
            .deflate_block(&self.pending[..cut], &mut compressed)
            .map_err(|e| IoError::Io(format!("bgzf deflate: {e}")))?;
        let boundary = self
            .records
            .partition_point(|record| record.offset_in_block < cut32);
        let mut records = std::mem::take(&mut self.records);
        let mut continuation = records.split_off(boundary);
        for record in &mut continuation {
            record.offset_in_block -= cut32;
        }
        self.records = continuation;
        self.blocks.push(OutBlock {
            compressed,
            uncompressed_len: cut32,
            records,
        });
        self.pending.drain(..cut);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgzf::inflate::Inflater;

    fn key(start1: u64) -> IndexKey {
        IndexKey { contig: 0, start1, ref_len: 1 }
    }

    /// Inflates every block back and checks the bytes and record anchors.
    fn round_trip(out: &ChunkOutput, expected: &[u8]) {
        let mut inflater = Inflater::new();
        let mut text = Vec::new();
        for block in &out.blocks {
            let mut scratch = Vec::new();
            inflater
                .inflate_block(&block.compressed, &mut scratch)
                .expect("invariant: block inflates");
            assert_eq!(scratch.len(), block.uncompressed_len as usize);
            text.extend_from_slice(&scratch);
        }
        assert_eq!(text, expected);
    }

    #[test]
    fn cuts_at_the_member_limit_and_anchors_records() {
        let mut writer = BlockWriter::new(6, true).expect("invariant: level 6");
        let a = vec![b'a'; 40_000];
        let b = vec![b'b'; 40_000];
        writer.push_record(&a, key(1)).expect("invariant: push");
        writer.push_record(&b, key(2)).expect("invariant: push");
        let out = writer.finish_chunk().expect("invariant: finish");
        // 80,000 bytes: one cut at the member limit, the rest in a second block.
        assert_eq!(out.blocks.len(), 2);
        assert_eq!(out.blocks[0].uncompressed_len, BGZF_BLOCK_SIZE as u32);
        assert_eq!(out.blocks[1].uncompressed_len, 80_000 - BGZF_BLOCK_SIZE as u32);
        // Both records start before the cut, so both anchor in block 0.
        let anchors: Vec<u32> = out.blocks[0].records.iter().map(|r| r.offset_in_block).collect();
        assert_eq!(anchors, vec![0, 40_000]);
        assert!(out.blocks[1].records.is_empty());
        let mut expected = a.clone();
        expected.extend_from_slice(&b);
        round_trip(&out, &expected);
    }

    #[test]
    fn a_record_wider_than_a_block_spans_blocks() {
        let mut writer = BlockWriter::new(6, true).expect("invariant: level 6");
        let wide = vec![b'x'; BGZF_BLOCK_SIZE + 1_000];
        writer.push_record(&wide, key(1)).expect("invariant: push");
        let out = writer.finish_chunk().expect("invariant: finish");
        assert_eq!(out.blocks.len(), 2);
        assert_eq!(out.blocks[0].records.len(), 1, "the record anchors in the first block");
        assert_eq!(out.blocks[0].records[0].offset_in_block, 0);
        assert!(out.blocks[1].records.is_empty());
        round_trip(&out, &wide);
    }

    #[test]
    fn records_after_a_cut_are_rebased() {
        let mut writer = BlockWriter::new(6, true).expect("invariant: level 6");
        let wide = vec![b'x'; BGZF_BLOCK_SIZE + 1_000];
        writer.push_record(&wide, key(1)).expect("invariant: push");
        writer.push_record(b"tail", key(2)).expect("invariant: push");
        let out = writer.finish_chunk().expect("invariant: finish");
        // `tail` starts after the cut → it is re-based into block 1.
        assert_eq!(out.blocks[1].records.len(), 1);
        assert_eq!(out.blocks[1].records[0].offset_in_block, 1_000);
        let mut expected = wide.clone();
        expected.extend_from_slice(b"tail");
        round_trip(&out, &expected);
    }

    #[test]
    fn plain_path_keeps_bytes_and_no_blocks() {
        let mut writer = BlockWriter::new(6, false).expect("invariant: plain");
        writer.push_record(b"record\n", key(1)).expect("invariant: push");
        let out = writer.finish_chunk().expect("invariant: finish");
        assert!(out.blocks.is_empty());
        assert_eq!(out.raw, b"record\n");
    }
}
