//! Reducer-side ordered output (`docs/rounds/v2/03-architecture.md` §5).
//!
//! Appends worker blocks in chunk order, tracks the compressed offset, turns
//! each record's block-local anchor into a virtual offset, and feeds the
//! indexer while writing — no second read for `--write-index` (§6.6).
//!
//! Index chunks follow the *record* stream, not the block stream: a record's
//! chunk ends where the next record starts, even when that sits in a later
//! block (records may span blocks).

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::bgzf::block::EOF_BLOCK;
use crate::bgzf::deflate::Deflater;
use crate::error::IoError;
use crate::write::block_writer::{BGZF_BLOCK_SIZE, ChunkOutput, IndexRecord};
use crate::write::index::Indexer;

/// A record whose index chunk is still open: its key plus the virtual offset
/// where it starts.
type PendingRecord = (IndexRecord, u64);

/// Ordered sink for the filter's and `reed index`'s output streams.
pub struct OrderedOutput {
    writer: Box<dyn Write + Send>,
    compressed: bool,
    level: u8,
    indexer: Option<Indexer>,
    index_path: Option<PathBuf>,
    contig_names: Arc<[Vec<u8>]>,
    /// Compressed bytes written so far.
    coffset: u64,
    /// The record whose index chunk is still open.
    pending_record: Option<PendingRecord>,
    /// Virtual offset just past the last written block.
    last_block_end: u64,
}

impl OrderedOutput {
    /// Creates the sink. `indexer`/`index_path` are both present or both
    /// absent (`-W` with a `.vcf.gz` output).
    #[must_use]
    pub fn new(
        writer: Box<dyn Write + Send>,
        compressed: bool,
        level: u8,
        indexer: Option<Indexer>,
        index_path: Option<PathBuf>,
        contig_names: Arc<[Vec<u8>]>,
    ) -> Self {
        Self {
            writer,
            compressed,
            level,
            indexer,
            index_path,
            contig_names,
            coffset: 0,
            pending_record: None,
            last_block_end: 0,
        }
    }

    /// Writes the (already rendered) header bytes as BGZF block(s).
    ///
    /// # Errors
    ///
    /// [`IoError`] from the deflater or the underlying writer.
    pub fn write_header(&mut self, header: &[u8]) -> Result<(), IoError> {
        if !self.compressed {
            self.writer.write_all(header)?;
            return Ok(());
        }
        let mut deflater =
            Deflater::new(i32::from(self.level)).map_err(|e| IoError::Io(e.to_string()))?;
        for chunk in header.chunks(BGZF_BLOCK_SIZE) {
            let mut block = Vec::with_capacity(chunk.len() / 2 + 64);
            deflater
                .deflate_block(chunk, &mut block)
                .map_err(|e| IoError::Io(format!("bgzf deflate: {e}")))?;
            self.coffset += block.len() as u64;
            self.writer.write_all(&block)?;
        }
        Ok(())
    }

    /// Appends one chunk's output in chunk order.
    ///
    /// # Errors
    ///
    /// [`IoError`] from the writer or the indexer.
    pub fn accept_chunk(&mut self, chunk: &ChunkOutput) -> Result<(), IoError> {
        if self.compressed {
            for block in &chunk.blocks {
                self.accept_block(block)?;
            }
        } else {
            self.writer.write_all(&chunk.raw)?;
        }
        Ok(())
    }

    /// Finishes the stream: EOF marker, flush, the last record's index chunk,
    /// and the index file itself.
    ///
    /// # Errors
    ///
    /// [`IoError`] from the writer or the indexer.
    pub fn finish(mut self) -> Result<(), IoError> {
        if self.compressed {
            self.writer.write_all(&EOF_BLOCK)?;
        }
        self.writer.flush()?;
        if let Some(pending) = self.pending_record.take() {
            // The last record's chunk runs to the end of the last block.
            self.close_index_record(pending, self.last_block_end)?;
        }
        if let (Some(indexer), Some(path)) = (self.indexer.take(), self.index_path.take()) {
            indexer.write(&path)?;
        }
        Ok(())
    }

    fn accept_block(
        &mut self,
        block: &crate::write::block_writer::OutBlock,
    ) -> Result<(), IoError> {
        let coffset_before = self.coffset;
        self.writer.write_all(&block.compressed)?;
        self.coffset += block.compressed.len() as u64;
        self.last_block_end = (coffset_before << 16) | u64::from(block.uncompressed_len);
        if self.indexer.is_some() {
            for record in &block.records {
                let voff_start = (coffset_before << 16) | u64::from(record.offset_in_block);
                if let Some(pending) = self.pending_record.take() {
                    self.close_index_record(pending, voff_start)?;
                }
                self.pending_record = Some((*record, voff_start));
            }
        }
        Ok(())
    }

    /// Closes one record's index chunk `[start, end)` against the indexer.
    fn close_index_record(&mut self, pending: PendingRecord, voff_end: u64) -> Result<(), IoError> {
        let Some(indexer) = &mut self.indexer else {
            return Ok(());
        };
        let (record, voff_start) = pending;
        // The CSI path (the default) keys by contig id and never touches the
        // name, so no per-record String is built (A6-12).
        let name = self
            .contig_names
            .get(record.key.contig.get() as usize)
            .map_or(&[][..], |v| v.as_slice());
        indexer.add(
            name,
            record.key.contig,
            record.key.start1,
            record.key.ref_len,
            voff_start,
            voff_end,
        )
    }

    /// The index path for an output path and flavor.
    #[must_use]
    pub fn index_path_for(output: &Path, flavor: crate::write::index::WriteIndex) -> PathBuf {
        flavor.path_for(output)
    }
}
