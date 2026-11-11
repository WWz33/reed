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

use crate::bgzf::address::VirtualOffset;
use crate::bgzf::block::EOF_BLOCK;
use crate::bgzf::deflate::Deflater;
use crate::error::IoError;
use crate::write::atomic::AtomicFile;
use crate::write::block_writer::{BGZF_BLOCK_SIZE, ChunkOutput, IndexRecord};
use crate::write::index::Indexer;

/// A record whose index chunk is still open: its key plus the virtual offset
/// where it starts.
type PendingRecord = (IndexRecord, u64);

/// Narrows a block- or record-relative offset to the low 16 bits a virtual
/// offset carries, with the reason a violation is impossible in well-formed
/// BGZF: a block's inflated payload is at most 65,536 bytes.
fn in_block(offset: u32, what: &'static str) -> Result<u16, IoError> {
    u16::try_from(offset).map_err(|_| {
        IoError::Io(format!(
            "{what} ({offset}) does not fit a virtual offset's 16-bit half"
        ))
    })
}

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
    last_block_end: VirtualOffset,
    /// Deflate time observed on the chunks handed in (A6-4), reported to the
    /// engine through `OrderedSink::deflate_time`.
    deflate_nanos: u64,
    /// The atomic owner of the output file when the path was opened here
    /// (doc_v26 V26-08). `None` for a caller-supplied writer (stdout, a
    /// temp spill, a test buffer), where there is no file to publish.
    atomic: Option<AtomicFile>,
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
            last_block_end: VirtualOffset::from_raw(0),
            deflate_nanos: 0,
            atomic: None,
        }
    }

    /// Makes this sink the owner of `atomic`: [`Self::finish`] publishes it.
    ///
    /// Used when the caller opened the file itself (the filter clones the
    /// handle for its `BufWriter`).
    #[must_use]
    pub fn adopting(mut self, atomic: AtomicFile) -> Self {
        self.atomic = Some(atomic);
        self
    }

    /// Creates the sink for a **file** output: the bytes go to a sibling
    /// temporary and the file only appears on [`Self::finish`]
    /// (doc_v26 V26-08).
    ///
    /// The index sidecar is written to its own temporary and published in
    /// `finish` too, so a failure before `finish` leaves neither a truncated
    /// VCF nor an index that describes a file that never appeared. The two
    /// renames are still individual: an interrupt between them leaves the new
    /// VCF with the old index — the per-file promise documented on
    /// [`crate::write::atomic`].
    ///
    /// # Errors
    ///
    /// I/O failures creating the temporary file.
    pub fn create(
        output: &Path,
        compressed: bool,
        level: u8,
        indexer: Option<Indexer>,
        index_path: Option<PathBuf>,
        contig_names: Arc<[Vec<u8>]>,
    ) -> Result<Self, IoError> {
        let (atomic, file) = AtomicFile::create(output)?;
        let writer: Box<dyn Write + Send> =
            Box::new(std::io::BufWriter::with_capacity(1 << 20, file));
        let mut out = Self::new(writer, compressed, level, indexer, index_path, contig_names);
        out.atomic = Some(atomic);
        Ok(out)
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
        self.deflate_nanos += chunk.deflate_nanos;
        if self.compressed {
            for block in &chunk.blocks {
                self.accept_block(block)?;
            }
        } else {
            self.writer.write_all(&chunk.raw)?;
        }
        Ok(())
    }

    /// Cumulative worker-side deflate time observed (A6-4).
    #[must_use]
    pub const fn deflate_time(&self) -> std::time::Duration {
        std::time::Duration::from_nanos(self.deflate_nanos)
    }

    /// Finishes the stream: EOF marker, flush, the last record's index chunk,
    /// and the index file itself.
    ///
    /// # Errors
    ///
    /// [`IoError`] from the writer or the indexer.
    pub fn finish(self) -> Result<(), IoError> {
        self.finish_pending()?.commit()
    }

    /// Finishes writing but does **not** publish the file or its index.
    ///
    /// The caller owns the final integrity check (for filter, `job.finish`
    /// rejects a BCF truncated at its last record) and calls
    /// [`PendingOrderedOutput::commit`] only after it passes. Dropping the
    /// pending value removes the temporary and leaves the old output and
    /// index untouched (doc_v28 V28-06).
    ///
    /// # Errors
    ///
    /// [`IoError`] from the writer or the indexer.
    pub fn finish_pending(mut self) -> Result<PendingOrderedOutput, IoError> {
        if self.compressed {
            self.writer.write_all(&EOF_BLOCK)?;
        }
        self.writer.flush()?;
        if let Some(pending) = self.pending_record.take() {
            // The last record's chunk runs to the end of the last block.
            self.close_index_record(pending, self.last_block_end.raw())?;
        }
        Ok(PendingOrderedOutput {
            atomic: self.atomic.take(),
            indexer: self.indexer.take(),
            index_path: self.index_path.take(),
        })
    }

    fn accept_block(
        &mut self,
        block: &crate::write::block_writer::OutBlock,
    ) -> Result<(), IoError> {
        let coffset_before = self.coffset;
        self.writer.write_all(&block.compressed)?;
        self.coffset += block.compressed.len() as u64;
        self.last_block_end = VirtualOffset::new(
            coffset_before,
            in_block(block.uncompressed_len, "a block's inflated length")?,
        );
        if self.indexer.is_some() {
            for record in &block.records {
                let voff_start = VirtualOffset::new(
                    coffset_before,
                    in_block(record.offset_in_block, "a record's in-block offset")?,
                )
                .raw();
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

/// A finished output whose publication is still pending an external check.
///
/// Dropping it removes the data temporary and writes no index. `commit`
/// publishes the data file first, then its index — a new VCF with a stale
/// index is a visible mismatch; a new index over an old VCF would answer
/// region queries with the wrong records silently. The two renames are
/// individual, not a group transaction.
pub struct PendingOrderedOutput {
    atomic: Option<AtomicFile>,
    indexer: Option<Indexer>,
    index_path: Option<PathBuf>,
}

impl PendingOrderedOutput {
    /// Publishes the data file and, when requested, its index.
    ///
    /// # Errors
    ///
    /// The commit or index-write failure.
    pub fn commit(self) -> Result<(), IoError> {
        if let Some(atomic) = self.atomic {
            atomic.commit()?;
        }
        if let (Some(indexer), Some(path)) = (self.indexer, self.index_path) {
            indexer.write(&path)?;
        }
        Ok(())
    }
}
