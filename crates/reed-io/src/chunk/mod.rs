//! Block-aligned parallel parsing inputs (`docs/spec/architecture.md` §6.2).
//!
//! A chunk is a run of consecutive BGZF blocks cut at a block boundary once
//! its compressed size reaches [`CHUNK_TARGET_COMPRESSED`]. Chunk boundaries
//! depend only on the input bytes — never on thread counts — which is what
//! makes reed's output deterministic (§2 principle 5).

pub mod bgzf;
pub mod framing;
pub mod plain;
pub mod region;

use std::ops::Range;
use std::sync::Arc;

use reed_core::ChunkId;

use crate::error::IoError;

/// Compressed-byte target per chunk. `// determinism-relevant`: changing it
/// changes chunk boundaries, so output could shift at the last float digit;
/// bump the version and record an ADR before touching it
/// (`docs/spec/architecture.md` §6.2).
pub const CHUNK_TARGET_COMPRESSED: u64 = 1 << 20; // 1 MiB

/// A chunk's compressed bytes. Owned when read from a stream; a range of a
/// memory-mapped input otherwise, in which case the reader copies nothing
/// and issues no read syscall (doc_v4 03 §3.4 step 3).
#[derive(Clone, Debug)]
pub enum ChunkPayload {
    /// Bytes read into a heap buffer (files without mmap, stdin, plain
    /// input, region slices).
    Owned(Arc<[u8]>),
    /// A byte range of the mapped input file, shared by every chunk.
    Mapped {
        /// The mapping (one per run).
        map: Arc<memmap2::Mmap>,
        /// The chunk's range within the mapping.
        range: Range<usize>,
    },
}

impl ChunkPayload {
    /// The chunk's bytes.
    #[must_use]
    pub fn as_slice(&self) -> &[u8] {
        match self {
            Self::Owned(bytes) => bytes,
            Self::Mapped { map, range } => &map[range.clone()],
        }
    }

    /// The chunk's compressed length.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Owned(bytes) => bytes.len(),
            Self::Mapped { range, .. } => range.len(),
        }
    }

    /// Whether the chunk has no bytes.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl std::ops::Deref for ChunkPayload {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        self.as_slice()
    }
}

impl AsRef<[u8]> for ChunkPayload {
    fn as_ref(&self) -> &[u8] {
        self.as_slice()
    }
}

impl From<Arc<[u8]>> for ChunkPayload {
    fn from(bytes: Arc<[u8]>) -> Self {
        Self::Owned(bytes)
    }
}

impl From<Vec<u8>> for ChunkPayload {
    fn from(bytes: Vec<u8>) -> Self {
        Self::Owned(Arc::from(bytes.into_boxed_slice()))
    }
}

impl From<Box<[u8]>> for ChunkPayload {
    fn from(bytes: Box<[u8]>) -> Self {
        Self::Owned(Arc::from(bytes))
    }
}

/// One chunk of compressed input, ready for a worker to inflate.
#[derive(Clone, Debug)]
pub struct ChunkBytes {
    /// Input-order index; the only basis of reduction order.
    pub id: ChunkId,
    /// Whole BGZF blocks, concatenated.
    pub compressed: ChunkPayload,
    /// Inflated-text offset where records start (nonzero only for the chunk
    /// that carries the end of the header).
    pub first_in_block_offset: u64,
    /// The chunk holds plain VCF text, not BGZF blocks (uncompressed input).
    pub plain: bool,
    /// Compressed offset of the chunk's first block within the input
    /// (parallel indexing needs it to rebuild virtual offsets).
    pub compressed_start: u64,
    /// Per-block `(compressed length, inflated length)` pairs, in order
    /// (empty for plain chunks). Walks the chunk's inflated text back to
    /// `(compressed offset, in-block offset)` without inflating twice.
    pub blocks: Vec<(u32, u32)>,
    /// Inflated-text prefix length to keep. Set only by the region source on
    /// the chunk that carries an index range's end: the index's end virtual
    /// offset may fall mid-block, so the last block read can hold records
    /// past the region, and the engine must drop them (they otherwise look
    /// like an unterminated final record — doc_v21 M3). `None` everywhere
    /// else, where the stream ends at a real record boundary.
    pub truncate_to: Option<u64>,
}

/// The inflated content of one chunk, split by line ownership.
///
/// `head` (chunks after the first) is the partial line at the start — it
/// belongs to the previous chunk's record and goes back to the reducer for
/// seam stitching. `body` is the byte range of complete lines this chunk
/// owns. `tail` is the trailing partial line after the last `\n` (None when
/// the chunk ends exactly at a line end). All ranges index into `text`.
#[derive(Debug)]
pub struct InflatedChunk<'a> {
    /// Full inflated text.
    pub text: &'a [u8],
    /// Byte range of complete lines this chunk owns (for the first chunk
    /// this starts at the end of the header).
    pub body: Range<usize>,
    /// Partial line at the start (`None` for the first chunk).
    pub head: Option<Range<usize>>,
    /// Partial line after the last `\n` (`None` when empty).
    pub tail: Option<Range<usize>>,
}

impl InflatedChunk<'_> {
    /// The complete lines this chunk owns, as a byte slice.
    #[must_use]
    pub fn body_bytes<'a>(&self, text: &'a [u8]) -> &'a [u8] {
        &text[self.body.clone()]
    }

    /// The head fragment, as a byte slice.
    #[must_use]
    pub fn head_bytes<'a>(&self, text: &'a [u8]) -> Option<&'a [u8]> {
        self.head.as_ref().map(|r| &text[r.clone()])
    }

    /// The tail fragment, as a byte slice.
    #[must_use]
    pub fn tail_bytes<'a>(&self, text: &'a [u8]) -> Option<&'a [u8]> {
        self.tail.as_ref().map(|r| &text[r.clone()])
    }
}

/// Inflates every block of `chunk` into `scratch` and splits the text into
/// head/body/tail line ownership.
///
/// The chunk's `blocks` metadata (built by the reader) drives the walk: no
/// second header parse, and `scratch` is grown once to the chunk's total
/// inflated size — every block inflates straight into its exact-size window,
/// so after the first chunk a worker does no zero-fill at all (A6-6).
///
/// # Errors
///
/// [`IoError::MalformedBgzf`] when the metadata disagrees with the payload
/// or a block is corrupt.
pub fn inflate_chunk<'a>(
    chunk: &ChunkBytes,
    inflater: &mut crate::bgzf::inflate::Inflater,
    scratch: &'a mut Vec<u8>,
) -> Result<InflatedChunk<'a>, IoError> {
    if chunk.plain {
        scratch.clear();
        // Uncompressed input: the "compressed" bytes are the text itself.
        scratch.extend_from_slice(&chunk.compressed);
        let text: &'a [u8] = truncate(chunk, scratch.as_slice());
        return split_text(chunk, text);
    }
    let total_inflated: usize = chunk.blocks.iter().map(|&(_, ilen)| ilen as usize).sum();
    let total_compressed: usize = chunk.blocks.iter().map(|&(clen, _)| clen as usize).sum();
    if total_compressed != chunk.compressed.len() {
        return Err(IoError::MalformedBgzf {
            locus: None,
            reason: "chunk block table disagrees with the payload length",
        });
    }
    if scratch.len() < total_inflated {
        scratch.resize(total_inflated, 0);
    }
    let mut compressed: &[u8] = chunk.compressed.as_ref();
    let mut pos = 0_usize;
    for &(clen, ilen) in &chunk.blocks {
        let block = &compressed[..clen as usize];
        compressed = &compressed[clen as usize..];
        if ilen == 0 {
            continue; // Empty member (EOF marker or concatenation join).
        }
        pos += inflater.inflate_window(block, &mut scratch[pos..pos + ilen as usize])?;
    }
    let text: &'a [u8] = truncate(chunk, &scratch[..pos]);
    split_text(chunk, text)
}

/// The chunk's inflated text, cut to `truncate_to` when the region source set
/// it. The cut lands on a record boundary (the index's end virtual offset), so
/// the shortened text never ends mid-record (doc_v21 M3).
fn truncate<'a>(chunk: &ChunkBytes, text: &'a [u8]) -> &'a [u8] {
    match chunk.truncate_to {
        Some(limit) => {
            let len = usize::try_from(limit).unwrap_or(usize::MAX).min(text.len());
            &text[..len]
        }
        None => text,
    }
}

/// Splits already-inflated text into head/body/tail line ownership.
///
/// # Errors
///
/// [`IoError::MalformedBgzf`] when the chunk's start offset — an index
/// virtual offset's low 16 bits (region.rs) — is beyond the inflated text:
/// a corrupt or hand-edited .csi/.tbi can name one, and the slice used to
/// panic instead of surfacing the error. A within-bounds offset keeps the
/// documented behavior (records before it are dropped) (doc_v20 S40).
fn split_text<'a>(chunk: &ChunkBytes, text: &'a [u8]) -> Result<InflatedChunk<'a>, IoError> {
    let scan_start = if chunk.first_in_block_offset > 0 {
        usize::try_from(chunk.first_in_block_offset).unwrap_or(0)
    } else {
        0
    };
    if scan_start > text.len() {
        return Err(IoError::MalformedBgzf {
            locus: None,
            reason: "chunk start offset is beyond the inflated text (corrupt index virtual offset)",
        });
    }

    if scan_start > 0 || chunk.id.get() == 0 {
        return Ok(InflatedChunk {
            text,
            body: owned_body_range(text, scan_start),
            head: None,
            tail: tail_range(text),
        });
    }

    let Some(first_nl) = memchr::memchr(b'\n', &text[scan_start..]) else {
        return Ok(InflatedChunk {
            text,
            body: scan_start..scan_start,
            head: Some(scan_start..text.len()),
            tail: None,
        });
    };
    let head_end = scan_start + first_nl + 1;
    Ok(InflatedChunk {
        text,
        body: owned_body_range(text, head_end),
        head: Some(scan_start..head_end),
        tail: tail_range(text),
    })
}

/// Range of complete lines: from `from` to one past the last `\n`.
fn owned_body_range(text: &[u8], from: usize) -> Range<usize> {
    match memchr::memrchr(b'\n', &text[from..]) {
        Some(rel) => from..from + rel + 1,
        None => from..from,
    }
}

/// Tail fragment range (after the last `\n`), or None when empty.
fn tail_range(text: &[u8]) -> Option<Range<usize>> {
    match memchr::memrchr(b'\n', text) {
        Some(rel) if rel + 1 < text.len() => Some(rel + 1..text.len()),
        Some(_) => None,
        None => (!text.is_empty()).then_some(0..text.len()),
    }
}

/// The chunk stream the engine consumes (`docs/spec/architecture.md` §7.1).
pub trait ChunkSource: Send {
    /// The next chunk, or `None` at end of stream.
    ///
    /// # Errors
    ///
    /// Transport failures (abort is handled by the engine, not the source).
    fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, IoError>;
}

/// A BGZF file's chunk source; construction performs the header scan, so the
/// first chunk starts exactly at the first record.
pub struct BgzfChunkSource<R: std::io::Read + Send> {
    reader: bgzf::BgzfChunkReader<R>,
}

impl<R: std::io::Read + Send> BgzfChunkSource<R> {
    /// Opens the source and scans the header.
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when the header is missing or unparsable.
    pub fn open(inner: R) -> Result<(Self, crate::header::HeaderFacts), IoError> {
        let mut reader = bgzf::BgzfChunkReader::new(inner);
        let scan = reader.scan_header(&mut crate::bgzf::inflate::Inflater::new())?;
        let facts = crate::header::HeaderFacts::parse(&scan.text)?;
        Ok((Self { reader }, facts))
    }
}

impl<R: std::io::Read + Send> ChunkSource for BgzfChunkSource<R> {
    fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, IoError> {
        self.reader.next_chunk()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chunk(first_in_block_offset: u64, id: u64) -> ChunkBytes {
        ChunkBytes {
            id: ChunkId(id),
            compressed: ChunkPayload::from(Vec::new()),
            first_in_block_offset,
            plain: true,
            compressed_start: 0,
            blocks: Vec::new(),
            truncate_to: None,
        }
    }

    /// doc_v20 S40: `first_in_block_offset` comes straight from an index
    /// file's chunk-start virtual offset (the low 16 bits, region.rs), so a
    /// corrupt or hand-edited .csi/.tbi can name an offset past the end of
    /// the inflated block. That used to panic in the `&text[scan_start..]`
    /// slice ("range start index 64 out of range for slice of length 5")
    /// instead of surfacing as a typed error.
    #[test]
    fn s40_start_offset_beyond_the_text_is_a_typed_error() {
        let text = b"ab\ncd";
        for (offset, id) in [(6, 0), (64, 0), (6, 1), (64, 1)] {
            let err = split_text(&chunk(offset, id), text)
                .expect_err("an offset beyond the inflated text must be MalformedBgzf");
            assert!(
                matches!(err, IoError::MalformedBgzf { .. }),
                "offset {offset} id {id}: {err:?}"
            );
        }
    }

    /// Control (doc_v20 S40): a within-bounds offset keeps the documented
    /// behavior — the bytes before it are dropped, the rest is split by line
    /// ownership.
    #[test]
    fn s40_start_offset_within_the_text_still_splits() {
        let view = split_text(&chunk(3, 0), b"ab\ncd").expect("invariant: within bounds");
        assert_eq!(view.body, 3..3);
        assert_eq!(&view.text[view.tail.expect("invariant: tail")], b"cd");
        assert!(view.head.is_none());
    }
}
