//! Block-aligned parallel parsing inputs (`docs/architecture-spec.md` §6.2).
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
/// (`docs/architecture-spec.md` §6.2).
pub const CHUNK_TARGET_COMPRESSED: u64 = 1 << 20; // 1 MiB

/// One chunk of compressed input, ready for a worker to inflate.
#[derive(Clone, Debug)]
pub struct ChunkBytes {
    /// Input-order index; the only basis of reduction order.
    pub id: ChunkId,
    /// Whole BGZF blocks, concatenated.
    pub compressed: Arc<[u8]>,
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

/// Inflates every block of `chunk` into `scratch` (replacing its content)
/// and splits the text into head/body/tail line ownership.
///
/// # Errors
///
/// [`IoError::MalformedBgzf`] on corrupt blocks.
pub fn inflate_chunk<'a>(
    chunk: &ChunkBytes,
    inflater: &mut crate::bgzf::inflate::Inflater,
    scratch: &'a mut Vec<u8>,
) -> Result<InflatedChunk<'a>, IoError> {
    scratch.clear();
    if chunk.plain {
        // Uncompressed input: the "compressed" bytes are the text itself.
        scratch.extend_from_slice(&chunk.compressed);
        return Ok(split_text(chunk, scratch));
    }
    let mut compressed: &[u8] = chunk.compressed.as_ref();
    while !compressed.is_empty() {
        let Some(total) = crate::bgzf::block::parse_block_size(compressed) else {
            return Err(IoError::MalformedBgzf {
                locus: None,
                reason: "chunk contains a non-BGZF block header",
            });
        };
        inflater.inflate_block(&compressed[..total], scratch)?;
        compressed = &compressed[total..];
    }
    Ok(split_text(chunk, scratch))
}

/// Splits already-inflated text into head/body/tail line ownership.
fn split_text<'a>(chunk: &ChunkBytes, scratch: &'a mut Vec<u8>) -> InflatedChunk<'a> {
    let text: &'a [u8] = scratch.as_slice();

    let scan_start = if chunk.first_in_block_offset > 0 {
        usize::try_from(chunk.first_in_block_offset).unwrap_or(0)
    } else {
        0
    };

    if scan_start > 0 || chunk.id.get() == 0 {
        return InflatedChunk {
            text,
            body: owned_body_range(text, scan_start),
            head: None,
            tail: tail_range(text),
        };
    }

    let Some(first_nl) = memchr::memchr(b'\n', &text[scan_start..]) else {
        return InflatedChunk {
            text,
            body: scan_start..scan_start,
            head: Some(scan_start..text.len()),
            tail: None,
        };
    };
    let head_end = scan_start + first_nl + 1;
    InflatedChunk {
        text,
        body: owned_body_range(text, head_end),
        head: Some(scan_start..head_end),
        tail: tail_range(text),
    }
}

/// Range of complete lines: from `from` to one past the last `\n`.
fn owned_body_range(text: &[u8], from: usize) -> Range<usize> {
    match text[from..].iter().rposition(|&b| b == b'\n') {
        Some(rel) => from..from + rel + 1,
        None => from..from,
    }
}

/// Tail fragment range (after the last `\n`), or None when empty.
fn tail_range(text: &[u8]) -> Option<Range<usize>> {
    match text.iter().rposition(|&b| b == b'\n') {
        Some(rel) if rel + 1 < text.len() => Some(rel + 1..text.len()),
        Some(_) => None,
        None if !text.is_empty() => Some(0..text.len()),
        None => None,
    }
}

/// The chunk stream the engine consumes (`docs/architecture-spec.md` §7.1).
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
    pub fn open(
        inner: R,
    ) -> Result<(Self, crate::header::HeaderFacts), IoError> {
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
