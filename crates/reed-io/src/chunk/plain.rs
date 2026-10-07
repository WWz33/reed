//! Plain VCF and gzip VCF chunk source (03-architecture.md §6). Plain text
//! is read in fixed 1 MiB pieces; gzip is inflated serially and then cut
//! the same way. Both produce `ChunkBytes { plain: true }` that
//! `inflate_chunk` passes through without deflation.

use std::io::{BufRead, Read};
use std::sync::Arc;

use reed_core::ChunkId;

use crate::chunk::{ChunkBytes, CHUNK_TARGET_COMPRESSED};
use crate::error::IoError;
use crate::header::HeaderFacts;

/// Fixed text-byte target per chunk (plain input is already text).
const PLAIN_CHUNK_BYTES: usize = CHUNK_TARGET_COMPRESSED as usize;

/// Reads a plain (uncompressed) VCF or a gzip VCF as plain-text chunks.
pub struct PlainChunkSource {
    reader: Box<dyn BufRead + Send>,
    next_chunk_id: u64,
    header_done: bool,
    first_offset: u64,
    buf: Vec<u8>,
}

impl PlainChunkSource {
    /// Opens and scans the header.
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when no `#CHROM` line is found.
    pub fn open_vcf(reader: impl Read + Send + 'static) -> Result<(Self, HeaderFacts), IoError> {
        Self::open_inner(Box::new(std::io::BufReader::with_capacity(1 << 20, reader)))
    }

    /// Opens a gzip stream (serial inflate) and scans the header.
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when no `#CHROM` line is found.
    pub fn open_gzip(reader: impl Read + Send + 'static) -> Result<(Self, HeaderFacts), IoError> {
        let gz = flate2::read::GzDecoder::new(reader);
        Self::open_inner(Box::new(std::io::BufReader::with_capacity(1 << 20, gz)))
    }

    fn open_inner(mut reader: Box<dyn BufRead + Send>) -> Result<(Self, HeaderFacts), IoError> {
        let mut header_text = Vec::new();
        loop {
            let mut line = Vec::new();
            let n = reader.read_until(b'\n', &mut line)?;
            if n == 0 {
                return Err(IoError::Header("no #CHROM line found in the stream".into()));
            }
            header_text.extend_from_slice(&line);
            if line.starts_with(b"#CHROM") {
                break;
            }
        }
        let facts = HeaderFacts::parse(&header_text)?;
        Ok((
            Self {
                reader,
                next_chunk_id: 0,
                header_done: true,
                first_offset: 0,
                buf: Vec::with_capacity(PLAIN_CHUNK_BYTES + (1 << 14)),
            },
            facts,
        ))
    }
}

impl crate::chunk::ChunkSource for PlainChunkSource {
    fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, IoError> {
        if !self.header_done {
            return Err(IoError::Header("header not scanned".into()));
        }
        self.buf.clear();
        while self.buf.len() < PLAIN_CHUNK_BYTES {
            let available = self.reader.fill_buf()?;
            if available.is_empty() {
                break;
            }
            let take = available.len().min(PLAIN_CHUNK_BYTES - self.buf.len());
            self.buf.extend_from_slice(&available[..take]);
            self.reader.consume(take);
        }
        if self.buf.is_empty() {
            return Ok(None);
        }
        let id = ChunkId(self.next_chunk_id);
        self.next_chunk_id += 1;
        let offset = std::mem::take(&mut self.first_offset);
        Ok(Some(ChunkBytes {
            id,
            compressed: Arc::from(self.buf.as_slice()),
            first_in_block_offset: offset,
            plain: true,
        }))
    }
}
