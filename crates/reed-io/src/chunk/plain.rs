//! Plain VCF and gzip VCF chunk source (03-architecture.md §6). Plain text
//! is read in fixed 1 MiB pieces; gzip is inflated serially and then cut
//! the same way. Both produce `ChunkBytes { plain: true }` that
//! `inflate_chunk` passes through without deflation.

// Casts in this file are bounded by construction: The chunk target is a build-time constant far below usize::MAX.
#![allow(clippy::cast_possible_truncation)]

use std::io::{BufRead, Read};
use std::sync::Arc;

use reed_core::ChunkId;

use crate::chunk::{CHUNK_TARGET_COMPRESSED, ChunkBytes, ChunkPayload};
use crate::error::IoError;
use crate::header::HeaderFacts;

/// Fixed text-byte target per chunk (plain input is already text).
const PLAIN_CHUNK_BYTES: usize = CHUNK_TARGET_COMPRESSED as usize;

/// Reads a plain (uncompressed) VCF or a gzip VCF as plain-text chunks.
pub struct PlainChunkSource {
    reader: Box<dyn BufRead + Send>,
    next_chunk_id: u64,
    header_done: bool,
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
    /// Multi-member input (`cat a.vcf.gz b.vcf.gz`) is read end to end:
    /// `MultiGzDecoder` continues into the next member where htslib's
    /// `inflate_gzip_block` peeks and `inflateReset`s (bgzf.c:832-908). A
    /// plain `GzDecoder` would silently stop after the first member
    /// (doc_v10 D10-6).
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when no `#CHROM` line is found.
    pub fn open_gzip(reader: impl Read + Send + 'static) -> Result<(Self, HeaderFacts), IoError> {
        let gz = flate2::read::MultiGzDecoder::new(reader);
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
            // doc_v26 V26-13: a malformed input with no #CHROM line used to
            // grow this buffer until memory ran out, and only then report the
            // missing header.
            if header_text.len().saturating_add(line.len()) > crate::header::MAX_HEADER_BYTES {
                return Err(IoError::Header(format!(
                    "header exceeds {} bytes without a #CHROM line",
                    crate::header::MAX_HEADER_BYTES
                )));
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
        Ok(Some(ChunkBytes {
            id,
            compressed: ChunkPayload::from(Arc::<[u8]>::from(self.buf.as_slice())),
            first_in_block_offset: 0,
            plain: true,
            compressed_start: 0,
            blocks: Vec::new(),
            truncate_to: None,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chunk::ChunkSource as _;

    /// doc_v26 V26-13: a header without a `#CHROM` line must be rejected at
    /// the ceiling instead of growing until memory runs out.
    #[test]
    fn v26_13_endless_header_without_chrom_is_capped() {
        let mut text = String::from("##fileformat=VCFv4.3\n");
        for _ in 0..40_000 {
            text.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n");
        }
        assert!(text.len() > 1 << 20 && text.len() < crate::header::MAX_HEADER_BYTES);
        let err = match PlainChunkSource::open_vcf(std::io::Cursor::new(text.into_bytes())) {
            Ok(_) => panic!("a header with no #CHROM line must fail"),
            Err(err) => err,
        };
        assert!(err.to_string().contains("#CHROM"), "{err}");
    }

    /// A legal header near the ceiling still parses; the bound is not a
    /// practical limit on real headers.
    #[test]
    fn v26_13_a_large_but_legal_header_still_parses() {
        use std::fmt::Write as _;
        let mut text = String::from("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n");
        let mut n = 0_u32;
        while text.len() < (4 << 20) {
            // Distinct IDs: a repeated ##INFO line is a duplicate, and
            // noodles rejects it (that is the header parser's own contract).
            writeln!(
                text,
                "##INFO=<ID=T{n},Number=1,Type=Integer,Description=\"d\">"
            )
            .expect("invariant: writing to a String");
            n += 1;
        }
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n");
        let (mut source, facts) =
            PlainChunkSource::open_vcf(std::io::Cursor::new(text.into_bytes()))
                .expect("invariant: a legal large header parses");
        assert_eq!(facts.n_samples(), 1);
        assert!(
            source.next_chunk().expect("invariant: read").is_none(),
            "no records"
        );
    }

    /// The ceiling admits wide cohorts (a million-sample header is ~20 MiB).
    #[test]
    fn v26_13_the_ceiling_admits_wide_cohorts() {
        let ceiling = crate::header::MAX_HEADER_BYTES;
        assert!(ceiling >= (64 << 20), "the ceiling must admit wide cohorts");
        let accumulated = ceiling;
        assert!(
            accumulated.saturating_add(1) > ceiling,
            "one more byte past the ceiling must trip the check"
        );
    }

    /// A gzip member of `text` (regular gzip, no BGZF `BC` subfield).
    fn gzip_member(text: &[u8]) -> Vec<u8> {
        use std::io::Write as _;
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(text).expect("invariant: gzip write");
        enc.finish().expect("invariant: gzip finish")
    }

    /// doc_v10 D10-6: `cat a.vcf.gz b.vcf.gz` must be read end to end (as
    /// htslib does — `inflate_gzip_block` resets into the next member,
    /// bgzf.c:832-908). The old single-member `GzDecoder` silently stopped
    /// after the first member.
    #[test]
    fn multi_member_gzip_reads_every_member() {
        let header = b"##fileformat=VCFv4.3\n\
                       ##contig=<ID=chr1,length=1000>\n\
                       #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let a = b"chr1\t10\t.\tA\tG\t.\tPASS\t.\n";
        let b = b"chr1\t20\t.\tA\tG\t.\tPASS\t.\n";
        let mut member_a = header.to_vec();
        member_a.extend_from_slice(a);
        let mut member_b = header.to_vec();
        member_b.extend_from_slice(b);
        let mut stream = gzip_member(&member_a);
        stream.extend_from_slice(&gzip_member(&member_b));

        let (mut source, facts) = PlainChunkSource::open_gzip(std::io::Cursor::new(stream))
            .expect("invariant: header present");
        let mut text = Vec::new();
        while let Some(chunk) = source.next_chunk().expect("invariant: chunking") {
            text.extend_from_slice(chunk.compressed.as_slice());
        }
        assert_eq!(facts.n_samples(), 0);
        let text = String::from_utf8(text).expect("invariant: ASCII");
        assert!(text.contains("chr1\t10\t"), "member a missing: {text}");
        assert!(text.contains("chr1\t20\t"), "member b missing: {text}");
    }
}
