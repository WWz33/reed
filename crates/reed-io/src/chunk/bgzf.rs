//! BGZF chunk reader: the serial reader of BAPP (`docs/architecture-spec.md`
//! §6.2). It only parses 18-byte block headers — no inflation — and cuts
//! chunks at block boundaries once the compressed target is reached.

use std::io::Read;
use std::sync::Arc;

use reed_core::ChunkId;

use crate::bgzf::block::{is_empty_member, is_eof_block, parse_block_size, BGZF_HEADER_LEN};
use crate::chunk::{ChunkBytes, CHUNK_TARGET_COMPRESSED};
use crate::error::IoError;

/// Locates the end of the VCF header and hands out chunks from there.
pub struct BgzfChunkReader<R: Read> {
    inner: R,
    chunk_target: u64,
    next_chunk_id: u64,
    /// Compressed bytes of the block where the header ends (chunk 0 starts
    /// here, so it is not read twice).
    pending_block: Option<Vec<u8>>,
    /// Compressed offset of `pending_block`.
    pending_block_offset: u64,
    /// Compressed offset of the next byte to read (bytes consumed so far,
    /// empty members included).
    compressed_offset: u64,
    /// Offset of the first record inside that block's inflated text.
    first_in_block_offset: u64,
    header_done: bool,
    /// Whether the canonical EOF marker was the last block read (htslib
    /// warns when a stream ends without one).
    last_was_eof_marker: bool,
    /// Whether the missing-EOF-marker warning already fired (once per run).
    warned_missing_eof: bool,
}

/// Result of the header scan: the accumulated header text plus where
/// records begin.
#[derive(Debug)]
pub struct HeaderScan {
    /// Full header text including the `#CHROM` line and its newline.
    pub text: Vec<u8>,
    /// Number of sample columns on the `#CHROM` line.
    pub n_samples: usize,
}

impl<R: Read> BgzfChunkReader<R> {
    /// Opens a reader with the canonical chunk target.
    #[must_use]
    pub fn new(inner: R) -> Self {
        Self::with_chunk_target(inner, CHUNK_TARGET_COMPRESSED)
    }

    /// Opens a reader with an explicit chunk target (tests use small values
    /// to force seams; production code uses the determinism-relevant
    /// constant).
    #[must_use]
    pub fn with_chunk_target(inner: R, chunk_target: u64) -> Self {
        Self {
            inner,
            chunk_target,
            next_chunk_id: 0,
            pending_block: None,
            pending_block_offset: 0,
            compressed_offset: 0,
            first_in_block_offset: 0,
            header_done: false,
            last_was_eof_marker: false,
            warned_missing_eof: false,
        }
    }

    /// Reads one BGZF block (header + payload) and its compressed offset;
    /// `None` at the end of the underlying stream.
    ///
    /// Empty members (ISIZE = 0, including the canonical EOF marker) are
    /// skipped — concatenated BGZF files carry one at every join, and the
    /// stream ends only at the underlying EOF (01-review D2). A clean end
    /// without a trailing EOF marker warns once (htslib behavior).
    fn read_block(&mut self) -> Result<Option<(Vec<u8>, u64)>, IoError> {
        loop {
            let block_offset = self.compressed_offset;
            let mut header = [0_u8; BGZF_HEADER_LEN];
            let mut filled = 0;
            while filled < BGZF_HEADER_LEN {
                match self.inner.read(&mut header[filled..]) {
                    Ok(0) => {
                        if filled == 0 {
                            self.warn_if_missing_eof_marker();
                            return Ok(None);
                        }
                        return Err(IoError::MalformedBgzf {
                            locus: None,
                            reason: "stream ends inside a block header",
                        });
                    }
                    Ok(n) => filled += n,
                    Err(err) if err.kind() == std::io::ErrorKind::Interrupted => {}
                    Err(err) => return Err(err.into()),
                }
            }
            let Some(total) = parse_block_size(&header) else {
                return Err(IoError::MalformedBgzf {
                    locus: None,
                    reason: "not a BGZF block header (trailing garbage after the last block?)",
                });
            };
            let mut block = Vec::with_capacity(total);
            block.extend_from_slice(&header);
            block.resize(total, 0);
            self.inner.read_exact(&mut block[BGZF_HEADER_LEN..])?;
            self.compressed_offset += total as u64;
            if is_empty_member(&block) {
                self.last_was_eof_marker = is_eof_block(&block);
                continue;
            }
            self.last_was_eof_marker = false;
            return Ok(Some((block, block_offset)));
        }
    }

    /// Warns once when the stream ended without the canonical EOF marker.
    fn warn_if_missing_eof_marker(&mut self) {
        if !self.warned_missing_eof && self.header_done && !self.last_was_eof_marker {
            self.warned_missing_eof = true;
            eprintln!("reed: warning: BGZF stream ends without an EOF marker");
        }
    }

    /// Reads blocks, inflating as it goes, until the `#CHROM` line ends.
    /// Remembers the block containing that position for chunk 0.
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when the stream ends before a `#CHROM` line.
    pub fn scan_header(
        &mut self,
        inflater: &mut crate::bgzf::inflate::Inflater,
    ) -> Result<HeaderScan, IoError> {
        let mut text = Vec::new();
        loop {
            let Some((block, block_offset)) = self.read_block()? else {
                return Err(IoError::Header("no #CHROM line found in the stream".into()));
            };
            let text_len_before = text.len();
            inflater.inflate_block(&block, &mut text)?;
            if let Some(end) = find_chrom_line_end(&text) {
                self.pending_block = Some(block);
                self.pending_block_offset = block_offset;
                // The offset is *within this block's* inflated text, not in
                // the accumulated header text.
                self.first_in_block_offset = (end - text_len_before) as u64;
                self.header_done = true;
                text.truncate(end);
                let n_samples = count_samples(&text).ok_or_else(|| {
                    IoError::Header("#CHROM line found but has no sample columns".into())
                })?;
                return Ok(HeaderScan { text, n_samples });
            }
        }
    }

    /// Yields the next chunk of compressed blocks, or `None` at end of
    /// stream. Requires [`BgzfChunkReader::scan_header`] first.
    pub fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, IoError> {
        if !self.header_done {
            return Err(IoError::Header("scan_header must run before chunks".into()));
        }
        let (first, compressed_start) = match self.pending_block.take() {
            Some(block) => (block, self.pending_block_offset),
            None => match self.read_block()? {
                Some(pair) => pair,
                None => return Ok(None),
            },
        };
        let mut blocks = vec![block_pair(&first)];
        let mut compressed = first;
        while (compressed.len() as u64) < self.chunk_target {
            match self.read_block()? {
                Some((block, _offset)) => {
                    blocks.push(block_pair(&block));
                    compressed.extend_from_slice(&block);
                }
                None => break,
            }
        }
        let id = ChunkId(self.next_chunk_id);
        self.next_chunk_id += 1;
        Ok(Some(ChunkBytes {
            id,
            compressed: Arc::from(compressed.into_boxed_slice()),
            first_in_block_offset: std::mem::take(&mut self.first_in_block_offset),
            plain: false,
            compressed_start,
            blocks,
        }))
    }
}

/// One block's `(compressed length, inflated length)`.
fn block_pair(block: &[u8]) -> (u32, u32) {
    (
        u32::try_from(block.len()).unwrap_or(u32::MAX),
        crate::bgzf::block::block_inflated_size(block),
    )
}

/// Offset just past the newline ending the `#CHROM` line, scanning text
/// from `from`. The caller re-scans a small overlap, because the line may
/// span block boundaries; matching `\n#CHROM` via memmem stays correct for
/// any overlap >= the needle length.
fn find_chrom_line_end(text: &[u8]) -> Option<usize> {
    let start = if text.starts_with(b"#CHROM") {
        0
    } else {
        // Headers are a handful of blocks; a whole-text search per block is
        // negligible and immune to lines spanning block boundaries.
        memchr::memmem::find(text, b"\n#CHROM")? + 1
    };
    let nl = memchr::memchr(b'\n', &text[start..])? + start;
    Some(nl + 1)
}

/// Sample column count on the `#CHROM` line (tabs after FORMAT).
fn count_samples(header_text: &[u8]) -> Option<usize> {
    let line = header_text
        .split(|&b| b == b'\n')
        .find(|l| l.starts_with(b"#CHROM"))?;
    Some(line.iter().filter(|&&b| b == b'\t').count().saturating_sub(8))
}

#[cfg(test)]
mod tests {
    use super::super::inflate_chunk;
    use super::*;
    use crate::bgzf::deflate::Deflater;

    /// Builds a BGZF stream by cutting the concatenated text into blocks of
    /// `block_data` uncompressed bytes — block boundaries may fall mid-line,
    /// which is exactly what forces seams.
    fn compress_lines(lines: &[&[u8]], block_data: usize) -> Vec<u8> {
        let mut deflater = Deflater::new(6).expect("invariant: level 6 is valid");
        let mut out = Vec::new();
        let mut text: Vec<u8> = Vec::new();
        for line in lines {
            text.extend_from_slice(line);
        }
        for block in text.chunks(block_data) {
            deflater.deflate_block(block, &mut out).expect("invariant: test block fits");
        }
        out.extend_from_slice(&crate::bgzf::block::EOF_BLOCK);
        out
    }

    fn header_text(n_samples: usize) -> Vec<u8> {
        let mut h = String::from("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100000>\n");
        h.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for i in 0..n_samples {
            h.push_str(&format!("\ts{i}"));
        }
        h.push('\n');
        h.into_bytes()
    }

    #[test]
    fn header_scan_counts_samples() {
        let mut text = header_text(3);
        text.extend_from_slice(b"chr1\t10\t.\tA\tG\t1\tPASS\t.\tGT\t0/0\t0/1\t1/1\n");
        let stream = compress_lines(&[&text], 64);
        let mut reader = BgzfChunkReader::with_chunk_target(Cursor(stream), 64);
        let mut inflater = crate::bgzf::inflate::Inflater::new();
        let scan = reader.scan_header(&mut inflater).expect("invariant: header present");
        assert_eq!(scan.n_samples, 3);
        assert!(scan.text.ends_with(b"\ts2\n"));
    }

    struct Cursor(Vec<u8>);
    impl Read for Cursor {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            let n = buf.len().min(self.0.len());
            buf[..n].copy_from_slice(&self.0[..n]);
            self.0.drain(..n);
            Ok(n)
        }
    }

    #[test]
    fn chunks_cover_all_records_with_many_seams() {
        // 100 tiny records, 32-byte blocks, 40-byte chunks: seam-rich.
        let mut text = header_text(2);
        for i in 0..100 {
            text.extend_from_slice(
                format!("chr1\t{}\t.\tA\tG\t1\tPASS\t.\tGT\t0/0\t0/1\n", 100 + i * 7).as_bytes(),
            );
        }
        let stream = compress_lines(&[&text], 32);
        let mut reader = BgzfChunkReader::with_chunk_target(Cursor(stream), 40);
        let mut inflater = crate::bgzf::inflate::Inflater::new();
        let _ = reader.scan_header(&mut inflater).expect("invariant: header present");

        let mut n_record_lines = 0_usize;
        let mut n_seam_bytes = 0_usize;
        // Chained fragments of the record spanning the seam(s): the previous
        // chunk's tail, plus any number of whole-chunk head fragments, until
        // a newline completes the record (spec §6.2 step 4, chained case).
        let mut fragment: Option<Vec<u8>> = None;
        let mut scratch = Vec::new();
        while let Some(chunk) = reader.next_chunk().expect("invariant: chunking") {
            let view = inflate_chunk(&chunk, &mut inflater, &mut scratch).expect("invariant: inflate");
            if let Some(head) = view.head_bytes(view.text) {
                let mut seam = fragment.take().unwrap_or_default();
                n_seam_bytes += head.len();
                seam.extend_from_slice(head);
                if seam.ends_with(b"\n") {
                    assert_eq!(seam.iter().filter(|&&b| b == b'\n').count(), 1, "seam is one line");
                    assert!(seam.starts_with(b"chr1\t"));
                    n_record_lines += 1;
                } else {
                    fragment = Some(seam);
                }
            }
            for line in view.body_bytes(view.text).split(|&b| b == b'\n') {
                if line.is_empty() {
                    continue;
                }
                assert!(line.starts_with(b"chr1\t"));
                n_record_lines += 1;
            }
            fragment = match (fragment, view.tail_bytes(view.text)) {
                (Some(mut acc), Some(tail)) => {
                    acc.extend_from_slice(tail);
                    Some(acc)
                }
                (Some(acc), None) => Some(acc),
                (None, Some(tail)) => Some(tail.to_vec()),
                (None, None) => None,
            };
        }
        // The stream ends with a newline, so no fragment may remain.
        assert!(fragment.is_none());
        assert!(n_seam_bytes > 0, "test must actually exercise seams");
        assert_eq!(n_record_lines, 100);
    }

    /// Compresses `text` into whole blocks of `block_data` bytes, each cut
    /// at block boundaries, plus a trailing EOF marker. Returns the block
    /// list *without* the marker so tests can splice it anywhere.
    fn compress_blocks(text: &[u8], block_data: usize) -> Vec<Vec<u8>> {
        let mut deflater = Deflater::new(6).expect("invariant: level 6 is valid");
        let mut blocks = Vec::new();
        for piece in text.chunks(block_data) {
            let mut block = Vec::new();
            deflater.deflate_block(piece, &mut block).expect("invariant: test block fits");
            blocks.push(block);
        }
        blocks
    }

    fn record_text(n_records: usize) -> Vec<u8> {
        let mut text = header_text(2);
        for i in 0..n_records {
            text.extend_from_slice(
                format!("chr1\t{}\t.\tA\tG\t1\tPASS\t.\tGT\t0/0\t0/1\n", 100 + i * 7).as_bytes(),
            );
        }
        text
    }

    /// Counts complete `chr1` record lines a full scan of `stream` yields.
    fn count_records_in_stream(stream: &[u8], chunk_target: u64) -> usize {
        let mut reader = BgzfChunkReader::with_chunk_target(Cursor(stream.to_vec()), chunk_target);
        let mut inflater = crate::bgzf::inflate::Inflater::new();
        let _ = reader.scan_header(&mut inflater).expect("invariant: header present");
        let mut n = 0_usize;
        let mut scratch = Vec::new();
        let mut fragment: Vec<u8> = Vec::new();
        while let Some(chunk) = reader.next_chunk().expect("invariant: chunking") {
            let view = inflate_chunk(&chunk, &mut inflater, &mut scratch).expect("invariant: inflate");
            let text = view.text;
            let mut all = std::mem::take(&mut fragment);
            all.extend_from_slice(text);
            let mut lines = all.split(|&b| b == b'\n').collect::<Vec<_>>();
            // The final piece may be an incomplete record; keep it pending
            // unless it ends with a newline.
            fragment = if text.last() == Some(&b'\n') {
                Vec::new()
            } else {
                lines.pop().unwrap_or_default().to_vec()
            };
            n += lines.iter().filter(|l| l.starts_with(b"chr1\t")).count();
        }
        assert!(fragment.is_empty(), "stream ends mid-record");
        n
    }

    /// Regression (01-review D2): the EOF marker and empty members may sit
    /// anywhere in a concatenated stream — mid-chunk or at a chunk start —
    /// and must never terminate the read. Every splice position of the EOF
    /// marker across the whole stream must still yield every record.
    #[test]
    fn empty_members_anywhere_never_end_the_stream() {
        let text = record_text(60);
        let blocks = compress_blocks(&text, 24);
        assert!(blocks.len() >= 4, "test needs several blocks");
        let n_records = 60;

        for splice in 0..=blocks.len() {
            let mut stream: Vec<u8> = Vec::new();
            for (i, block) in blocks.iter().enumerate() {
                if i == splice {
                    stream.extend_from_slice(&crate::bgzf::block::EOF_BLOCK);
                    // An extra non-canonical empty member at the same spot.
                    let mut empty = Vec::new();
                    Deflater::new(6)
                        .expect("invariant: level 6 is valid")
                        .deflate_block(b"", &mut empty)
                        .expect("invariant: empty member");
                    stream.extend_from_slice(&empty);
                }
                stream.extend_from_slice(block);
            }
            if splice == blocks.len() {
                stream.extend_from_slice(&crate::bgzf::block::EOF_BLOCK);
            }
            // Two chunk targets: block-aligned and not, so the empty block
            // also lands at a chunk *start* (the case that used to end the
            // whole stream).
            assert_eq!(
                count_records_in_stream(&stream, 40),
                n_records,
                "EOF marker spliced before block {splice} (chunk target 40)"
            );
            assert_eq!(
                count_records_in_stream(&stream, 72),
                n_records,
                "EOF marker spliced before block {splice} (chunk target 72)"
            );
        }
    }

    /// Regression (01-review D2): `cat a.vcf.gz b.vcf.gz` — two complete
    /// BGZF files joined — must be read end to end.
    #[test]
    fn concatenated_files_read_through_every_join() {
        let a = record_text(30);
        let b = record_text(30);
        let mut stream = compress_lines(&[&a], 48);
        stream.extend_from_slice(&compress_lines(&[&b], 48));
        assert_eq!(count_records_in_stream(&stream, 64), 60);
    }
}
