//! BGZF chunk reader: the serial reader of BAPP (`docs/spec/architecture.md`
//! §6.2). It only parses 18-byte block headers — no inflation — and cuts
//! chunks at block boundaries once the compressed target is reached.

use std::io::Read;
use std::sync::Arc;

use reed_core::ChunkId;

use crate::bgzf::block::{BGZF_HEADER_LEN, is_empty_member, is_eof_block, parse_block_size};
use crate::chunk::{CHUNK_TARGET_COMPRESSED, ChunkBytes, ChunkPayload};
use crate::error::IoError;

/// Locates the end of the VCF header and hands out chunks from there.
pub struct BgzfChunkReader<R: Read> {
    inner: R,
    chunk_target: u64,
    next_chunk_id: u64,
    /// Compressed bytes read but not yet handed to a chunk: the block where
    /// the header ends (chunk 0 starts here, so it is not read twice), and
    /// the tail of the last slab read that belongs to the next chunk.
    carry: Vec<u8>,
    /// Compressed offset of `carry[0]`.
    carry_offset: u64,
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
            carry: Vec::new(),
            carry_offset: 0,
            compressed_offset: 0,
            first_in_block_offset: 0,
            header_done: false,
            last_was_eof_marker: false,
            warned_missing_eof: false,
        }
    }

    /// Reads up to `want` bytes into `out`, directly into the chunk buffer
    /// (one syscall per call in practice). Returns how many bytes arrived;
    /// `0` means end of stream.
    ///
    /// The old path read through a `BufReader`, which copied every byte
    /// twice (kernel → buffer, buffer → chunk); `--timing` showed the reader
    /// at 87% of the stats wall, so the copy is gone (doc_v4 03 §3.4 step 2).
    fn read_into(&mut self, out: &mut Vec<u8>, want: usize) -> Result<usize, IoError> {
        let start = out.len();
        out.resize(start + want, 0);
        let mut filled = 0;
        while filled < want {
            match self.inner.read(&mut out[start + filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(err) if err.kind() == std::io::ErrorKind::Interrupted => {}
                Err(err) => {
                    out.truncate(start);
                    return Err(err.into());
                }
            }
        }
        out.truncate(start + filled);
        self.compressed_offset += filled as u64;
        Ok(filled)
    }

    /// Grows `out` until it holds at least `chunk_target` bytes of whole
    /// BGZF blocks, and returns the offset just past the last whole block
    /// (bytes beyond it belong to the next chunk).
    ///
    /// Empty members (ISIZE = 0: the EOF marker and concatenation joins) are
    /// removed in place and never terminate the stream — the stream ends
    /// only at the underlying EOF (01-review D2).
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedBgzf`] for trailing garbage, a truncated header
    /// or a truncated payload.
    fn fill_to_target(
        &mut self,
        out: &mut Vec<u8>,
        blocks: &mut Vec<(u32, u32)>,
    ) -> Result<usize, IoError> {
        let mut pos = 0_usize;
        loop {
            // Parse every whole block currently in `out`.
            while out.len() >= pos + BGZF_HEADER_LEN {
                let Some(total) = parse_block_size(&out[pos..pos + BGZF_HEADER_LEN]) else {
                    return Err(IoError::MalformedBgzf {
                        locus: None,
                        reason: "not a BGZF block header (trailing garbage after the last block?)",
                    });
                };
                if out.len() < pos + total {
                    break; // The block spans further: read more bytes.
                }
                if is_empty_member(&out[pos..pos + total]) {
                    self.last_was_eof_marker = is_eof_block(&out[pos..pos + total]);
                    out.drain(pos..pos + total);
                    continue; // The next block now starts at `pos`.
                }
                self.last_was_eof_marker = false;
                blocks.push((
                    u32::try_from(total).unwrap_or(u32::MAX),
                    crate::bgzf::block::block_inflated_size(&out[pos..pos + total]),
                ));
                pos += total;
            }
            // At least one whole block and the target reached: the chunk is
            // complete (the cut is always a block boundary).
            if pos > 0 && out.len() as u64 >= self.chunk_target {
                return Ok(pos);
            }
            // Read what is still missing. When a block announced more bytes
            // than the target, read exactly that block's remainder so a
            // chunk is never empty and the loop always progresses.
            let mut want = self.chunk_target.saturating_sub(out.len() as u64);
            if want == 0 {
                want = if out.len() >= BGZF_HEADER_LEN {
                    parse_block_size(&out[..BGZF_HEADER_LEN])
                        .map_or(1, |total| (total as u64).saturating_sub(out.len() as u64))
                } else {
                    (BGZF_HEADER_LEN - out.len()) as u64
                };
            }
            let want = usize::try_from(want.max(1)).unwrap_or(usize::MAX);
            if self.read_into(out, want)? == 0 {
                self.warn_if_missing_eof_marker();
                if pos == out.len() {
                    return Ok(pos); // Clean end at a block boundary.
                }
                return Err(IoError::MalformedBgzf {
                    locus: None,
                    reason: if out.len() - pos < BGZF_HEADER_LEN {
                        "stream ends inside a block header"
                    } else {
                        "stream ends inside a block payload"
                    },
                });
            }
        }
    }

    /// Reads one BGZF block (header + payload) and its compressed offset;
    /// `None` at the end of the underlying stream.
    ///
    /// Reads *exactly* one block: the header scan uses this before any chunk
    /// exists, so buffered read-ahead would drop bytes that the chunk walk
    /// still needs.
    ///
    /// Empty members (ISIZE = 0, including the canonical EOF marker) are
    /// skipped — concatenated BGZF files carry one at every join, and the
    /// stream ends only at the underlying EOF (01-review D2). A clean end
    /// without a trailing EOF marker warns once (htslib behavior).
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedBgzf`] for a truncated header, trailing garbage or
    /// a short payload; transport errors pass through.
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
            tracing::warn!("BGZF stream ends without an EOF marker");
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
                self.carry = block;
                self.carry_offset = block_offset;
                // The offset is *within this block's* inflated text, not in
                // the accumulated header text. A BCF stream NUL-terminates
                // the header, so its records start one byte later.
                let is_bcf = text.starts_with(b"BCF\x02");
                let record_start = if is_bcf && text.get(end) == Some(&0) {
                    end + 1
                } else {
                    end
                };
                self.first_in_block_offset = (record_start - text_len_before) as u64;
                self.header_done = true;
                text.truncate(end);
                // A BCF stream prefixes the VCF header with the magic and
                // `l_text`; the header facts start at the VCF text
                // (doc_v3 02-bcf.md §2.1).
                let start = if text.starts_with(b"BCF\x02") { 9 } else { 0 };
                let header_text = text[start..].to_vec();
                let n_samples = count_samples(&header_text).ok_or_else(|| {
                    IoError::Header("#CHROM line found but has no sample columns".into())
                })?;
                return Ok(HeaderScan {
                    text: header_text,
                    n_samples,
                });
            }
        }
    }

    /// Yields the next chunk of compressed blocks, or `None` at end of
    /// stream. Requires [`BgzfChunkReader::scan_header`] first.
    pub fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, IoError> {
        if !self.header_done {
            return Err(IoError::Header("scan_header must run before chunks".into()));
        }
        // Reserve the target up front: growth reallocs would copy the whole
        // chunk several times over.
        let mut compressed = Vec::with_capacity(
            usize::try_from(self.chunk_target).unwrap_or(usize::MAX - (1 << 17)) + (1 << 17),
        );
        let mut blocks = Vec::with_capacity(24);
        // The chunk starts with the bytes carried over from the previous
        // read (the header block for chunk 0, else a partial block tail).
        let compressed_start = self.carry_offset;
        compressed.extend_from_slice(&self.carry);
        self.carry.clear();
        let cut = self.fill_to_target(&mut compressed, &mut blocks)?;
        if compressed.is_empty() {
            return Ok(None);
        }
        // Bytes beyond the last whole block belong to the next chunk; the
        // overflow is bounded by one block (<= 64 KiB + 26).
        self.carry_offset = compressed_start + cut as u64;
        self.carry.extend_from_slice(&compressed[cut..]);
        compressed.truncate(cut);
        let id = ChunkId(self.next_chunk_id);
        self.next_chunk_id += 1;
        Ok(Some(ChunkBytes {
            id,
            compressed: ChunkPayload::from(compressed),
            first_in_block_offset: std::mem::take(&mut self.first_in_block_offset),
            plain: false,
            compressed_start,
            blocks,
        }))
    }
}

/// A memory-mapped BGZF file's chunk source (doc_v4 03 §3.4 step 3): a chunk
/// is a byte range of the mapping, so the reader copies nothing and issues no
/// read syscall. Chunk boundaries follow the same rule as the streaming
/// reader — whole blocks until the compressed target is reached — so the
/// two readers produce identical chunks and identical output.
pub struct MmapChunkSource {
    map: Arc<memmap2::Mmap>,
    chunk_target: u64,
    next_chunk_id: u64,
    /// Compressed offset of the next chunk's first block.
    next_offset: usize,
    /// Records' start offset inside that block's inflated text.
    first_in_block_offset: u64,
    /// The canonical EOF marker was the last empty member seen.
    last_was_eof_marker: bool,
    warned_missing_eof: bool,
}

impl MmapChunkSource {
    /// Maps `file` and scans its header. The mapping is read-only and lives
    /// exactly as long as the source.
    ///
    /// # Errors
    ///
    /// Mapping failures (the caller falls back to the streaming reader) and
    /// [`IoError::Header`] when the header is missing or unparsable.
    pub fn open(file: std::fs::File) -> Result<(Self, crate::header::HeaderFacts), IoError> {
        // SAFETY: `Mmap::map` is unsafe because a concurrent writer could
        // change the file under the mapping (the worst case is a SIGBUS on
        // truncation). reed opens the input read-only, never writes it, and
        // the mapping outlives every chunk; the callers fall back to the
        // streaming reader when the input is not a regular file. This is the
        // single audited `unsafe` site in the workspace (doc_v6 06 Y0 deleted
        // the AVX-512 count kernel).
        #[allow(unsafe_code)]
        let map = unsafe { memmap2::Mmap::map(&file) }?;
        let map = Arc::new(map);
        // The header scan reuses the streaming reader over the mapped bytes:
        // it copies a handful of header blocks, which is negligible.
        let mut probe = BgzfChunkReader::new(std::io::Cursor::new(&map[..]));
        let scan = probe.scan_header(&mut crate::bgzf::inflate::Inflater::new())?;
        let facts = crate::header::HeaderFacts::parse(&scan.text)?;
        // The probe borrows the mapping; read what it found before moving
        // the `Arc` into the source.
        let first_offset = usize::try_from(probe.carry_offset).unwrap_or(0);
        let first_in_block_offset = probe.first_in_block_offset;
        drop(probe);
        let source = Self {
            map,
            chunk_target: CHUNK_TARGET_COMPRESSED,
            next_chunk_id: 0,
            next_offset: first_offset,
            first_in_block_offset,
            last_was_eof_marker: false,
            warned_missing_eof: false,
        };
        Ok((source, facts))
    }

    /// Warns once when the stream ended without the canonical EOF marker.
    fn warn_if_missing_eof_marker(&mut self) {
        if !self.warned_missing_eof && !self.last_was_eof_marker {
            self.warned_missing_eof = true;
            tracing::warn!("BGZF stream ends without an EOF marker");
        }
    }
}

impl crate::chunk::ChunkSource for MmapChunkSource {
    fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, IoError> {
        let start = self.next_offset;
        let mut pos = start;
        let mut blocks: Vec<(u32, u32)> = Vec::with_capacity(24);
        let mut total = 0_u64;
        loop {
            if pos >= self.map.len() {
                break;
            }
            if pos + BGZF_HEADER_LEN > self.map.len() {
                return Err(IoError::MalformedBgzf {
                    locus: None,
                    reason: "stream ends inside a block header",
                });
            }
            let Some(size) = parse_block_size(&self.map[pos..pos + BGZF_HEADER_LEN]) else {
                return Err(IoError::MalformedBgzf {
                    locus: None,
                    reason: "not a BGZF block header (trailing garbage after the last block?)",
                });
            };
            if pos + size > self.map.len() {
                return Err(IoError::MalformedBgzf {
                    locus: None,
                    reason: "stream ends inside a block payload",
                });
            }
            let block = &self.map[pos..pos + size];
            if is_empty_member(block) {
                // Empty members stay in the chunk's byte range (their
                // inflated length is zero) so the block offsets remain
                // contiguous for the index.
                self.last_was_eof_marker = is_eof_block(block);
                blocks.push((u32::try_from(size).unwrap_or(u32::MAX), 0));
                pos += size;
                continue;
            }
            self.last_was_eof_marker = false;
            blocks.push((
                u32::try_from(size).unwrap_or(u32::MAX),
                crate::bgzf::block::block_inflated_size(block),
            ));
            pos += size;
            total += size as u64;
            if total >= self.chunk_target {
                break;
            }
        }
        if blocks.is_empty() {
            self.warn_if_missing_eof_marker();
            return Ok(None);
        }
        let id = ChunkId(self.next_chunk_id);
        self.next_chunk_id += 1;
        let first_in_block_offset = if self.next_chunk_id == 1 {
            self.first_in_block_offset
        } else {
            0
        };
        self.next_offset = pos;
        Ok(Some(ChunkBytes {
            id,
            compressed: ChunkPayload::Mapped {
                map: Arc::clone(&self.map),
                range: start..pos,
            },
            first_in_block_offset,
            plain: false,
            compressed_start: start as u64,
            blocks,
        }))
    }
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
    Some(memchr::memchr_iter(b'\t', line).count().saturating_sub(8))
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
            deflater
                .deflate_block(block, &mut out)
                .expect("invariant: test block fits");
        }
        out.extend_from_slice(&crate::bgzf::block::EOF_BLOCK);
        out
    }

    fn header_text(n_samples: usize) -> Vec<u8> {
        let mut h = String::from("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100000>\n");
        h.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for i in 0..n_samples {
            use std::fmt::Write as _;
            write!(h, "\ts{i}").expect("invariant: writing to String");
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
        let scan = reader
            .scan_header(&mut inflater)
            .expect("invariant: header present");
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
        let _ = reader
            .scan_header(&mut inflater)
            .expect("invariant: header present");

        let mut n_record_lines = 0_usize;
        let mut n_seam_bytes = 0_usize;
        // Chained fragments of the record spanning the seam(s): the previous
        // chunk's tail, plus any number of whole-chunk head fragments, until
        // a newline completes the record (spec §6.2 step 4, chained case).
        let mut fragment: Option<Vec<u8>> = None;
        let mut scratch = Vec::new();
        while let Some(chunk) = reader.next_chunk().expect("invariant: chunking") {
            let view =
                inflate_chunk(&chunk, &mut inflater, &mut scratch).expect("invariant: inflate");
            if let Some(head) = view.head_bytes(view.text) {
                let mut seam = fragment.take().unwrap_or_default();
                n_seam_bytes += head.len();
                seam.extend_from_slice(head);
                if seam.ends_with(b"\n") {
                    assert_eq!(
                        memchr::memchr_iter(b'\n', &seam).count(),
                        1,
                        "seam is one line"
                    );
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
            deflater
                .deflate_block(piece, &mut block)
                .expect("invariant: test block fits");
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
        let _ = reader
            .scan_header(&mut inflater)
            .expect("invariant: header present");
        let mut n = 0_usize;
        let mut scratch = Vec::new();
        let mut fragment: Vec<u8> = Vec::new();
        while let Some(chunk) = reader.next_chunk().expect("invariant: chunking") {
            let view =
                inflate_chunk(&chunk, &mut inflater, &mut scratch).expect("invariant: inflate");
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
