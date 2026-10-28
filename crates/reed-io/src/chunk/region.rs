//! Region-restricted BGZF chunk source (02-cli.md §2, `-r/-R`): the index
//! maps the requested regions to a compressed byte range, the reader seeks
//! there and yields chunks only for that range. Records outside the regions
//! are dropped by the caller's region predicate (the same rule `-t` uses),
//! so the index only saves I/O.
//!
//! The first chunk of the range starts mid-record by construction: the
//! source inflates its first block once to find the first line boundary and
//! hands that offset to the pipeline as `first_in_block_offset` (the same
//! mechanism the header-carrying first chunk uses).

use std::io::{Read, Seek, SeekFrom};

use reed_core::ChunkId;

use crate::bgzf::block::{BGZF_HEADER_LEN, parse_block_size};
use crate::chunk::{CHUNK_TARGET_COMPRESSED, ChunkBytes};
use crate::error::IoError;

/// One parsed region: contig, 1-based inclusive bounds (`None` = open).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    /// Contig name.
    pub contig: String,
    /// 1-based start (`None` = from the beginning).
    pub start1: Option<u64>,
    /// 1-based end (`None` = to the end).
    pub end1: Option<u64>,
}

/// Reads only the compressed byte range covering `regions`.
pub struct RegionChunkSource<R: Read + Seek> {
    inner: R,
    end_coffset: u64,
    /// In-block offset of the range end (low 16 bits of the index's end
    /// virtual offset): a record boundary. The block read at `end_coffset`
    /// may extend past it, so the chunk carrying that block truncates its
    /// inflated text here (doc_v21 M3).
    end_uoffset: u64,
    chunk_target: u64,
    next_chunk_id: u64,
    pending_block: Option<Vec<u8>>,
    pending_block_offset: u64,
    first_in_block_offset: u64,
    done: bool,
}

impl<R: Read + Seek> RegionChunkSource<R> {
    /// Computes the byte range for `regions` from `index` and positions the
    /// reader at its first block.
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when a region names an unknown contig, and
    /// [`IoError::UnsupportedInput`] when the index holds no chunks for it.
    pub fn open(
        mut inner: R,
        index: &dyn IndexQuery,
        regions: &[Region],
        contig_names: &[String],
    ) -> Result<Self, IoError> {
        let Some((start, end)) = index.query_range(regions, contig_names)? else {
            return Ok(Self::open_empty(inner));
        };
        // The index's virtual offset carries both halves: the compressed
        // offset of the record's block and its offset inside that block's
        // inflated text — exactly the chunk's `first_in_block_offset`.
        inner.seek(SeekFrom::Start(start >> 16))?;
        Ok(Self {
            inner,
            end_coffset: end >> 16,
            end_uoffset: end & 0xFFFF,
            chunk_target: CHUNK_TARGET_COMPRESSED,
            next_chunk_id: 0,
            pending_block: None,
            pending_block_offset: 0,
            first_in_block_offset: start & 0xFFFF,
            done: false,
        })
    }

    /// Reads one block; `None` once the region's range is exhausted or the
    /// underlying stream ends.
    fn read_block(&mut self) -> Result<Option<(Vec<u8>, u64)>, IoError> {
        if self.done {
            return Ok(None);
        }
        let position = self.inner.stream_position()?;
        // doc_v21 M3: a range whose end offset sits exactly on a block
        // boundary needs nothing from the block starting there.
        if position > self.end_coffset || (position == self.end_coffset && self.end_uoffset == 0) {
            self.done = true;
            return Ok(None);
        }
        let mut header = [0_u8; BGZF_HEADER_LEN];
        let mut filled = 0;
        while filled < BGZF_HEADER_LEN {
            match self.inner.read(&mut header[filled..]) {
                Ok(0) => {
                    self.done = true;
                    return Ok(None);
                }
                Ok(n) => filled += n,
                Err(err) if err.kind() == std::io::ErrorKind::Interrupted => {}
                Err(err) => return Err(err.into()),
            }
        }
        let Some(total) = parse_block_size(&header) else {
            return Err(IoError::MalformedBgzf {
                locus: None,
                reason: "not a BGZF block header in the region range",
            });
        };
        let mut block = Vec::with_capacity(total);
        block.extend_from_slice(&header);
        block.resize(total, 0);
        let block_offset = position;
        self.inner.read_exact(&mut block[BGZF_HEADER_LEN..])?;
        Ok(Some((block, block_offset)))
    }

    /// The next chunk of the range, or `None` at its end.
    ///
    /// # Errors
    ///
    /// Transport failures.
    pub fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, IoError> {
        if self.done && self.pending_block.is_none() {
            return Ok(None);
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
        let mut last_block_offset = compressed_start;
        while (compressed.len() as u64) < self.chunk_target {
            match self.read_block()? {
                Some((block, offset)) => {
                    last_block_offset = offset;
                    blocks.push(block_pair(&block));
                    compressed.extend_from_slice(&block);
                }
                None => break,
            }
        }
        // doc_v21 M3: when this chunk ends on the block the index range ends
        // in, that block's inflated text runs past the range (BGZF blocks are
        // cut at 64 KB, not at record boundaries). Keep only the prefix up to
        // the end virtual offset's in-block position — a record boundary — so
        // the engine never sees an unterminated trailing record.
                let truncate_to = (last_block_offset == self.end_coffset).then(|| {
            let prior: u64 = blocks[..blocks.len() - 1]
                .iter()
                .map(|&(_, ilen)| u64::from(ilen))
                .sum();
            prior + self.end_uoffset
        });
        let id = ChunkId(self.next_chunk_id);
        self.next_chunk_id += 1;
        Ok(Some(ChunkBytes {
            id,
            compressed: crate::chunk::ChunkPayload::from(compressed),
            first_in_block_offset: std::mem::take(&mut self.first_in_block_offset),
            plain: false,
            compressed_start,
            blocks,
            truncate_to,
        }))
    }
}

impl<R: Read + Seek + Send> crate::chunk::ChunkSource for RegionChunkSource<R> {
    fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, IoError> {
        RegionChunkSource::next_chunk(self)
    }
}

impl<R: Read + Seek> RegionChunkSource<R> {
    /// A source that yields nothing (a region with no indexed chunks).
    #[must_use]
    pub fn open_empty(inner: R) -> Self {
        Self {
            inner,
            end_coffset: 0,
            end_uoffset: 0,
            chunk_target: CHUNK_TARGET_COMPRESSED,
            next_chunk_id: 0,
            pending_block: None,
            pending_block_offset: 0,
            first_in_block_offset: 0,
            done: true,
        }
    }
}

/// One block's `(compressed length, inflated length)`.
fn block_pair(block: &[u8]) -> (u32, u32) {
    (
        u32::try_from(block.len()).unwrap_or(u32::MAX),
        crate::bgzf::block::block_inflated_size(block),
    )
}

/// The slice of the index API the region source needs (keeps reed-io free of
/// a hard noodles-index dependency in its public surface).
pub trait IndexQuery {
    /// Chunks covering `tid`'s `[start1, end1]` (1-based, `None` = open),
    /// as `(compressed_offset << 16 | uncompressed_offset)` pairs.
    ///
    /// # Errors
    ///
    /// Index read failures.
    fn query(
        &self,
        tid: usize,
        start1: Option<u64>,
        end1: Option<u64>,
    ) -> Result<Vec<(u64, u64)>, IoError>;

    /// Whether the index carries a reference-sequence name table. A TBI
    /// always does (first-record-appearance order); a CSI does when htslib
    /// wrote one. When present, [`Self::tid_of`] is authoritative — an id
    /// resolved from the VCF header order would be wrong (doc_v21 M2).
    fn has_reference_names(&self) -> bool {
        false
    }

    /// The index's own id for `name`, or `None` when the index has no name
    /// table or does not name it.
    fn tid_of(&self, name: &str) -> Option<usize> {
        let _ = name;
        None
    }

    /// The virtual-offset range covering every region, or `None` when no
    /// chunk matches.
    ///
    /// # Errors
    ///
    /// Unknown contigs and index failures.
    fn query_range(
        &self,
        regions: &[Region],
        contig_names: &[String],
    ) -> Result<Option<(u64, u64)>, IoError> {
        let unknown = |name: &str| {
            IoError::Header(format!("region contig {name:?} is not in the header"))
        };
        let mut start = u64::MAX;
        let mut end = 0_u64;
        for region in regions {
            // doc_v21 M2: the id must be the index's own. A TBI's reference
            // table (and an htslib CSI's) is first-record-appearance order,
            // not the VCF header order, so a header position would pick the
            // wrong bins — or, when the index has fewer refs than the header,
            // fail with "invalid reference sequence ID". Fall back to the
            // header position only for a CSI that names nothing (contig-id
            // keyed: there the header position *is* the id).
            let tid = match self.tid_of(&region.contig) {
                Some(tid) => tid,
                None if self.has_reference_names() => {
                    // Named contigs and no match: the contig has no records.
                    // Still reject one the header never declared.
                    if contig_names.iter().any(|name| name == &region.contig) {
                        continue;
                    }
                    return Err(unknown(&region.contig));
                }
                None => contig_names
                    .iter()
                    .position(|name| name == &region.contig)
                    .ok_or_else(|| unknown(&region.contig))?,
            };
            for (s, e) in self.query(tid, region.start1, region.end1)? {
                start = start.min(s);
                end = end.max(e);
            }
        }
        Ok((start != u64::MAX).then_some((start, end)))
    }
}

/// A noodles-backed index file (CSI or TBI), sniffed by magic.
pub enum FileIndex {
    /// Coordinate-sorted index (`.csi`).
    Csi(Box<noodles::csi::Index>),
    /// Tabix index (`.tbi`).
    Tbi(Box<noodles::tabix::Index>),
}

impl FileIndex {
    /// Opens an index file, detecting CSI vs TBI from its magic.
    ///
    /// # Errors
    ///
    /// [`IoError::Io`] when the file cannot be read, and
    /// [`IoError::UnsupportedInput`] for an unknown magic.
    pub fn open(path: &std::path::Path) -> Result<Self, IoError> {
        // htslib writes indexes as BGZF streams, so the magic is read
        // through a BGZF reader.
        let mut magic = [0_u8; 4];
        {
            let file = std::fs::File::open(path)?;
            let mut reader =
                noodles::bgzf::io::Reader::new(std::io::BufReader::with_capacity(1 << 16, file));
            reader
                .read_exact(&mut magic)
                .map_err(|_| IoError::UnsupportedInput("index file is shorter than its magic"))?;
        }
        match &magic {
            b"CSI\x01" => Ok(Self::Csi(Box::new(noodles::csi::fs::read(path)?))),
            b"TBI\x01" => Ok(Self::Tbi(Box::new(noodles::tabix::fs::read(path)?))),
            _ => Err(IoError::UnsupportedInput("index is neither CSI nor TBI")),
        }
    }
}

impl FileIndex {
    /// The index's reference-sequence name table, when it has one.
    fn reference_names(
        &self,
    ) -> Option<&noodles::csi::binning_index::index::header::ReferenceSequenceNames> {
        use noodles::csi::BinningIndex as _;
        match self {
            Self::Csi(index) => index.header().map(|h| h.reference_sequence_names()),
            Self::Tbi(index) => index.header().map(|h| h.reference_sequence_names()),
        }
    }
}

impl IndexQuery for FileIndex {
    fn has_reference_names(&self) -> bool {
        self.reference_names().is_some_and(|names| !names.is_empty())
    }

    fn tid_of(&self, name: &str) -> Option<usize> {
        self.reference_names()?
            .iter()
            .position(|entry| <bstr::BString as AsRef<[u8]>>::as_ref(entry) == name.as_bytes())
    }

    fn query(
        &self,
        tid: usize,
        start1: Option<u64>,
        end1: Option<u64>,
    ) -> Result<Vec<(u64, u64)>, IoError> {
        use noodles::core::{Position, region::Interval};
        let position = |v: u64| {
            usize::try_from(v)
                .ok()
                .and_then(|v| Position::try_from(v).ok())
                .ok_or_else(|| IoError::Header(format!("region coordinate {v} is out of range")))
        };
        let interval = match (start1, end1) {
            (Some(s), Some(e)) => Interval::from(position(s)?..=position(e)?),
            (Some(s), None) => Interval::from(position(s)?..),
            (None, Some(e)) => Interval::from(..=position(e)?),
            (None, None) => Interval::from(..),
        };
        use noodles::csi::BinningIndex as _;
        let chunks = match self {
            Self::Csi(index) => index.query(tid, interval),
            Self::Tbi(index) => index.query(tid, interval),
        }
        .map_err(|e| IoError::Io(e.to_string()))?;
        Ok(chunks
            .iter()
            .map(|c| (c.start().into(), c.end().into()))
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A one-region index whose single chunk spans `[start_voff, end_voff)`.
    /// The region source only consumes the virtual offsets, so a stub is
    /// enough to drive the end-of-range logic exactly.
    struct StubIndex {
        start: u64,
        end: u64,
        names: Vec<&'static str>,
    }

    impl IndexQuery for StubIndex {
        fn query(
            &self,
            tid: usize,
            _start1: Option<u64>,
            _end1: Option<u64>,
        ) -> Result<Vec<(u64, u64)>, IoError> {
            Ok(if tid == 0 {
                vec![(self.start, self.end)]
            } else {
                Vec::new()
            })
        }

        fn has_reference_names(&self) -> bool {
            !self.names.is_empty()
        }

        fn tid_of(&self, name: &str) -> Option<usize> {
            self.names.iter().position(|n| *n == name)
        }
    }

    fn bgzf(blocks: &[&[u8]]) -> Vec<u8> {
        let mut deflater =
            crate::bgzf::deflate::Deflater::new(6).expect("invariant: level 6 is valid");
        let mut out = Vec::new();
        for block in blocks {
            deflater
                .deflate_block(block, &mut out)
                .expect("invariant: block fits");
        }
        out.extend_from_slice(&crate::bgzf::block::EOF_BLOCK);
        out
    }

    fn voff(coffset: u64, uoffset: u64) -> u64 {
        (coffset << 16) | uoffset
    }

    /// doc_v21 M3: the index's end virtual offset is a record boundary inside
    /// a block; everything after it in that block belongs to later records and
    /// must be dropped. The region source used to hand on the whole block, and
    /// the engine then reported an unterminated trailing record.
    #[test]
    fn m3_chunk_truncates_to_the_range_end_inside_the_last_block() {
        // Block 1: three records; the range ends after the first one.
        let b1: &[u8] = b"chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\nchr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\nchr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
        let b2: &[u8] = b"chr1\t4\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
        let bytes = bgzf(&[b1, b2]);
        let block1_len = {
            let mut deflater =
                crate::bgzf::deflate::Deflater::new(6).expect("invariant: level 6");
            let mut one = Vec::new();
            deflater.deflate_block(b1, &mut one).expect("invariant: fits");
            one.len() as u64
        };
        let end_uoffset = b1
            .iter()
            .position(|&c| c == b'\n')
            .expect("invariant: first record ends")
            as u64
            + 1;
        // The chunk starts at block 0 and ends inside block 0, after record 1.
        let index = StubIndex {
            start: voff(0, 0),
            end: voff(0, end_uoffset),
            names: vec!["chr1"],
        };
        let mut source = RegionChunkSource::open(
            std::io::Cursor::new(bytes),
            &index,
            &[Region {
                contig: "chr1".into(),
                start1: Some(1),
                end1: Some(1),
            }],
            &["chr1".to_string()],
        )
        .expect("invariant: source opens");

        let chunk = source
            .next_chunk()
            .expect("invariant: read")
            .expect("invariant: one chunk");
        // The end offset lies inside block 0, so the chunk must carry a
        // truncation point (just past record 1), not the whole blocks.
        assert_eq!(
            chunk.truncate_to,
            Some(end_uoffset),
            "the range end must cut the inflated text at a record boundary"
        );
        assert_eq!(block1_len, chunk.blocks[0].0 as u64);
        // Inflating with the truncation applied drops records 2 and 3.
        let mut inflater = crate::bgzf::inflate::Inflater::new();
        let mut scratch = Vec::new();
        let view = crate::chunk::inflate_chunk(&chunk, &mut inflater, &mut scratch)
            .expect("invariant: inflates");
        assert_eq!(view.text, b"chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n");
        assert!(view.tail.is_none(), "truncation lands on a record boundary");
    }

    /// The control: when the range ends exactly on a block boundary there is
    /// nothing to truncate, and the chunk ends at a record boundary anyway.
    #[test]
    fn m3_block_aligned_range_end_needs_no_truncation() {
        let b1: &[u8] = b"chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
        let b2: &[u8] = b"chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
        let bytes = bgzf(&[b1, b2]);
        let block1_len = {
            let mut deflater =
                crate::bgzf::deflate::Deflater::new(6).expect("invariant: level 6");
            let mut one = Vec::new();
            deflater.deflate_block(b1, &mut one).expect("invariant: fits");
            one.len() as u64
        };
        let index = StubIndex {
            start: voff(0, 0),
            end: voff(block1_len, 0),
            names: vec!["chr1"],
        };
        let mut source = RegionChunkSource::open(
            std::io::Cursor::new(bytes),
            &index,
            &[Region {
                contig: "chr1".into(),
                start1: Some(1),
                end1: Some(1),
            }],
            &["chr1".to_string()],
        )
        .expect("invariant: source opens");
        let chunk = source
            .next_chunk()
            .expect("invariant: read")
            .expect("invariant: one chunk");
        assert_eq!(chunk.truncate_to, None, "block-aligned end needs no cut");
    }

    /// doc_v21 M2: the contig id must come from the index's own name table.
    /// A TBI's order is first-record-appearance, not header order, so the
    /// header position picked the wrong bins (or failed with "invalid
    /// reference sequence ID" when the index had fewer refs).
    #[test]
    fn m2_query_range_resolves_tids_through_the_index_names() {
        // Index lists only chr2 at tid 0; the header declares chr1, chr2, chr3.
        let index = StubIndex {
            start: voff(0, 0),
            end: voff(0, 10),
            names: vec!["chr2"],
        };
        let header = ["chr1".to_string(), "chr2".to_string(), "chr3".to_string()];
        let got = index
            .query_range(
                &[Region {
                    contig: "chr2".into(),
                    start1: None,
                    end1: None,
                }],
                &header,
            )
            .expect("chr2 must resolve through the index");
        assert_eq!(got, Some((voff(0, 0), voff(0, 10))));

        // A header contig the index never saw has no records — not an error
        // (routine after `bcftools view -r`).
        let empty = index
            .query_range(
                &[Region {
                    contig: "chr1".into(),
                    start1: None,
                    end1: None,
                }],
                &header,
            )
            .expect("a declared contig with no records is empty, not an error");
        assert_eq!(empty, None);

        // A contig the header never declared is still a usage error.
        let err = index
            .query_range(
                &[Region {
                    contig: "chrZ".into(),
                    start1: None,
                    end1: None,
                }],
                &header,
            )
            .expect_err("unknown contig");
        assert!(err.to_string().contains("not in the header"), "{err}");
    }
}
