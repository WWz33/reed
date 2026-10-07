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

use crate::bgzf::block::{parse_block_size, BGZF_HEADER_LEN};
use crate::chunk::{ChunkBytes, CHUNK_TARGET_COMPRESSED};
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
    chunk_target: u64,
    next_chunk_id: u64,
    pending_block: Option<Vec<u8>>,
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
            chunk_target: CHUNK_TARGET_COMPRESSED,
            next_chunk_id: 0,
            pending_block: None,
            first_in_block_offset: start & 0xFFFF,
            done: false,
        })
    }

    /// Reads one block; `None` once the region's range is exhausted or the
    /// underlying stream ends.
    fn read_block(&mut self) -> Result<Option<Vec<u8>>, IoError> {
        if self.done {
            return Ok(None);
        }
        let position = self.inner.stream_position()?;
        if position > self.end_coffset {
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
        self.inner.read_exact(&mut block[BGZF_HEADER_LEN..])?;
        Ok(Some(block))
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
        let first = match self.pending_block.take() {
            Some(block) => block,
            None => match self.read_block()? {
                Some(block) => block,
                None => return Ok(None),
            },
        };
        let mut compressed = first;
        while (compressed.len() as u64) < self.chunk_target {
            match self.read_block()? {
                Some(block) => compressed.extend_from_slice(&block),
                None => break,
            }
        }
        let id = ChunkId(self.next_chunk_id);
        self.next_chunk_id += 1;
        Ok(Some(ChunkBytes {
            id,
            compressed: std::sync::Arc::from(compressed.into_boxed_slice()),
            first_in_block_offset: std::mem::take(&mut self.first_in_block_offset),
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
            chunk_target: CHUNK_TARGET_COMPRESSED,
            next_chunk_id: 0,
            pending_block: None,
            first_in_block_offset: 0,
            done: true,
        }
    }
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
        let mut start = u64::MAX;
        let mut end = 0_u64;
        for region in regions {
            let tid = contig_names
                .iter()
                .position(|name| name == &region.contig)
                .ok_or_else(|| {
                    IoError::Header(format!(
                        "region contig {:?} is not in the header",
                        region.contig
                    ))
                })?;
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
            reader.read_exact(&mut magic).map_err(|_| {
                IoError::UnsupportedInput("index file is shorter than its magic")
            })?;
        }
        match &magic {
            b"CSI\x01" => Ok(Self::Csi(Box::new(noodles::csi::fs::read(path)?))),
            b"TBI\x01" => Ok(Self::Tbi(Box::new(noodles::tabix::fs::read(path)?))),
            _ => Err(IoError::UnsupportedInput("index is neither CSI nor TBI")),
        }
    }
}

impl IndexQuery for FileIndex {
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
        Ok(chunks.iter().map(|c| (c.start().into(), c.end().into())).collect())
    }
}
