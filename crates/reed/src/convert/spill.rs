//! Convert's rendered-row spill blobs (doc_v25 V25-09 ③): the hapmap row
//! stream, the MSA frame stream and their readers. Every on-disk file goes
//! through [`crate::spill::TempSpill`] (O_EXCL + unlink on drop).

use std::fs::{File, OpenOptions};
use std::io::{Read, Write};

use reed_engine::error::EngineError;

/// Rendered hapmap rows: kept in memory for a small run, spilled to an
/// append-only temp file past [`HMP_SPILL_THRESHOLD`] so a whole-genome run
/// holds only one chunk's text at a time (doc_v15 D15-7, bounded form). The
/// temp file is removed on drop — including on failure or abort, since the
/// engine drops the sink.
#[derive(Debug)]
pub(crate) enum RowsBlob {
    Memory(Vec<u8>),
    /// The rows past the memory threshold, in a [`TempSpill`] that removes
    /// the file on drop (doc_v25 V25-09: the same O_EXCL machinery as the
    /// stdin replay, so a pre-planted symlink cannot redirect the spill).
    Spilled {
        spill: crate::spill::TempSpill,
        len: u64,
    },
}

impl Default for RowsBlob {
    fn default() -> Self {
        Self::Memory(Vec::new())
    }
}

/// Spill each buffered blob past this many bytes. A whole-genome run is far
/// larger; a small run never touches the disk.
pub(crate) const SPILL_THRESHOLD: usize = 4 << 20;

/// The sink's total in-memory blob budget: per-sample blobs are small
/// individually (one byte per site), so a whole-genome, many-sample run
/// would only cross the per-blob threshold on the largest ones. Past this
/// total the sink spills whole buffers in a fixed order (bytes unchanged).
pub(crate) const SPILL_TOTAL_BUDGET: usize = 8 << 20;

impl RowsBlob {
    /// Absorbs one chunk's rendered rows. The ordered sink appends in chunk
    /// order, so the bytes stay identical to the in-memory form.
    pub(crate) fn append(
        &mut self,
        rows: &mut Vec<u8>,
        threshold: usize,
    ) -> Result<(), EngineError> {
        if let Self::Spilled { spill, len } = self {
            let mut file = OpenOptions::new()
                .append(true)
                .open(spill.path())
                .map_err(EngineError::Io)?;
            file.write_all(rows).map_err(EngineError::Io)?;
            *len = len.saturating_add(u64::try_from(rows.len()).unwrap_or(u64::MAX));
            return Ok(());
        }
        // A `Drop` type cannot be destructured by move: fill the buffer,
        // then flip the variant once the borrow has ended.
        let mut spilled: Option<(crate::spill::TempSpill, u64)> = None;
        if let Self::Memory(buf) = self {
            if buf.len() + rows.len() <= threshold {
                buf.append(rows);
                return Ok(());
            }
            buf.append(rows);
            let (mut file, spill) =
                crate::spill::TempSpill::create(None, "convert", "spill")
                    .map_err(EngineError::Io)?;
            file.write_all(buf).map_err(EngineError::Io)?;
            file.flush().map_err(EngineError::Io)?;
            let len = u64::try_from(buf.len()).unwrap_or(u64::MAX);
            spilled = Some((spill, len));
        }
        if let Some((spill, len)) = spilled {
            *self = Self::Spilled { spill, len };
        }
        Ok(())
    }

    /// Streams the bytes through `f` in read-sized chunks (the transformed
    /// renderers: `--012` maps each byte back to a cell token).
    pub(crate) fn visit(
        &self,
        mut f: impl FnMut(&[u8]) -> std::io::Result<()>,
    ) -> Result<(), EngineError> {
        match self {
            Self::Memory(bytes) => f(bytes).map_err(EngineError::Io),
            Self::Spilled { spill, .. } => {
                let mut src =
                    std::io::BufReader::new(File::open(spill.path()).map_err(EngineError::Io)?);
                let mut buf = vec![0_u8; 64 << 10];
                loop {
                    let n = src.read(&mut buf).map_err(EngineError::Io)?;
                    if n == 0 {
                        return Ok(());
                    }
                    f(&buf[..n]).map_err(EngineError::Io)?;
                }
            }
        }
    }

    /// In-memory bytes (`0` once spilled).
    #[must_use]
    pub(crate) fn mem_len(&self) -> usize {
        match self {
            Self::Memory(bytes) => bytes.len(),
            Self::Spilled { .. } => 0,
        }
    }

    /// Forces the in-memory buffer to a spill file (no-op when empty or
    /// already spilled).
    pub(crate) fn spill_now(&mut self) -> Result<(), EngineError> {
        let mut spilled: Option<(crate::spill::TempSpill, u64)> = None;
        if let Self::Memory(buf) = self
            && !buf.is_empty()
        {
            let (mut file, spill) =
                crate::spill::TempSpill::create(None, "convert", "spill")
                    .map_err(EngineError::Io)?;
            file.write_all(buf).map_err(EngineError::Io)?;
            file.flush().map_err(EngineError::Io)?;
            let len = u64::try_from(buf.len()).unwrap_or(u64::MAX);
            spilled = Some((spill, len));
        }
        if let Some((spill, len)) = spilled {
            *self = Self::Spilled { spill, len };
        }
        Ok(())
    }

    /// Streams the rows into `out` (the render path).
    pub(crate) fn write_to(&self, out: &mut impl std::io::Write) -> Result<(), EngineError> {
        match self {
            Self::Memory(bytes) => out.write_all(bytes).map_err(EngineError::Io),
            Self::Spilled { spill, len } => {
                let mut src =
                    std::io::BufReader::new(File::open(spill.path()).map_err(EngineError::Io)?);
                let copied = std::io::copy(&mut src, out).map_err(EngineError::Io)?;
                debug_assert_eq!(copied, *len, "invariant: the spill length tracks writes");
                Ok(())
            }
        }
    }
}


/// The MSA sites, one frame stream per contig run (the `runs` merge rule is
/// preserved: a run continues while the contig does not change). Each frame
/// is `pos0: u64 LE || width: u32 LE || vals (n_samples x width)` — the
/// render walks the frames once per sample, so nothing region-sized is ever
/// buffered (doc_v15 D15-7).
pub(crate) struct SpilledMsa {
    pub(crate) runs: Vec<(Vec<u8>, RowsBlob)>,
}

/// A pull-based reader over one contig's frames. The render re-walks the
/// frames per sample, so the reader rewinds: memory slices by index, spill
/// files by `seek`.
pub(crate) struct MsaFrameReader<'a> {
    source: MsaFrameSource<'a>,
    n_samples: usize,
}

enum MsaFrameSource<'a> {
    Memory { bytes: &'a [u8], cursor: usize },
    Spilled { file: std::io::BufReader<File> },
}

impl<'a> MsaFrameReader<'a> {
    pub(crate) fn open(blob: &'a RowsBlob, n_samples: usize) -> Result<Self, EngineError> {
        let source = match blob {
            RowsBlob::Memory(bytes) => MsaFrameSource::Memory { bytes, cursor: 0 },
            RowsBlob::Spilled { spill, .. } => MsaFrameSource::Spilled {
                file: std::io::BufReader::new(
                    File::open(spill.path()).map_err(EngineError::Io)?,
                ),
            },
        };
        Ok(Self { source, n_samples })
    }

    /// Reads the next frame and copies only `sample`'s slice into `buf`.
    /// `Ok(None)` at the end of the stream.
    ///
    /// The render walks the frames once per sample; reading the whole frame
    /// every pass would be the per-sample repeated-IO shape the convert
    /// design forbids (doc_v13 D13-2), so the spill path seeks past the
    /// other samples' bytes (doc_v15 D15-7).
    pub(crate) fn next_frame_for_sample(
        &mut self,
        sample: usize,
        buf: &mut Vec<u8>,
    ) -> Result<Option<u64>, EngineError> {
        let mut header = [0_u8; 12];
        match &mut self.source {
            MsaFrameSource::Memory { bytes, cursor } => {
                if *cursor + 12 > bytes.len() {
                    return if *cursor == bytes.len() {
                        Ok(None)
                    } else {
                        Err(EngineError::job("MSA spill: truncated frame header"))
                    };
                }
                header.copy_from_slice(&bytes[*cursor..*cursor + 12]);
                *cursor += 12;
            }
            MsaFrameSource::Spilled { file } => match file.read_exact(&mut header) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
                Err(e) => return Err(EngineError::Io(e)),
            },
        }
        let pos0 = u64::from_le_bytes(header[..8].try_into().expect("invariant: 8 bytes"));
        let width = u32::from_le_bytes(header[8..].try_into().expect("invariant: 4 bytes")) as usize;
        let len = self.n_samples * width;
        let lo = sample * width;
        buf.clear();
        buf.resize(width, 0);
        match &mut self.source {
            MsaFrameSource::Memory { bytes, cursor } => {
                if *cursor + len > bytes.len() {
                    return Err(EngineError::job("MSA spill: truncated frame values"));
                }
                buf.copy_from_slice(&bytes[*cursor + lo..*cursor + lo + width]);
                *cursor += len;
            }
            MsaFrameSource::Spilled { file } => {
                file.seek_relative(lo as i64).map_err(EngineError::Io)?;
                file.read_exact(buf).map_err(EngineError::Io)?;
                let tail = u64::try_from(len - lo - width).unwrap_or(u64::MAX);
                file.seek_relative(tail as i64).map_err(EngineError::Io)?;
            }
        }
        Ok(Some(pos0))
    }
}

