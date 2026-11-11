//! Unified input entry (`docs/rounds/v2/03-architecture.md` §6): one open yields
//! the chunk source *and* the header — there is no second header-only
//! reopen. `-` (or `None`) reads stdin.
//!
//! stdin streams straight into the chunk sources (the format probe's prefix
//! is chained back), so a single-pass command starts computing before the
//! producer finishes (03-performance.md O5). Commands that must read the
//! input more than once (filter's passes) spill stdin once through
//! [`StdinSpill`], which deletes the temporary file on drop — on error and
//! early return as well.

use std::io::Read;
use std::path::Path;

use reed_io::IoError;
use reed_io::chunk::plain::PlainChunkSource;
use reed_io::chunk::region::{FileIndex, Region, RegionChunkSource};
use reed_io::chunk::{BgzfChunkSource, ChunkSource};
use reed_io::format::Format;
use reed_io::header::HeaderFacts;

/// An opened input: the chunk source plus the header facts that were parsed
/// while locating the first record.
///
/// The fields are crate-internal (doc_v11 N9): the transport types behind
/// them belong to the lower crates, so exposing them would freeze those
/// crates into the facade's stability promise. Callers read the header
/// through [`OpenedInput::facts`] and hand the value to a `run_*_opened`
/// entry point.
pub struct OpenedInput {
    /// Chunks from the first record on.
    pub(crate) source: Box<dyn ChunkSource>,
    /// Header facts (noodles header, samples, contigs).
    pub(crate) facts: HeaderFacts,
    /// Present when the input is BCF: the tag dictionaries the renderer
    /// needs (doc_v3 02-bcf.md B1).
    pub(crate) bcf: Option<std::sync::Arc<reed_io::bcf::BcfDict>>,
    /// The opened path (`None` for stdin). Output planning compares it with
    /// every path the run would write (v29 V29-03).
    pub(crate) path: Option<std::path::PathBuf>,
}

impl OpenedInput {
    /// The opened file's path, or `None` when the input is stdin.
    #[must_use]
    pub fn path(&self) -> Option<&Path> {
        self.path.as_deref()
    }

    /// Header facts (noodles header, samples, contigs) of the opened input.
    #[must_use]
    pub const fn facts(&self) -> &HeaderFacts {
        &self.facts
    }

    /// The tag dictionaries when the input is BCF (doc_v3 02-bcf.md B1).
    #[must_use]
    pub fn bcf_dict(&self) -> Option<&std::sync::Arc<reed_io::bcf::BcfDict>> {
        self.bcf.as_ref()
    }

    /// The underlying chunk source, for the crate's own integration tests
    /// that read raw chunks. Not part of the stable surface (N9).
    #[doc(hidden)]
    pub fn source_mut(&mut self) -> &mut dyn ChunkSource {
        self.source.as_mut()
    }
}

/// Opens a one-shot input: file at `path`, or stdin when `None`. stdin is
/// streamed, not spilled — the caller must read it exactly once.
///
/// # Errors
///
/// [`IoError`] when the stream ends before a valid `#CHROM` header.
pub fn open(path: Option<&Path>) -> Result<OpenedInput, IoError> {
    open_regions(path, &[])
}

/// Opens `p`, tagging a failure with the path so the operator sees which
/// file failed (doc_v7 D7-24: the plain `io::Error` conversion dropped it).
fn open_path(p: &Path) -> std::io::Result<std::fs::File> {
    std::fs::File::open(p)
        .map_err(|e| std::io::Error::new(e.kind(), format!("cannot open {}: {e}", p.display())))
}

/// Like `open`, but with an index-accelerated region restriction
/// (`-r/-R`): the reader seeks to the compressed range the index maps the
/// regions to. Records outside the regions are the caller's predicate to
/// drop (the same rule `-t` uses).
///
/// # Errors
///
/// [`IoError`] when the stream has no valid header, and region errors when
/// an index is required but missing.
pub fn open_regions(path: Option<&Path>, regions: &[Region]) -> Result<OpenedInput, IoError> {
    if regions.is_empty() {
        return open_all(path);
    }
    let Some(p) = path else {
        return Err(IoError::UnsupportedInput(
            "regions need a seekable file, not stdin",
        ));
    };
    use std::io::{Seek, SeekFrom};
    // doc_v28 V28-26: one open serves the sniff, the header scan and the
    // region reader (it used to open the same path three times, leaving two
    // windows for the file to change between them).
    let mut file = open_path(p)?;
    let (prefix, inflated) = sniff(&mut file)?;
    let format = Format::detect(&prefix, inflated.as_deref())?;
    file.seek(SeekFrom::Start(0))?;
    // The header scan runs on a dup of the handle, so `file`'s own offset
    // stays at 0 for the region reader's first seek.
    let (_, facts) = BgzfChunkSource::open(file.try_clone()?)?;
    // The index sits next to the input as `.csi` or `.tbi`.
    let index = sidecar_index(p)?;
    let source = RegionChunkSource::open(
        std::io::BufReader::with_capacity(1 << 20, file),
        &index,
        regions,
        &facts.contig_names,
    )?;
    let bcf = if format == Format::Bcf {
        Some(std::sync::Arc::new(reed_io::bcf::BcfDict::from_header(
            &facts.text,
        )?))
    } else {
        None
    };
    Ok(OpenedInput {
        source: Box::new(source),
        facts,
        bcf,
        path: Some(p.to_path_buf()),
    })
}

/// Opens the `.csi` or `.tbi` next to `p` — the one place the sidecar rule
/// lives (doc_v28 V28-26).
fn sidecar_index(p: &Path) -> Result<FileIndex, IoError> {
    for suffix in [".csi", ".tbi"] {
        let mut sidecar = p.as_os_str().to_os_string();
        sidecar.push(suffix);
        let sidecar = Path::new(&sidecar);
        if sidecar.exists() {
            return FileIndex::open(sidecar);
        }
    }
    Err(IoError::UnsupportedInput(
        "regions need an index next to the input (.csi or .tbi; build one with `reed index`)",
    ))
}

/// Detects a file's input format from its magic bytes without opening the
/// whole stream (used to reject unindexable inputs up front, doc_v7 D7-3).
///
/// # Errors
///
/// [`IoError`] when the file cannot be read or no known magic matches.
pub fn detect_format(path: &Path) -> Result<Format, IoError> {
    let mut probe = open_path(path)?;
    let (prefix, inflated) = sniff(&mut probe)?;
    Format::detect(&prefix, inflated.as_deref())
}

/// Reads enough of the stream to detect its format: the whole first BGZF
/// block when the bytes are BGZF (the inflated content tells BGZF-VCF from
/// BGZF-BCF), else a short prefix. The raw bytes are returned so the caller
/// can chain them back onto the stream.
fn sniff(reader: &mut impl Read) -> Result<(Vec<u8>, Option<Vec<u8>>), IoError> {
    let mut header = [0_u8; reed_io::bgzf::block::BGZF_HEADER_LEN];
    let mut filled = 0;
    while filled < header.len() {
        let n = reader.read(&mut header[filled..])?;
        if n == 0 {
            break;
        }
        filled += n;
    }
    if filled == 0 {
        return Err(IoError::UnsupportedInput("input is empty"));
    }
    if let Some(total) = reed_io::bgzf_total_size(&header[..filled])
        && total > filled
    {
        let mut block = header[..filled].to_vec();
        block.resize(total, 0);
        reader.read_exact(&mut block[filled..])?;
        let mut inflater = reed_io::bgzf::inflate::Inflater::new();
        let mut text = Vec::new();
        inflater.inflate_block(&block, &mut text)?;
        return Ok((block, Some(text)));
    }
    // Not (or not yet known to be) a full BGZF member: a short text prefix
    // is enough to tell VCF text from gzip.
    let mut prefix = header[..filled].to_vec();
    let mut extra = [0_u8; 64];
    let n = reader.read(&mut extra)?;
    prefix.extend_from_slice(&extra[..n]);
    Ok((prefix, None))
}

/// A stdin stream spilled to a temporary file so a multi-pass driver can
/// re-read the same bytes. The file is removed when this guard drops —
/// A stdin replay file: the stream is copied into an exclusive temp file so
/// the driver can open it twice (regions/replay); the [`TempSpill`] removes
/// it on every exit path (03-performance.md O5/04 §C6).
pub struct StdinSpill {
    spill: crate::spill::TempSpill,
}

impl StdinSpill {
    /// Spills stdin into `tmp_dir` (default: the platform temp directory).
    ///
    /// # Errors
    ///
    /// I/O failures creating or filling the temporary file.
    pub fn new(tmp_dir: Option<&Path>) -> Result<Self, IoError> {
        let (mut out, spill) = crate::spill::TempSpill::create(tmp_dir, "stdin", "vcf")?;
        std::io::copy(&mut std::io::stdin(), &mut out)?;
        out.sync_all()?;
        Ok(Self { spill })
    }

    /// The temporary file holding the spilled bytes.
    #[must_use]
    pub fn path(&self) -> &Path {
        self.spill.path()
    }
}

/// Opens the whole stream (no region restriction), detecting the format by
/// magic: BGZF VCF, plain VCF, plain gzip, BCF (03-architecture.md §6).
/// stdin is wrapped directly, with the probe's bytes chained back.
fn open_all(path: Option<&Path>) -> Result<OpenedInput, IoError> {
    match path {
        Some(p) => open_file(p),
        None => open_stdin(),
    }
}

/// Opens a file once: sniffs its format, rewinds, and hands the same handle
/// to the source the format names (doc_v28 V28-26).
///
/// The previous version opened the path twice — once to sniff and once for
/// the source — which also left a window for the file to change between the
/// two opens. One handle serves both.
fn open_file(p: &Path) -> Result<OpenedInput, IoError> {
    use std::io::{Seek, SeekFrom};
    let mut file = open_path(p)?;
    let (prefix, inflated) = sniff(&mut file)?;
    let format = Format::detect(&prefix, inflated.as_deref())?;
    file.seek(SeekFrom::Start(0))?;
    build_source(format, file, Some(p.to_path_buf()))
}

/// The four format arms, written once for files and stdin alike.
///
/// `reader` yields the stream from its first byte; the caller owns the sniff
/// and the rewind (a file seeks, stdin chains the probed prefix back).
fn build_source(
    format: Format,
    reader: impl Read + Send + 'static,
    path: Option<std::path::PathBuf>,
) -> Result<OpenedInput, IoError> {
    match format {
        Format::VcfBgzf => {
            let (source, facts) = BgzfChunkSource::open(reader)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: None,
                path,
            })
        }
        Format::Vcf => {
            let (source, facts) = PlainChunkSource::open_vcf(reader)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: None,
                path,
            })
        }
        Format::VcfGzip => {
            warn_serial_gzip();
            let (source, facts) = PlainChunkSource::open_gzip(reader)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: None,
                path,
            })
        }
        Format::Bcf => {
            // BCF shares the BGZF chunk source; the dictionaries come from
            // the header text embedded in the file.
            let (source, facts) = BgzfChunkSource::open(reader)?;
            let dict = reed_io::bcf::BcfDict::from_header(&facts.text)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: Some(std::sync::Arc::new(dict)),
                path,
            })
        }
    }
}

fn open_stdin() -> Result<OpenedInput, IoError> {
    let mut stdin = std::io::stdin();
    let (prefix, inflated) = sniff(&mut stdin)?;
    let format = Format::detect(&prefix, inflated.as_deref())?;
    // stdin cannot seek, so the probed bytes are chained back on the front.
    build_source(format, std::io::Cursor::new(prefix).chain(stdin), None)
}

/// Opens a BGZF file with the streaming reader, not a memory map.
///
/// This used to prefer the memory-mapped source (doc_v4 03 §3.4 step 3 — no
/// copies, no read syscalls). `memmap2::Mmap::map` is only sound while the
/// backing file is not modified by *any* process for the lifetime of the
/// mapping, and "reed opened it read-only" does not establish that: a shared
/// working directory, an NFS export or a concurrent in-place rewrite all
/// violate it, with silent corruption or SIGBUS as the failure mode
/// (doc_v26 V26-00). The default is therefore the streaming reader, which
/// reads into a private buffer that cannot be invalidated underneath the run.
///
/// [`reed_io::chunk::bgzf::MmapChunkSource`] remains available to callers
/// that can prove an immutable snapshot; `reed`'s own commands do not use it.
fn warn_serial_gzip() {
    tracing::warn!(
        "the input is plain gzip, not BGZF; decompression is serial \
         (repack with `bgzip` for parallel reads)"
    );
}
