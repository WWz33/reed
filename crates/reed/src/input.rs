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
use std::path::{Path, PathBuf};

use reed_io::IoError;
use reed_io::chunk::plain::PlainChunkSource;
use reed_io::chunk::region::{FileIndex, Region, RegionChunkSource};
use reed_io::chunk::{BgzfChunkSource, ChunkSource};
use reed_io::format::Format;
use reed_io::header::HeaderFacts;

/// An opened input: the chunk source plus the header facts that were parsed
/// while locating the first record.
pub struct OpenedInput {
    /// Chunks from the first record on.
    pub source: Box<dyn ChunkSource>,
    /// Header facts (noodles header, samples, contigs).
    pub facts: HeaderFacts,
    /// Present when the input is BCF: the tag dictionaries the renderer
    /// needs (doc_v3 02-bcf.md B1).
    pub bcf: Option<std::sync::Arc<reed_io::bcf::BcfDict>>,
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
    let format = {
        let mut probe = std::fs::File::open(p)?;
        let (prefix, inflated) = sniff(&mut probe)?;
        Format::detect(&prefix, inflated.as_deref())?
    };
    let file = std::fs::File::open(p)?;
    let (_, facts) = BgzfChunkSource::open(file)?;
    // The index sits next to the input as `.csi` or `.tbi`.
    let mut csi = p.as_os_str().to_os_string();
    csi.push(".csi");
    let mut tbi = p.as_os_str().to_os_string();
    tbi.push(".tbi");
    let index = if Path::new(&csi).exists() {
        FileIndex::open(Path::new(&csi))?
    } else if Path::new(&tbi).exists() {
        FileIndex::open(Path::new(&tbi))?
    } else {
        return Err(IoError::UnsupportedInput(
            "regions need an index next to the input (.csi or .tbi; build one with `reed index`)",
        ));
    };
    let file = std::fs::File::open(p)?;
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
    })
}

/// Detects a file's input format from its magic bytes without opening the
/// whole stream (used to reject unindexable inputs up front, doc_v7 D7-3).
///
/// # Errors
///
/// [`IoError`] when the file cannot be read or no known magic matches.
pub fn detect_format(path: &Path) -> Result<Format, IoError> {
    let mut probe = std::fs::File::open(path)?;
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
/// normal exit, error returns and early returns included (03-performance.md
/// O5/04 §C6).
pub struct StdinSpill {
    path: PathBuf,
}

impl StdinSpill {
    /// Spills stdin into `tmp_dir` (default: the platform temp directory).
    ///
    /// # Errors
    ///
    /// I/O failures creating or filling the temporary file.
    pub fn new(tmp_dir: Option<&Path>) -> Result<Self, IoError> {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let dir = tmp_dir.map_or_else(std::env::temp_dir, Path::to_path_buf);
        let path = dir.join(format!(
            "reed-stdin-{}-{}.vcf",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        let mut out = std::fs::File::create(&path)?;
        std::io::copy(&mut std::io::stdin(), &mut out)?;
        out.sync_all()?;
        Ok(Self { path })
    }

    /// The temporary file holding the spilled bytes.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for StdinSpill {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
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

fn open_file(p: &Path) -> Result<OpenedInput, IoError> {
    let (prefix, inflated) = {
        let mut file = std::fs::File::open(p)?;
        sniff(&mut file)?
    };
    match Format::detect(&prefix, inflated.as_deref())? {
        Format::VcfBgzf => {
            let (source, facts) = open_bgzf(p)?;
            Ok(OpenedInput {
                source,
                facts,
                bcf: None,
            })
        }
        Format::Vcf => {
            let file = std::fs::File::open(p)?;
            let (source, facts) = PlainChunkSource::open_vcf(file)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: None,
            })
        }
        Format::VcfGzip => {
            let file = std::fs::File::open(p)?;
            warn_serial_gzip();
            let (source, facts) = PlainChunkSource::open_gzip(file)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: None,
            })
        }
        Format::Bcf => {
            // BCF shares the BGZF chunk source; the dictionaries come from
            // the header text embedded in the file.
            let (source, facts) = open_bgzf(p)?;
            let dict = reed_io::bcf::BcfDict::from_header(&facts.text)?;
            Ok(OpenedInput {
                source,
                facts,
                bcf: Some(std::sync::Arc::new(dict)),
            })
        }
    }
}

/// Opens a BGZF file: memory-mapped when possible (doc_v4 03 §3.4 step 3 —
/// no copies, no read syscalls), else the streaming reader. Both cut chunks
/// by the same rule, so the choice cannot change the output.
fn open_bgzf(p: &Path) -> Result<(Box<dyn ChunkSource>, HeaderFacts), IoError> {
    match reed_io::chunk::bgzf::MmapChunkSource::open(std::fs::File::open(p)?) {
        Ok((source, facts)) => Ok((Box::new(source), facts)),
        Err(_) => {
            let (source, facts) = BgzfChunkSource::open(std::fs::File::open(p)?)?;
            Ok((Box::new(source), facts))
        }
    }
}

/// Streams stdin: the format probe's prefix is chained back so the source
/// sees every byte exactly once.
fn open_stdin() -> Result<OpenedInput, IoError> {
    let mut stdin = std::io::stdin();
    let (prefix, inflated) = sniff(&mut stdin)?;
    let format = Format::detect(&prefix, inflated.as_deref())?;
    let reader = std::io::Cursor::new(prefix).chain(stdin);
    match format {
        Format::VcfBgzf => {
            let (source, facts) = BgzfChunkSource::open(reader)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: None,
            })
        }
        Format::Vcf => {
            let (source, facts) = PlainChunkSource::open_vcf(reader)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: None,
            })
        }
        Format::VcfGzip => {
            warn_serial_gzip();
            let (source, facts) = PlainChunkSource::open_gzip(reader)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: None,
            })
        }
        Format::Bcf => {
            let (source, facts) = BgzfChunkSource::open(reader)?;
            let dict = reed_io::bcf::BcfDict::from_header(&facts.text)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: Some(std::sync::Arc::new(dict)),
            })
        }
    }
}

fn warn_serial_gzip() {
    tracing::warn!(
        "the input is plain gzip, not BGZF; decompression is serial \
         (repack with `bgzip` for parallel reads)"
    );
}
