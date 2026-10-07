//! Unified input entry (`doc_v2/03-architecture.md` §6): one open yields
//! the chunk source *and* the header — there is no second header-only
//! reopen. `-` (or `None`) reads BGZF VCF from stdin.

use std::path::{Path, PathBuf};

use reed_io::chunk::plain::PlainChunkSource;
use reed_io::chunk::region::{FileIndex, IndexQuery as _, Region, RegionChunkSource};
use reed_io::chunk::{BgzfChunkSource, ChunkSource};
use reed_io::format::Format;
use reed_io::header::HeaderFacts;
use reed_io::IoError;

/// Opens a BGZF VCF stream: file at `path`, or stdin when `None`. The
/// returned source yields record chunks from the first record on; the facts
/// carry the header text and the parsed noodles header.
///
/// # Errors
///
/// [`IoError`] when the stream ends before a valid `#CHROM` header.
pub fn open(path: Option<&Path>) -> Result<(Box<dyn ChunkSource>, HeaderFacts), IoError> {
    open_regions(path, &[])
}

/// Like [`open`], but with an index-accelerated region restriction
/// (`-r/-R`): the reader seeks to the compressed range the index maps the
/// regions to. Records outside the regions are the caller's predicate to
/// drop (the same rule `-t` uses).
///
/// # Errors
///
/// [`IoError`] when the stream has no valid header, and region errors when
/// an index is required but missing.
pub fn open_regions(
    path: Option<&Path>,
    regions: &[Region],
) -> Result<(Box<dyn ChunkSource>, HeaderFacts), IoError> {
    if regions.is_empty() {
        return open_all(path);
    }
    let Some(p) = path else {
        return Err(IoError::UnsupportedInput("regions need a seekable file, not stdin"));
    };
    let file = std::fs::File::open(p)?;
    let (_, facts) =
        BgzfChunkSource::open(std::io::BufReader::with_capacity(1 << 20, file))?;
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
    Ok((Box::new(source), facts))
}

/// Reads enough of the stream to detect its format: the whole first BGZF
/// block when the bytes are BGZF (the inflated content tells BGZF-VCF from
/// BGZF-BCF), else a short prefix. The bytes are returned so the caller can
/// chain them back onto the stream.
fn sniff(reader: &mut impl std::io::Read) -> Result<(Vec<u8>, Option<Vec<u8>>), IoError> {
    use std::io::Read as _;
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

/// The temp file a stdin stream is spilled to, so the drivers' several
/// opens (header, then one per pass) all see the same bytes. Stdin cannot
/// be re-read; a file can.
fn stdin_spill() -> Result<&'static Path, IoError> {
    use std::sync::OnceLock;
    static SPILL: OnceLock<PathBuf> = OnceLock::new();
    if let Some(path) = SPILL.get() {
        return Ok(path);
    }
    let mut path = std::env::temp_dir();
    path.push(format!("reed-stdin-{}.vcf", std::process::id()));
    {
        use std::io::Write as _;
        let mut out = std::fs::File::create(&path)?;
        std::io::copy(&mut std::io::stdin(), &mut out)?;
    }
    Ok(SPILL.get_or_init(|| path))
}

/// Opens the whole stream (no region restriction), detecting the format by
/// magic: BGZF VCF, plain VCF, plain gzip, BCF (03-architecture.md §6).
fn open_all(path: Option<&Path>) -> Result<(Box<dyn ChunkSource>, HeaderFacts), IoError> {
    let spilled;
    let path = match path {
        Some(p) => Some(p),
        None => {
            spilled = stdin_spill()?;
            Some(spilled)
        }
    };
    match path {
        Some(p) => {
            let (prefix, inflated) = {
                let mut file = std::fs::File::open(p)?;
                sniff(&mut file)?
            };
            match Format::detect(&prefix, inflated.as_deref())? {
                Format::VcfBgzf => {
                    let file = std::fs::File::open(p)?;
                    let (source, facts) = BgzfChunkSource::open(
                        std::io::BufReader::with_capacity(1 << 20, file),
                    )?;
                    Ok((Box::new(source), facts))
                }
                Format::Vcf => {
                    let file = std::fs::File::open(p)?;
                    let (source, facts) = PlainChunkSource::open_vcf(file)?;
                    Ok((Box::new(source), facts))
                }
                Format::VcfGzip => {
                    let file = std::fs::File::open(p)?;
                    eprintln!(
                        "reed: warning: the input is plain gzip, not BGZF; decompression is \
                         serial (repack with `bgzip` for parallel reads)"
                    );
                    let (source, facts) = PlainChunkSource::open_gzip(file)?;
                    Ok((Box::new(source), facts))
                }
                Format::Bcf => Err(IoError::UnsupportedInput(
                    "BCF input is not supported yet; convert with `bcftools view -Oz`",
                )),
            }
        }
        // `path` is always `Some` here (stdin was spilled above).
        None => unreachable!("invariant: stdin is spilled to a file first"),
    }
}
