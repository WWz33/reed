//! Unified input entry (`doc_v2/03-architecture.md` §6): one open yields
//! the chunk source *and* the header — there is no second header-only
//! reopen. `-` (or `None`) reads BGZF VCF from stdin.

use std::path::Path;

use reed_io::chunk::region::{FileIndex, IndexQuery as _, Region, RegionChunkSource};
use reed_io::chunk::{BgzfChunkSource, ChunkSource};
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

/// Opens the whole stream (no region restriction).
fn open_all(path: Option<&Path>) -> Result<(Box<dyn ChunkSource>, HeaderFacts), IoError> {
    match path {
        Some(p) => {
            let file = std::fs::File::open(p)?;
            let (source, facts) = BgzfChunkSource::open(std::io::BufReader::with_capacity(
                1 << 20,
                file,
            ))?;
            Ok((Box::new(source), facts))
        }
        // `Stdin` is `Read + Send` (internally synchronized) and `Clone`.
        None => {
            let (source, facts) = BgzfChunkSource::open(std::io::BufReader::with_capacity(
                1 << 20,
                std::io::stdin(),
            ))?;
            Ok((Box::new(source), facts))
        }
    }
}
