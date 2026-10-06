//! Unified input entry (`doc_v2/03-architecture.md` §6): one open yields
//! the chunk source *and* the header — there is no second header-only
//! reopen. `-` (or `None`) reads BGZF VCF from stdin.

use std::path::Path;

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
