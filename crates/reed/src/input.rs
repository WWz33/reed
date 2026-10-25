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
use std::sync::atomic::{AtomicU64, Ordering};

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
}

impl OpenedInput {
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
    let format = {
        let mut probe = open_path(p)?;
        let (prefix, inflated) = sniff(&mut probe)?;
        Format::detect(&prefix, inflated.as_deref())?
    };
    let file = open_path(p)?;
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
    let file = open_path(p)?;
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
/// normal exit, error returns and early returns included (03-performance.md
/// O5/04 §C6).
pub struct StdinSpill {
    path: PathBuf,
}

/// Process-wide spill-name counter (doc_v20 SEC-1): one sequence per
/// process, starting at 0, so spill names are predictable for operators.
static SPILL_COUNTER: AtomicU64 = AtomicU64::new(0);

/// The candidate spill path `counter` generations on in `dir`
/// (`reed-stdin-{pid}-{counter}.vcf`, doc_v20 SEC-1).
fn spill_path(dir: &Path, counter: u64) -> PathBuf {
    dir.join(format!(
        "reed-stdin-{}-{}.vcf",
        std::process::id(),
        counter
    ))
}

/// How many candidate names one spill may burn through before giving up
/// (doc_v20 SEC-1).
const SPILL_ATTEMPTS: usize = 16;

/// Creates the next spill file in `dir`, exclusively (doc_v20 SEC-1).
///
/// The name stays predictable (`pid` + process counter) so operators can
/// recognise spill files, but creation uses O_EXCL (`create_new`), which
/// fails on an existing file *and* on a symlink — no follow, no truncate.
/// Exclusivity alone defeats a pre-planted name on a shared host; random
/// names are not needed. On `AlreadyExists` the counter is bumped and the
/// next candidate tried, up to [`SPILL_ATTEMPTS`] names in all.
fn create_spill_file(dir: &Path) -> std::io::Result<(std::fs::File, PathBuf)> {
    for _ in 0..SPILL_ATTEMPTS {
        let path = spill_path(dir, SPILL_COUNTER.fetch_add(1, Ordering::Relaxed));
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
        {
            Ok(file) => return Ok((file, path)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(e) => return Err(e),
        }
    }
    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        format!(
            "cannot create a spill file in {}: all {} candidate names are taken",
            dir.display(),
            SPILL_ATTEMPTS
        ),
    ))
}

impl StdinSpill {
    /// Spills stdin into `tmp_dir` (default: the platform temp directory).
    ///
    /// # Errors
    ///
    /// I/O failures creating or filling the temporary file.
    pub fn new(tmp_dir: Option<&Path>) -> Result<Self, IoError> {
        let dir = tmp_dir.map_or_else(std::env::temp_dir, Path::to_path_buf);
        let (mut out, path) = create_spill_file(&dir)?;
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
        // doc_v20 SEC-1: `self.path` is the name `create_new` made ours.
        // `remove_file` is unlink(2), which never follows symlinks — even if
        // the name were swapped for a link after creation, this removes the
        // link entry, not a victim target.
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
        let mut file = open_path(p)?;
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
            let file = open_path(p)?;
            let (source, facts) = PlainChunkSource::open_vcf(file)?;
            Ok(OpenedInput {
                source: Box::new(source),
                facts,
                bcf: None,
            })
        }
        Format::VcfGzip => {
            let file = open_path(p)?;
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
    match reed_io::chunk::bgzf::MmapChunkSource::open(open_path(p)?) {
        Ok((source, facts)) => Ok((Box::new(source), facts)),
        Err(_) => {
            let (source, facts) = BgzfChunkSource::open(open_path(p)?)?;
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

#[cfg(test)]
mod tests {
    use super::{SPILL_COUNTER, create_spill_file, spill_path};
    use std::path::PathBuf;
    use std::sync::atomic::Ordering;

    /// The sec1 tests predict the next [`SPILL_COUNTER`] value, so they must
    /// not allocate concurrently (cargo test runs tests on parallel threads).
    static SEC1_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    /// A private scratch dir, mirroring the stdin_streaming.rs tests' style.
    fn sec1_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("reed-sec1-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("invariant: scratch dir");
        dir
    }

    /// doc_v20 SEC-1: the spill name is predictable (`temp_dir` + pid +
    /// process counter), so on a shared host another user can pre-create a
    /// file at the exact path the spill will use. The spill must refuse the
    /// name and leave the victim's bytes untouched (it used to truncate).
    #[test]
    fn sec1_precreated_file_is_not_truncated() {
        let _serialized = SEC1_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let dir = sec1_dir("file");
        let victim = spill_path(&dir, SPILL_COUNTER.load(Ordering::Relaxed));
        std::fs::write(&victim, b"victim bytes").expect("invariant: plant victim");
        let (_spill, path) = create_spill_file(&dir).expect("invariant: create spill");
        assert_eq!(
            std::fs::read(&victim).expect("invariant: re-read victim"),
            b"victim bytes".as_slice(),
            "the pre-created victim must be untouched"
        );
        assert_ne!(
            path, victim,
            "the spill must refuse the pre-created name and pick another"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// doc_v20 SEC-1: same attack with a symlink planted at the predicted
    /// path, pointing at a victim file elsewhere in the dir. The spill must
    /// not follow the link; the victim's bytes must survive.
    #[cfg(unix)]
    #[test]
    fn sec1_precreated_symlink_is_not_followed() {
        let _serialized = SEC1_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let dir = sec1_dir("symlink");
        let victim = dir.join("victim.vcf");
        std::fs::write(&victim, b"victim bytes").expect("invariant: plant victim");
        let predicted = spill_path(&dir, SPILL_COUNTER.load(Ordering::Relaxed));
        std::os::unix::fs::symlink(&victim, &predicted).expect("invariant: plant symlink");
        let (_spill, path) = create_spill_file(&dir).expect("invariant: create spill");
        assert_eq!(
            std::fs::read(&victim).expect("invariant: re-read victim"),
            b"victim bytes".as_slice(),
            "the symlink's target must be untouched"
        );
        assert_ne!(path, predicted, "the spill must not write through the link");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// doc_v20 SEC-1: refusals retry with bumped counters, bounded — when
    /// every candidate name is taken the spill is an error, not a clobber,
    /// and none of the planted files was touched.
    #[test]
    fn sec1_exhausted_candidates_error_without_clobbering() {
        let _serialized = SEC1_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let dir = sec1_dir("exhausted");
        let first = SPILL_COUNTER.load(Ordering::Relaxed);
        for n in first..first + super::SPILL_ATTEMPTS as u64 {
            std::fs::write(spill_path(&dir, n), b"victim").expect("invariant: plant candidate");
        }
        let err = create_spill_file(&dir).expect_err("every candidate name is taken");
        assert_eq!(err.kind(), std::io::ErrorKind::AlreadyExists);
        for n in first..first + super::SPILL_ATTEMPTS as u64 {
            assert_eq!(
                std::fs::read(spill_path(&dir, n)).expect("invariant: re-read candidate"),
                b"victim".as_slice(),
                "candidate {n} must not have been touched"
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}
