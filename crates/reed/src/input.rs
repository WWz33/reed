//! Unified input entry (`docs/rounds/v2/03-architecture.md` §6): one open yields
//! the chunk source *and* the header — there is no second header-only
//! reopen. `-` (or `None`) reads stdin.
//!
//! Capacities are types (doc_v30 R05–R07): a [`StreamedInput`] is a stream
//! that can be consumed once; a [`ReplayableInput`] is a file (or a stdin
//! spill) that can be opened again; a [`SeekableInput`] is a replayable input
//! whose index and BGZF framing were verified, so `-r/-R` can seek it. The
//! driver asks these capabilities instead of remembering a bare
//! `path: Option<PathBuf>`, and the one `open` count is visible in
//! [`PendingInput`], which holds the already-opened file between the header
//! probe and the source construction.
//!
//! stdin streams straight into the chunk sources (the format probe's prefix
//! is chained back), so a single-pass command starts computing before the
//! producer finishes (03-performance.md O5). Commands that must read the
//! input more than once ask [`ReplayableInput::from_stdin`] for a spill; the
//! [`TempSpill`](crate::spill::TempSpill) removes it on drop — on error and
//! early return as well.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use reed_io::IoError;
use reed_io::bcf::BcfDict;
use reed_io::chunk::plain::PlainChunkSource;
use reed_io::chunk::region::{FileIndex, Region, RegionChunkSource};
use reed_io::chunk::{BgzfChunkSource, ChunkSource};
use reed_io::format::Format;
use reed_io::header::HeaderFacts;

/// How the bytes are framed on the way in (the transport half of the
/// capability, doc_v30 R07).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Framing {
    /// BGZF-compressed VCF or BCF.
    Bgzf,
    /// A regular gzip stream (serial decompression).
    Gzip,
    /// Uncompressed VCF text.
    Plain,
}

impl Framing {
    /// The framing a detected [`Format`] uses.
    #[must_use]
    pub const fn of(format: Format) -> Self {
        match format {
            Format::VcfBgzf | Format::Bcf => Self::Bgzf,
            Format::VcfGzip => Self::Gzip,
            Format::Vcf => Self::Plain,
        }
    }
}

/// Format, framing and the native dictionaries, carried as one value
/// (doc_v30 R07): consumers ask [`Self::is_bcf`] instead of inferring BCF
/// from "the dictionary field happens to be `Some`".
#[derive(Clone)]
pub struct InputMeta {
    format: Format,
    framing: Framing,
    bcf_dict: Option<Arc<BcfDict>>,
}

impl InputMeta {
    fn new(format: Format, header_text: &[u8]) -> Result<Self, IoError> {
        let bcf_dict = if format == Format::Bcf {
            Some(Arc::new(BcfDict::from_header(header_text)?))
        } else {
            None
        };
        Ok(Self {
            format,
            framing: Framing::of(format),
            bcf_dict,
        })
    }

    /// The detected input format.
    #[must_use]
    pub const fn format(&self) -> Format {
        self.format
    }

    /// How the input is framed.
    #[must_use]
    pub const fn framing(&self) -> Framing {
        self.framing
    }

    /// Whether the input is BCF.
    #[must_use]
    pub const fn is_bcf(&self) -> bool {
        matches!(self.format, Format::Bcf)
    }

    /// The tag dictionaries a BCF input carries.
    #[must_use]
    pub fn bcf_dict(&self) -> Option<&Arc<BcfDict>> {
        self.bcf_dict.as_ref()
    }
}

/// A stream that can be consumed exactly once (stdin, or a file read once).
///
/// The fields stay crate-internal (doc_v11 N9): the transport types behind
/// them belong to the lower crates, so exposing them would freeze those
/// crates into the facade's stability promise. [`Self::into_parts`] is the
/// one way to take the source out.
pub struct StreamedInput {
    source: Box<dyn ChunkSource>,
    facts: HeaderFacts,
    meta: InputMeta,
    path: Option<PathBuf>,
    capability: InputCapability,
}

/// The historical name of [`StreamedInput`] (doc_v30 R05 keeps the public
/// entry-point type name stable while the capability types are new).
pub type OpenedInput = StreamedInput;

impl StreamedInput {
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

    /// Format, framing and native dictionaries.
    #[must_use]
    pub const fn meta(&self) -> &InputMeta {
        &self.meta
    }

    /// The detected format.
    #[must_use]
    pub const fn format(&self) -> Format {
        self.meta.format()
    }

    /// Whether the input is BCF (one predicate, not a field probe).
    #[must_use]
    pub const fn is_bcf(&self) -> bool {
        self.meta.is_bcf()
    }

    /// The tag dictionaries when the input is BCF (doc_v3 02-bcf.md B1).
    #[must_use]
    pub fn bcf_dict(&self) -> Option<&Arc<BcfDict>> {
        self.meta.bcf_dict()
    }

    /// Whether the input can be opened again (file, or stdin spilled by the
    /// caller before opening).
    #[must_use]
    pub const fn can_replay(&self) -> bool {
        !matches!(self.capability, InputCapability::OneShot)
    }

    /// Takes the replay/seek capability out, leaving a one-shot stream.
    ///
    /// Callers take it before consuming the stream so a second pass can be
    /// opened without remembering a path of their own.
    pub fn take_capability(&mut self) -> InputCapability {
        std::mem::take(&mut self.capability)
    }

    /// The one way to take the source (doc_v30 R07).
    #[must_use]
    pub fn into_parts(self) -> (Box<dyn ChunkSource>, HeaderFacts, InputMeta) {
        (self.source, self.facts, self.meta)
    }
}

/// A file-backed input that can be opened again, possibly after a stdin
/// spill. Holding this value keeps a spilled temporary file alive.
pub struct ReplayableInput {
    path: PathBuf,
    /// Keeps a spilled stdin file alive for the capability's lifetime.
    _spill: Option<StdinSpill>,
}

impl ReplayableInput {
    /// The file every replay opens.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Spills stdin into `tmp_dir` (default: the platform temp directory)
    /// and returns the replay capability for it.
    ///
    /// The caller must not have opened stdin already: this copies the raw
    /// stream, which is why the spill happens before the first open.
    ///
    /// # Errors
    ///
    /// I/O failures creating or filling the temporary file.
    pub fn from_stdin(tmp_dir: Option<&Path>) -> Result<Self, IoError> {
        let spill = StdinSpill::new(tmp_dir)?;
        Ok(Self {
            path: spill.path().to_path_buf(),
            _spill: Some(spill),
        })
    }

    /// Opens one stream over the same bytes, restricted to `regions` when
    /// non-empty. The returned stream is one-shot; the caller keeps this
    /// capability for the next pass.
    ///
    /// # Errors
    ///
    /// Transport failures opening the file, and [`IoError::Capability`] when
    /// regions need an index the input does not have.
    pub fn replay(&self, regions: &[Region]) -> Result<StreamedInput, IoError> {
        let mut opened = open_regions(Some(&self.path), regions)?;
        opened.capability = InputCapability::OneShot;
        Ok(opened)
    }
}

/// A replayable input whose BGZF/BCF framing and sidecar index were verified
/// when [`open_regions`] built it, so `-r/-R` can seek it (doc_v30 R05/R06).
pub struct SeekableInput {
    replay: ReplayableInput,
}

impl SeekableInput {
    /// The file every seek opens.
    #[must_use]
    pub fn path(&self) -> &Path {
        self.replay.path()
    }

    /// The replay capability behind the seek.
    #[must_use]
    pub const fn replayable(&self) -> &ReplayableInput {
        &self.replay
    }

    /// Opens one region-restricted stream over the same bytes.
    ///
    /// # Errors
    ///
    /// Transport failures opening the file.
    pub fn seek(&self, regions: &[Region]) -> Result<StreamedInput, IoError> {
        self.replay.replay(regions)
    }
}

/// What the opened input can do besides being read once (doc_v30 R05/R07).
#[derive(Default)]
pub enum InputCapability {
    /// One stream; a second pass is a usage error.
    #[default]
    OneShot,
    /// A plain file (or stdin spill) that can be opened again.
    Replayable(ReplayableInput),
    /// A file with an index and BGZF framing, seekable by region.
    Seekable(SeekableInput),
}

impl InputCapability {
    /// Whether a second pass can be opened.
    #[must_use]
    pub const fn can_replay(&self) -> bool {
        !matches!(self, Self::OneShot)
    }

    /// The replay capability, when there is one.
    #[must_use]
    pub const fn replayable(&self) -> Option<&ReplayableInput> {
        match self {
            Self::OneShot => None,
            Self::Replayable(replay) => Some(replay),
            Self::Seekable(seek) => Some(seek.replayable()),
        }
    }

    /// Opens one stream, restricted to `regions` when non-empty.
    ///
    /// # Errors
    ///
    /// [`IoError::Capability`] with [`reed_io::Needed::Replay`] when the
    /// input is one-shot; transport failures otherwise.
    pub fn replay(&self, regions: &[Region]) -> Result<StreamedInput, IoError> {
        match self {
            Self::OneShot => Err(IoError::Capability {
                needed: reed_io::Needed::Replay,
                fix: "write the stream to a file and pass that path, or read stdin once",
            }),
            Self::Replayable(replay) => replay.replay(regions),
            Self::Seekable(seek) => seek.seek(regions),
        }
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
/// The format and the sidecar are checked before the region reader is
/// built, so a plain VCF, a regular gzip stream or a missing index is a
/// capability error ([`reed_io::Needed`]) rather than a malformed-BGZF
/// surprise (doc_v30 R06).
///
/// # Errors
///
/// [`IoError`] when the stream has no valid header, and
/// [`IoError::Capability`] when regions need a capability the input lacks.
pub fn open_regions(path: Option<&Path>, regions: &[Region]) -> Result<OpenedInput, IoError> {
    match path {
        Some(p) => PendingInput::open(p)?.into_opened(regions),
        None if regions.is_empty() => open_stdin(),
        None => Err(IoError::Capability {
            needed: reed_io::Needed::Seek,
            fix: "regions need a file path; stdin cannot be seeked",
        }),
    }
}

/// An opened file that has been probed (format plus header) but not yet
/// framed into a chunk source.
///
/// The driver needs the header's contig vocabulary before it can parse
/// `-r/-R` text, and the parsed regions decide how the same bytes are
/// framed. Keeping the file handle here means the probe and the source
/// share one `open` (doc_v30 R05/R06); the previous shape opened the path
/// once to read the header and again for the seek.
pub struct PendingInput {
    file: std::fs::File,
    format: Format,
    facts: HeaderFacts,
    path: PathBuf,
}

impl PendingInput {
    /// Opens `path` once, detects the format and reads the header.
    ///
    /// # Errors
    ///
    /// [`IoError`] when the file cannot be read, the format is unknown or
    /// the header is missing.
    pub fn open(path: &Path) -> Result<Self, IoError> {
        use std::io::{Seek, SeekFrom};
        let mut file = open_path(path)?;
        let (prefix, inflated) = sniff(&mut file)?;
        let format = Format::detect(&prefix, inflated.as_deref())?;
        file.seek(SeekFrom::Start(0))?;
        let facts = read_facts(format, file.try_clone()?)?;
        // `try_clone` shares the file offset: the header scan leaves the
        // original at EOF, so rewind before the source reads it.
        file.seek(SeekFrom::Start(0))?;
        Ok(Self {
            file,
            format,
            facts,
            path: path.to_path_buf(),
        })
    }

    /// Header facts of the probed input.
    #[must_use]
    pub const fn facts(&self) -> &HeaderFacts {
        &self.facts
    }

    /// Frames the already-opened file into a stream, seeking to `regions`
    /// when non-empty.
    ///
    /// # Errors
    ///
    /// [`IoError::Capability`] when `regions` needs BGZF framing or an index
    /// the input lacks; transport failures building the source.
    pub fn into_opened(self, regions: &[Region]) -> Result<StreamedInput, IoError> {
        if regions.is_empty() {
            let (source, facts) = build_source(self.format, self.file)?;
            let meta = InputMeta::new(self.format, &facts.text)?;
            return Ok(StreamedInput {
                source,
                facts,
                meta,
                path: Some(self.path.clone()),
                capability: InputCapability::Replayable(ReplayableInput {
                    path: self.path,
                    _spill: None,
                }),
            });
        }
        // The capability check reads the framing, not the format: the two agree
        // today, but the framing is the property regions actually need (BGZF
        // blocks a seek can land on), and routing the check through it gives
        // `Framing` a production reader instead of a frozen unused API item
        // (doc_v30 §4.7).
        if Framing::of(self.format) != Framing::Bgzf {
            return Err(IoError::Capability {
                needed: reed_io::Needed::Bgzf,
                fix: "recompress plain VCF/regular gzip with bgzip (or drop -r/-R)",
            });
        }
        let index = sidecar_index(&self.path)?;
        let source = RegionChunkSource::open(
            std::io::BufReader::with_capacity(1 << 20, self.file),
            &index,
            regions,
            &self.facts.contig_names,
        )?;
        let meta = InputMeta::new(self.format, &self.facts.text)?;
        Ok(StreamedInput {
            source: Box::new(source),
            facts: self.facts,
            meta,
            path: Some(self.path.clone()),
            capability: InputCapability::Seekable(SeekableInput {
                replay: ReplayableInput {
                    path: self.path,
                    _spill: None,
                },
            }),
        })
    }
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
    Err(IoError::Capability {
        needed: reed_io::Needed::Index,
        fix: "build an index with `reed index` (.csi or .tbi)",
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

/// Reads the header of `reader` under the source the format names, without
/// keeping the source.
fn read_facts(format: Format, reader: impl Read + Send + 'static) -> Result<HeaderFacts, IoError> {
    match format {
        Format::VcfBgzf | Format::Bcf => BgzfChunkSource::open(reader).map(|(_, facts)| facts),
        Format::Vcf => PlainChunkSource::open_vcf(reader).map(|(_, facts)| facts),
        Format::VcfGzip => {
            warn_serial_gzip();
            PlainChunkSource::open_gzip(reader).map(|(_, facts)| facts)
        }
    }
}

/// Opens the whole stream from the first byte under the source the format
/// names: BGZF VCF, plain VCF, plain gzip, BCF (03-architecture.md §6).
fn build_source(
    format: Format,
    reader: impl Read + Send + 'static,
) -> Result<(Box<dyn ChunkSource>, HeaderFacts), IoError> {
    match format {
        Format::VcfBgzf | Format::Bcf => {
            let (source, facts) = BgzfChunkSource::open(reader)?;
            Ok((Box::new(source), facts))
        }
        Format::Vcf => {
            let (source, facts) = PlainChunkSource::open_vcf(reader)?;
            Ok((Box::new(source), facts))
        }
        Format::VcfGzip => {
            warn_serial_gzip();
            let (source, facts) = PlainChunkSource::open_gzip(reader)?;
            Ok((Box::new(source), facts))
        }
    }
}

fn open_stdin() -> Result<OpenedInput, IoError> {
    let mut stdin = std::io::stdin();
    let (prefix, inflated) = sniff(&mut stdin)?;
    let format = Format::detect(&prefix, inflated.as_deref())?;
    // stdin cannot seek, so the probed bytes are chained back on the front.
    let (source, facts) = build_source(format, std::io::Cursor::new(prefix).chain(stdin))?;
    let meta = InputMeta::new(format, &facts.text)?;
    Ok(StreamedInput {
        source,
        facts,
        meta,
        path: None,
        capability: InputCapability::OneShot,
    })
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

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("reed-input-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        dir
    }

    const VCF: &[u8] = b"##fileformat=VCFv4.2\n\
        ##contig=<ID=chr1,length=100>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
        chr1\t1\t.\tA\tG\t50\tPASS\t.\n";

    fn bgzf(bytes: &[u8]) -> Vec<u8> {
        let mut deflater =
            reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6 is valid");
        let mut out = Vec::new();
        deflater
            .deflate_block(bytes, &mut out)
            .expect("invariant: fixture fits one block");
        out.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
        out
    }

    /// doc_v30 R07: the format is a field, not an inference from whether a
    /// dictionary happens to exist.
    #[test]
    fn format_comes_from_the_probe_not_the_dictionary() {
        let dir = temp_dir("format");
        let vcf = dir.join("in.vcf");
        std::fs::write(&vcf, VCF).expect("invariant: fixture");
        let pending = PendingInput::open(&vcf).expect("invariant: probes");
        let opened = pending.into_opened(&[]).expect("invariant: opens");
        assert_eq!(opened.format(), Format::Vcf);
        assert!(!opened.is_bcf());
        assert!(opened.bcf_dict().is_none());
        assert_eq!(opened.meta().framing(), Framing::Plain);
        drop(opened);

        let bcf =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../reed-io/tests/fixtures/bcf/fix.bcf");
        let opened = open(Some(&bcf)).expect("invariant: fixture opens");
        assert_eq!(opened.format(), Format::Bcf);
        assert!(opened.is_bcf());
        assert!(opened.bcf_dict().is_some());
        assert_eq!(opened.meta().framing(), Framing::Bgzf);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// doc_v30 R05: stdin is one-shot; a file can be replayed; a `-r` open
    /// advertises a seek capability.
    #[test]
    fn capabilities_describe_what_the_input_can_do() {
        let dir = temp_dir("capability");
        let vcf = dir.join("in.vcf");
        std::fs::write(&vcf, VCF).expect("invariant: fixture");
        let mut opened = open(Some(&vcf)).expect("invariant: opens");
        assert!(opened.can_replay());
        assert!(matches!(
            opened.take_capability(),
            InputCapability::Replayable(_)
        ));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// doc_v30 R06: regions on a plain VCF fail as a capability before the
    /// region reader sees a "malformed BGZF" block.
    #[test]
    fn regions_on_plain_vcf_are_a_capability_error() {
        let dir = temp_dir("cap-bgzf");
        let vcf = dir.join("in.vcf");
        std::fs::write(&vcf, VCF).expect("invariant: fixture");
        let regions = [Region {
            contig: "chr1".to_string(),
            start1: Some(1),
            end1: Some(10),
        }];
        let err = open_regions(Some(&vcf), &regions)
            .err()
            .expect("plain VCF cannot seek");
        assert!(matches!(
            err,
            IoError::Capability {
                needed: reed_io::Needed::Bgzf,
                ..
            }
        ));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// doc_v30 R06: BGZF without a sidecar names the missing index.
    #[test]
    fn regions_without_an_index_are_a_capability_error() {
        let dir = temp_dir("cap-index");
        let vcf = dir.join("in.vcf.gz");
        std::fs::write(&vcf, bgzf(VCF)).expect("invariant: fixture");
        let regions = [Region {
            contig: "chr1".to_string(),
            start1: Some(1),
            end1: Some(10),
        }];
        let err = open_regions(Some(&vcf), &regions)
            .err()
            .expect("no sidecar");
        assert!(matches!(
            err,
            IoError::Capability {
                needed: reed_io::Needed::Index,
                ..
            }
        ));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
