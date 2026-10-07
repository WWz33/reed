//! `reed index`: build a CSI (default) or TBI index for an existing BGZF
//! VCF, like `bcftools index` (02-cli.md §1).
//!
//! The scan is sequential: each record's start is the virtual position
//! before its bytes, and a record's chunk runs to the next record's start
//! (the same convention the filter writer uses). CSI is built with
//! min_shift 14 / depth 6 so it covers every VCF-legal position; TBI
//! refuses positions ≥ 2^29.

use std::io::BufRead;
use std::path::{Path, PathBuf};

use reed_engine::error::EngineError;
use reed_io::header::HeaderFacts;
use reed_io::vcf::fixed::parse_pos1;
use reed_io::vcf::line::LineRef;
use reed_io::Locus;

/// TBI coordinates are limited to 29 bits (htslib); anything wider needs CSI.
pub const TBI_MAX_POS1: u64 = (1 << 29) - 1;

/// Which index flavor to write. CSI is the default (bcftools parity).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WriteIndex {
    /// No index.
    None,
    /// Coordinate-sorted index (`.csi`), the default.
    #[default]
    Csi,
    /// Tabix index (`.tbi`); positions ≥ 2^29 are refused.
    Tbi,
}

impl WriteIndex {
    /// Parses the optional `--write-index[=csi|tbi]` value.
    ///
    /// # Errors
    ///
    /// On an unknown flavor name.
    pub fn parse(value: Option<&str>) -> Result<Self, EngineError> {
        match value {
            None | Some("") | Some("csi") => Ok(Self::Csi),
            Some("tbi") => Ok(Self::Tbi),
            Some(other) => Err(EngineError::Job(format!(
                "unknown index flavor {other:?} (expected csi or tbi)"
            ))),
        }
    }

    /// Whether an index is requested at all.
    #[must_use]
    pub const fn is_some(self) -> bool {
        !matches!(self, Self::None)
    }

    /// File suffix for this flavor.
    #[must_use]
    pub const fn suffix(self) -> &'static str {
        match self {
            Self::Tbi => ".tbi",
            Self::Csi | Self::None => ".csi",
        }
    }
}

/// The index writer, per flavor. CSI is keyed by header contig id; TBI maps
/// names itself and is limited to 29-bit coordinates.
pub enum Indexer {
    /// Coordinate-sorted index (contig-id keyed, with the contig count).
    Csi(
        Box<noodles::csi::binning_index::indexer::Indexer<noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex>>,
        usize,
    ),
    /// Tabix index (name keyed).
    Tbi(Box<noodles::tabix::index::Indexer>),
}

impl Indexer {
    /// Builds an indexer for `flavor` with `contig_names` in header order.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] when the flavor is `None` or the contig count
    /// overflows CSI's reference-sequence ids.
    pub fn new(flavor: WriteIndex, contig_names: &[Vec<u8>]) -> Result<Self, EngineError> {
        use noodles::csi::binning_index::index::header::ReferenceSequenceNames;
        let names: ReferenceSequenceNames = contig_names
            .iter()
            .map(|v| bstr::BString::from(v.as_slice()))
            .collect();
        let header = noodles::csi::binning_index::index::Header::builder()
            .set_format(noodles::csi::binning_index::index::header::Format::Vcf)
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(None)
            .set_line_comment_prefix(b'#')
            .set_reference_sequence_names(names)
            .build();
        match flavor {
            WriteIndex::Tbi => {
                let mut indexer = noodles::tabix::index::Indexer::default();
                indexer.set_header(header);
                Ok(Self::Tbi(Box::new(indexer)))
            }
            WriteIndex::Csi | WriteIndex::None => {
                use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;
                let indexer =
                    noodles::csi::binning_index::indexer::Indexer::<BinnedIndex>::new(14, 6)
                        .expect("invariant: min_shift/depth are valid")
                        .set_header(header);
                Ok(Self::Csi(Box::new(indexer), contig_names.len()))
            }
        }
    }

    /// Adds one record's chunk `[voff_start, voff_end)`.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] for out-of-range positions or an over-long TBI
    /// coordinate.
    pub fn add(
        &mut self,
        name: &str,
        tid: u32,
        pos1: u64,
        ref_len: u32,
        voff_start: u64,
        voff_end: u64,
    ) -> Result<(), EngineError> {
        use noodles::bgzf::VirtualPosition;
        use noodles::core::Position as CorePosition;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        let chunk = Chunk::new(
            VirtualPosition::from(voff_start),
            VirtualPosition::from(voff_end),
        );
        let start = usize::try_from(pos1)
            .ok()
            .and_then(|p| CorePosition::try_from(p).ok())
            .ok_or_else(|| EngineError::Job(format!("index position {pos1} out of range")))?;
        let end_pos1 = pos1.saturating_add(u64::from(ref_len.saturating_sub(1)));
        let end = usize::try_from(end_pos1)
            .ok()
            .and_then(|p| CorePosition::try_from(p).ok())
            .ok_or_else(|| EngineError::Job(format!("index end {end_pos1} out of range")))?;
        match self {
            Self::Csi(indexer, _) => indexer
                .add_record(
                    Some((
                        usize::try_from(tid).expect("invariant: tid < u32::MAX"),
                        start,
                        end,
                        true,
                    )),
                    chunk,
                )
                .map_err(|e| EngineError::Job(format!("index add_record: {e}"))),
            Self::Tbi(indexer) => {
                if end_pos1 > TBI_MAX_POS1 {
                    return Err(EngineError::Job(format!(
                        "{end_pos1}: TBI coordinates are limited to {TBI_MAX_POS1}; use --write-index=csi"
                    )));
                }
                indexer
                    .add_record(name, start, end, chunk)
                    .map_err(|e| EngineError::Job(format!("index add_record: {e}")))
            }
        }
    }

    /// Writes the index to `path`.
    ///
    /// # Errors
    ///
    /// [`EngineError::Io`] on write failure.
    pub fn write(self, path: &Path) -> Result<(), EngineError> {
        match self {
            Self::Csi(indexer, n_ref) => noodles::csi::fs::write(path, &indexer.build(n_ref))
                .map_err(|e| EngineError::Io(e.to_string())),
            Self::Tbi(indexer) => noodles::tabix::fs::write(path, &indexer.build())
                .map_err(|e| EngineError::Io(e.to_string())),
        }
    }
}

/// Result of `reed index`.
#[derive(Clone, Debug)]
pub struct IndexOutcome {
    /// Index written.
    pub path: PathBuf,
    /// Records indexed.
    pub n_records: u64,
}

/// Builds an index for `input`, writing it next to the input (or at
/// `output`).
///
/// # Errors
///
/// Transport failures, malformed records, unsorted input (positions must not
/// go backwards within a contig), and output failures.
pub fn run_index(
    input: &Path,
    flavor: WriteIndex,
    output: Option<&Path>,
) -> Result<IndexOutcome, EngineError> {
    if !flavor.is_some() {
        return Err(EngineError::Job("no index flavor selected".into()));
    }
    let file = std::fs::File::open(input).map_err(|e| EngineError::Io(e.to_string()))?;
    let mut reader =
        noodles::bgzf::io::Reader::new(std::io::BufReader::with_capacity(1 << 20, file));

    // Header: lines up to and including #CHROM.
    let mut header_text = Vec::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        let n = reader
            .read_until(b'\n', &mut line)
            .map_err(|e| EngineError::Io(e.to_string()))?;
        if n == 0 {
            return Err(EngineError::Job("input ends before the #CHROM line".into()));
        }
        header_text.extend_from_slice(&line);
        if line.starts_with(b"#CHROM") {
            break;
        }
    }
    let facts = HeaderFacts::parse(&header_text).map_err(EngineError::from)?;
    let contig_names: Vec<Vec<u8>> = facts
        .contig_names
        .iter()
        .map(|n| n.as_bytes().to_vec())
        .collect();
    let mut indexer = Indexer::new(flavor, &contig_names)?;

    let mut pending: Option<(u32, u64, u32, u64)> = None;
    let mut last_pos1: Option<(u32, u64)> = None;
    let mut n_records = 0_u64;
    loop {
        let start: u64 = reader.virtual_position().into();
        line.clear();
        let n = reader
            .read_until(b'\n', &mut line)
            .map_err(|e| EngineError::Io(e.to_string()))?;
        if n == 0 {
            break;
        }
        let record = line.strip_suffix(b"\n").unwrap_or(&line);
        if record.is_empty() {
            continue;
        }
        let Some(lr) = LineRef::parse(record) else {
            return Err(EngineError::Job("record line has fewer than 9 columns".into()));
        };
        let chrom = lr.chrom().to_vec();
        let locus = Locus { contig_name: String::from_utf8_lossy(&chrom).into_owned(), pos1: 0 };
        let pos1 = parse_pos1(lr.pos_text(), &locus).map_err(|e| EngineError::Job(e.to_string()))?;
        let Some(tid) = contig_names
            .iter()
            .position(|name| name == &chrom)
            .and_then(|i| u32::try_from(i).ok())
        else {
            return Err(EngineError::Job(format!(
                "{}: contig is not declared in the header",
                String::from_utf8_lossy(&chrom)
            )));
        };
        if let Some((prev_tid, prev_pos1)) = last_pos1
            && prev_tid == tid
            && pos1 < prev_pos1
        {
            return Err(EngineError::Job(format!(
                "{locus}: positions are not ascending (indexing needs sorted input)"
            )));
        }
        if let Some((p_tid, p_pos1, p_ref_len, p_start)) = pending.take() {
            let name = String::from_utf8_lossy(&contig_names[p_tid as usize]).into_owned();
            indexer.add(&name, p_tid, p_pos1, p_ref_len, p_start, start)?;
        }
        let ref_len = u32::try_from(lr.reference().len()).unwrap_or(1).max(1);
        pending = Some((tid, pos1, ref_len, start));
        last_pos1 = Some((tid, pos1));
        n_records += 1;
    }
    if let Some((tid, pos1, ref_len, voff)) = pending {
        let end: u64 = reader.virtual_position().into();
        let name = String::from_utf8_lossy(&contig_names[tid as usize]).into_owned();
        indexer.add(&name, tid, pos1, ref_len, voff, end)?;
    }
    let path = match output {
        Some(p) => p.to_path_buf(),
        None => {
            let mut s = input.as_os_str().to_os_string();
            s.push(flavor.suffix());
            PathBuf::from(s)
        }
    };
    indexer.write(&path)?;
    Ok(IndexOutcome { path, n_records })
}
