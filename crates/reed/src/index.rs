//! `reed index`: build a CSI (default) or TBI index for an existing BGZF
//! VCF, like `bcftools index` (02-cli.md §1).
//!
//! The scan is sequential: each record's start is the virtual position
//! before its bytes, and a record's chunk runs to the next record's start
//! (the same convention the filter writer uses). The indexer itself lives in
//! `reed_io::write` — one implementation for both this command and the
//! filter's `--write-index`.

use std::io::BufRead;
use std::path::{Path, PathBuf};

use reed_engine::error::EngineError;
use reed_io::header::HeaderFacts;
use reed_io::vcf::fixed::parse_pos1;
use reed_io::vcf::line::LineRef;
use reed_io::write::index::Indexer;
use reed_io::Locus;

pub use reed_io::write::index::WriteIndex;

/// A transport-layer failure inside the index driver.
fn job(err: reed_io::IoError) -> EngineError {
    EngineError::Job(err.to_string())
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
    let mut indexer = Indexer::new(flavor, &contig_names).map_err(job)?;

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
            indexer.add(&name, p_tid, p_pos1, p_ref_len, p_start, start).map_err(job)?;
        }
        let ref_len = u32::try_from(lr.reference().len()).unwrap_or(1).max(1);
        pending = Some((tid, pos1, ref_len, start));
        last_pos1 = Some((tid, pos1));
        n_records += 1;
    }
    if let Some((tid, pos1, ref_len, voff)) = pending {
        let end: u64 = reader.virtual_position().into();
        let name = String::from_utf8_lossy(&contig_names[tid as usize]).into_owned();
        indexer
            .add(&name, tid, pos1, ref_len, voff, end)
            .map_err(job)?;
    }
    let path = match output {
        Some(p) => p.to_path_buf(),
        None => {
            let mut s = input.as_os_str().to_os_string();
            s.push(flavor.suffix());
            PathBuf::from(s)
        }
    };
    indexer.write(&path).map_err(job)?;
    Ok(IndexOutcome { path, n_records })
}
