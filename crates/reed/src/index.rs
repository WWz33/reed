//! `reed index`: build a CSI (default) or TBI index for an existing BGZF
//! VCF, like `bcftools index` (02-cli.md §1).
//!
//! Parallel (03-performance.md O4): the reader cuts BGZF chunks and records
//! each block's `(compressed, inflated)` sizes; workers inflate the chunk and
//! walk its record starts, converting each to a virtual offset; the ordered
//! sink adds the records in chunk order (`Indexer::add` stays serial but
//! cheap). A record's chunk runs from its start to the next record's start —
//! the convention the sequential scan and the filter writer share.

use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use reed_core::ChunkId;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_io::chunk::BgzfChunkSource;
use reed_io::chunk::{inflate_chunk, ChunkBytes};
use reed_io::vcf::fixed::parse_pos1;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::ContigTable;
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

/// One record start (in chunk order).
#[derive(Clone, Copy, Debug)]
struct RecordStart {
    tid: u32,
    pos1: u64,
    ref_len: u32,
    /// Virtual offset of the record's first byte.
    voff: u64,
}

/// A chunk's slice of the record-start sequence.
#[derive(Default)]
struct IndexPartial {
    /// The chunk's first line, when it is a fresh record (the reducer drops
    /// it when the previous chunk left a partial line to complete).
    head: Option<RecordStart>,
    /// Records that start in the body.
    body: Vec<RecordStart>,
    /// The trailing partial line, when it starts a record.
    tail: Option<RecordStart>,
    /// Whether the chunk's text ends mid-line (the next chunk continues it).
    continues: bool,
    /// A seam-record partial: the record's start was already recorded by the
    /// chunk that held its tail, so this carries nothing and changes nothing.
    seam: bool,
}

/// Worker scratch.
struct IndexScratch {
    inflater: reed_io::bgzf::inflate::Inflater,
    text: Vec<u8>,
    contigs: ContigTable,
    starts: Vec<RecordStart>,
}

/// Maps inflated-text offsets back to virtual offsets through the chunk's
/// block table.
struct VoMap<'a> {
    chunk_start: u64,
    blocks: &'a [(u32, u32)],
}

impl VoMap<'_> {
    /// The virtual offset of `text_off` within the chunk's inflated text.
    fn voff(&self, text_off: usize) -> u64 {
        let mut remaining = text_off as u64;
        let mut coff = self.chunk_start;
        for &(clen, ilen) in self.blocks {
            if remaining < u64::from(ilen) {
                return (coff << 16) | remaining;
            }
            remaining -= u64::from(ilen);
            coff += u64::from(clen);
        }
        // Past the last block (empty chunk): clamp to the chunk's end.
        coff << 16
    }
}

/// The chunk job: inflate → walk record starts → virtual offsets.
struct IndexJob {
    contig_names: Arc<[Vec<u8>]>,
}

impl IndexJob {
    fn parse_record(
        &self,
        record: &[u8],
        text_off: usize,
        map: &VoMap<'_>,
        contigs: &mut ContigTable,
    ) -> Result<RecordStart, EngineError> {
        let Some(lr) = LineRef::parse(record) else {
            return Err(EngineError::Job("record line has fewer than 9 columns".into()));
        };
        let Some(tid) = contigs.id_of(lr.chrom()) else {
            return Err(EngineError::Job(format!(
                "{}: contig is not declared in the header",
                String::from_utf8_lossy(lr.chrom())
            )));
        };
        let locus = Locus {
            contig_name: String::from_utf8_lossy(lr.chrom()).into_owned(),
            pos1: 0,
        };
        let pos1 = parse_pos1(lr.pos_text(), &locus).map_err(|e| EngineError::Job(e.to_string()))?;
        let ref_len = u32::try_from(lr.reference().len()).unwrap_or(1).max(1);
        Ok(RecordStart { tid, pos1, ref_len, voff: map.voff(text_off) })
    }
}

impl ChunkJob for IndexJob {
    type Partial = IndexPartial;
    type Scratch = Option<IndexScratch>;

    fn run_chunk(
        &self,
        chunk: &ChunkBytes,
        _body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| IndexScratch {
            inflater: reed_io::bgzf::inflate::Inflater::new(),
            text: Vec::new(),
            contigs: ContigTable::new(self.contig_names.iter().cloned()),
            starts: Vec::new(),
        });
        let view = {
            let IndexScratch { inflater, text, .. } = scratch;
            inflate_chunk(chunk, inflater, text).map_err(job)?
        };
        let map = VoMap { chunk_start: chunk.compressed_start, blocks: &chunk.blocks };
        let contigs = &mut scratch.contigs;
        let mut partial = IndexPartial::default();

        // The first line: a fresh record only when the previous chunk did
        // not leave a tail (the reducer decides; it owns the ordering).
        if let Some(head) = view.head.clone()
            && !view.text[head.clone()].is_empty()
        {
            partial.head = Some(self.parse_record(&view.text[head.clone()], head.start, &map, contigs)?);
        }
        let body = view.body.clone();
        // Body lines, one at a time.
        let mut off = body.start;
        while off < body.end {
            let rest = &view.text[off..body.end];
            let nl = rest.iter().position(|&b| b == b'\n').unwrap_or(rest.len());
            let record = &rest[..nl];
            if !record.is_empty() {
                partial.body.push(self.parse_record(record, off, &map, contigs)?);
            }
            off += nl + 1;
        }
        // The trailing partial line always follows a newline (or is the
        // whole first chunk); when the whole chunk is one partial line the
        // framing reports it as the `head`, never as the tail.
        if let Some(tail) = view.tail.clone()
            && !view.text[tail.clone()].is_empty()
        {
            partial.tail = Some(self.parse_record(&view.text[tail.clone()], tail.start, &map, contigs)?);
        }
        partial.continues = view.text.last() != Some(&b'\n');
        Ok(partial)
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        _line: &[u8],
        _scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        // The seam record's start lives in the previous chunk's tail (or in
        // this chunk's head when the previous chunk ended line-aligned); the
        // ordered sink already has it.
        Ok(IndexPartial { seam: true, ..IndexPartial::default() })
    }
}

/// The ordered sink: pairs consecutive starts into chunks and feeds the
/// indexer; the last record ends at `eof_voff`.
struct IndexSink {
    indexer: Indexer,
    path: PathBuf,
    contig_names: Arc<[Vec<u8>]>,
    pending: Option<RecordStart>,
    prev_continues: bool,
    last_pos1: Option<(u32, u64)>,
    n_records: u64,
    eof_voff: u64,
    out: Arc<Mutex<Option<u64>>>,
}

impl IndexSink {
    fn push(&mut self, entry: RecordStart) -> Result<(), EngineError> {
        if let Some((tid, pos1)) = self.last_pos1
            && tid == entry.tid
            && entry.pos1 < pos1
        {
            let name = String::from_utf8_lossy(&self.contig_names[entry.tid as usize]).into_owned();
            return Err(EngineError::Job(format!(
                "{name}:{}: positions are not ascending (indexing needs sorted input)",
                entry.pos1
            )));
        }
        self.last_pos1 = Some((entry.tid, entry.pos1));
        if let Some(prev) = self.pending.replace(entry) {
            self.add(prev, entry.voff)?;
        }
        self.n_records += 1;
        Ok(())
    }

    fn add(&mut self, entry: RecordStart, end: u64) -> Result<(), EngineError> {
        let name = String::from_utf8_lossy(&self.contig_names[entry.tid as usize]).into_owned();
        self.indexer
            .add(&name, entry.tid, entry.pos1, entry.ref_len, entry.voff, end)
            .map_err(job)
    }
}

impl OrderedSink<IndexPartial> for IndexSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut IndexPartial) -> Result<(), EngineError> {
        if partial.seam {
            return Ok(());
        }
        // A head is a record start unless the previous chunk left a partial
        // line for it to complete.
        if let Some(head) = partial.head.take()
            && !self.prev_continues
        {
            self.push(head)?;
        }
        for entry in &partial.body {
            self.push(*entry)?;
        }
        if let Some(tail) = partial.tail.take() {
            self.push(tail)?;
        }
        self.prev_continues = partial.continues;
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        if let Some(last) = self.pending.take() {
            self.add(last, self.eof_voff)?;
        }
        let IndexSink { indexer, path, n_records, out, .. } = self;
        indexer.write(&path).map_err(job)?;
        *out.lock().expect("invariant: index slot not poisoned") = Some(n_records);
        Ok(n_records)
    }
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
    threads: usize,
) -> Result<IndexOutcome, EngineError> {
    if !flavor.is_some() {
        return Err(EngineError::Job("no index flavor selected".into()));
    }
    let path = match output {
        Some(p) => p.to_path_buf(),
        None => {
            let mut s = input.as_os_str().to_os_string();
            s.push(flavor.suffix());
            PathBuf::from(s)
        }
    };
    let file = std::fs::File::open(input).map_err(|e| EngineError::Io(e.to_string()))?;
    // The compressed file size bounds the last record's chunk.
    let eof_voff = file.metadata().map_err(|e| EngineError::Io(e.to_string()))?.len() << 16;
    let (source, facts) =
        BgzfChunkSource::open(std::io::BufReader::with_capacity(1 << 20, file)).map_err(job)?;
    let contig_names: Arc<[Vec<u8>]> = Arc::from(
        facts
            .contig_names
            .iter()
            .map(|n| n.as_bytes().to_vec())
            .collect::<Vec<_>>()
            .into_boxed_slice(),
    );
    let indexer = Indexer::new(flavor, &contig_names).map_err(job)?;
    let out: Arc<Mutex<Option<u64>>> = Arc::new(Mutex::new(None));
    let sink = IndexSink {
        indexer,
        path: path.clone(),
        contig_names: Arc::clone(&contig_names),
        pending: None,
        prev_continues: false,
        last_pos1: None,
        n_records: 0,
        eof_voff,
        out: Arc::clone(&out),
    };
    let budget = ThreadBudget {
        threads: threads.max(1),
        mem: reed_engine::budget::MemBudget::default_half_memory(),
    };
    let abort = reed_engine::abort::abort();
    let job = IndexJob { contig_names: Arc::clone(&contig_names) };
    let mut source = source;
    run(&mut source, &job, sink, &budget, &abort)?;
    let n_records = out
        .lock()
        .expect("invariant: index slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the count");
    Ok(IndexOutcome { path, n_records })
}
