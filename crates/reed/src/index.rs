//! `reed index`: build a CSI (default) or TBI index for an existing BGZF
//! VCF, like `bcftools index` (02-cli.md §1).
//!
//! Parallel (03-performance.md O4): the reader cuts BGZF chunks and records
//! each block's `(compressed, inflated)` sizes; workers inflate the chunk and
//! walk its record starts, converting each to a virtual offset; the ordered
//! sink adds the records in chunk order (`Indexer::add` stays serial but
//! cheap). A record's chunk runs from its start to the next record's start —
//! the convention the sequential scan and the filter writer share.

// Casts in this file are bounded by construction: Carried record prefixes are bounded by the 1 GiB record ceiling.
#![allow(clippy::cast_possible_truncation)]

use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use reed_core::{ChunkId, ContigId};
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::frame_chain::{FrameChain, SharedChain};
use reed_engine::pipeline::{ChunkJob, OrderedSink, run};
use reed_io::bcf::{VirtualMap, record_len, record_site};
use reed_io::chunk::framing::RecordFraming;
use reed_io::chunk::{ChunkBytes, inflate_chunk};
use reed_io::vcf::fixed::parse_pos1;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::ContigTable;
use reed_io::write::index::Indexer;

pub use reed_io::write::index::WriteIndex;

/// A transport-layer failure inside the index driver.
fn job(err: reed_io::IoError) -> EngineError {
    EngineError::from(err)
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
pub(super) struct RecordStart {
    tid: ContigId,
    pos1: u64,
    ref_len: u32,
    /// Virtual offset of the record's first byte.
    voff: u64,
}

/// A chunk's slice of the record-start sequence.
#[derive(Default)]
pub(super) struct IndexPartial {
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
            return Err(EngineError::job("record line has fewer than 9 columns"));
        };
        let Some(tid) = contigs.id_of(lr.chrom()) else {
            return Err(EngineError::job(format!(
                "{}: contig is not declared in the header",
                String::from_utf8_lossy(lr.chrom())
            )));
        };
        let pos1 = parse_pos1(lr.pos_text(), lr.chrom()).map_err(EngineError::from)?;
        let ref_len = u32::try_from(lr.reference().len()).unwrap_or(1).max(1);
        Ok(RecordStart {
            tid,
            pos1,
            ref_len,
            voff: map.voff(text_off),
        })
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
        });
        let view = {
            let IndexScratch { inflater, text, .. } = scratch;
            inflate_chunk(chunk, inflater, text).map_err(job)?
        };
        let map = VoMap {
            chunk_start: chunk.compressed_start,
            blocks: &chunk.blocks,
        };
        let contigs = &mut scratch.contigs;
        let mut partial = IndexPartial::default();

        // The first line: a fresh record only when the previous chunk did
        // not leave a tail (the reducer decides; it owns the ordering).
        if let Some(head) = view.head.clone()
            && !view.text[head.clone()].is_empty()
        {
            partial.head =
                Some(self.parse_record(&view.text[head.clone()], head.start, &map, contigs)?);
        }
        let body = view.body.clone();
        // Body lines, one at a time.
        let mut off = body.start;
        while off < body.end {
            let rest = &view.text[off..body.end];
            let nl = rest.iter().position(|&b| b == b'\n').unwrap_or(rest.len());
            let record = &rest[..nl];
            if !record.is_empty() {
                partial
                    .body
                    .push(self.parse_record(record, off, &map, contigs)?);
            }
            off += nl + 1;
        }
        // The trailing partial line always follows a newline (or is the
        // whole first chunk); when the whole chunk is one partial line the
        // framing reports it as the `head`, never as the tail.
        if let Some(tail) = view.tail.clone()
            && !view.text[tail.clone()].is_empty()
        {
            partial.tail =
                Some(self.parse_record(&view.text[tail.clone()], tail.start, &map, contigs)?);
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
        Ok(IndexPartial {
            seam: true,
            ..IndexPartial::default()
        })
    }
}

/// The BCF index job: same entries as the text job, but the framing walks
/// length prefixes and the chain carries the crossing record's start.
pub(super) struct BcfIndexJob {
    chain: SharedChain,
}

impl ChunkJob for BcfIndexJob {
    type Partial = IndexPartial;
    type Scratch = ();

    fn framing(&self) -> RecordFraming {
        RecordFraming::Raw
    }

    fn run_chunk(
        &self,
        chunk: &ChunkBytes,
        text: &[u8],
        _scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let chunk_id = chunk.id.get();
        let skip = if chunk_id == 0 {
            usize::try_from(chunk.first_in_block_offset).unwrap_or(0)
        } else {
            0
        };
        let map = VirtualMap {
            chunk_start: chunk.compressed_start,
            blocks: &chunk.blocks,
        };
        let mut carried = self.chain.wait(chunk_id)?;
        let mut partial = IndexPartial::default();
        let mut guard = self.chain.publish_guard(chunk_id + 1);
        if skip >= text.len() {
            guard.publish(reed_engine::frame_chain::Carried::default());
            return Ok(partial);
        }
        let body_start;
        let body_end;
        if carried.is_empty() {
            body_start = skip;
            body_end = walk_bcf_records(text, skip);
        } else {
            let total = declared_total(&carried.head, text)
                .ok_or_else(|| EngineError::job("BCF record length prefix is truncated"))?;
            let need = total.saturating_sub(carried.prefix_len as usize);
            if need > text.len() {
                // The record still spans this chunk: keep the header bytes.
                let mut head = carried.head;
                let extra = 24_usize.saturating_sub(head.len()).min(text.len());
                head.extend_from_slice(&text[..extra]);
                carried.prefix_len += text.len() as u64;
                guard.publish(reed_engine::frame_chain::Carried {
                    head,
                    prefix_len: carried.prefix_len,
                    voff: carried.voff,
                });
                return Ok(partial);
            }
            body_start = need;
            body_end = walk_bcf_records(text, need);
            partial.head = Some(entry_at(&carried.head, text, body_start, carried.voff)?);
        }
        // Publish the crossing record of the next chunk (its first 24 bytes
        // carry everything the index decodes).
        let mut head = text[body_end..].to_vec();
        head.truncate(24);
        guard.publish(reed_engine::frame_chain::Carried {
            head,
            prefix_len: (text.len() - body_end) as u64,
            voff: map.voff(body_end),
        });
        let mut q = body_start;
        while q < body_end {
            let len = record_len(&text[q..]).expect("invariant: walked records are complete");
            // The shared block starts after the 8-byte length prefix.
            let (tid, pos0, rlen) = record_site(&text[q + 8..]).map_err(EngineError::from)?;
            partial.body.push(RecordStart {
                tid: ContigId(
                    u32::try_from(tid)
                        .map_err(|_| EngineError::job("negative BCF contig index"))?,
                ),
                pos1: u64::try_from(pos0).unwrap_or(0) + 1,
                ref_len: u32::try_from(rlen).unwrap_or(1).max(1),
                voff: map.voff(q),
            });
            q += len;
        }
        Ok(partial)
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        _line: &[u8],
        _scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        Err(EngineError::job(
            "internal: the BCF index never receives seam lines",
        ))
    }
}

/// Walks length prefixes; the returned offset is where the crossing record
/// starts (== `text.len()` when the chunk ends record-aligned).
fn walk_bcf_records(text: &[u8], from: usize) -> usize {
    let mut p = from;
    while let Some(len) = record_len(&text[p..]) {
        p += len;
    }
    p
}

/// The total length declared by the record's first eight bytes, taken from
/// `head` and (when `head` is short) the start of `text`.
fn declared_total(head: &[u8], text: &[u8]) -> Option<usize> {
    if head.len() >= 8 {
        return reed_io::bcf::declared_len(head);
    }
    let need = 8 - head.len();
    let mut buf = [0_u8; 8];
    buf[..head.len()].copy_from_slice(head);
    buf[head.len()..].copy_from_slice(&text[..need]);
    reed_io::bcf::declared_len(&buf)
}

/// Decodes a crossing record's entry: `head` holds the record's leading
/// bytes (length prefix included), completed by the start of `text` when
/// the previous chunk held fewer than 20 bytes.
fn entry_at(
    head: &[u8],
    text: &[u8],
    body_start: usize,
    voff: u64,
) -> Result<RecordStart, EngineError> {
    // chrom/pos/rlen sit at record offsets 8..20.
    let mut buf = [0_u8; 20];
    let held = head.len().min(20);
    buf[..held].copy_from_slice(&head[..held]);
    if held < 20 {
        // The record's own bytes continue after the carried prefix; the
        // missing bytes are the *first* ones of this chunk.
        let missing = 20 - held;
        if text.len() < missing {
            return Err(EngineError::job("BCF record header spans multiple chunks"));
        }
        buf[held..].copy_from_slice(&text[..missing]);
    }
    let _ = body_start;
    let (tid, pos0, rlen) = record_site(&buf[8..]).map_err(EngineError::from)?;
    Ok(RecordStart {
        tid: ContigId(
            u32::try_from(tid).map_err(|_| EngineError::job("negative BCF contig index"))?,
        ),
        pos1: u64::try_from(pos0).unwrap_or(0) + 1,
        ref_len: u32::try_from(rlen).unwrap_or(1).max(1),
        voff,
    })
}

/// The ordered sink: pairs consecutive starts into chunks and feeds the
/// indexer; the last record ends at `eof_voff`.
struct IndexSink {
    indexer: Indexer,
    path: PathBuf,
    contig_names: Arc<[Vec<u8>]>,
    pending: Option<RecordStart>,
    prev_continues: bool,
    last_pos1: Option<(ContigId, u64)>,
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
            let name =
                String::from_utf8_lossy(&self.contig_names[entry.tid.get() as usize]).into_owned();
            return Err(EngineError::job(format!(
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
        // The CSI path never touches the name: no per-record String (A6-12).
        let name = &self.contig_names[entry.tid.get() as usize];
        self.indexer
            .add(name, entry.tid, entry.pos1, entry.ref_len, entry.voff, end)
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
        let IndexSink {
            indexer,
            path,
            n_records,
            out,
            ..
        } = self;
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
        return Err(EngineError::job("no index flavor selected"));
    }
    let path = match output {
        Some(p) => p.to_path_buf(),
        None => {
            let mut s = input.as_os_str().to_os_string();
            s.push(flavor.suffix());
            PathBuf::from(s)
        }
    };
    // The compressed file size bounds the last record's chunk.
    let eof_voff = std::fs::metadata(input).map_err(EngineError::Io)?.len() << 16;
    let mut opened = crate::input::open(Some(input)).map_err(job)?;
    let is_bcf = opened.bcf.is_some();
    let contig_names: Arc<[Vec<u8>]> = Arc::from(
        opened
            .facts
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
    let report = if is_bcf {
        let chain = Arc::new(FrameChain::new(abort.signal.clone()));
        chain.publish(0, reed_engine::frame_chain::Carried::default());
        let job = BcfIndexJob {
            chain: Arc::clone(&chain),
        };
        let report = run(opened.source.as_mut(), &job, sink, &budget, &abort)?;
        chain.finish();
        if chain
            .peek(report.n_chunks)
            .is_some_and(|carried| !carried.is_empty())
        {
            return Err(EngineError::Job(
                "BCF input ends inside a record (truncated file)".into(),
            ));
        }
        report
    } else {
        let job = IndexJob {
            contig_names: Arc::clone(&contig_names),
        };
        run(opened.source.as_mut(), &job, sink, &budget, &abort)?
    };
    let _ = report;
    let n_records = out
        .lock()
        .expect("invariant: index slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the count");
    Ok(IndexOutcome { path, n_records })
}
