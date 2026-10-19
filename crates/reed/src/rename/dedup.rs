//! `--rm-dup-ids` two-phase output (doc_v15 D15-28).
//!
//! The deduplicated text is ordered by the sink — the suffix must follow the
//! global record order, which only the ordered stage knows — but deflating
//! it there too made the whole path serial: on the 3M-record repro the
//! sink's deflater alone spent 0.6 s of 1.4 s. The two phases split that:
//!
//! 1. **serial**: the sink assigns suffixes and fills byte-driven batches
//!    (fixed target, whole records only);
//! 2. **parallel**: sealed batches are deflated by a small pool and handed
//!    back in batch order.
//!
//! Batch seals are driven by accumulated bytes, never by timing or thread
//! count, so the block layout — and with it the whole compressed stream —
//! is identical for every `--threads` value.

// Casts in this file are bounded by construction: elapsed nanoseconds are
// clamped to u64 in practice (584 years).
#![allow(clippy::cast_possible_truncation)]

use std::collections::VecDeque;
use std::sync::mpsc;

use reed_engine::error::EngineError;
use reed_io::IoError;
use reed_io::bgzf::deflate::Deflater;
use reed_io::write::block_writer::{BGZF_BLOCK_SIZE, ChunkOutput, IndexKey, IndexRecord, OutBlock};
use reed_io::write::ordered_output::OrderedOutput;

/// Uncompressed bytes per sealed batch. Whole records only, so a batch is
/// this plus at most one record.
const BATCH_BYTES: usize = 1 << 20;

/// In-flight batches beyond the pool size: the seal blocks on the oldest
/// result once the window is full, which bounds the batch memory.
const INFLIGHT_EXTRA: usize = 2;

/// The memory budget's share per in-flight batch: the window never exceeds
/// `mem / BATCH_COST_DIVISOR`, so `--memory` also bounds the pool's batches
/// (they are output buffers the engine's own accounting does not see).
const BATCH_COST_DIVISOR: u64 = 4;

/// One sealed batch on its way through the pool.
struct Job {
    bytes: Vec<u8>,
    anchors: Vec<(u32, IndexKey)>,
    done: mpsc::Sender<Result<ChunkOutput, IoError>>,
}

/// The deflate pool: `threads` workers behind one shared job queue. Each
/// job carries its own result channel, so the sink's in-order drain is
/// independent of the order the workers happen to finish in.
struct Pool {
    /// `None` once dropped: disconnecting the queue is what ends the
    /// workers' loop.
    jobs: Option<mpsc::Sender<Job>>,
    handles: Vec<std::thread::JoinHandle<()>>,
}

impl Pool {
    fn new(threads: usize, level: u8) -> Result<Self, EngineError> {
        let (tx, rx) = mpsc::channel::<Job>();
        let rx = std::sync::Arc::new(std::sync::Mutex::new(rx));
        let mut handles = Vec::with_capacity(threads);
        for _ in 0..threads {
            let rx = std::sync::Arc::clone(&rx);
            let handle = std::thread::Builder::new()
                .name("reed-rename-deflate".into())
                .spawn(move || {
                    let mut deflater = Deflater::new(i32::from(level))
                        .expect("invariant: the level comes from the CLI clamp");
                    loop {
                        // Holding the queue lock across `recv` is the point:
                        // there is no other work to hand out meanwhile, and
                        // the lock is released before the batch is deflated.
                        let Ok(job) = rx
                            .lock()
                            .expect("invariant: the job queue lock is not poisoned")
                            .recv()
                        else {
                            break; // the sink finished or aborted
                        };
                        let result = deflate_batch(&job.bytes, &job.anchors, &mut deflater);
                        // A vanished receiver means the sink aborted; the
                        // worker just moves on to the disconnect check.
                        let _ = job.done.send(result);
                    }
                })
                .map_err(EngineError::from)?;
            handles.push(handle);
        }
        Ok(Self {
            jobs: Some(tx),
            handles,
        })
    }

    fn submit(&self, job: Job) -> Result<(), EngineError> {
        self.jobs
            .as_ref()
            .expect("invariant: the pool is live until its Drop")
            .send(job)
            .map_err(|_| EngineError::job("rename: the deflate pool is gone"))
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        // Disconnect the queue first: idle workers see it and exit, a worker
        // mid-batch finishes the batch and exits on its next round.
        self.jobs.take();
        for handle in self.handles.drain(..) {
            let _ = handle.join();
        }
    }
}

/// The sink-side state: the serial batch builder plus the pool bookkeeping.
pub(super) struct BatchWriter {
    level: u8,
    compressed: bool,
    /// An index is being written: per-record anchors are kept (and the
    /// caller passes real keys; without an index the keys are never read,
    /// exactly as on the non-dedup path).
    need_key: bool,
    /// Pool size.
    threads: usize,
    /// Seals block on the oldest result once this many batches are in
    /// flight: the pool size plus a couple, capped by the memory budget.
    window: usize,
    batch: Vec<u8>,
    anchors: Vec<(u32, IndexKey)>,
    pool: Option<Pool>,
    inflight: VecDeque<mpsc::Receiver<Result<ChunkOutput, IoError>>>,
}

impl BatchWriter {
    pub(super) fn new(
        level: u8,
        compressed: bool,
        need_key: bool,
        threads: usize,
        mem_bytes: u64,
    ) -> Self {
        let threads = threads.max(1);
        Self {
            level,
            compressed,
            need_key,
            threads,
            window: inflight_window(threads, mem_bytes),
            batch: Vec::with_capacity(BATCH_BYTES + (1 << 16)),
            anchors: Vec::new(),
            pool: None,
            inflight: VecDeque::new(),
        }
    }

    /// Appends one record and seals the batch once it reached the target.
    /// The caller passes the bare line; the terminator is added here (a
    /// line that already carries one is not doubled).
    pub(super) fn push_record(
        &mut self,
        out: &mut OrderedOutput,
        bytes: &[u8],
        key: IndexKey,
    ) -> Result<(), EngineError> {
        if self.need_key {
            let offset = u32::try_from(self.batch.len())
                .map_err(|_| EngineError::job("rename: batch exceeds u32::MAX bytes"))?;
            self.anchors.push((offset, key));
        }
        self.batch.extend_from_slice(bytes);
        if !bytes.ends_with(b"\n") {
            self.batch.push(b'\n');
        }
        if self.batch.len() >= BATCH_BYTES {
            self.seal(out)?;
        }
        Ok(())
    }

    /// Seals the tail batch and hands every outstanding result to `out`.
    pub(super) fn finish(&mut self, out: &mut OrderedOutput) -> Result<(), EngineError> {
        if !self.batch.is_empty() {
            if self.compressed && self.pool.is_none() {
                // Output that never filled a batch: deflate the tail inline
                // instead of spawning a pool for one small job.
                let bytes = std::mem::take(&mut self.batch);
                let anchors = std::mem::take(&mut self.anchors);
                let mut deflater =
                    Deflater::new(i32::from(self.level)).map_err(EngineError::from)?;
                let chunk =
                    deflate_batch(&bytes, &anchors, &mut deflater).map_err(EngineError::from)?;
                out.accept_chunk(&chunk).map_err(EngineError::from)?;
            } else {
                self.seal(out)?;
            }
        }
        while !self.inflight.is_empty() {
            self.drain_front(out)?;
        }
        Ok(())
    }

    /// Seals the pending batch: the plain path writes the bytes straight
    /// through, the compressed path hands a job to the pool.
    fn seal(&mut self, out: &mut OrderedOutput) -> Result<(), EngineError> {
        if self.batch.is_empty() {
            return Ok(());
        }
        let bytes = std::mem::take(&mut self.batch);
        let anchors = std::mem::take(&mut self.anchors);
        if !self.compressed {
            let chunk = ChunkOutput {
                raw: bytes,
                ..ChunkOutput::default()
            };
            return out.accept_chunk(&chunk).map_err(EngineError::from);
        }
        if self.pool.is_none() {
            self.pool = Some(Pool::new(self.threads, self.level)?);
        }
        let (tx, rx) = mpsc::channel();
        self.pool
            .as_ref()
            .expect("invariant: inserted above")
            .submit(Job {
                bytes,
                anchors,
                done: tx,
            })?;
        self.inflight.push_back(rx);
        while self.inflight.len() > self.window {
            self.drain_front(out)?;
        }
        Ok(())
    }

    /// Blocks on the oldest outstanding batch and appends its blocks.
    fn drain_front(&mut self, out: &mut OrderedOutput) -> Result<(), EngineError> {
        let rx = self
            .inflight
            .pop_front()
            .expect("invariant: the caller checked for an in-flight batch");
        let chunk = rx
            .recv()
            .map_err(|_| EngineError::job("rename: a deflate pool worker panicked"))?
            .map_err(EngineError::from)?;
        out.accept_chunk(&chunk).map_err(EngineError::from)
    }
}

/// The in-flight window: the pool size plus a couple of batches, capped by
/// the run's memory budget (a quarter of it, at least two batches — the
/// floor keeps a seal able to make progress under any budget).
fn inflight_window(threads: usize, mem_bytes: u64) -> usize {
    let cap = usize::try_from(mem_bytes / (BATCH_COST_DIVISOR * BATCH_BYTES as u64))
        .unwrap_or(usize::MAX)
        .max(2);
    (threads + INFLIGHT_EXTRA).min(cap)
}

/// Compresses one batch into BGZF blocks. Cut points are the fixed
/// [`BGZF_BLOCK_SIZE`] boundaries — the same layout
/// [`reed_io::write::block_writer::BlockWriter`] produces — and the record
/// anchors are re-based into the block holding each record's start.
fn deflate_batch(
    bytes: &[u8],
    anchors: &[(u32, IndexKey)],
    deflater: &mut Deflater,
) -> Result<ChunkOutput, IoError> {
    let mut out = ChunkOutput::default();
    let mut cursor = 0_usize;
    let mut base = 0_u32;
    for chunk in bytes.chunks(BGZF_BLOCK_SIZE) {
        let len = u32::try_from(chunk.len())
            .map_err(|_| IoError::Io("block exceeds u32::MAX bytes".into()))?;
        let mut compressed = Vec::with_capacity(chunk.len() / 2 + 64);
        let started = std::time::Instant::now();
        deflater.deflate_block(chunk, &mut compressed)?;
        out.deflate_nanos += started.elapsed().as_nanos() as u64;
        let hi = base + len;
        let end = anchors.partition_point(|(offset, _)| *offset < hi);
        let records = anchors[cursor..end]
            .iter()
            .map(|(offset, key)| IndexRecord {
                key: *key,
                offset_in_block: offset - base,
            })
            .collect();
        cursor = end;
        base = hi;
        out.blocks.push(OutBlock {
            compressed,
            uncompressed_len: len,
            records,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_io::bgzf::inflate::Inflater;

    /// A record's key, distinguishable per record.
    fn key(start1: u64) -> IndexKey {
        IndexKey {
            contig: reed_core::ContigId(3),
            start1,
            ref_len: u16::try_from(start1 % 7 + 1).map_or(1, u32::from),
        }
    }

    /// doc_v15 D15-28: the in-flight window follows the pool size, and the
    /// memory budget caps it once that is the smaller bound.
    #[test]
    fn the_inflight_window_follows_the_budget() {
        assert_eq!(inflight_window(8, u64::MAX), 10);
        // A quarter of a 64 MiB budget holds 16 batches of 1 MiB.
        assert_eq!(inflight_window(54, 64 << 20), 16);
        // Below two batches the floor wins: progress beats the budget.
        assert_eq!(inflight_window(54, 1 << 20), 2);
        assert_eq!(inflight_window(1, 0), 2);
    }

    /// doc_v15 D15-28: dropping the writer with a job still in flight (the
    /// engine's abort path drops the sink without `finish`) must join the
    /// pool instead of hanging on it.
    #[test]
    fn dropping_with_a_job_in_flight_joins_the_pool() {
        let mut writer = BatchWriter::new(6, true, false, 2, u64::MAX);
        let dummy = std::io::sink();
        let mut out = OrderedOutput::new(
            Box::new(dummy),
            true,
            6,
            None,
            None,
            std::sync::Arc::from(Vec::new()),
        );
        let key = key(1);
        // Two seals: the second waits behind a scheduled batch, so at least
        // one result is outstanding when the writer goes away.
        for i in 0..40_000_u64 {
            let line = format!("chr1\t{i}\trs{i}\tA\tG\t50\t.\t.\tGT\t0/1\n");
            writer
                .push_record(&mut out, line.as_bytes(), key)
                .expect("invariant: push");
        }
        drop(writer);
    }

    /// doc_v15 D15-28: a batch cuts at the fixed member boundaries and every
    /// anchor lands in the block holding the record's *start* — including a
    /// record that starts before a cut and runs past it.
    #[test]
    fn batch_blocks_cut_at_the_member_limit_and_anchor_records() {
        // Records of varying width, with the 4th one straddling the first
        // cut: 65,280 is not a record boundary here.
        let mut bytes = Vec::new();
        let mut anchors = Vec::new();
        let mut start1 = 1_u64;
        while bytes.len() < 200_000 {
            let offset = u32::try_from(bytes.len()).expect("invariant: fits u32");
            anchors.push((offset, key(start1)));
            let width = if start1 == 4 {
                BGZF_BLOCK_SIZE + 37
            } else {
                40 + (start1 as usize * 7) % 90
            };
            bytes.extend(std::iter::repeat_n(b'a', width));
            bytes.push(b'\n');
            start1 += 1;
        }
        let mut deflater = Deflater::new(6).expect("invariant: level 6");
        let out = deflate_batch(&bytes, &anchors, &mut deflater).expect("invariant: deflate");

        // Cuts are at the member limit: every full block is exactly
        // BGZF_BLOCK_SIZE bytes, the tail is the remainder.
        let mut burned = 0_usize;
        for block in &out.blocks {
            let expect_len = (bytes.len() - burned).min(BGZF_BLOCK_SIZE);
            assert_eq!(block.uncompressed_len as usize, expect_len);
            burned += expect_len;
        }
        assert_eq!(burned, bytes.len());

        // Every record anchors in the block that holds its start offset.
        let mut seen = 0_usize;
        for block in &out.blocks {
            for record in &block.records {
                let (offset, want_key) = anchors[seen];
                let base = u32::try_from(
                    out.blocks
                        .iter()
                        .take_while(|b| !std::ptr::eq(*b, block))
                        .map(|b| b.uncompressed_len as usize)
                        .sum::<usize>(),
                )
                .expect("invariant: fits u32");
                assert_eq!(record.offset_in_block, offset - base);
                assert_eq!(record.key, want_key);
                assert!((base..base + block.uncompressed_len).contains(&offset));
                seen += 1;
            }
        }
        assert_eq!(seen, anchors.len());

        // The decompressed stream is exactly the batch that went in.
        let mut inflater = Inflater::new();
        let mut text = Vec::new();
        for block in &out.blocks {
            let mut scratch = Vec::new();
            inflater
                .inflate_block(&block.compressed, &mut scratch)
                .expect("invariant: block inflates");
            assert_eq!(scratch.len(), block.uncompressed_len as usize);
            text.extend_from_slice(&scratch);
        }
        assert_eq!(text, bytes);
        assert!(out.deflate_nanos > 0);
    }
}
