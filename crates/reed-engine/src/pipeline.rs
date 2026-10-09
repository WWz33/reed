//! The chunk pipeline skeleton (`docs/architecture-spec.md` §7.1–§7.6):
//! reader → bounded chunk channel (ticket window) → worker pool (fused
//! inflate + job) → ordered reducer with seam chain-stitching.
//!
//! Threads never exchange per-record heap objects — only compressed chunks
//! and small partials (§7.2). Every blocking point selects on the abort
//! signal (§7.5), and in-flight compressed bytes are accounted against the
//! memory budget (§7.7).

// Casts in this file are bounded by construction: Elapsed nanoseconds are clamped to u64 in practice (584 years).
#![allow(clippy::cast_possible_truncation)]

use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

use crossbeam_channel::bounded;
use reed_core::ChunkId;
use reed_io::chunk::framing::RecordFraming;
use reed_io::chunk::{inflate_chunk, ChunkBytes, ChunkSource};

pub use reed_io::chunk::ChunkBytes as ChunkBytesExport;

use crate::abort::{Abort, AbortSignal, FirstError};
use crate::budget::{BytePermits, ThreadBudget};
use crate::error::EngineError;
use crate::metrics::StageMetrics;
use crate::reorder::ReorderBuffer;

/// One chunk's work output: the reduced partial plus the seam fragments the
/// reducer needs (`head` belongs to the previous chunk's last record).
pub struct ChunkOutput<P> {
    /// Input-order chunk id.
    pub id: ChunkId,
    /// Partial result of this chunk's complete records.
    pub partial: P,
    /// Partial line at the chunk start (`None` for the first chunk).
    pub head: Option<Vec<u8>>,
    /// Partial line after the last newline (`None` when line-aligned).
    pub tail: Option<Vec<u8>>,
}

/// A fused chunk computation (`docs/architecture-spec.md` §7.1).
///
/// `run_chunk` receives the chunk's complete-record bytes (head fragment
/// excluded); `run_seam` receives one chain-stitched seam record. Partials
/// merge at the sink in chunk order.
pub trait ChunkJob: Sync {
    /// Additive partial result (`Default` is the identity element).
    type Partial: Send + Default;
    /// Worker-local reusable buffers.
    type Scratch: Send + Default;

    /// How the chunk's inflated text is delimited into records. With
    /// [`RecordFraming::Raw`] the engine hands the job the whole text and
    /// leaves seam stitching to it (the job owns the frame chain).
    fn framing(&self) -> RecordFraming {
        RecordFraming::Newline
    }

    /// Processes the chunk's complete records (newline-terminated lines).
    ///
    /// # Errors
    ///
    /// Whatever the job considers a run-stopping failure.
    fn run_chunk(
        &self,
        chunk: &ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError>;

    /// Cumulative time workers spent waiting on job-internal
    /// synchronization (the BCF frame chain). Summed into
    /// [`StageMetrics::chain_wait`]; zero for jobs without a chain.
    fn wait_time(&self) -> Duration {
        Duration::ZERO
    }

    /// Processes one seam record (a complete line ending in `\n`).
    ///
    /// `id` is the chunk whose head fragment completed the seam — the seam
    /// record itself belongs to the *previous* chunk's stream (jobs that
    /// key state by chunk use `id - 1`).
    ///
    /// # Errors
    ///
    /// Whatever the job considers a run-stopping failure.
    fn run_seam(
        &self,
        id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError>;
}

/// Ordered consumer of partials: `accept` merges in chunk order; seam
/// records arrive before the chunk they precede.
pub trait OrderedSink<P>: Send {
    /// Accepts one partial in order.
    ///
    /// # Errors
    ///
    /// Output failures (the engine aborts on them).
    fn accept(&mut self, id: ChunkId, partial: &mut P) -> Result<(), EngineError>;

    /// Finalizes the sink after the last partial; returns its own report
    /// value (e.g. the total number of accepted records).
    ///
    /// # Errors
    ///
    /// Whatever flushing can produce.
    fn finish(self) -> Result<u64, EngineError>;
}

/// Result of one [`run`] call.
#[derive(Clone, Debug)]
pub struct RunReport {
    /// Chunks processed and delivered in order.
    pub n_chunks: u64,
    /// Value returned by the sink's `finish`.
    pub sink_report: u64,
    /// Total wall time.
    pub elapsed: Duration,
    /// Aggregate stage metrics (§7.10).
    pub metrics: StageMetrics,
}

/// Runs the pipeline: source → workers → ordered sink, honoring budgets and
/// abort (`docs/architecture-spec.md` §7.1).
///
/// Determinism contract (§7.6): chunk boundaries come from the input only;
/// partials merge strictly in chunk order; seam records are stitched in
/// between. Output is independent of thread count — locked by tests.
///
/// # Errors
///
/// The first error recorded by any stage, [`EngineError::Aborted`] on
/// external abort, [`EngineError::MemBudgetExceeded`] when a single chunk
/// exceeds the whole memory budget.
pub fn run<J, S>(
    source: &mut dyn ChunkSource,
    job: &J,
    sink: S,
    budget: &ThreadBudget,
    abort: &Abort,
) -> Result<RunReport, EngineError>
where
    J: ChunkJob + 'static,
    S: OrderedSink<J::Partial> + 'static,
{
    let started = Instant::now();
    // Caller and engine share one channel: an internal error triggers the
    // caller's handle, an external trigger aborts the stages — no watcher.
    let first = Arc::new(FirstError::new(abort.handle.clone()));
    let outcome = run_parallel(source, job, sink, budget, &first, &abort.signal);
    // The recorded cause outranks the Aborted consequence it produced.
    if let Some(err) = FirstError::take_shared(&first) {
        return Err(err);
    }
    let (n_chunks, sink_report, metrics) = outcome?;
    Ok(RunReport {
        n_chunks,
        sink_report,
        elapsed: started.elapsed(),
        metrics,
    })
}

/// Shared state handed to the pipeline threads.
struct Shared {
    in_flight: AtomicU64,
    in_flight_peak: AtomicU64,
}

impl Shared {
    fn record_in_flight(&self, bytes: u64) {
        let now = self.in_flight.fetch_add(bytes, Ordering::AcqRel) + bytes;
        self.in_flight_peak.fetch_max(now, Ordering::AcqRel);
    }

    fn release_in_flight(&self, bytes: u64) {
        self.in_flight.fetch_sub(bytes, Ordering::AcqRel);
    }
}

/// Parallel mode: reader thread, `workers()` compute workers, the caller's
/// thread reduces in chunk order (§7.3 thread layout).
fn run_parallel<J, S>(
    source: &mut dyn ChunkSource,
    job: &J,
    sink: S,
    budget: &ThreadBudget,
    first: &Arc<FirstError>,
    abort: &AbortSignal,
) -> Result<(u64, u64, StageMetrics), EngineError>
where
    J: ChunkJob + 'static,
    S: OrderedSink<J::Partial> + 'static,
{
    let n_workers = budget.workers();
    let tickets = budget.tickets();
    let (chunk_tx, chunk_rx) = bounded::<InFlightChunk>(tickets);
    let (result_tx, result_rx) = bounded::<WorkerMsg<J::Partial>>(tickets);
    // Memory backpressure (03 §4.2): in-flight chunks hold permits sized by
    // their estimated footprint; the reader blocks while the pool is empty.
    let permits = Arc::new(BytePermits::new(budget.mem.bytes));
    // Real tickets (03 §4.1, fixes 01-review D1): the reader acquires one
    // permit *before* cutting a chunk, and the reducer returns it only after
    // the chunk is delivered in order. The number of chunks cut but not yet
    // delivered is therefore bounded by `tickets`, so every in-flight id
    // stays inside `[expected, expected + tickets)` and the reorder ring of
    // the same capacity can never overflow.
    let (ticket_tx, ticket_rx) = bounded::<()>(tickets);
    for _ in 0..tickets {
        // Capacity equals the token count, so these never block.
        let _ = ticket_tx.send(());
    }
    let shared = Arc::new(Shared {
        in_flight: AtomicU64::new(0),
        in_flight_peak: AtomicU64::new(0),
    });
    let window = tickets;
    let worker_busy_ns = Arc::new(AtomicU64::new(0));
    // doc_v4 03 A4-5: reader and reducer cumulative times, so "reader-bound"
    // is measured instead of inferred.
    let reader_ns = Arc::new(AtomicU64::new(0));
    let reduce_ns = Arc::new(AtomicU64::new(0));
    let chunks_read = Arc::new(AtomicU64::new(0));

    let scope_result = std::thread::scope(|scope| -> Result<_, EngineError> {
        // Reader: acquire a ticket, then cut one chunk; sends block when
        // the ticket window is full. Every blocking point races abort, and
        // an abort/error breaks the loop so the sender disconnects and
        // unblocks everyone downstream.
        scope.spawn({
            let chunk_tx = chunk_tx.clone();
            let first = Arc::clone(first);
            let shared = Arc::clone(&shared);
            let ticket_rx = ticket_rx.clone();
            let permits = Arc::clone(&permits);
            let reader_ns = Arc::clone(&reader_ns);
            let chunks_read = Arc::clone(&chunks_read);
            move || {
                loop {
                    if first.failed() {
                        break;
                    }
                    let ticket = crossbeam_channel::select! {
                        recv(abort.receiver()) -> _ => false,
                        recv(ticket_rx) -> token => token.is_ok(),
                    };
                    if !ticket || first.failed() {
                        break;
                    }
                    let read_started = Instant::now();
                    let chunk = source.next_chunk();
                    reader_ns.fetch_add(
                        read_started.elapsed().as_nanos() as u64,
                        Ordering::AcqRel,
                    );
                    let chunk = match chunk {
                        Ok(Some(chunk)) => {
                            chunks_read.fetch_add(1, Ordering::AcqRel);
                            chunk
                        }
                        Ok(None) => break,
                        Err(err) => {
                            first.record(EngineError::from(err));
                            break;
                        }
                    };
                    // Memory backpressure: acquire the footprint before the
                    // chunk enters the pipeline. One chunk larger than the
                    // whole budget cannot continue (03 §4.2).
                    let footprint = budget.mem.estimate(chunk.compressed.len() as u64);
                    if footprint > budget.mem.bytes {
                        first.record(EngineError::MemBudgetExceeded {
                            peak_bytes: footprint,
                            budget_bytes: budget.mem.bytes,
                        });
                        break;
                    }
                    if !permits.acquire(footprint, abort) {
                        break;
                    }
                    if first.failed() {
                        permits.release(footprint);
                        break;
                    }
                    shared.record_in_flight(footprint);
                    let sent = crossbeam_channel::select! {
                        recv(abort.receiver()) -> _ => {
                            permits.release(footprint);
                            false
                        }
                        send(chunk_tx, InFlightChunk { chunk, footprint }) -> _ => true,
                    };
                    if !sent || first.failed() {
                        break;
                    }
                }
            }
        });
        drop(chunk_tx);

        // Workers: inflate + job, fused per chunk (§7.2). Both blocking
        // points (chunk recv, result send) race abort.
        for _ in 0..n_workers {
            scope.spawn({
                let chunk_rx = chunk_rx.clone();
                let result_tx = result_tx.clone();
                let first = Arc::clone(first);
                let worker_busy_ns = Arc::clone(&worker_busy_ns);
                move || {
                    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
                    let mut inflate_scratch = Vec::new();
                    let mut scratch = J::Scratch::default();
                    loop {
                        if first.failed() {
                            break;
                        }
                        let in_flight = crossbeam_channel::select! {
                            recv(abort.receiver()) -> _ => break,
                            recv(chunk_rx) -> msg => match msg {
                                Ok(in_flight) => in_flight,
                                Err(_) => break,
                            },
                        };
                        let chunk = in_flight.chunk;
                        let started = Instant::now();
                        let framing = job.framing();
                        let outcome = catch_unwind(AssertUnwindSafe(|| {
                            let view =
                                inflate_chunk(&chunk, &mut inflater, &mut inflate_scratch)?;
                            let (body, head, tail) = match framing {
                                RecordFraming::Newline => (
                                    view.body_bytes(view.text),
                                    view.head_bytes(view.text).map(Vec::from),
                                    view.tail_bytes(view.text).map(Vec::from),
                                ),
                                // The job walks the length prefixes itself.
                                RecordFraming::Raw => (view.text, None, None),
                            };
                            let partial = job.run_chunk(&chunk, body, &mut scratch)?;
                            Ok::<_, EngineError>((partial, head, tail))
                        }));
                        worker_busy_ns.fetch_add(
                            started.elapsed().as_nanos() as u64,
                            Ordering::AcqRel,
                        );
                        let (partial, head, tail, failure) = match outcome {
                            Ok(Ok((partial, head, tail))) => (partial, head, tail, None),
                            Ok(Err(err)) => (J::Partial::default(), None, None, Some(err)),
                            Err(panic) => (
                                J::Partial::default(),
                                None,
                                None,
                                Some(EngineError::WorkerPanicked {
                                    chunk: chunk.id,
                                    message: panic_message(&panic),
                                }),
                            ),
                        };
                        if let Some(err) = failure {
                            first.record(err);
                            break;
                        }
                        let sent = crossbeam_channel::select! {
                            recv(abort.receiver()) -> _ => false,
                            send(result_tx, WorkerMsg {
                                id: chunk.id,
                                footprint: in_flight.footprint,
                                partial,
                                head,
                                tail,
                            }) -> _ => true,
                        };
                        if !sent {
                            break;
                        }
                    }
                }
            });
        }
        drop(result_tx);
        // Only worker clones receive chunks now: when the last worker exits
        // (error, abort, or drain), the chunk channel disconnects and the
        // reader's pending send fails instead of blocking forever.
        drop(chunk_rx);

        // Reducer (this thread): deliver in chunk order, chain-stitching
        // seam records between chunks (§6.2 step 4, §7.6 rule 3). Errors
        // return immediately — the dropped senders/receivers unblock every
        // other stage, and `run` reports the recorded cause.
        let mut ring: ReorderBuffer<WorkerMsg<J::Partial>> = ReorderBuffer::new(window);
        let mut fragment: Option<Vec<u8>> = None;
        let mut seam_scratch = J::Scratch::default();
        let mut delivered: u64 = 0;
        let mut sink = sink;
        'outer: loop {
            if first.failed() {
                drop(result_rx);
                return Err(EngineError::Aborted);
            }
            while ring.ready() {
                let Some((id, mut msg)) = ring.take_ready() else {
                    break;
                };
                fragment = match stitch_seam(
                    fragment,
                    msg.head.take(),
                    job,
                    &mut seam_scratch,
                    &mut sink,
                    id,
                ) {
                    Ok(fragment) => fragment,
                    Err(err) => {
                        drop(result_rx);
                        return Err(err);
                    }
                };
                let accept_started = Instant::now();
                let accepted = sink.accept(id, &mut msg.partial);
                reduce_ns.fetch_add(
                    accept_started.elapsed().as_nanos() as u64,
                    Ordering::AcqRel,
                );
                if let Err(err) = accepted {
                    drop(result_rx);
                    return Err(err);
                }
                shared.release_in_flight(msg.footprint);
                permits.release(msg.footprint);
                // Delivered in order: the chunk's ticket returns to the
                // reader (03 §4.1). Never blocks — returns ≤ acquisitions.
                let _ = ticket_tx.send(());
                fragment = finish_fragment(fragment, msg.tail.take());
                delivered += 1;
            }
            if first.failed() {
                drop(result_rx);
                return Err(EngineError::Aborted);
            }
            let msg = crossbeam_channel::select! {
                recv(abort.receiver()) -> _ => None,
                recv(result_rx) -> msg => match msg {
                    Ok(msg) => Some(msg),
                    Err(_) => {
                        // All senders gone: deliver whatever is already in
                        // the ring, then finish.
                        while let Some((id, mut msg)) = ring.take_ready() {
                            fragment = stitch_seam(
                                fragment,
                                msg.head.take(),
                                job,
                                &mut seam_scratch,
                                &mut sink,
                                id,
                            )?;
                            let accept_started = Instant::now();
                            let accepted = sink.accept(id, &mut msg.partial);
                            reduce_ns.fetch_add(
                                accept_started.elapsed().as_nanos() as u64,
                                Ordering::AcqRel,
                            );
                            accepted?;
                            shared.release_in_flight(msg.footprint);
                            permits.release(msg.footprint);
                            let _ = ticket_tx.send(());
                            fragment = finish_fragment(fragment, msg.tail.take());
                            delivered += 1;
                        }
                        break 'outer;
                    }
                },
            };
            match msg {
                Some(msg) => {
                    ring.insert(msg.id, msg);
                }
                None => {
                    drop(result_rx);
                    return Err(EngineError::Aborted);
                }
            }
        }
        // Draining finished; dropping the receiver unblocks any worker
        // still sending, so the scope join below cannot deadlock.
        drop(result_rx);
        if first.failed() {
            return Err(EngineError::Aborted);
        }
        let sink_report = sink.finish()?;
        Ok((delivered, sink_report))
    });

    let (delivered, sink_report) = scope_result?;
    Ok((
        delivered,
        sink_report,
        StageMetrics {
            n_chunks_read: chunks_read.load(Ordering::Acquire),
            worker_busy: Duration::from_nanos(worker_busy_ns.load(Ordering::Acquire)),
            in_flight_peak_bytes: shared.in_flight_peak.load(Ordering::Acquire),
            reader_time: Duration::from_nanos(reader_ns.load(Ordering::Acquire)),
            reduce_time: Duration::from_nanos(reduce_ns.load(Ordering::Acquire)),
            chain_wait: job.wait_time(),
        },
    ))
}

/// A worker result on its way to the reducer.
struct WorkerMsg<P> {
    id: ChunkId,
    /// The memory permits this chunk holds; the reducer releases them.
    footprint: u64,
    partial: P,
    head: Option<Vec<u8>>,
    tail: Option<Vec<u8>>,
}

/// A chunk in the pipeline together with its memory footprint.
struct InFlightChunk {
    chunk: ChunkBytes,
    footprint: u64,
}

/// Chain-stitches the seam record between the previous fragment and this
/// chunk's head, running `job.run_seam` when the line completes (§6.2).
fn stitch_seam<J, S>(
    fragment: Option<Vec<u8>>,
    head: Option<Vec<u8>>,
    job: &J,
    seam_scratch: &mut J::Scratch,
    sink: &mut S,
    id: ChunkId,
) -> Result<Option<Vec<u8>>, EngineError>
where
    J: ChunkJob,
    S: OrderedSink<J::Partial>,
{
    let Some(head) = head else {
        // First chunk: nothing to stitch (its head is complete by
        // construction); any prior fragment would be a bug.
        debug_assert!(fragment.is_none(), "invariant: fragment before the first chunk");
        return Ok(None);
    };
    let mut seam = fragment.unwrap_or_default();
    seam.extend_from_slice(&head);
    if seam.ends_with(b"\n") {
        let mut partial = job.run_seam(id, &seam, seam_scratch)?;
        sink.accept(id, &mut partial)?;
        Ok(None)
    } else {
        // The record spans further chunks: keep chaining.
        Ok(Some(seam))
    }
}

/// Merges a chunk's tail into the running fragment.
fn finish_fragment(fragment: Option<Vec<u8>>, tail: Option<Vec<u8>>) -> Option<Vec<u8>> {
    match (fragment, tail) {
        (Some(mut acc), Some(tail)) => {
            acc.extend_from_slice(&tail);
            Some(acc)
        }
        (Some(acc), None) => Some(acc),
        (None, Some(tail)) => Some(tail),
        (None, None) => None,
    }
}

/// Extracts a printable message from a panic payload.
fn panic_message(panic: &Box<dyn std::any::Any + Send>) -> String {
    if let Some(s) = panic.downcast_ref::<&str>() {
        (*s).to_string()
    } else if let Some(s) = panic.downcast_ref::<String>() {
        s.clone()
    } else {
        "unknown panic".to_string()
    }
}
