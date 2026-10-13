//! The chunk pipeline skeleton (`docs/spec/architecture.md` §7.1–§7.6,
//! `docs/spec/v6_y5/README.md` §2): reader → bounded chunk channel (ticket
//! window) → worker pool (fused inflate + job) → ordered reducer with seam
//! chain-stitching; `threads <= 2` takes the sequential path (S4).
//!
//! Threads never exchange per-record heap objects — only compressed chunks
//! and small partials (§7.2). Every blocking point selects on the abort
//! signal (§7.5), and in-flight compressed bytes are accounted against the
//! memory budget (§7.7); the permits live in a drop guard, so no exit path
//! can leak them (S2).

// Casts in this file are bounded by construction: Elapsed nanoseconds are clamped to u64 in practice (584 years).
#![allow(clippy::cast_possible_truncation)]

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use crossbeam_channel::bounded;
use reed_core::ChunkId;
use reed_io::chunk::framing::RecordFraming;
use reed_io::chunk::{ChunkBytes, ChunkSource, inflate_chunk};

use crate::abort::{Abort, AbortSignal, FirstError};
use crate::budget::{BytePermits, ThreadBudget};
use crate::error::EngineError;
use crate::metrics::StageMetrics;
use crate::reorder::ReorderBuffer;

/// The largest record the engine reassembles across chunk seams. Shared with
/// the BCF frame chain, which enforces the same ceiling on carried prefixes
/// (`docs/rounds/v6/03-architecture.md` A6-2).
pub const MAX_RECORD_BYTES: usize = 1 << 30;

/// Message for an input whose last record has no terminating newline. Shared
/// with the BCF path (`reed/src/bcf_input.rs`) so every input format reports
/// the same failure (doc_v9 D9-3).
pub const TRUNCATED_INPUT: &str = "input ends inside a record (truncated file)";

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

/// A fused chunk computation (`docs/spec/architecture.md` §7.1).
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

    /// Cumulative output-compression time this sink observed (A6-4). The
    /// engine reads it once, before `finish` consumes the sink; sinks that
    /// do not compress leave it at zero.
    fn deflate_time(&self) -> Duration {
        Duration::ZERO
    }
}

/// Reassembles the record that spans a chunk seam from its fragments, with
/// a hard ceiling: a pathological newline-free input must fail the run, not
/// grow the reducer's buffer to the file size (A6-2).
#[derive(Default)]
pub struct NewlineStitcher {
    fragment: Vec<u8>,
}

impl NewlineStitcher {
    fn grow(&mut self, bytes: &[u8]) -> Result<(), EngineError> {
        if self.fragment.len() + bytes.len() > MAX_RECORD_BYTES {
            return Err(EngineError::job(format!(
                "record exceeds {MAX_RECORD_BYTES} bytes across chunk seams"
            )));
        }
        self.fragment.extend_from_slice(bytes);
        Ok(())
    }

    /// Absorbs a chunk's head fragment; `Some(line)` when the boundary line
    /// closes here. A `None` head means the first chunk (or none pending).
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] when the accumulated fragment passes
    /// [`MAX_RECORD_BYTES`].
    pub fn push_head(&mut self, head: Option<&[u8]>) -> Result<Option<Vec<u8>>, EngineError> {
        let Some(head) = head else {
            debug_assert!(
                self.fragment.is_empty(),
                "invariant: fragment before the first chunk"
            );
            return Ok(None);
        };
        self.grow(head)?;
        Ok(self
            .fragment
            .ends_with(b"\n")
            .then(|| std::mem::take(&mut self.fragment)))
    }

    /// Merges a chunk's tail fragment into the running fragment.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] when the accumulated fragment passes
    /// [`MAX_RECORD_BYTES`].
    pub fn push_tail(&mut self, tail: Option<&[u8]>) -> Result<(), EngineError> {
        if let Some(tail) = tail {
            self.grow(tail)?;
        }
        Ok(())
    }

    /// Consumes the stitcher and returns the fragment left after the last
    /// chunk, or `None` when the stream ended on a record boundary.
    ///
    /// A `Some` result is an unterminated final record: the caller rejects it
    /// (doc_v9 D9-3). This is the stream's end — it is **not**
    /// [`push_head(None)`](Self::push_head), which only means "this chunk has
    /// no head fragment".
    #[must_use]
    pub fn finish(self) -> Option<Vec<u8>> {
        (!self.fragment.is_empty()).then_some(self.fragment)
    }
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
/// abort (`docs/spec/architecture.md` §7.1). `threads <= 2` takes the
/// single-threaded sequential path, whose output is byte-identical (§7.3).
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
    let outcome = if budget.is_sequential() {
        run_sequential(source, job, sink, budget, &abort.signal)
    } else {
        run_parallel(source, job, sink, budget, &first, abort)
    };
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

    fn in_flight_final(&self) -> u64 {
        self.in_flight.load(Ordering::Acquire)
    }
}

/// Owns one in-flight chunk's memory permit and its in-flight accounting.
/// Dropped, never explicitly released: on in-order delivery, on worker
/// failure, or when an aborted pipeline drops the pending message (S2).
struct PermitGuard {
    permits: Arc<BytePermits>,
    shared: Arc<Shared>,
    bytes: u64,
}

impl PermitGuard {
    fn new(permits: Arc<BytePermits>, shared: Arc<Shared>, bytes: u64) -> Self {
        shared.record_in_flight(bytes);
        Self {
            permits,
            shared,
            bytes,
        }
    }
}

impl Drop for PermitGuard {
    fn drop(&mut self) {
        self.permits.release(self.bytes);
        self.shared.release_in_flight(self.bytes);
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
    abort: &Abort,
) -> Result<(u64, u64, StageMetrics), EngineError>
where
    J: ChunkJob + 'static,
    S: OrderedSink<J::Partial> + 'static,
{
    let signal = &abort.signal;
    let n_workers = budget.workers();
    let tickets = budget.tickets();
    let (chunk_tx, chunk_rx) = bounded::<InFlightChunk>(tickets);
    let (result_tx, result_rx) = bounded::<WorkerMsg<J::Partial>>(tickets);
    // Memory backpressure (03 §4.2): in-flight chunks hold permits sized by
    // their estimated footprint; the reader blocks while the pool is empty.
    let permits = Arc::new(BytePermits::new(budget.mem.bytes));
    // The abort handle wakes these waiters directly (no polling, §7.5).
    abort.handle.register(&permits);
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
    // is measured instead of inferred. A6-4 splits the wait and inflate/job
    // shares out of the busy totals.
    let reader_ns = Arc::new(AtomicU64::new(0));
    let reduce_ns = Arc::new(AtomicU64::new(0));
    let reader_wait_ns = Arc::new(AtomicU64::new(0));
    let reducer_wait_ns = Arc::new(AtomicU64::new(0));
    let inflate_ns = Arc::new(AtomicU64::new(0));
    let job_ns = Arc::new(AtomicU64::new(0));
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
            let reader_wait_ns = Arc::clone(&reader_wait_ns);
            let chunks_read = Arc::clone(&chunks_read);
            move || {
                loop {
                    if first.failed() {
                        break;
                    }
                    let wait_started = Instant::now();
                    let ticket = crossbeam_channel::select! {
                        recv(signal.receiver()) -> _ => false,
                        recv(ticket_rx) -> token => token.is_ok(),
                    };
                    reader_wait_ns
                        .fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::AcqRel);
                    if !ticket || first.failed() {
                        break;
                    }
                    let read_started = Instant::now();
                    let chunk = source.next_chunk();
                    reader_ns.fetch_add(read_started.elapsed().as_nanos() as u64, Ordering::AcqRel);
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
                    let wait_started = Instant::now();
                    let acquired = permits.acquire(footprint, signal);
                    reader_wait_ns
                        .fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::AcqRel);
                    if !acquired || first.failed() {
                        break;
                    }
                    let guard =
                        PermitGuard::new(Arc::clone(&permits), Arc::clone(&shared), footprint);
                    let wait_started = Instant::now();
                    let sent = crossbeam_channel::select! {
                        recv(signal.receiver()) -> _ => false,
                        send(chunk_tx, InFlightChunk { chunk, guard }) -> res => res.is_ok(),
                    };
                    reader_wait_ns
                        .fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::AcqRel);
                    if !sent || first.failed() {
                        // Abort won, or no worker is left to receive: the
                        // guard is a local and drops here either way (D9-2).
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
                let inflate_ns = Arc::clone(&inflate_ns);
                let job_ns = Arc::clone(&job_ns);
                move || {
                    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
                    let mut inflate_scratch = Vec::new();
                    let mut scratch = J::Scratch::default();
                    loop {
                        if first.failed() {
                            break;
                        }
                        let in_flight = crossbeam_channel::select! {
                            recv(signal.receiver()) -> _ => break,
                            recv(chunk_rx) -> msg => match msg {
                                Ok(in_flight) => in_flight,
                                Err(_) => break,
                            },
                        };
                        let InFlightChunk { chunk, guard } = in_flight;
                        let started = Instant::now();
                        let framing = job.framing();
                        let outcome = catch_unwind(AssertUnwindSafe(|| {
                            let inflate_started = Instant::now();
                            let view = inflate_chunk(&chunk, &mut inflater, &mut inflate_scratch)?;
                            let inflate_elapsed = inflate_started.elapsed();
                            let (body, head, tail) = match framing {
                                RecordFraming::Newline => (
                                    view.body_bytes(view.text),
                                    view.head_bytes(view.text).map(Vec::from),
                                    view.tail_bytes(view.text).map(Vec::from),
                                ),
                                // The job walks the length prefixes itself.
                                RecordFraming::Raw => (view.text, None, None),
                            };
                            let job_started = Instant::now();
                            let partial = job.run_chunk(&chunk, body, &mut scratch)?;
                            let job_elapsed = job_started.elapsed();
                            Ok::<_, EngineError>((
                                partial,
                                head,
                                tail,
                                inflate_elapsed,
                                job_elapsed,
                            ))
                        }));
                        worker_busy_ns
                            .fetch_add(started.elapsed().as_nanos() as u64, Ordering::AcqRel);
                        let (partial, head, tail, failure) = match outcome {
                            Ok(Ok((partial, head, tail, inflate_elapsed, job_elapsed))) => {
                                inflate_ns
                                    .fetch_add(inflate_elapsed.as_nanos() as u64, Ordering::AcqRel);
                                job_ns.fetch_add(job_elapsed.as_nanos() as u64, Ordering::AcqRel);
                                (partial, head, tail, None)
                            }
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
                            break; // `guard` drops: permit + in_flight returned.
                        }
                        let sent = crossbeam_channel::select! {
                            recv(signal.receiver()) -> _ => false,
                            send(result_tx, WorkerMsg {
                                id: chunk.id,
                                guard,
                                partial,
                                head,
                                tail,
                            }) -> res => res.is_ok(),
                        };
                        if !sent {
                            // Aborted, or the reducer dropped its receiver: the
                            // unsent message drops with its guard (doc_v9 D9-2).
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
        // seam records between chunks (§6.2 step 4, §7.6 rule 3).
        //
        // A reducer-side error or panic must exit through `bail`: recording
        // into `first` is what fires the abort that unblocks the reader's
        // ticket wait — without it the scope join deadlocks instead of
        // failing the run (doc_v9 D9-1/D9-1c).
        let mut ring: ReorderBuffer<WorkerMsg<J::Partial>> = ReorderBuffer::new(window);
        let mut stitcher = NewlineStitcher::default();
        let mut seam_scratch = J::Scratch::default();
        let mut delivered: u64 = 0;
        let mut sink = sink;
        // One delivery step, shared by the hot loop and the drain path.
        macro_rules! deliver {
            ($id:expr, $msg:expr) => {{
                let mut msg = $msg;
                if let Some(seam) = stitcher.push_head(msg.head.take().as_deref())? {
                    // The seam runs on the reducer, not a worker: count it in
                    // job_time *and* worker_busy (so `worker_busy >= inflate +
                    // job` still holds) or the work is invisible to the
                    // metrics (doc_v9 D9-6).
                    let seam_job_started = Instant::now();
                    let mut partial = job.run_seam($id, &seam, &mut seam_scratch)?;
                    let seam_ns = seam_job_started.elapsed().as_nanos() as u64;
                    job_ns.fetch_add(seam_ns, Ordering::AcqRel);
                    worker_busy_ns.fetch_add(seam_ns, Ordering::AcqRel);
                    let seam_accept_started = Instant::now();
                    sink.accept($id, &mut partial)?;
                    reduce_ns.fetch_add(
                        seam_accept_started.elapsed().as_nanos() as u64,
                        Ordering::AcqRel,
                    );
                }
                let accept_started = Instant::now();
                let accepted = sink.accept($id, &mut msg.partial);
                reduce_ns.fetch_add(accept_started.elapsed().as_nanos() as u64, Ordering::AcqRel);
                accepted?;
                let tail = msg.tail.take();
                drop(msg); // permits + in_flight return here (S2)
                let _ = ticket_tx.send(());
                stitcher.push_tail(tail.as_deref())?;
                delivered += 1;
            }};
        }
        // The whole reducer body runs under `catch_unwind` (as the workers'
        // per-chunk work does): a panic anywhere — `accept`, `run_seam`,
        // `finish` — becomes `ReducerPanicked` instead of unwinding out of
        // the scope, which would deadlock the join like any reducer error.
        let reducer = catch_unwind(AssertUnwindSafe(|| -> Result<_, EngineError> {
            'outer: loop {
                if first.failed() {
                    return Err(EngineError::Aborted);
                }
                while ring.ready() {
                    let Some((id, msg)) = ring.take_ready() else {
                        break;
                    };
                    deliver!(id, msg);
                }
                if first.failed() {
                    return Err(EngineError::Aborted);
                }
                let wait_started = Instant::now();
                let msg = crossbeam_channel::select! {
                    recv(signal.receiver()) -> _ => None,
                    recv(result_rx) -> msg => match msg {
                        Ok(msg) => Some(msg),
                        Err(_) => {
                            // All senders gone: deliver whatever is already in
                            // the ring, then finish.
                            reducer_wait_ns
                                .fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::AcqRel);
                            while let Some((id, msg)) = ring.take_ready() {
                                deliver!(id, msg);
                            }
                            break 'outer;
                        }
                    },
                };
                reducer_wait_ns
                    .fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::AcqRel);
                match msg {
                    Some(msg) => {
                        ring.insert(msg.id, msg);
                    }
                    None => return Err(EngineError::Aborted),
                }
            }
            // Draining finished; dropping the receiver unblocks any worker
            // still sending, so the scope join below cannot deadlock.
            drop(result_rx);
            if first.failed() {
                return Err(EngineError::Aborted);
            }
            let deflate_time = sink.deflate_time();
            // A tail fragment left after the last chunk is an unterminated
            // final record: reject it exactly as the BCF path does (D9-3).
            if stitcher.finish().is_some() {
                return Err(EngineError::job(TRUNCATED_INPUT));
            }
            let sink_report = sink.finish()?;
            Ok((delivered, sink_report, deflate_time))
        }));
        match reducer {
            Ok(Ok(report)) => Ok(report),
            Ok(Err(err)) => Err(bail(first, err)),
            Err(panic) => Err(bail(
                first,
                EngineError::ReducerPanicked {
                    // `delivered` is how many chunks landed before the panic.
                    chunk: ChunkId(delivered),
                    message: panic_message(&panic),
                },
            )),
        }
    });

    let (delivered, sink_report, deflate_time) = scope_result?;
    Ok((
        delivered,
        sink_report,
        StageMetrics {
            n_chunks_read: chunks_read.load(Ordering::Acquire),
            worker_busy: Duration::from_nanos(worker_busy_ns.load(Ordering::Acquire)),
            in_flight_peak_bytes: shared.in_flight_peak.load(Ordering::Acquire),
            in_flight_final_bytes: shared.in_flight_final(),
            reader_time: Duration::from_nanos(reader_ns.load(Ordering::Acquire)),
            reduce_time: Duration::from_nanos(reduce_ns.load(Ordering::Acquire)),
            chain_wait: job.wait_time(),
            reader_wait: Duration::from_nanos(reader_wait_ns.load(Ordering::Acquire)),
            reducer_wait: Duration::from_nanos(reducer_wait_ns.load(Ordering::Acquire)),
            inflate_time: Duration::from_nanos(inflate_ns.load(Ordering::Acquire)),
            job_time: Duration::from_nanos(job_ns.load(Ordering::Acquire)),
            deflate_time,
        },
    ))
}

/// Sequential mode (§7.3): one thread reads, inflates, computes and writes.
/// Output must equal the parallel path byte for byte — the same stitcher and
/// the same sink, no channels in between (S4).
fn run_sequential<J, S>(
    source: &mut dyn ChunkSource,
    job: &J,
    sink: S,
    budget: &ThreadBudget,
    abort: &AbortSignal,
) -> Result<(u64, u64, StageMetrics), EngineError>
where
    J: ChunkJob + 'static,
    S: OrderedSink<J::Partial> + 'static,
{
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut inflate_scratch = Vec::new();
    let mut scratch = J::Scratch::default();
    let mut seam_scratch = J::Scratch::default();
    let mut stitcher = NewlineStitcher::default();
    let mut sink = sink;
    let mut reader_time = Duration::ZERO;
    let mut inflate_time = Duration::ZERO;
    let mut job_time = Duration::ZERO;
    let mut reduce_time = Duration::ZERO;
    let mut chunks_read = 0_u64;
    let mut delivered = 0_u64;
    let mut in_flight_peak = 0_u64;

    loop {
        if abort.is_aborted() {
            return Err(EngineError::Aborted);
        }
        let started = Instant::now();
        let chunk = source.next_chunk()?;
        reader_time += started.elapsed();
        let Some(chunk) = chunk else {
            break;
        };
        chunks_read += 1;

        let footprint = budget.mem.estimate(chunk.compressed.len() as u64);
        if footprint > budget.mem.bytes {
            return Err(EngineError::MemBudgetExceeded {
                peak_bytes: footprint,
                budget_bytes: budget.mem.bytes,
            });
        }
        in_flight_peak = in_flight_peak.max(footprint);

        let started = Instant::now();
        let view = inflate_chunk(&chunk, &mut inflater, &mut inflate_scratch)?;
        inflate_time += started.elapsed();

        match job.framing() {
            RecordFraming::Newline => {
                if let Some(seam) = stitcher.push_head(view.head_bytes(view.text))? {
                    let started = Instant::now();
                    let mut partial = job.run_seam(chunk.id, &seam, &mut seam_scratch)?;
                    job_time += started.elapsed();
                    let started = Instant::now();
                    sink.accept(chunk.id, &mut partial)?;
                    reduce_time += started.elapsed();
                }
                let started = Instant::now();
                let mut partial =
                    job.run_chunk(&chunk, view.body_bytes(view.text), &mut scratch)?;
                job_time += started.elapsed();
                let started = Instant::now();
                sink.accept(chunk.id, &mut partial)?;
                reduce_time += started.elapsed();
                stitcher.push_tail(view.tail_bytes(view.text))?;
            }
            // The job walks the length prefixes; no seam assembly here.
            RecordFraming::Raw => {
                let started = Instant::now();
                let mut partial = job.run_chunk(&chunk, view.text, &mut scratch)?;
                job_time += started.elapsed();
                let started = Instant::now();
                sink.accept(chunk.id, &mut partial)?;
                reduce_time += started.elapsed();
            }
        }
        delivered += 1;
    }

    // A tail fragment left after the last chunk is an unterminated final
    // record: reject it as the BCF path does (doc_v9 D9-3). Only the Newline
    // framing produces fragments; Raw jobs own their own framing.
    if stitcher.finish().is_some() {
        return Err(EngineError::job(TRUNCATED_INPUT));
    }
    let deflate_time = sink.deflate_time();
    let sink_report = sink.finish()?;
    Ok((
        delivered,
        sink_report,
        StageMetrics {
            n_chunks_read: chunks_read,
            worker_busy: inflate_time + job_time,
            // One chunk at a time here: the peak is the largest footprint
            // that was in flight (doc_v9 D9-5).
            in_flight_peak_bytes: in_flight_peak,
            in_flight_final_bytes: 0,
            reader_time,
            reduce_time,
            chain_wait: job.wait_time(),
            reader_wait: Duration::ZERO,
            reducer_wait: Duration::ZERO,
            inflate_time,
            job_time,
            deflate_time,
        },
    ))
}

/// A worker result on its way to the reducer.
struct WorkerMsg<P> {
    id: ChunkId,
    /// The memory permits this chunk holds; the reducer drops them after
    /// in-order delivery. Never read — held for its `Drop` (S2).
    #[allow(dead_code, reason = "held for its Drop impl")]
    guard: PermitGuard,
    partial: P,
    head: Option<Vec<u8>>,
    tail: Option<Vec<u8>>,
}

/// A chunk in the pipeline together with its memory permit.
struct InFlightChunk {
    chunk: ChunkBytes,
    guard: PermitGuard,
}

/// The reducer's single error exit (doc_v9 §2.2, 规约 W2): record `err` so
/// `FirstError` triggers the abort, then return the `Aborted` consequence.
///
/// The reader blocks on its *ticket* wait, which only the abort signal can
/// unblock (the ticket sender lives until `run_parallel` returns, i.e. until
/// the scope join completes) — so recording is what prevents the deadlock;
/// `run` then reports the recorded cause over this `Aborted` (`run`).
fn bail(first: &FirstError, err: EngineError) -> EngineError {
    first.record(err);
    EngineError::Aborted
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

#[cfg(test)]
mod tests {
    use super::*;

    /// The stitcher reassembles a record split across many chunks and
    /// enforces the ceiling instead of growing without bound (A6-2).
    #[test]
    fn stitcher_reassembles_and_caps_fragments() {
        let mut stitcher = NewlineStitcher::default();
        // No fragment before the first chunk.
        assert!(stitcher.push_head(None).unwrap().is_none());
        // "abc" tail, then "def\n" head: the line completes in the head.
        stitcher.push_tail(Some(b"abc")).unwrap();
        let line = stitcher.push_head(Some(b"def\n")).unwrap().unwrap();
        assert_eq!(line, b"abcdef\n");
        // An unterminated head keeps chaining.
        stitcher.push_tail(Some(b"x")).unwrap();
        assert!(stitcher.push_head(Some(b"y")).unwrap().is_none());
        // `finish` (not `push_head(None)`) is the stream's end: the leftover
        // fragment is an unterminated final record (doc_v9 D9-3/D9-8).
        assert_eq!(stitcher.finish().unwrap(), b"xy");
    }

    #[test]
    fn stitcher_rejects_a_record_past_the_ceiling() {
        let mut stitcher = NewlineStitcher {
            fragment: vec![0; MAX_RECORD_BYTES],
        };
        let err = stitcher
            .push_head(Some(b"x"))
            .expect_err("a fragment past the ceiling must fail");
        assert!(matches!(err, EngineError::Job(_)), "got {err:?}");
        assert!(err.to_string().contains("across chunk seams"), "{err}");
    }

    /// S2 gate: the guard is the only release path — dropping it returns the
    /// bytes to the pool and zeroes the in-flight accounting.
    #[test]
    fn permit_guard_returns_its_bytes_on_drop() {
        let permits = Arc::new(BytePermits::new(100));
        let shared = Arc::new(Shared {
            in_flight: AtomicU64::new(0),
            in_flight_peak: AtomicU64::new(0),
        });
        let guard = PermitGuard::new(Arc::clone(&permits), Arc::clone(&shared), 40);
        assert_eq!(shared.in_flight.load(Ordering::Acquire), 40);
        assert_eq!(shared.in_flight_peak.load(Ordering::Acquire), 40);
        drop(guard);
        assert_eq!(shared.in_flight.load(Ordering::Acquire), 0);
        // The 40 bytes are back: a second acquisition of the full 100 works.
        assert!(permits.acquire(100, &crate::abort::abort().signal));
    }
}
