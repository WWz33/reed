//! Engine acceptance tests (`docs/architecture-spec.md` §19 M3): ordered
//! reduction under arbitrary worker completion order, panic containment,
//! abort, memory-budget enforcement, seam chain-stitching, and thread-count
//! invariance.

use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Duration;

use reed_core::ChunkId;
use reed_engine::abort::abort;
use reed_engine::budget::{MemBudget, ThreadBudget};
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_io::bgzf::deflate::Deflater;
use reed_io::chunk::{ChunkBytes, ChunkSource};

/// Collects the exact byte sequence of every record line, in arrival order.
/// The sink re-folds each line exactly as a serial scan would, so any line
/// dropped, duplicated, or reordered changes the digest.
#[derive(Default)]
struct CountJob;

impl ChunkJob for CountJob {
    type Partial = Vec<Vec<u8>>;
    type Scratch = ();

    fn run_chunk(
        &self,
        _chunk: &ChunkBytes,
        body: &[u8],
        _scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        Ok(lines_of(body))
    }

    fn run_seam(
        &self,
        _id: reed_core::ChunkId,
        line: &[u8],
        _scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        Ok(lines_of(line))
    }
}

fn lines_of(bytes: &[u8]) -> Vec<Vec<u8>> {
    bytes
        .split(|&b| b == b'\n')
        .filter(|l| !l.is_empty())
        .map(Vec::from)
        .collect()
}

/// Folds lines exactly as `fold_lines` below, in the order the engine
/// delivers partials.
#[derive(Default)]
struct VecSink {
    fold: u64,
    n_lines: u64,
}



fn fold_line(mut fold: u64, line: &[u8]) -> u64 {
    for &b in line {
        fold = (fold ^ u64::from(b)).wrapping_mul(0x100000001b3);
    }
    fold.wrapping_add(1)
}

fn fold_lines(bytes: &[u8]) -> u64 {
    // Same zero seed as VecSink::fold: both sides fold identically.
    let mut fold = 0_u64;
    for line in bytes.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        fold = fold_line(fold, line);
    }
    fold
}

impl OrderedSink<Vec<Vec<u8>>> for VecSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut Vec<Vec<u8>>) -> Result<(), EngineError> {
        for line in partial.drain(..) {
            self.fold = fold_line(self.fold, &line);
            self.n_lines += 1;
        }
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        Ok(self.fold)
    }
}

/// Sink for jobs with `Partial = ()`.
#[derive(Default)]
struct UnitSink;

impl OrderedSink<()> for UnitSink {
    fn accept(&mut self, _id: ChunkId, _partial: &mut ()) -> Result<(), EngineError> {
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        Ok(0)
    }
}

/// In-memory chunk source.
struct VecSource {
    chunks: std::vec::IntoIter<ChunkBytes>,
}

impl ChunkSource for VecSource {
    fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, reed_io::IoError> {
        Ok(self.chunks.next())
    }
}

/// Compresses `text` into chunks of `chunk_bytes` uncompressed text each,
/// cutting block boundaries mid-line to force seams.
fn make_chunks(text: &[u8], chunk_bytes: usize) -> Vec<ChunkBytes> {
    let mut deflater = Deflater::new(6).expect("invariant: level 6 is valid");
    let mut chunks = Vec::new();
    for (id, piece) in text.chunks(chunk_bytes).enumerate() {
        let mut compressed = Vec::new();
        deflater.deflate_block(piece, &mut compressed).expect("invariant: test block fits");
        compressed.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
        chunks.push(ChunkBytes {
            id: ChunkId(id as u64),
            compressed: compressed.into_boxed_slice().into(),
            first_in_block_offset: 0,
            plain: false,
        });
    }
    chunks
}

fn fixture_text(n_records: usize) -> Vec<u8> {
    let mut text = Vec::new();
    for i in 0..n_records {
        writeln!(
            &mut text,
            "chr1\t{}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1|0\t./.",
            100 + i * 3
        )
        .expect("invariant: write to Vec");
    }
    text
}

fn budget(threads: usize) -> ThreadBudget {
    ThreadBudget {
        threads,
        mem: MemBudget::from_bytes(u64::MAX),
    }
}

#[test]
fn output_is_invariant_across_thread_counts() {
    let text = fixture_text(2_000);
    let chunks = make_chunks(&text, 700); // seams + multiple chunks
    let expected = {
        let mut source = VecSource { chunks: chunks.clone().into_iter() };
        let abort = abort();
        run(&mut source, &CountJob, VecSink::default(), &budget(1), &abort).expect("invariant: run")
    };
    for threads in [2_usize, 3, 7, 16] {
        let mut source = VecSource { chunks: chunks.clone().into_iter() };
        let abort = abort();
        let report = run(&mut source, &CountJob, VecSink::default(), &budget(threads), &abort)
            .expect("invariant: parallel run matches sequential");
        assert_eq!(report.sink_report, expected.sink_report, "threads={threads}");
        assert_eq!(report.n_chunks, expected.n_chunks, "threads={threads}");
    }
    // Sanity: the fixture must have produced records and multiple chunks.
    assert!(expected.n_chunks > 2, "test must span multiple chunks");
}

/// Regression (01-review D1): a slow chunk must not let later chunk ids
/// overflow the reorder window. The reducer used to pull finished results
/// into the ring while chunk 0 was still running, freeing result-channel
/// slots until ids exceeded the ring capacity and the `insert` assert
/// panicked. With real tickets the reader blocks at the token budget, so
/// the run completes instead.
#[test]
fn slow_first_chunk_does_not_overflow_the_reorder_window() {
    struct DelayFirstJob;
    impl ChunkJob for DelayFirstJob {
        type Partial = Vec<Vec<u8>>;
        type Scratch = ();
        fn run_chunk(
            &self,
            chunk: &ChunkBytes,
            body: &[u8],
            _scratch: &mut Self::Scratch,
        ) -> Result<Self::Partial, EngineError> {
            if chunk.id.get() == 0 {
                std::thread::sleep(Duration::from_millis(500));
            }
            Ok(lines_of(body))
        }
        fn run_seam(
            &self,
            _id: ChunkId,
            line: &[u8],
            _scratch: &mut Self::Scratch,
        ) -> Result<Self::Partial, EngineError> {
            Ok(lines_of(line))
        }
    }
    // Far more chunks than the old window (tickets + workers), so the
    // fast chunks would overflow the ring while chunk 0 sleeps.
    let text = fixture_text(8_000);
    let chunks = make_chunks(&text, 500);
    assert!(chunks.len() >= 60, "test needs many chunks, got {}", chunks.len());
    let n_chunks = chunks.len() as u64;
    let mut source = VecSource { chunks: chunks.into_iter() };
    let abort = abort();
    let report = run(&mut source, &DelayFirstJob, VecSink::default(), &budget(8), &abort)
        .expect("invariant: slow chunk must not overflow the reorder window");
    assert_eq!(report.n_chunks, n_chunks);
}

#[test]
fn seam_records_are_counted_exactly_once() {
    // 50 records compressed as one block each is too coarse; 100-byte
    // chunks guarantee mid-line cuts and thus seams.
    let text = fixture_text(50);
    let all_chunks = make_chunks(&text, 100);
    let n_chunks = all_chunks.len() as u64;
    assert!(n_chunks > 10, "test needs many seam-rich chunks, got {n_chunks}");
    let mut seq_source = VecSource { chunks: all_chunks.clone().into_iter() };
    let abort2 = abort();
    let seq = run(&mut seq_source, &CountJob, VecSink::default(), &budget(1), &abort2)
        .expect("invariant: sequential reference");
    let mut source = VecSource { chunks: all_chunks.into_iter() };
    let abort = abort();
    let report = run(&mut source, &CountJob, VecSink::default(), &budget(4), &abort)
        .expect("invariant: run with seams");
    assert_eq!(report.sink_report, seq.sink_report, "parallel must equal sequential");
    assert_eq!(seq.sink_report, fold_lines(&text), "sequential must equal a serial scan");
    assert_eq!(report.n_chunks, n_chunks);
}

#[test]
fn worker_panic_becomes_worker_panicked() {
    struct PanickyJob;
    impl ChunkJob for PanickyJob {
        type Partial = ();
        type Scratch = ();
        fn run_chunk(&self, chunk: &ChunkBytes, _body: &[u8], _: &mut ()) -> Result<(), EngineError> {
            if chunk.id.get() == 1 {
                panic!("boom on chunk 1");
            }
            Ok(())
        }
        fn run_seam(&self, _id: reed_core::ChunkId, _line: &[u8], _: &mut ()) -> Result<(), EngineError> {
            Ok(())
        }
    }
    let text = fixture_text(200);
    let chunks = make_chunks(&text, 500);
    let mut source = VecSource { chunks: chunks.into_iter() };
    let abort = abort();
    let err = run(&mut source, &PanickyJob, UnitSink, &budget(3), &abort).unwrap_err();
    match err {
        EngineError::WorkerPanicked { chunk, message } => {
            assert_eq!(chunk.get(), 1);
            assert_eq!(message, "boom on chunk 1");
        }
        other => panic!("expected WorkerPanicked, got {other:?}"),
    }
}

#[test]
fn external_abort_stops_the_run() {
    struct SlowJob;
    impl ChunkJob for SlowJob {
        type Partial = ();
        type Scratch = ();
        fn run_chunk(&self, _chunk: &ChunkBytes, _body: &[u8], _: &mut ()) -> Result<(), EngineError> {
            std::thread::sleep(Duration::from_millis(5));
            Ok(())
        }
        fn run_seam(&self, _id: reed_core::ChunkId, _line: &[u8], _: &mut ()) -> Result<(), EngineError> {
            Ok(())
        }
    }
    let text = fixture_text(5_000);
    let chunks = make_chunks(&text, 500);
    let mut source = VecSource { chunks: chunks.into_iter() };
    let abort = abort();
    // Abort shortly after the run starts; the pipeline must observe it at
    // its next blocking point and return Aborted.
    {
        let handle = abort.handle.clone();
        std::thread::spawn(move || {
            std::thread::sleep(Duration::from_millis(50));
            handle.trigger();
        });
    }
    let result = run(&mut source, &SlowJob, UnitSink, &budget(4), &abort);
    assert!(matches!(result, Err(EngineError::Aborted)));
}

#[test]
fn memory_budget_stops_the_run() {
    // Chunks compress to a few dozen bytes; a 64 B budget trips as soon as
    // more than one chunk is in flight (tickets = 2 × workers).
    let text = fixture_text(2_000);
    let chunks = make_chunks(&text, 2_000);
    assert!(chunks.len() > 3, "test needs several chunks");
    let mut source = VecSource { chunks: chunks.into_iter() };
    let abort = abort();
    let b = ThreadBudget {
        threads: 4,
        mem: MemBudget::from_bytes(64),
    };
    let result = run(&mut source, &CountJob, VecSink::default(), &b, &abort);
    assert!(
        matches!(result, Err(EngineError::MemBudgetExceeded { .. })),
        "expected MemBudgetExceeded, got {result:?}"
    );
}

#[test]
fn source_errors_propagate() {
    struct FailingSource;
    impl ChunkSource for FailingSource {
        fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, reed_io::IoError> {
            Err(reed_io::IoError::Io("disk exploded".into()))
        }
    }
    let mut source = FailingSource;
    let abort = abort();
    let err = run(&mut source, &CountJob, VecSink::default(), &budget(3), &abort).unwrap_err();
    assert!(matches!(err, EngineError::Io(ref m) if m.contains("disk exploded")));
}

#[test]
fn partials_arrive_in_chunk_order_at_the_sink() {
    // A job that records its own chunk id; the sink asserts ids increase.
    struct IdJob {
        seen: Arc<AtomicU64>,
    }
    impl ChunkJob for IdJob {
        type Partial = u64;
        type Scratch = ();
        fn run_chunk(&self, chunk: &ChunkBytes, _body: &[u8], _: &mut ()) -> Result<u64, EngineError> {
            Ok(chunk.id.get())
        }
        fn run_seam(&self, _id: reed_core::ChunkId, _line: &[u8], _: &mut ()) -> Result<u64, EngineError> {
            Ok(u64::MAX) // seam marker; excluded from the order assertion below
        }
    }
    struct IdSink {
        last_id: u64,
        n_seams: u64,
    }
    impl OrderedSink<u64> for IdSink {
        fn accept(&mut self, _id: ChunkId, partial: &mut u64) -> Result<(), EngineError> {
            if *partial == u64::MAX {
                self.n_seams += 1;
                return Ok(());
            }
            // Chunk ids arrive strictly ascending (seam markers excluded).
            assert!(
                *partial > self.last_id || self.last_id == 0,
                "chunk ids must arrive in order: {} after {}",
                partial,
                self.last_id
            );
            self.last_id = *partial;
            Ok(())
        }
        fn finish(self) -> Result<u64, EngineError> {
            Ok(self.last_id)
        }
    }
    // Sleep per chunk id parity forces out-of-order worker completion.
    struct SleepJob;
    impl ChunkJob for SleepJob {
        type Partial = u64;
        type Scratch = ();
        fn run_chunk(&self, chunk: &ChunkBytes, _body: &[u8], _: &mut ()) -> Result<u64, EngineError> {
            if chunk.id.get().is_multiple_of(2) {
                std::thread::sleep(Duration::from_millis(2));
            }
            Ok(chunk.id.get())
        }
        fn run_seam(&self, _id: reed_core::ChunkId, _line: &[u8], _: &mut ()) -> Result<u64, EngineError> {
            Ok(u64::MAX)
        }
    }
    let text = fixture_text(3_000);
    let chunks = make_chunks(&text, 600);
    let n_chunks = chunks.len();
    let mut source = VecSource { chunks: chunks.into_iter() };
    let abort = abort();
    let report = run(
        &mut source,
        &SleepJob,
        IdSink { last_id: 0, n_seams: 0 },
        &budget(6),
        &abort,
    )
    .expect("invariant: ordered delivery");
    assert_eq!(report.sink_report, (n_chunks - 1) as u64);
}


/// R3 gate (03-architecture.md §4.2): the memory budget is backpressure —
/// the run *completes* with a small pool, and the observed in-flight peak
/// stays within the budget.
#[test]
fn memory_budget_bounds_in_flight_bytes() {
    // Chunks of ~1 KB compressed (footprint ≈ 8 KB each): a 64 KB pool
    // admits at most a handful in flight, and the run still finishes.
    let text = fixture_text(20_000);
    let chunks = make_chunks(&text, 1_000);
    assert!(chunks.len() > 20, "test needs many chunks");
    let budget_bytes: u64 = 64 * 1024;
    let b = ThreadBudget {
        threads: 4,
        mem: MemBudget::from_bytes(budget_bytes),
    };
    let mut source = VecSource { chunks: chunks.into_iter() };
    let abort = abort();
    let report = run(&mut source, &CountJob, VecSink::default(), &b, &abort)
        .expect("invariant: backpressure completes the run");
    assert!(
        report.metrics.in_flight_peak_bytes <= budget_bytes,
        "in-flight peak {} must respect the budget {budget_bytes}",
        report.metrics.in_flight_peak_bytes
    );
    assert!(report.metrics.in_flight_peak_bytes > 0, "peak must be observed");
}
