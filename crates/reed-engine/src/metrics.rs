//! Stage metrics (`docs/architecture-spec.md` §7.10): minimal M3 form —
//! aggregate counters and timings, enough for the §17.3 baseline log.

use std::time::Duration;

/// Aggregate metrics of one engine run.
#[derive(Clone, Copy, Debug, Default)]
pub struct StageMetrics {
    /// Chunks read from the source.
    pub n_chunks_read: u64,
    /// Sum of worker busy time (inflate + job) across all workers.
    pub worker_busy: Duration,
    /// Peak in-flight compressed bytes (ticket-window occupancy).
    pub in_flight_peak_bytes: u64,
    /// Cumulative time the reader spent in `ChunkSource::next_chunk`
    /// (transport + I/O). doc_v4 03 A4-5: without it, "reader-bound" is an
    /// inference.
    pub reader_time: Duration,
    /// Cumulative time the reducer spent in `OrderedSink::accept`.
    pub reduce_time: Duration,
    /// Cumulative time workers spent waiting on a job-internal chain
    /// ([`crate::pipeline::ChunkJob::wait_time`]; the BCF frame chain).
    pub chain_wait: Duration,
}

impl StageMetrics {
    /// Adds another run's metrics (the multi-pass drivers sum pass0 and
    /// pass1; the in-flight peak is the maximum, not a sum).
    pub fn merge(&mut self, other: &Self) {
        self.n_chunks_read += other.n_chunks_read;
        self.worker_busy += other.worker_busy;
        self.reader_time += other.reader_time;
        self.reduce_time += other.reduce_time;
        self.chain_wait += other.chain_wait;
        self.in_flight_peak_bytes = self.in_flight_peak_bytes.max(other.in_flight_peak_bytes);
    }
}
