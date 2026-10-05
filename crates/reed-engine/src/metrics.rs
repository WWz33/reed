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
}
