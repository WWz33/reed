//! Stage metrics (`docs/spec/architecture.md` §7.10; `docs/spec/v6_y5/README.md`
//! §2 S5): aggregate counters and timings, enough for the §17.3 baseline log.

use std::time::Duration;

/// Aggregate metrics of one engine run.
#[derive(Clone, Copy, Debug, Default)]
pub struct StageMetrics {
    /// Chunks read from the source.
    pub n_chunks_read: u64,
    /// Sum of busy time (inflate + job) across all workers, plus the seam
    /// work the reducer performs between chunks (doc_v9 D9-6), so it always
    /// covers `inflate_time + job_time`.
    pub worker_busy: Duration,
    /// Peak in-flight *estimated footprint*: the largest sum of
    /// `MemBudget::estimate` over the chunks concurrently in flight, i.e.
    /// compressed bytes × [`crate::budget::MEM_FACTOR`] (doc_v9 D9-5). Not
    /// raw compressed bytes; the sequential path reports the largest single
    /// chunk's footprint.
    pub in_flight_peak_bytes: u64,
    /// In-flight bytes still held when the run ended (a panic or abort must
    /// leave this at zero; the Y5 gate asserts it).
    pub in_flight_final_bytes: u64,
    /// Cumulative time the reader spent in `ChunkSource::next_chunk`
    /// (transport + I/O). doc_v4 03 A4-5: without it, "reader-bound" is an
    /// inference.
    pub reader_time: Duration,
    /// Cumulative time in `OrderedSink::accept`, including the seam record
    /// the reducer delivers between chunks (doc_v9 D9-6).
    pub reduce_time: Duration,
    /// Cumulative time workers spent waiting on a job-internal chain
    /// ([`crate::pipeline::ChunkJob::wait_time`]; the BCF frame chain).
    pub chain_wait: Duration,
    /// Cumulative time the reader spent blocked on the ticket window, the
    /// memory permits, or the chunk send (A6-4).
    pub reader_wait: Duration,
    /// Cumulative time the reducer spent blocked on the result channel.
    pub reducer_wait: Duration,
    /// Cumulative `inflate_chunk` time across all workers (A6-4).
    pub inflate_time: Duration,
    /// Cumulative `ChunkJob::run_chunk`/`run_seam` time. Seams run on the
    /// reducer, so this spans workers **and** the reducer (doc_v9 D9-6);
    /// output deflate is a subset ([`StageMetrics::deflate_time`]).
    pub job_time: Duration,
    /// Cumulative output-compression time, reported by the sink (A6-4).
    pub deflate_time: Duration,
}

impl StageMetrics {
    /// Adds another run's metrics (the multi-pass drivers sum pass0 and
    /// pass1; the peaks are maxima, not sums).
    pub fn merge(&mut self, other: &Self) {
        self.n_chunks_read += other.n_chunks_read;
        self.worker_busy += other.worker_busy;
        self.reader_time += other.reader_time;
        self.reduce_time += other.reduce_time;
        self.chain_wait += other.chain_wait;
        self.reader_wait += other.reader_wait;
        self.reducer_wait += other.reducer_wait;
        self.inflate_time += other.inflate_time;
        self.job_time += other.job_time;
        self.deflate_time += other.deflate_time;
        self.in_flight_peak_bytes = self.in_flight_peak_bytes.max(other.in_flight_peak_bytes);
        self.in_flight_final_bytes = self.in_flight_final_bytes.max(other.in_flight_final_bytes);
    }

    /// The canonical `--timing` key/value block, one `key<TAB>value` line per
    /// metric (doc_v25 V25-08b).
    ///
    /// Every command that prints engine timings renders through this method:
    /// filter, popg, stats, convert and index used to spell the same metrics
    /// with four different key names (`engine.`-prefixed reader/reduce times
    /// vs the bare names, `chunks` vs `chunks_read`, missing
    /// `worker_busy_ms`/`chain_wait_ms`),
    /// so a performance-regression parser needed one table per subcommand.
    #[must_use]
    pub fn render_tsv(&self) -> String {
        use std::fmt::Write as _;
        let mut text = String::with_capacity(256);
        let ms = |d: Duration| d.as_secs_f64() * 1e3;
        let _ = writeln!(text, "chunks_read\t{}", self.n_chunks_read);
        let _ = writeln!(text, "worker_busy_ms\t{:.3}", ms(self.worker_busy));
        let _ = writeln!(text, "reader_ms\t{:.3}", ms(self.reader_time));
        let _ = writeln!(text, "reduce_ms\t{:.3}", ms(self.reduce_time));
        let _ = writeln!(text, "chain_wait_ms\t{:.3}", ms(self.chain_wait));
        let _ = writeln!(text, "reader_wait_ms\t{:.3}", ms(self.reader_wait));
        let _ = writeln!(text, "reduce_wait_ms\t{:.3}", ms(self.reducer_wait));
        let _ = writeln!(text, "inflate_ms\t{:.3}", ms(self.inflate_time));
        let _ = writeln!(text, "job_ms\t{:.3}", ms(self.job_time));
        let _ = writeln!(text, "deflate_ms\t{:.3}", ms(self.deflate_time));
        let _ = writeln!(text, "in_flight_peak_bytes\t{}", self.in_flight_peak_bytes);
        text
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn render_tsv_names_every_engine_metric_once() {
        let metrics = StageMetrics {
            n_chunks_read: 7,
            worker_busy: Duration::from_millis(2),
            reader_time: Duration::from_millis(3),
            reduce_time: Duration::from_millis(4),
            chain_wait: Duration::from_millis(5),
            reader_wait: Duration::from_millis(6),
            reducer_wait: Duration::from_millis(7),
            inflate_time: Duration::from_millis(8),
            job_time: Duration::from_millis(9),
            deflate_time: Duration::from_millis(10),
            in_flight_peak_bytes: 11,
            in_flight_final_bytes: 0,
        };
        let text = metrics.render_tsv();
        let keys: Vec<&str> = text
            .lines()
            .map(|line| line.split('\t').next().unwrap_or_default())
            .collect();
        assert_eq!(
            keys,
            [
                "chunks_read",
                "worker_busy_ms",
                "reader_ms",
                "reduce_ms",
                "chain_wait_ms",
                "reader_wait_ms",
                "reduce_wait_ms",
                "inflate_ms",
                "job_ms",
                "deflate_ms",
                "in_flight_peak_bytes",
            ]
        );
        assert!(text.contains("chunks_read\t7"));
        assert!(text.contains("job_ms\t9.000"));
        assert!(text.contains("in_flight_peak_bytes\t11"));
    }
}
