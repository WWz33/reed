//! Resource budgets (`docs/spec/architecture.md` §7.3, §7.7).

use std::path::Path;

/// In-flight chunks per worker (§7.4): the ticket window capacity is
/// `TICKETS_PER_WORKER × workers`.
pub const TICKETS_PER_WORKER: usize = 2;

/// The single compute-resource knob: total threads reed may use.
///
/// Layout: 1 reader + 1 reducer/writer + `workers()` compute workers.
/// `threads <= 2` is the intended sequential-mode selector (§7.3), whose
/// output must be byte-identical to the parallel mode; the sequential path
/// itself is not implemented yet (`is_sequential` has no caller — F6-6,
/// decided in the Y5 engine work, docs/rounds/v6/06-roadmap.md).
#[derive(Clone, Copy, Debug)]
pub struct ThreadBudget {
    /// Total threads (default: all available cores).
    pub threads: usize,
    /// In-memory upper bound for in-flight chunk data (§7.7).
    pub mem: MemBudget,
}

impl ThreadBudget {
    /// Budget using all available cores and the default memory budget.
    #[must_use]
    pub fn all() -> Self {
        Self {
            threads: std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get),
            mem: MemBudget::default_half_memory(),
        }
    }

    /// Compute-worker count (the rest of the budget goes to reader/reducer).
    #[must_use]
    pub fn workers(&self) -> usize {
        self.threads.saturating_sub(2).max(1)
    }

    /// Ticket-window capacity: the hard upper bound on in-flight chunks.
    #[must_use]
    pub fn tickets(&self) -> usize {
        TICKETS_PER_WORKER * self.workers()
    }

    /// Whether the run takes the single-threaded sequential path.
    #[must_use]
    pub const fn is_sequential(&self) -> bool {
        self.threads <= 2
    }
}

impl Default for ThreadBudget {
    fn default() -> Self {
        Self::all()
    }
}

/// In-memory upper bound for in-flight chunk data (§7.7, §2 principle 11).
///
/// The budget is enforced as *backpressure*: the reader blocks until the
/// reducer returns permits for delivered chunks. It fails the run only when
/// a single chunk alone exceeds the whole budget.
#[derive(Clone, Copy, Debug)]
pub struct MemBudget {
    /// Maximum in-flight (estimated) bytes.
    pub bytes: u64,
}

/// Estimated in-memory footprint per compressed input byte (compressed ≈
/// 5× smaller than inflated, plus output blocks and partial results).
pub const MEM_FACTOR: u64 = 8;

impl MemBudget {
    /// An explicit byte budget.
    #[must_use]
    pub const fn from_bytes(bytes: u64) -> Self {
        Self { bytes }
    }

    /// Default: 50% of physical memory (Linux; other platforms: unlimited).
    #[must_use]
    pub fn default_half_memory() -> Self {
        Self {
            bytes: physical_memory_bytes().map_or(u64::MAX, |total| total / 2),
        }
    }

    /// Estimated in-memory footprint of one chunk.
    #[must_use]
    pub const fn estimate(&self, compressed_len: u64) -> u64 {
        compressed_len * MEM_FACTOR
    }
}

/// Byte permits for in-flight chunks (03-architecture.md §4.2): the reader
/// acquires a chunk's estimated footprint before handing it to a worker and
/// the reducer releases it after in-order delivery. Over budget → the
/// reader *blocks*, not fails.
#[derive(Debug)]
pub struct BytePermits {
    available: std::sync::Mutex<u64>,
    cv: std::sync::Condvar,
}

impl BytePermits {
    /// A pool of `bytes` permits.
    #[must_use]
    pub fn new(bytes: u64) -> Self {
        Self {
            available: std::sync::Mutex::new(bytes),
            cv: std::sync::Condvar::new(),
        }
    }

    /// Blocks until `bytes` permits are available (or abort fires); returns
    /// `false` on abort.
    pub fn acquire(&self, bytes: u64, abort: &crate::abort::AbortSignal) -> bool {
        let mut guard = self
            .available
            .lock()
            .expect("invariant: permits not poisoned");
        loop {
            if abort.is_aborted() {
                return false;
            }
            if *guard >= bytes {
                *guard -= bytes;
                return true;
            }
            // A short poll keeps abort latency low without a watcher thread.
            let (g, _) = self
                .cv
                .wait_timeout(guard, std::time::Duration::from_millis(10))
                .expect("invariant: permits not poisoned");
            guard = g;
        }
    }

    /// Returns `bytes` permits to the pool.
    pub fn release(&self, bytes: u64) {
        let mut guard = self
            .available
            .lock()
            .expect("invariant: permits not poisoned");
        *guard = guard.saturating_add(bytes);
        drop(guard);
        self.cv.notify_all();
    }
}

/// Physical memory in bytes (Linux `/proc/meminfo`), `None` if unknown.
#[must_use]
fn physical_memory_bytes() -> Option<u64> {
    let info = std::fs::read_to_string(Path::new("/proc/meminfo")).ok()?;
    for line in info.lines() {
        if let Some(rest) = line.strip_prefix("MemTotal:") {
            let kib: u64 = rest.trim().trim_end_matches(" kB").trim().parse().ok()?;
            return Some(kib * 1024);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn threads(n: usize) -> ThreadBudget {
        ThreadBudget {
            threads: n,
            mem: MemBudget::from_bytes(u64::MAX),
        }
    }

    #[test]
    fn thread_layout_follows_the_spec() {
        assert!(threads(2).is_sequential());
        assert!(threads(1).is_sequential());
        assert_eq!(threads(2).workers(), 1);
        assert_eq!(threads(56).workers(), 54);
        assert_eq!(threads(56).tickets(), 108);
    }

    #[test]
    fn mem_budget_estimates_footprints() {
        let budget = MemBudget::from_bytes(1000);
        assert_eq!(budget.estimate(10), 80);
        // Backpressure, not error: over-budget chunks block the reader; only
        // a single chunk larger than the whole budget fails the run (03 §4.2).
        let permits = BytePermits::new(80);
        assert!(permits.acquire(80, &crate::abort::abort().signal));
        permits.release(80);
        assert!(*permits.available.lock().unwrap() >= 80);
    }
}
