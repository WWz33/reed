//! Resource budgets (`docs/spec/architecture.md` §7.3, §7.7).

use std::path::Path;

/// In-flight chunks per worker (§7.4): the ticket window capacity is
/// `TICKETS_PER_WORKER × workers`.
pub const TICKETS_PER_WORKER: usize = 2;

/// The single compute-resource knob: total threads reed may use.
///
/// Layout: 1 reader + 1 reducer/writer + `workers()` compute workers.
/// `threads <= 2` selects the sequential path (§7.3), whose output is
/// byte-identical to the parallel mode (locked by tests).
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

    /// Thread count for the auxiliary rayon stages and per-contig caps
    /// (`docs/rounds/v8` D8-1): the mid-stage decisions, the LD prune and
    /// the LD decay take their parallelism from this one place, so no
    /// stage can outgrow the budget and land on the rayon global pool.
    /// Sequential budgets stay at 1: the auxiliary stages must run serial
    /// when the pipeline does.
    #[must_use]
    pub fn auxiliary(&self) -> usize {
        if self.is_sequential() {
            1
        } else {
            self.workers()
        }
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

/// Runs `f` on the auxiliary stage's parallelism: a throwaway rayon pool of
/// exactly `budget.auxiliary()` threads (doc_v8 D8-1) — never the rayon
/// global pool, whose thread count would escape the budget. Sequential
/// budgets skip rayon entirely: `f` receives `parallel == false` and must
/// run plain iterators (a par_iter inside a serial f would still lazily
/// spawn the global pool).
///
/// # Errors
///
/// `String` when the pool cannot be built (thread creation refused).
///
/// # Panics
///
/// Propagates `f`'s panics.
pub fn run_auxiliary<T, F>(budget: &ThreadBudget, f: F) -> Result<T, String>
where
    F: FnOnce(bool) -> T + Send,
    T: Send,
{
    let n = budget.auxiliary();
    if n <= 1 {
        // Serial: f(false) must run plain iterators — a par_iter in the
        // closure would lazily spawn the rayon *global* pool (all cores),
        // the exact escape doc_v8 D8-1 closes.
        return Ok(f(false));
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .map_err(|e| e.to_string())?;
    Ok(pool.install(|| f(true)))
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
    ///
    /// The wait is a plain condvar wait, not a poll: the abort handle
    /// notifies this pool on `trigger` (see
    /// [`crate::abort::AbortHandle::register`]), so abort latency is the
    /// wake-up latency (§7.5 "中止即时生效").
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
            guard = self
                .cv
                .wait(guard)
                .expect("invariant: permits not poisoned");
        }
    }

    /// Wakes every waiter so it re-checks the abort flag. Called by
    /// [`crate::abort::AbortHandle::trigger`].
    pub fn notify_abort(&self) {
        self.cv.notify_all();
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
    fn auxiliary_stages_share_the_worker_budget() {
        // doc_v8 D8-1a: the rayon stages and per-contig caps take their
        // thread count from one place, so `--threads 8` bounds them at 6,
        // not at the rayon global pool's core count.
        assert_eq!(threads(8).auxiliary(), 6);
        assert_eq!(threads(1).auxiliary(), 1);
        assert_eq!(threads(2).auxiliary(), 1);
        assert_eq!(threads(56).auxiliary(), 54);
    }

    #[test]
    fn auxiliary_pool_stays_inside_the_budget() {
        // doc_v8 D8-1: the throwaway pool must have exactly auxiliary()
        // threads; the sequential budget runs on the calling thread, where
        // rayon reports the ambient global pool (whose size is irrelevant
        // because no par_iter executes on that path).
        let seen = |b: &ThreadBudget| {
            run_auxiliary(b, |_| rayon::current_num_threads()).expect("invariant: pool builds")
        };
        assert_eq!(seen(&threads(8)), 6);
        assert_eq!(seen(&threads(56)), 54);
        assert_eq!(seen(&threads(1)), rayon::current_num_threads());
        assert_eq!(seen(&threads(2)), rayon::current_num_threads());
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
