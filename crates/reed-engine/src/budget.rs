//! Resource budgets (`docs/spec/architecture.md` §7.3, §7.7).

use std::path::Path;

/// In-flight chunks per worker (§7.4): the ticket window capacity is
/// `TICKETS_PER_WORKER × workers`.
pub const TICKETS_PER_WORKER: usize = 2;

/// The largest `--threads` value reed accepts (doc_v26 V26-03).
///
/// The old parser rejected only `0`, so an arbitrary `usize` reached the
/// derived capacities — `threads + 4` abort tokens, `TICKETS_PER_WORKER ×
/// workers` channel slots and reorder cells, `threads` thread spawns. A
/// four-digit cap keeps every derivation inside a sane resource envelope
/// (tickets ≤ 2048 slots), is machine-independent so CI and production
/// agree, and is far above any supported deployment (the reference grid is
/// 348 samples on 56 cores).
pub const MAX_THREADS: usize = 1024;

/// A validated `--threads` value (doc_v25 V25-06).
///
/// `0` used to mean "all cores" on the distance path (rayon's convention)
/// and "sequential" on every `ThreadBudget` path — a silent resource split.
/// The CLI now rejects it once, at parse time, so every command agrees.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Threads(usize);

impl Threads {
    /// The thread count (≥ 1).
    #[must_use]
    pub const fn get(self) -> usize {
        self.0
    }
}

impl std::fmt::Display for Threads {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::str::FromStr for Threads {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let value: usize = text
            .trim()
            .parse()
            .map_err(|_| format!("--threads: {text:?} is not a thread count"))?;
        if value == 0 {
            return Err("--threads: must be >= 1 (0 does not mean \"all cores\")".to_string());
        }
        if value > MAX_THREADS {
            return Err(format!(
                "--threads: {value} exceeds the supported maximum of {MAX_THREADS}"
            ));
        }
        Ok(Self(value))
    }
}

/// A validated `--memory` value in MiB (doc_v25 V25-06).
///
/// The old `saturating_mul(1 << 20)` silently turned an overflowing value
/// into an unlimited budget; the newtype rejects `0` and anything that does
/// not fit the byte budget at parse time.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MemMiB(u32);

impl MemMiB {
    /// The value in MiB (≥ 1).
    #[must_use]
    pub const fn get(self) -> u32 {
        self.0
    }

    /// The value in bytes (exact: `u32` MiB always fits `u64`).
    #[must_use]
    pub const fn bytes(self) -> u64 {
        (self.0 as u64) << 20
    }
}

/// Clap `value_parser` for `--threads` (the `FromStr` grammar, one copy).
///
/// # Errors
///
/// The `FromStr` message.
pub fn parse_threads(text: &str) -> Result<Threads, String> {
    text.parse()
}

/// Clap `value_parser` for `--memory`.
///
/// # Errors
///
/// The `FromStr` message.
pub fn parse_mem_mib(text: &str) -> Result<MemMiB, String> {
    text.parse()
}

impl std::fmt::Display for MemMiB {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::str::FromStr for MemMiB {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let value: u32 = text
            .trim()
            .parse()
            .map_err(|_| format!("--memory: {text:?} is not a MiB count"))?;
        if value == 0 {
            return Err("--memory: must be >= 1 MiB".to_string());
        }
        Ok(Self(value))
    }
}

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

    /// The CLI's single budget constructor (doc_v25 V25-06): `None` threads
    /// means all available cores (1 when the platform refuses to report —
    /// the same fallback as [`Self::all`], never a panic), `None` memory the
    /// default half-memory budget. Both values are already validated by the
    /// CLI's value parsers ([`Threads`], [`MemMiB`]).
    #[must_use]
    pub fn from_cli(threads: Option<Threads>, mem_mib: Option<MemMiB>) -> Self {
        Self {
            threads: threads.map_or_else(
                || {
                    std::thread::available_parallelism()
                        .map_or(1, std::num::NonZeroUsize::get)
                },
                Threads::get,
            ),
            mem: mem_mib.map_or_else(MemBudget::default_half_memory, |m| {
                MemBudget::from_bytes(m.bytes())
            }),
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
    ///
    /// `--threads` is capped at [`MAX_THREADS`], so the product cannot
    /// overflow; the saturating form keeps the function total for budgets
    /// built directly in library code (doc_v26 V26-03).
    #[must_use]
    pub fn tickets(&self) -> usize {
        TICKETS_PER_WORKER.saturating_mul(self.workers())
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
    ///
    /// The notify holds the same mutex the waiter checks the flag under, so
    /// a notify that lands between the waiter's check and its `wait` cannot
    /// be lost (doc_v15 D15-32).
    pub fn notify_abort(&self) {
        let _guard = self
            .available
            .lock()
            .expect("invariant: permits not poisoned");
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

    /// doc_v26 V26-03: the parser rejects 0 and anything above the cap, and
    /// accepts the cap itself.
    #[test]
    fn the_thread_parser_has_a_documented_upper_bound() {
        assert_eq!(MAX_THREADS, 1024, "the CLI contract test pins this value");
        assert_eq!("1".parse::<Threads>().expect("invariant: 1").get(), 1);
        assert_eq!(
            MAX_THREADS.to_string().parse::<Threads>().expect("invariant: cap").get(),
            MAX_THREADS
        );
        for bad in ["0", "0 ", &(MAX_THREADS + 1).to_string(), "18446744073709551615"] {
            let err = bad.parse::<Threads>().expect_err("must be refused");
            assert!(
                err.contains(">= 1") || err.contains("exceeds the supported maximum"),
                "{bad:?}: {err}"
            );
        }
        // A budget built directly in library code cannot overflow the derived
        // capacities (the CLI would refuse it, but the library is public).
        let huge = threads(usize::MAX);
        assert_eq!(huge.tickets(), usize::MAX);
        assert_eq!(huge.workers(), usize::MAX - 2);
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
