//! Length-prefix chains (doc_v3 02-bcf.md §2.2): a BCF record has no
//! synchronisation marker, so a chunk cannot find its first record start on
//! its own. Only one small piece has to travel along the chunk sequence —
//! the bytes of the record that started in a previous chunk — and walking a
//! chunk's length prefixes is O(records) with no extra inflation, so the
//! hand-off is microseconds.
//!
//! Cell `k` holds the carried bytes of the record that started before chunk
//! `k` (empty = chunk `k` begins at a record boundary). The reader publishes
//! cell 0 before the run; worker `k` waits for cell `k`, publishes cell
//! `k + 1` as soon as its own walk is done, and only then decodes — a
//! predecessor never waits for a successor, so the chain cannot deadlock
//! (chunks are dispatched in order and chunk 0 never waits).

// Casts in this file are bounded by construction: Chunk ids are counted in memory and fit `usize` on every supported target.
#![allow(clippy::cast_possible_truncation)]

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant};

use crate::abort::AbortSignal;
use crate::error::EngineError;

/// The production watchdog: a chain that publishes nothing for this long is
/// treated as stalled (doc_v4 01 F4-14).
const DEFAULT_STALL_LIMIT: Duration = Duration::from_secs(120);

/// One carried record fragment: the tail bytes of the previous chunk plus
/// the virtual offset where that record starts (the index needs it; the
/// renderer ignores it).
#[derive(Clone, Debug, Default)]
pub struct Carried {
    /// Bytes of the record that started before this chunk. The renderer
    /// keeps all of them; the index keeps only the shared header it decodes.
    pub head: Vec<u8>,
    /// How many bytes of that record earlier chunks already hold.
    pub prefix_len: u64,
    /// Virtual offset of that record's first byte.
    pub voff: u64,
}

impl Carried {
    /// Whether there is no carried record.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.prefix_len == 0
    }
}

/// One cell: `None` until published.
type Cell = Option<Carried>;

/// The per-run chain state.
struct State {
    /// Cells by chunk id; `None` = not yet published.
    cells: Vec<Cell>,
    /// Bit per chunk id: a waiter is parked on that cell. Publishes notify
    /// only when the cell's own waiter is present — a blanket `notify_all`
    /// with a full ticket window costs a thundering herd per chunk.
    waiting: Vec<u64>,
    /// A worker failed before publishing: every waiter errors out.
    poisoned: bool,
}

/// Sets bit `id` in a word-per-64 bitmap, growing it as needed.
fn set_waiter(bits: &mut Vec<u64>, id: u64) {
    let word = (id >> 6) as usize;
    if bits.len() <= word {
        bits.resize(word + 1, 0);
    }
    bits[word] |= 1 << (id & 63);
}

/// Clears bit `id`, when set.
fn clear_waiter(bits: &mut [u64], id: u64) {
    if let Some(word) = bits.get_mut((id >> 6) as usize) {
        *word &= !(1 << (id & 63));
    }
}

/// Whether bit `id` is set.
fn has_waiter(bits: &[u64], id: u64) -> bool {
    bits.get((id >> 6) as usize).is_some_and(|word| word >> (id & 63) & 1 == 1)
}

/// The shared chain.
pub struct FrameChain {
    state: Mutex<State>,
    cv: Condvar,
    abort: AbortSignal,
    /// Set by the driver when the run ends, so stragglers do not block.
    finished: AtomicBool,
    /// Cumulative nanoseconds spent in [`FrameChain::wait`] (doc_v4 03
    /// A4-5: the chain wait must be measurable, not inferred).
    wait_ns: AtomicU64,
    /// Reference point for the publish timestamps below.
    epoch: Instant,
    /// Nanoseconds (since `epoch`) of the last publish of any cell, and the
    /// id it published — the watchdog is a *global* silence, not a per-cell
    /// wait (doc_v4 01 F4-14: a stalled NFS blocks every waiter even while
    /// some cell is being retried).
    last_publish_ns: AtomicU64,
    last_publish_id: AtomicU64,
    /// The watchdog's silence limit (production: 120 s; tests shorten it).
    stall_limit: Duration,
}

impl FrameChain {
    /// Creates the chain bound to the run's abort signal.
    #[must_use]
    pub fn new(abort: AbortSignal) -> Self {
        Self {
            state: Mutex::new(State {
                cells: vec![Some(Carried::default())],
                waiting: Vec::new(),
                poisoned: false,
            }),
            cv: Condvar::new(),
            abort,
            finished: AtomicBool::new(false),
            wait_ns: AtomicU64::new(0),
            epoch: Instant::now(),
            last_publish_ns: AtomicU64::new(0),
            last_publish_id: AtomicU64::new(0),
            stall_limit: DEFAULT_STALL_LIMIT,
        }
    }

    /// Takes the carried bytes for `id`, waiting until the previous chunk
    /// published them.
    ///
    /// # Errors
    ///
    /// [`EngineError::Aborted`] when the run aborts or a worker poisons the
    /// chain before publishing this cell.
    pub fn wait(&self, id: u64) -> Result<Carried, EngineError> {
        let started = Instant::now();
        let result = self.wait_inner(id, started);
        self.wait_ns
            .fetch_add(started.elapsed().as_nanos() as u64, Ordering::AcqRel);
        result
    }

    /// The cumulative wait time over all callers (engine metrics).
    #[must_use]
    pub fn wait_time(&self) -> Duration {
        Duration::from_nanos(self.wait_ns.load(Ordering::Acquire))
    }

    fn wait_inner(&self, id: u64, started: Instant) -> Result<Carried, EngineError> {
        // A predecessor's walk is microseconds (one length-prefix pass over
        // one chunk); anything past this safety net means another stage
        // failed without poisoning the chain, and blocking forever would
        // hang the scope join. The clock is the *last publish of any cell*
        // (01 F4-14): a slow storage that still makes progress must not trip
        // the watchdog.
        let mut state = self.state.lock().expect("invariant: chain not poisoned");
        loop {
            let now_ns = self.epoch.elapsed().as_nanos() as u64;
            let silence_ns = now_ns.saturating_sub(self.last_publish_ns.load(Ordering::Acquire));
            if u128::from(silence_ns) > self.stall_limit.as_nanos() {
                let waited = started.elapsed().as_secs_f64();
                return Err(EngineError::job(format!(
                    "chain stalled: chunk {id} waited {waited:.1}s; no cell published since \
                     chunk {} ({:.1}s of silence)",
                    self.last_publish_id.load(Ordering::Acquire),
                    silence_ns as f64 / 1e9,
                )));
            }
            if let Some(slot) = state.cells.get_mut(id as usize)
                && let Some(cell) = slot.take()
            {
                return Ok(cell);
            }
            if state.poisoned || self.finished.load(Ordering::Acquire) || self.abort.is_aborted() {
                return Err(EngineError::Aborted);
            }
            set_waiter(&mut state.waiting, id);
            let (guard, _timeout) = self
                .cv
                .wait_timeout(state, Duration::from_millis(10))
                .expect("invariant: chain not poisoned");
            state = guard;
            clear_waiter(&mut state.waiting, id);
        }
    }

    /// Publishes the carried record for `id` (chunk `id`'s head fragment).
    pub fn publish(&self, id: u64, carried: Carried) {
        self.last_publish_ns
            .store(self.epoch.elapsed().as_nanos() as u64, Ordering::Release);
        self.last_publish_id.store(id, Ordering::Release);
        let mut state = self.state.lock().expect("invariant: chain not poisoned");
        let idx = id as usize;
        if state.cells.len() <= idx {
            state.cells.resize_with(idx + 1, || None);
        }
        state.cells[idx] = Some(carried);
        // Only the cell's own waiter needs waking; the set/check pair is
        // under one lock, so no wakeup can be lost.
        let notify = has_waiter(&state.waiting, id);
        drop(state);
        if notify {
            self.cv.notify_all();
        }
    }

    /// Marks the chain failed: every current and future waiter errors out
    /// instead of blocking on a cell that will never arrive.
    pub fn poison(&self) {
        let mut state = self.state.lock().expect("invariant: chain not poisoned");
        state.poisoned = true;
        drop(state);
        self.cv.notify_all();
    }

    /// Wakes and releases every waiter at the end of the run.
    pub fn finish(&self) {
        self.finished.store(true, Ordering::Release);
        self.cv.notify_all();
    }

    /// The carried bytes still held for `id`, if any (the driver uses it to
    /// detect a record truncated by the end of the file).
    #[must_use]
    pub fn peek(&self, id: u64) -> Option<Carried> {
        let state = self.state.lock().expect("invariant: chain not poisoned");
        state.cells.get(id as usize).and_then(Clone::clone)
    }
}

/// Publishes one cell on the happy path; if the worker unwinds (a panic in
/// `run_chunk`, caught by the engine's `catch_unwind`) or returns before
/// publishing, the chain is poisoned so no waiter blocks forever. The abort
/// signal alone cannot do this: a blocked `select!` may consume its token.
pub struct PublishGuard<'a> {
    chain: &'a FrameChain,
    id: u64,
    done: bool,
}

impl FrameChain {
    /// A guard that publishes cell `id` (or poisons the chain instead).
    #[must_use]
    pub const fn publish_guard(&self, id: u64) -> PublishGuard<'_> {
        PublishGuard { chain: self, id, done: false }
    }
}

impl PublishGuard<'_> {
    /// Publishes the cell and disarms the guard.
    pub fn publish(&mut self, carried: Carried) {
        self.chain.publish(self.id, carried);
        self.done = true;
    }
}

impl Drop for PublishGuard<'_> {
    fn drop(&mut self) {
        if !self.done {
            self.chain.poison();
        }
    }
}

/// One chain per run, shared by the job.
pub type SharedChain = Arc<FrameChain>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::Ordering as _;

    /// F4-14: the watchdog fires on *global* silence and names the waiter,
    /// the waited time and the last publish.
    #[test]
    fn stall_watchdog_reports_the_global_silence() {
        let abort = crate::abort::abort();
        let mut chain = FrameChain::new(abort.signal.clone());
        chain.stall_limit = Duration::from_millis(20);
        // Cell 1 is never published: only the construction-time cell 0 exists.
        let err = chain.wait(1).expect_err("a silent chain must trip the watchdog");
        let text = err.to_string();
        assert!(text.contains("chain stalled"), "{text}");
        assert!(text.contains("chunk 1"), "{text}");
        assert!(text.contains("since chunk 0"), "{text}");

        // A publish after the wait starts resets the silence clock: the next
        // waiter is served instead of timing out.
        chain.publish(1, Carried::default());
        assert!(chain.wait(1).is_ok());
        assert!(chain.wait_time() > Duration::ZERO);
        assert_eq!(chain.last_publish_id.load(Ordering::Acquire), 1);
    }
}
