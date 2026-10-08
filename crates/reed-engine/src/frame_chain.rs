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

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant};

use crate::abort::AbortSignal;
use crate::error::EngineError;

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
    /// A worker failed before publishing: every waiter errors out.
    poisoned: bool,
}

/// The shared chain.
pub struct FrameChain {
    state: Mutex<State>,
    cv: Condvar,
    abort: AbortSignal,
    /// Set by the driver when the run ends, so stragglers do not block.
    finished: AtomicBool,
}

impl FrameChain {
    /// Creates the chain bound to the run's abort signal.
    #[must_use]
    pub fn new(abort: AbortSignal) -> Self {
        Self {
            state: Mutex::new(State { cells: vec![Some(Carried::default())], poisoned: false }),
            cv: Condvar::new(),
            abort,
            finished: AtomicBool::new(false),
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
        // A predecessor's walk is microseconds (one length-prefix pass over
        // one chunk); anything past this safety net means another stage
        // failed without poisoning the chain, and blocking forever would
        // hang the scope join.
        const STALL_LIMIT: Duration = Duration::from_secs(120);
        let started = Instant::now();
        let mut state = self.state.lock().expect("invariant: chain not poisoned");
        loop {
            if started.elapsed() > STALL_LIMIT {
                return Err(EngineError::Job(
                    "timed out waiting for the previous chunk's framing".into(),
                ));
            }
            if let Some(slot) = state.cells.get_mut(id as usize)
                && let Some(cell) = slot.take()
            {
                return Ok(cell);
            }
            if state.poisoned || self.finished.load(Ordering::Acquire) || self.abort.is_aborted() {
                return Err(EngineError::Aborted);
            }
            let (guard, _timeout) = self
                .cv
                .wait_timeout(state, Duration::from_millis(10))
                .expect("invariant: chain not poisoned");
            state = guard;
        }
    }

    /// Publishes the carried record for `id` (chunk `id`'s head fragment).
    pub fn publish(&self, id: u64, carried: Carried) {
        let mut state = self.state.lock().expect("invariant: chain not poisoned");
        let idx = id as usize;
        if state.cells.len() <= idx {
            state.cells.resize_with(idx + 1, || None);
        }
        state.cells[idx] = Some(carried);
        drop(state);
        self.cv.notify_all();
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
