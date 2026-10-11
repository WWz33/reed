//! Abortion and first-error capture (`docs/spec/architecture.md` §7.5,
//! `docs/spec/v6_y5/README.md` §2 S1).
//!
//! One close-on-drop channel plus an atomic flag: the flag is the source of
//! truth (`is_aborted` never flips back when a waiter consumes a token), the
//! channel wakes every blocked `select!`, and registered byte-permit pools
//! are notified so their waiters re-check the flag immediately instead of
//! polling.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, Weak};

use crossbeam_channel::{Receiver, Sender};

use crate::budget::BytePermits;
use crate::error::EngineError;

/// The sending side: trigger abort manually, or let [`FirstError`] do it.
#[derive(Clone, Debug)]
pub struct AbortHandle {
    tx: Sender<()>,
    flag: Arc<AtomicBool>,
    /// Byte-permit pools blocked in `acquire`; notified on trigger so the
    /// abort is observed without a polling timeout (§7.5).
    permits: Arc<Mutex<Vec<Weak<BytePermits>>>>,
}

impl AbortHandle {
    /// Registers a permit pool to wake on abort. Call once per run, when the
    /// pool is created.
    pub fn register(&self, permits: &Arc<BytePermits>) {
        self.permits
            .lock()
            .expect("invariant: permit registry not poisoned")
            .push(Arc::downgrade(permits));
    }

    /// Triggers the abort signal: sets the flag, fills the channel so every
    /// blocked `select!` observes it, and wakes the registered permit
    /// waiters. A full or closed channel is fine.
    pub fn trigger(&self) {
        self.flag.store(true, Ordering::Release);
        while self.tx.try_send(()).is_ok() {}
        let mut registry = self
            .permits
            .lock()
            .expect("invariant: permit registry not poisoned");
        registry.retain(|weak| {
            let Some(permits) = weak.upgrade() else {
                return false;
            };
            permits.notify_abort();
            true
        });
    }
}

/// The receiving side, held by the engine.
#[derive(Clone, Debug)]
pub struct AbortSignal {
    rx: Receiver<()>,
    flag: Arc<AtomicBool>,
}

impl AbortSignal {
    /// Whether abort has been triggered (non-blocking check). The flag is
    /// the source of truth: consuming a channel token never clears it.
    #[must_use]
    pub fn is_aborted(&self) -> bool {
        self.flag.load(Ordering::Acquire)
    }

    /// Waits for abort, racing `other` — returns true if abort won.
    #[must_use]
    pub fn raced_by_abort(&self, other: &Receiver<()>) -> bool {
        crossbeam_channel::select! {
            recv(self.rx) -> _ => true,
            recv(other) -> _ => false,
        }
    }

    /// The raw receiver for `select!` use inside the pipeline.
    #[must_use]
    pub const fn receiver(&self) -> &Receiver<()> {
        &self.rx
    }
}

/// Creates an abort channel whose wake-up capacity covers `threads` waiters
/// (reader + workers + reducer + margin): `trigger` fills the channel so
/// *each* blocked `select!` observes the abort rather than one consumer
/// eating the only token (§7.5). `--threads` may exceed the core count, so
/// the capacity follows the budget, not `available_parallelism`.
#[must_use]
pub fn abort_channel_for(threads: usize) -> (AbortHandle, AbortSignal) {
    let (tx, rx) = crossbeam_channel::bounded::<()>(threads + 4);
    let flag = Arc::new(AtomicBool::new(false));
    (
        AbortHandle {
            tx,
            flag: Arc::clone(&flag),
            permits: Arc::new(Mutex::new(Vec::new())),
        },
        AbortSignal { rx, flag },
    )
}

/// Creates an abort channel sized for the default parallelism (tests and
/// callers without a [`crate::budget::ThreadBudget`]).
#[must_use]
pub fn abort_channel() -> (AbortHandle, AbortSignal) {
    let threads = std::thread::available_parallelism().map_or(8, std::num::NonZeroUsize::get);
    abort_channel_for(threads)
}

/// One abort channel shared by the caller and the engine: internal errors
/// trigger the same handle the caller holds, so there is no watcher thread
/// forwarding an external signal (03-architecture.md §4.3).
#[derive(Clone, Debug)]
pub struct Abort {
    /// Trigger side (external abort and internal first-error both use it).
    pub handle: AbortHandle,
    /// Receive side (every engine blocking point selects on it).
    pub signal: AbortSignal,
}

/// Creates one shared abort channel.
#[must_use]
pub fn abort() -> Abort {
    let (handle, signal) = abort_channel();
    Abort { handle, signal }
}

/// Creates one shared abort channel sized by the run's thread budget.
#[must_use]
pub fn abort_for(budget: &crate::budget::ThreadBudget) -> Abort {
    let (handle, signal) = abort_channel_for(budget.threads);
    Abort { handle, signal }
}

/// First error wins: records the first error and triggers abort; later
/// errors are dropped (they are consequences, not causes).
#[derive(Debug)]
pub struct FirstError {
    err: Mutex<Option<EngineError>>,
    abort: AbortHandle,
}

impl FirstError {
    /// Creates the slot bound to an abort handle.
    #[must_use]
    pub fn new(abort: AbortHandle) -> Self {
        Self {
            err: Mutex::new(None),
            abort,
        }
    }

    /// Records `err` if it is the first, and triggers abort.
    pub fn record(&self, err: EngineError) {
        let mut slot = self.err.lock().expect("invariant: error slot not poisoned");
        if slot.is_none() {
            *slot = Some(err);
        }
        self.abort.trigger();
    }

    /// Whether an error has been recorded.
    #[must_use]
    pub fn failed(&self) -> bool {
        self.err
            .lock()
            .expect("invariant: error slot not poisoned")
            .is_some()
    }

    /// Takes the recorded error, if any.
    #[must_use]
    pub fn take(self) -> Option<EngineError> {
        self.err
            .into_inner()
            .expect("invariant: error slot not poisoned")
    }

    /// Same as [`FirstError::take`] through an `Arc` (for the pipeline, which
    /// shares the slot across threads).
    #[must_use]
    pub fn take_shared(this: &Arc<Self>) -> Option<EngineError> {
        // The run has joined all threads, so the slot is uncontended.
        this.err
            .lock()
            .expect("invariant: error slot not poisoned")
            .take()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_error_wins_and_aborts() {
        let (handle, signal) = abort_channel();
        let first = FirstError::new(handle);
        assert!(!first.failed());
        assert!(!signal.is_aborted());
        first.record(EngineError::Aborted);
        first.record(EngineError::job("later"));
        assert!(first.failed());
        assert!(signal.is_aborted());
        assert!(matches!(first.take(), Some(EngineError::Aborted)));
    }

    #[test]
    fn manual_abort_is_observed() {
        let (handle, signal) = abort_channel();
        handle.trigger();
        assert!(signal.is_aborted());
    }

    /// The flag survives token consumption: `is_aborted` must not flip back
    /// when a `select!` waiter eats the channel token (A6-1).
    #[test]
    fn consuming_the_token_does_not_clear_the_flag() {
        let (handle, signal) = abort_channel();
        handle.trigger();
        while signal.receiver().try_recv().is_ok() {}
        assert!(
            signal.is_aborted(),
            "flag must remain set after the channel drains"
        );
    }

    /// A registered permit pool is woken by `trigger`, so `acquire` returns
    /// without waiting out any polling timeout (A6-1, §7.5).
    #[test]
    fn trigger_wakes_registered_permit_waiters() {
        let (handle, signal) = abort_channel();
        let permits = Arc::new(BytePermits::new(0));
        handle.register(&permits);
        let waiter = {
            let permits = Arc::clone(&permits);
            let signal = signal.clone();
            std::thread::spawn(move || permits.acquire(1, &signal))
        };
        // Give the waiter time to block inside `acquire`, then abort.
        std::thread::sleep(std::time::Duration::from_millis(20));
        let started = std::time::Instant::now();
        handle.trigger();
        assert!(!waiter.join().expect("invariant: waiter does not panic"));
        assert!(
            started.elapsed() < std::time::Duration::from_secs(1),
            "acquire must wake on trigger, took {:?}",
            started.elapsed()
        );
    }
}
