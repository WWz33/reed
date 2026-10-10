//! Abortion and first-error capture (`docs/spec/architecture.md` §7.5).
//!
//! One close-on-drop channel: recording the first error triggers abort, and
//! every blocking point in the engine selects on the abort signal, so a
//! failure stops the whole run immediately (no polling).

use std::sync::{Arc, Mutex};

use crossbeam_channel::{Receiver, Sender};

use crate::error::EngineError;

/// The sending side: trigger abort manually, or let [`FirstError`] do it.
#[derive(Clone, Debug)]
pub struct AbortHandle {
    tx: Sender<()>,
}

impl AbortHandle {
    /// Triggers the abort signal: fills the channel so every blocked
    /// waiter observes it; a full or closed channel is fine.
    pub fn trigger(&self) {
        while self.tx.try_send(()).is_ok() {}
    }
}

/// The receiving side, held by the engine.
#[derive(Clone, Debug)]
pub struct AbortSignal {
    rx: Receiver<()>,
}

impl AbortSignal {
    /// Whether abort has been triggered (non-blocking check).
    #[must_use]
    pub fn is_aborted(&self) -> bool {
        !self.rx.is_empty()
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

/// Creates an abort channel: the handle may be cloned to multiple owners.
///
/// The capacity covers every thread that may block on the signal inside one
/// run (reader + workers + reducer + a safety margin): `trigger` fills the
/// channel so *each* blocked `select!` observes the abort rather than one
/// consumer eating the only token (§7.5 "中止即时生效").
#[must_use]
pub fn abort_channel() -> (AbortHandle, AbortSignal) {
    let waiters = std::thread::available_parallelism().map_or(8, std::num::NonZeroUsize::get) + 4;
    let (tx, rx) = crossbeam_channel::bounded::<()>(waiters);
    (AbortHandle { tx }, AbortSignal { rx })
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
}
