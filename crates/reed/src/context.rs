//! The per-run context every driver receives (doc_v30 P2 §2).
//!
//! Before this, each driver took its own `&ThreadBudget` and built its own
//! abort channel, so a caller could not cancel a run and every driver picked
//! its own channel capacity. One value now carries both: the caller owns the
//! [`Abort`](reed_engine::abort::Abort) (and therefore the handle) and passes
//! the receive side down.

use reed_engine::abort::AbortSignal;
use reed_engine::budget::ThreadBudget;

/// The thread/memory budget and abort signal of one run.
///
/// Construct with [`RunContext::new`] from a caller-owned
/// [`AbortSignal`](reed_engine::abort::AbortSignal), or with
/// [`RunContext::idle`] when the caller has no cancellation source (tests and
/// one-shot library calls).
#[derive(Clone, Copy, Debug)]
pub struct RunContext<'a> {
    /// Thread count and memory ceiling for this run.
    pub budget: &'a ThreadBudget,
    /// Cancellation signal; an engine error triggers the caller's handle.
    pub abort: &'a AbortSignal,
}

impl<'a> RunContext<'a> {
    /// Binds a budget and an abort signal.
    #[must_use]
    pub const fn new(budget: &'a ThreadBudget, abort: &'a AbortSignal) -> Self {
        Self { budget, abort }
    }

    /// A context whose signal never aborts.
    ///
    /// Use this when the caller has no trigger side. The pipeline detects a
    /// signal without a trigger handle and substitutes a local channel for
    /// the run, so an internal first error still wakes every waiter; the
    /// error itself is returned through the engine's first-error slot.
    #[must_use]
    pub fn idle(budget: &'a ThreadBudget) -> Self {
        Self {
            budget,
            abort: AbortSignal::idle(),
        }
    }
}
