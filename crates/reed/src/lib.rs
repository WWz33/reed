//! `reed` (L4): the facade library. Wires the pure policy layers
//! (reed-stats; reed-filter) onto the parallel engine
//! (`docs/architecture-spec.md` §3.1, §7.1).
//!
//! **Embedded API**: the typed requests (`FilterSpec`, `StatsRequest`) with
//! the `run_*` entry points and their report/error types, plus the output
//! writers reachable as `reed::pixy` / `reed::reports` / `reed::summary`.
//! Everything else is implementation detail: `bcf_input` and `filtered` are
//! crate-private, `input` is hidden from the docs. The binary is the
//! product; the stability promise is on the CLI flags and output formats
//! (`docs/spec/cli.md` §0), not on these Rust items.

pub mod filter;
pub mod index;
pub mod stats;

pub(crate) mod bcf_input;
pub(crate) mod filtered;

/// Input opening stays reachable for the binary's one-shot opens
/// (`run_*` does its own); hidden from the docs.
#[doc(hidden)]
pub mod input;

/// Implementation types the integration tests drive directly (not part of
/// the embedded API).
#[doc(hidden)]
pub mod internal {
    pub use crate::bcf_input::{BcfRecordFn, InputJob, InputScratch};
    pub use crate::filtered::{FilteredJob, FilteredScratch, ScannedRecord, ScannedRecordFn};
    pub use crate::input::{open, open_regions, OpenedInput};
}

// The writers that belong to the embedded API: the pixy-layout statistics,
// the per-site/per-sample reports and the bcftools-stats summary renderer.
pub use stats::writers::{pixy, reports, summary};
