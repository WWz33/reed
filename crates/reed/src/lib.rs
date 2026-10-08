//! `reed` (L4): the facade library. Wires the pure policy layers
//! (reed-qc, reed-stats; reed-filter as it migrates) onto the
//! parallel engine (`docs/architecture-spec.md` §3.1, §7.1).

pub mod bcf_input;
pub mod filter;
pub mod filtered;
pub mod index;
pub mod input;
pub mod stats;

// Historical public paths (doc_v2 §8): the writers moved under `stats::writers`.
pub use stats::writers::{pixy, qc, qc_report, reports};
