//! `reed` (L4): the facade library. Wires the pure policy layers
//! (reed-qc, reed-stats; reed-filter as it migrates) onto the
//! parallel engine (`docs/architecture-spec.md` §3.1, §7.1).

pub mod filter;
pub mod qc;
pub mod qc_report;
pub mod stats;
