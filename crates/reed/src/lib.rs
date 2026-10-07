//! `reed` (L4): the facade library. Wires the pure policy layers
//! (reed-qc, reed-stats; reed-filter as it migrates) onto the
//! parallel engine (`docs/architecture-spec.md` §3.1, §7.1).

pub mod bcf_input;
pub mod filter;
pub mod index;
pub mod input;
pub mod pixy;
pub mod qc;
pub mod reports;
pub mod qc_report;
pub mod stats;
