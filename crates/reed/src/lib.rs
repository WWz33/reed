//! `reed` (L4): the facade library. Wires the pure policy layers
//! (reed-qc today; reed-stats and reed-filter as they migrate) onto the
//! parallel engine (`docs/architecture-spec.md` §3.1, §7.1).

pub mod qc;
