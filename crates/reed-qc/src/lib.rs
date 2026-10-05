//! `reed-qc`: QC summary statistics (L3), pure — no I/O, no threads
//! (`docs/architecture-spec.md` §10). Everything here satisfies the
//! accumulator laws (identity, associativity, split invariance), so QC runs
//! on the same chunk-parallel driver as statistics.

pub mod summary;
