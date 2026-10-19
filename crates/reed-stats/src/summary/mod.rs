//! QC summary statistics (spec §10), merged into `reed-stats` per
//! 04-cleanup §C4: pure — no I/O, no threads — and satisfying the
//! accumulator laws (identity, associativity, split invariance), so QC runs
//! on the same chunk-parallel driver as the population statistics.

//! The `summary` section (§10.1): record and variant-class counts, the
//! bcftools-stats `SN` analogue.

pub mod histogram;
pub mod report;

pub use report::SummaryAcc;
