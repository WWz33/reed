//! The statistics drivers (`docs/rounds/v2/03-architecture.md` §8): the chunk-parallel
//! job plus the output writers (pixy layout, plink2/vcftools reports, the
//! bcftools-stats QC summary).
//!
//! The writers are reachable as `reed::pixy` / `reed::reports` /
//! `reed::summary` through the crate-root re-exports; splitting them
//! further (per-format files) is pending.

pub mod job;
pub mod writers;

pub use job::*;
