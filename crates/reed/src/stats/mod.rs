//! The statistics drivers (`doc_v2/03-architecture.md` §8): the chunk-parallel
//! job plus the output writers (pixy layout, plink2/vcftools reports, the
//! bcftools-stats QC summary).
//!
//! The writers keep their historical module names and public paths
//! (`reed::reports`, `reed::pixy`, `reed::qc_report`) through the crate-root
//! re-exports; splitting them further (per-format files) is pending.

pub mod job;
pub mod writers;

pub use job::*;
