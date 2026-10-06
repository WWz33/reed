//! reed-io: formats, BGZF, chunk sources, VCF text decoding, editing,
//! output and indexing (L1). The single noodles boundary of the workspace;
//! downstream crates take `noodles_vcf` from here.

pub mod bgzf;
pub mod chunk;
pub mod error;
pub mod field_plan;
pub mod format;
pub mod header;
pub mod vcf;

// Public re-export: reed-io is the single noodles boundary, and downstream
// crates (reed-filter) have no direct noodles dependency.
pub use noodles::vcf as noodles_vcf;
pub use error::{IoError, Locus};
