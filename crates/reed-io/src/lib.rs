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
pub mod write;

// Public re-export: reed-io is the single noodles boundary, and downstream
// crates (reed-filter) have no direct noodles dependency.
pub use noodles::vcf as noodles_vcf;
pub use error::{IoError, Locus};

/// The total byte length of a BGZF block from its 18-byte header, when the
/// header is present in `bytes`.
#[must_use]
pub fn bgzf_total_size(bytes: &[u8]) -> Option<usize> {
    crate::bgzf::block::parse_block_size(bytes)
}
