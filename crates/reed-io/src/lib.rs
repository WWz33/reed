//! reed-io: formats, BGZF, chunk sources, VCF text decoding, editing,
//! output and indexing (L1). The single noodles boundary of the workspace;
//! crates that need noodles types declare the dependency themselves (08 §4
//! L-3: no noodles re-export).

pub mod bcf;
pub mod bgzf;
pub mod chunk;
pub mod error;
pub mod field_plan;
pub mod format;
pub mod header;
pub mod vcf;
pub mod write;

pub use error::{IoError, Locus};

/// The total byte length of a BGZF block from its 18-byte header, when the
/// header is present in `bytes`.
#[must_use]
pub fn bgzf_total_size(bytes: &[u8]) -> Option<usize> {
    crate::bgzf::block::parse_block_size(bytes)
}
