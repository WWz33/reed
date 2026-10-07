//! BCFv2.2 decoding and VCF text rendering (doc_v3 02-bcf.md).
//!
//! The behaviour follows hts-specs BCFv2.2 and htslib's `vcf.c`
//! (`bcf_fmt_array`, `bcf_format_gt_v2`, `kputd`, `vcf_format1`); reed does
//! not invent encoding rules. Records are decoded lazily, field by field,
//! and rendered straight into a byte buffer — no per-record allocation and
//! no intermediate objects.

pub mod dict;
pub mod record;
pub mod render;

pub use dict::BcfDict;
pub use record::{declared_len, record_len, record_ref, record_site, BcfRecordRef, VirtualMap};
pub use render::render_record;
