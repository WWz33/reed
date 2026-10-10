//! BCFv2.2 decoding and VCF text rendering (doc_v3 02-bcf.md).
//!
//! The behaviour follows hts-specs BCFv2.2 and htslib's `vcf.c`
//! (`bcf_fmt_array`, `bcf_format_gt_v2`, `kputd`, `vcf_format1`); reed does
//! not invent encoding rules. Records are decoded lazily, field by field,
//! and rendered straight into a byte buffer — no per-record allocation and
//! no intermediate objects.

pub mod dict;
pub mod encode;
pub mod gt;
pub mod native_encode;
pub mod record;
pub mod render;

pub use dict::BcfDict;
pub use encode::{BcfEncoder, bcf_header_text};
pub use gt::{GtQuery, GtTable, SiteHeader, decode_record_gt, format_key, site_header};
pub use native_encode::NativeEncoder;
pub use record::{BcfRecordRef, VirtualMap, declared_len, record_len, record_ref, record_site};
pub use render::{
    decode_format_columns, render_record, render_sample_part, render_site_line, render_site_part,
};
