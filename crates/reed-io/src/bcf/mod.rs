//! BCFv2.2 codecs: decode to VCF text or typed columns, and encode from VCF
//! text or a framed BCF record (doc_v3 02-bcf.md).
//!
//! The behaviour follows hts-specs BCFv2.2 and htslib's `vcf.c`
//! (`bcf_fmt_array`, `bcf_format_gt_v2`, `kputd`, `vcf_format1`); reed does
//! not invent encoding rules. Records are decoded lazily, field by field,
//! and rendered straight into a byte buffer — no per-record allocation and
//! no intermediate objects.
//!
//! # Module map (doc_v22 R22-2)
//!
//! The text side (`crate::vcf`) is a line scanner; BCF is not symmetric
//! because it is length-framed, has a dictionary header, and needs real
//! encoders. The files therefore follow *operations*, not the text module's
//! layers:
//!
//! ```text
//! read/decode
//!   record     framing (l_shared/l_indiv), Cursor, sentinels
//!   dict       BCF header dictionaries (contig/tag ids)
//!   gt         site header + GT lookup table (built from the text decoder)
//!   render     BCF -> VCF text (byte-compatible with bcftools view)
//!   columns    native FORMAT columns for the filter mask (no text pass)
//!
//! write/encode
//!   encode_noodles   text -> BCF through noodles (generic path)
//!   encode_native    edit-as-you-encode fast path (`-Ob`; ~800 µs/record
//!                    cheaper than re-parsing the formatted values)
//!
//! engine adapter (outside this module)
//!   `reed::bcf_input`  chunk framing -> rendered text -> line jobs
//! ```
//!
//! Two codecs are intentional: reed avoids materializing a `bcf1_t`-style
//! record, so the text scanner and the native BCF decode are separate
//! implementations of the same semantics. `tests/bcf_vs_noodles.rs`,
//! `tests/bcf_encode.rs` and `tests/format_columns_vs_bcf.rs` guard the
//! pairing (doc_v22 R22-4); there is no 1:1 file mirror on purpose.

pub mod columns;
pub mod dict;
pub mod encode_native;
pub mod encode_noodles;
pub mod gt;
pub(crate) mod gt_decode;
pub mod record;
pub mod render;

pub use columns::decode_format_columns;
pub use dict::BcfDict;
pub use encode_native::NativeEncoder;
pub use encode_noodles::{BcfEncoder, bcf_file_header, bcf_header_text};
pub use gt::{GtQuery, GtTable, SiteHeader, decode_record_gt, format_key, site_header};
pub use record::{BcfRecordRef, VirtualMap, declared_len, record_len, record_ref, record_site};
pub use render::{render_record, render_sample_part, render_site_line, render_site_part};
