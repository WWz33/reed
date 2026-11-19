//! Pure VCF text semantics shared by the transport and policy layers
//! (doc_v30 R12).
//!
//! These modules hold no I/O and no buffer ownership: the column view
//! ([`mod@line`]), the INFO tag index ([`info`]), the FORMAT layout cache
//! ([`format_layout`]), the fixed-column value parser ([`fixed`]) and the
//! byte-level edit values ([`line_edit`]). `reed-io` re-exports them and adds
//! the transport: located fixed-column errors, site views and the scanner.

pub mod fixed;
pub mod format_layout;
pub mod info;
pub mod line;
pub mod line_edit;
