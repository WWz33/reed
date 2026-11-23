//! VCF text fast paths: line column views, fixed-column parsing, site views,
//! INFO fragment iteration, the FORMAT layout cache and the single-pass
//! genotype scanner (`docs/rounds/v2/03-architecture.md` §3).
//!
//! The pure text semantics moved to `reed-core::vcf` (doc_v30 R12) so the
//! filter layer can depend on L0 alone; this module re-exports them and adds
//! the transport pieces that need located errors ([`fixed`]), site views
//! ([`site`]) and the scanner ([`genotype_scan`]).
pub use reed_core::vcf::{format_layout, info, line, line_edit};

pub mod fixed;
pub mod genotype_scan;
pub mod site;
