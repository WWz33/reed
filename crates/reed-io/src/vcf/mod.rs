//! VCF text fast paths: line column views, fixed-column parsing, site views,
//! INFO fragment iteration, the FORMAT layout cache and the single-pass
//! genotype scanner (`docs/rounds/v2/03-architecture.md` §3).

pub mod fixed;
pub mod format_layout;
pub mod genotype_scan;
pub mod info;
pub mod line;
pub mod line_edit;
pub mod site;
