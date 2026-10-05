//! VCF text fast paths: line column views, fixed-column parsing, INFO
//! fragment iteration and the single-pass sample-column scanner
//! (`docs/architecture-spec.md` §6.3–§6.4).

pub mod fixed;
pub mod gt;
pub mod info;
pub mod line;
pub mod sample_scan;
