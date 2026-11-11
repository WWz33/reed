//! BGZF block primitives: header parsing, libdeflate inflate/deflate
//! (`docs/spec/architecture.md` §6.2, §6.6), and virtual-offset addressing
//! (doc_v28 V28-27).

pub mod address;
pub mod block;
pub mod deflate;
pub mod inflate;
