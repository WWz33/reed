//! BGZF block primitives: header parsing, libdeflate inflate/deflate
//! (`docs/spec/architecture.md` §6.2, §6.6).

pub mod block;
pub mod deflate;
pub mod inflate;
