//! Output layer (`doc_v2/03-architecture.md` §5): the worker-local
//! [`BlockWriter`], the reducer-side [`OrderedOutput`] and the index writer,
//! shared by every command that emits VCF/BGZF.

pub mod block_writer;
pub mod index;
pub mod ordered_output;

pub use block_writer::{BGZF_BLOCK_SIZE, BlockWriter, ChunkOutput, IndexKey, IndexRecord, OutBlock};
pub use index::{Indexer, WriteIndex};
pub use ordered_output::OrderedOutput;
