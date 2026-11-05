//! Output layer (`docs/rounds/v2/03-architecture.md` §5): the worker-local
//! [`BlockWriter`], the reducer-side [`OrderedOutput`], the index writer and
//! the atomic output file, shared by every command that emits VCF/BGZF or a
//! report.

pub mod atomic;
pub mod block_writer;
pub mod index;
pub mod ordered_output;

pub use block_writer::{
    BGZF_BLOCK_SIZE, BlockWriter, ChunkOutput, IndexKey, IndexRecord, OutBlock,
};
pub use index::{Indexer, WriteIndex};
pub use atomic::{AtomicFile, write_atomic};
pub use ordered_output::OrderedOutput;
