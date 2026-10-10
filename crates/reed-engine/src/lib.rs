//! `reed-engine` (L2): the parallel execution base
//! (`docs/spec/architecture.md` §7). Policy crates never touch threads; the
//! facade wraps their pure logic into [`pipeline::ChunkJob`] and
//! [`pipeline::OrderedSink`] and calls [`pipeline::run`].

pub mod abort;
pub mod budget;
pub mod error;
pub mod frame_chain;
pub mod metrics;
pub mod pipeline;
pub mod reorder;

// The crate's central types, reachable at the root (naming.md §2.9: `pub mod`
// per concept plus a short `pub use` list).
pub use budget::ThreadBudget;
pub use error::EngineError;
pub use pipeline::{ChunkJob, OrderedSink, run};
