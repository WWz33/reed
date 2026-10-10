//! Streaming hard-filter engine for reed (policy layer; transport and I/O
//! live in reed-io).
//!
//! Architecture: `docs/spec/filter.md`. Gold standard and parity gate: `gold/`.
//!
//! Layers:
//! - [`spec`]: the sealed filter spec (value type; `compat` mirrors the gold pipeline config)
//! - [`plan`]: fixed schedule expansion + ascertainment derivation + provenance hash
//!
//! The caller declares *what* (spec); the engine decides *how* (auto strategy);
//! the sink decides *where* (output stream). Every backend touch point sits in
//! reed-io; this crate depends on reed-io and reed-core only.

pub mod error;
pub mod eval;
pub mod expr;
pub mod funnel;
pub mod plan;
pub mod spec;

pub use error::FilterError;
pub use spec::FilterSpec;
