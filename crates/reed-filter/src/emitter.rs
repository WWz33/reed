//! Output layer of the filter engine: a thin re-export of reed-io's sink.
//!
//! reed-io owns the writer stack: plain-text VCF on stdout is
//! the default (Unix filter semantics; SIGPIPE is left to the caller), `*.bcf`
//! and `*.gz` compress on a multithreaded BGZF pool sized by `threads`, and
//! anything else is a plain-text VCF file (docs/filter-architecture.md §9).
//! reed-filter adds no I/O of its own; it only hands the opened sink to
//! [`crate::exec::execute`].

pub use reed_io::sink::open as open_sink;
pub use reed_io::sink::VariantSink;
