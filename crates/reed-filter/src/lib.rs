//! Streaming hard-filter engine for reed (policy layer; transport and I/O
//! live in reed-io).
//!
//! Architecture: `docs/filter-architecture.md`. Gold standard and parity gate: `gold/`.
//!
//! Layers:
//! - [`spec`]: the sealed filter spec (value type; `from_conf` mirrors the gold pipeline config)
//! - [`plan`]: fixed schedule expansion + ascertainment derivation + provenance hash
//! - [`exec`]: the engine (two streaming passes over reed-io's record transport:
//!   reduce = genome-wide sample missing-rate reduction, scan = filter + edit + emit)
//! - [`emitter`]: thin re-export of reed-io's output sink (stdout plain-text VCF is
//!   the default; .vcf.gz / .bcf are dispatched by extension)
//!
//! The caller declares *what* (spec); the engine decides *how* (auto strategy);
//! the sink decides *where* (output stream). Every backend touch point sits in
//! reed-io; this crate depends on reed-io and reed-core only.

pub mod emitter;
pub mod exec;
pub mod plan;
pub mod spec;

pub use exec::{execute, ExecOptions, RunReport};
pub use spec::FilterSpec;

use std::fmt;

use reed_io::ReadError;

/// Unified error type. Construction-time validation reports [`FilterError::Spec`];
/// runtime failures report Io/Format/Decode.
#[derive(Debug)]
pub enum FilterError {
    /// An I/O failure while reading or writing.
    Io(std::io::Error),
    /// Spec construction/parsing failed. Validation happens here, never as a runtime panic.
    Spec(String),
    /// Input format unsupported or unrecognized.
    Format(String),
    /// Per-record decode failure, located by site.
    Decode {
        /// Contig name from the record.
        contig: String,
        /// 1-based record position.
        pos: u64,
        /// The underlying decode failure.
        msg: String,
    },
}

impl fmt::Display for FilterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FilterError::Io(e) => write!(f, "io error: {e}"),
            FilterError::Spec(msg) => write!(f, "invalid filter spec: {msg}"),
            FilterError::Format(msg) => write!(f, "unsupported input format: {msg}"),
            FilterError::Decode { contig, pos, msg } => {
                write!(f, "decode failed at {contig}:{pos}: {msg}")
            }
        }
    }
}

impl std::error::Error for FilterError {}

impl From<std::io::Error> for FilterError {
    fn from(e: std::io::Error) -> Self {
        FilterError::Io(e)
    }
}

/// Transport failures arrive through reed-io's ReadError: I/O and per-record
/// decode keep their payload, every other failure mode (unsupported format,
/// header contradictions) becomes a Format error.
impl From<ReadError> for FilterError {
    fn from(e: ReadError) -> Self {
        match e {
            ReadError::Io(err) => FilterError::Io(err),
            ReadError::Decode { contig, pos, msg } => FilterError::Decode { contig, pos, msg },
            ReadError::UnsupportedFormat(name) => FilterError::Format(format!(
                "unsupported input {name:?} (v0.1: .vcf, .vcf.gz, .vcf.bgz, .bcf)"
            )),
            other => FilterError::Format(other.to_string()),
        }
    }
}
