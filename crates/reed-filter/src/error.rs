//! The crate's single error type (`docs/spec/naming.md` §2.7: one `*Error`
//! per crate, defined in `error.rs` like the other layers).

use std::fmt;

use reed_io::IoError;

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

/// Transport failures arrive through reed-io's [`IoError`].
impl From<IoError> for FilterError {
    fn from(e: IoError) -> Self {
        FilterError::Io(std::io::Error::other(e.to_string()))
    }
}
