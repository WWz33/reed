//! Transport-layer errors with locus context (`docs/spec/architecture.md` §6,
//! §15). Supersedes the legacy `ReadError` once the old window path
//! retires.

use std::fmt;

/// A genomic position for error context (1-based, VCF convention).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Locus {
    /// Contig name.
    pub contig_name: String,
    /// 1-based position.
    pub pos1: u64,
}

impl fmt::Display for Locus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.contig_name, self.pos1)
    }
}

/// A capability an input lacks but the command needs (doc_v30 R06).
///
/// This is a command-line problem, not a data problem: the caller chose a
/// mode (regions, two passes) the input cannot serve. The facade maps it to
/// its usage class (exit 2) instead of a runtime failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Needed {
    /// Random access to the input; stdin cannot seek.
    Seek,
    /// A `.csi`/`.tbi` sidecar next to the input.
    Index,
    /// BGZF/BCF framing; plain VCF and regular gzip have no block offsets.
    Bgzf,
    /// A second pass over the same bytes; stdin is read once.
    Replay,
}

impl fmt::Display for Needed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Seek => "a seekable input",
            Self::Index => "an index (.csi or .tbi)",
            Self::Bgzf => "BGZF/BCF framing",
            Self::Replay => "a replayable input",
        })
    }
}

/// Errors from the reed-io transport layer.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum IoError {
    /// The input is not a supported format.
    #[error("unsupported input: {0}")]
    UnsupportedInput(
        /// What was wrong with the bytes.
        &'static str,
    ),
    /// The input lacks a capability the command needs (doc_v30 R06).
    #[error("unsupported input: needs {needed}; {fix}")]
    Capability {
        /// Which capability is missing.
        needed: Needed,
        /// How the operator can provide it.
        fix: &'static str,
    },
    /// Reading from the underlying stream failed.
    #[error("io error: {0}")]
    Io(String),
    /// A BGZF block or chunk is malformed.
    #[error("malformed BGZF data{}: {reason}", locus.as_ref().map_or(String::new(), |l| format!(" at {l}")))]
    MalformedBgzf {
        /// Where the failure was noticed.
        locus: Option<Locus>,
        /// What is wrong.
        reason: String,
    },
    /// A record line is malformed.
    #[error("malformed record at {locus}: {reason}")]
    MalformedRecord {
        /// Where the failure was noticed.
        locus: Locus,
        /// What is wrong.
        reason: String,
    },
    /// A genotype violates the declared input policy.
    #[error("policy rejected genotype at {locus}: {reason}")]
    PolicyRejected {
        /// Where the failure was noticed.
        locus: Locus,
        /// Which policy and why.
        reason: &'static str,
    },
    /// A header problem (missing #CHROM line, duplicated sample, ...).
    #[error("header error: {0}")]
    Header(String),
}

impl From<std::io::Error> for IoError {
    fn from(err: std::io::Error) -> Self {
        Self::Io(err.to_string())
    }
}

impl From<libdeflater::DecompressionError> for IoError {
    fn from(err: libdeflater::DecompressionError) -> Self {
        Self::MalformedBgzf {
            locus: None,
            reason: match err {
                libdeflater::DecompressionError::InsufficientSpace => {
                    "inflate output buffer too small (block exceeds BGZF limits)".to_string()
                }
                libdeflater::DecompressionError::BadData => "corrupt deflate stream".to_string(),
            },
        }
    }
}
