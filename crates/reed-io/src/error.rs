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
    /// Reading from the underlying stream failed.
    #[error("io error: {0}")]
    Io(String),
    /// A BGZF block or chunk is malformed.
    #[error("malformed BGZF data{}: {reason}", locus.as_ref().map_or(String::new(), |l| format!(" at {l}")))]
    MalformedBgzf {
        /// Where the failure was noticed.
        locus: Option<Locus>,
        /// What is wrong.
        reason: &'static str,
    },
    /// A record line is malformed.
    #[error("malformed record at {locus}: {reason}")]
    MalformedRecord {
        /// Where the failure was noticed.
        locus: Locus,
        /// What is wrong.
        reason: &'static str,
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
                    "inflate output buffer too small (block exceeds BGZF limits)"
                }
                libdeflater::DecompressionError::BadData => "corrupt deflate stream",
            },
        }
    }
}
