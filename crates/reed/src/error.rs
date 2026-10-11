//! The facade's single error type (doc_v5 08 §6.2, doc_v6 Y6): wraps the
//! lower crates' errors and classifies them for the front end's exit-code
//! mapping (docs/spec/cli.md §9: Usage = 2, Runtime = 1). Nothing here
//! prints or decides process exits — that is the CLI's job.

use reed_core::error::CoreError;
use reed_engine::error::EngineError;
use reed_filter::error::FilterError;
use reed_io::error::IoError;

/// One error type across every `run_*` entry point.
#[derive(Debug)]
pub enum Error {
    /// A model-layer failure (group/interval/sample construction).
    Core(CoreError),
    /// A transport failure (open, read, decode).
    Io(IoError),
    /// An engine failure (worker panic, memory budget, abort).
    Engine(EngineError),
    /// A filter-spec failure.
    Spec(FilterError),
    /// A usage failure the operator can fix on the command line: a bad flag
    /// combination, a missing population, a failed invariant-site
    /// declaration (cli.md §9, exit 2).
    Usage(String),
}

/// The exit-code class of an [`Error`] (cli.md §9: 2 = usage, 1 = runtime).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ErrorKind {
    /// Exit 2: the command line or input declaration is at fault.
    Usage,
    /// Exit 1: the run failed for reasons the command line cannot fix.
    Runtime,
}

impl Error {
    /// Wraps a message as a usage failure.
    #[must_use]
    pub fn usage(msg: impl Into<String>) -> Self {
        Self::Usage(msg.into())
    }

    /// The exit-code class (cli.md §9). Request-construction failures are
    /// usage: the `run_*` validation stages throw `Error::usage` there, and
    /// the model-layer `CoreError` variants are all request-shaped, so both
    /// are `Usage`. Runtime failures (I/O, decode with its `contig:pos`, the
    /// engine stages) are not — a `Spec` carrying a transport or per-record
    /// decode error rides the runtime class with them.
    #[must_use]
    pub fn kind(&self) -> ErrorKind {
        match self {
            Self::Usage(_)
            | Self::Core(_)
            | Self::Spec(FilterError::Spec(_) | FilterError::Format(_)) => ErrorKind::Usage,
            Self::Io(_)
            | Self::Spec(FilterError::Io(_) | FilterError::Decode { .. })
            | Self::Engine(_) => ErrorKind::Runtime,
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Core(e) => write!(f, "{e}"),
            Self::Io(e) => write!(f, "{e}"),
            Self::Engine(e) => write!(f, "{e}"),
            Self::Spec(e) => write!(f, "{e}"),
            Self::Usage(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Core(e) => Some(e),
            Self::Io(e) => Some(e),
            Self::Engine(e) => Some(e),
            Self::Spec(e) => Some(e),
            Self::Usage(_) => None,
        }
    }
}

impl From<CoreError> for Error {
    fn from(e: CoreError) -> Self {
        Self::Core(e)
    }
}

impl From<IoError> for Error {
    fn from(e: IoError) -> Self {
        Self::Io(e)
    }
}

impl From<EngineError> for Error {
    fn from(e: EngineError) -> Self {
        Self::Engine(e)
    }
}

impl From<FilterError> for Error {
    fn from(e: FilterError) -> Self {
        Self::Spec(e)
    }
}

/// Raw I/O (`?` in the CLI translation layer) enters the same `Io` class as
/// the typed transport errors.
impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Self::Io(IoError::from(e))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn usage_and_spec_are_usage_the_rest_runtime() {
        assert_eq!(Error::usage("x").kind(), ErrorKind::Usage);
        assert_eq!(
            Error::Core(CoreError::UnknownSample {
                sample_name: "s".into()
            })
            .kind(),
            ErrorKind::Usage
        );
        assert_eq!(
            Error::Spec(FilterError::Spec("s".into())).kind(),
            ErrorKind::Usage
        );
        assert_eq!(
            Error::Engine(EngineError::Aborted).kind(),
            ErrorKind::Runtime
        );
        assert_eq!(
            Error::Spec(FilterError::Io(std::io::Error::other("x"))).kind(),
            ErrorKind::Runtime,
            "a transport failure is runtime, not a usage error"
        );
        assert_eq!(
            Error::Spec(FilterError::Decode {
                contig: "chr1".into(),
                pos: 7,
                msg: "x".into()
            })
            .kind(),
            ErrorKind::Runtime,
            "decode carries contig:pos and exits 1 (cli.md §9)"
        );
        assert_eq!(
            Error::from(std::io::Error::other("x")).kind(),
            ErrorKind::Runtime
        );
    }
}
