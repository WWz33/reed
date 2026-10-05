//! Engine errors (`docs/architecture-spec.md` §7.5).

/// Errors from the parallel execution engine.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum EngineError {
    /// A transport/source error.
    #[error("io: {0}")]
    Io(String),
    /// A worker panicked while processing a chunk.
    #[error("worker panicked on chunk {chunk}: {message}")]
    WorkerPanicked {
        /// The chunk being processed.
        chunk: reed_core::ChunkId,
        /// The panic payload, when it carried a string.
        message: String,
    },
    /// In-flight bytes exceeded the memory budget (§7.7); never a silent OOM.
    #[error("memory budget exceeded: in-flight {peak_bytes} bytes > budget {budget_bytes} bytes")]
    MemBudgetExceeded {
        /// The observed peak.
        peak_bytes: u64,
        /// The configured budget.
        budget_bytes: u64,
    },
    /// A job-level (policy/decode) failure.
    #[error("job failed: {0}")]
    Job(String),
    /// The run was aborted before completion.
    #[error("aborted")]
    Aborted,
}

impl From<reed_io::IoError> for EngineError {
    fn from(err: reed_io::IoError) -> Self {
        Self::Io(err.to_string())
    }
}
