//! Engine errors (`doc_v2/03-architecture.md` §4.4): typed, with the original
//! error preserved so the `source` chain survives to the caller.

/// Errors from the parallel execution engine.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum EngineError {
    /// A transport or output failure (the original `io::Error`).
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    /// A typed transport error from `reed-io` (carries the record locus).
    #[error("{0}")]
    Input(#[from] reed_io::IoError),
    /// A job-level (policy/decode) failure; the boxed error keeps its type
    /// and `source` chain.
    #[error("job failed: {0}")]
    Job(Box<dyn std::error::Error + Send + Sync>),
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
    /// The run was aborted before completion.
    #[error("aborted")]
    Aborted,
}

impl From<reed_core::error::CoreError> for EngineError {
    fn from(err: reed_core::error::CoreError) -> Self {
        Self::Job(Box::new(err))
    }
}

impl EngineError {
    /// A job failure from a message (the common constructor for policy and
    /// decode rejections).
    #[must_use]
    pub fn job(message: impl Into<String>) -> Self {
        Self::Job(message.into().into())
    }
}
