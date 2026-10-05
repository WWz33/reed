//! Error type for the reed-core model and math modules.

/// Errors from constructing or using reed-core model types.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum CoreError {
    /// A mask or table length does not match the sample count it must cover.
    #[error("length mismatch: expected {expected}, got {got}")]
    LengthMismatch {
        /// Required length.
        expected: usize,
        /// Supplied length.
        got: usize,
    },
    /// A sample name is not present in the sample table.
    #[error("unknown sample name: {sample_name}")]
    UnknownSample {
        /// The offending name.
        sample_name: String,
    },
    /// A sample name appeared twice in a table that requires unique names.
    #[error("duplicate sample name: {sample_name}")]
    DuplicateSample {
        /// The duplicated name.
        sample_name: String,
    },
    /// A group name is not present in the group table.
    #[error("unknown group name: {group_name}")]
    UnknownGroup {
        /// The offending name.
        group_name: String,
    },
    /// A group has no samples, so no statistic can be computed over it.
    #[error("group {group_name} has no samples")]
    EmptyGroup {
        /// The empty group's name.
        group_name: String,
    },
    /// A contig is not declared in the VCF header.
    #[error("contig {contig_name} is not declared in the header")]
    UndeclaredContig {
        /// The undeclared contig name.
        contig_name: String,
    },
    /// A contig name was declared twice in the header.
    #[error("contig {contig_name} declared twice")]
    DuplicateContig {
        /// The duplicated contig name.
        contig_name: String,
    },
    /// A group name was declared twice.
    #[error("group {group_name} declared twice")]
    DuplicateGroup {
        /// The duplicated group name.
        group_name: String,
    },
    /// Reading a group file failed.
    #[error("group file error: {message}")]
    GroupFile {
        /// What went wrong and where.
        message: String,
    },
    /// Positions must be strictly ascending within a contig.
    #[error("positions must be strictly ascending: {prev_pos0} followed by {pos0}")]
    PositionNotAscending {
        /// Previous position (0-based).
        prev_pos0: u64,
        /// Offending position (0-based).
        pos0: u64,
    },
    /// A position lies beyond the declared contig length.
    #[error("position {pos1} (1-based) exceeds contig {contig_name} length {len_bp}")]
    PositionBeyondContigLength {
        /// Contig name.
        contig_name: String,
        /// Offending position (1-based, as in VCF).
        pos1: u64,
        /// Declared contig length in bp.
        len_bp: u64,
    },
    /// An interval is empty or has end before start.
    #[error("invalid interval [{start0}, {end0}): must be non-empty and ordered")]
    InvalidInterval {
        /// Interval start (0-based, inclusive).
        start0: u64,
        /// Interval end (0-based, exclusive).
        end0: u64,
    },
    /// A window grid definition is invalid (zero size or step).
    #[error("invalid window grid: size and step must both be positive")]
    InvalidWindow,
    /// An explicit input policy rejected the input as it forbids.
    #[error("policy rejected the input: {0}")]
    PolicyRejected(
        /// Human-readable rejection reason.
        &'static str,
    ),
    /// A declared policy or code path is not implemented yet.
    #[error("not implemented: {0}")]
    Unimplemented(
        /// What is missing.
        &'static str,
    ),
}
