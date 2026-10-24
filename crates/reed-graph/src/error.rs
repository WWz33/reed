//! Errors from tree construction and serialization.

/// Errors from [`crate::upgma`], [`crate::nj`] and [`crate::to_newick`].
#[derive(Clone, Debug, PartialEq)]
pub enum GraphError {
    /// The distance matrix was not square.
    NotSquare {
        /// Rows in the supplied matrix.
        rows: usize,
        /// Columns in the supplied matrix.
        cols: usize,
    },
    /// Fewer than the minimum number of taxa were supplied.
    TooFewTaxa {
        /// Taxa present.
        got: usize,
        /// Taxa required.
        need: usize,
    },
    /// An entry was `NaN` or infinite.
    NonFiniteDistance {
        /// Row of the offending entry.
        row: usize,
        /// Column of the offending entry.
        col: usize,
    },
    /// `d[i][j]` and `d[j][i]` disagreed beyond tolerance.
    Asymmetric {
        /// Row of the offending entry.
        row: usize,
        /// Column of the offending entry.
        col: usize,
    },
    /// A diagonal entry was not zero.
    NonZeroDiagonal {
        /// Index of the offending entry.
        index: usize,
    },
    /// A branch length was `NaN` or infinite when serializing.
    BadBranchLength {
        /// Index of the offending edge.
        edge: usize,
    },
    /// The tree nesting exceeded the serializer's depth budget.
    TreeTooDeep {
        /// The depth that was reached.
        depth: usize,
    },
}

impl std::fmt::Display for GraphError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotSquare { rows, cols } => {
                write!(f, "the distance matrix is {rows}x{cols}, expected square")
            }
            Self::TooFewTaxa { got, need } => {
                write!(f, "need at least {need} taxa, got {got}")
            }
            Self::NonFiniteDistance { row, col } => {
                write!(f, "distance ({row}, {col}) is not finite")
            }
            Self::Asymmetric { row, col } => {
                write!(f, "distance ({row}, {col}) differs from its transpose")
            }
            Self::NonZeroDiagonal { index } => {
                write!(f, "distance ({index}, {index}) is not zero")
            }
            Self::BadBranchLength { edge } => {
                write!(f, "branch length {edge} is not finite")
            }
            Self::TreeTooDeep { depth } => {
                write!(f, "the tree nests {depth} levels deep (limit 10000)")
            }
        }
    }
}

impl std::error::Error for GraphError {}
