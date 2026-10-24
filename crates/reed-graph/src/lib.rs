//! Trees over genetic distance matrices.
//!
//! Two constructions are provided, each matching a reference implementation
//! bit-for-bit:
//!
//! - [`upgma`] — average linkage, the same Lance–Williams recurrence
//!   `stats::hclust(method = "average")` uses (which is what `poppr::upgma`
//!   calls; `poppr` does not implement UPGMA itself).
//! - [`nj`] — neighbor-joining, the same agglomeration `ape::nj` performs.
//!
//! The [`Tree`] type follows `ape`'s `phylo` shape so that trees round-trip
//! through newick without a translation layer, and [`to_newick`] emits the
//! `ape::write.tree` format.
//!
//! # Scope
//!
//! This crate consumes a distance matrix and produces a tree. It deliberately
//! does not depend on `reed-distance`: how the distances were computed is the
//! caller's business, which keeps every algorithm here testable with a
//! hand-written matrix.
//!
//! Bootstrap support is in [`bootstrap`]: loci are resampled with replacement,
//! the tree is rebuilt per replicate, and each branch is counted
//! (`ape::boot.phylo`'s procedure).
//!
//! Not here (yet): minimum spanning networks and distance-threshold
//! clustering (`docs/rounds/v17b` §3).

#![allow(clippy::needless_range_loop)]
// Node ids are bounded by the sample count, which is far below u32::MAX.
#![allow(clippy::cast_possible_truncation)]

pub mod bootstrap;
pub mod error;
pub mod newick;
pub mod tree;

pub use bootstrap::{support, BootstrapConfig, BootstrapReport, TreeMethod};
pub use error::GraphError;
pub use newick::to_newick;
pub use tree::{Tree, nj, upgma};
