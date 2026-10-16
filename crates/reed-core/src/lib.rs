//! `reed-core`: model + kernels + math (L0).
//!
//! Identifiers, coordinates, bit sets, sample masks, group tables, genotype
//! bit planes, the per-site counting primitives (`count`), window grids and
//! the math helpers. No I/O except the one sanctioned group-file reader in
//! [`group`]. Estimator definitions live in reed-stats.

pub mod args;
pub mod bitset;
pub mod coord;
pub mod count;
pub mod error;
pub mod genotype;
pub mod grid;
pub mod group;
pub mod id;
pub mod iupac;
pub mod math;
pub mod sample;
pub mod site;
pub mod tag_table;

pub use bitset::Bitset;
pub use coord::{Interval0, Interval1, Pos0, Pos1};
pub use count::SiteGroupCounts;
pub use error::CoreError;
pub use genotype::{
    GenotypePolicies, GtCode, HalfCallPolicy, MultiallelicPolicy, Ploidy, PloidyPolicy,
    PolicyCounters,
};
pub use grid::{BinId, BpWindow, RegionSet, SitesWindow, WindowGrid};
pub use group::{Group, GroupFileOptions, GroupTable};
pub use id::{ChunkId, ContigId, GroupId, RecordOrdinal, SampleId, TagId};
pub use sample::{SampleMask, SampleTable};
pub use site::{SiteFlags, SiteKey, VariantClass};
