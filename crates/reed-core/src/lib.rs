//! `reed-core`: model + kernels (L0) and the legacy estimator library.
//!
//! Two layers live here during the architecture migration
//! (`docs/architecture-spec.md` §4):
//!
//! - The **target model** (L0): identifiers, coordinates, bit sets, sample
//!   masks, group tables, genotype bit planes, the popcount counting kernel,
//!   manifest v2, window grids and the math helpers. No I/O except the one
//!   sanctioned group-file reader in [`group`].
//! - The **legacy estimators** (v0.2, frozen semantics): [`window`] +
//!   [`estimators`], consumed by reed-io/reed-filter today and replaced by
//!   reed-stats in M4. [`ContigId`] at the crate root is still the legacy
//!   first-appearance id; the header-order id lives at [`id::ContigId`].
//!
//! Design (see `reed/DESIGN.md`):
//!
//! - [`WindowBuffer`] materialises decoded genotypes for one window.
//! - estimators are sealed definitions (formula + denominator + missing-data
//!   semantics + citation); availability is derived from a
//!   [`SiteCompleteness`] manifest; results always carry exact integer
//!   counts, never bare floats.

pub mod bits;
pub mod bitset;
pub mod coord;
pub mod count;
pub mod error;
pub mod estimators;
pub mod genotype;
pub mod grid;
pub mod group;
pub mod id;
pub mod manifest;
pub mod math;
pub mod population;
pub mod result;
pub mod sample;
pub mod site;
pub mod window;

pub use bits::Bits;
pub use bitset::Bitset;
pub use coord::{Interval0, Interval1, Pos0, Pos1};
pub use count::{CountTable, SiteGroupCounts};
pub use error::CoreError;
pub use genotype::{
    GenotypePolicies, HalfCallPolicy, MultiallelicPolicy, Ploidy, PloidyPolicy, PolicyCounters,
};
pub use grid::{AccessibilityMask, BinId, BpWindow, RegionSet, SitesWindow, WindowGrid};
pub use group::{Group, GroupFileOptions, GroupTable};
pub use id::{ChunkId, GroupId, RecordOrdinal, SampleId, TagId};
pub use manifest::{Ascertainment, Availability, SiteCompleteness};
pub use population::{PopError, PopulationView};
pub use result::{DxyResult, FstResult, PairCounts, PiResult};
pub use sample::{SampleMask, SampleTable};
pub use site::{SiteFlags, SiteKey, VariantClass};
pub use window::{AlleleCounts, ContigId, GenotypeCounts, WindowBuffer, WindowError};
