//! `reed-core`: canonical population-genetics statistics.
//!
//! Design (see `reed/docs/rounds/v1/design.md`):
//!
//! - L2: [`WindowBuffer`] materialises decoded genotypes for one window.
//! - L3: estimators are sealed definitions (formula + denominator + missing-data
//!   semantics + citation); availability is derived from a [`SiteCompleteness`]
//!   manifest; results always carry exact integer counts, never bare floats.

pub mod bits;
pub mod estimators;
pub mod manifest;
pub mod population;
pub mod result;
pub mod window;

pub use bits::Bits;
pub use manifest::{Ascertainment, Availability, SiteCompleteness};
pub use population::{PopError, PopulationView};
pub use result::{DxyResult, FstResult, PairCounts, PiResult};
pub use window::{AlleleCounts, ContigId, GenotypeCounts, WindowBuffer, WindowError};
