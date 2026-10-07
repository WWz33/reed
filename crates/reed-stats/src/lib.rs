//! `reed-stats` (L3): population-genetics statistics as additive
//! accumulators over per-site group counts (`docs/architecture-spec.md`
//! §8–§9). Every estimator is a sealed definition: the accumulator carries
//! its formula, denominator and missing-data semantics, and `merge` is
//! component-wise addition (exact integers, [`reed_core::math::sum::NeumaierSum`]
//! for floats), so per-chunk partials reduce deterministically.
//!
//! The single input view is a row of [`reed_core::count::SiteGroupCounts`]
//! (one entry per group) plus a site-wide "any ALT called" flag — the same
//! integer quantities the frozen v0.2 estimators derived from
//! `WindowBuffer::allele_counts`, so parity with the oracle
//! (`reed-core-v02-oracle`) is exact by construction.

// Count-to-float conversions are the domain here: every estimator turns
// exact integer tallies into f64 once per site or once per window. The
// results are covered by the 1e-12 oracle-parity gate, not by exactness.
#![allow(clippy::cast_precision_loss)]

pub mod dxy;
pub mod fst;
pub mod pi;
pub mod individual;
pub mod jackknife;
pub mod ld;
pub mod patterson;
pub mod pbs;
pub mod request;
pub mod sfs;
pub mod tajima;
pub mod theta;
pub mod window;

pub use dxy::DxyAcc;
pub use fst::{FstAcc, FstKind, FstMatrixAcc, FstMatrixEntry, HudsonFiniteAcc, NeiMultiAcc, WcMultiAcc};
pub use pi::PiAcc;
pub use individual::{IndividualFAcc, JointSfsAcc};
pub use jackknife::{BlockJackknife, BlockPartial};
pub use ld::{LdDecayAcc, LdPruner, PairSums, pair_sums};
pub use patterson::{F2Acc, F3Acc, F4Acc, PattersonKind, f4_ratio};
pub use pbs::PbsAcc;
pub use request::{StatAcc, StatKey, StatKind, build_accs};
pub use sfs::{ProjectionCache, ProjectedSfsAcc, SfsAcc, projection_weights};
pub use tajima::TajimaAcc;
pub use theta::ThetaAcc;
pub use window::{BpBinning, ContigBins};
