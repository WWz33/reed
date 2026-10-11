//! `reed` (L4): the facade library — the backend API boundary of the
//! front/back-end split (doc_v5 08 §6). The `reed` binary is the only front
//! end and depends on this crate alone. Nothing here prints to stderr or
//! knows about clap or exit codes.
//!
//! This file is a re-export list (08 §6.2), not a module tree: the typed
//! requests, the `run_*` entry points and their reports, the output writers,
//! the expression entry point and the parallel budgets. The implementation
//! modules (`input`, `bcf_input`, `filtered`, `stats`, `filter`, `index`)
//! are private; their items are `pub` but only reachable through the
//! re-exports below, so the crate-wide `unreachable_pub` lint is allowed
//! here by construction.
//!
//! # Examples
//!
//! Build a request from types and render it back to the CLI spelling
//! (`to_args` is the text contract the provenance header records):
//!
//! ```
//! use reed::spec::{FilterSpec, MissingRate};
//! let mut spec = FilterSpec::default();
//! spec.site_after_mask.max_missing = Some(MissingRate::parse("0.1").expect("valid"));
//! let args = spec.to_args();
//! assert!(args.windows(2).any(|w| w[0] == "--geno" && w[1] == "0.1"), "{args:?}");
//! ```
//!
//! A statistics request is a plain struct; unnamed fields take their default:
//!
//! ```
//! use reed::{StatKind, StatsRequest};
//! let request = StatsRequest {
//!     kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
//!     ..StatsRequest::default()
//! };
//! assert_eq!(request.kinds.len(), 1);
//! ```
//!
//! Drive a run (the caller opens the input path; `run_*` never prints):
//!
//! ```no_run
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let report = reed::run_filter(
//!     Some(std::path::Path::new("in.vcf.gz")),
//!     &reed::spec::FilterSpec::default(),
//!     &reed::FilterOptions::default(),
//! )?;
//! let _ = report;
//! # Ok(())
//! # }
//! ```

#![allow(unreachable_pub)]

mod bcf_input;
mod error;
mod filter;
mod filtered;
mod index;
mod input;
mod stats;

// Requests, drivers and their reports (08 §6.2).
pub use crate::error::{Error, ErrorKind};
pub use crate::filter::{FilterOptions, WriteIndex, run_filter};
pub use crate::index::run_index;
pub use crate::stats::{
    FStatTest, FStatsReport, FStatsRequest, FstKind, InvariantMode, JointSfsReport, LdReport,
    LdRequest, Levels, PbsKind, Populations, SfsReport, SfsRequest, StatKind, StatsReport,
    StatsRequest, WindowLayoutOptions, assign_groups, default_kinds, run_stats, run_stats_opened,
    window_kinds,
};
pub use reed_stats::ld::LdWindow;

// The output writers (08 §6.2 `reed::writers`): window layout, per-site reports
// and the bcftools-stats summary renderer.
pub use crate::stats::writers::{fstats, ld, reports, sfs, summary, windows};

// The request sub-modules the CLI builds against (08 §6.2): `spec`, `expr`,
// `plan`, `funnel`.
pub use reed_filter::{expr, funnel, plan, spec};

// Input opening (the binary's one-shot opens); hidden from the docs.
#[doc(hidden)]
pub use crate::input::open_regions;

// Front-end re-exports (08 §6.2): the parallel budgets, the region type and
// the BGZF helpers the tests write fixtures with.
pub use reed_core::coord::{Interval0, Pos0};
pub use reed_core::grid::{ContigRegions, WindowGrid};
pub use reed_engine::budget;
pub use reed_engine::metrics;
pub use reed_io::{bgzf, chunk, header};
