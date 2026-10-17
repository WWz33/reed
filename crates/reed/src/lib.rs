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
//! use reed::{StatKind, PopgenRequest};
//! let request = PopgenRequest {
//!     kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
//!     ..PopgenRequest::default()
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
//!
//! A statistics run borrows the request, so the caller keeps it afterwards
//! (the README's example, kept compiling here — doc_v11 A11-2):
//!
//! ```no_run
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! use reed::{Populations, StatKind, PopgenRequest};
//! let request = PopgenRequest {
//!     kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
//!     populations: Populations::All,
//!     ..PopgenRequest::default()
//! };
//! let report = reed::run_popgen(
//!     Some(std::path::Path::new("in.vcf.gz")),
//!     &request,
//!     &reed::budget::ThreadBudget::default(),
//! )?;
//! let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false);
//! let _ = tsv;
//! # Ok(())
//! # }
//! ```

#![allow(unreachable_pub)]

mod bcf_input;
mod convert;
mod error;
mod filter;
mod filtered;
mod index;
mod input;
mod rename;
mod stats;

// Requests, drivers and their reports (08 §6.2).
pub use crate::error::{Error, ErrorKind};
pub use crate::filter::{FilterOptions, WriteIndex, run_filter};
pub use crate::index::{IndexOutcome, run_index};
pub use crate::input::OpenedInput;
pub use crate::stats::{
    DecayBin, FStatRow, FStatTest, FStatsReport, FStatsRequest, FstKind, InvariantMode,
    JointSfsReport, LdReport, LdRequest, Levels, PbsKind, Populations, SfsReport, SfsRequest,
    StatKind, PopgenReport, PopgenRequest, WindowLayoutOptions, WindowRow, assign_groups,
    assign_groups_in_memory, default_kinds, run_popgen, run_popgen_opened, window_kinds,
};
// The convert command's request surface (doc_v13 §8): the typed targets, the
// driver entry points and their reports.
pub use crate::convert::{
    ConvertCounters, ConvertReport, ConvertRequest, ConvertTarget, MsaRegion, PhyloFormats,
    PhyloOptions, run_convert, run_convert_opened,
};
// The rename command's request surface (doc_v14 §4): the typed spec, the
// output target and the driver.
pub use crate::rename::{
    ChrStyle, RenameOutput, RenameReport, RenameSpec, run_rename, run_rename_opened,
};
// The report and driver of the per-site/per-sample report family (doc_v11
// N6): `run_filter`/`run_popgen`/`run_index` are already at the root, so
// `run_stats` and the types in its signature belong beside them rather
// than only under the writers.
pub use crate::stats::writers::reports::{
    StatsReportRequest, StatsOutcome, run_stats, run_stats_opened,
};
pub use reed_filter::funnel::FilterReport;
// The accumulator and report types that appear in the facade's public
// signatures (doc_v7 D7-22: they were `pub` in private modules but never
// re-exported, so a library caller could not name them).
pub use reed_stats::request::StatAcc;
pub use reed_stats::summary::report::QcReportAcc;

// The output writers (08 §6.2 `reed::writers`): window layout, per-site reports
// and the bcftools-stats summary renderer.
pub use crate::stats::writers::{fstats, histogram, ld, reports, sfs, summary, windows};

// The request sub-modules the CLI builds against (08 §6.2): `spec`, `expr`,
// `plan`, `funnel`.
pub use reed_filter::{expr, funnel, plan, spec};

// Input opening (the binary's one-shot opens); hidden from the docs.
#[doc(hidden)]
pub use crate::input::open_regions;

// Front-end re-exports (08 §6.2): the parallel budgets, the region type and
// the BGZF chunk helpers the front end and its integration tests seek with
// (`bgzf` writes the test fixtures; `chunk` is the seek/region API).
/// The header's declared tag table (doc_v14 §3.3: `--hist` resolves its
/// `INFO/`/`FMT/` requests against it before the first record).
pub use reed_core::tag_table::{TagNumber, TagTable, TagType};
/// Builds the tag table from an opened input's header.
pub use reed_io::header::tag_table;
pub use reed_core::coord::{Interval0, Pos0};
pub use reed_core::grid::{ContigRegions, WindowGrid};
pub use reed_engine::budget;
pub use reed_engine::metrics;
pub use reed_io::{bgzf, chunk};
