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
//! let budget = reed::ThreadBudget::default();
//! let abort = reed::abort_for(&budget);
//! let ctx = reed::RunContext::new(&budget, &abort.signal);
//! let report = reed::run_filter(
//!     Some(std::path::Path::new("in.vcf.gz")),
//!     &reed::spec::FilterSpec::default(),
//!     &reed::FilterOptions::default(),
//!     &ctx,
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
//! let budget = reed::ThreadBudget::default();
//! let abort = reed::abort_for(&budget);
//! let ctx = reed::RunContext::new(&budget, &abort.signal);
//! let report = reed::run_popgen(
//!     Some(std::path::Path::new("in.vcf.gz")),
//!     &request,
//!     &ctx,
//! )?;
//! let tsv = reed::render_tsv(&report.windows, &request.kinds, false);
//! let _ = tsv;
//! # Ok(())
//! # }
//! ```

//! # Stability
//!
//! The stability promise covers the root items and the `run_*` entry points
//! (`docs/spec/api-contract.md`). The `frontend` module is `#[doc(hidden)]`
//! and excluded from semver.

#![allow(unreachable_pub)]

mod bcf_input;
mod context;
mod convert;
mod error;
mod filter;
mod filtered;
mod index;
mod input;
mod rename;
mod spill;
mod stats;

// Requests, drivers and their reports (08 §6.2).
mod dapc;
mod genotype;
mod pca;

/// The curated one-line import set (`use reed::prelude::*;`).
pub mod prelude;

/// The per-run budget and abort signal every driver takes (doc_v30 P2 §2).
pub use crate::context::RunContext;
pub use crate::dapc::{
    DapcReport, DapcRequest, XvalOutcome, XvalRequest, run_dapc, run_dapc_opened,
};
pub use crate::error::{Error, ErrorKind};
pub use crate::filter::{FilterOptions, WriteIndex, run_filter};
pub use crate::index::{IndexOutcome, run_index};
pub use crate::input::{
    Framing, InputCapability, InputMeta, OpenedInput, ReplayableInput, SeekableInput, StreamedInput,
};
pub use crate::pca::{
    MatrixUse, PcaReport, PcaRequest, read_genotype_matrix, run_pca, run_pca_opened,
};
pub use crate::stats::{
    DecayBin, FStatRow, FStatTest, FStatsReport, FStatsRequest, FstKind, InvariantMode,
    JointSfsReport, LdReport, LdRequest, Levels, PbsKind, PopgenReport, PopgenRequest, Populations,
    SfsReport, SfsRequest, StatKind, WindowLayoutOptions, WindowOutputConfig, WindowRow,
    assign_groups, assign_groups_in_memory, default_kinds, run_popgen, run_popgen_opened,
    window_kinds,
};
/// Dense row-major matrix, shared with `reed-reduce` and `reed-graph`.
pub use reed_core::linalg::Mat;
pub use reed_engine::abort::{Abort, AbortHandle, AbortSignal, abort_for};
pub use reed_engine::budget::{MemBudget, ThreadBudget};
// The convert command's request surface (doc_v13 §8): the typed targets, the
// driver entry points and their reports.
pub use crate::convert::{
    ConvertCounters, ConvertReport, ConvertRequest, ConvertTarget, PhyloFormats, PhyloOptions,
    run_convert, run_convert_opened,
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
    StatsOutcome, StatsReportRequest, run_stats, run_stats_opened, validate_report_tags,
};
pub use reed_filter::funnel::FilterReport;
/// The filter request types that appear in `run_filter`'s signature and in
/// the typed request structs (doc_v30 P3: every signature type is nameable).
pub use reed_filter::spec::FilterSpec;
pub use reed_filter::spec::region::SiteSelector;
// The accumulator and report types that appear in the facade's public
// signatures (doc_v7 D7-22: they were `pub` in private modules but never
// re-exported, so a library caller could not name them).
pub use reed_stats::request::{StatAcc, StatOutput};
pub use reed_stats::summary::report::QcReportAcc;
// Budgets and errors the drivers take or return (doc_v30 P3 nameability).
pub use reed_engine::error::EngineError;
pub use reed_engine::metrics::StageMetrics;
// Input-side types in public methods: `format()`, `replay(&[Region])` and
// `into_parts` return them, so they must be nameable (api-contract §2).
pub use reed_io::bcf::BcfDict;
pub use reed_io::chunk::ChunkSource;
pub use reed_io::chunk::region::Region;
pub use reed_io::error::IoError;
pub use reed_io::format::Format;

// The request vocabulary (doc_v30 P3): `FilterSpec` and the rule types a
// caller builds. The build-time plan (`plan`), expression parser (`expr`)
// and funnel internals are not part of the designed API and live in
// `frontend` below.
/// The documented consumer of a windowed-statistics report. The other
/// writer modules are front-end-only (`frontend`); this renderer is what the
/// report's public example uses.
pub use crate::stats::writers::windows::render_tsv;
pub use reed_filter::spec;

// Front-end re-exports (08 §6.2): the parallel budgets, the region type and
pub use reed_core::coord::{Interval0, Pos0};
pub use reed_core::grid::{ContigRegions, WindowGrid};
/// The header's declared tag table (doc_v14 §3.3: `--hist` resolves its
/// `INFO/`/`FMT/` requests against it before the first record).
pub use reed_core::tag_table::{TagNumber, TagTable, TagType};
/// The VCF header as an opaque value (doc_v26 V26-10): the operations a
/// downstream user needs, without a third-party type in the surface.
pub use reed_io::header::{HeaderFacts, VcfHeader};

/// The in-repo front end's escape hatch (doc_v30 P3): everything `reed-cli`
/// and the integration tests need that is not part of the designed library
/// API — the output writers, the build-time plan, the expression parser, the
/// distance/graph/reduce crates, BGZF fixture helpers, the one-shot input
/// openers and the atomic-output primitives.
///
/// `#[doc(hidden)]`: not part of the documented surface and not covered by
/// the semver promise (ADR `api-contract` §2). External callers should use
/// the `run_*` entry points and the request/report types at the root.
#[doc(hidden)]
pub mod frontend {
    pub use crate::input::{PendingInput, open_regions};
    pub use crate::stats::writers::{fstats, histogram, ld, reports, sfs, summary, windows};
    pub use reed_distance as distance;
    pub use reed_engine::{budget, metrics};
    pub use reed_filter::{expr, funnel, plan};
    pub use reed_graph as graph;
    pub use reed_io::header;
    pub use reed_io::header::tag_table;
    pub use reed_io::write::atomic as atomic_file;
    pub use reed_io::write::atomic::write_atomic;
    pub use reed_io::{bgzf, chunk};
    pub use reed_reduce as reduce;
}
