//! The facade's one-line import set (doc_v30 P3).
//!
//! `use reed::prelude::*;` brings in the `run_*` entry points and the
//! vocabulary their requests and reports are built from: request/options
//! types, reports, errors, the run context with its budget/abort types, and
//! the already-opened input. The root re-exports remain the canonical paths;
//! this module only curates them for callers that do not want to list each
//! one.
//!
//! ```
//! use reed::prelude::*;
//!
//! let budget = ThreadBudget::default();
//! let ctx = RunContext::idle(&budget);
//! assert_eq!(ctx.budget.threads, budget.threads);
//! ```

pub use crate::{
    Abort, AbortHandle, AbortSignal, ChrStyle, ConvertCounters, ConvertReport, ConvertRequest,
    ConvertTarget, DapcReport, DapcRequest, DecayBin, DistanceFormat, DistanceKind, DistanceReport,
    DistanceRequest, EngineError, Error, ErrorKind, FStatRow, FStatTest, FStatsReport,
    FStatsRequest, FilterOptions, FilterReport, FilterSpec, FstKind, IndexOutcome, InvariantMode,
    IoError, JointSfsReport, LdReport, LdRequest, Levels, MatrixUse, MemBudget, MissingPolicy,
    OpenedInput, PbsKind, PcaReport, PcaRequest, PcaSolver, PhyloFormats, PhyloOptions,
    PopgenReport, PopgenRequest, Populations, QcReportAcc, RenameOutput, RenameReport, RenameSpec,
    RunContext, SfsReport, SfsRequest, Stand, StatAcc, StatKind, StatOutput, StatsOutcome,
    StatsReportRequest, ThreadBudget, TreeMethod, WindowLayoutOptions, WindowOutputConfig,
    WindowRow, WriteIndex, XvalOutcome, XvalRequest, XvalResultMode, abort_for, run_convert,
    run_convert_opened, run_dapc, run_dapc_opened, run_distance, run_distance_opened, run_filter,
    run_filter_opened, run_index, run_pca, run_pca_opened, run_popgen, run_popgen_opened,
    run_rename, run_rename_opened, run_stats, run_stats_opened,
};
