//! Nameability gate (doc_v30 P3, `docs/spec/api-contract.md` §2): every type
//! that appears in a `reed` public signature must be nameable from the
//! facade. The lower crates are internal, so a `reed_core::…`/`reed_io::…`
//! path in a signature is a contract leak even though it compiles inside the
//! workspace.
//!
//! `ci/check-nameable.sh` keeps the list complete: it extracts the terminal
//! type names from the lower-crate paths in `ci/api/reed.txt` and fails when
//! one is missing here. The imports below are the compile half of the gate —
//! if an export disappears, this file stops compiling.

#![allow(unused_imports)]

use reed::{
    Abort, AbortHandle, AbortSignal, BcfDict, ChunkSource, EngineError, FilterReport, FilterSpec,
    Format, FstKind, HeaderFacts, IoError, Mat, MemBudget, PbsKind, QcReportAcc, Region,
    RunContext, SiteSelector, StageMetrics, StatAcc, StatKind, StatOutput, TagTable, ThreadBudget,
    WindowGrid, WriteIndex,
};

/// Referencing every imported type in one signature keeps the list honest:
/// an import that is not exported (or a type that stops existing) fails the
/// build. The tuple is never constructed.
#[allow(dead_code)]
type Nameable<'a> = (
    Abort,
    AbortHandle,
    AbortSignal,
    BcfDict,
    Box<dyn ChunkSource>,
    EngineError,
    FilterReport,
    FilterSpec,
    Format,
    FstKind,
    HeaderFacts,
    IoError,
    Mat,
    MemBudget,
    PbsKind,
    QcReportAcc,
    Region,
    RunContext<'a>,
    SiteSelector,
    StageMetrics,
    StatAcc,
    StatKind,
    StatOutput,
    TagTable,
    ThreadBudget,
    WindowGrid,
    WriteIndex,
);
