---
title: Library use
---

The `reed` crate (`crates/reed`) is the public API; the command line is a wrapper around it.

```toml
[dependencies]
reed = { path = "crates/reed" }
```

The stability promise covers the crate root only; crates below `crates/` are internal implementation whose types and signatures may change. Errors are `reed::Error`, classified by `Error::kind()`, with the message read through `Display`/`source()`.

## Stability

The root items and `run_*` entry points follow semver (0.x: breaking changes need a minor bump); `#[doc(hidden)] reed::frontend` and every crate under `crates/` are internal. See `docs/spec/api-contract.md`.

## Entry points

The crate root exports one `run_*` per subcommand (`_opened` variants take an
already-opened `OpenedInput` and are listed in `ci/api/reed.txt`):

| Entry point | Request | Report |
|---|---|---|
| `run_filter` | `FilterSpec` | `FilterReport` |
| `run_popgen` | `PopgenRequest` | `PopgenReport` |
| `run_stats` | `StatsReportRequest` | `StatsOutcome` |
| `run_rename` | `RenameSpec` | `RenameReport` |
| `run_index` | `WriteIndex` (argument) | `IndexOutcome` |
| `run_convert` | `ConvertRequest` | `ConvertReport` |
| `run_pca` | `PcaRequest` | `PcaReport` |
| `run_dapc` | `DapcRequest` | `DapcReport` |
| `run_distance` | `DistanceRequest` | `DistanceReport` |

The signature pattern is uniform: `fn run_*(input: Option<&Path>, request: &XxxRequest, ctx: &RunContext<'_>) -> Result<XxxReport>`. `RunContext` (`reed::RunContext`) borrows the run's `ThreadBudget` and `AbortSignal`; build one with `RunContext::new(&budget, &abort.signal)` from an `Abort` (`reed::abort_for(&budget)`), or `RunContext::idle(&budget)` when the caller has no cancellation source. `ThreadBudget::default()` means "all cores + half of physical memory", the same as the CLI defaults. The same names are available as `use reed::prelude::*;`; the prelude only re-exports root items (the `run_*` entry points, requests, reports, errors, budgets and `RunContext`).

## Population statistics

```rust
use reed::{PopgenRequest, Populations, StatKind};

let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::File("pops.txt".into()),
    ..PopgenRequest::default()
};
let budget = reed::ThreadBudget::default();
let ctx = reed::RunContext::idle(&budget);
let report = reed::run_popgen(Some("in.vcf.gz".as_ref()), &request, &ctx)?;
let tsv = reed::render_tsv(&report.windows, &request.kinds, false);
```

The population mapping does not have to touch disk:

```rust
let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::Assignment(vec![
        ("sample1".into(), "popA".into()),
        ("sample2".into(), "popB".into()),
    ]),
    ..PopgenRequest::default()
};
```

`StatKind` variants cover every CLI statistic (`Pi`, `Dxy`, `Fst`, `WattersonTheta`, `TajimaD`, `Pbs`, the f-statistic family, `Sfs`, the LD family); requests the CLI cannot spell (per-bp statistics, a `Regions` grid) fail in `to_args()` rendering but work when called directly.

## Format conversion

```rust
use reed::{ConvertRequest, ConvertTarget, PhyloFormats, PhyloOptions};

let request = ConvertRequest {
    targets: vec![ConvertTarget::Phylo {
        formats: PhyloFormats { phylip: true, fasta: true, nexus: false },
        options: PhyloOptions::default(),
    }],
    prefix: "matrices".into(),
    ..ConvertRequest::default()
};
let budget = reed::ThreadBudget::default();
let ctx = reed::RunContext::idle(&budget);
let report = reed::run_convert(Some("in.vcf.gz".as_ref()), &request, &ctx)?;
```

## Requests and reports

| Request type | Report type |
|---|---|
| `PopgenRequest` | `PopgenReport` |
| `FilterSpec` | `FilterReport` |
| `StatsReportRequest` | `StatsOutcome` |
| `RenameSpec` | `RenameReport` |
| `ConvertRequest` | `ConvertReport` |
| (no request structure for index) | `IndexOutcome` |

`render_tsv` renders a windowed-statistics report; the remaining writer modules are front-end-only (behind `#[doc(hidden)] reed::frontend`, not part of the stability promise).

## Round-trip rendering

`PopgenRequest::to_args(&group_names)` and `FilterSpec::to_args()` render a request back to CLI spelling, and `parse(render(x)) == x` (locked by tests). Uses: recording a library-built request in a log, or building it in code and handing it to a CLI process.
