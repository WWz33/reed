# Library use

The `reed` crate (`crates/reed`) is the public API; the command line is a wrapper around it.

```toml
[dependencies]
reed = { path = "crates/reed" }
```

The stability promise covers the crate root only; crates below `crates/` are internal implementation whose types and signatures may change. Errors are `reed::Error`, classified by `Error::kind()`, with the message read through `Display`/`source()`.

## Entry points

The crate root exports six entry points, one per subcommand:

```rust
run_filter, run_popgen, run_stats, run_rename, run_index, run_convert
```

The signature pattern is uniform: `fn run_*(input: Option<&Path>, request: &XxxRequest, budget: &ThreadBudget) -> Result<XxxReport>`. `ThreadBudget` (`reed::budget`) holds the thread count and memory budget; `Default` means "all cores + half of physical memory", the same as the CLI defaults.

## Population statistics

```rust
use reed::{PopgenRequest, Populations, StatKind};

let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::File("pops.txt".into()),
    ..PopgenRequest::default()
};
let report = reed::run_popgen(Some("in.vcf.gz".as_ref()), &request,
                              &reed::budget::ThreadBudget::default())?;
let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false);
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
let report = reed::run_convert(Some("in.vcf.gz".as_ref()), &request,
                               &reed::budget::ThreadBudget::default())?;
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

The writer modules are usable directly: `reed::windows` (window TSV), `reed::reports`, `reed::histogram`, `reed::sfs`, `reed::fstats`, `reed::ld`, `reed::summary`.

## Round-trip rendering

`PopgenRequest::to_args(&group_names)` and `FilterSpec::to_args()` render a request back to CLI spelling, and `parse(render(x)) == x` (locked by tests). Uses: recording a library-built request in a log, or building it in code and handing it to a CLI process.
