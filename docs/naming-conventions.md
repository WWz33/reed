# reed Naming Conventions

Owner ruling: one unified architecture, one unified naming scheme across all crates. This document
describes the **target** scheme for the end state of the in-flight migration (record-level transport
moves into reed-io; reed-filter becomes policy-only; reed-core gains `estimators::site_counts`).
Current names that must change are listed only in the rename ledger (§3); nothing here blesses
existing accidents. Authority: `DESIGN.md` §4.5, `docs/filter-architecture.md` §3.

Applies to every `pub` item in reed-core / reed-io / reed-filter, and to doc comments that name
concepts. Published-API status: reed-core 0.2.0 and reed-io 0.2.0 are published surfaces; reed-filter
0.1.0 is still free to rename. §6 extends the scheme from identifiers to file names (Rust modules,
gold scripts) and bioinformatics artifacts (VCF/BCF outputs, index sidecars, logs, reports).

## 1. Concept → canonical name

| Concept | Canonical name | Home | Retired names / anti-names |
|---|---|---|---|
| Transport record, full fidelity, surgically editable | `TransparentRecord` | reed-io | `RecordBuf` (noodles' vendor name) in reed signatures, `buf`, `rec` |
| Decoded per-site/record view (FieldSet-driven) | `RecordView` | reed-io | `DecodedSite` (reed-filter), `PendingSite` (reed-io private) |
| Window materialisation (L2) | `WindowBuffer` | reed-core | — (already canonical) |
| Window grid + manifest declaration | `ReadOptions` | reed-io | — |
| Input format | `Format::{Vcf, VcfGz, Bcf}` + `Format::detect(path)` | reed-io | `SourceKind` (exec.rs), three ad-hoc suffix matchers |
| Record iteration callback | `for_each_record(path, opts, sink)`; dyn alias `RecordSink<'a>` | reed-io | `for_each_buf` (exec.rs), `reader.record_bufs` leakage |
| Window iteration callback | `for_each_window(path, opts, sink)`; dyn alias `WindowSink<'a>` | reed-io | — (alias must appear in a pub signature or become private) |
| Convenience materialisation | `read_records`, `read_windows` | reed-io | — |
| Sample-retention mask (pass0 output) | `SampleMask` (Bits-backed; predicate `keep(i)`) | reed-core | `keep: &[bool]` bare slices, `miss: Vec<u64>` naming |
| Genotype mask result (GT set to `./.`) | field `masked` + verb `mask` | `RecordView` | — (distinct concept from retention; never swap the words) |
| Per-site exact counts (count-first rule) | `estimators::site_counts::SiteCounts` + `site_counts()` | reed-core | `record::SiteCounts` (reed-filter), `DecodedSite::counts` |
| Per-site tallies over a sample subset | `window::AlleleCounts`, `window::GenotypeCounts` | reed-core | — |
| Window-aggregate group tallies | `AlleleCountsResult` | reed-core | watch: name is close to `SiteCounts`; do not add more `*CountsResult` types |
| Pair-level audit trail | `PairCounts` | reed-core | — |
| Filter declaration | `FilterSpec` (+ `CohortRule`, `RegionRule`, `SiteRules`, `InfoRules`, `GtMask`, `Recompute`) | reed-filter | — |
| Expanded plan | `ExecutionPlan` | reed-filter | `Plan` |
| Evaluation passes | `reduce` (= pass0, genome-wide sample missing-rate reduction), `scan` (= pass1, filter + edit + emit), `ordered_carry` (future LD) | reed-filter | `pass0`/`pass1` in identifiers (`read_header_pass1`); the numbers survive only in prose mapping to the gold pipeline |
| Output sink trait | `Emitter` (`write_header`, `write_record`, `finish`) | reed-filter | `write_record_buf`; private enum `Sink`/struct `Writer` naming |
| Emitter variants | `StdoutVcf`, `VcfGzWriter`, `BcfWriter`, `WindowFeeder` (reserved, v0.2) | reed-filter | — |
| Header manipulation | `HeaderInfo`, `header_info()` (reed-io); `retain_samples(header, &SampleMask)`, `ensure_info_defs()`, `STALE_INFO_TAGS` (reed-filter `header` module) | reed-io / reed-filter | `filter_header`, `ensure_recompute_defs`, header fns living under `record::vcf` |
| Errors | `WindowError`, `PopError`, `ReadError`, `FilterError`, `Unavailable` | per crate | variants `User`, `Core`, `Format` (see §3) |
| Input knobs | `*Options` (`ReadOptions`, `ExecOptions`), value types (`Clone + Debug`) | — | — |
| Output reports | `*Report` (`RunReport`, `StreamReport`) | — | `StreamSummary` |
| Contig identity | `ContigId` (stream-local numeric id) in reed-core results; contig **name** (`String`) at the transport/policy boundary (`RecordView.contig`, `Region.contig`) | — | never use bare `contig` for a `ContigId` field; field is `contig` only when it is the name |
| Estimator definition enums | `PiEstimator`, `DxyEstimator`, `FstEstimator`, `PbsEstimator`, … with `id() / citation() / availability()` | reed-core | — (pattern for all future estimator families) |
| Test helpers | `fixture()` = the canonical domain object of the module; `fixture_*` / `make_*` constructors; `fixture_path(name)` for file paths | all crates | `site()`, `test_records()`, `sample_record_buf()`, `fixture(name)` as a path helper |

## 2. Rules

### 2.1 Casing
- Rust API guidelines: types/traits `UpperCamelCase`, functions/methods `snake_case`, consts
  `SCREAMING_SNAKE_CASE`, modules lower-case singular. Acronyms are camel-cased inside type names:
  `VcfGz`, `HweSiteResult`, `LdPair`, `PbsResult` — never `VCFGZ`/`HWE`.
- Counts are `n_<unit>` (`n_sites`, `n_samples`, `n_masked_gt`); booleans are `is_/has_/require_`
  (`is_variant`, `is_pass`, `require_pass`); thresholds are `min_/max_` (`min_maf`, `dp_max`).

### 2.2 Trait vs struct suffixes
- A trait is a **role noun** with no suffix: `Emitter`. A concrete type is named by format + role:
  `VcfGzWriter`, `BcfWriter`. Never both a trait and an impl named `Writer`.
- **`Source` is banned for I/O types.** A type that reads is `*Reader`; a free function that streams
  is `for_each_*`. (`VcfSource`/`BcfSource` are retired with the transport migration.)
- `Sink` is reserved for the reader-side callback aliases (`RecordSink`, `WindowSink`); the emitter's
  private dispatch enum must not be called `Sink`.
- `open(path, …)` is the fallible path constructor on both sides (`Format::open` for input,
  `emitter::open` for output).

### 2.3 Verb vocabulary (one verb per action)
| Verb | Meaning | Canonical bearers |
|---|---|---|
| `read` | pull bytes/records/headers from a path or reader | `read_header`, `read_records`, `read_windows`, `header_info` |
| `decode` | bytes → typed view (FieldSet-driven) | `RecordView` construction, `decode_gt` (private) |
| `stream`/`for_each` | push every item to a callback in order | `for_each_record`, `for_each_window` |
| `write` | serialize an item | `Emitter::write_header`, `Emitter::write_record` |
| `emit` | the sink side of write (trait-level role) | `Emitter` |
| `keep` | retention predicate over samples | `SampleMask::keep(i)`, `CohortRule` |
| `prune` | drop sample columns / header names from output | `retain_samples`, sample-column pruning in edits |
| `mask` | set a genotype to missing (`./.`) | `mask`, `GtMask`, `masked` |
| `recompute` | rewrite derived INFO tags | `Recompute`, `ensure_info_defs` |
| `detect` | sniff input format from extension/magic | `Format::detect` |
| `compile` | spec → plan expansion | `FilterSpec::compile`, `plan::compile` |
| `retain` | reserved for `Vec::retain` semantics; never a public API verb here | — |

### 2.4 One concept, one name; one name, one concept
- A concept may not change names across crates (`keep` means sample retention everywhere; `mask`
  means set-to-missing everywhere).
- A name may not mean two things: today `contig` is a `ContigId` in reed-core and a `String` in
  reed-filter, and `pos` changes base between layers (see 2.6). Both are allowed only under the
  explicit domain rules in §1/§2.6, and every such field carries a doc comment stating which.

### 2.5 No vendor vocabulary in public APIs
- noodles' `RecordBuf` is an implementation detail; reed's public surface names the reed concept
  (`TransparentRecord`, `RecordView`). Vendor method names (`record_bufs`, `write_variant_record`)
  do not appear in reed signatures or docs.

### 2.6 Coordinates
- reed's own coordinate space is **0-based, half-open** (`WindowBuffer.start/end/pos`,
  `*SiteResult.pos`, window grids).
- **VCF-native record space is 1-based** (`RecordView.pos` = VCF POS; `Region` is 1-based inclusive,
  mirroring `Chr01:100-200` text). Conversion happens exactly once per direction, in one helper each
  (`pos0()` / `pos1()` in reed-io), never inline arithmetic at call sites.
- Error messages always report 1-based VCF POS. Every `pos`/`start`/`end` field is documented with
  its base. (Verified: noodles-core 0.21 `Position` is 1-based; reed-io and reed-filter currently
  convert with two different ad-hoc idioms.)

### 2.7 Errors
- One error enum per crate, `*Error` suffixed; variants name the **failure mode**, never an actor or
  a crate (`User`, `Core` are banned). Structured payloads over strings where location exists
  (`Decode { contig, pos, msg }` pattern). `From` chains must keep the underlying error reachable
  via `Error::source`.
- `Unavailable` (reed-core) stays as-is: it is a declared refusal, not an I/O failure.

### 2.8 Options / Reports
- Inputs: `*Options`, value types, no I/O side effects at construction. Outputs: `*Report` (a run
  with timing/funnel) — `StreamSummary` becomes `StreamReport`. Report counters use `n_*`.

### 2.9 Modules
- Module names singular; a function must not duplicate its module name (`hwe::hwe_exact`,
  `pi::pi`, `fst::fst` stutter). Tolerated in reed-core 0.2.0; the fix is re-exporting functions at
  `estimators::` level when the published API next breaks — not renaming call sites.

### 2.10 Abbreviations
- Allowed (domain canon from VCF/popgen, also the INFO/FORMAT tag names): `GT`, `AC`, `AN`, `NS`,
  `MAF`, `HWE`, `LD`, `SFS`, `Fst`, `PBS`, `Dxy`, `VCF`, `BCF`, `BGZF`, `DP`, `GQ`, `QUAL`, `INFO`.
- Banned in public APIs, doc comments naming concepts, and type/field/function names: `buf`, `rec`,
  `iv`, `idx`, `src`, `tmp`, `cnt`. Fine as locals.

### 2.11 Tests
- Every test module has at most one `fixture()` (the canonical domain object); parameterized
  constructors are `fixture_*`/`make_*`; file-path helpers are `fixture_path(name)`. Test names read
  as behaviour sentences (`dot_sample_column_decodes_as_missing`).

## 3. Rename ledger

Risk key: **internal** = private item, free; **rf-0.1** = reed-filter published surface but
pre-release, rename now; **rio-0.2** = breaks reed-io 0.2.0; **rc-0.2** = breaks reed-core 0.2.0.

| Current | Target | Location | Reason | Risk |
|---|---|---|---|---|
| `VcfSource`, `BcfSource` | deleted (transport moves to reed-io); interim rename `VcfReader`/`BcfReader` | reed-filter/src/exec.rs | `Source` banned for I/O; I/O semantics leave reed-filter | internal |
| `SourceKind` | `reed_io::Format` | reed-filter/src/exec.rs | third copy of format dispatch | internal |
| `for_each_buf` | `for_each_record` (reed-io) | reed-filter/src/exec.rs | `buf` anti-name; verb vocabulary | internal |
| `DecodedSite` | `reed_io::RecordView` | reed-filter/src/record/mod.rs | one decoded-view type, owned by transport | rf-0.1 |
| `record::SiteCounts` (+ `DecodedSite::counts`) | `reed_core::estimators::site_counts::SiteCounts` + `site_counts()` | reed-filter/src/record/mod.rs → reed-core | count-first rule belongs to core; reused by estimators | rf-0.1 (move), additive to rc |
| `keep: &[bool]` (in `mask`, `counts`, `apply`, `filter_header`) | `&SampleMask` | reed-filter/src/record/*.rs | named type for the retention mask; no bare bool-slice signatures | rf-0.1 |
| `record::vcf::filter_header` | `header::retain_samples(header, &SampleMask)` | reed-filter/src/record/vcf.rs | header fns out of the format module; `prune`/`retain` vocabulary | rf-0.1 |
| `record::vcf::apply` | `record::apply_edits` | reed-filter/src/record/vcf.rs | `apply` is uninformative; edits vocabulary | rf-0.1 |
| `record::vcf::count_missing` | `record::count_sample_missing` | reed-filter/src/record/vcf.rs | it is the pass0 per-sample reduction input, not a site count | rf-0.1 |
| `record::vcf::decode` | reed-io FieldSet decode; filter-side helper dies | reed-filter/src/record/vcf.rs | decode is transport | rf-0.1 |
| `STALE_INFO` | `STALE_INFO_TAGS` (stays with policy) | reed-filter/src/record/vcf.rs | name the thing (tags), live in the policy module | rf-0.1 |
| `Emitter::write_record_buf` | `Emitter::write_record` | reed-filter/src/emitter.rs | `buf` anti-name; parameter becomes `TransparentRecord`/`&RecordView` | rf-0.1 |
| private `enum Sink`, `struct Writer` | `enum Backend`, `struct EmitterImpl` (or any non-reserved names) | reed-filter/src/emitter.rs | `Sink`/`Writer` are reserved vocabulary (§2.2) | internal |
| `Plan` | `ExecutionPlan` | reed-filter/src/plan.rs | matches filter-architecture.md §3; `Plan` is over-generic | rf-0.1 |
| `read_header_pass1` | `read_header_for_scan` | reed-filter/src/exec.rs | pass numbers out of identifiers (§1 passes row) | internal |
| `ensure_recompute_defs` | `ensure_info_defs` | reed-filter/src/exec.rs | it adds INFO definitions, not "recompute" | internal |
| `FilterError::Format` | `FilterError::UnsupportedInput` | reed-filter/src/lib.rs | "Format" reads as formatting; name the failure mode | rf-0.1 |
| `RunReport.funnel: Vec<(String, u64)>` | `Vec<FunnelStage { rule: RuleId, passed: u64 }>` | reed-filter/src/exec.rs | stringly-typed rule names in a public report | rf-0.1 (defer, additive) |
| `site()`, `test_records()`, `sample_record_buf()`, `fixture(name)` | `fixture_site()`, `fixture_records()`, `fixture_record()`, `fixture_path(name)` | reed-filter tests, reed-io tests | §2.11 | internal |
| `ReadError::User` | split: `ReadError::Options(String)` (validation) + `ReadError::Callback(String)` (sink errors) | reed-io/src/lib.rs (+ examples/scan.rs) | variant named after an actor; two failure modes conflated | rio-0.2 |
| `ReadError::Core` | `ReadError::Window` | reed-io/src/lib.rs | variant named after a crate, not the payload concept | rio-0.2 |
| `StreamSummary` | `StreamReport` | reed-io/src/lib.rs | `*Report` rule (§2.8) | rio-0.2 |
| `WindowSink` (pub, used by no pub fn) | used in `for_each_window`'s dyn overload or made private | reed-io/src/lib.rs | pub item that is not part of any pub signature | rio-0.2 (signature change) |
| (new) `Format`, `RecordView`, `TransparentRecord`, `for_each_record`, `RecordSink`, `header::retain_samples` | — | reed-io | transport migration targets | additive, rio-0.3 |
| `WattersonResult.value` naming, `theta_per_bp`/`total_theta` | keep | reed-core/src/estimators/spectrum.rs | already conform | — |
| `AlleleCountsResult` | keep, but freeze the `*CountsResult` slot | reed-core/src/estimators/allele_counts.rs | too close to the new `SiteCounts`; rename only if reed-core breaks anyway | rc-0.2 (deferred) |
| estimator fns re-exported at `estimators::` (kills `pi::pi`, `fst::fst`, `hwe::hwe_exact` stutter) | do at next breaking reed-core release | reed-core/src/estimators/* | §2.9 | rc-0.2 (deferred) |
| `d_stat` fn vs `DStatResult` | keep (acronym rule covers `DStat`); consider `DResult` only at the breaking release | reed-core/src/estimators/fstats.rs | cosmetic | rc-0.2 (deferred) |

Deliberately **not** renamed (already canonical): `FilterSpec` + rule types, `GtMask` (GT is the VCF
column name; conf keys `GT_*` depend on it), `Emitter` trait name, `WindowBuffer`, `ContigId`,
`Bits`, `PopulationView`, `SiteCompleteness`/`Ascertainment`/`Availability`, `PairCounts`,
`*Result` types, `*Estimator` enums, `ExecOptions`/`ReadOptions`, `WindowError`/`PopError`/
`FilterError`, `Unavailable`.

## 4. Migration order

| Phase | Renames | Why now / constraint |
|---|---|---|
| **0 — safe today** (no published API touched) | All **internal** rows: `VcfSource`/`BcfSource` → `*Reader` (or leave for deletion), `read_header_pass1`, `ensure_recompute_defs`, private emitter `Sink`/`Writer`, test helpers; all **rf-0.1** rows that stay inside reed-filter (`FilterError::Format`, `Plan` → `ExecutionPlan`, `apply` → `apply_edits`, `count_missing` → `count_sample_missing`, `STALE_INFO` → `STALE_INFO_TAGS`) | reed-filter 0.1.0 surface is free; parallel refactors only touch these files. Do **not** start `keep → SampleMask` or `DecodedSite` renames in reed-filter while subagents are mid-edit — sequence them with phase 1 |
| **1 — with the record-transport migration** (reed-io 0.3.0 + reed-filter 0.2.0, one breaking batch) | `for_each_record` + `RecordSink`; `RecordView`/`TransparentRecord` land and `DecodedSite`/`record` module retire; `Format::detect` replaces all suffix matching (fixes the `.bgz` divergence); `SampleMask` introduced and all `keep: &[bool]` signatures flip; `SiteCounts` moves to reed-core `estimators::site_counts` (additive); `Emitter::write_record_buf` → `write_record`; `retain_samples`/`header` module; reed-io rows: `ReadError::User`/`Core`, `StreamSummary` → `StreamReport`, `WindowSink` resolution | These all share the transport seam — splitting them across releases would churn the same signatures twice. reed-core gets only additive changes here (0.3.0 minor) |
| **2 — deferred breaking** (next reed-core breaking release, 0.4.0) | estimator fns re-exported at `estimators::`; `AlleleCountsResult` rename decision; `DStatResult` cosmetic; `RunReport.funnel` typing if not done in phase 1 | reed-core 0.2.0 is published; these are cosmetic-or-breaking with no functional pull |

Conflict-risk notes for parallel agents: (a) `record::vcf.rs` is a frozen-contract file owned by
another subagent — its renames (decode/apply/count_missing/filter_header) must land with the
transport migration, not ad hoc; (b) `SampleMask` touches every filter signature at once — one
commit, no gradual migration; (c) do not introduce `RecordView`/`Format` names in reed-filter before
reed-io owns them, or the migration will fork the vocabulary.

## 5. Reserved names (do not squat)

`WindowFeeder` (emitter inline mode), `ordered_carry` (LD pass), `ExecutionPlan`, `HaplotypeStream`,
`MatrixAccumulator`, `DetectedFacts`, `is_accessible` (v0.3 site mask), `RecordSink`, `SampleMask`,
`site_counts`.

## 6. File and artifact naming

Extends §1–§5 from identifiers to file names and bioinformatics artifacts. Status tags:

| Tag | Meaning |
|---|---|
| **APPLIED** | rule in force now; renamed or codified in this change set (2026-10-05) |
| **PROPOSED** | agreed target, blocked on in-flight work — do not execute yet |
| **FROZEN** | historical oracle record — never rename; paths are pinned by existing reports |
| **KEEP** | audited 2026-10-05, already conforms |

### 6.1 Rust module files (APPLIED as rules; per-file verdicts below)

- One concept per file. The file is named for its module; the module for its concept
  (§2.1: lower-case, singular). A file whose name needs "and" is two files.
- File names track module names, not type names: `plan.rs` holds `ExecutionPlan`, `spec.rs` holds
  `FilterSpec`. Renaming the file to match the type would only duplicate it (file-level form of the
  §2.9 stutter rule).
- Use a directory module when a concept has per-format or per-pass backends: `<concept>/mod.rs`
  plus one file per backend (`record/{vcf,bcf}.rs` is the model).
- `lib.rs` (crate root) and `main.rs` (bin target) are Rust canon, exempt.
- Integration tests: `tests/<behaviour-or-contract>.rs` at crate level
  (`reed-core/tests/{contracts,core,families,unequal_groups}.rs` conform).

| File | Primary content | Verdict |
|---|---|---|
| reed-filter/src/spec.rs | `FilterSpec` + builder + `from_conf` | KEEP |
| reed-filter/src/plan.rs | plan compile + ascertainment + provenance hash | KEEP (split `provenance.rs` out when provenance outgrows the plan — PROPOSED) |
| reed-filter/src/exec.rs | engine (reduce + scan passes in one file) | PROPOSED → `exec/{mod,scan,reduce,ld}.rs` per filter-architecture.md §3 — **WAIT**: exec.rs is under a parallel edit; split only after it lands |
| reed-filter/src/emitter.rs | output-sink re-export (reed-io) | KEEP |
| reed-filter/src/main.rs | CLI | KEEP |
| reed-io/src/record/mod.rs | record transport (`for_each_record`, `RecordView`, surgical edits) | KEEP |
| reed-io/src/record/vcf.rs, record/bcf.rs | per-format codecs | KEEP |
| reed-io/src/sink.rs | output sinks (`open` dispatch on output path) | KEEP |
| reed-core/src/{bits,window,population,result,manifest}.rs | one concept each | KEEP |
| reed-core/src/estimators/*.rs | one estimator family per file; `mod.rs` re-export hub | KEEP (`site_counts.rs` already present) |
| reed-filter/src/header.rs | header policy (`retain_samples`, `STALE_INFO_TAGS`; §1 header row) | PROPOSED — **WAIT**: land with the transport migration, not while `record::vcf` still owns the functions |

(filter-architecture.md §3's `reed-filter/src/record.rs` row landed as reed-io `record/` instead;
the `exec/` split is the remaining delta in that tree.)

### 6.2 gold/ scripts and config (APPLIED)

- Action scripts are `<verb>_<object>.<sh|py>`; one branch script per executor
  (`run_<executor>.sh`). `lib.sh` and config are role files, exempt from the verb rule.
- Config: `config/filters.conf` — singular domain noun, `KEY=VALUE` with SCREAMING_SNAKE_CASE keys
  (file-level mirror of §2.1 consts); sourced by `lib.sh`, parsed by `compare_tools.py` /
  `check_parity.py`.

| File | Verb + object | Verdict |
|---|---|---|
| run_bcftools.sh / run_vcftools.sh / run_plink2.sh | run + executor | KEEP |
| make_subset.sh | make + subset | KEEP |
| profile_input.py | profile + input | KEEP (audited: conforms — it profiles the *input*; `profile.py` would stutter with its own `profile.{md,json}` output) |
| compare_tools.py | compare + tools | KEEP |
| check_parity.py (was parity.py) | check + parity | **APPLIED** — the only bare-noun deviator; renamed 2026-10-05, docstring + references updated in gold/README.md and docs/filter-architecture.md |
| lib.sh | role file | KEEP (exempt) |
| config/filters.conf | config | KEEP |

### 6.3 Bioinformatics artifacts (APPLIED as rules)

General grammar for every artifact reed tooling mints outside the frozen parity layout (§6.4):

```
<base>.<stage>[.<qualifier>][.<executor>].<format>
```

| Slot | Values | Rule |
|---|---|---|
| base | dataset lineage = input stem minus format (`clean` from `clean.vcf.gz`) | no dots inside a token (use `_`) |
| stage | `subset` / `raw` / `common` / `filtered` / `gold` | ordered pipeline position |
| qualifier | region/scope tokens (`Chr01-03`) | data-derived names verbatim — `Chr01`, never `chr01` |
| executor | `bcftools` / `vcftools` / `plink2` / `reed` | omitted for the oracle and when unambiguous |
| format | `vcf.gz` / `bcf` / `txt` / `tsv` / `md` / `json` | — |

- Tokens are dot-separated, lower-case, snake inside a token (§2.1 casing, file level).
- Index sidecars: same stem + `.tbi` (tabix) or `.csi` (CSI). Never rename an artifact without its
  sidecar; never orphan an index.
- Logs: the log of artifact `X.<format>` is `X.log` (stem minus format suffix).
- Reports: three-way concordance = `concordance.{txt,tsv,md}`; the parity-gate report for candidate
  `Y` is `Y.parity.{txt,json}` (via `check_parity.py --report`).
- Population/config files: `<role>.<action>.txt` (`samples.keep.txt`).
- Scratch/temp: `work/` only, never `out/`. Free-form but producer-prefixed snake_case
  (`bcftools.tagged.bcf`, `vcftools.raw.vcf`, `sample_stats.psc`); exit-code markers are
  `<script-stem>.rc` (`run_bcftools.rc`). `work/` is ephemeral and exempt from the grammar; `out/`
  holds only finished artifacts.
- Run dirs: `out/<scope>/` with scope `full` or `subset.<ContigRange>` (`out/full/`,
  `out/subset.Chr01-03/`). Parity candidates (outputs under test, e.g. reed's) sit directly under
  `out/`, beside their logs.

Examples:
- `clean.filtered.reed.vcf.gz` + `clean.filtered.reed.vcf.gz.tbi` + `clean.filtered.reed.log` —
  reed's full-data filter output (APPLIED 2026-10-05; was `reed-full.{vcf.gz,log}`).
- `clean.filtered.Chr01-03.reed.vcf.gz` — region-restricted run (naming only; not yet produced).
- `gold.vcf.gz` is the oracle's fixed name — never mint `clean.filtered.gold.*`.

### 6.4 Frozen parity layout (FROZEN)

The `run_all.sh` product family inside a run dir is **frozen**: the records under `out/full/` and
`out/subset.Chr01-03/` are oracle evidence (three-way bit-identical agreement, 1,768,759 sites)
referenced by existing reports, and `run_*.sh` must keep emitting these exact names so a re-run
reproduces the layout. The executor-first spelling is a grandfathered exception to §6.3, not a
second convention.

| Artifact | Meaning |
|---|---|
| `gold.vcf.gz` (+`.tbi`) | the oracle (= bcftools full layer) |
| `<executor>.common.vcf.gz` (+`.tbi`) | common layer, three-way comparable |
| `<executor>.full.vcf.gz` (+`.tbi`) | executor full layer (bcftools = gold source) |
| `<executor>.stats.json` | site/sample counts + spec snapshot |
| `samples.keep.txt` | shared sample retention list |
| `concordance.{txt,tsv,md}` | three-way concordance report |
| `profile.{md,json}` (out/ root) | input profile (`profile_input.py --outdir`) |
| `out/logs/<branch>.err`, `out/logs/run_<branch>.out` | per-branch stderr/stdout |
| `out/full_run.log` | historical run console log (frozen) |

### 6.5 File rename ledger (this change set)

| From | To | Status | Why |
|---|---|---|---|
| gold/out/reed-full.vcf.gz | gold/out/clean.filtered.reed.vcf.gz | APPLIED | §6.3 grammar; run had completed (final funnel line in the log) and no process held the path; zero textual references existed to update |
| gold/out/reed-full.vcf.gz.tbi | gold/out/clean.filtered.reed.vcf.gz.tbi | APPLIED | sidecar follows its artifact (§6.3) |
| gold/out/reed-full.log | gold/out/clean.filtered.reed.log | APPLIED | §6.3 log rule |
| gold/scripts/parity.py | gold/scripts/check_parity.py | APPLIED | §6.2 verb_object rule; self-docstring + gold/README.md + docs/filter-architecture.md updated |
| reed-filter/src/exec.rs | reed-filter/src/exec/{mod,scan,reduce,ld}.rs | PROPOSED — **WAIT** | three pass concepts in one file (filter-architecture.md §3); blocked by the parallel exec.rs edit |
| (new) reed-filter/src/provenance.rs | — | PROPOSED | split from plan.rs when provenance outgrows the plan |
| (new) reed-filter/src/header.rs | — | PROPOSED — **WAIT** | §1 header row; land with the transport migration |

Verified 2026-10-05: `bash -n gold/scripts/*.sh` and `python3 -m py_compile gold/scripts/*.py` clean;
`grep -rn "reed-full"` and `grep -rn "parity.py"` clean over gold/ + docs/ outside the rename notes.
