#!/usr/bin/env bash
# Structural gates that keep the doc_v25 sinks from drifting (00 P3/P4):
# a semantic that was unified must not grow a second implementation.
#
# Run from the repo root (ci/check-structure.sh); exits 1 on the first drift.
set -uo pipefail
cd "$(dirname "$0")/.."

fail=0
check() {
    local name="$1"
    shift
    if "$@"; then
        printf 'ok   %s\n' "$name"
    else
        printf 'DRIFT %s\n' "$name" >&2
        fail=1
    fi
}

# V25-06: one budget constructor. `available_parallelism` may only appear in
# the engine's budget (ThreadBudget) and abort (channel capacity) — the CLI
# resolves through ThreadBudget::from_cli.
check "the CLI has no available_parallelism fallback of its own (doc_v25 V25-06)" \
    bash -c '! grep -rq "available_parallelism" crates/reed-cli/src --include=*.rs'

# V25-06: `-O` is parsed once. The level/letter grammar lives in
# reed-cli/src/params.rs; no translate file may re-derive it.
check "-O parsing lives only in reed-cli/src/params.rs (doc_v25 V25-06)" \
    bash -c '[ "$(grep -rn "invalid -O level" crates/reed-cli/src --include=*.rs | grep -v "src/params.rs" | wc -l)" -eq 0 ]'

# V25-06: `--ld-decay` is parsed by its clap value_parser; the translate
# layers never re-parse the raw text.
check "--ld-decay is not re-parsed outside its value_parser (doc_v25 V25-06)" \
    bash -c '! grep -rn "ld_decay" crates/reed-cli/src/translate --include=*.rs | grep -q "parse_size"'

# V25-01/V25-02: GT code has one home for syntax and one for position.
check "GT syntax helpers live only in reed-io/src/bcf/gt_decode.rs (doc_v25 V25-02)" \
    bash -c '[ "$(grep -rn "fn gt_codes_into\|fn gt_shape" crates --include=*.rs | grep -v "gt_decode.rs" | wc -l)" -eq 0 ]'

# V25-10: the dead linalg copy stays deleted and the env back door is gone.
check "the dead reed-reduce linalg file is gone (doc_v25 V25-10, R12)" \
    bash -c 'test ! -e crates/reed-reduce/src/linalg.rs && ! grep -rq "mod linalg\\|crate::linalg" crates/reed-reduce/src --include=*.rs'

check "no STREAM_DEBUG env back door remains (doc_v25 V25-10)" \
    bash -c '! grep -rq "STREAM_DEBUG" crates --include=*.rs'

# V25-08b: one StageMetrics renderer; the CLI no longer spells the keys.
check "StageMetrics has the one timing renderer (doc_v25 V25-08b)" \
    bash -c 'grep -q "fn render_tsv" crates/reed-engine/src/metrics.rs && ! grep -rq "engine.reader_ms\\|engine.reduce_ms\\|engine.chain_wait_ms" crates --include=*.rs'

# V25-08c: provenance lives in core and argv is captured (and quoted) once.
check "provenance lives in reed-core (doc_v25 V25-08c)" \
    bash -c 'grep -q "pub struct Provenance" crates/reed-core/src/provenance.rs && test ! -e crates/reed-filter/src/plan/provenance.rs'

check "argv is captured once, shell-quoted (doc_v25 V25-08c)" \
    bash -c '[ "$(grep -rn "std::env::args()" crates/reed-cli/src --include=*.rs | grep -v "src/params.rs" | wc -l)" -eq 0 ]'

# doc_v25 08 ④: compile() is four private section modules, not one function.
check "compile() is split into section modules (doc_v25 08 ④)" \
    bash -c 'test -e crates/reed-filter/src/plan/compile/head.rs && test -e crates/reed-filter/src/plan/compile/mask.rs && test -e crates/reed-filter/src/plan/compile/tail.rs && test -e crates/reed-filter/src/plan/compile/finish.rs && ! grep -q "pub fn compile(" crates/reed-filter/src/plan/mod.rs'

# doc_v25 08 ③: convert's spill blobs live in their own module.
check "convert spill types live in convert/spill.rs (doc_v25 08 ③)" \
    bash -c 'grep -q "enum RowsBlob" crates/reed/src/convert/spill.rs && ! grep -rq "enum RowsBlob" crates/reed/src/convert/job'

# V25-09: every on-disk spill goes through the exclusive TempSpill; no
# private name builder may grow back.
check "spills use the shared TempSpill (doc_v25 V25-09)" \
    bash -c 'grep -q "TempSpill::create" crates/reed/src/input.rs && grep -q "TempSpill::create" crates/reed/src/convert/job/mod.rs && ! grep -rq "fn spill_path()" crates/reed/src --include=*.rs'

# V25-08: the CLI depends on the facade only; distance/reduce/graph are
# reachable through reed's re-exports.
check "reed-cli depends on the facade, not the model crates (doc_v25 V25-08)" \
    bash -c '! grep -qE "^reed-(distance|reduce|graph)" crates/reed-cli/Cargo.toml'

# V25-07a: one counting implementation; the frozen oracle keeps its own
# parity reference on purpose.
check "counting semantics live in reed-core only (doc_v25 V25-07a)" \
    bash -c '[ "$(grep -rn "fn count_site\b" crates --include=*.rs | grep -v "/oracle/" | wc -l)" -eq 0 ]'

# V25-07b: no implicit short-mask default survives.
check "no short-keep default remains (doc_v25 V25-07b)" \
    bash -c '! grep -rEn "keep\.get\([^)]*\)\.copied\(\)\.unwrap_or\(true\)" crates --include=*.rs'

# V25-07c: the ghost trigger is gone and the trigger bits come from the
# canonical field table.
check "the ghost haploid-as-missing trigger is gone (doc_v25 V25-07c)" \
    bash -c '! grep -rq "n_haploid_as_missing\|TRIG_HAPLOID_MISSING" crates --include=*.rs'

# doc_v30 R12: the filter layer is L3 over reed-core. Record-text semantics
# (`LineRef`, `InfoIndex`, INFO fragments, the POS/QUAL grammar) live in L0,
# and the `reed-filter -> reed-io` edge must not grow back.
check "reed-filter depends on reed-core only, not reed-io (doc_v30 R12)" \
    bash -c '! grep -q "^reed-io" crates/reed-filter/Cargo.toml && ! grep -rq "reed_io" crates/reed-filter/src --include=*.rs'

check "record-text semantics live in reed-core (doc_v30 R12)" \
    bash -c 'grep -q "pub struct LineRef" crates/reed-core/src/vcf/line.rs && grep -q "pub struct InfoIndex" crates/reed-core/src/vcf/info.rs && test ! -e crates/reed-io/src/vcf/line.rs && test ! -e crates/reed-io/src/vcf/info.rs && test ! -e crates/reed-io/src/vcf/line_edit.rs && test ! -e crates/reed-io/src/vcf/format_layout.rs'

# doc_v30 R05–R07: input capabilities are types, the source has one way
# out, and BCF is asked of `InputMeta` instead of probing an `Option`.
check "input capabilities are typed, not path probes (doc_v30 R05-R07)" \
    bash -c 'grep -q "pub enum InputCapability" crates/reed/src/input.rs && grep -q "pub struct ReplayableInput" crates/reed/src/input.rs && grep -q "pub struct SeekableInput" crates/reed/src/input.rs && grep -q "pub struct InputMeta" crates/reed/src/input.rs && ! grep -rq "source_mut()" crates --include=*.rs'

check "BCF capability comes from InputMeta, not an Option probe (doc_v30 R07)" \
    bash -c '! grep -rq "bcf\.is_some()" crates --include=*.rs'

# doc_v30 P2 §2: one run context carries the budget and the abort signal into
# every driver. The per-driver `run_budget` constructor is gone, so a stage
# cannot derive a second budget (and `FilterOptions` cannot carry one).
check "the run context owns budget and abort (doc_v30 P2 §2)" \
    bash -c 'test -e crates/reed/src/context.rs &&
        grep -q "pub struct RunContext" crates/reed/src/context.rs &&
        test ! -e crates/reed/src/run_budget.rs &&
        ! grep -rq "run_budget::for_run" crates/reed/src --include=*.rs &&
        [ "$(grep -rn "fn for_run" crates/reed/src --include=*.rs | wc -l)" -eq 0 ]'

check "every run driver takes the run context (doc_v30 P2 §2)" \
    bash -c 'grep -q "ctx: &RunContext" crates/reed/src/filter/mod.rs &&
        grep -q "ctx: &RunContext" crates/reed/src/stats/driver.rs &&
        grep -q "ctx: &RunContext" crates/reed/src/stats/writers/reports.rs &&
        grep -q "ctx: &RunContext" crates/reed/src/convert/job/mod.rs &&
        grep -q "ctx: &RunContext" crates/reed/src/rename/mod.rs &&
        grep -q "ctx: &RunContext" crates/reed/src/pca/mod.rs &&
        grep -q "ctx: &RunContext" crates/reed/src/dapc/mod.rs &&
        grep -q "ctx: &RunContext" crates/reed/src/index.rs'

check "FilterOptions carries no second budget source (doc_v30 P2 §2)" \
    bash -c '! grep -q "pub mem_bytes" crates/reed/src/filter/mod.rs &&
        ! grep -q "pub threads" crates/reed/src/filter/mod.rs'

# doc_v30 R13: one `--fill-tags` registry. Each tag's name/type/number/
# description/compute lives in reed-core; the filter spec, the fill loop and
# the `##INFO` renderer derive from it. The CLI help is checked against the
# same registry by `fill_tags_help_names_every_registry_tag`.
check "the --fill-tags registry is defined once (doc_v30 R13)" \
    bash -c '[ "$(grep -rn "pub const RECOMPUTE_TAGS" crates --include=*.rs | wc -l)" -eq 1 ] &&
        grep -q "pub const RECOMPUTE_TAGS" crates/reed-core/src/count/recompute.rs &&
        grep -q "RECOMPUTE_TAGS.iter()" crates/reed-filter/src/spec/recompute.rs'

check "the INFO recompute definitions come from the registry (doc_v30 R13)" \
    bash -c 'grep -q "tag.name" crates/reed/src/filter/mod.rs &&
        ! grep -q "id: \"AC\"" crates/reed/src/filter/mod.rs &&
        ! grep -rEq "recompute\.(ac|an|ns|maf|hwe)\b" crates --include=*.rs'

# doc_v30 R14: the survivor's edit is a pure value with one constructor.
# The three emit paths (single pass, sketch emit, filtered passthrough) used
# to spell `EditPlan { .. }` themselves; a fourth must go through the value.
check "EditPlan is assembled only through SurvivorEdit (doc_v30 R14)" \
    bash -c 'grep -q "pub struct SurvivorEdit" crates/reed-filter/src/eval/edit_plan.rs &&
        [ "$(grep -rn "EditPlan {" crates/reed/src --include=*.rs | wc -l)" -eq 0 ]'

# doc_v30 P3: the facade is curated. Requests, reports, errors and the
# `run_*` entry points are public; the build-time plan, writers, distance/
# graph/reduce aliases and the fixture/open primitives live behind the one
# hidden `frontend` module (`#[doc(hidden)]` is allowed only there).
check "the prelude only re-exports facade items (doc_v30 P3)" \
    bash -c 'grep -q "pub use crate::{" crates/reed/src/prelude.rs &&
        ! grep -q "pub use reed_" crates/reed/src/prelude.rs'

check "front-end-only API lives behind the hidden frontend module (doc_v30 P3)" \
    bash -c '[ "$(grep -c "^#\[doc(hidden)\]" crates/reed/src/lib.rs)" -eq 1 ] &&
        grep -q "^pub mod prelude" crates/reed/src/lib.rs &&
        grep -q "^pub mod frontend" crates/reed/src/lib.rs &&
        ! grep -q "^pub use reed_filter::{expr" crates/reed/src/lib.rs &&
        ! grep -q "^pub use reed_distance as distance" crates/reed/src/lib.rs &&
        ! grep -q "^pub use reed_graph as graph" crates/reed/src/lib.rs &&
        ! grep -q "^pub use reed_reduce as reduce" crates/reed/src/lib.rs &&
        ! grep -q "^pub use crate::stats::writers::{" crates/reed/src/lib.rs'

# doc_v30 P4: the convert job is a directory module (types/sink, the job
# impl, the MSA helpers) and no source file grows past the audit ceiling.
check "convert job is split into job/{mod,convert_job,msa}.rs (doc_v30 P4)" \
    bash -c 'test -e crates/reed/src/convert/job/mod.rs &&
        test -e crates/reed/src/convert/job/convert_job.rs &&
        test -e crates/reed/src/convert/job/msa.rs &&
        test ! -e crates/reed/src/convert/job.rs'

check "no source file exceeds 2000 lines (doc_v30 P4)" \
    bash -c '[ "$(find crates -path "*/src/*" -name "*.rs" -exec wc -l {} + | awk "\$1 > 2000 && \$2 != \"total\"" | wc -l)" -eq 0 ]'

# doc_v30 R10: spills share AtomicFile's live table (so a sweep can never
# delete a spill that is still being written) and the spill subsystem sweeps
# its directory once at startup, removing only provably orphaned names.
check "spills register in the shared live table (doc_v30 R10)" \
    bash -c 'grep -q "pub fn register_live_temp" crates/reed-io/src/write/atomic.rs &&
        grep -q "register_live_temp(&path)" crates/reed/src/spill.rs &&
        grep -q "unregister_live_temp" crates/reed/src/spill.rs &&
        grep -q "fn sweep_stale_spills" crates/reed/src/spill.rs &&
        grep -q "fn a_live_spill_survives_the_sweep" crates/reed/src/spill.rs &&
        grep -q "fn a_dead_owners_spill_is_swept" crates/reed/src/spill.rs'

exit "$fail"
