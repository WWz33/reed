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
    bash -c '! grep -rq "fn gt_codes_into\|fn gt_shape" crates --include=*.rs | grep -v "gt_decode.rs"'

# V25-10: the dead linalg copy stays deleted and the env back door is gone.
check "the dead reed-reduce linalg file is gone (doc_v25 V25-10, R12)" \
    bash -c 'test ! -e crates/reed-reduce/src/linalg.rs && ! grep -rq "mod linalg\\|crate::linalg" crates/reed-reduce/src --include=*.rs'

check "no STREAM_DEBUG env back door remains (doc_v25 V25-10)" \
    bash -c '! grep -rq "STREAM_DEBUG" crates --include=*.rs'

# V25-09: every on-disk spill goes through the exclusive TempSpill; no
# private name builder may grow back.
check "spills use the shared TempSpill (doc_v25 V25-09)" \
    bash -c 'grep -q "TempSpill::create" crates/reed/src/input.rs && grep -q "TempSpill::create" crates/reed/src/convert/job.rs && ! grep -rq "fn spill_path()" crates/reed/src --include=*.rs'

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

exit "$fail"
