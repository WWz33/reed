#!/usr/bin/env bash
# scripts/check-doc-drift.sh — doc/code drift guards (doc_v4 03-architecture.md §3.6).
#
# Every guard asserts a claim a STATUS file makes about the code, so a
# reversal of the code fails here instead of silently contradicting the
# documentation (the 01 F4-2 class of error: "C3 not done" while the code
# was done).
#
# Usage: bash scripts/check-doc-drift.sh   (exit 1 on drift)
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

fail=0
check() {
    local label="$1"; shift
    if "$@"; then
        printf 'ok   %s\n' "$label"
    else
        printf 'DRIFT %s\n' "$label" >&2
        fail=1
    fi
}

# doc_v3/STATUS.md: C3 (reed-filter is noodles-free) is done.
check "reed-filter has no noodles dependency (doc_v3 C3)" \
    bash -c '! grep -q noodles crates/reed-filter/Cargo.toml'

# doc_v3/STATUS.md: C4 (reed-qc merged into reed-stats) is done.
check "the reed-qc crate is gone (doc_v3 C4)" \
    bash -c '[ ! -d crates/reed-qc ]'

# doc_v2/05-naming.md §4: no version suffixes on types.
check "no FilterSpecV2/SiteRulesV2 remain (doc_v2 05 §4)" \
    bash -c '! grep -rq "FilterSpecV2\|SiteRulesV2" crates/*/src'

# doc_v4 W0: the keep-direction spellings exist and the drop-direction flags
# are parse-time errors (never silent aliases).
check "the --max-missing family exists (doc_v4 W0)" \
    grep -q -- '--max-missing-sample' crates/reed-cli/src/cli.rs
check "removed --geno errors with the conversion (doc_v4 W0)" \
    grep -q 'was removed: it stated the fraction to drop' crates/reed-cli/src/cli/filter_flags.rs

# doc_v4 W1: index documents and enforces its file-only input.
check "reed index rejects stdin (doc_v4 W1/F4-4)" \
    grep -q 'stdin is not supported' crates/reed-cli/src/translate/index.rs

# doc_v4 03 §3.3: the verbosity contract lives in one place.
check "the shared Verbosity enum exists (doc_v4 03 §3.3)" \
    grep -q 'pub(crate) enum Verbosity' crates/reed-cli/src/cli/verbosity.rs

exit "$fail"
