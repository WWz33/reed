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

# doc_v4 02 §5 (W0b): the plink2 spellings are the real flags; the withdrawn
# --max-missing family survives only as parse-time errors (never as aliases).
# NOTE: these three guards report DRIFT until W0b lands -- expected, since the
# documents are reverted and the code is not yet.
check "the withdrawn --max-missing flags error (doc_v4 02 §5.2)" \
    grep -q 'fn removed_max_missing' crates/reed-cli/src/cli/filter_flags.rs
check "no text points at the withdrawn --max-missing name (doc_v4 02 §5.2)" \
    bash -c '! grep -rq "use --max-missing" crates/*/src'
check "the funnel ids are after_mask.missing / sample.mind (doc_v4 02 §5.3)" \
    grep -q 'after_mask.missing' crates/reed-filter/src/funnel.rs

# doc_v4 W1: index documents and enforces its file-only input.
check "reed index rejects stdin (doc_v4 W1/F4-4)" \
    grep -q 'stdin is not supported' crates/reed-cli/src/translate/index.rs

# doc_v4 03 §3.3: the verbosity contract lives in one place.
check "the shared Verbosity enum exists (doc_v4 03 §3.3)" \
    grep -q 'pub(crate) enum Verbosity' crates/reed-cli/src/cli/verbosity.rs

# doc_v5 03 §2 / A5-2: exactly two audited `unsafe` sites (the AVX-512 count
# kernel and the mmap). A new one must change the docs and this guard first.
check "exactly two allow(unsafe_code) sites (doc_v5 03 §2)" \
    bash -c '[ "$(grep -rl "allow(unsafe_code)" crates --include=*.rs | wc -l)" -eq 2 ]'

# doc_v5 03 §6 / A5-3: the stale names are gone from the code.
check "no code says sample.max_missing (doc_v5 A5-3)" \
    bash -c '! grep -rq "sample\.max_missing" crates/*/src'
check "no code says --fst-pops (doc_v5 A5-3)" \
    bash -c '! grep -rq -e "--fst-pops" crates/*/src'

# doc_v5 03 §6 / C-2: the deprecated `reed qc` subcommand is gone.
check "reed qc is gone (doc_v5 03 C-2)" \
    bash -c '! grep -q "Command::Qc" crates/reed-cli/src/main.rs'

# doc_v5 08 §6.4: the front/back-end boundary. The front end depends on the
# facade alone; the backend never prints; the facade is a re-export list.
check "reed-cli depends only on the facade (doc_v5 08 §6.4.1)" \
    bash -c '! grep -qE "reed-(core|io|engine|stats|filter)" crates/reed-cli/Cargo.toml'
check "the backend never prints to stdout/stderr (doc_v5 08 §6.4.2)" \
    bash -c '! grep -rqE "eprintln!|println!" crates/reed/src crates/reed-core/src crates/reed-io/src crates/reed-engine/src crates/reed-stats/src crates/reed-filter/src --include=*.rs'
check "the facade has no pub mod (doc_v5 08 §6.4.3)" \
    bash -c '[ "$(grep -c "^pub mod" crates/reed/src/lib.rs)" -eq 0 ]'

# doc_v5 08 §4 L-1: the README describes v0.3, not the v0.2 API.
check "the README has no v0.2 API examples (doc_v5 08 §4 L-1)" \
    bash -c '! grep -qE "WindowBuffer|exec::execute|--groups" README.md'

# doc_v5 08 §4 L-3: reed-io no longer re-exports noodles.
check "no noodles re-export remains (doc_v5 08 §4 L-3)" \
    bash -c '! grep -rq "noodles_vcf" crates'

exit "$fail"
