#!/usr/bin/env bash
# ci/check-doc-drift.sh — doc/code drift guards (doc_v4 03-architecture.md §3.6).
#
# Every guard asserts a claim a STATUS file makes about the code, so a
# reversal of the code fails here instead of silently contradicting the
# documentation (the 01 F4-2 class of error: "C3 not done" while the code
# was done).
#
# Usage: bash ci/check-doc-drift.sh   (exit 1 on drift)
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

# docs/rounds/v3/STATUS.md: C3 (reed-filter is noodles-free) is done.
check "reed-filter has no noodles dependency (doc_v3 C3)" \
    bash -c '! grep -q noodles crates/reed-filter/Cargo.toml'

# docs/rounds/v3/STATUS.md: C4 (reed-qc merged into reed-stats) is done.
check "the reed-qc crate is gone (doc_v3 C4)" \
    bash -c '[ ! -d crates/reed-qc ]'

# docs/rounds/v2/05-naming.md §4: no version suffixes on types.
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

# doc_v6 06 Y0 / A6-22: the unsound AVX-512 count kernel was deleted, leaving
# the mmap as the only audited `unsafe` site (doc_v6 03 §5). A new one must
# change the docs and this guard first.
check "exactly one allow(unsafe_code) site (doc_v6 06 Y0)" \
    bash -c '[ "$(grep -rl "allow(unsafe_code)" crates --include=*.rs | wc -l)" -eq 1 ]'

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

# doc_v5 06-layout / 07 X2: the directory tree moved; no top-level doc_vN/
# or scripts/ remains, and the docs index uses the new paths.
check "the top-level doc_vN/ dirs are gone (doc_v5 X2)" \
    bash -c '! ls -d doc_v[0-9] 2>/dev/null | grep -q .'
check "scripts/ is now ci/ (doc_v5 X2)" \
    bash -c '[ ! -d scripts ] && [ -f ci/check-doc-drift.sh ]'
check "the docs index has no doc_vN/ path (doc_v5 X2)" \
    bash -c '! grep -rq "doc_v[0-9]/" docs/README.md'
check "the gold products live under target/ (doc_v5 X2)" \
    bash -c '[ ! -d gold/work ] && [ ! -d gold/out ]'

# doc_v5 03 §5 / 07 X3: the summary consumes the reports scan; it has no
# scanner of its own any more.
check "the summary has no second scanner (doc_v5 X3)" \
    bash -c '! grep -qE "GenotypeScanner|QcScratch" crates/reed/src/stats/writers/summary.rs'

# doc_v5 05 §3.4 / X4.4: the LD path is chunk-boundary independent — the
# worker carry resets per chunk and the sink owns every cross-chunk pair.
check "the LD carry resets per chunk (doc_v5 X4.4)" \
    grep -q 'The LD carry is per chunk' crates/reed/src/stats/job.rs

# doc_v6 06 Y0 / E6-13: vendor vocabulary is out of the production code
# (test directories keep the reference-tool names, naming.md §2.5).
check "no pixy vocabulary in production src (doc_v6 E6-13)" \
    bash -c '! grep -rniE "\bpixy" crates/reed/src crates/reed-cli/src crates/reed-core/src crates/reed-io/src crates/reed-engine/src crates/reed-stats/src crates/reed-filter/src --include=*.rs'

# doc_v6 06 Y0 / E6-14: the non-standard `ponytail` markers are gone.
check "no ponytail markers remain (doc_v6 E6-14)" \
    bash -c '! grep -rq ponytail crates'

# doc_v6 06 Y0 / E6-16: one module style — no `foo.rs` + `foo/` sibling pairs.
check "no sibling foo.rs + foo/ module pairs (doc_v6 E6-16)" \
    bash -c 'for d in crates/*/src; do for f in "$d"/*.rs; do [ -d "${f%.rs}" ] && exit 1; done; done; exit 0'

# doc_v6 06 Y0 / E6-16: every product crate defines its error in error.rs.
check "the four product crates have error.rs (doc_v6 E6-16)" \
    bash -c '[ "$(ls crates/reed-core/src/error.rs crates/reed-io/src/error.rs crates/reed-engine/src/error.rs crates/reed-filter/src/error.rs 2>/dev/null | wc -l)" -eq 4 ]'

# doc_v6 06 Y0 / E6-11: the binary carries its build provenance.
check "reed-cli has a build.rs for --version provenance (doc_v6 E6-11)" \
    grep -q 'REED_BUILD_COMMIT' crates/reed-cli/build.rs

# doc_v6 06 Y0 / E6-10: the declared license has its texts.
check "LICENSE-MIT and LICENSE-APACHE exist (doc_v6 E6-10)" \
    bash -c '[ -f LICENSE-MIT ] && [ -f LICENSE-APACHE ]'

exit "$fail"
