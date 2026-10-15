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
check "the front end never compiles plans or assigns groups (doc_v5 08 §6.2)" \
    bash -c '! grep -rqE "assign_groups|plan::compile" crates/reed-cli/src --include=*.rs'
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
    bash -c '[ ! -f docs/README.md ] || ! grep -rq "doc_v[0-9]/" docs/README.md'
check "the gold products live under target/ (doc_v5 X2)" \
    bash -c '[ ! -d gold/work ] && [ ! -d gold/out ]'

# doc_v5 03 §5 / 07 X3: the summary consumes the reports scan; it has no
# scanner of its own any more.
check "the summary has no second scanner (doc_v5 X3)" \
    bash -c '! grep -qE "GenotypeScanner|QcScratch" crates/reed/src/stats/writers/summary.rs'

# doc_v5 05 §3.4 / X4.4: the LD path is chunk-boundary independent — the
# worker carry resets per chunk and the sink owns every cross-chunk pair.
# doc_v6 02 §2.2 / 06 Y1: the stats LD path collects the site-major matrix
# (the reducer only concatenates runs); the head/tail carry is retired.
check "the stats LD path is matrix-based, carry retired (doc_v6 06 Y1)" \
    bash -c '! grep -rq "The LD carry is per chunk" crates/reed/src/stats/ && grep -rq "ld_matrix" crates/reed/src/stats/'

# doc_v6 Y6 S4: job.rs is split into request/run/sink, and the flat table
# renders from structured rows (no text buffer in the sink, no tsv field).
check "stats/job.rs is split into request/run/sink/driver (doc_v6 Y6 S4)" \
    bash -c '[ ! -f crates/reed/src/stats/job.rs ] && [ -f crates/reed/src/stats/run.rs ] && [ -f crates/reed/src/stats/sink.rs ] && [ -f crates/reed/src/stats/request.rs ] && [ -f crates/reed/src/stats/driver.rs ]'
check "the stats report carries no pre-rendered tsv (doc_v6 Y6 S4)" \
    bash -c '! grep -q "pub tsv" crates/reed/src/stats/request.rs'

# doc_v6 06 Y0 / E6-13: vendor vocabulary is out of the production code
# (test directories keep the reference-tool names, naming.md §2.5).
check "no pixy vocabulary in production src (doc_v6 E6-13)" \
    bash -c '! grep -rniE "\bpixy" crates/*/src --include=*.rs'

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

# doc_v8 01 §1 / D8-5b + D8-6: every rayon par_iter/into_par_iter in
# production src sits inside a run_auxiliary closure's *parallel branch*
# (`if parallel {`), never on the rayon global pool. D8-6 taught the lesson:
# a par_iter in the serial path still lazily spawns the 56-thread global
# pool. New par_* outside an `if parallel {` arm fails here and must be
# registered in docs/rounds/v8/01-review.md §1 first.
check "rayon parallelism is behind the parallel branch (doc_v8 D8-1/D8-6)" \
    bash -c 'for f in $(grep -rl "par_iter\|into_par_iter" crates/reed/src --include=*.rs); do
        n_par=$(grep -c "par_iter\|into_par_iter" "$f")
        n_guard=$(grep -cE "if parallel \{" "$f")
        [ "$n_par" -le "$n_guard" ] || { echo "unguarded par_iter in $f"; exit 1; }
    done; exit 0'

# doc_v8 D8-4: the per-contig caps read budget.auxiliary(), not the raw CLI
# thread count.
check "per-contig caps use budget.auxiliary (doc_v8 D8-4)" \
    bash -c '! grep -rn "max_prune: opts.threads\|max_ld_jobs: opts.threads" crates/reed/src'

# doc_v9 D9-7: spec §7.7 describes backpressure, not "over budget → error".
# The implementation follows v2/03 §4.2; the spec line was stale.
check "spec §7.7 describes backpressure, not an error (doc_v9 D9-7)" \
    bash -c '! sed -n "/### 7.7/,/^### 7.8/p" docs/spec/architecture.md | grep -q "超出时报 \`EngineError::MemBudgetExceeded\`"'

# doc_v9 D9-2: the pipeline channels must honour the send result
# (`res.is_ok()`), never bind it to `_` and report success on a disconnect.
check "pipeline sends check their result (doc_v9 D9-2)" \
    bash -c '! grep -q "send(.*) -> _ => true" crates/reed-engine/src/pipeline.rs'

# doc_v9 D9-1: the reducer records its own errors before returning (`bail`):
# the helper exists and both error exits (the loop's `Err` and the panic arm)
# call it.
check "reducer errors go through bail (doc_v9 D9-1)" \
    bash -c 'grep -q "fn bail" crates/reed-engine/src/pipeline.rs && [ "$(grep -c "bail(" crates/reed-engine/src/pipeline.rs)" -ge 3 ]'

# doc_v9 D9-3: the newline path finalizes the stitcher on both the parallel
# and sequential paths, so an unterminated final record cannot be dropped.
check "NewlineStitcher is finalized at stream end (doc_v9 D9-3)" \
    bash -c '[ "$(grep -c "stitcher.finish()" crates/reed-engine/src/pipeline.rs)" -ge 2 ]'

# doc_v10 D10-1: the GT width follows the encoded value domain — INT8 ends
# at allele 62, INT16 at 16382 (2a+3 <= 127 / 32767). A regression to the
# old allele-domain thresholds (32766) truncates codes with `as i16`.
check "GT width follows the encoded value domain (doc_v10 D10-1)" \
    bash -c 'grep -q "fn gt_width" crates/reed-io/src/bcf/native_encode.rs &&
        grep -q "max_allele <= 16382" crates/reed-io/src/bcf/native_encode.rs &&
        ! grep -q "max_allele <= 32766" crates/reed-io/src/bcf/native_encode.rs'

# doc_v10 D10-2: a (NULL, count>0) descriptor is rejected at the one read
# gate (`Cursor::typed`), so the empty payload never reaches the renderers.
check "null descriptors with a nonzero count are rejected (doc_v10 D10-2)" \
    grep -q "null-typed descriptor has a nonzero count" crates/reed-io/src/bcf/record.rs

# doc_v10 D10-3: the dictionary is sparse (a BTreeMap), so an attacker
# `IDX=` cannot drive an allocation; the htslib ceiling is enforced.
check "the BCF dictionary is sparse and IDX-bounded (doc_v10 D10-3)" \
    bash -c 'grep -q "strings: BTreeMap<u32, Vec<u8>>" crates/reed-io/src/bcf/dict.rs &&
        grep -q "IDX_MAX" crates/reed-io/src/bcf/dict.rs'

# doc_v10 D10-4: BSIZE below the 18-byte header length is rejected in
# `parse_block_size`, before any caller slices or underflows.
check "undersized BGZF blocks are rejected (doc_v10 D10-4)" \
    bash -c 'grep -q "bsize + 1 < BGZF_HEADER_LEN" crates/reed-io/src/bgzf/block.rs'

# doc_v10 D10-5: the text index carries possibly-incomplete lines to the
# reducer and decodes each record once, whole (`feed`/`CarriedLine`); the
# worker no longer parses head/tail fragments.
check "the text index decodes carried fragments on the reducer (doc_v10 D10-5)" \
    bash -c 'grep -q "fn feed" crates/reed/src/index.rs &&
        grep -q "CarriedLine" crates/reed/src/index.rs &&
        grep -q "HeadEntry::Carried" crates/reed/src/index.rs'

# doc_v10 D10-6: plain gzip is read as a multi-member stream (htslib resets
# into the next member), not a single-member decoder.
check "plain gzip is decoded multi-member (doc_v10 D10-6)" \
    bash -c 'grep -q "MultiGzDecoder" crates/reed-io/src/chunk/plain.rs &&
        ! grep -qE "[^i]GzDecoder::new" crates/reed-io/src/chunk/plain.rs'

# doc_v10 D10-7: the output header declares PASS even when the input does
# not (htslib injects it on both sides).
check "the output header declares PASS (doc_v10 D10-7)" \
    bash -c 'grep -q "fn ensure_pass_filter" crates/reed/src/filter/mod.rs &&
        grep -q "ensure_pass_filter(out_header)" crates/reed/src/filter/mod.rs'

# doc_v10 D10-8: record field sizes go through the checked `field_bytes`
# helper, never a raw `type_width(ty) * count` product.
check "BCF field sizes are computed with checked multiplication (doc_v10 D10-8)" \
    bash -c 'grep -q "pub fn field_bytes" crates/reed-io/src/bcf/record.rs &&
        ! grep -rn "type_width(ty) \* count" crates/reed-io/src/bcf/render.rs crates/reed-io/src/bcf/native_encode.rs'

# ---- doc_v11 (facade API) --------------------------------------------------
# N1: `Levels`'s default is the window level, not the all-false bool tuple.
check "the default Levels is the window level (doc_v11 A11-1)" \
    bash -c 'grep -q "impl Default for Levels" crates/reed/src/stats/request.rs &&
        ! grep -q "derive(Clone, Copy, Debug, Default, PartialEq, Eq)" crates/reed/src/stats/request.rs'

# N2: the stats driver borrows the request like its sibling drivers.
check "run_stats borrows its request (doc_v11 A11-2)" \
    bash -c 'grep -q "request: &StatsRequest" crates/reed/src/stats/driver.rs'

# N3: `to_args` reports unspellable requests as usage errors, never panics.
check "StatsRequest::to_args errors instead of panicking (doc_v11 A11-3)" \
    bash -c 'grep -q "pub fn to_args(&self, group_names: &\[String\]) -> Result<Vec<String>, crate::error::Error>" crates/reed/src/stats/render.rs &&
        ! grep -q "panic!" crates/reed/src/stats/render.rs &&
        ! grep -q "assert!" crates/reed/src/stats/render.rs'

# N5: `Error` is a sealed struct; the variant payload types stay out of the
# public surface.
check "the facade Error is sealed (doc_v11 A11-5)" \
    bash -c 'grep -q "pub struct Error(Box<ErrorInner>)" crates/reed/src/error.rs &&
        grep -q "^enum ErrorInner" crates/reed/src/error.rs'

# N6: the four drivers and their reports are all at the root.
check "run_reports and its types are at the facade root (doc_v11 A11-6)" \
    bash -c 'grep -q "run_reports" crates/reed/src/lib.rs &&
        grep -q "ReportsOutcome" crates/reed/src/lib.rs &&
        grep -q "FilterReport" crates/reed/src/lib.rs'

# N8: the dead re-exports stay gone (bgzf/chunk are live front-end helpers).
check "no dead LdWindow/header re-exports (doc_v11 A11-8)" \
    bash -c '! grep -q "LdWindow" crates/reed/src/lib.rs &&
        ! grep -q "reed_io::{bgzf, chunk, header}" crates/reed/src/lib.rs'

# N9: OpenedInput keeps its transport fields crate-internal.
check "OpenedInput fields are crate-internal (doc_v11 A11-9)" \
    bash -c '! grep -q "pub source: Box<dyn ChunkSource>" crates/reed/src/input.rs &&
        ! grep -q "pub facts: HeaderFacts" crates/reed/src/input.rs'

# N13: one `--flag value` renderer shared by both to_args contracts.
check "the to_args flag renderer is defined once (doc_v11 A11-13)" \
    bash -c 'grep -q "pub fn push_flag" crates/reed-core/src/args.rs &&
        [ "$(grep -rn "fn flag(" crates --include=*.rs | wc -l)" -eq 0 ]'

# The API-surface gate: a test that imports only the facade.
check "the facade API-surface test exists (doc_v11)" \
    bash -c 'grep -q "use reed::{" crates/reed/tests/api_surface.rs'

exit "$fail"
