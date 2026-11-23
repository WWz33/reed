#!/usr/bin/env bash
# ci/check-api.sh — the public API contract snapshots (docs/spec/api-contract.md).
#
# Every published crate is a lockstep stable surface (doc_v30 P1-7): `reed`
# plus its eight library dependencies. `cargo public-api` lists each crate's
# reachable items; the snapshot in ci/api/<crate>.txt is the reviewed surface:
# a change is a contract change, not a refactor, so it must be regenerated
# deliberately.
#
# Module re-exports are not expanded by the tool, so this is a lower bound on
# the reachable surface; the nameability gate closes that gap.
#
# Usage: bash ci/check-api.sh   (exit 1 on drift)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

crates=(
    reed-core reed-io reed-engine reed-stats reed-filter
    reed-graph reed-distance reed-reduce reed
)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
# The toolchain is spelled out rather than left as the floating name
# "nightly": `cargo public-api` silently switches a stable toolchain to
# `nightly`, and rustdoc JSON differs between nightlies, so an unpinned
# runner can fail this gate for reasons that are not in the code (doc_v30
# §4.2). This is the same pinned nightly the external-types job uses.
status=0
for crate in "${crates[@]}"; do
    cargo +nightly-2026-03-20 public-api -p "$crate" -s >"$tmp/$crate.txt"
    if ! diff -u "ci/api/$crate.txt" "$tmp/$crate.txt"; then
        echo "public API changed for $crate: review it and update ci/api/$crate.txt deliberately" >&2
        status=1
    fi
done
if [[ $status -ne 0 ]]; then
    echo "see docs/spec/api-contract.md before accepting a breaking change" >&2
    exit 1
fi
echo "public API snapshots match for ${#crates[@]} published crates"
