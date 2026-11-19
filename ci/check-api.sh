#!/usr/bin/env bash
# ci/check-api.sh — the facade contract snapshot (docs/spec/api-contract.md §6.1).
#
# `cargo public-api` lists the reachable items of the `reed` facade. The
# snapshot in ci/api/reed.txt is the reviewed surface: a change is a contract
# change, not a refactor, so it must be regenerated deliberately:
#
# Module re-exports are not expanded by the tool, so this is a lower bound on
# the reachable surface (V30-T03); the nameability gate in P3 closes that gap.
#
# Usage: bash ci/check-api.sh   (exit 1 on drift)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
# The toolchain is spelled out rather than left as the floating name
# "nightly": `cargo public-api` silently switches a stable toolchain to
# `nightly`, and rustdoc JSON differs between nightlies, so an unpinned
# runner can fail this gate for reasons that are not in the code (doc_v30
# §4.2). This is the same pinned nightly the external-types job uses.
#
#   cargo +nightly-2026-03-20 public-api -p reed -s > ci/api/reed.txt
cargo +nightly-2026-03-20 public-api -p reed -s >"$tmp"
if ! diff -u ci/api/reed.txt "$tmp"; then
    echo "public API changed: review it and update ci/api/reed.txt deliberately" >&2
    echo "see docs/spec/api-contract.md (P0) before accepting a breaking change" >&2
    exit 1
fi
