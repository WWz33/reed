#!/usr/bin/env bash
# ci/check-features.sh — ADR api-contract §4 / doc_v30 P3/P2-3.
#
# Features are for optional backends only. Every feature combination of every
# workspace crate must compile, including the empty set. The package set comes
# from the workspace rather than a hand-written list, so a new feature-bearing
# crate cannot escape the gate (doc_v30 P2-3).
set -euo pipefail
cd "$(dirname "$0")/.."

cargo hack --each-feature check --workspace
cargo check --workspace --no-default-features
