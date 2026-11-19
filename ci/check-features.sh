#!/usr/bin/env bash
# ci/check-features.sh — ADR api-contract §4 / doc_v30 P3.
#
# Features are for optional backends only. Every feature combination of the
# one feature-bearing crate must compile, including the empty set.
set -euo pipefail
cd "$(dirname "$0")/.."

cargo hack --each-feature check -p reed-reduce
cargo check -p reed-reduce --no-default-features
