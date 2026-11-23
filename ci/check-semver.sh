#!/usr/bin/env bash
# ci/check-semver.sh — ADR api-contract §5 / doc_v30 P1-7.
#
# Every published crate is a lockstep stable surface, so every one is checked
# against the latest release tag. A breaking change must come with a version
# bump: for 0.x, a minor bump is the breaking release, which is exactly what
# cargo-semver-checks enforces (it skips the per-item lints only when the
# detected change is already major).
#
# CI checks out the full history (`fetch-depth: 0`) so the tags exist. Run
# from the repo root; exits 1 on a breaking change without a bump.
#
# Usage: bash ci/check-semver.sh [baseline-rev]
set -euo pipefail
cd "$(dirname "$0")/.."

if [[ $# -gt 0 ]]; then
    baseline="$1"
else
    baseline=$(git tag --list 'v[0-9]*' --sort=-version:refname | sed -n '1p')
fi
if [[ -z "${baseline:-}" ]]; then
    echo "no versioned release tag found: cannot run semver-checks" >&2
    exit 1
fi
echo "semver baseline: $baseline"
# `cargo doc` (the documentation gate) leaves HTML in the target directory.
# cargo-semver-checks then sees a fresh doc artifact and does not regenerate
# rustdoc JSON, failing with "could not find expected rustdoc output". Clear
# the doc artifacts so each crate is genuinely rebuilt for the check.
target_dir=$(cargo metadata --format-version 1 --no-deps |
    python3 -c 'import json, sys; print(json.load(sys.stdin)["target_directory"])')
rm -rf "$target_dir/doc"
crates=(
    reed-core reed-io reed-engine reed-stats reed-filter
    reed-graph reed-distance reed-reduce reed
)
for crate in "${crates[@]}"; do
    echo "== $crate =="
    cargo semver-checks -p "$crate" --baseline-rev "$baseline"
done
