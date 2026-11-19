#!/usr/bin/env bash
# ci/check-semver.sh — ADR api-contract §5 / doc_v30 P3.
#
# Compares the facade (`reed`, the only crate with a semver promise) against
# the latest release tag. A breaking change must come with a version bump:
# for 0.x, a minor bump is the breaking release, which is exactly what
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
cargo semver-checks -p reed --baseline-rev "$baseline"
