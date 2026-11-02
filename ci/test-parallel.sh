#!/usr/bin/env bash
# Run the workspace test binaries concurrently.
#
# `cargo test --workspace` runs one test binary at a time (each binary's
# tests are parallel inside libtest, but the binaries themselves are
# serial), so the wall time is the SUM of every integration suite. This
# runner builds everything once, then runs each test executable with
# `xargs -P`, which brought the local suite from ~50 s to ~10 s
# (doc_v22 R22-7).
#
# Usage:
#   ci/test-parallel.sh              # JOBS=$(nproc), one test per binary
#   JOBS=8 TEST_THREADS=4 ci/...     # tune oversubscription
#
# Doctests are not emitted as executables by `cargo test --no-run`; they run
# afterwards through the normal driver.
set -uo pipefail
cd "$(dirname "$0")/.."

JOBS="${JOBS:-$(nproc)}"
TEST_THREADS="${TEST_THREADS:-1}"
# `cargo +stable` because this checkout pins a toolchain channel that may be
# incomplete on a dev host; override with e.g. CARGO=cargo in CI.
CARGO="${CARGO:-cargo +stable}"
logdir="$(mktemp -d)"
export logdir TEST_THREADS

# Build all test harnesses and list only the *test* targets (the same JSON
# stream also names bins and examples). `pipefail` + the explicit status
# check make a compile failure fatal: without it a broken integration-test
# target left the other ~55 binaries running and the gate reported green.
build_err="$logdir/build.stderr"
if ! $CARGO test --workspace --no-run --message-format=json 2>"$build_err" \
  | python3 -c '
import json, sys
for line in sys.stdin:
    try:
        msg = json.loads(line)
    except json.JSONDecodeError:
        continue
    profile = msg.get("profile") or {}
    exe = msg.get("executable")
    if profile.get("test") and exe:
        print(exe)
' | sort -u >"$logdir/executables"; then
  echo "test build failed:" >&2
  cat "$build_err" >&2
  exit 1
fi

if [ ! -s "$logdir/executables" ]; then
  echo "no test executables found (build failed?)" >&2
  exit 1
fi
expected=$(wc -l <"$logdir/executables")

start=$(date +%s)
# `--test-threads=1` per binary: the parallelism comes from running many
# binaries at once, so the total stays near $JOBS even on a large host.
if ! xargs -P "$JOBS" -I{} bash -c '
  exe="$1"
  name=$(basename "$exe")
  "$exe" --test-threads="$TEST_THREADS" >"$logdir/$name.log" 2>&1 \
    || echo "FAIL $exe" >>"$logdir/failures"
' _ {} <"$logdir/executables"; then
  echo "test runner failed (JOBS=$JOBS)" >&2
  exit 1
fi
end=$(date +%s)
ran=$(find "$logdir" -maxdepth 1 -name '*.log' | wc -l)
if [ "$ran" -ne "$expected" ]; then
  echo "ran $ran of $expected test binaries" >&2
  exit 1
fi

if [ -s "$logdir/failures" ]; then
  echo "failed binaries:" >&2
  sed 's/^/  /' "$logdir/failures" >&2
  grep -l 'FAILED' "$logdir"/*.log 2>/dev/null | while read -r log; do
    echo "--- $log" >&2
    grep -A20 '^failures:' "$log" >&2 || true
  done
  echo "logs: $logdir" >&2
  exit 1
fi

passed=$(grep -h '^test result' "$logdir"/*.log \
  | sed -E 's/.*ok\. ([0-9]+) passed; ([0-9]+) failed.*/\1/' \
  | awk '{p+=$1} END{print p+0}')
failed=$(grep -h '^test result' "$logdir"/*.log \
  | sed -E 's/.*ok\. ([0-9]+) passed; ([0-9]+) failed.*/\2/' \
  | awk '{f+=$1} END{print f+0}')
echo "test binaries: $passed passed, $failed failed, $((end - start))s wall, $JOBS jobs"

# Doctests are not harness executables.
if ! $CARGO test --workspace --doc 2>&1 | tail -5; then
  echo "doctests failed" >&2
  exit 1
fi
