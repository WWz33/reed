#!/usr/bin/env bash
# Reproducible comparison protocol (doc_v24 V24-02).
#
# Every (binary, input) pair is warmed up once and then measured REPS times.
# The tool order alternates per repetition (ABBA...) so neither binary always
# pays the cold-cache/first-run cost, stderr is kept per run, and a non-zero
# exit or an n_rec mismatch fails the harness instead of silently recording
# secs=0. Results are JSONL; the median per (input, binary) is printed as
# JSON at the end.
#
# Usage:
#   bench/run_matrix.sh [INPUT...]
#   REPS=7 WARMUP=1 OUT=/tmp/b.jsonl bench/run_matrix.sh a.bcf b.vcf.gz
#
# Env:
#   R        release target dir (default /tmp/bench-target/release)
#   REPS     measured repetitions (default 5, minimum 3)
#   WARMUP   discarded warmup runs per (binary, input) (default 1)
#   OUT      JSONL path (default /tmp/reed-bench-<timestamp>.jsonl)
set -euo pipefail

R="${R:-/tmp/bench-target/release}"
REPS="${REPS:-5}"
WARMUP="${WARMUP:-1}"
OUT="${OUT:-/tmp/reed-bench-$(date +%Y%m%d-%H%M%S).jsonl}"
DEFAULT_INPUTS=(
  "/sri/home/sri2025201067/reed/clean.bcf"
  "/sri/home/sri2025201067/reed/clean.vcf.gz"
  "/sri/home/sri2025201067/reed/clean_bg.vcf.gz"
)
if [ "$#" -gt 0 ]; then
  INPUTS=("$@")
else
  INPUTS=("${DEFAULT_INPUTS[@]}")
fi
if [ "$REPS" -lt 3 ]; then
  echo "REPS must be >= 3 (got $REPS): a median of two is not a protocol" >&2
  exit 2
fi

for bin in noo hts; do
  if [ ! -x "$R/bench-$bin" ]; then
    echo "missing $R/bench-$bin; build with:" >&2
    echo "  cargo +stable build --release --manifest-path bench/Cargo.toml -p bench-$bin" >&2
    exit 2
  fi
done

LOGDIR="$(mktemp -d)"
rustc_version="$(rustc +stable -V 2>/dev/null || echo unknown)"
printf '{"type":"meta","date":"%s","host":"%s","nproc":%s,"rustc":"%s","reps":%d,"warmup":%d}\n' \
  "$(date -Iseconds)" "$(hostname)" "$(nproc)" "$rustc_version" "$REPS" "$WARMUP" >"$OUT"

# One timed run; echoes "secs n_rec" and keeps stdout+stderr in $LOGDIR.
run_one() {
  local bin="$1" input="$2" tag="$3"
  local log="$LOGDIR/$tag.log"
  local out rc secs n_rec
  set +e
  out="$("$R/bench-$bin" "$input" 2>"$log")"
  rc=$?
  set -e
  if [ "$rc" -ne 0 ]; then
    echo "bench-$bin $input failed (exit $rc); stderr kept at $log" >&2
    exit 1
  fi
  secs="$(printf '%s\n' "$out" | awk -F= '/^secs=/{print $2}')"
  n_rec="$(printf '%s\n' "$out" | awk -F= '/^n_rec=/{print $2}')"
  if [ -z "$secs" ] || [ -z "$n_rec" ] || [ "$n_rec" -eq 0 ]; then
    echo "bench-$bin $input produced no usable output (secs='$secs' n_rec='$n_rec'); log $log" >&2
    exit 1
  fi
  printf '%s %s\n' "$secs" "$n_rec"
}

for input in "${INPUTS[@]}"; do
  if [ ! -f "$input" ]; then
    echo "input not found: $input" >&2
    exit 2
  fi
  for warm in $(seq 1 "$WARMUP"); do
    for bin in noo hts; do
      read -r secs n_rec < <(run_one "$bin" "$input" "warm-$bin-$warm")
      echo "warmup $input bench-$bin ${secs}s n_rec=$n_rec (discarded)" >&2
    done
  done

  for rep in $(seq 1 "$REPS"); do
    # Alternate the order: even reps run noo first, odd reps hts first.
    if [ $((rep % 2)) -eq 0 ]; then order=(noo hts); else order=(hts noo); fi
    for bin in "${order[@]}"; do
      read -r secs n_rec < <(run_one "$bin" "$input" "rep-$bin-$rep")
      printf '{"type":"run","input":"%s","binary":"bench-%s","rep":%d,"secs":%s,"n_rec":%s}\n' \
        "$input" "$bin" "$rep" "$secs" "$n_rec" >>"$OUT"
    done
  done
done

python3 - "$OUT" <<'PY'
import json, statistics, sys
from collections import defaultdict

runs = defaultdict(list)
records = defaultdict(set)
with open(sys.argv[1], encoding="utf-8") as fh:
    for line in fh:
        row = json.loads(line)
        if row.get("type") != "run":
            continue
        key = (row["input"], row["binary"])
        runs[key].append(float(row["secs"]))
        records[row["input"]].add(int(row["n_rec"]))

summary = []
for (input_path, binary), samples in sorted(runs.items()):
    summary.append({
        "input": input_path,
        "binary": binary,
        "reps": len(samples),
        "median_secs": round(statistics.median(samples), 6),
        "min_secs": round(min(samples), 6),
        "max_secs": round(max(samples), 6),
    })
out = {
    "jsonl": sys.argv[1],
    "combined_n_rec_consistent": all(len(v) == 1 for v in records.values()),
    "n_rec_by_input": {k: sorted(v) for k, v in records.items()},
    "medians": summary,
}
print(json.dumps(out, indent=2))
if not out["combined_n_rec_consistent"]:
    print("ERROR: the two binaries did not process the same record count", file=sys.stderr)
    sys.exit(1)
PY

echo "jsonl: $OUT" >&2
echo "stderr logs: $LOGDIR" >&2
