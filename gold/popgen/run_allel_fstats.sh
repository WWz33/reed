#!/usr/bin/env bash
# gold/popgen/run_allel_fstats.sh — Patterson f-statistics parity against
# scikit-allel (docs/rounds/v5/05-algorithms.md §3.3).
#
# Usage:
#   bash gold/popgen/run_allel_fstats.sh [VCF] [POPS] [OUTDIR]
#
# Defaults to the shared fixture (3 populations, 2 samples each). On real
# data, pass a `--pop` file with at least 4 populations and give the tests
# explicitly through the FSTATS_* variables, e.g.:
#   FSTATS_F2="pop1,pop2" FSTATS_F3="pop3;pop1,pop2" \
#   FSTATS_F4="pop1,pop2;pop3,pop4" FSTATS_D="pop1,pop2;pop3,pop4" \
#   bash gold/popgen/run_allel_fstats.sh clean.vcf.gz clean.3pop.txt out
#
# Gate: estimates within rtol 1e-12 on the fixture (05 §3.3). The SE is
# compared only when the block grids coincide (reed blocks by bp, allel by
# variant count) — hence `--no-se` here.
set -euo pipefail

GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}"
REED="${REED:-/tmp/reed-target/release/reed}"
PY="${ALLEL_PYTHON:-/tmp/allelvenv/bin/python}"

VCF="${1:-$GOLD_DIR/fixtures/fstats.vcf}"
POPS="${2:-$GOLD_DIR/fixtures/fstats.pops}"
OUT="${3:-$SCRATCH/out/allel_fstats}"

[[ -x "$REED" ]] || { echo "reed binary not executable: $REED" >&2; exit 1; }
[[ -x "$PY" ]] || { echo "allel python not executable: $PY" >&2; exit 1; }
mkdir -p "$OUT"

F2="${FSTATS_F2:-popA,popB}"
F3="${FSTATS_F3:-popC;popA,popB}"
F4="${FSTATS_F4:-popA,popB;popC,popC}"
D="${FSTATS_D:-popA,popB;popA,popC}"
RATIO="${FSTATS_RATIO:-}"

args=(popg "$VCF" --pop "$POPS" -o "$OUT/reed" --verbosity 0)
[[ -n "$F2" ]] && args+=(--f2 "$F2")
[[ -n "$F3" ]] && args+=(--f3 "$F3")
[[ -n "$F4" ]] && args+=(--f4 "$F4")
[[ -n "$D" ]] && args+=(--dstat "$D")
[[ -n "$RATIO" ]] && args+=(--f4-ratio "$RATIO")
[[ -n "${FSTATS_BLOCK_BP:-}" ]] && args+=(--block-size "$FSTATS_BLOCK_BP")
# The fixture's sites sit within a few bp: 1 bp blocks make the jackknife
# meaningful (the production default is 5M).
[[ -z "${FSTATS_BLOCK_BP:-}" && "$VCF" == "$GOLD_DIR/fixtures/fstats.vcf" ]] && args+=(--block-size 1)

echo "== reed =="
"$REED" "${args[@]}"
cat "$OUT/reed.fstats"

echo "== scikit-allel =="
cmp_args=("$OUT/reed.fstats" "$VCF" "$POPS" --no-se)
[[ -n "${FSTATS_RTOL:-}" ]] && cmp_args+=(--rtol "$FSTATS_RTOL")
"$PY" "$GOLD_DIR/popgen/compare_fstats.py" "${cmp_args[@]}"
