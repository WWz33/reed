#!/usr/bin/env bash
# gold/popgen/run_plink2_ld.sh — LD parity against plink2
# (docs/rounds/v5/05-algorithms.md §3.4).
#
# Usage:
#   bash gold/popgen/run_plink2_ld.sh [VCF] [POPS] [OUTDIR]
#
# Checks two things on the shared fixture (6 samples, 6 sites):
#   1. `--r2-unphased` per-pair r² against reed's `--ld-decay --ld-pairs`;
#   2. `--indep-pairwise W S R2` prune.in/prune.out against reed's.
# plink2 needs unique variant ids and ≥ 50 samples for LD, so the fixture
# run passes `--set-all-var-ids '@_#'` and `--bad-ld`. Two plink2
# conventions the comparison normalizes:
#   * `.vcor` prints six significant digits and drops pairs below
#     `--ld-window-r2` (0.2 by default), so the run sets that to 0 and the
#     gate is 1e-6 on the printed values;
#   * chromosome names are recoded to numeric codes (`Chr01` → 1), so the
#     prune keys are compared on the trailing position.
#
# Gate: per-pair r² within 2e-6 (plink2 prints %g: rounding up to 5e-7); identical prune
# sets (05 §3.4).
set -euo pipefail

GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}"
REED="${REED:-/tmp/reed-target/release/reed}"
PLINK2="${PLINK2:-$HOME/biosoft/plink2}"

VCF="${1:-$GOLD_DIR/fixtures/ld.vcf}"
POPS="${2:-$GOLD_DIR/fixtures/ld.pops}"
OUT="${3:-$SCRATCH/out/plink2_ld}"
W="${LD_W:-10}"
S="${LD_S:-1}"
R2="${LD_R2:-0.5}"

[[ -x "$REED" ]] || { echo "reed binary not executable: $REED" >&2; exit 1; }
[[ -x "$PLINK2" ]] || { echo "plink2 not executable: $PLINK2" >&2; exit 1; }
mkdir -p "$OUT"

echo "== plink2 =="
"$PLINK2" --vcf "$VCF" --set-all-var-ids '@_#' --r2-unphased \
    --ld-window-kb 100000 --ld-window-r2 0 --out "$OUT/plink2_r2" >/dev/null 2>&1
"$PLINK2" --vcf "$VCF" --set-all-var-ids '@_#' --indep-pairwise "$W" "$S" "$R2" \
    --bad-ld --out "$OUT/plink2_prune" >/dev/null 2>&1

echo "== reed =="
"$REED" stats "$VCF" --pop "$POPS" --variant-only --ld-decay 100k --ld-pairs \
    --indep-pairwise "$W $S $R2" -o "$OUT/reed" --verbosity 0

python3 - "$OUT" <<'PY'
import sys

out = sys.argv[1]

def pairs_plink(path):
    rows = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            a, b = int(c[1]), int(c[4])
            rows[(min(a, b), max(a, b))] = float(c[6])
    return rows

def pairs_reed(path):
    rows = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            a, b = int(c[1]), int(c[4])
            rows[(min(a, b), max(a, b))] = float(c[6])
    return rows

want = pairs_plink(f"{out}/plink2_r2.vcor")
got = pairs_reed(f"{out}/reed.ld")
bad = 0
worst = 0.0
for key, w in want.items():
    g = got.get(key)
    if g is None:
        print(f"FAIL pair {key}: missing in reed", file=sys.stderr); bad += 1; continue
    rel = abs(g - w) / abs(w) if w != 0 else abs(g - w)
    worst = max(worst, rel)
    if rel > 2e-6:
        print(f"FAIL pair {key}: reed={g} plink2={w} rel={rel:.3e}", file=sys.stderr); bad += 1
extra = set(got) - set(want)
if extra:
    print(f"FAIL: reed has {len(extra)} pair(s) plink2 did not report: {sorted(extra)[:5]}", file=sys.stderr)
    bad += len(extra)
print(f"{'PASS' if bad == 0 else 'FAIL'} r2: {len(want)} pairs, max rel {worst:.3e} (rtol 2e-6, plink2 %g)")

def prune(path):
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}

def numeric(keys):
    """plink2 recodes chromosome names; compare the trailing position."""
    return {k.replace(":", "_").rsplit("_", 1)[-1] for k in keys}

want_in = numeric(prune(f"{out}/plink2_prune.prune.in"))
want_out = numeric(prune(f"{out}/plink2_prune.prune.out"))
got_in = numeric(prune(f"{out}/reed.prune.in"))
got_out = numeric(prune(f"{out}/reed.prune.out"))
ok = want_in == got_in and want_out == got_out
if not ok:
    print(f"FAIL prune: in reed={sorted(got_in)} plink2={sorted(want_in)}", file=sys.stderr)
    print(f"             out reed={sorted(got_out)} plink2={sorted(want_out)}", file=sys.stderr)
    bad += 1
print(f"{'PASS' if ok else 'FAIL'} prune: {len(got_in)} kept / {len(got_out)} pruned")

sys.exit(1 if bad else 0)
PY
