#!/usr/bin/env bash
# gold/popgen/run_pixy_bed.sh — `--bed-file` windows vs pixy `--bed_file`
# (docs/rounds/v6/06-roadmap.md Y6, docs/spec/cli.md §7.3).
#
# Usage:
#   bash gold/popgen/run_pixy_bed.sh [VCF] [POPS] [BED] [OUTDIR]
#
# Defaults to the shared 4-sample fixture split two ways (`sfs2.pops`) and a
# BED with an overlapping pair of rows, so the overlap rule (each covering
# region accumulates the site) is exercised. pixy needs a bgzipped, tabix-
# indexed VCF; the script indexes a copy when handed a plain one.
#
# Gate: window coordinates and avg_pi/no_sites/count_diffs/count_comparisons
# column-for-column equal (pixy prints 14 significant digits; the comparison
# rounds reed's f64 to the same width).
set -euo pipefail

GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}"
REED="${REED:-/tmp/reed-target/release/reed}"
PIXY="${PIXY:-python3 -m pixy}"

VCF="${1:-$GOLD_DIR/fixtures/sfs.vcf}"
POPS="${2:-$GOLD_DIR/fixtures/sfs2.pops}"
BED="${3:-$GOLD_DIR/fixtures/sfs.bed}"
OUT="${4:-$SCRATCH/out/pixy_bed}"

[[ -x "$REED" ]] || { echo "reed binary not executable: $REED" >&2; exit 1; }
[[ -f "$VCF" ]] || { echo "VCF not found: $VCF" >&2; exit 1; }
mkdir -p "$OUT"

# pixy wants bgzip + tabix. This repo's filesystem (dpc) lags mtimes
# several minutes behind writes, which makes htslib refuse any index as
# "older than the data file" — so the indexed copy and pixy's output both
# live on tmpfs; only the final comparison reads back from $OUT.
STAGE="$(mktemp -d /dev/shm/reed-pixy-bed.XXXXXX)"
trap 'rm -rf "$STAGE"' EXIT
GZ="$STAGE/in.vcf.gz"
if [[ "$VCF" == *.gz ]]; then
    cp "$VCF" "$GZ"
    cp "$VCF.tbi" "$GZ.tbi" 2>/dev/null || cp "$VCF.csi" "$GZ.csi" 2>/dev/null || true
else
    cp "$VCF" "$STAGE/in.vcf"
    bgzip -f "$STAGE/in.vcf"
fi
tabix -f -p vcf "$GZ"

echo "== pixy --bed_file =="
$PIXY --stats pi --vcf "$GZ" --populations "$POPS" --bed_file "$BED" \
    --bypass_invariant_check --output_folder "$STAGE/pixy" --n_cores 2 >/dev/null

echo "== reed --bed-file =="
"$REED" popg "$VCF" --pop "$POPS" --pi --bypass-invariant-check \
    --bed-file "$BED" -o "$OUT/reed" --verbosity 0

python3 - "$STAGE/pixy/pixy_pi.txt" "$OUT/reed_pi.txt" << 'PY'
import sys

def rows(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        out = []
        for line in fh:
            if line.strip():
                out.append(line.rstrip("\n").split("\t"))
    return header, out

h1, ref = rows(sys.argv[1])
h2, got = rows(sys.argv[2])
cols = ["pop", "chromosome", "window_pos_1", "window_pos_2", "avg_pi",
        "no_sites", "count_diffs", "count_comparisons", "count_missing"]
i1 = [h1.index(c) for c in cols]
i2 = [h2.index(c) for c in cols]
if len(ref) != len(got):
    print(f"FAIL: {len(got)} reed rows vs {len(ref)} pixy rows", file=sys.stderr)
    sys.exit(1)
bad = 0
for r, g in zip(ref, got):
    a = [r[i] for i in i1]
    b = [g[i] for i in i2]
    # avg_pi: pixy prints %.14g; compare numerically at that width.
    if float(a[4]) != float(f"{float(b[4]):.14g}") or a[:4] != b[:4] or a[5:] != b[5:]:
        bad += 1
        print(f"FAIL: reed={b} pixy={a}", file=sys.stderr)
print(f"{'PASS' if bad == 0 else 'FAIL'} pixy_bed: {len(ref)} rows, "
      f"{'all columns equal' if bad == 0 else str(bad) + ' mismatched'}")
sys.exit(0 if bad == 0 else 1)
PY
