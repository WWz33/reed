#!/usr/bin/env bash
# gold/popgen/run_vcftools_fst.sh — global (r-population) WC84 Fst parity
# against vcftools (docs/rounds/v5/05-algorithms.md §3.2).
#
# Usage:
#   bash gold/popgen/run_vcftools_fst.sh [VCF] [POPS] [OUTDIR]
#
# vcftools takes one --weir-fst-pop per population and reports the *global*
# (r-population) WC84 estimator: the per-site WEIR_AND_COCKERHAM_FST column
# and the log's weighted Fst (= Σa/Σ(a+b+c)). reed's `--fst --fst-global`
# computes the same quantity.
#
# Gate: the genome-level value within 1e-9 (05 §3.2); the per-site values
# within 1e-6 because vcftools prints them with six significant digits, which
# is the resolution of the reference file itself.
set -euo pipefail

GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}"
REED="${REED:-/tmp/reed-target/release/reed}"
VCFTOOLS="${VCFTOOLS:-$HOME/biosoft/vcftools}"

VCF="${1:-$GOLD_DIR/fixtures/fstats.vcf}"
POPS="${2:-$GOLD_DIR/fixtures/fstats.pops}"
OUT="${3:-$SCRATCH/out/vcftools_fst}"

[[ -x "$REED" ]] || { echo "reed binary not executable: $REED" >&2; exit 1; }
[[ -x "$VCFTOOLS" ]] || { echo "vcftools not executable: $VCFTOOLS" >&2; exit 1; }
mkdir -p "$OUT"

# One --weir-fst-pop file per population (vcftools wants sample lists).
pop_files=()
while read -r pop; do
    f="$OUT/pop_$pop.txt"
    awk -F'\t' -v p="$pop" '$2 == p {print $1}' "$POPS" > "$f"
    pop_files+=(--weir-fst-pop "$f")
done < <(awk -F'\t' '!/^#/ && NF>0 {print $2}' "$POPS" | sort -u)

echo "== vcftools =="
"$VCFTOOLS" --vcf "$VCF" "${pop_files[@]}" --out "$OUT/vcftools" 2> "$OUT/vcftools.log"
grep -i "weighted Fst" "$OUT/vcftools.log" || true

echo "== reed =="
# Per-site windows (size 1) plus the genome level.
"$REED" stats "$VCF" --pop "$POPS" --variant-only --fst --fst-global \
    --window 1 --levels window -o "$OUT/reed_site" --verbosity 0
"$REED" stats "$VCF" --pop "$POPS" --variant-only --fst --fst-global \
    --window 100k --levels genome -o "$OUT/reed_genome" --verbosity 0

python3 - "$OUT" <<'PY'
import re
import sys

out = sys.argv[1]

def read_reed(path):
    rows = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            rows[int(c[2])] = float(c[4])
    return rows

def read_vcftools(path):
    rows = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            rows[int(c[1])] = float(c[2])
    return rows

rtol_site = 2e-6  # vcftools prints %g: six significant digits (rounding up to 5e-7)
rtol = 1e-9
reed = read_reed(f"{out}/reed_site_fst_global.txt")
vcft = read_vcftools(f"{out}/vcftools.weir.fst")
bad = 0
worst = 0.0
for pos, want in vcft.items():
    got = reed.get(pos)
    if got is None:
        print(f"FAIL per-site pos {pos}: missing in reed", file=sys.stderr)
        bad += 1
        continue
    rel = abs(got - want) / abs(want) if want != 0 else abs(got - want)
    worst = max(worst, rel)
    if rel > rtol_site:
        print(f"FAIL per-site pos {pos}: reed={got} vcftools={want} rel={rel:.3e}", file=sys.stderr)
        bad += 1
print(f"{'PASS' if bad == 0 else 'FAIL'} per-site: {len(vcft)} sites, max rel {worst:.3e} (rtol {rtol_site:g})")

# Genome-level: reed's genome row vs vcftools' logged weighted Fst.
with open(f"{out}/reed_genome_fst_global.txt") as fh:
    fh.readline()
    genome = float(fh.readline().split("\t")[4])
log = open(f"{out}/vcftools.log").read()
m = re.search(r"weighted Fst estimate:\s*([-\d.eE+]+)", log)
if not m:
    print("FAIL: no weighted Fst in the vcftools log", file=sys.stderr)
    sys.exit(1)
want = float(m.group(1))
rel = abs(genome - want) / abs(want) if want != 0 else abs(genome - want)
ok = rel <= rtol
print(f"{'PASS' if ok else 'FAIL'} genome: reed={genome} vcftools={want} rel={rel:.3e}")
sys.exit(0 if (ok and bad == 0) else 1)
PY
