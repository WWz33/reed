#!/bin/bash
# pixy reference for the popgen parity gate (04-algorithms.md §5).
# Usage: run_pixy.sh VCF POPFILE OUTDIR [WINDOW]
set -euo pipefail
VCF="${1:?vcf}"; POPS="${2:?populations file}"; OUT="${3:?output dir}"; WINDOW="${4:-100000}"
mkdir -p "$OUT"
python3 -m pixy --stats pi dxy watterson_theta tajima_d --vcf "$VCF" --populations "$POPS" \
  --window_size "$WINDOW" --bypass_invariant_check --output_folder "$OUT/pixy_bypass" --n_cores 8
python3 -m pixy --stats pi dxy --vcf "$VCF" --populations "$POPS" \
  --window_size "$WINDOW" --output_folder "$OUT/pixy_default" --n_cores 8
for method in hudson wc; do
  python3 -m pixy --stats fst --fst_type "$method" --vcf "$VCF" --populations "$POPS" \
    --window_size "$WINDOW" --bypass_invariant_check \
    --output_folder "$OUT/pixy_fst_$method" --n_cores 8
done
