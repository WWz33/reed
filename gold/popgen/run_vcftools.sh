#!/bin/bash
# vcftools --window-pi reference: the --variant-only π definition (04 §5).
# Usage: run_vcftools.sh VCF POPFILE OUTDIR [WINDOW]
set -euo pipefail
VCF="${1:?vcf}"; POPS="${2:?populations file}"; OUT="${3:?output dir}"; WINDOW="${4:-100000}"
mkdir -p "$OUT"
cut -f1 "$POPS" | sort -u > "$OUT/samples.txt"
while IFS=$'\t' read -r pop; do :; done < <(cut -f2 "$POPS" | sort -u)
for pop in $(cut -f2 "$POPS" | sort -u); do
  awk -v p="$pop" -F'\t' '$2==p {print $1}' "$POPS" > "$OUT/$pop.keep"
  vcftools --gzvcf "$VCF" --keep "$OUT/$pop.keep" --window-pi "$WINDOW" --out "$OUT/vcftools_pi_$pop"
done
