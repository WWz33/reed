#!/bin/bash
# pixy reference for the popgen parity gate (04-algorithms.md §5).
#
# The four pixy invocations write to distinct folders, so they run
# concurrently (doc_v22 R22-7). JOBS bounds the process count and each pixy
# gets `nproc/JOBS` cores (override with PIXY_CORES). Defaults to the local
# subset callset, e.g.:
#   ./run_pixy.sh ../../../subvcf.vcf.gz ../../../clean.3pop.txt /tmp/popgen_sub
#
# `pixy_default` (without --bypass_invariant_check) needs an all-sites input
# with ALT="." records. The documented default input is the variant-only
# subset, where that case can only fail ("no invariant sites"), so it is
# added only when the VCF actually carries invariant sites. The other three
# runs (bypass/fst) are the meaningful comparisons on variant-only input.
#
# Usage: run_pixy.sh VCF POPFILE OUTDIR [WINDOW]
set -euo pipefail
VCF="${1:?vcf}"; POPS="${2:?populations file}"; OUT="${3:?output dir}"; WINDOW="${4:-100000}"
GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JOBS="${JOBS:-4}"
cores=$(nproc)
PER_JOB="${PIXY_CORES:-$(( cores / JOBS ))}"
[[ "$PER_JOB" -ge 1 ]] || PER_JOB=1

mkdir -p "$OUT"
cmds=(
  "python3 -m pixy --stats pi dxy watterson_theta tajima_d --vcf \"$VCF\" --populations \"$POPS\" --window_size \"$WINDOW\" --bypass_invariant_check --output_folder \"$OUT/pixy_bypass\" --n_cores $PER_JOB"
)

# `pixy_default` only runs when the input has at least one ALT="." record.
if bcftools view -H "$VCF" 2>/dev/null \
    | awk -F'\t' '$5 == "." { found = 1 } END { exit !found }'; then
  cmds+=(
    "python3 -m pixy --stats pi dxy --vcf \"$VCF\" --populations \"$POPS\" --window_size \"$WINDOW\" --output_folder \"$OUT/pixy_default\" --n_cores $PER_JOB"
  )
else
  echo "run_pixy: $VCF has no invariant sites (ALT=\".\"); skipping pixy_default" >&2
fi

cmds+=(
  "python3 -m pixy --stats fst --fst_type hudson --vcf \"$VCF\" --populations \"$POPS\" --window_size \"$WINDOW\" --bypass_invariant_check --output_folder \"$OUT/pixy_fst_hudson\" --n_cores $PER_JOB"
  "python3 -m pixy --stats fst --fst_type wc --vcf \"$VCF\" --populations \"$POPS\" --window_size \"$WINDOW\" --bypass_invariant_check --output_folder \"$OUT/pixy_fst_wc\" --n_cores $PER_JOB"
)

printf '%s\n' "${cmds[@]}" | python3 "$GOLD_DIR/scripts/run_parallel.py" "$JOBS"
