#!/usr/bin/env bash
# gold/popgen/run_bcftools_af.sh — `-q/-Q` allele-frequency parity against
# bcftools (docs/rounds/v5/02-cli.md §3.2): all five `:type` suffixes
# (nref, alt1, minor, major, nonmajor) at several thresholds, compared as
# site sets (the same CHROM+POS must survive on both sides).
#
# The 40 cases are independent, so they run concurrently through
# gold/scripts/run_parallel.py (doc_v22 R22-7). Each case uses one reed
# thread; parallelism comes from the case count. JOBS/REED_JOBS bounds it.
#
# Usage:
#   bash gold/popgen/run_bcftools_af.sh [VCF] [OUTDIR]
#
# Defaults to the local subset callset (`subvcf.vcf.gz`, 3,000 records x
# 348 samples). Pass the full `clean.vcf.gz` explicitly for a full-data run.
#
# Gate: every (threshold, type) pair has a zero symmetric difference.
set -euo pipefail

GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}"
REED="${REED:-/tmp/reed-target/release/reed}"
BCFTOOLS="${BCFTOOLS:-$HOME/biosoft/bcftools}"
SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

# One case: run both sides, compare the site sets.
run_case() {
    local flag="$1" thr="$2" type="$3"
    local tag="${flag}_${thr}_${type}"
    local b="$OUT/bcftools.$tag.txt" r="$OUT/reed.$tag.txt" v="$OUT/reed.$tag.vcf"
    "$BCFTOOLS" view -H "-$flag" "$thr:$type" "$VCF" 2>/dev/null \
        | cut -f1,2 | sort > "$b"
    "$REED" filter "$VCF" "-$flag" "$thr:$type" -Ov -o "$v" \
        --verbosity 0 --threads 1 2>/dev/null
    # `grep -v` exits 1 when every record was filtered out; that is a valid
    # case (the extreme thresholds keep zero sites), not a failure.
    { grep -v '^#' "$v" || true; } | cut -f1,2 | sort > "$r"
    local n_b n_r
    n_b=$(wc -l < "$b")
    n_r=$(wc -l < "$r")
    if cmp -s "$b" "$r"; then
        printf 'PASS -%s %s:%s (%s sites)\n' "$flag" "$thr" "$type" "$n_b"
    else
        printf 'FAIL -%s %s:%s (bcftools %s vs reed %s)\n' \
            "$flag" "$thr" "$type" "$n_b" "$n_r" >&2
        return 1
    fi
}

# Re-exec entry point used by the parallel driver: the parent exports
# REED/BCFTOOLS/VCF/OUT, so the child needs no positional arguments.
if [[ "${1:-}" == "--run-case" ]]; then
    VCF="${VCF:?--run-case needs VCF in the environment}"
    OUT="${OUT:?--run-case needs OUT in the environment}"
    run_case "$2" "$3" "$4"
    exit $?
fi

VCF="${1:-$REPO_ROOT/../subvcf.vcf.gz}"
OUT="${2:-$SCRATCH/out/bcftools_af}"

[[ -x "$REED" ]] || { echo "reed binary not executable: $REED" >&2; exit 1; }
[[ -x "$BCFTOOLS" ]] || { echo "bcftools not executable: $BCFTOOLS" >&2; exit 1; }
[[ -f "$VCF" ]] || { echo "input not found: $VCF" >&2; exit 1; }
mkdir -p "$OUT"
export REED BCFTOOLS VCF OUT

JOBS="${JOBS:-${REED_JOBS:-$(nproc)}}"
cmds=()
for flag in q Q; do
    for thr in 0.05 0.2 0.5 0.95; do
        for type in nref alt1 minor major nonmajor; do
            cmds+=("bash \"$SCRIPT\" --run-case $flag $thr $type")
        done
    done
done

printf '%s\n' "${cmds[@]}" | python3 "$GOLD_DIR/scripts/run_parallel.py" "$JOBS"

echo "AF-type parity: PASS (all -q/-Q thresholds and :types)"
