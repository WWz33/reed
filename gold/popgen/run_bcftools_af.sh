#!/usr/bin/env bash
# gold/popgen/run_bcftools_af.sh — `-q/-Q` allele-frequency parity against
# bcftools (docs/rounds/v5/02-cli.md §3.2): all five `:type` suffixes
# (nref, alt1, minor, major, nonmajor) at several thresholds, compared as
# site sets (the same CHROM+POS must survive on both sides).
#
# Usage:
#   bash gold/popgen/run_bcftools_af.sh [VCF] [OUTDIR]
#
# Defaults to the Chr01:1-2,000,000 region of the gold subset (8,364 sites,
# enough boundary cases to exercise the f32/double comparison and the INFO
# AC/AN preference).
#
# Gate: every (threshold, type) pair has a zero symmetric difference.
set -euo pipefail

GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}"
REED="${REED:-/tmp/reed-target/release/reed}"
BCFTOOLS="${BCFTOOLS:-$HOME/biosoft/bcftools}"

VCF="${1:-$SCRATCH/out/ld_prune_sub/region.vcf.gz}"
OUT="${2:-$SCRATCH/out/bcftools_af}"

[[ -x "$REED" ]] || { echo "reed binary not executable: $REED" >&2; exit 1; }
[[ -x "$BCFTOOLS" ]] || { echo "bcftools not executable: $BCFTOOLS" >&2; exit 1; }
[[ -f "$VCF" ]] || { echo "input not found: $VCF" >&2; exit 1; }
mkdir -p "$OUT"

fail=0
for flag in q Q; do
    for thr in 0.05 0.2 0.5 0.95; do
        for type in nref alt1 minor major nonmajor; do
            # -Q above 0.5 keeps almost everything; -q at 0.95 nearly nothing.
            "$BCFTOOLS" view -H "-$flag" "$thr:$type" "$VCF" 2>/dev/null \
                | cut -f1,2 | sort > "$OUT/bcftools.txt"
            "$REED" filter "$VCF" "-$flag" "$thr:$type" -Ov -o "$OUT/reed.vcf" \
                --verbosity 0 2>/dev/null
            grep -v '^#' "$OUT/reed.vcf" | cut -f1,2 | sort > "$OUT/reed.txt"
            n_b=$(wc -l < "$OUT/bcftools.txt")
            n_r=$(wc -l < "$OUT/reed.txt")
            if cmp -s "$OUT/bcftools.txt" "$OUT/reed.txt"; then
                printf 'PASS -%s %s:%s (%s sites)\n' "$flag" "$thr" "$type" "$n_b"
            else
                printf 'FAIL -%s %s:%s (bcftools %s vs reed %s)\n' \
                    "$flag" "$thr" "$type" "$n_b" "$n_r" >&2
                fail=1
            fi
        done
    done
done

[[ "$fail" -eq 0 ]] || { echo "AF-type parity: FAIL" >&2; exit 1; }
echo "AF-type parity: PASS (all -q/-Q thresholds and :types)"
