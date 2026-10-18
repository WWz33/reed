#!/usr/bin/env bash
# gold/popgen/run_popqc.sh — `reed stats --pop` stratified QC vs the oracle
# (docs/rounds/v14/01-review.md §7, C3-d). There is no single reference
# command for stratified QC: vcftools' --singletons has no POP column and
# plink2 --within only tags .smiss with CATPHENO. The oracle is therefore a
# per-population run: cut each population's samples out of the VCF, run
# plink2 --freq/--missing (and reed's own cohort mode) on it, and compare
# the stratified rows group by group. The cohort files must stay
# byte-identical to a run without --pop.
#
# Usage: bash gold/popgen/run_popqc.sh [VCF] [POPFILE] [OUTDIR]
#   VCF     default gold/fixtures/singletons.vcf (the subset works too:
#           bash run_popqc.sh target/gold/work/subset.Chr01-03.vcf.gz \
#                ../../clean.3pop.txt)
#   POPFILE sample<TAB>population rows, no header
set -euo pipefail
GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}/popqc"
REED="${REED:-/tmp/reed-target/release/reed}"
PLINK2="${PLINK2:-$HOME/biosoft/plink2}"

VCF="${1:-$GOLD_DIR/fixtures/singletons.vcf}"
if [[ -n "${2:-}" ]]; then
    POPS="$2"
else
    POPS="$SCRATCH/pops.tsv"
    mkdir -p "$SCRATCH"
    printf 's0\tpopA\ns1\tpopA\ns2\tpopB\n' > "$POPS"   # s3 stays unassigned
fi
OUT="${3:-$SCRATCH}"
mkdir -p "$OUT"

[[ -x "$REED" ]] || { echo "reed binary not found: $REED" >&2; exit 1; }
[[ -x "$PLINK2" ]] || { echo "plink2 not found: $PLINK2" >&2; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo "bcftools not found on PATH" >&2; exit 1; }

THREADS="${THREADS:-4}"
FAIL=0

# ---- reed: stratified + cohort runs ------------------------------------
"$REED" stats "$VCF" --pop "$POPS" --force-samples \
    --freq --missing --singletons -o "$OUT/strat" --verbosity 0 --threads "$THREADS"
"$REED" stats "$VCF" --freq --missing --singletons \
    -o "$OUT/cohort" --verbosity 0 --threads "$THREADS"

# Cohort invariance (doc_v14 02 §2.2): --pop must not touch the cohort
# files. .singletons gains the POP column by design, so its cohort check is
# the first five columns. A missing file is a failure, not an empty match
# (two absent files used to compare equal).
for suffix in .afreq .vmiss .singletons; do
    [[ -s "$OUT/strat$suffix" && -s "$OUT/cohort$suffix" ]] || {
        echo "missing output: $OUT/{strat,cohort}$suffix" >&2; FAIL=1; continue; }
done
for suffix in .afreq .vmiss; do
    if cmp -s "$OUT/cohort$suffix" "$OUT/strat$suffix"; then
        echo "cohort$suffix unchanged: SAME"
    else
        echo "cohort$suffix unchanged: DIFF" >&2; FAIL=1
    fi
done
if [[ -f "$OUT/cohort.singletons" && -f "$OUT/strat.singletons" ]]; then
    if cmp -s <(cut -f1-5 "$OUT/cohort.singletons") \
              <(cut -f1-5 "$OUT/strat.singletons"); then
        echo "cohort.singletons columns 1-5 unchanged: SAME"
    else
        echo "cohort.singletons columns 1-5 unchanged: DIFF" >&2; FAIL=1
    fi
fi

# ---- oracle: per-population plink2 runs ---------------------------------
# The group's samples are cut down to the ones the VCF carries: a popfile
# entry naming a sample the input lacks is legal for reed (--force-samples)
# but `bcftools view -S` refuses it, which used to abort the whole gate.
# `cut -f10-` keeps the sample columns on one line; split them into one
# array entry each, or every comparison below reads the whole list as a
# single name (which silently skipped every group).
mapfile -t NAMES < <(gzip -dc "$VCF" 2>/dev/null | grep "^#CHROM" | tail -1 | cut -f10- | tr '\t' '\n')
if [[ ${#NAMES[@]} -eq 0 ]]; then
    # plain (uncompressed) VCF
    mapfile -t NAMES < <(grep "^#CHROM" "$VCF" | tail -1 | cut -f10- | tr '\t' '\n')
fi
[[ ${#NAMES[@]} -gt 0 ]] || { echo "no samples in $VCF" >&2; exit 1; }
in_header() {
    local want="$1" name
    for name in "${NAMES[@]}"; do
        [[ "$name" == "$want" ]] && return 0
    done
    return 1
}

for GROUP in $(awk -F'\t' '{print $2}' "$POPS" | awk '!seen[$0]++'); do
    # The group's samples, as a keep-list file for bcftools.
    KEEP="$OUT/group_$GROUP.samples"
    : > "$KEEP"
    while IFS=$'\t' read -r sample group; do
        [[ "$group" == "$GROUP" ]] || continue
        in_header "$sample" && printf '%s\n' "$sample" >> "$KEEP"
    done < "$POPS"
    if [[ ! -s "$KEEP" ]]; then
        echo "afreq.pop $GROUP: no sample of the group is in the VCF (skipped)"
        continue
    fi
    GVCF="$OUT/group_$GROUP.vcf.gz"
    bcftools view -S "$KEEP" -Oz -o "$GVCF" "$VCF" 2>/dev/null
    bcftools index -t "$GVCF" 2>/dev/null
    "$PLINK2" --vcf "$GVCF" --freq --missing --out "$OUT/p2_$GROUP" \
        >/dev/null 2>&1

    # reed's rows for this group: strip the POP column to recover the
    # cohort layout, then compare the value columns (REF ALT ALT_FREQS
    # OBS_CT). reed collapses multiallelic ALT into one frequency (the
    # cohort .afreq convention, `CollapseAlt`); plink2 prints one frequency
    # per ALT — the per-ALT values are summed before the comparison, and
    # the collapsed-vs-split difference is the only accepted transform
    # (registered with the pop stratification, doc_v14 02 §2.2).
    awk -F'\t' -v g="$GROUP" 'NR>1 && $3==g {
        n = split($6, f, ",");
        s = 0; for (i = 1; i <= n; i++) s += f[i];
        print $4"\t"$5"\t"s"\t"$7
    }' "$OUT/strat.afreq.pop" > "$OUT/reed_$GROUP.afreq.vals"
    awk 'NR>1 {
        n = split($5, f, ",");
        s = 0; for (i = 1; i <= n; i++) s += f[i];
        print $3"\t"$4"\t"s"\t"$6
    }' "$OUT/p2_$GROUP.afreq" > "$OUT/p2_$GROUP.afreq.vals"
    # Numeric, not byte-wise: plink2 prints 6 significant digits, so summing
    # its per-ALT frequencies can land one ulp below reed's collapsed value
    # (0.333333+0.333333 = 0.666666 vs 0.666667).
    if python3 "$GOLD_DIR/popgen/compare_numeric.py" \
        "$OUT/reed_$GROUP.afreq.vals" "$OUT/p2_$GROUP.afreq.vals"; then
        echo "afreq.pop $GROUP vs plink2: SAME"
    else
        echo "afreq.pop $GROUP vs plink2: DIFF" >&2; FAIL=1
        diff "$OUT/reed_$GROUP.afreq.vals" "$OUT/p2_$GROUP.afreq.vals" | head -4 >&2
    fi

    awk -F'\t' -v g="$GROUP" 'NR>1 && $3==g {print $4"\t"$5"\t"$6}' \
        "$OUT/strat.vmiss.pop" > "$OUT/reed_$GROUP.vmiss.vals"
    awk 'NR>1 {print $3"\t"$4"\t"$5}' \
        "$OUT/p2_$GROUP.vmiss" > "$OUT/p2_$GROUP.vmiss.vals"
    if cmp -s "$OUT/reed_$GROUP.vmiss.vals" "$OUT/p2_$GROUP.vmiss.vals"; then
        echo "vmiss.pop $GROUP vs plink2: SAME"
    else
        echo "vmiss.pop $GROUP vs plink2: DIFF" >&2; FAIL=1
        diff "$OUT/reed_$GROUP.vmiss.vals" "$OUT/p2_$GROUP.vmiss.vals" | head -4 >&2
    fi
done

# The singletons POP column must name the carrier's population, resolved
# from the popfile (`.` for a carrier outside it) — a name whitelist used to
# pass any of five hardcoded strings and fail any other legitimate popfile.
rows=0
[[ -f "$OUT/strat.singletons" ]] && rows=$(($(wc -l < "$OUT/strat.singletons") - 1))
if [[ $rows -le 0 ]]; then
    echo "singletons: no rows on this input (nothing to check)"
else
    if awk -F'\t' 'NR == FNR { pop[$1] = $2; next }
        FNR > 1 {
            want = ($5 in pop) ? pop[$5] : "."
            if ($6 != want) {
                printf "row %d: carrier %s is in %s, the POP column says %s\n", FNR, $5, want, $6
                bad = 1
            }
        }
        END { exit bad }' "$POPS" "$OUT/strat.singletons"; then
        echo "singletons POP column: OK ($rows rows)"
    else
        echo "singletons POP column: DIFF" >&2; FAIL=1
    fi
fi

if [[ $FAIL -eq 0 ]]; then
    echo "popqc parity: PASS"
else
    echo "popqc parity: FAIL" >&2
    exit 1
fi
