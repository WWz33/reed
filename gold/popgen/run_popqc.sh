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

THREADS="${THREADS:-4}"
FAIL=0

# ---- reed: stratified + cohort runs ------------------------------------
"$REED" stats "$VCF" --pop "$POPS" --force-samples \
    --freq --missing --singletons -o "$OUT/strat" --verbosity 0 --threads "$THREADS"
"$REED" stats "$VCF" --freq --missing --singletons \
    -o "$OUT/cohort" --verbosity 0 --threads "$THREADS"

# Cohort invariance (doc_v14 02 §2.2): --pop must not touch the cohort
# files. .singletons gains the POP column by design, so its cohort check is
# the first five columns.
for suffix in .afreq .vmiss; do
    if cmp -s "$OUT/cohort$suffix" "$OUT/strat$suffix"; then
        echo "cohort$suffix unchanged: SAME"
    else
        echo "cohort$suffix unchanged: DIFF" >&2; FAIL=1
    fi
done
if cmp -s <(cut -f1-5 "$OUT/cohort.singletons") \
          <(cut -f1-5 "$OUT/strat.singletons"); then
    echo "cohort.singletons columns 1-5 unchanged: SAME"
else
    echo "cohort.singletons columns 1-5 unchanged: DIFF" >&2; FAIL=1
fi

# ---- oracle: per-population plink2 runs ---------------------------------
# The header sample list drives the cut; populations resolve per group.
SAMPLES=$("$REED" stats "$VCF" --freq -o "$OUT/hdr" --verbosity 0 --threads 1 >/dev/null 2>&1; true)
# Read the sample names from the VCF header line (last #CHROM line).
mapfile -t NAMES < <(gzip -dc "$VCF" 2>/dev/null | grep "^#CHROM" | tail -1 | cut -f10-)
if [[ ${#NAMES[@]} -eq 0 ]]; then
    # plain (uncompressed) VCF
    mapfile -t NAMES < <(grep "^#CHROM" "$VCF" | tail -1 | cut -f10-)
fi

for GROUP in $(awk -F'\t' '{print $2}' "$POPS" | awk '!seen[$0]++'); do
    # The group's samples, as a keep-list file for bcftools.
    KEEP="$OUT/group_$GROUP.samples"
    awk -F'\t' -v g="$GROUP" '$2==g {print $1}' "$POPS" > "$KEEP"
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
    if cmp -s "$OUT/reed_$GROUP.afreq.vals" "$OUT/p2_$GROUP.afreq.vals"; then
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

# The singletons rows: the POP column must name the carrier's group and the
# first five columns must equal the cohort run's rows (counting runs over
# the whole cohort, doc_v14 02 §2.2).
if [[ $(wc -l < "$OUT/strat.singletons") -gt 1 ]]; then
    if awk -F'\t' 'NR>1 {print $6}' "$OUT/strat.singletons" | \
       grep -qvE '^(popA|popB|pop1|pop2|pop3|\.)$'; then
        echo "singletons POP column: DIFF (unexpected names)" >&2; FAIL=1
    else
        echo "singletons POP column: OK"
    fi
fi

if [[ $FAIL -eq 0 ]]; then
    echo "popqc parity: PASS"
else
    echo "popqc parity: FAIL" >&2
    exit 1
fi
