#!/usr/bin/env bash
# gold/popgen/run_hist.sh — `reed stats --hist` vs bcftools stats sections
# (docs/rounds/v14/01-review.md §7, C3-c). The built-in metrics are compared
# on their non-zero bins: the QUAL section (SNP+indel merged, iqual bins),
# the DP section (per-genotype and per-site columns) and the IDD section
# (signed lengths).
#
# Usage: bash gold/popgen/run_hist.sh [VCF] [OUTDIR]
set -euo pipefail
GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}/hist"
REED="${REED:-/tmp/reed-target/release/reed}"
BCFTOOLS="${BCFTOOLS:-$HOME/biosoft/bcftools}"

VCF="${1:-$GOLD_DIR/fixtures/hist.vcf}"
OUT="${2:-$SCRATCH}"
mkdir -p "$OUT"

[[ -x "$REED" ]] || { echo "reed binary not found: $REED" >&2; exit 1; }
[[ -x "$BCFTOOLS" ]] || { echo "bcftools not found: $BCFTOOLS" >&2; exit 1; }

"$REED" stats "$VCF" --hist qual,dp,dp-site,idl,af -o "$OUT/reed" --verbosity 0
# `-s -` adds the per-genotype (FMT/DP) rows to the DP section; without it
# only the per-site (INFO/DP) rows carry counts.
"$BCFTOOLS" stats -s - "$VCF" > "$OUT/bcftools.stats" 2>/dev/null

fail=0
# reed .hist -> "key<TAB>count", non-zero rows only.
reed_bins() { awk -F'\t' '$3+0 > 0 { print $1"\t"$3 }' "$1" | sort; }

check() {
    if diff -q "$1" "$2" >/dev/null; then
        echo "  SAME  $3"
    else
        echo "  DIFF  $3" >&2; fail=1
    fi
}

# QUAL: reed merges SNP+indel and labels bins in QUAL units
# (`iqual = 1 + qual·10`, so bin lower = quality). bcftools prints the
# quality in $3 (a `.` row is its missing bucket: reed counts those in
# `#missing` instead) and the counts as nsnp ($4) and nindels ($7). Only the
# exact region (iqual <= 10^5, 0.1-wide bins) compares bin by bin: above it
# bcftools widens the bins and prints a rounded float that cannot recover
# the exact bin edge, so the tails are compared by total instead.
awk -F'\t' '$1=="QUAL" && $3 != "." && int($3*10+0.5)+1 <= 100001 { printf "%.6f\t%d\n", $3, $4+$7 }' \
    "$OUT/bcftools.stats" | sort > "$OUT/qual.bcf"
awk -F'\t' '$3+0>0 && $1 <= 10000.0 { printf "%.6f\t%d\n", $1, $3 }' \
    "$OUT/reed.qual.hist" | sort > "$OUT/qual.reed"
check "$OUT/qual.bcf" "$OUT/qual.reed" "QUAL exact region (0.1 bins)"
awk -F'\t' '$1=="QUAL" && $3 != "." && int($3*10+0.5)+1 > 100001 { s += $4+$7 } END { print s+0 }' \
    "$OUT/bcftools.stats" > "$OUT/qual.tail.bcf"
awk -F'\t' '$3+0>0 && $1 > 10000.0 { s += $3 } END { print s+0 }' "$OUT/reed.qual.hist" > "$OUT/qual.tail.reed"
check "$OUT/qual.tail.bcf" "$OUT/qual.tail.reed" "QUAL widened tail (total)"

# DP: per-genotype ($4) and per-site ($6). reed labels the overflow bucket
# `>max` exactly like bcftools, so the rows compare as they are.
awk -F'\t' '$1=="DP" && $4+0>0 { print $3"\t"$4 }' "$OUT/bcftools.stats" | sort > "$OUT/dp.bcf"
reed_bins "$OUT/reed.dp.hist" > "$OUT/dp.reed"
check "$OUT/dp.bcf" "$OUT/dp.reed" "DP per genotype"
awk -F'\t' '$1=="DP" && $6+0>0 { print $3"\t"$6 }' "$OUT/bcftools.stats" | sort > "$OUT/dpsite.bcf"
reed_bins "$OUT/reed.dp-site.hist" > "$OUT/dpsite.reed"
check "$OUT/dpsite.bcf" "$OUT/dpsite.reed" "DP per site"

# IDD: signed length ($3) -> number of sites ($4).
awk -F'\t' '$1=="IDD" { print $3"\t"$4 }' "$OUT/bcftools.stats" | sort > "$OUT/idd.bcf"
reed_bins "$OUT/reed.idl.hist" > "$OUT/idd.reed"
check "$OUT/idd.bcf" "$OUT/idd.reed" "IDD (indel length)"

# AF: bcftools prints the frequency (`(i−1)/(m_af−1)`, singletons folded
# into the first row) in $3 and the SNP ($4) and indel ($7) counts; reed's
# one count column sums both.
awk -F'\t' '$1=="AF" && ($4+$7)+0>0 { printf "%.6f\t%d\n", $3, $4+$7 }' \
    "$OUT/bcftools.stats" | sort > "$OUT/af.bcf"
awk -F'\t' '$3+0>0 { printf "%.6f\t%d\n", $1, $3 }' "$OUT/reed.af.hist" | sort > "$OUT/af.reed"
check "$OUT/af.bcf" "$OUT/af.reed" "AF bins (snps+indels)"

if [[ $fail -eq 0 ]]; then
    echo "hist parity: PASS"
else
    echo "hist parity: FAIL" >&2
    exit 1
fi
