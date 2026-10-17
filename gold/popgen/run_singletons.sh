#!/usr/bin/env bash
# gold/popgen/run_singletons.sh — `reed stats --singletons` vs vcftools
# (docs/rounds/v14/01-review.md §7). vcftools `--singletons` is the layout
# reference: CHROM POS SINGLETON/DOUBLETON ALLELE INDV.
#
# Usage: bash gold/popgen/run_singletons.sh [VCF] [OUTDIR]
#        bash gold/popgen/run_singletons.sh --multi    # multi-chunk gate
#
# `--multi` generates a 50 000-record fixture (about a dozen chunks) so the
# cross-chunk path is exercised: worker partials are merged by the ordered
# sink, and a single-chunk fixture never reaches that code.
#
# The fixture carries a singleton, a private doubleton and a multiallelic
# singleton so the per-allele semantics (including REF) are exercised.
set -euo pipefail
GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}/singletons"
REED="${REED:-/tmp/reed-target/release/reed}"
VCFTOOLS="${VCFTOOLS:-$HOME/biosoft/vcftools}"

MULTI=0
if [[ "${1:-}" == "--multi" ]]; then
    MULTI=1; shift
fi
VCF="${1:-$GOLD_DIR/fixtures/singletons.vcf}"
OUT="${2:-$SCRATCH}"
mkdir -p "$OUT"

if [[ $MULTI -eq 1 ]]; then
    VCF="$OUT/multi.vcf"
    awk 'BEGIN{
        print "##fileformat=VCFv4.2";
        print "##contig=<ID=chr1,length=1000000>";
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">";
        printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n";
        for (i = 0; i < 50000; i++) {
            k = i % 3;
            if (k == 0)      printf "chr1\t%d\t.\tA\tG\t50\t.\t.\tGT\t0/1\t0/0\t0/0\t0/0\n", i+1;
            else if (k == 1) printf "chr1\t%d\t.\tA\tG\t50\t.\t.\tGT\t1/1\t0/0\t0/0\t0/0\n", i+1;
            else             printf "chr1\t%d\t.\tA\tG,T\t50\t.\t.\tGT\t0/2\t0/0\t0/0\t0/0\n", i+1;
        }
    }' > "$VCF"
fi

[[ -x "$REED" ]] || { echo "reed binary not found: $REED" >&2; exit 1; }
[[ -x "$VCFTOOLS" ]] || { echo "vcftools not found: $VCFTOOLS" >&2; exit 1; }

"$REED" stats "$VCF" --singletons -o "$OUT/reed" --verbosity 0 --threads "${THREADS:-4}"
"$VCFTOOLS" --vcf "$VCF" --singletons --out "$OUT/vcftools" >/dev/null 2>&1

if diff -u "$OUT/vcftools.singletons" "$OUT/reed.singletons"; then
    echo "singletons parity: PASS (byte-equal)"
else
    echo "singletons parity: FAIL" >&2
    exit 1
fi
