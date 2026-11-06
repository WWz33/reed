#!/usr/bin/env bash
# Runs on Linux/macOS and inside the old-glibc compatibility container.
set -euo pipefail

binary=$1
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
"$binary" --version
"$binary" --help > "$scratch/help.txt"
# Keep a real read/write/compute path in the gate, not just dynamic loading.
printf '%s\n' \
    '##fileformat=VCFv4.2' \
    '##contig=<ID=chr1,length=1000>' \
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
    > "$scratch/in.vcf"
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n' >> "$scratch/in.vcf"
printf 'chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\n' >> "$scratch/in.vcf"
printf 'chr1\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t1/1\n' >> "$scratch/in.vcf"
printf 'chr1\t3\t.\tG\tA\t50\tPASS\t.\tGT\t1/1\t0/1\t0/0\t0/1\n' >> "$scratch/in.vcf"

"$binary" filter "$scratch/in.vcf" --threads 1 --no-version -Oz -W -o "$scratch/out.vcf.gz"
test -s "$scratch/out.vcf.gz.csi"
"$binary" filter "$scratch/out.vcf.gz" --threads 4 --no-version -Ov -o "$scratch/roundtrip.vcf"
grep -v '^#' "$scratch/in.vcf" > "$scratch/expected.txt"
grep -v '^#' "$scratch/roundtrip.vcf" > "$scratch/actual.txt"
diff -u "$scratch/expected.txt" "$scratch/actual.txt"
# Exercise region/index lookup and native BCF decoding as well.
"$binary" filter "$scratch/out.vcf.gz" -r chr1:2-2 --threads 1 --no-version -Ob -o "$scratch/region.bcf"
"$binary" filter "$scratch/region.bcf" --threads 4 --no-version -Ov -o "$scratch/region.vcf"
grep -v '^#' "$scratch/region.vcf" > "$scratch/region.txt"
sed -n '2p' "$scratch/expected.txt" > "$scratch/one.txt"
diff -u "$scratch/one.txt" "$scratch/region.txt"
# Dense PCA reaches faer's CPU dispatch, not only reed's CLI parser.
"$binary" pca "$scratch/out.vcf.gz" --solver dense --missing meanimpute --pcs 2 --threads 1 -o "$scratch/pcs"
test -s "$scratch/pcs.eigenvec"
test "$(wc -l < "$scratch/pcs.eigenval")" -eq 2
printf 'release smoke: PASS\n'
