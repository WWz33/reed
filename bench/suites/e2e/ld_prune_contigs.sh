#!/bin/bash
# doc_v15 D15-17 scenario: `--indep-pairwise` over the same 300k sites laid
# out as one contig vs 10,000 contigs. The pass-1 keep-set lookup used to
# scan the whole contig list per surviving record (linear); it now
# binary-searches, so the many-contig file must stay in the same time class.
#
# Historical reference (2026-10-18, 56-core host, --threads 1):
#   one contig   0.36 s
#   10k contigs  8.71 s   (before the fix — ~24x)
#   after the fix the 10k-contig run is ~0.8 s (~2.2x)
#
# Usage: bench/suites/e2e/ld_prune_contigs.sh [REED_BINARY]
set -euo pipefail

REED="${1:-${REED:-/tmp/reed-target/release/reed}}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

python3 - "$TMP" <<'PY'
import random, sys
out = sys.argv[1]
random.seed(7)

def rec(chrom, pos):
    gts = "\t".join("0/1" if random.random() < 0.5 else "0/0" for _ in range(4))
    return f"{chrom}\t{pos}\t.\tA\tT\t50\t.\t.\tGT\t{gts}\n"

with open(f"{out}/one.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=c0,length=1000000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n")
    for i in range(300_000):
        f.write(rec("c0", i + 1))

with open(f"{out}/many.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    for c in range(10_000):
        f.write(f"##contig=<ID=c{c:05d},length=1000>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n")
    for c in range(10_000):
        for i in range(30):
            f.write(rec(f"c{c:05d}", i + 1))
PY

run() {
  local name="$1" file="$2" start end
  start=$(date +%s.%N)
  "$REED" filter "$TMP/$file" --indep-pairwise 50 5 0.999 -Oz -o /dev/null --threads 1 2>/dev/null
  end=$(date +%s.%N)
  echo "$name $(echo "$end - $start" | bc)"
}

echo "== reed: $REED"
run one-contig  one.vcf
run 10k-contigs many.vcf
echo "(doc_v15 D15-17: the 10k-contig time must stay in the one-contig class;"
echo " before the fix it was ~24x slower)"
