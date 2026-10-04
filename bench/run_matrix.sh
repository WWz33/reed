#!/bin/bash
R=/tmp/bench-target/release
D=/sri/home/sri2025201067/reed
echo "=== PHASE A: clean.bcf x3 ==="
for round in 1 2 3; do
  for bin in noo hts; do
    echo "--- round=$round bin=$bin file=clean.bcf ---"
    $R/bench-$bin $D/clean.bcf 2>/dev/null
  done
done
echo "=== PHASE B: clean.vcf.gz x1 ==="
for bin in noo hts; do
  echo "--- bin=$bin file=clean.vcf.gz ---"
  $R/bench-$bin $D/clean.vcf.gz 2>/dev/null
done
echo "=== PHASE C: clean_bg.vcf.gz x2 more ==="
for round in 2 3; do
  for bin in noo hts; do
    echo "--- round=$round bin=$bin file=clean_bg.vcf.gz ---"
    $R/bench-$bin $D/clean_bg.vcf.gz 2>/dev/null
  done
done
echo BENCH-DONE
