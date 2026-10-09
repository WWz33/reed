#!/usr/bin/env bash
# gold/scripts/make_subset.sh — 切染色体子集，便于快速迭代（大文件全量跑之前先在小集上定型参数）。
#
# 用法:
#   REED_SCRATCH=target/gold bash gold/scripts/make_subset.sh Chr01,Chr02,Chr03 "$REED_SCRATCH/work/subset.Chr01-03.vcf.gz"
# 之后把 INPUT 指向子集:
#   INPUT=$REED_SCRATCH/work/subset.Chr01-03.vcf.gz bash gold/scripts/run_all.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_tools BCFTOOLS

CHROMS="${1:?用法: make_subset.sh CHROMS OUTFILE}"
OUT="${2:?用法: make_subset.sh CHROMS OUTFILE}"
mkdir -p "$(dirname "$OUT")"
log "切子集: $CHROMS -> $OUT (from $INPUT)"
"$BCFTOOLS" view -r "$CHROMS" "$INPUT" -Oz -o "$OUT"
"$BCFTOOLS" index -t "$OUT"
log "完成: $("$BCFTOOLS" view -H "$OUT" | wc -l) sites"
