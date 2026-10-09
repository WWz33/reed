#!/usr/bin/env bash
# gold/scripts/run_bcftools.sh — bcftools 分支。
# 产出：
#   out/bcftools.common.vcf.gz  common 层 (与 vcftools/plink2 可对拍)
#   out/bcftools.full.vcf.gz    common + INFO 层 = **黄金标准 gold**
#
# 关键实现要点 (详见 README「踩坑记录」)：
#   1. 输入 INFO/AC,AN,AF,MLEAC,MLEAF 是陈旧的 (348 样本是大队列子集、未重算)，
#      必须在 fill-tags 之前剔除，否则 MAF/HWE 会用错值。
#   2. bcftools 表达式里 FMT 字段上的 `||` 有 bug，必须用按位 `|` 连接失败条件。
#   3. 位点级 `&&` 正常。
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_tools BCFTOOLS BGZIP

log "[bcftools] 开始 $(tool_versions | tr '\n' ' ')"

[[ -s "$OUTDIR/samples.keep.txt" ]] || build_sample_keep

VIEW_ARGS=(-m2 -M2 -v snps -i "QUAL>=$MIN_QUAL")
[[ "$REQUIRE_PASS" == "1" ]] && VIEW_ARGS+=(-f PASS)
[[ -n "${REGION:-}" ]] && VIEW_ARGS+=(-r "$REGION")
[[ -s "$OUTDIR/samples.keep.txt" ]] && VIEW_ARGS+=(-S "$OUTDIR/samples.keep.txt")

STALE_INFO="INFO/AC,INFO/AN,INFO/AF,INFO/MLEAC,INFO/MLEAF,INFO/ExcessHet,INFO/InbreedingCoeff,INFO/DP"
MASK_EXPR="FMT/DP<$GT_DP_MIN | FMT/DP>$GT_DP_MAX | FMT/GQ<$GT_GQ_MIN"
SITE_EXPR="F_MISSING<=$SITE_MISS_MAX && MAF>=$MAF_MIN && HWE>=$HWE_MIN"
INFO_EXPR="QD>=$INFO_QD_MIN && FS<=$INFO_FS_MAX && MQ>=$INFO_MQ_MIN && SOR<=$INFO_SOR_MAX"

TAGGED="$WORKDIR/bcftools.tagged.bcf"
log "[bcftools] pass1: type/qual/pass -> strip stale INFO -> mask GT -> fill-tags"
"$BCFTOOLS" view --threads "$THREADS" "${VIEW_ARGS[@]}" "$INPUT" -Ou 2> "$LOGDIR/bcftools.err" \
    | "$BCFTOOLS" annotate -x "$STALE_INFO" -Ou 2>> "$LOGDIR/bcftools.err" \
    | "$BCFTOOLS" filter -e "$MASK_EXPR" -S . -Ou 2>> "$LOGDIR/bcftools.err" \
    | "$BCFTOOLS" +fill-tags -Ou -- -t AC,AN,NS,MAF,HWE 2>> "$LOGDIR/bcftools.err" \
    > "$TAGGED"

log "[bcftools] pass2a: common 层 (site: $SITE_EXPR)"
"$BCFTOOLS" view -i "$SITE_EXPR" "$TAGGED" -Oz -o "$OUTDIR/bcftools.common.vcf.gz" 2>> "$LOGDIR/bcftools.err"
"$BCFTOOLS" index -t "$OUTDIR/bcftools.common.vcf.gz" 2>> "$LOGDIR/bcftools.err"

log "[bcftools] pass2b: full 层 (+INFO: $INFO_EXPR)"
"$BCFTOOLS" view -i "$SITE_EXPR && $INFO_EXPR" "$TAGGED" -Oz -o "$OUTDIR/bcftools.full.vcf.gz" 2>> "$LOGDIR/bcftools.err"
"$BCFTOOLS" index -t "$OUTDIR/bcftools.full.vcf.gz" 2>> "$LOGDIR/bcftools.err"

write_stats bcftools "$OUTDIR/bcftools.common.vcf.gz" \
  "\"common_sites\": $(bcftools view -H "$OUTDIR/bcftools.common.vcf.gz" 2>/dev/null | wc -l), \"full_sites\": $("$BCFTOOLS" view -H "$OUTDIR/bcftools.full.vcf.gz" 2>/dev/null | wc -l)"
log "[bcftools] 完成"
