#!/usr/bin/env bash
# gold/scripts/run_reed.sh — reed 分支。
# 与 run_bcftools.sh 使用同一份 filters.conf 和 run_all.sh 统一算出的
# samples.keep.txt；02-cli.md §6 的等价命令由本脚本生成（无 TOML/conf 输入）。
#
# 产出：
#   out/reed.vcf.gz   common + INFO 层 = 与 bcftools.full 对拍的 reed 输出
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_tools REED BGZIP

log "[reed] 开始 $(tool_versions | tr '\n' ' ')"

[[ -s "$OUTDIR/samples.keep.txt" ]] || build_sample_keep

STALE_INFO="INFO/AC,INFO/AN,INFO/AF,INFO/MLEAC,INFO/MLEAF,INFO/ExcessHet,INFO/InbreedingCoeff,INFO/DP"

ARGS=()
[[ -n "${REGION:-}" ]] && ARGS+=(-r "$REGION")
[[ -s "$OUTDIR/samples.keep.txt" ]] && ARGS+=(-S "$OUTDIR/samples.keep.txt")
[[ "$KEEP_BIALLELIC_SNPS" == "1" ]] && ARGS+=(--snps-only=just-acgt -M 2)
[[ "$REQUIRE_PASS" == "1" ]] && ARGS+=(-f PASS)
ARGS+=(-i "QUAL>=$MIN_QUAL")
ARGS+=(-e "QD<$INFO_QD_MIN || FS>$INFO_FS_MAX || MQ<$INFO_MQ_MIN || SOR>$INFO_SOR_MAX")
ARGS+=(--min-dp "$GT_DP_MIN" --max-dp "$GT_DP_MAX" --min-gq "$GT_GQ_MIN")
# doc_v4 02-cli-semantics.md 回退版 §6: reed takes the missing-rate spelling
# directly, like run_plink2.sh — `SITE_MISS_MAX` is the spec's own value.
ARGS+=(--geno "$SITE_MISS_MAX" --maf "$MAF_MIN")
[[ "$HWE_MIN" != "0" ]] && ARGS+=(--hwe "$HWE_MIN")
ARGS+=(-x "$STALE_INFO")
ARGS+=(--fill-tags AC,AN,NS,MAF,HWE)
ARGS+=(-Oz -W -o "$OUTDIR/reed.vcf.gz" --threads "$THREADS")

log "[reed] ${REED} filter ${INPUT} ${ARGS[*]}"
"$REED" filter "$INPUT" "${ARGS[@]}" 2> "$LOGDIR/reed.err" \
    || { log "[reed] 失败，见 $LOGDIR/reed.err"; exit 1; }
log "[reed] 完成：$OUTDIR/reed.vcf.gz"
