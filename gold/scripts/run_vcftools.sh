#!/usr/bin/env bash
# gold/scripts/run_vcftools.sh — vcftools 分支 (common 层)。
# vcftools 无法读取 INFO 字段 (QD/FS/MQ/SOR)，故只跑 common 层。
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_tools VCFTOOLS BGZIP TABIX

log "[vcftools] 开始 $(tool_versions | tr '\n' ' ')"
[[ -s "$OUTDIR/samples.keep.txt" ]] || build_sample_keep

# REGION=Chr01:1-2000000 -> --chr Chr01 --from-bp 1 --to-bp 2000000
REG_ARGS=()
if [[ -n "${REGION:-}" ]]; then
    chr="${REGION%%:*}"; rest="${REGION#*:}"; from="${rest%%-*}"; to="${rest#*-}"
    REG_ARGS=(--chr "$chr" --from-bp "$from" --to-bp "$to")
fi

# --max-missing 语义：0=全缺失允许，1=不允许缺失 => 1 - 缺失率上限
MAX_MISSING=$(awk -v m="$SITE_MISS_MAX" 'BEGIN{printf "%.10g", 1-m}')

FILT_ARGS=()
[[ "$REQUIRE_PASS" == "1" ]] && FILT_ARGS+=(--remove-filtered-all)
[[ "$KEEP_BIALLELIC_SNPS" == "1" ]] && FILT_ARGS+=(--remove-indels --min-alleles 2 --max-alleles 2)
[[ "$HWE_MIN" != "0" ]] && FILT_ARGS+=(--hwe "$HWE_MIN")

OUT="$OUTDIR/vcftools.vcf.gz"
log "[vcftools] 运行 (keep=$(wc -l < "$OUTDIR/samples.keep.txt") samples, max-missing=$MAX_MISSING, maf=$MAF_MIN, hwe=$HWE_MIN)"
"$VCFTOOLS" --gzvcf "$INPUT" \
    "${REG_ARGS[@]}" \
    --keep "$OUTDIR/samples.keep.txt" \
    "${FILT_ARGS[@]}" \
    --minQ "$MIN_QUAL" \
    --minDP "$GT_DP_MIN" --maxDP "$GT_DP_MAX" --minGQ "$GT_GQ_MIN" \
    --max-missing "$MAX_MISSING" \
    --maf "$MAF_MIN" \
    --recode --recode-INFO-all --stdout \
    > "$WORKDIR/vcftools.raw.vcf" 2> "$LOGDIR/vcftools.err"

grep -iE 'kept|Individuals|Sites' "$LOGDIR/vcftools.err" | tail -4 | sed 's/^/[vcftools] /' >&2

"$BGZIP" -c "$WORKDIR/vcftools.raw.vcf" > "$OUT"
"$TABIX" -f -p vcf "$OUT"
write_stats vcftools "$OUT"
log "[vcftools] 完成"
