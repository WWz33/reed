#!/usr/bin/env bash
# gold/scripts/run_plink2.sh — plink2 分支 (common 层)。
# plink2 不读 INFO 字段，故只跑 common 层；输出回写为 VCF.gz。
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_tools PLINK2 BGZIP TABIX

log "[plink2] 开始 $(tool_versions | tr '\n' ' ')"
[[ -s "$OUTDIR/samples.keep.txt" ]] || build_sample_keep

REG_ARGS=()
if [[ -n "${REGION:-}" ]]; then
    chr="${REGION%%:*}"; rest="${REGION#*:}"; from="${rest%%-*}"; to="${rest#*-}"
    REG_ARGS=(--chr "$chr" --from-bp "$from" --to-bp "$to")
fi

TYPE_ARGS=()
[[ "$KEEP_BIALLELIC_SNPS" == "1" ]] && TYPE_ARGS+=(--snps-only just-acgt --min-alleles 2 --max-alleles 2)
PASS_ARGS=()
[[ "$REQUIRE_PASS" == "1" ]] && PASS_ARGS+=(--vcf-filter)

PREFIX="$WORKDIR/plink2"
log "[plink2] 运行 (geno=$SITE_MISS_MAX, maf=$MAF_MIN, hwe=$HWE_MIN)"
"$PLINK2" --vcf "$INPUT" \
    "${REG_ARGS[@]}" \
    --keep "$OUTDIR/samples.keep.txt" \
    "${PASS_ARGS[@]}" "${TYPE_ARGS[@]}" \
    --var-min-qual "$MIN_QUAL" \
    --vcf-min-dp "$GT_DP_MIN" --vcf-max-dp "$GT_DP_MAX" --vcf-min-gq "$GT_GQ_MIN" \
    --geno "$SITE_MISS_MAX" \
    --maf "$MAF_MIN" \
    --hwe "$HWE_MIN" \
    --threads "$THREADS" \
    --export vcf bgz --out "$PREFIX" \
    > "$LOGDIR/plink2.err" 2>&1

grep -iE 'remaining after main filters|variants removed|samples remaining' "$LOGDIR/plink2.err" | tail -8 | sed 's/^/[plink2] /' >&2

# plink2 --export vcf 会把 Chr01 写成数字码 1；用输入 header 的 contig 名还原，
# 使三工具输出可用同一 (CHROM,POS,REF,ALT) 键直接对拍。
RENAME="$WORKDIR/chrom_rename.txt"
"$BCFTOOLS" view -h "$INPUT" 2>/dev/null \
    | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' \
    | awk '{n=$0; gsub(/[^0-9]/,"",n); if(n!=""){n=n+0; print n"\t"$0}}' > "$RENAME"
if [[ -s "$RENAME" ]]; then
    log "[plink2] 还原染色体名 ($(wc -l < "$RENAME") contigs)"
    "$BCFTOOLS" annotate --rename-chrs "$RENAME" "$PREFIX.vcf.gz" -Oz -o "$OUTDIR/plink2.vcf.gz" 2>> "$LOGDIR/plink2.err"
else
    mv "$PREFIX.vcf.gz" "$OUTDIR/plink2.vcf.gz"
fi
"$TABIX" -f -p vcf "$OUTDIR/plink2.vcf.gz"
write_stats plink2 "$OUTDIR/plink2.vcf.gz"
log "[plink2] 完成"
