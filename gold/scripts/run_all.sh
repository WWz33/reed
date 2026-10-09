#!/usr/bin/env bash
# gold/scripts/run_all.sh — 并行跑三个工具，然后做三方一致性对拍。
#   1. 在全量输入上统一计算样本保留名单 (三工具共享，保证样本集一致)；
#   2. bcftools / vcftools / plink2 三个分支并行执行；
#   3. 汇总 -> out/gold.vcf.gz (= bcftools.full) + out/concordance.{txt,tsv,md}
#
# 用法:
#   bash gold/scripts/run_all.sh
#   REGION=Chr01:1-2000000 OUTDIR=/tmp/goldtest bash gold/scripts/run_all.sh   # 烟测
#   HWE_MIN=0 bash gold/scripts/run_all.sh                                     # 关闭 HWE 层
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
require_tools BCFTOOLS VCFTOOLS PLINK2 BGZIP TABIX

t0=$(date +%s)
log "==================== gold filter run ===================="
log "input : $INPUT"
log "region: ${REGION:-<全量>}"
log "config: $CONFIG_FILE"
tool_versions | sed 's/^/        /' >&2
log "spec  : qual>=$MIN_QUAL pass=$REQUIRE_PASS dp=[$GT_DP_MIN,$GT_DP_MAX] gq>=$GT_GQ_MIN site_miss<=$SITE_MISS_MAX maf>=$MAF_MIN hwe>=$HWE_MIN"
log "        INFO(bcftools only): QD>=$INFO_QD_MIN FS<=$INFO_FS_MAX MQ>=$INFO_MQ_MIN SOR<=$INFO_SOR_MAX"

# 1) 共享样本名单 (全量)
build_sample_keep

# 2) 三分支并行
declare -a names=(bcftools vcftools plink2)
declare -a pids=()
for n in "${names[@]}"; do
    ( bash "$GOLD_DIR/scripts/run_$n.sh" > "$LOGDIR/run_$n.out" 2>&1; echo $? > "$WORKDIR/run_$n.rc" ) &
    pids+=($!)
    log "launched run_$n.sh (pid $!)"
done
for p in "${pids[@]}"; do wait "$p" || true; done

# 3) 检查
fail=0
for n in "${names[@]}"; do
    rc=$(cat "$WORKDIR/run_$n.rc" 2>/dev/null || echo 1)
    if [[ "$rc" != "0" ]]; then
        log "FAIL run_$n.sh rc=$rc"; tail -25 "$LOGDIR/run_$n.out" >&2; fail=1
    else
        log "OK   run_$n.sh"
    fi
done
[[ "$fail" == "0" ]] || die "至少一个工具分支失败，见 $LOGDIR/"

# 4) gold = bcftools.full
cp -f "$OUTDIR/bcftools.full.vcf.gz" "$OUTDIR/gold.vcf.gz"
cp -f "$OUTDIR/bcftools.full.vcf.gz.tbi" "$OUTDIR/gold.vcf.gz.tbi" 2>/dev/null || true
log "gold  -> $OUTDIR/gold.vcf.gz ($("$BCFTOOLS" view -H "$OUTDIR/gold.vcf.gz" | wc -l) sites)"

# 5) 三方对拍
python3 "$GOLD_DIR/scripts/compare_tools.py" --outdir "$OUTDIR" 2>&1 | tee "$OUTDIR/concordance.txt"

t1=$(date +%s)
log "==================== done in $((t1-t0))s ===================="
