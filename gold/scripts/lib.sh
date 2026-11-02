#!/usr/bin/env bash
# gold/scripts/lib.sh — 公共函数库。被 run_*.sh source。
# 负责：解析路径、加载配置、检查工具、日志、样本保留名单、结果 JSON。
set -euo pipefail

# ---- 路径解析 --------------------------------------------------------------
GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${CONFIG_FILE:-$GOLD_DIR/config/filters.conf}"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"          # reed-rs/
# 生成物一律进 target/（或 $REED_SCRATCH）；仓库里不堆大文件（docs/rounds/v5/06-layout.md）。
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}"
OUTDIR="${OUTDIR:-$SCRATCH/out}"
WORKDIR="${WORKDIR:-$SCRATCH/work}"
LOGDIR="${LOGDIR:-$OUTDIR/logs}"

mkdir -p "$OUTDIR" "$WORKDIR" "$LOGDIR"

# ---- 加载配置 --------------------------------------------------------------
# shellcheck source=/dev/null
source "$CONFIG_FILE"

# 输入默认值：仓库上级的 subvcf（小集，3000 记录 / 348 样本）。
# 全量 clean.* 只在显式 INPUT=... 时使用；测试脚本默认都跑小集。
if [[ -z "${INPUT:-}" ]]; then
    for cand in "$REPO_ROOT/../subvcf.vcf.gz" "$REPO_ROOT/subvcf.vcf.gz" \
                "$REPO_ROOT/../subvcf.bcf" "$REPO_ROOT/subvcf.bcf"; do
        [[ -e "$cand" ]] && INPUT="$cand" && break
    done
fi
[[ -n "${INPUT:-}" ]] || die "INPUT 未设置且未找到默认输入 (subvcf.vcf.gz / subvcf.bcf)"

# PATH：把 biosoft 及 htslib 放前面
export PATH="$BIOSOFT:$BIOSOFT/htslib-1.23:$PATH"
export BCFTOOLS_PLUGINS

# ---- 日志 ------------------------------------------------------------------
_ts() { date '+%Y-%m-%d %H:%M:%S'; }
log()  { printf '[%s] %s\n' "$(_ts)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(_ts)" "$*" >&2; exit 1; }

# ---- 工具检查 --------------------------------------------------------------
require_tools() {
    local t
    for t in "$@"; do
        local bin="${!t:-}"
        [[ -n "$bin" && -x "$bin" ]] || die "工具不可执行: $t=$bin"
    done
}

tool_versions() {
    {
        echo "bcftools: $("$BCFTOOLS" --version 2>&1 | head -1)"
        echo "vcftools: $("$VCFTOOLS" --version 2>&1 | head -1)"
        echo "plink2:   $("$PLINK2" --version 2>&1 | head -1)"
    }
}

# ---- 样本保留名单 (全量输入上统一计算一次，三工具共享) ----------------------
# 用 bcftools stats -s - 的 PSC 行：$3=sample $4=nRefHom $5=nNonRefHom
# $6=nHets $14=nMissing。缺失率 = nMissing / (called + missing)。
build_sample_keep() {
    local out="$OUTDIR/samples.keep.txt"
    if [[ -s "$out" ]]; then log "复用已有样本名单: $out"; return; fi
    log "计算样本缺失率 (阈值 $SAMPLE_MISS_MAX) ..."
    local stats="$WORKDIR/sample_stats.psc"
    # shellcheck disable=SC2086
    "$BCFTOOLS" stats ${REGION:+-r "$REGION"} -s - "$INPUT" 2>/dev/null \
        | awk -v FS='\t' -v max="$SAMPLE_MISS_MAX" \
              '$1=="PSC"{called=$4+$5+$6; tot=called+$14;
                        if(tot>0 && $14/tot<=max) print $3 > "'"$out"'"}'
    local n; n=$(wc -l < "$out")
    log "样本保留: $n 个 (阈值缺失率<=$SAMPLE_MISS_MAX)"
}

# ---- 结果 JSON -------------------------------------------------------------
# write_stats <tool> <vcf> <extra_json_fragment>
write_stats() {
    local tool="$1" vcf="$2" extra="${3:-}"
    local nsites=0 nsamples=0
    if [[ -s "$vcf" ]]; then
        nsites=$("$BCFTOOLS" view -H "$vcf" 2>/dev/null | wc -l)
        nsamples=$("$BCFTOOLS" query -l "$vcf" 2>/dev/null | wc -l)
    fi
    cat > "$OUTDIR/$tool.stats.json" <<EOF
{
  "tool": "$tool",
  "vcf": "$(basename "$vcf")",
  "sites": $nsites,
  "samples": $nsamples,
  "input": "$INPUT",
  "region": "${REGION:-}",
  "spec": {
    "sample_miss_max": $SAMPLE_MISS_MAX,
    "min_qual": $MIN_QUAL,
    "require_pass": $REQUIRE_PASS,
    "gt_dp_min": $GT_DP_MIN,
    "gt_dp_max": $GT_DP_MAX,
    "gt_gq_min": $GT_GQ_MIN,
    "site_miss_max": $SITE_MISS_MAX,
    "maf_min": $MAF_MIN,
    "hwe_min": $HWE_MIN
  }${extra:+,$extra}
}
EOF
    log "[$tool] sites=$nsites samples=$nsamples -> $OUTDIR/$tool.stats.json"
}
