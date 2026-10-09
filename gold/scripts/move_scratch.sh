#!/usr/bin/env bash
# gold/scripts/move_scratch.sh — 把仓库内的历史生成物搬到 $REED_SCRATCH。
#
# 一次性迁移（docs/rounds/v5/06-layout.md §4）：`gold/work/`（中间件）与
# `gold/out/`（产物）移到 `$REED_SCRATCH`（默认 `target/gold`），仓库里不再
# 堆 79 G。`mv` 在同一文件系统上是改名，不复制数据。
#
# 用法： bash gold/scripts/move_scratch.sh [--dry-run]
set -euo pipefail

GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}"
DRY="${1:-}"

move_dir() {
    local src="$1" dst="$2"
    [[ -d "$src" ]] || { echo "skip (absent): $src"; return; }
    if [[ "$DRY" == "--dry-run" ]]; then
        echo "would move: $src -> $dst"
        return
    fi
    mkdir -p "$dst"
    # 逐项移动，目标已存在时保留（重复运行安全）。
    shopt -s dotglob nullglob
    local entry base
    for entry in "$src"/*; do
        base="$(basename "$entry")"
        if [[ -e "$dst/$base" ]]; then
            echo "keep existing: $dst/$base (source left in $src)"
            continue
        fi
        mv "$entry" "$dst/$base"
    done
    shopt -u dotglob nullglob
    rmdir "$src" 2>/dev/null || true
    echo "moved: $src -> $dst"
}

move_dir "$GOLD_DIR/work" "$SCRATCH/work"
# gold/out/work 是历史套娃，合并进 scratch/work。
if [[ -d "$GOLD_DIR/out/work" ]]; then
    move_dir "$GOLD_DIR/out/work" "$SCRATCH/work"
fi
move_dir "$GOLD_DIR/out" "$SCRATCH/out"

echo "scratch root: $SCRATCH"
