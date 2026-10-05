#!/bin/bash
# Publish the public subset of reed to origin/main.
#
# Internal paths (DESIGN.md, docs/, bench/) are commit-only — never pushed.
# This script builds a scrubbed tree via a temporary index (no branch
# switching, working tree untouched), then force-updates origin/main.
#
# Usage: scripts/publish.sh ["commit message"]

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

PUB_PATHS=".gitignore Cargo.toml Cargo.lock README.md README.zh-CN.md crates scripts"

# Safety net: never let an internal path into the published tree.
for p in DESIGN.md docs bench; do
    if [[ "$PUB_PATHS" == *"$p"* ]]; then
        echo "publish.sh: internal path '$p' found in PUB_PATHS" >&2
        exit 1
    fi
done

# Stage the public subset into a temporary index (leave the real index alone).
idx="$(mktemp -u /tmp/reed-pub-idx.XXXXXX)"
export GIT_INDEX_FILE="$idx"
trap 'rm -f "$GIT_INDEX_FILE"' EXIT
git read-tree --empty
git add $PUB_PATHS
tree="$(git write-tree)"

if git ls-tree -r --name-only "$tree" | grep -qE '^(DESIGN\.md|docs/|bench/)'; then
    echo "publish.sh: internal path leaked into publish tree" >&2
    exit 1
fi

parent="$(git ls-remote origin refs/heads/main | cut -f1)"
if [[ -n "$parent" ]]; then
    commit="$(git commit-tree "$tree" -p "$parent" -m "${1:-publish: sync public subset}")"
else
    commit="$(git commit-tree "$tree" -m "${1:-reed v0.1: initial public import}")"
fi

git push --force origin "$commit:refs/heads/main"
echo "published: $commit -> origin/main"
