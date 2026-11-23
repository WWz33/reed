#!/bin/bash
# Publish the public subset of the exact commit at origin/dev to origin/main.
#
# Internal paths (docs/, docsite/, gold/, bench/, and crate tests) are excluded.
# Usage: ci/publish.sh ["commit message"]

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

branch="$(git symbolic-ref --quiet --short HEAD || true)"
if [[ "$branch" != dev ]]; then
    echo "publish.sh: run from the dev branch" >&2
    exit 1
fi

status="$(git status --short --untracked-files=all)"
if [[ -n "$status" ]]; then
    printf 'publish.sh: worktree must be clean before publishing:\n%s\n' "$status" >&2
    exit 1
fi

source_commit="$(git rev-parse --verify 'HEAD^{commit}')"
remote_dev="$(git ls-remote origin refs/heads/dev | awk '$2 == "refs/heads/dev" { print $1; exit }')"
if [[ -z "$remote_dev" || "$source_commit" != "$remote_dev" ]]; then
    echo "publish.sh: push and sync this exact commit to origin/dev before publishing" >&2
    exit 1
fi

PUB_PATHS=(
    .gitignore .github Cargo.toml Cargo.lock LICENSE-APACHE LICENSE-MIT
    README.md README.zh-CN.md clippy.toml deny.toml rustfmt.toml crates ci
)

# Safety net: never allow an internal root into the public path list.
for p in docs docsite gold bench; do
    for public_path in "${PUB_PATHS[@]}"; do
        if [[ "$public_path" == "$p" ]]; then
            echo "publish.sh: internal path '$p' found in PUB_PATHS" >&2
            exit 1
        fi
    done
done

# Build the index from committed objects, never from files in the worktree.
tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/reed-pub.XXXXXX")"
export GIT_INDEX_FILE="$tmp_dir/index"
trap 'rm -rf -- "$tmp_dir"' EXIT
git read-tree --empty
git reset --quiet "$source_commit" -- "${PUB_PATHS[@]}"
while IFS= read -r -d '' path; do
    git update-index --force-remove -- "$path"
done < <(git ls-files -z -- ':(glob)crates/**/tests/**' .github/workflows/docs.yml)
tree="$(git write-tree)"

if git ls-tree -r --name-only "$tree" | grep -E '^(DESIGN\.md|docs/|docsite/|gold/|bench/|tests/)|/tests/|^\.github/workflows/docs\.yml$' >/dev/null; then
    echo "publish.sh: internal path leaked into publish tree" >&2
    exit 1
fi

parent="$(git ls-remote origin refs/heads/main | awk '$2 == "refs/heads/main" { print $1; exit }')"
message="${1:-publish: sync public subset}"
if [[ -n "$parent" ]]; then
    commit="$(git commit-tree "$tree" -p "$parent" -m "$message" -m "Source-Commit: $source_commit")"
else
    commit="$(git commit-tree "$tree" -m "${1:-reed v0.1: initial public import}" -m "Source-Commit: $source_commit")"
fi

# The new snapshot is a child of the observed main tip. A concurrent update
# makes this a non-fast-forward push, which Git rejects without overwriting it.
git push origin "$commit:refs/heads/main"
echo "published: $commit -> origin/main (source $source_commit)"
