#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d "${TMPDIR:-/tmp}/reed-publish-test.XXXXXX")"
trap 'rm -rf -- "$tmp"' EXIT
real_git="$(command -v git)"
remote="$tmp/remote.git"
source_repo="$tmp/source"
git init --bare --initial-branch=main "$remote" >/dev/null
git init --initial-branch=dev "$source_repo" >/dev/null
git -C "$source_repo" config user.name "Publish test"
git -C "$source_repo" config user.email "publish-test@example.invalid"
git -C "$source_repo" config commit.gpgsign false
git -C "$source_repo" remote add origin "$remote"

write_file() {
    mkdir -p "$(dirname "$1")"
    printf '%s\n' "$2" > "$1"
}

write_file "$source_repo/.gitignore" "target/"
write_file "$source_repo/.github/workflows/ci.yml" "public workflow"
write_file "$source_repo/.github/workflows/docs.yml" "private workflow"
write_file "$source_repo/Cargo.toml" "public manifest"
write_file "$source_repo/Cargo.lock" "public lockfile"
write_file "$source_repo/LICENSE-APACHE" "license"
write_file "$source_repo/LICENSE-MIT" "license"
write_file "$source_repo/README.md" "committed README"
write_file "$source_repo/README.zh-CN.md" "committed README zh"
write_file "$source_repo/clippy.toml" "clippy"
write_file "$source_repo/deny.toml" "deny"
write_file "$source_repo/rustfmt.toml" "rustfmt"
write_file "$source_repo/crates/reed/src/lib.rs" "pub fn public() {}"
write_file "$source_repo/crates/reed/tests/private.rs" "private test"
mkdir -p "$source_repo/ci"
cp "$repo_root/ci/publish.sh" "$source_repo/ci/publish.sh"
write_file "$source_repo/docs/rounds/v30/06-final-audit.md" "internal docs"
write_file "$source_repo/docsite/index.md" "internal docs"
write_file "$source_repo/gold/fixture.txt" "internal fixture"
write_file "$source_repo/bench/bench.rs" "internal benchmark"

git -C "$source_repo" add -A
git -C "$source_repo" commit -qm "fixture source"
source_commit="$(git -C "$source_repo" rev-parse HEAD)"
git -C "$source_repo" push -q origin HEAD:refs/heads/dev
git -C "$source_repo" push -q origin HEAD:refs/heads/main

before_main="$(git --git-dir="$remote" rev-parse refs/heads/main)"
printf 'uncommitted change\n' >> "$source_repo/README.md"
if (cd "$source_repo" && bash ci/publish.sh) >"$tmp/dirty.log" 2>&1; then
    echo "publish test: dirty worktree was accepted" >&2
    exit 1
fi
grep -q 'worktree must be clean' "$tmp/dirty.log"
test "$(git --git-dir="$remote" rev-parse refs/heads/main)" = "$before_main"
git -C "$source_repo" restore -- README.md

(cd "$source_repo" && bash ci/publish.sh "test publish")
published="$(git --git-dir="$remote" rev-parse refs/heads/main)"
test "$(git --git-dir="$remote" show "$published:README.md")" = "committed README"
git --git-dir="$remote" show -s --format=%B "$published" | grep -q "Source-Commit: $source_commit"

assert_absent() {
    if "$real_git" --git-dir="$remote" cat-file -e "$published:$1" 2>/dev/null; then
        printf 'publish test: internal path leaked: %s\n' "$1" >&2
        exit 1
    fi
}
assert_absent "docs/rounds/v30/06-final-audit.md"
assert_absent "docsite/index.md"
assert_absent "gold/fixture.txt"
assert_absent "bench/bench.rs"
assert_absent "crates/reed/tests/private.rs"
assert_absent ".github/workflows/docs.yml"

racer="$tmp/racer"
git clone -q "$remote" "$racer"
git -C "$racer" config user.name "Concurrent publisher"
git -C "$racer" config user.email "concurrent@example.invalid"
write_file "$racer/concurrent.txt" "concurrent update"
git -C "$racer" add concurrent.txt
git -C "$racer" commit -qm "concurrent main update"
race_commit="$(git -C "$racer" rev-parse HEAD)"
shim_dir="$tmp/git-shim"
mkdir -p "$shim_dir"
cat > "$shim_dir/git" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == push && -n "${PUBLISH_TEST_RACE_REPO:-}" ]]; then
    race_repo="$PUBLISH_TEST_RACE_REPO"
    unset PUBLISH_TEST_RACE_REPO
    GIT_INDEX_FILE="$race_repo/.git/index" \
        "$PUBLISH_TEST_REAL_GIT" -C "$race_repo" push origin HEAD:refs/heads/main >/dev/null
fi
exec "$PUBLISH_TEST_REAL_GIT" "$@"
SH
chmod +x "$shim_dir/git"

if (cd "$source_repo" && PATH="$shim_dir:$PATH" \
    PUBLISH_TEST_REAL_GIT="$real_git" PUBLISH_TEST_RACE_REPO="$racer" \
    bash ci/publish.sh "race test") >"$tmp/race.log" 2>&1; then
    echo "publish test: concurrent main update was overwritten" >&2
    exit 1
fi
test "$(git --git-dir="$remote" rev-parse refs/heads/main)" = "$race_commit"
printf 'publish safety: PASS\n'
