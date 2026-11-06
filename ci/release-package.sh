#!/usr/bin/env bash
# Package one already-built release binary with its source/build identity.
set -euo pipefail

binary=$1
version=$2
target=$3
compatibility=$4
out=$5
case "$version" in
    v[0-9]*) ;;
    *) echo 'Expected a version beginning with v and a digit' >&2; exit 1 ;;
esac
[[ "$version" =~ ^[a-zA-Z0-9.+-]+$ ]]
case "$target:$compatibility" in
    x86_64-unknown-linux-gnu:glibc217|aarch64-unknown-linux-gnu:glibc217|aarch64-apple-darwin:macos11|x86_64-apple-darwin:macos11) ;;
    *) echo "Unsupported release target/baseline: $target:$compatibility" >&2; exit 1 ;;
esac

name="reed-$version-$target-$compatibility"
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
mkdir -p "$stage/$name" "$out"
cp "$binary" "$stage/$name/reed"
chmod 755 "$stage/$name/reed"
cp LICENSE-MIT LICENSE-APACHE README.md README.zh-CN.md "$stage/$name/"
{
    printf 'Version: %s\nTarget: %s\nCompatibility: %s\n' "$version" "$target" "$compatibility"
    printf 'Source: %s\n' "$(git rev-parse HEAD)"
    rustc -Vv
    "$binary" --version
} > "$stage/$name/BUILD-INFO.txt"
tar -czf "$out/$name.tar.gz" -C "$stage" "$name"
