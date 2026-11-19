#!/usr/bin/env bash
# ci/check-package.sh — validate the crates.io package graph for `reed`.
#
# `reed` has normal dependencies on the workspace library crates. Cargo turns
# those path dependencies into registry dependencies in a published tarball,
# so every one of them must be publishable, versioned, and released first.
# This gate validates that local contract. It does not pretend that a clean
# `cargo package --no-verify` can pass before those crates exist on crates.io.
set -euo pipefail
cd "$(dirname "$0")/.."

list=$(mktemp)
trap 'rm -f "$list"' EXIT

cargo package -p reed --list --allow-dirty >"$list"
grep -qx 'src/lib.rs' "$list" || {
    echo "package gate: src/lib.rs is not in the published file list" >&2
    exit 1
}

python3 - <<'PY2'
import json
import subprocess
import sys

meta = json.loads(
    subprocess.check_output(["cargo", "metadata", "--format-version", "1", "--no-deps"])
)
packages = {package["name"]: package for package in meta["packages"]}
reed = packages.get("reed")
if reed is None:
    sys.exit("package gate: reed is not in the workspace")

if reed.get("publish") == []:
    sys.exit("package gate: reed must be a publish target")

publish_order = [
    "reed-core",
    "reed-io",
    "reed-engine",
    "reed-stats",
    "reed-filter",
    "reed-graph",
    "reed-distance",
    "reed-reduce",
    "reed",
]

missing = [name for name in publish_order if name not in packages]
if missing:
    sys.exit("package gate: missing workspace package(s): " + ", ".join(missing))

facade_version = reed["version"]
for name in publish_order:
    package = packages[name]
    if package.get("publish") == []:
        sys.exit(f"package gate: {name} is publish = false but is in reed's release graph")
    if package["version"] != facade_version:
        sys.exit(
            f"package gate: {name} version {package['version']} does not match reed "
            f"version {facade_version}"
        )

    for dep in package["dependencies"]:
        if dep.get("kind") is not None or not dep.get("path"):
            continue
        dependency = packages.get(dep["name"])
        if dependency is None:
            sys.exit(
                f"package gate: {name} has a path dependency outside the workspace: "
                f"{dep['name']}"
            )
        if dep.get("req") in (None, "*"):
            sys.exit(
                f"package gate: {name} -> {dep['name']} has no registry version requirement"
            )
        if dependency.get("publish") == []:
            sys.exit(
                f"package gate: {name} -> {dep['name']} depends on publish = false crate"
            )

positions = {name: index for index, name in enumerate(publish_order)}
for name in publish_order:
    for dep in packages[name]["dependencies"]:
        if dep.get("kind") is None and dep.get("path") and dep["name"] in positions:
            if positions[dep["name"]] >= positions[name]:
                sys.exit(
                    f"package gate: release order is invalid for {name} -> {dep['name']}"
                )

print("ok   reed package file selection is valid")
print("ok   all normal workspace dependencies are publishable and versioned")
print("release order: " + " -> ".join(publish_order))
print(
    "note: cargo package -p reed --no-verify requires the preceding crates "
    "to be published to the selected registry first"
)
PY2
