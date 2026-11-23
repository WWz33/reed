#!/usr/bin/env bash
# ci/check-docs.sh — doc_v30 P4: documentation derived from code, not
# hand-maintained in parallel.
#
#   1. rustdoc has zero warnings (`-D warnings`);
#   2. every library crate's doctests run (each crate has at least one);
#   3. no module is orphaned (`cargo modules orphans`);
#   4. the library guide's entry table matches the facade's `run_*` exports.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== rustdoc"
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --workspace --quiet

echo "== doctests"
cargo test --doc --workspace --quiet

echo "== orphan modules"
# `cargo modules` is installed by the CI gates job (pinned 0.28.0); locally
# the bare error "no such subcommand" reads like a repo bug, so name the
# install command instead (doc_v30 §4.2).
command -v cargo-modules >/dev/null 2>&1 || {
    echo "cargo-modules is not installed; run: cargo install cargo-modules@0.28.0" >&2
    exit 1
}
for crate in reed-core reed-io reed-engine reed-stats reed-filter \
             reed-graph reed-distance reed-reduce reed-oracle; do
    cargo modules orphans --lib -p "$crate"
done

echo "== entry table"
python3 - <<'PY'
import re
import sys
from pathlib import Path

lib = Path("crates/reed/src/lib.rs").read_text()
docs = [
    Path("docsite/src/content/docs/library.md").read_text(),
    Path("docsite/src/content/docs/zh-cn/library.md").read_text(),
]

# `run_*` names the facade actually exports (inside `pub use ...;` items).
exported = set()
for item in re.findall(r"pub use [^;]+;", lib):
    for name in re.findall(r"\brun_[a-z_]+\b", item):
        exported.add(name.removesuffix("_opened"))

for path, text in zip(["library.md", "zh-cn/library.md"], docs):
    documented = set(re.findall(r"`(run_[a-z_]+)`", text))
    if documented != exported:
        print(
            f"{path}: entry table drift; missing={sorted(exported - documented)} "
            f"extra={sorted(documented - exported)}",
            file=sys.stderr,
        )
        sys.exit(1)
print(f"ok   entry table lists all {len(exported)} entry points")
PY
