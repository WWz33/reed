#!/usr/bin/env bash
# ci/check-nameable.sh — the nameability gate (docs/spec/api-contract.md §2,
# doc_v30 P3). Every type in a facade signature must be nameable from `reed`.
#
# `cargo public-api` prints the defining crate's canonical path for types
# that are not defined in `reed`, so any `reed_core::…`/`reed_io::…` path in
# the snapshot is a nameability candidate. The compile probe in
# crates/reed/tests/nameability.rs imports each candidate from `reed`; this
# script keeps the probe complete and compiles it.
#
# Usage: bash ci/check-nameable.sh   (exit 1 on drift)
set -euo pipefail
cd "$(dirname "$0")/.."

python3 - <<'PY'
import re
import sys
from pathlib import Path

snapshot = Path("ci/api/reed.txt").read_text().splitlines()
probe = Path("crates/reed/tests/nameability.rs").read_text()


def terminal(path: str) -> str:
    return path.rstrip(">,) ").split("::")[-1]


names = set()
for line in snapshot:
    # `From<lower::Error> for reed::Error` is error plumbing, not a signature.
    if line.startswith("impl core::convert::From<") or "reed::Error::from(" in line:
        continue
    for match in re.finditer(
        r"\breed_(?:core|io|filter|stats|engine|distance|graph|reduce)::[^ ,;()\[\]{}]+",
        line,
    ):
        name = terminal(match.group(0))
        if re.fullmatch(r"[A-Z][A-Za-z0-9_]*", name):
            names.add(name)

missing = sorted(n for n in names if not re.search(rf"\b{re.escape(n)}\b", probe))
if missing:
    print(
        "nameability: these types appear in the facade API but are not in "
        "crates/reed/tests/nameability.rs: " + ", ".join(missing),
        file=sys.stderr,
    )
    sys.exit(1)
print(f"ok   every facade signature type is in the probe ({len(names)})")
PY

cargo test -p reed --test nameability --no-run
