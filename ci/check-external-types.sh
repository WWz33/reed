#!/usr/bin/env bash
# ci/check-external-types.sh — ADR api-contract §3, doc_v30 P3.
#
# Third-party types are opaque in the facade: noodles, serde, rayon and
# memmap2 must not appear in a `reed` public signature. The workspace's own
# model crates are re-exported through `reed` and are listed in the
# allowlist; anything else fails.
#
# The pinned nightly is part of the contract (the tool reads rustdoc JSON).
set -euo pipefail
cd "$(dirname "$0")/.."

cargo +nightly-2026-03-20 check-external-types \
    --manifest-path crates/reed/Cargo.toml \
    --config ci/external-types-allowlist.toml
