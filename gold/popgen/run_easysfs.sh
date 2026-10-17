#!/usr/bin/env bash
# gold/popgen/run_easysfs.sh — projected-SFS parity against easySFS
# (docs/rounds/v5/05-algorithms.md §3.1).
#
# Usage:
#   bash gold/popgen/run_easysfs.sh [VCF] [POPS] [PROJ] [OUTDIR]
#
# Defaults to the shared fixture (4 diploid samples, one population,
# projection 6). For the full-data run pass clean.vcf.gz, clean.3pop.txt and
# `PROJ = 2·min(n_i)` (the smallest group's chromosome count) so every group
# projects to the same size.
#
# Gate: every bin within rtol 1e-6 of easySFS's, before rounding (the int
# form can differ by 1 at a half-integer; see compare_sfs.py). Both SFS
# entries run — `reed popg --sfs` and `reed convert --sfs` — and their
# trees must be byte-identical (doc_v13 §4.1, D13-6).
set -euo pipefail

GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}"
EASYSFS="${EASYSFS:-$HOME/reed/easySFS/easySFS.py}"
REED="${REED:-/tmp/reed-target/release/reed}"

VCF="${1:-$GOLD_DIR/fixtures/sfs.vcf}"
POPS="${2:-$GOLD_DIR/fixtures/sfs.pops}"
PROJ="${3:-6}"
OUT="${4:-$SCRATCH/out/easysfs}"

[[ -x "$REED" ]] || { echo "reed binary not executable: $REED" >&2; exit 1; }
[[ -f "$EASYSFS" ]] || { echo "easySFS.py not found: $EASYSFS" >&2; exit 1; }

# easySFS wants one projection per population; reed broadcasts a single
# value. Repeat the value NPOPS times.
NPOPS=$(awk -F'\t' '!/^#/ && NF>0 {print $2}' "$POPS" | sort -u | wc -l)
EASY_PROJ=""
for ((i = 0; i < NPOPS; i++)); do EASY_PROJ+="${EASY_PROJ:+,}$PROJ"; done

fail=0
for fold in unfolded folded; do
    mode=(); [[ "$fold" == "folded" ]] && mode=(folded)
    easy="$OUT/easysfs_$fold"
    # The two entries write sibling directories whose prefixes share the
    # basename `out`: the joint fastsimcoal2 file names embed the prefix, so
    # only equal basenames make the byte-level `diff -r` meaningful.
    stat_dir="$OUT/reed_stat_$fold"
    conv_dir="$OUT/reed_conv_$fold"
    reed="$stat_dir/out"
    mkdir -p "$easy" "$stat_dir" "$conv_dir"
    echo "== easySFS ($fold) =="
    easy_args=(-i "$VCF" -p "$POPS" --proj "$EASY_PROJ" -a -f -o "$easy")
    [[ "$fold" == "unfolded" ]] && easy_args+=(--unfolded)
    python3 "$EASYSFS" "${easy_args[@]}" >/dev/null
    echo "== reed ($fold) =="
    "$REED" popg "$VCF" --pop "$POPS" --sfs "${mode[@]}" --proj "$PROJ" -o "$reed" --verbosity 0
    "$REED" convert "$VCF" --pop "$POPS" --sfs "${mode[@]}" --proj "$PROJ" -o "$conv_dir/out" --verbosity 0
    if diff -r "$stat_dir/out_sfs" "$conv_dir/out_sfs" >/dev/null 2>&1; then
        echo "  convert == popg ($fold)"
    else
        echo "  convert != popg ($fold)" >&2
        fail=1
    fi

    # easySFS also writes each population's unprojected `<pop>.sfs`; the
    # projected one is the comparison target.
    for ref in "$easy"/dadi/*-"$PROJ".sfs; do
        base="$(basename "$ref")"
        got="$reed""_sfs/dadi/$base"
        if [[ ! -f "$got" ]]; then
            echo "missing reed output: $got" >&2
            fail=1
            continue
        fi
        python3 "$GOLD_DIR/popgen/compare_sfs.py" "$got" "$ref" --label "$fold/$base" || fail=1
    done
    # The pairwise joint spectra: easySFS writes `<A>-<B>.sfs` (no projection
    # suffix) per pair; the values line is the comparison target (the header
    # carries population labels reed omits, like the 1D writer).
    for ref in "$easy"/dadi/*-*.sfs; do
        base="$(basename "$ref")"
        case "$base" in *-"$PROJ".sfs) continue ;; esac
        got="$reed""_sfs/dadi/$base"
        if [[ ! -f "$got" ]]; then
            echo "missing reed joint output: $got" >&2
            fail=1
            continue
        fi
        python3 "$GOLD_DIR/popgen/compare_sfs.py" "$got" "$ref" --label "$fold/joint/$base" || fail=1
    done
done

[[ "$fail" -eq 0 ]] || { echo "SFS parity: FAIL" >&2; exit 1; }
echo "SFS parity: PASS (rtol 1e-6, before rounding)"
