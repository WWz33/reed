#!/usr/bin/env bash
# gold/popgen/run_split_parity.sh — the doc_v14 split parity gate
# (docs/rounds/v14/01-review.md §7).
#
# The 0.3.x `reed stats` command is the truth. After the v14 split
# (`stats` keeps the QC reports, `popg` takes the windowed statistics), every
# moved invocation must produce byte-identical output trees. This script
# freezes the 0.3.x outputs once (the "record" mode) and then compares the
# split binary against them (the "check" mode).
#
# Usage:
#   bash gold/popgen/run_split_parity.sh record [VCF POPS REGION]
#   bash gold/popgen/run_split_parity.sh check  [VCF POPS REGION]
#   bash gold/popgen/run_split_parity.sh        # check on the fixture
#
# The frozen tree is keyed by the input's basename, so the fixture and the
# subset (target/gold/work/subset.Chr01-03.vcf.gz, 348 samples, 3 contigs)
# keep separate references. Running it on the subset is what exercises the
# multi-chunk paths (partial merges, the ordered sink) that a 5-record
# fixture cannot: one fixture chunk never reaches `HistAcc::merge` or the
# singletons sink aggregation.
#
# THREADS (default 1) may be raised for the subset runs: outputs are
# thread-invariant by contract, and the parity gate compares bytes.
#
# Environment overrides:
#   REED_REF   0.3.x reference binary (record mode; default the frozen copy
#              under $REED_SCRATCH/v14/)
#   REED       0.4.0 binary under test (check mode; default
#              /tmp/reed-target/release/reed)
#   REED_SCRATCH  scratch root (default $REPO_ROOT/target/gold)
#   FIXTURES   fixture directory (default gold/fixtures)
#
# Gates (each row of CASES below):
#   1. QC family: --freq/--hardy/--het/--missing/--depth/--summary, alone and
#      combined, stay on `stats` with byte-equal outputs.
#   2. Windowed family: pi/dxy/theta/tajima/fst(wc84|hudson)/pbs/sfs/ld,
#      alone and combined, move to `popg` with byte-equal outputs.
#   3. Selection flags (-r/-s/--keep/--geno/--maf/--min-dp) behave the same
#      on both sides.
#   4. Window parameters (--window/--step/--levels/--skip-empty) unchanged.
set -euo pipefail

MODE="${1:-check}"
VCF_ARG="${2:-}"
POPS_ARG="${3:-}"
REGION_ARG="${4:-}"
GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}/v14"
FIXTURES="${FIXTURES:-$GOLD_DIR/fixtures}"
REED_REF="${REED_REF:-$SCRATCH/reed-0.3.0-ref}"
REED="${REED:-/tmp/reed-target/release/reed}"

VCF="${VCF_ARG:-$FIXTURES/fstats.vcf}"
POPS="${POPS_ARG:-$FIXTURES/fstats.pops}"
# The region case needs a contig of the input; the fixture knows `chr1`,
# the subset `Chr01`.
REGION="${REGION_ARG:-chr1}"
THREADS="${THREADS:-1}"
# The three --pbs populations, in first-appearance order.
POPS_LIST="$(cut -f2 "$POPS" | awk '!seen[$0]++' | tr '\n' ' ')"
read -r P1 P2 P3 _ <<< "$POPS_LIST"
[[ -n "${P3:-}" ]] || { echo "$POPS declares fewer than three populations" >&2; exit 1; }

[[ -f "$VCF" ]] || { echo "missing fixture: $VCF" >&2; exit 1; }
[[ -f "$POPS" ]] || { echo "missing fixture: $POPS" >&2; exit 1; }

mkdir -p "$SCRATCH"

# Each case: name | stats flags (0.3.x spelling) | popg flags (0.4.0 spelling)
# The record mode runs the stats spelling against the 0.3.x binary; the check
# mode runs the popg spelling against the binary under test. QC cases keep
# the stats spelling on both sides.
CASES=(
    # --- QC family (stays on `stats`) ---
    "qc_freq|--freq|--freq"
    "qc_hardy|--hardy|--hardy"
    "qc_het|--het|--het"
    "qc_missing|--missing|--missing"
    "qc_depth|--depth|--depth"
    "qc_summary|--summary|--summary"
    "qc_all|--freq --hardy --het --missing --depth --summary|--freq --hardy --het --missing --depth --summary"
    # --- Windowed family (moves to `popg`) ---
    # The fixture is variant-only: --variant-only declares it (the pixy
    # convention), otherwise the invariant-sites check rejects the run.
    "pi|--variant-only --pi|--variant-only --pi"
    "dxy|--variant-only --dxy|--variant-only --dxy"
    "theta|--variant-only --theta|--variant-only --theta"
    "tajima|--variant-only --tajima-d|--variant-only --tajima-d"
    "fst_wc84|--variant-only --fst|--variant-only --fst"
    "fst_hudson|--variant-only --fst hudson|--variant-only --fst hudson"
    "pbs|--variant-only --pbs $P1,$P2,$P3|--variant-only --pbs $P1,$P2,$P3"
    "sfs|--variant-only --sfs --proj 4,4,4|--variant-only --sfs --proj 4,4,4"
    "ld_decay|--variant-only --ld-decay|--variant-only --ld-decay"
    "ld_pairs|--variant-only --ld-decay --ld-pairs|--variant-only --ld-decay --ld-pairs"
    "popg_all|--variant-only --pi --dxy --theta --tajima-d --fst --pbs $P1,$P2,$P3|--variant-only --pi --dxy --theta --tajima-d --fst --pbs $P1,$P2,$P3"
    # --- Selection + window parameters ---
    "window_10k|--variant-only --pi --window 10k --step 5k --levels window,genome|--variant-only --pi --window 10k --step 5k --levels window,genome"
    "skip_empty|--variant-only --pi --skip-empty|--variant-only --pi --skip-empty"
    "sel_geno_maf|--variant-only --pi --geno 0.5 --maf 0.01|--variant-only --pi --geno 0.5 --maf 0.01"
    "sel_region|--variant-only --pi -t $REGION|--variant-only --pi -t $REGION"
)

# The `popg`-spelled cases: which subcommand carries them in 0.4.0.
popg_cases="pi dxy theta tajima fst_wc84 fst_hudson pbs sfs ld_decay ld_pairs popg_all window_10k skip_empty sel_geno_maf sel_region"

is_popg() { [[ " $popg_cases " == *" $1 "* ]]; }

run_case() {
    local bin="$1" sub="$2" flags="$3" outdir="$4"
    local popflag=""
    mkdir -p "$outdir"
    # `--pop` is needed by the windowed family; the QC family ignores it but
    # 0.3.x tolerates it, and the split `stats` keeps it optional, so pass it
    # whenever the case is a popg case.
    if is_popg "$5"; then
        popflag="--pop $POPS"
    fi
    # shellcheck disable=SC2086
    "$bin" "$sub" "$VCF" $popflag $flags -o "$outdir/out" --threads "$THREADS" \
        > "$outdir/stdout.txt" 2> "$outdir/stderr.txt"
    echo "exit=0" > "$outdir/exit.txt"
}

record() {
    [[ -x "$REED_REF" ]] || { echo "reference binary not found: $REED_REF" >&2; exit 1; }
    # The frozen reference is pinned by hash when the digest file is next to
    # it: a rebuilt or replaced 0.3.x binary must not silently define truth.
    local digest="$REED_REF.sha256"
    if [[ -f "$digest" ]]; then
        local want have
        want="$(awk '{print $1}' "$digest")"
        have="$(sha256sum "$REED_REF" | awk '{print $1}')"
        [[ "$want" == "$have" ]] || {
            echo "reference binary $REED_REF does not match $digest" >&2; exit 1; }
    else
        echo "warning: no digest at $digest; the reference is not hash-pinned" >&2
    fi
    local root="$SCRATCH/ref-0.3.0/$(basename "$VCF")"
    rm -rf "$root"
    [[ ! -e "$root" ]] || { echo "cannot clear the stale reference tree $root" >&2; exit 1; }
    mkdir -p "$root"
    for case in "${CASES[@]}"; do
        IFS='|' read -r name stats_flags _ <<< "$case"
        echo "record $name"
        run_case "$REED_REF" "stats" "$stats_flags" "$root/$name" "$name"
    done
    echo "frozen tree: $root"
}

check() {
    [[ -x "$REED" ]] || { echo "binary under test not found: $REED" >&2; exit 1; }
    local ref="$SCRATCH/ref-0.3.0/$(basename "$VCF")"
    [[ -d "$ref" ]] || { echo "run 'record' first (no tree at $ref)" >&2; exit 1; }
    local root="$SCRATCH/new/$(basename "$VCF")"
    rm -rf "$root"
    # A leftover tree would be compared as if freshly produced (false PASS),
    # so a failed clear is fatal rather than ignored.
    [[ ! -e "$root" ]] || { echo "cannot clear the stale output tree $root" >&2; exit 1; }
    mkdir -p "$root"
    local fail=0
    for case in "${CASES[@]}"; do
        IFS='|' read -r name stats_flags popg_flags <<< "$case"
        local sub flags
        if is_popg "$name"; then
            sub="popg"; flags="$popg_flags"
        else
            sub="stats"; flags="$stats_flags"
        fi
        echo "check  $name ($sub)"
        # A missing reference case would otherwise compare nothing and pass:
        # an empty truth tree must fail the gate, not satisfy it.
        if [[ ! -d "$ref/$name" ]]; then
            echo "  FAIL $name: no reference case at $ref/$name"; fail=1; continue
        fi
        if ! run_case "$REED" "$sub" "$flags" "$root/$name" "$name" 2>/dev/null; then
            echo "  FAIL $name: run failed"; fail=1; continue
        fi
        # Compare every produced file except the log (it carries the command
        # line, which legitimately changed: stats -> popg) and stderr.
        local drift=0 compared=0
        while IFS= read -r rel; do
            if [[ "$rel" == *".log" || "$rel" == *"stdout.txt" || "$rel" == *"stderr.txt" || "$rel" == *"exit.txt" ]]; then
                continue
            fi
            compared=$((compared + 1))
            if ! cmp -s "$ref/$name/$rel" "$root/$name/$rel"; then
                echo "  FAIL $name: differs: $rel"; drift=1
            fi
        done < <(cd "$ref/$name" && find . -type f | sed 's|^\./||' | sort)
        # Also assert no file the reference produced is missing.
        while IFS= read -r rel; do
            [[ -e "$root/$name/$rel" ]] || { echo "  FAIL $name: missing: $rel"; drift=1; }
        done < <(cd "$ref/$name" && find . -type f | sed 's|^\./||' | sort)
        # The reference case must carry at least one output to compare (the
        # run's own stdout/stderr/exit files are excluded above).
        if [[ $compared -eq 0 ]]; then
            echo "  FAIL $name: the reference case has no output to compare"; drift=1
        fi
        [[ $drift -eq 0 ]] || fail=1
    done
    if [[ $fail -eq 0 ]]; then
        echo "split parity: PASS (${#CASES[@]} cases, byte-equal)"
    else
        echo "split parity: FAIL" >&2
        exit 1
    fi
}

case "$MODE" in
    record) record ;;
    check)  check ;;
    *) echo "usage: $0 [record|check]" >&2; exit 2 ;;
esac
