#!/usr/bin/env python3
"""Compares two projected-SFS files (reed vs easySFS), numerically.

Both layouts carry the spectrum on a single line: dadi writes
`<dim> <folded>` / values / mask, fastsimcoal2 writes `1 observation` /
`d0_…` / values. The comparison reads the values line, so either layout is
accepted, and checks `|reed - ref| <= rtol * |ref|` (default rtol 1e-6, the
05 §3.1 gate: before rounding — easySFS's own float accumulation differs
from reed's compensated sums at ~1e-15).

Known, documented difference: the integer form (`--sfs-int` vs
`--dtype int`) can differ by 1 in a bin whose exact value is a half-integer,
because easySFS's float noise can land just below the half (e.g. its
1.499999999999999 rounds to 1 while reed's 1.5 rounds to 2).

Usage: compare_sfs.py REED_FILE EASYSFS_FILE [--rtol R] [--label NAME]
Exit 1 on a mismatch beyond the tolerance.
"""
import sys


def values(path):
    """The spectrum line: the first line whose fields all parse as floats."""
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if not fields:
                continue
            try:
                return [float(x) for x in fields]
            except ValueError:
                continue
    raise SystemExit(f"compare_sfs: no value line in {path}")


def main() -> int:
    args = sys.argv[1:]
    rtol = 1e-6
    label = ""
    files = []
    i = 0
    while i < len(args):
        if args[i] == "--rtol":
            rtol = float(args[i + 1])
            i += 2
        elif args[i] == "--label":
            label = args[i + 1]
            i += 2
        else:
            files.append(args[i])
            i += 1
    if len(files) != 2:
        raise SystemExit(__doc__)
    reed, ref = files
    a, b = values(reed), values(ref)
    if len(a) != len(b):
        print(f"FAIL {label}: {len(a)} bins vs {len(b)}", file=sys.stderr)
        return 1
    worst = 0.0
    bad = 0
    for i, (x, y) in enumerate(zip(a, b)):
        diff = abs(x - y)
        rel = diff / abs(y) if y != 0 else diff
        worst = max(worst, rel)
        if rel > rtol:
            bad += 1
            print(f"FAIL {label}: bin {i}: reed={x} easySFS={y} rel={rel:.3e}", file=sys.stderr)
    status = "PASS" if bad == 0 else "FAIL"
    print(f"{status} {label}: {len(a)} bins, max rel {worst:.3e} (rtol {rtol:g})")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
