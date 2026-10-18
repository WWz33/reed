#!/usr/bin/env python3
"""Numeric row comparison for the parity gates.

`cmp` is the rule for byte-equality gates, but two tools printing the same
value at 6 significant digits can differ by one unit in the last digit —
plink2's per-ALT frequencies summed to recover reed's collapsed frequency
are the case this exists for (0.333333 + 0.333333 = 0.666666 against
0.666667).

Usage: compare_numeric.py LEFT RIGHT [RTOL]

Rows must be tab-separated and aligned (same order); numeric cells compare
with `abs(a-b) <= 1e-6 + rtol*|b|`, other cells must match exactly. Exit 0
when every row agrees, 1 otherwise (the offending rows go to stderr).
"""
import sys


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    left_path, right_path = sys.argv[1], sys.argv[2]
    rtol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-5
    with open(left_path) as fh:
        left = [line.rstrip("\n").split("\t") for line in fh]
    with open(right_path) as fh:
        right = [line.rstrip("\n").split("\t") for line in fh]
    if len(left) != len(right):
        print(f"row counts differ: {len(left)} vs {len(right)}", file=sys.stderr)
        return 1
    bad = 0
    for n, (a, b) in enumerate(zip(left, right), 1):
        if len(a) != len(b):
            print(f"row {n}: {len(a)} vs {len(b)} columns", file=sys.stderr)
            bad += 1
            continue
        for ca, cb in zip(a, b):
            if ca == cb:
                continue
            try:
                fa, fb = float(ca), float(cb)
            except ValueError:
                print(f"row {n}: {ca!r} != {cb!r}", file=sys.stderr)
                bad += 1
                break
            if abs(fa - fb) > 1e-6 + rtol * abs(fb):
                print(f"row {n}: {ca} != {cb}", file=sys.stderr)
                bad += 1
                break
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
