#!/usr/bin/env python3
"""Per-site/per-sample report parity against plink2 (04-algorithms.md §5).

Known, documented layout differences (values are compared keyed by
chromosome+position):
- plink2 recodes chromosome names to its numeric codes (Chr01 -> 1);
- with a missing VCF ID both tools synthesize `chrom_pos`, so the position is
  recoverable from the ID column.

Usage: compare_reports.py REED_PREFIX PLINK2_PREFIX
"""
import sys


def rows(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return header, [l.rstrip("\n").split("\t") for l in fh]


def key_from_id(id_field):
    """`Chr01_36` / `1:36` -> the trailing integer position."""
    return int(id_field.replace(":", "_").rsplit("_", 1)[-1])


def key(row):
    """(chromosome number, position): plink2 recodes chromosome names, so the
    numeric part identifies the contig on both sides."""
    digits = "".join(c for c in row[0] if c.isdigit())
    return (int(digits) if digits else 0, key_from_id(row[1]))


def cell_ok(a, b, abs_tol, rel_tol):
    """Values agree within `abs_tol + rel_tol*|b|`; NA must match NA."""
    if a == b:
        return True, 0.0, False
    if a == "NA" or b == "NA":
        return False, float("inf"), True
    fa, fb = float(a), float(b)
    diff = abs(fa - fb)
    return diff <= abs_tol + rel_tol * abs(fb), diff, False


def compare(name, reed_path, plink_path, cols, rel_tol=1e-5, abs_tol=1e-6):
    rh, rr = rows(reed_path)
    ph, pr = rows(plink_path)
    want = {}
    for row in pr:
        want[key(row)] = row
    worst_abs, worst_row, checked, bad = 0.0, None, 0, 0
    for row in rr:
        other = want.get(key(row))
        if other is None:
            continue
        for ci, pi in cols:
            ok, diff, na = cell_ok(row[ci], other[pi], abs_tol, rel_tol)
            if not ok:
                bad += 1
                if diff > worst_abs:
                    worst_abs, worst_row = diff, (key(row), rh[ci], row[ci], other[pi], na)
        checked += 1
    ok = bad == 0
    print(f"{'OK' if ok else 'FAIL'} {name}: {checked} rows, {bad} mismatching cells, "
          f"worst |diff| {worst_abs:.3e} (abs {abs_tol:.0e} + rel {rel_tol:.0e})")
    if worst_row:
        print(f"     worst: pos {worst_row[0]} {worst_row[1]}: reed {worst_row[2]!r} "
              f"vs plink2 {worst_row[3]!r} (NA mismatch: {worst_row[4]})")
    return ok


def guarded(fn, *args, **kwargs):
    """Runs a comparison section, skipping it when its files are absent."""
    try:
        return fn(*args, **kwargs)
    except FileNotFoundError:
        return True


def main():
    reed, plink = sys.argv[1], sys.argv[2]
    ok = True
    # plink2 prints ~6 significant digits, so two renderings of the same
    # value can differ by one unit in the last printed digit (1e-5 relative).
    ok &= guarded(compare, "afreq", f"{reed}.afreq", f"{plink}.afreq", [(4, 4), (5, 5)])
    ok &= guarded(compare, "hardy", f"{reed}.hardy", f"{plink}.hardy",
                  [(4, 4), (5, 5), (6, 6), (7, 7), (8, 8), (9, 9)])
    ok &= guarded(compare, "vmiss", f"{reed}.vmiss", f"{plink}.vmiss", [(2, 2), (3, 3), (4, 4)])
    # per-sample files: same order, keyed by IID
    for name, cols in (("het", [(1, 1), (2, 2), (3, 3), (4, 4)]),
                       ("smiss", [(1, 1), (2, 2), (3, 3)])):
        try:
            rh, rr = rows(f"{reed}.{name}")
            ph, pr = rows(f"{plink}.{name}")
        except FileNotFoundError:
            continue
        want = {row[0]: row for row in pr}
        worst, bad = 0.0, 0
        for row in rr:
            other = want.get(row[0])
            if other is None:
                continue
            for ci, pi in cols:
                good, diff, _ = cell_ok(row[ci], other[pi], 1e-6, 1e-5)
                if not good:
                    bad += 1
                    worst = max(worst, diff)
        status = "OK" if bad == 0 else "FAIL"
        ok &= status == "OK"
        print(f"{status} {name}: {len(rr)} samples, {bad} mismatching cells, worst |diff| {worst:.3e}")
    # vcftools depth layouts: idepth (per sample) and ldepth.mean (per site).
    try:
        rh, rr = rows(f"{reed}.idepth")
        ph, pr = rows(f"{plink}.idepth")
        want = {row[0]: row for row in pr}
        worst, bad = 0.0, 0
        for row in rr:
            other = want.get(row[0])
            if other is None:
                continue
            for ci, pi in ((1, 1), (2, 2)):
                good, diff, _ = cell_ok(row[ci], other[pi], 1e-6, 1e-5)
                if not good:
                    bad += 1
                    worst = max(worst, diff)
        status = "OK" if bad == 0 else "FAIL"
        ok &= status == "OK"
        print(f"{status} idepth: {len(rr)} samples, {bad} mismatching cells, worst |diff| {worst:.3e}")
    except FileNotFoundError:
        pass
    try:
        _, rr = rows(f"{reed}.ldepth.mean")
        _, pr = rows(f"{plink}.ldepth.mean")
        want = {(row[0], row[1]): row for row in pr}
        worst, bad, checked = 0.0, 0, 0
        for row in rr:
            other = want.get((row[0], row[1]))
            if other is None:
                continue
            checked += 1
            for ci, pi in ((2, 2), (3, 3)):
                good, diff, _ = cell_ok(row[ci], other[pi], 1e-4, 1e-5)
                if not good:
                    bad += 1
                    worst = max(worst, diff)
        status = "OK" if bad == 0 else "FAIL"
        ok &= status == "OK"
        print(f"{status} ldepth.mean: {checked} sites, {bad} mismatching cells, worst |diff| {worst:.3e}")
    except FileNotFoundError:
        pass
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
