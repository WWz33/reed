#!/usr/bin/env python3
"""Compares reed's `.fstats` against scikit-allel's admixture.py
(docs/rounds/v5/05-algorithms.md §3.3).

Usage: compare_fstats.py REED_FSTATS VCF POPS [--block-bp N]

The VCF is read with allel's own reader; the populations come from a
two-column `sample<TAB>population` file (reed's `--pop`). Only the tests
reed wrote are checked, keyed by `(stat, populations)`.

Gate: fixture estimates and SEs within 1e-12 relative; full-data estimates
within 1e-9 (05 §3.3). The SE comparison is fixture-only: allel blocks by
*variant count* (`blen`), reed blocks by *bp* (Patterson/admixtools style),
so on real data the two block grids differ and the SEs are not comparable
by construction (the estimates still are).
"""
import sys

import numpy as np
import allel
from allel.stats.admixture import patterson_f2, patterson_f3, patterson_d


def load(path):
    rows = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            rows[(c[0], c[1])] = (c[2], c[3])
    return rows


def main() -> int:
    args = sys.argv[1:]
    rtol = 1e-12
    check_se = True
    files = []
    i = 0
    while i < len(args):
        if args[i] == "--rtol":
            rtol = float(args[i + 1]); i += 2
        elif args[i] == "--no-se":
            check_se = False; i += 1
        else:
            files.append(args[i]); i += 1
    if len(files) != 3:
        raise SystemExit(__doc__)
    reed_path, vcf_path, pops_path = files

    pops = {}
    with open(pops_path) as fh:
        for line in fh:
            if line.strip():
                s, p = line.rstrip("\n").split("\t")[:2]
                pops.setdefault(p, []).append(s)

    callset = allel.read_vcf(vcf_path, fields=["samples", "GT"])
    samples = list(callset["samples"])
    gt = allel.GenotypeArray(callset["calldata/GT"])

    def ac(pop):
        idx = [samples.index(s) for s in pops[pop]]
        sub = gt[:, idx, :]
        return np.column_stack([
            (sub == 0).sum(axis=(1, 2)),
            (sub == 1).sum(axis=(1, 2)),
        ])

    counts = {p: ac(p) for p in pops}
    got = load(reed_path)
    bad = 0
    worst = 0.0

    def compare(key, want_est, want_se=None):
        nonlocal bad, worst
        if key not in got:
            print(f"FAIL {key}: missing in reed", file=sys.stderr); bad += 1; return
        est, se = got[key]
        est = float(est)
        rel = abs(est - want_est) / abs(want_est) if want_est != 0 else abs(est - want_est)
        worst = max(worst, rel)
        if rel > rtol:
            print(f"FAIL {key}: reed={est} allel={want_est} rel={rel:.3e}", file=sys.stderr); bad += 1
        if check_se and want_se is not None and se != "NA":
            rel_se = abs(float(se) - want_se) / abs(want_se) if want_se != 0 else abs(float(se) - want_se)
            if rel_se > rtol:
                print(f"FAIL {key} se: reed={se} allel={want_se} rel={rel_se:.3e}", file=sys.stderr); bad += 1

    for (stat, popkey), (est, se) in got.items():
        names = popkey.split(",")
        if stat == "f2":
            a, b = names
            compare((stat, popkey), float(np.nansum(patterson_f2(counts[a], counts[b]))))
        elif stat == "f3":
            c, a, b = names
            t, _ = patterson_f3(counts[c], counts[a], counts[b])
            compare((stat, popkey), float(np.nansum(t)))
        elif stat in ("f4", "dstat"):
            a, b, c, d = names
            num, den = patterson_d(counts[a], counts[b], counts[c], counts[d])
            want = float(np.nansum(num)) if stat == "f4" else float(np.nansum(num) / np.nansum(den))
            compare((stat, popkey), want)
        else:
            print(f"note: {stat} has no allel counterpart, skipped")

    status = "PASS" if bad == 0 else "FAIL"
    print(f"{status} f-stats: {len(got)} test(s), max rel {worst:.3e} (rtol {rtol:g})")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
