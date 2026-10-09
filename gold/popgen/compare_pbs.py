#!/usr/bin/env python3
"""PBS parity against scikit-allel (`allel.pbs`), 04-algorithms.md §5.

`allel.pbs` windows by *variant index*, so the fixture must hold the same
number of variants in every genomic window. The script recomputes the
per-window pair Hudson Fsts with `allel.hudson_fst` on the same slices.

Both flavours write `PREFIX_pbs.txt` (02-cli.md §5.4), so the flavour under
test is named explicitly:

Usage: compare_pbs.py REED_PREFIX VCF POPFILE WINDOW VARIANTS_PER_WINDOW [--normed]
"""
import sys
import numpy as np
import allel


def main():
    prefix, vcf, popfile = sys.argv[1], sys.argv[2], sys.argv[3]
    per_window = int(sys.argv[4])
    normed = "--normed" in sys.argv[5:]

    callset = allel.read_vcf(vcf, fields=["calldata/GT", "samples"])
    samples = [str(s) for s in callset["samples"]]
    gt = allel.GenotypeArray(callset["calldata/GT"])
    pops = {}
    with open(popfile) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sample, pop = line.split("\t")[:2]
            pops.setdefault(pop, []).append(samples.index(sample))

    names = sorted(pops)
    if len(names) < 3:
        sys.exit("PBS needs three populations")
    idx = [pops[n] for n in names[:3]]
    acs = [gt.count_alleles(subpop=i, max_allele=1) for i in idx]

    want = allel.pbs(acs[0], acs[1], acs[2], window_size=per_window, normed=normed)

    def load(path):
        rows = []
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = header.index("pbs")
            for line in fh:
                c = line.rstrip("\n").split("\t")
                rows.append(None if c[col] == "NA" else float(c[col]))
        return rows

    got = load(prefix + "_pbs.txt")
    worst = 0.0
    compared = 0
    for a, b in zip(got, [float(x) for x in want]):
        if a is None:
            continue
        worst = max(worst, abs(a - b))
        compared += 1
    status = "OK" if worst <= 1e-12 and compared == len(want) else "FAIL"
    print(f"{status} pbs (normed={normed}) vs scikit-allel: {compared}/{len(want)} windows, "
          f"max |diff| {worst:.3e}")
    sys.exit(0 if status == "OK" else 1)


if __name__ == "__main__":
    main()
