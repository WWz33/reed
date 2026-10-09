#!/usr/bin/env python3
"""Compares reed pixy-layout output against pixy / vcftools references.

Window identity is `(population(s), chromosome, window_pos_1)`: reed clips the
last window of a contig to the contig length (04-algorithms.md §4), while
pixy prints the full window size there. The values are unaffected; only the
printed `window_pos_2` of that one window per contig differs, so the end
coordinate is not part of the key.

Usage: compare.py REED_PREFIX REFDIR [--vcftools]

- pixy reference:   REFDIR/pixy_bypass/pixy_pi.txt, pixy_dxy.txt
- vcftools:         REFDIR/vcftools_pi_<pop>.windowed.pi (π only)
Exit 1 on any mismatch beyond the stated tolerance.
"""
import sys, os, glob

def load_tsv(path, key_cols, val_col):
    """Value column as float, or None for `NA` (empty windows)."""
    rows = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            v = c[val_col]
            rows[tuple(c[i] for i in key_cols)] = None if v == "NA" else float(v)
    return rows

def compare(name, expected, actual, tol, abs_tol=1e-12):
    """Values agree when |a - v| <= abs_tol + tol*|v|.

    The absolute floor keeps near-zero statistics (Fst can be ~1e-8) from
    failing on f64 round-off while still catching real algorithmic errors.
    """
    worst, worst_abs, n, missing, na_mismatch = 0.0, 0.0, 0, 0, 0
    for k, v in expected.items():
        if k not in actual:
            missing += 1
            continue
        a = actual[k]
        if v is None or a is None:
            if v is not a:
                na_mismatch += 1
            continue
        denom = max(abs(v), 1e-300)
        diff = abs(a - v)
        worst_abs = max(worst_abs, diff)
        worst = max(worst, max(0.0, diff - abs_tol) / denom)
        n += 1
    ok = worst <= tol and missing == 0 and na_mismatch == 0
    status = "OK" if ok else "FAIL"
    print(f"{status} {name}: {n} numeric rows, max |diff| {worst_abs:.3e}, "
          f"max rel err (over {abs_tol:.0e} floor) {worst:.3e}, "
          f"missing {missing}, NA mismatches {na_mismatch} (tol {tol:.0e})")
    return ok

def main():
    prefix, refdir = sys.argv[1], sys.argv[2]
    vcftools = "--vcftools" in sys.argv
    vcftools_only = "--vcftools-only" in sys.argv
    ok = True
    pixy_pi = os.path.join(refdir, "pixy_bypass", "pixy_pi.txt")
    pixy_dxy = os.path.join(refdir, "pixy_bypass", "pixy_dxy.txt")
    if vcftools_only:
        pixy_pi = pixy_dxy = ""
    if os.path.exists(pixy_pi) and os.path.exists(prefix + "_pi.txt"):
        ok &= compare("pi vs pixy --bypass_invariant_check",
                      load_tsv(pixy_pi, [0, 1, 2], 4),
                      load_tsv(prefix + "_pi.txt", [0, 1, 2], 4), 1e-9)
    if os.path.exists(pixy_dxy) and os.path.exists(prefix + "_dxy.txt"):
        ok &= compare("dxy vs pixy --bypass_invariant_check",
                      load_tsv(pixy_dxy, [0, 1, 2, 3], 5),
                      load_tsv(prefix + "_dxy.txt", [0, 1, 2, 3], 5), 1e-9)
    pixy_theta = os.path.join(refdir, "pixy_bypass", "pixy_watterson_theta.txt")
    pixy_tajima = os.path.join(refdir, "pixy_bypass", "pixy_tajima_d.txt")
    if os.path.exists(pixy_theta) and os.path.exists(prefix + "_watterson_theta.txt"):
        ok &= compare("watterson_theta vs pixy",
                      load_tsv(pixy_theta, [0, 1, 2], 4),
                      load_tsv(prefix + "_watterson_theta.txt", [0, 1, 2], 4), 1e-9)
    if os.path.exists(pixy_tajima) and os.path.exists(prefix + "_tajima_d.txt"):
        ok &= compare("tajima_d vs pixy",
                      load_tsv(pixy_tajima, [0, 1, 2], 4),
                      load_tsv(prefix + "_tajima_d.txt", [0, 1, 2], 4), 1e-9)
    for method in ("hudson", "wc"):
        ref = os.path.join(refdir, f"pixy_fst_{method}", "pixy_fst.txt")
        reed_path = prefix + "_fst.txt"
        if not (os.path.exists(ref) and os.path.exists(reed_path)):
            continue
        # Only compare rows of the matching method (reed writes one method
        # per run; the header records which).
        with open(reed_path) as fh:
            header = fh.readline()
        if f"avg_{method}_fst" not in header:
            continue
        ok &= compare(f"fst ({method}) vs pixy",
                      load_tsv(ref, [0, 1, 2, 3], 5),
                      load_tsv(reed_path, [0, 1, 2, 3], 5), 1e-9)
    if vcftools:
        reed = load_tsv(prefix + "_pi.txt", [0, 1, 2], 4)
        paths = sorted(glob.glob(os.path.join(refdir, "vcftools_pi_*.windowed.pi")))
        paths += sorted(glob.glob(os.path.join(refdir, "vc_pi_*.windowed.pi")))
        for path in paths:
            base = os.path.basename(path)
            pop = base.split("_pi_", 1)[1][:-len(".windowed.pi")]
            expected = {}
            with open(path) as fh:
                fh.readline()
                for line in fh:
                    c = line.split()
                    expected[(pop, c[0], c[1])] = float(c[5])
            # vcftools prints ~6 significant digits.
            rounded = {k: (None if v is None else float(f"{v:.6g}")) for k, v in reed.items()}
            # Documented convention difference (04-algorithms.md §4): reed
            # clips the last window of a contig to the contig length, while
            # vcftools/pixy extend it to the full window size and count the
            # positions past the contig end as monomorphic. Those rows are
            # excluded here and reported.
            last_windows = {}
            for (p_, chrom, start) in rounded:
                if p_ == pop:
                    last_windows[chrom] = max(last_windows.get(chrom, 0), int(start))
            excluded = {k for k in rounded if k[0] == pop and int(k[2]) == last_windows.get(k[1])}
            expected = {k: v for k, v in expected.items() if k not in excluded}
            rounded = {k: v for k, v in rounded.items() if k not in excluded}
            print(f"note: excluded {len(excluded)} clipped last-window rows for {pop} "
                  f"(reed clips to the contig length; vcftools extends to the full window)")
            # vcftools applies a per-site finite-sample correction
            # (n/(n-1)) while reed/pixy use ratio-of-sums, so windows with
            # missing data differ by up to ~1e-5 relative (measured max
            # 1.05e-5 on the full file); vcftools also prints limited
            # precision. 1e-4 keeps the gate meaningful (a real algorithmic
            # error is orders of magnitude larger) without false failures.
            ok &= compare(f"pi vs vcftools --window-pi ({pop})", expected, rounded, 1e-4)
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
