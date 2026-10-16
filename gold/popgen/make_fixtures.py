#!/usr/bin/env python3
"""Deterministic fixtures shared by the parity scripts (doc_v4 05 §4).

Every fixture is small, hand-checkable and fixed (no RNG): the parity
scripts (`run_easysfs.sh`, `run_plink2_ld.sh`, `compare_fstats.py`, the
missing-rate boundary gates) all read these instead of each generating its
own.

Usage:
    python3 make_fixtures.py [--out DIR]      # default: gold/popgen/fixtures
"""
from __future__ import annotations

import argparse
import os

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
"""


def write(path: str, header: str, rows: list[tuple[int, str, str, list[str]]]) -> None:
    """Writes one fixture: (pos, ref, alt, genotypes) rows."""
    with open(path, "w") as fh:
        fh.write(header)
        for pos, ref, alt, gts in rows:
            fh.write(f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")


def samples_header(names: list[str]) -> str:
    return HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n"


def write_pops(path: str, assignments: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for sample, pop in assignments:
            fh.write(f"{sample}\t{pop}\n")


def missing_rate(out: str) -> None:
    """10 samples; record `pos` has `pos−1` fully-missing samples (call
    rates 1.0, 0.9, 0.8, 0.7, 0.5, 0.1, 0.0), plus half calls (doc_v4 02
    §7.1: `--geno X` == plink2 `--geno X` == vcftools `--max-missing 1−X`)."""
    names = [f"s{i}" for i in range(10)]
    rows = []
    for pos, missing in [(1, 0), (2, 1), (3, 2), (4, 3), (5, 5), (6, 9), (7, 10)]:
        gts = ["0/0"] * (10 - missing) + ["./."] * missing
        rows.append((pos, "A", "G", gts))
    # Half calls: 0/. contributes one called allele (vcftools semantics).
    rows.append((8, "A", "G", ["0/.", "0/.", "./."] + ["0/0"] * 7))
    rows.append((9, "A", "G", ["0/."] + ["./."] * 9))
    rows.append((10, "A", "G", ["0/."] * 10))
    write(os.path.join(out, "missing_rate.vcf"), samples_header(names), rows)


def sfs(out: str) -> None:
    """4 diploid samples, one population; sites with ALT counts 1,1,2,2,3,3,
    4,4,5,6,7 of 8 chromosomes (the doc_v4 05 §2 projected-SFS fixture)."""
    names = ["s0", "s1", "s2", "s3"]
    rows = [
        (1, "A", "G", ["0/0", "0/0", "0/0", "0/1"]),
        (2, "A", "G", ["0/1", "0/0", "0/0", "0/0"]),
        (3, "A", "G", ["0/0", "0/0", "0/1", "0/1"]),
        (4, "A", "G", ["0/1", "0/0", "0/1", "0/0"]),
        (5, "A", "G", ["0/1", "0/1", "0/1", "0/0"]),
        (6, "A", "G", ["1/1", "0/1", "0/0", "0/0"]),
        (7, "A", "G", ["1/1", "0/1", "0/1", "0/0"]),
        (8, "A", "G", ["1/1", "1/1", "0/0", "0/0"]),
        (9, "A", "G", ["1/1", "1/1", "1/0", "0/0"]),
        (10, "A", "G", ["1/1", "1/1", "1/1", "0/0"]),
        (11, "A", "G", ["1/1", "1/1", "1/1", "1/0"]),
    ]
    write(os.path.join(out, "sfs.vcf"), samples_header(names), rows)
    write_pops(os.path.join(out, "sfs.pops"), [(n, "pop1") for n in names])
    # BED windows with an overlapping pair (rows 2 and 3 share position 10),
    # for the pixy --bed_file parity (`run_pixy_bed.sh`).
    with open(os.path.join(out, "sfs.bed"), "w") as fh:
        fh.write("chr1\t0\t5\nchr1\t5\t10\nchr1\t9\t13\n")
    # The joint (2D) case: the same sites split s0,s1 | s2,s3, projection 4,4.
    write_pops(
        os.path.join(out, "sfs2.pops"),
        [("s0", "popA"), ("s1", "popA"), ("s2", "popB"), ("s3", "popB")],
    )


def multiallelic(out: str) -> None:
    """Multi-allelic sites for the per-ALT recompute parity (doc_v4 05 §3 /
    W6): a biallelic control, two ALTs with different counts, a three-ALT
    site, and one with a missing genotype. `AC`/`AN` in INFO are deliberately
    *stale* (the point of `--fill-tags`)."""
    names = ["s0", "s1", "s2", "s3"]
    header = samples_header(names).replace(
        "##FORMAT=<ID=GQ", "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred likelihoods\">\n##FORMAT=<ID=GQ"
    )
    with open(os.path.join(out, "multiallelic.vcf"), "w") as fh:
        fh.write(header)
        rows = [
            (1, "A", "G", ["0/1", "0/0", "1/1", "0/0"]),
            (2, "A", "G,T", ["0/1", "0/2", "1/2", "0/0"]),
            (3, "A", "G,T", ["1/1", "1/2", "2/2", "./."]),
            (4, "A", "C,G,T", ["0/1", "2/3", "1/1", "0/2"]),
        ]
        for pos, ref, alt, gts in rows:
            fh.write(
                f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tAC=99;AN=999\tGT\t"
                + "\t".join(gts)
                + "\n"
            )


def fstats(out: str) -> None:
    """3 populations (2 samples each) for the Patterson f-stat parity against
    scikit-allel (`admixture.py`): fixed, hand-checkable allele patterns."""
    names = [f"s{i}" for i in range(6)]
    rows = [
        (1, "A", "G", ["0/0", "0/0", "0/1", "0/1", "1/1", "1/1"]),
        (2, "A", "G", ["0/1", "0/1", "0/0", "0/1", "1/1", "0/1"]),
        (3, "A", "G", ["0/0", "0/1", "0/1", "0/1", "0/0", "0/1"]),
        (4, "A", "G", ["1/1", "0/1", "0/0", "0/0", "0/1", "0/0"]),
        (5, "A", "G", ["0/1", "0/0", "0/1", "0/0", "0/1", "0/1"]),
    ]
    write(os.path.join(out, "fstats.vcf"), samples_header(names), rows)
    write_pops(
        os.path.join(out, "fstats.pops"),
        [("s0", "popA"), ("s1", "popA"), ("s2", "popB"), ("s3", "popB"), ("s4", "popC"), ("s5", "popC")],
    )


def ld(out: str) -> None:
    """6 samples, 6 sites in two blocks of linked sites (for the plink2
    `--r2-unphased` / `--indep-pairwise` parity)."""
    names = [f"s{i}" for i in range(6)]
    rows = [
        (1, "A", "G", ["0/1", "0/1", "0/0", "1/1", "0/1", "0/0"]),
        (2, "A", "G", ["0/1", "0/1", "0/0", "1/1", "0/1", "0/0"]),
        (3, "A", "G", ["0/0", "0/1", "0/1", "0/0", "0/1", "1/1"]),
        (4, "A", "G", ["0/0", "0/1", "0/1", "0/0", "0/1", "1/1"]),
        (5, "A", "G", ["0/0", "0/0", "0/1", "0/1", "0/0", "0/1"]),
        (6, "A", "G", ["1/1", "0/1", "0/0", "0/1", "0/0", "0/1"]),
    ]
    write(os.path.join(out, "ld.vcf"), samples_header(names), rows)
    write_pops(os.path.join(out, "ld.pops"), [(n, "pop1") for n in names])



def phylo(out: str) -> None:
    """The `reed convert` phylogenetic fixture (doc_v13 01 §7): 8 samples x
    20 records covering biallelic/multiallelic SNPs, MNPs, indels,
    `*`/`<NON_REF>`/`.` ALT values, phased and haploid genotypes, half
    calls, a fully missing sample and one malformed column-count line.
    Reference tools: vcf2phylip 2.9, vcftools 0.1.18, TASSEL 5, vcf2msa."""
    names = [f"s{i}" for i in range(8)]
    header = HEADER + "##contig=<ID=chr2,length=100000>\n"
    rows = [
        ("chr1", 1, ".", "A", "G", ["0/0", "0/0", "0/0", "0/1", "0/1", "1/1", "1/1", "./."]),
        ("chr1", 2, ".", "A", "G", ["0/1", "0/0", "0/0", "0/1", "0/1", "1/1", "1/1", "./."]),
        ("chr1", 3, ".", "A", "G", ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "1/1", "./."]),
        ("chr1", 4, ".", "C", "T", ["0/0", "0/0", "0/0", "0/0", "0/1", "1/1", "1/1", "./."]),
        ("chr1", 5, ".", "A", "G,T", ["0/0", "0/1", "2/2", "1/2", "0/2", "1/1", "2/2", "./."]),
        ("chr1", 6, ".", "AT", "GC", ["0/0", "0/1", "1/1", "0/1", "1/1", "0/0", "0/1", "./."]),
        ("chr1", 7, ".", "A", "GT", ["0/0", "0/1", "1/1", "0/1", "1/1", "0/0", "0/1", "./."]),
        ("chr1", 8, ".", "A", "G,*", ["0/0", "0/1", "1/1", "2/2", "0/2", "1/1", "2/2", "./."]),
        ("chr1", 9, ".", "A", "<NON_REF>", ["0/0", "0/0", "0/0", "0/1", "0/1", "1/1", "1/1", "./."]),
        ("chr1", 10, ".", "A", "G", ["0/0"] * 7 + ["./."]),
        ("chr1", 11, ".", "A", "G", ["0|1", "1|0", "0|0", "1|1", "0/0", "0/1", "1/1", "./."]),
        ("chr1", 12, ".", "A", "G", ["0", "0", "1", "1", "0", "1", "0", "1"]),
        ("chr1", 14, ".", "A", "G", ["0/.", "0/0", "0/0", "0/0", "0/1", "1/1", "1/1", "./."]),
        ("chr1", 15, ".", "A", ".", ["0/0"] * 7 + ["./."]),
        ("chr2", 1, ".", "A", "G", ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "./."]),
        ("chr2", 2, ".", "C", "A", ["0/0", "0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "./."]),
        ("chr2", 3, ".", "T", "G", ["0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1", "./."]),
        ("chr2", 4, ".", "A", "C", ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1", "0/0", "./."]),
        ("chr2", 5, ".", "G", "T", ["0/1", "0/1", "0/1", "0/1", "0/1", "0/1", "0/1", "./."]),
    ]
    with open(os.path.join(out, "phylo.vcf"), "w") as fh:
        fh.write(header)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
        for chrom, pos, ident, ref, alt, gts in rows:
            fh.write(f"{chrom}\t{pos}\t{ident}\t{ref}\t{alt}\t50\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
        # One malformed line: two samples instead of eight (vcf2phylip skips
        # it; reed counts it in the convert counters, B13-2).
        fh.write("chr1\t16\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\n")


def convert012(out: str) -> None:
    """The `--012` fixture: full diploid calls plus the B13-8 rows (haploid
    and half calls), a multiallelic site (skipped) and an `ALT='.'` site."""
    names = [f"s{i}" for i in range(4)]
    with open(os.path.join(out, "convert012.vcf"), "w") as fh:
        fh.write(HEADER)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
        rows = [
            (1, "A", "G", ["0/0", "0/1", "1/1", "./."]),
            (2, "A", "G", ["1/1", "1/0", "0/1", "0/0"]),
            (3, "A", "G,T", ["0/0", "0/1", "2/2", "./."]),
            (4, "A", ".", ["0/0", "0/0", "0/0", "0/0"]),
            (5, "A", "G", ["0", "1", "0/.", "./1"]),
        ]
        for pos, ref, alt, gts in rows:
            fh.write(f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")


def convert_msa(out: str) -> None:
    """The `--msa` fixture: a 10 bp reference and three unphased SNP sites
    (vcf2msa's supported domain)."""
    with open(os.path.join(out, "convert_msa.fa"), "w") as fh:
        fh.write(">chr1\nAACCGGTTAA\n")
    names = [f"s{i}" for i in range(3)]
    with open(os.path.join(out, "convert_msa.vcf"), "w") as fh:
        fh.write(HEADER)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
        for pos, ref, alt, gts in [
            (2, "A", "G", ["0/1", "0/0", "1/1"]),
            (5, "G", "T", ["0/0", "0/1", "0/0"]),
            (8, "T", "C", ["1/1", "0/1", "0/0"]),
        ]:
            fh.write(f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")

def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=os.path.join(here, "fixtures"))
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    missing_rate(args.out)
    sfs(args.out)
    multiallelic(args.out)
    fstats(args.out)
    ld(args.out)
    phylo(args.out)
    convert012(args.out)
    convert_msa(args.out)
    print(f"fixtures written to {args.out}")
    for name in sorted(os.listdir(args.out)):
        print(f"  {name}")


if __name__ == "__main__":
    main()
