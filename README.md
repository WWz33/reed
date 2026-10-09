# reed

[English](README.md) | [简体中文](README.zh-CN.md)

**reed** is a multithreaded VCF/BCF toolkit for population-genetics
statistics and hard filtering, in pure Rust. It is a command-line product:
the `reed` binary is the deliverable and the crates are its internal
structure. Only the CLI flags and the output file formats are a stability
promise (doc_v5 08 §6).

No htslib, no GPL: transport is built on noodles, plus reed's own BGZF and
BCF codecs.

## Commands

| Command | What it does | Reference |
|---|---|---|
| `reed filter` | Hard filtering, genotype masking, INFO recomputation; byte-faithful VCF/BCF output | bcftools `view` + `filter -S .` + `+fill-tags` + `annotate -x`; plink2 QC thresholds |
| `reed stats` | Every statistic in one command: π, dxy, Fst, Watterson θ, Tajima's D, PBS, SFS, per-site/per-sample reports, bcftools-stats summary | vcftools, pixy, plink2, easySFS, `bcftools stats` |
| `reed index` | CSI (default) or TBI | `bcftools index` |

`reed <command> --help` is authoritative for the flags; the full parameter
specification is `docs/spec/cli.md` (the merged CLI document).

```bash
# windowed statistics over an all-sites input, 56 threads
reed stats allsites.vcf.gz --pop pops.txt --pi --dxy --fst --tajima-d --window 10k -o out

# variant-only input; genome-wide and windowed Fst
reed stats clean.vcf.gz --pop pops.txt --variant-only --pi --fst --window 100k --step 50k -o out

# per-site/per-sample reports plus the bcftools-stats summary, one pass
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# hard filter with a post-mask site gate and an index
reed filter in.vcf.gz -o out.vcf.gz -Oz6 --geno 0.1 --maf 0.05 --min-dp 5 --min-gq 20 -W

# index
reed index in.vcf.gz
```

`-o` is a file for `filter` and `index`, and an output **prefix** for
`stats` (each statistic writes its own file). Statistics parameters are
named after the method (`--pi`, `--tajima-d`, `--fst wc84|hudson`); software
names appear only in the parity suites.

## Verification

Correctness is pinned to reference tools, not to self-consistency. The
numbers below are from the doc_v5 review (2026-10-09); each has a
reproducible script under `gold/`.

| Check | Reference | Result |
|---|---|---|
| π, dxy (with invariant sites) | pixy 2.2.4 | 28,470 windows, max relative error 4.8e-14 |
| Watterson θ, Tajima's D | pixy | 3.1e-14 / 8.9e-11 |
| Fst wc84 / hudson (pairwise) | pixy | 2.9e-15 / 2.2e-15 |
| PBS | scikit-allel | per-column exact |
| `--freq/--hardy/--het/--missing` | plink2 | row-for-row |
| `--depth` | vcftools | row-for-row |
| `--summary` | `bcftools stats` | nine sections row-for-row |
| gold filter rule | bcftools | 1,768,759 sites, `check_parity` PASS, body md5 `7aedc59a…` |
| `--geno/--mind` boundaries | plink2, vcftools | N=10/7 and half-calls agree; subsets 263,420 / 381,672 |
| BCF vs VCF, mmap vs streaming | bcftools | byte-identical bodies; identical CSI |

Determinism is a gate: output is byte-identical across thread counts
(single-pass / two-pass / sketch strategies too).

## Performance

Full `clean.vcf.gz` (7.0 GB, 3,315,604 records, 348 samples), 56 threads,
warm cache, median of three runs:

| Command | wall | 
|---|---|
| `filter` gold rule, VCF input, `-Oz6` | 8.79 s |
| `filter` gold rule, BCF input, `-Oz6` | 8.87 s |
| `stats` π/dxy/θ/D/Fst, VCF / BCF | 1.74 s / 1.51 s |
| `stats` reports + summary | 3.76 s |
| `index` VCF.gz / BCF | 1.60 s / 1.05 s |

The reader is memory-mapped for regular files; deflate dominates `filter`.

## Repository layout

```
crates/    the eight crates; the only binary is reed-cli
gold/      verification: fixtures and parity scripts against the reference tools
bench/     measurement suites and result logs
docs/      the specifications (spec/ = current, rounds/ = history)
ci/        drift guards and the CI entry point
```

## Build and test

```bash
cargo test --release --workspace
cargo build --release --workspace --all-targets
bash ci/check-doc-drift.sh
```

## License

MIT OR Apache-2.0.
