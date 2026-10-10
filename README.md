# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed computes population-genetics statistics and applies hard filters to
VCF, VCF.gz and BCF files, splitting the work over every core. It is written
in pure Rust: the transport comes from noodles, the BGZF and BCF codecs are
reed's own, and no GPL code is linked. The `reed` binary is the product; the
crates under `crates/` are its internal structure. Only the CLI flags and the
output file formats carry a stability promise.

## Installation

reed needs Rust 1.91 or newer. The repository pins 1.99 in
`rust-toolchain.toml`.

```bash
cargo install --path crates/reed-cli
reed --version   # reed 0.3.0 (64e5dce 2026-10-10, rustc 1.99.0)
```

Building in place works too:

```bash
cargo build --release -p reed-cli     # binary at target/release/reed
```

There are no system libraries to install. The BGZF and BCF codecs are
compiled from Rust sources, so the binary runs anywhere the toolchain does.

## Quick start

```bash
# Genome-wide and windowed π, dxy, Fst, θ and Tajima's D over an all-sites file
reed stats allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out

# Variant-only input, pairwise and global Fst
reed stats clean.vcf.gz --pop pops.txt --variant-only \
    --fst --window 100k --step 50k -o out

# Per-site and per-sample QC reports plus the bcftools-stats summary, one pass
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# Site and genotype filters, with an index on the output
reed filter in.vcf.gz -o out.vcf.gz -Oz6 --geno 0.1 --maf 0.05 --min-dp 5 --min-gq 20 -W

# LD pruning: one command writes the pruned VCF
reed filter in.vcf.gz --indep-pairwise 50kb 5 0.2 -Oz6 -o pruned.vcf.gz

# CSI index
reed index in.vcf.gz
```

## Commands

| Command | Purpose | Reference |
|---|---|---|
| `reed filter` | Hard filtering, genotype masking, INFO recomputation; byte-faithful VCF/BCF output | bcftools `view`, `filter -S .`, `+fill-tags`, `annotate -x`; plink2 QC thresholds |
| `reed stats` | π, dxy, Fst, Watterson θ, Tajima's D, PBS, SFS, f-statistics, LD, per-site and per-sample reports, bcftools-stats summary | vcftools, pixy, plink2, easySFS, scikit-allel, `bcftools stats` |
| `reed index` | CSI (default) or TBI | `bcftools index` |

`reed <command> --help` lists every option. The full parameter
specification is `docs/spec/cli.md`.

### reed stats

Statistics flags. Pass any combination; each one adds its output file.

| Flag | Output |
|---|---|
| `--pi` | π per population, `PREFIX_pi.txt` |
| `--dxy` | dxy per population pair, `PREFIX_dxy.txt` |
| `--fst [wc84\|hudson]` | Fst for all population pairs, `PREFIX_fst.txt` (the method shows in the column header) |
| `--fst-global` | Global WC84 Fst over all populations, `PREFIX_fst_global.txt` |
| `--theta` | Watterson's θ per population, `PREFIX_watterson_theta.txt` |
| `--tajima-d` | Tajima's D per population, `PREFIX_tajima_d.txt` |
| `--pbs A,B,C` | PBS for the named triple, `PREFIX_pbs.txt`; add `--pbs-normed` for PBSn1 |
| `--sfs` | Projected SFS under `PREFIX_sfs/dadi/` and `PREFIX_sfs/fastsimcoal2/`; needs `--proj N` |
| `--f2 A,B`, `--f3 C;A,B`, `--f4 A,B;C,D`, `--dstat A,B;C,D`, `--f4-ratio A,B;C,D;X` | Patterson f-statistics with block jackknife, `PREFIX.fstats` |
| `--ld-decay [MAX_BP]` | Pairwise r² in 1 kb bins, `PREFIX.ld.decay`; add `--ld-pairs` for the per-pair table `PREFIX.ld` |
| `--freq`, `--hardy`, `--het`, `--missing`, `--depth` | `PREFIX.afreq`, `PREFIX.hardy`, `PREFIX.het`, `PREFIX.smiss` + `PREFIX.vmiss`, `PREFIX.idepth` + `PREFIX.ldepth.mean` |
| `--summary` | QC summary in bcftools-stats layout, `PREFIX.stats` |

Window flags:

| Flag | Meaning |
|---|---|
| `--window <SIZE>` | Window size in bp; accepts `100k`, `1M` (default 100000) |
| `--step <STEP>` | Window step (default: the window size) |
| `--levels <LIST>` | Rows to report: `window`, `contig`, `genome` (default `window`) |
| `--skip-empty` | Drop windows that no site fell into |
| `--pop <FILE>` | `sample<TAB>population` rows, no header; required by the windowed modes and by every per-population statistic |
| `--variant-only` | Treat the input as all-variant: positions the file lacks count as monomorphic and fully called |
| `--bypass-invariant-check` | Compute over the sites present, accepting the π/dxy denominator bias |
| `--contig-lengths <FILE>` | `.fai`-style `name<TAB>length` pairs, for headers without `##contig` lines |

`--variant-only` and `--bypass-invariant-check` exist because π, dxy, θ and
Tajima's D need a denominator that includes monomorphic sites. By default
reed verifies that the input carries invariant records and stops with an
error when it does not.

Sample and region selection works the same in `stats` and `filter`: `-r/-R`
regions (need an index), `-t/-T` targets (streamed), `-s/-S`, `--keep`,
`--remove`. Site filters (`-v/-V`, `--snps-only`, `-m/-M`, `-f`, `-i/-e`,
`-q/-Q`, `--geno`, `--maf/--max-maf`, `--mac/--max-mac`, `--hwe`) select what
the statistics see.

### reed filter

Selection, from bcftools: `-r/-R`, `-t/-T`, `-s/-S`, `--keep`, `--remove`
(with `--force-samples` to downgrade unknown names to warnings), `-v/-V`,
`--snps-only`, `-m/-M`, `-f`, `-i/-e`.

Site filters, applied after masking: `-q/-Q/--min-af/--max-af` with an
optional `:nref|alt1|minor|major|nonmajor` suffix, `--geno`, `--maf`,
`--max-maf`, `--mac`, `--max-mac`, `--hwe <P [midp]>`.

Genotype masking: `--min-dp`, `--max-dp`, `--min-gq` mask from FMT/DP and
FMT/GQ; `-i/-e` FORMAT expressions with `--set-gts .` mask exactly the
matches, as `bcftools filter -S .` does. Masked genotypes count as missing
for the site filters that follow.

Sample filters: `--mind <FLOAT>` drops samples whose missing rate exceeds the
threshold. `--indep-pairwise <W[kb] [S] R2>` prunes sites by LD using the
plink2 algorithm and emits the kept sites; the window `W` may carry a `kb`
suffix, the step `S` defaults to 1, and `R2` is the r² ceiling.

INFO handling: reed updates `INFO/AC` and `INFO/AN` when the sample set
changes; `-I/--no-update` turns that off, `-x/--remove-info LIST` drops tags,
and `--fill-tags LIST` recomputes them.

Output: `-O v|z|b|u[0-9]` selects plain VCF, BGZF VCF, BCF or uncompressed
BCF, with an optional compression level. `-o` names the file (stdout by
default), `-W[=csi|tbi]` writes an index, `--report` writes a JSON run
report, `--timing` records per-phase timings, `--no-version` and
`--no-command` suppress the `##reed_*` provenance header lines.

### reed index

`-c/--csi` (default) or `-t/--tbi`, `-o` to write elsewhere, `-f` to
overwrite. The input must be a file: an index needs a seekable source.

## Input requirements

- VCF, VCF.gz or BCF, detected from the first bytes. `-` reads VCF text from
  stdin; `reed index` requires a file.
- The header must declare `##contig=<ID=…,length=…>`. Windowed statistics
  clip the last window of a contig to the declared length, and the index
  needs the contig order. Pass `--contig-lengths FILE` when the header lacks
  the lines.
- Records must be sorted by contig and position. Indexing and windowed
  statistics rely on it.
- `-r/-R` needs a `.csi` or `.tbi` next to the input, built with
  `reed index`. Plain scans do not.
- Genotypes must be diploid. reed rejects multiallelic records by default
  (`-M 2`), and treats half-calls such as `0/.` at the allele level, so
  `--geno` and `--mind` count them the way plink2 does.

## Output

`filter` and `index` write the file named by `-o`. `stats` treats `-o` as a
prefix and writes one file per statistic, listed in the tables above. Without
`-o`, window statistics print TSV to stdout and the report files are skipped.

Provenance header lines land in the output of `filter` (`##reed_version`,
`##reed_command`, `##reed_filterCommand`, `##reed_filter`) and inside
`PREFIX.log` for `stats`.

Windowed tables share a shape: one row per population (or population pair)
per window, with `chromosome`, `window_pos_1`, `window_pos_2`, the estimate,
and the counts behind it (`no_sites`, `count_diffs`, `count_comparisons`,
`count_missing` for π and dxy; `no_snps` for Fst, PBS and the f-statistics).
Estimates with a zero denominator print `NA`.

## Verification

Correctness is pinned to reference tools rather than to self-consistency.
Each row below has a script under `gold/`, and the numbers come from the
2026-10-10 review recorded in `docs/rounds/v6/`.

| Check | Reference | Result |
|---|---|---|
| π, dxy with invariant sites | pixy 2.2.4 | 28,470 windows, max relative error 4.8e-14 |
| Watterson θ, Tajima's D | pixy | 3.1e-14 / 8.9e-11 |
| Fst wc84 / hudson, pairwise | pixy | 2.9e-15 / 2.2e-15 |
| Global WC84 Fst | vcftools weighted Fst | 2.2e-16 |
| PBS | scikit-allel | per-column exact |
| Projected SFS | easySFS, dadi 2.4.4 | 1e-15 on the fixture, 5.7e-11 on the subset |
| f2/f3/f4/D estimates | scikit-allel | 1.3e-16 |
| LD r² and pruning sets | plink2 | identical kept sites on the full file (235,403) |
| `--freq/--hardy/--het/--missing` | plink2 | row for row |
| `--depth` | vcftools | row for row |
| `--summary` | `bcftools stats` | nine sections, row for row |
| gold filter rule | bcftools | 1,768,759 sites, `check_parity` PASS, body md5 `7aedc59a…` |
| `--geno/--mind` boundaries | plink2, vcftools | agree at N=10/7 and on half-calls; subsets 263,420 / 381,672 |
| BCF against VCF, mmap against streaming | bcftools | byte-identical bodies, identical CSI |

Thread count does not change the output: the test suite hashes the same
results at 1, 2, 3, 7 and 56 threads, across all three execution strategies.

## Performance

Full `clean.vcf.gz` (7.0 GB, 3,315,604 records, 348 samples), 56 threads,
warm cache, median of three runs:

| Command | Wall time |
|---|---|
| `filter` gold rule, VCF input, `-Oz6` | 8.85 s |
| `filter` gold rule, BCF input, `-Oz6` | 9.00 s |
| `filter` gold rule, `-Ob` output | 9.98 s |
| `stats` π/dxy/θ/D/Fst, VCF / BCF | 1.78 s / 1.65 s |
| `stats` reports + summary | 3.0 s to 5.4 s |
| `index` VCF.gz / BCF | 1.63 s / 1.12 s |

The reader memory-maps regular files and cuts compressed chunks on BGZF block
boundaries, so threads hand work over without resynchronising. Deflate
dominates `filter`. On a small machine, `--memory <MB>` caps the in-flight
buffers.

## Limitations

- Multiallelic records are rejected. reed does not recompute per-allele
  statistics when `-M` is raised past 2, so the default keeps them out.
- LD pruning lives in `filter` and needs its second pass. `reed stats`
  rejects `--indep-pairwise` and points at `reed filter`.
- `--ld-decay` asks for `--pop FILE` even though LD does not depend on
  population labels.
- `--window` currently requires an explicit `--step`; the documented default
  of "window size" is not accepted.
- f-statistics standard errors use block jackknife over base pairs, while
  scikit-allel uses variant counts. The estimates match; the errors are not
  comparable to allel's.
- Tajima's D divides by the variance at the window-average allele count, the
  pixy convention, which differs from vcftools' fixed 2N.

## Build and test

```bash
cargo test --release --workspace
cargo build --release --workspace --all-targets
cargo clippy --release --workspace --all-targets -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
cargo fmt --all --check
bash ci/check-doc-drift.sh
```

The lint gate covers the library and binary targets. Test targets carry
pedantic lints that are still open, so `-D warnings` over `--all-targets`
reports them.

## Citation

reed implements published estimators and checks its output against the
reference tools above. Cite the methods it uses:

- Nei & Li 1979 for π; Nei 1987 for dxy; Watterson 1975 for θ; Tajima 1989
  for D.
- Hudson 1992 and Bhatia et al. 2013 for Hudson Fst; Weir & Cockerham 1984
  for WC84 Fst.
- Patterson et al. 2012 for f-statistics; Yi et al. 2010 for PBS.
- Korunes & Samuk 2021 for the missing-data treatment of π, dxy, θ and D.

## Repository layout

```
crates/    eight crates; reed-cli builds the only binary
gold/      fixtures and parity scripts for the reference tools
bench/     measurement suites and result logs
docs/      specifications (spec/ current, rounds/ history)
ci/        drift guards and the CI entry point
```

## License

MIT OR Apache-2.0. See [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).
