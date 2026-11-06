# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed is a population-genetics command-line tool written in Rust. It reads VCF/BCF and computes π, dxy, Fst, θ, Tajima's D, PBS, SFS, f-statistics and LD decay by window or genome-wide. It also does hard filtering, LD pruning, QC reports, PCA/DAPC, between-individual matrices, format conversion and indexing.

## Installation

Download the archive matching your platform from [GitHub Releases](https://github.com/WWz33/reed/releases). These
binaries do not require Rust, Cargo, htslib or bcftools at runtime.

| Platform | Archive suffix | Minimum system |
|---|---|---|
| Linux x86_64 | `x86_64-unknown-linux-gnu-glibc217.tar.gz` | glibc 2.17+, kernel 3.2+ |
| Linux ARM64 | `aarch64-unknown-linux-gnu-glibc217.tar.gz` | glibc 2.17+, kernel 4.1+ |
| macOS Apple Silicon | `aarch64-apple-darwin-macos11.tar.gz` | macOS 11.0 |
| macOS Intel | `x86_64-apple-darwin-macos11.tar.gz` | macOS 11.0 |

Archives are named `reed-v<VERSION>-<SUFFIX>`. Download `SHA256SUMS` with the
archive, verify it with `sha256sum --ignore-missing --check SHA256SUMS` on Linux
or `shasum -a 256 --ignore-missing --check SHA256SUMS` on macOS, then extract it
and place `reed` on your PATH. Each archive includes the licenses and
`BUILD-INFO.txt` with the target, source commit and compiler version. The macOS
binaries are unsigned and are not notarized.

To install from a checkout instead:

```bash
cargo install --locked --path crates/reed-cli
reed --version
```

Source builds support Rust 1.91 or newer and use the Rust toolchain active
in your environment. A C compiler is needed for bundled libdeflate; no separately
installed libdeflate or htslib is needed. CI tests current stable and checks the
1.91.1 MSRV; tagged release binaries are built with 1.91.1.

`--threads N` accepts 1–1024. Release builds use baseline CPU instructions;
they do not require the build runner's CPU model.

## Quick start

```bash
# index the input
reed index raw.vcf.gz

# hard filter: mask shallow genotypes, then site thresholds
reed filter raw.vcf.gz -o clean.vcf.gz -Oz6 --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 -W

# windowed statistics on a variant-only file
reed popg clean.vcf.gz --pop pops.txt --variant-only --pi --fst --window 100k -o out

# QC report and histograms
mkdir -p qc
reed stats raw.vcf.gz --freq --missing --depth --hist qual,dp-site,gq -o qc/raw

# LD pruning, then PCA
reed filter clean.vcf.gz --indep-pairwise 50 5 0.2 -Oz6 -o pruned.vcf.gz
reed pca pruned.vcf.gz --pcs 20 -o pcs

# between-individual matrices and a tree
reed distance clean.vcf.gz --type grm -o rel
reed distance clean.vcf.gz --type 1-ibs --tree nj --bootstrap 100 -o tree
```

## Documentation

The full manual is in [`docsite/`](docsite/): installation, one page per command with every option, a step-by-step quick start, filter expressions, output conventions and the Rust API. It is built with Astro Starlight and published at <https://wwz33.github.io/reed/>. Preview it locally with:

```bash
cd docsite && npm install && npm run dev
```

## Citation

Cite the papers behind the statistics in use; the full list is in [`docsite/src/content/docs/citation.md`](docsite/src/content/docs/citation.md).

## License

MIT OR Apache-2.0. See [LICENSE-MIT](LICENSE-MIT) and [LICENSE-APACHE](LICENSE-APACHE).
