# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed is a population-genetics command-line tool written in Rust. It reads VCF/BCF and computes π, dxy, Fst, θ, Tajima's D, PBS, SFS, f-statistics and LD decay by window or genome-wide. It also does hard filtering, LD pruning, QC reports, PCA/DAPC, between-individual matrices, format conversion and indexing.

## Installation

```bash
cargo install --path crates/reed-cli
reed --version
```

Requires Rust 1.91 or newer.

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

The full manual is in [`docsify_doc/README.md`](docsify_doc/README.md): installation, one page per command with every option, a step-by-step quick start, filter expressions, output conventions and the Rust API. Preview it locally with:

```bash
docsify serve docsify_doc
```

## Citation

Cite the papers behind the statistics in use; the full list is in [`docsify_doc/citation.md`](docsify_doc/citation.md).

## License

MIT OR Apache-2.0. See [LICENSE-MIT](LICENSE-MIT) and [LICENSE-APACHE](LICENSE-APACHE).
