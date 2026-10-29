# reed

[English](/#/en/) | [简体中文](/#/)

reed is a population-genetics tool written in Rust. It reads VCF/BCF and performs hard filtering, QC reports, windowed population-genetics statistics, principal components and between-individual relationship matrices, format conversion and indexing. Current version 0.4.0.

This is the English translation. When the two versions differ, the Chinese pages are authoritative.

The binary is `reed`, a wrapper around the `reed` crate; the same library can also be called directly from Rust code, see [Library use](library.md).

## Commands

| Command | Alias | Purpose |
|---|---|---|
| `filter` | | Hard filtering: site and sample selection, genotype masking, QC thresholds, LD pruning |
| `stats` | | Site-level and sample-level QC reports: allele frequency, HWE, heterozygosity, missingness, depth, histograms |
| `popg` | `popgen` | Windowed population statistics: π, dxy, Fst, θ, Tajima's D, PBS, SFS, f-statistics, LD decay |
| `pca` | | Principal components of the sample relationship matrix (GRM); `--dapc` discriminant analysis, `--xval` cross-validation |
| `distance` | | Between-individual matrices: allele-difference distance, IBS, GRM, KING, PI_HAT, relatedness, discordance; `--tree` builds a tree |
| `convert` | | Phylogenetic matrices, genotype matrices, hapmap, reference MSA, projected SFS; hapmap back to VCF |
| `rename` | `reheader` | Identifier edits: sample names, contig names, variant IDs, header lines |
| `index` | | Build CSI/TBI indexes |

## Conventions

1. Parameter names and semantics are identical across commands: sample subsets (`-s/-S/--keep/--remove`), region selection (`-r/-R/-t/-T`) and output (`-o/-O`) are spelled the same everywhere.
2. Threshold parameters state the quantity to keep: lower bounds for what to retain (`--maf/--mac/--hwe`), upper bounds for what to suppress (`--max-maf/--max-mac/--geno/--mind`). No `1−x` conversion happens on the command line.
3. Missingness thresholds (`--geno`, `--mind`) keep the boundary: a site or sample exactly at the threshold stays.
4. Malformed records are skipped and counted, never fatal; the counts go to the log or the JSON report.
5. All parameters are validated before the first record; errors exit with code 2 and leave no partial output.

## Input and output conventions

Input: VCF, bgzip-compressed VCF, BCF, or stdin (`-`). Region queries (`-r/-R`) need a `.csi`/`.tbi` index next to the input (built by `reed index`); position-based selection (`-t/-T`) needs no index.

Output: `filter`, `rename`, `convert` and `index` treat `-o` as the output file name; `stats`, `popg`, `pca` and `distance` treat it as an output prefix and append suffixes (`.pi.txt`, `.eigenvec`, ...). Estimates with a zero denominator are written as `NA`, not 0. Coordinates are 1-based closed intervals (BED input excepted, read as 0-based half-open).

`--threads N` controls worker threads (default: all cores). `--memory MB` sets the in-flight memory budget in MiB (default: half of physical memory). Length parameters accept `100000`, `100k`, `100kb`, `1M` and `1Mb`; suffixes are case-insensitive.
