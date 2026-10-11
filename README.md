# reed

[English](README.md) | [简体中文](README.zh-CN.md)

**reed is a Popgen API library written in Rust.**

## Installation

```bash
cargo install --path crates/reed-cli     # requires Rust 1.91+
reed --version
```

## Usage

```bash
# Genome-wide and windowed pi, dxy, Fst, theta and Tajima's D
reed stats allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out

# Variant-only input
reed stats clean.vcf.gz --pop pops.txt --variant-only \
    --fst --window 100k --step 50k -o out

# Per-site and per-sample QC reports plus the bcftools-stats summary, one pass
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# Site and genotype filters, with an index on the output
reed filter in.vcf.gz -o out.vcf.gz -Oz6 --geno 0.1 --maf 0.05 \
    --min-dp 5 --min-gq 20 -W

# LD pruning
reed filter in.vcf.gz --indep-pairwise 50kb 5 0.2 -Oz6 -o pruned.vcf.gz

# CSI index
reed index in.vcf.gz
```

`reed <command> --help` lists every option.

## reed stats

Pass any combination of statistics flags; each one adds its output file.

| Flag | Output |
|---|---|
| `--pi` | pi per population, `PREFIX_pi.txt` |
| `--dxy` | dxy per population pair, `PREFIX_dxy.txt` |
| `--fst [wc84\|hudson]` | Fst for all population pairs, `PREFIX_fst.txt` |
| `--fst-global` | Global WC84 Fst over all populations, `PREFIX_fst_global.txt` |
| `--theta` | Watterson's theta per population, `PREFIX_watterson_theta.txt` |
| `--tajima-d` | Tajima's D per population, `PREFIX_tajima_d.txt` |
| `--pbs A,B,C` | PBS for the named triple, `PREFIX_pbs.txt`; `--pbs-normed` for PBSn1 |
| `--sfs` | Projected SFS under `PREFIX_sfs/dadi/` and `PREFIX_sfs/fastsimcoal2/`; needs `--proj N` |
| `--f2 A,B`, `--f3 C;A,B`, `--f4 A,B;C,D`, `--dstat A,B;C,D`, `--f4-ratio A,B;C,D;X` | Patterson f-statistics with block jackknife, `PREFIX.fstats` |
| `--ld-decay [MAX_BP]` | Pairwise r² in 1 kb bins, `PREFIX.ld.decay`; `--ld-pairs` for the per-pair table `PREFIX.ld` |
| `--freq`, `--hardy`, `--het`, `--missing`, `--depth` | `PREFIX.afreq`, `PREFIX.hardy`, `PREFIX.het`, `PREFIX.smiss` + `PREFIX.vmiss`, `PREFIX.idepth` + `PREFIX.ldepth.mean` |
| `--summary` | QC summary in bcftools-stats layout, `PREFIX.stats` |

Window options:

| Flag | Meaning |
|---|---|
| `--window <SIZE>` | Window size in bp; accepts `100k`, `1M` (default 100000) |
| `--step <STEP>` | Window step (default: the window size) |
| `--levels <LIST>` | Rows to report: `window`, `contig`, `genome` (default `window`) |
| `--skip-empty` | Drop windows that no site fell into |
| `--pop <FILE>` | `sample<TAB>population` rows, no header; required by the windowed modes |
| `--variant-only` | Treat the input as all-variant: missing positions count as monomorphic |
| `--bypass-invariant-check` | Compute over the sites present (biases pi/dxy denominators) |
| `--contig-lengths <FILE>` | `.fai`-style `name<TAB>length`, for headers without `##contig` lines |

Sample and region selection works the same in `stats` and `filter`: `-r/-R`
regions (need an index), `-t/-T` targets (streamed), `-s/-S`, `--keep`,
`--remove`. Site filters (`-v/-V`, `--snps-only`, `-m/-M`, `-f`, `-i/-e`,
`-q/-Q`, `--geno`, `--maf/--max-maf`, `--mac/--max-mac`, `--hwe`) select what
the statistics see.

## reed filter

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

## reed index

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
`count_missing` for pi and dxy; `no_snps` for Fst, PBS and the f-statistics).
Estimates with a zero denominator print `NA`.

## Library use

The `reed` crate is the public API: typed requests go in, structured reports
come out, and the writers turn reports into the same files the CLI writes.

```toml
[dependencies]
reed = { path = "crates/reed" }
```

```rust
use reed::{Populations, StatKind, StatsRequest};

let request = StatsRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::File("pops.txt".into()),
    ..StatsRequest::default()
};
let report = reed::run_stats(Some("in.vcf.gz".as_ref()), request, &reed::budget::ThreadBudget::default())?;
let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false);
```

The facade re-exports everything needed: `run_filter`, `run_stats`,
`run_index`, the request types (`StatsRequest`, `FilterSpec`), the reports,
the writers (`reed::windows`, `reed::reports`, `reed::sfs`, `reed::fstats`,
`reed::ld`, `reed::summary`), and a single `reed::Error` whose `kind()`
separates usage errors from runtime failures. `StatsRequest::to_args()` and
`FilterSpec::to_args()` render a request back to the CLI spelling, so a tool
that builds requests programmatically can also log the equivalent command
line. Only this facade carries a stability promise; the lower crates under
`crates/` are internal.

## Citation

reed implements published estimators. Cite the methods it uses:

- Nei & Li 1979 for pi; Nei 1987 for dxy; Watterson 1975 for theta; Tajima
  1989 for D.
- Hudson 1992 and Bhatia et al. 2013 for Hudson Fst; Weir & Cockerham 1984
  for WC84 Fst.
- Patterson et al. 2012 for f-statistics; Yi et al. 2010 for PBS.
- Korunes & Samuk 2021 for the missing-data treatment of pi, dxy, theta and
  D.

## License

MIT OR Apache-2.0. See [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).
