# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed computes population-genetics statistics on VCF/BCF files and applies
hard filters. No system libraries, multi-threaded, with the BGZF and BCF
codecs compiled in; also a library, crate name `reed`.

## Installation

```bash
cargo install --path crates/reed-cli     # requires Rust 1.91+
reed --version
```

## Usage

```bash
# Genome-wide and windowed pi, dxy, Fst, theta and Tajima's D
reed popg allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out

# Variant-only input
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --fst --window 100k --step 50k -o out

# Per-site and per-sample QC reports plus the bcftools-stats summary, one pass
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# Histograms in the same pass: built-ins and any numeric INFO/FORMAT tag
reed stats clean.vcf.gz --hist qual,dp-site,gq,INFO/DP -o qc

# Identifier edits: sample names, contig names, variant IDs, header lines
reed rename in.vcf.gz -s "s0 sampleA" --rename-chrs chrs.txt \
    --set-missing-var-ids '%CHROM:%POS' -Oz -o renamed.vcf.gz

# Site and genotype filters, with an index on the output
reed filter in.vcf.gz -o out.vcf.gz -Oz6 --geno 0.1 --maf 0.05 \
    --min-dp 5 --min-gq 20 -W

# LD pruning
reed filter in.vcf.gz --indep-pairwise 50 5 0.2 -Oz6 -o pruned.vcf.gz

# CSI index
reed index in.vcf.gz

# Phylogenetic matrices, the 0/1/2 matrix, hapmap and a reference-driven MSA
reed convert in.vcf.gz --phylip --fasta --nexus --nexus-binary --012 --hapmap -o matrices
reed convert in.vcf.gz --msa --reference ref.fa --region chr1:1-100000 -o aln

# ...and back: a hapmap into VCF (BGZF with -Oz; no -o streams to stdout)
reed convert panels.hmp.txt --vcf -Oz -o panels.vcf.gz

# The projected SFS
reed convert clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out
```

`reed <command> --help` lists every option.

## reed popg

The windowed population-genetics statistics. Pass any combination of
statistics flags; each one adds its output file.

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

`popgen` is a visible alias. The windowed flags were on `reed stats` in 0.3.x;
there they are now parse-time errors naming the rewrite.

## reed stats

The per-site and per-sample QC reports. Each flag adds its output file.

| Flag | Output |
|---|---|
| `--freq`, `--hardy`, `--het`, `--missing`, `--depth` | `PREFIX.afreq`, `PREFIX.hardy`, `PREFIX.het`, `PREFIX.smiss` + `PREFIX.vmiss`, `PREFIX.idepth` + `PREFIX.ldepth.mean` |
| `--summary` | QC summary in bcftools-stats layout, `PREFIX.stats` |
| `--hist <LIST>` | One `PREFIX.<m>.hist` TSV per metric |
| `--singletons` | Per-allele singleton/doubleton table, `PREFIX.singletons` (vcftools layout) |

`--hist` takes comma-separated built-ins (`qual`, `dp`, `dp-site`, `gq`,
`idl`, `af`) and any numeric tag as `INFO/<TAG>` or `FMT/<TAG>`. Every file
has `bin_lower`, `bin_upper`, `count`, `percent` columns; integer metrics also
carry `#mean`, `#median` and `#p95` lines. A tag that the header does not
declare, or one that is not Integer/Float, is a usage error.

`--pop FILE` is optional here (stratified QC); it is required by `popg`'s
windowed modes. With one of `--freq/--missing/--singletons` it stratifies:
`PREFIX.afreq.pop` and `PREFIX.vmiss.pop` carry one row per population
(first-appearance order), and `PREFIX.singletons` gains a POP column (`.` for
a carrier outside every group). The cohort files are unchanged, and `--pop`
without one of those three flags is a usage error. Multiallelic ALT
frequencies stay collapsed, as in `.afreq`.

## reed rename

Identifier edits: sample names, contig names, variant IDs and header lines.
Genotypes are never decoded; records pass through unchanged.

| Flag | Meaning |
|---|---|
| `-s "old new,…"`, `-N FILE` | Rename samples (`old<TAB>new` in the file) |
| `--replace-sample-name PATTERN` | Literal rewrite (`pre*` strips a prefix, `*suf` a suffix, else removes the text; regex metacharacters are an error) |
| `--rename-chrs FILE` | `old<TAB>new` contig names; the `##contig` lines and the CHROM column change together |
| `--chr-style plain\|ucsc\|mt-is-m` | Normalize contig names (no `chr`, `chr`, `chrM`) |
| `--set-id PATTERN` | Rewrite every ID from `%CHROM %POS %REF %ALT` |
| `--set-missing-var-ids PATTERN` | Fill in IDs that are `.` only |
| `--update-name FILE`, `--idspace-to SEP`, `--rm-dup-ids`, `--clear-ids` | ID mapping, separator swap, unique suffixes, all-dot |
| `-h FILE`, `-f FILE`, `--add-meta LINE`, `--remove-meta KEY` | Replace the header, set `##contig` lengths from a `.fai`, add or drop header lines |

Output: `-o FILE` (stdout by default), `-O v|z[level]`. `reheader` is a
visible alias. A rename onto an existing name, an unknown name and an empty
request are usage errors.

## reed filter

Selection, from bcftools: `-r/-R`, `-t/-T`, `-s/-S`, `--keep`, `--remove`
(with `--force-samples` to downgrade unknown names to warnings), `-v/-V`,
`--snps-only`, `-m/-M`, `-f`, `-i/-e`.

Site filters, applied after masking: `-q/-Q/--min-af/--max-af` with an
optional `:nref|alt1|minor|major|nonmajor` suffix, `--geno`, `--maf`,
`--max-maf`, `--mac`, `--max-mac`, `--hwe <P [midp]>`.

Genotype masking: `--min-dp`, `--max-dp`, `--min-gq` mask from FMT/DP and
FMT/GQ; `-i/-e` FORMAT expressions with `--set-gts .` mask the matches, as
`bcftools filter -S .` does. Masked genotypes count as missing for the site
filters that follow.

Sample filters: `--mind <FLOAT>` drops samples whose missing rate exceeds the
threshold. `--indep-pairwise <W[kb] [S] R2>` (filter only) prunes sites by LD
and emits the kept sites, aligned with plink2; the window `W` is a variant
count as a bare number and an interval with a `kb` suffix (the step is then
fixed at 1), `S` defaults to 1, and `R2` is the r² ceiling.

INFO handling: reed updates `INFO/AC` and `INFO/AN` when the sample set
changes (a record without a GT key keeps its input values); `-I/--no-update`
turns that off. `-x/--remove-info LIST` drops tags — by default
`AC,AN,AF,MLEAC,MLEAF,ExcessHet,InbreedingCoeff,DP`, so pass `-x ''` to keep
everything. `--fill-tags LIST` recomputes AC, AN, NS, MAF or HWE
(case-insensitive; an unsupported tag is an error, like bcftools).

Output: `-O v|z|b|u[0-9]` selects plain VCF, BGZF VCF, BCF or uncompressed
BCF, with an optional compression level. `-o` names the file (stdout by
default), `-W[=csi|tbi]` writes an index, `--report` writes a JSON run
report, `--timing` records per-phase timings, `--no-version` and
`--no-command` suppress the `##reed_*` provenance header lines.

## reed index

`-c/--csi` (default) or `-t/--tbi`, `-o` to write elsewhere, `-f` to
overwrite. The input must be a file.

## reed convert

Reads once and writes any number of targets at the same time. Supports
vcf2phylip, vcftools `--012`, TASSEL hapmap, MSA, and the reverse
hapmap → VCF.

| Flag | Output |
|---|---|
| `--phylip` / `--fasta` / `--nexus` | `PREFIX.phy` / `.fasta` / `.nexus` nucleotide matrices |
| `--nexus-binary` | `PREFIX.bin.nexus` (SNAPP) |
| `--012` | `PREFIX.012` + `.012.pos` + `.012.indv` |
| `--hapmap` | `PREFIX.hmp.txt` (TASSEL layout) |
| `--vcf [-O v\|z[level]]` | a hapmap input back to VCF |
| `--msa --reference ref.fa [--region CHR:B-E]` | `PREFIX_<contig>.fasta` |
| `--sfs --proj N[,N…] [--pop FILE]` | `PREFIX_sfs/dadi` and `PREFIX_sfs/fastsimcoal2` |
| `--outgroup NAME`, `--min-samples-locus N`, `--resolve-iupac`, `--write-used-sites` | controls, see `--help` |

The deliberate differences from the reference tools and how they are checked
are in `docs/convert-parity.md`.

## Input requirements

- VCF, VCF.gz or BCF, detected from the first bytes. `-` reads VCF text from
  stdin; `reed index` requires a file.
- The header must declare `##contig=<ID=…,length=…>`. Pass
  `--contig-lengths FILE` when the header lacks the lines.
- Records must be sorted by contig and position. Indexing and windowed
  statistics rely on it.
- `-r/-R` needs a `.csi` or `.tbi` next to the input, built with
  `reed index`. Plain scans do not.
- Genotypes must be diploid. reed rejects multiallelic records by default
  (`-M 2`), and treats half-calls such as `0/.` at the allele level, so
  `--geno` and `--mind` count them the way plink2 does.

## Output

`filter`, `convert` and `index` write the file named by `-o`. `stats` treats
`-o` as a prefix and writes one file per statistic, listed in the tables
above (`stats` has no stdout mode). `popg`'s window mode also needs `-o` — it
writes one file per statistic plus a `PREFIX.log`; its bare debug TSV (no
statistic requested) is the only output that goes to stdout.

Provenance header lines land in the output of `filter` (`##reed_version`,
`##reed_command`, `##reed_filterCommand`, `##reed_filter`) and in `popg`'s
`PREFIX.log`.

Windowed table format: one row per population (or population pair) per
window, with `chromosome`, `window_pos_1`, `window_pos_2`, the estimate, and
the counts behind it (`no_sites`, `count_diffs`, `count_comparisons`,
`count_missing` for pi and dxy; `no_snps` for Fst, PBS and the
f-statistics). Estimates with a zero denominator print `NA`.

## Library use

The `reed` crate is the public API; the writers produce the same files the
CLI writes.

```toml
[dependencies]
reed = { path = "crates/reed" }
```

```rust
use reed::{PopgenRequest, Populations, StatKind};

let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::File("pops.txt".into()),
    ..PopgenRequest::default()
};
let report = reed::run_popgen(Some("in.vcf.gz".as_ref()), &request, &reed::budget::ThreadBudget::default())?;
let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false);
```

A conversion run takes the same shape; the targets are typed and the report
lists every file written:

```rust
use reed::{ConvertRequest, ConvertTarget, PhyloFormats, PhyloOptions};

let request = ConvertRequest {
    targets: vec![ConvertTarget::Phylo {
        formats: PhyloFormats { phylip: true, fasta: true, nexus: false },
        options: PhyloOptions::default(),
    }],
    prefix: "matrices".into(),
    ..ConvertRequest::default()
};
let report = reed::run_convert(Some("in.vcf.gz".as_ref()), &request, &reed::budget::ThreadBudget::default())?;
```

A `sample -> population` mapping can be passed directly instead of writing a
populations file:

```rust
let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::Assignment(vec![
        ("sample1".into(), "popA".into()),
        ("sample2".into(), "popB".into()),
    ]),
    ..PopgenRequest::default()
};
```

The facade re-exports everything needed: `run_filter`, `run_popgen`,
`run_stats`, `run_rename`, `run_index`, `run_convert`, the request types
(`PopgenRequest`, `FilterSpec`, `StatsReportRequest`, `RenameSpec`,
`ConvertRequest`), the reports (`PopgenReport`, `FilterReport`,
`StatsOutcome`, `RenameReport`, `IndexOutcome`, `ConvertReport`), the writers
(`reed::windows`, `reed::reports`, `reed::histogram`, `reed::sfs`,
`reed::fstats`, `reed::ld`, `reed::summary`), and a single `reed::Error`.
`PopgenRequest::to_args(&group_names)` and `FilterSpec::to_args()` render a
request back to the CLI spelling; requests the CLI cannot spell (per-bp
kinds, a `Regions` grid) return a usage error. Only this facade carries a
stability promise; the lower crates under `crates/` are internal. Classify
errors with `Error::kind()`, read messages with `Display`/`source()`.

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