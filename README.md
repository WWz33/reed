# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed computes population-genetics statistics on VCF/BCF files and applies
hard filters with byte-faithful output. It is pure Rust, uses every core by
default, and needs no system libraries: the BGZF and BCF codecs are compiled
in, so the binary runs anywhere the toolchain does. It is also a Rust library
(`reed`), so other tools can drive the same statistics without shelling out.

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
reed filter in.vcf.gz --indep-pairwise 50kb 5 0.2 -Oz6 -o pruned.vcf.gz

# CSI index
reed index in.vcf.gz

# Phylogenetic matrices, the 0/1/2 matrix, hapmap and a reference-driven MSA
reed convert in.vcf.gz --phylip --fasta --nexus --nexus-binary --012 --hapmap -o matrices
reed convert in.vcf.gz --msa --reference ref.fa --region chr1:1-100000 -o aln

# ...and back: a hapmap into VCF (BGZF with -Oz; no -o streams to stdout)
reed convert panels.hmp.txt --vcf -Oz -o panels.vcf.gz

# The projected SFS, byte-identical to `reed popg --sfs`
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

`--hist` takes comma-separated built-ins (`qual`, `dp`, `dp-site`, `gq`,
`idl`, `af`) and any numeric tag as `INFO/<TAG>` or `FMT/<TAG>`. Every file
has `bin_lower`, `bin_upper`, `count`, `percent` columns; integer metrics also
carry exact `#mean`, `#median` and `#p95` lines. A tag that the header does not
declare, or one that is not Integer/Float, is a usage error.

`--pop FILE` is optional here (stratified QC); it is required by `popg`'s
windowed modes. With one of `--freq/--missing/--singletons` it stratifies:
`PREFIX.afreq.pop` and `PREFIX.vmiss.pop` carry one row per population
(first-appearance order; the values equal a per-population plink2 run), and
`PREFIX.singletons` gains a POP column (`.` for a carrier outside every
group). The cohort files are unchanged, and `--pop` without one of those
three flags is a usage error. Multiallelic ALT frequencies stay collapsed,
as in `.afreq`.

## reed rename

Identifier edits: sample names, contig names, variant IDs and header lines.
Genotypes are never decoded, so the records pass through byte-for-byte.

| Flag | Meaning |
|---|---|
| `-s "old new,…"`, `-N FILE` | Rename samples (`old<TAB>new` in the file) |
| `--replace-sample-name PATTERN` | Rewrite matching names (`pre*` strips a prefix, `*suf` a suffix, else removes the text) |
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
FMT/GQ; `-i/-e` FORMAT expressions with `--set-gts .` mask exactly the
matches, as `bcftools filter -S .` does. Masked genotypes count as missing
for the site filters that follow.

Sample filters: `--mind <FLOAT>` drops samples whose missing rate exceeds the
threshold. `--indep-pairwise <W[kb] [S] R2>` prunes sites by LD using the
plink2 algorithm and emits the kept sites; the window `W` may carry a `kb`
suffix, the step `S` defaults to 1, and `R2` is the r² ceiling.

INFO handling: reed updates `INFO/AC` and `INFO/AN` when the sample set
changes (a record without a GT key keeps its input values); `-I/--no-update`
turns that off. `-x/--remove-info LIST` drops tags — by default the gold
pipeline's stale list (`AC,AN,AF,MLEAC,MLEAF,ExcessHet,InbreedingCoeff,DP`),
so pass `-x ''` to keep everything. `--fill-tags LIST` recomputes AC, AN, NS,
MAF or HWE (case-insensitive; an unsupported tag is an error, like bcftools).

Output: `-O v|z|b|u[0-9]` selects plain VCF, BGZF VCF, BCF or uncompressed
BCF, with an optional compression level. `-o` names the file (stdout by
default), `-W[=csi|tbi]` writes an index, `--report` writes a JSON run
report, `--timing` records per-phase timings, `--no-version` and
`--no-command` suppress the `##reed_*` provenance header lines.

## reed convert

Converts the input into phylogenetic and population-genetics matrices in one
parallel pass: the requested targets share one worker pool, and the output is
independent of `--threads`. The semantics follow the reference tools exactly
(vcf2phylip 2.9 by Ortiz 2019, vcftools `--012`, TASSEL hapmap, vcf2msa);
`gold/popgen/run_convert_parity.sh` checks them byte for byte.

| Flag | Output |
|---|---|
| `--phylip` / `--fasta` / `--nexus` | `PREFIX.phy` / `.fasta` / `.nexus` nucleotide matrices (relaxed PHYLIP layout, IUPAC ambiguity codes; `--phylip` is opt-in, unlike vcf2phylip) |
| `--nexus-binary` | `PREFIX.bin.nexus` for SNAPP (diploid biallelic sites; other genotypes are `?`) |
| `--012` | `PREFIX.012` (0/1/2/-1 per sample), `.012.pos`, `.012.indv` |
| `--hapmap` | `PREFIX.hmp.txt` (TASSEL layout) |
| `--vcf [-O v\|z[level]]` | **Reverse direction**: a hapmap (`.hmp.txt`) input back to VCF — `ID` from `rs#`, REF/ALT in declared order, IUPAC cells expanded to `GT` (`R` → `0/1`); `-O z` writes BGZF and omitting `-o` streams to stdout |
| `--msa --reference ref.fa [--region CHR:B-E]` | `PREFIX_<contig>.fasta`: whole-site alignment against the reference |
| `--sfs --proj N[,N…] [--pop FILE]` | `PREFIX_sfs/dadi` and `PREFIX_sfs/fastsimcoal2`, identical to `reed popg --sfs` |
| `--outgroup NAME`, `--min-samples-locus N`, `--resolve-iupac`, `--write-used-sites` | Alignment controls (see `--help`) |

`--hapmap` and `--vcf` are the two directions of one pair: on its lossless
domain (biallelic diploid SNPs without missing calls) a VCF → hapmap → VCF
round trip reproduces every record, and TASSEL 5 reads reed's VCF back.

Deliberate differences from the reference tools (registered in
`docs/spec/cli.md` §8b): malformed records are skipped and counted instead of
crashing or flooding stdout; `--resolve-iupac` is deterministic; `--012`
writes `-1` for half calls and haploids (vcftools writes `1`); the hapmap
`alleles` column lists the declared REF/ALT alleles; the MSA pads indels
deterministically instead of calling an external aligner; the hapmap→VCF
header carries no `##contig` (TASSEL writes none either) and its `##reed_*`
provenance lines can be turned off.

On an 8 GB VCF (348 samples, 3.3M records) the five phylogenetic outputs
come out byte-identical to vcf2phylip 2.9 while taking **10 s** against its
27 min (28 cores), because the per-sample temporary-file transpose is
replaced by in-memory column blocks.

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
use reed::{PopgenRequest, Populations, StatKind};

let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::File("pops.txt".into()),
    ..PopgenRequest::default()
};
let report = reed::run_popgen(Some("in.vcf.gz".as_ref()), &request, &reed::budget::ThreadBudget::default())?;
let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false);
```

A caller that already holds a `sample -> population` mapping can pass it
directly instead of writing a populations file:

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
for file in &report.files {
    println!("{}", file.display());
}
```

The facade re-exports everything needed: `run_filter`, `run_popgen`,
`run_stats`, `run_rename`, `run_index`, `run_convert`, the request types
(`PopgenRequest`, `FilterSpec`, `StatsReportRequest`, `RenameSpec`,
`ConvertRequest`), the reports (`PopgenReport`, `FilterReport`,
`StatsOutcome`, `RenameReport`, `IndexOutcome`, `ConvertReport`), the writers
(`reed::windows`, `reed::reports`, `reed::histogram`, `reed::sfs`,
`reed::fstats`, `reed::ld`, `reed::summary`),
and a single `reed::Error` whose `kind()` separates usage errors from
runtime failures. `PopgenRequest::to_args(&group_names)` and
`FilterSpec::to_args()` render a request back to the CLI spelling, so a tool
that builds requests programmatically can also log the equivalent command
line; `PopgenRequest::to_args` returns a usage error for the requests the CLI
cannot spell (per-bp kinds, a `Regions` grid). Only this facade carries a
stability promise; the lower crates under `crates/` are internal, and
`reed::Error` keeps them out of its signature by sealing its variants —
classify with `kind()`, read with `Display`/`source()`.

## Limitations

- Multiallelic records are rejected. reed does not recompute per-allele
  statistics when `-M` is raised past 2, so the default keeps them out.
- LD pruning lives in `filter` and needs its second pass. `reed stats`
  rejects `--indep-pairwise` and points at `reed filter`.
- f-statistics standard errors use block jackknife over base pairs, while
  scikit-allel uses variant counts. The estimates match; the errors are not
  comparable to allel's.
- Tajima's D divides by the variance at the window-average allele count, the
  pixy convention, which differs from vcftools' fixed 2N.

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
