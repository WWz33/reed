# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed is a population-genetics API library written in Rust. It computes π,
dxy, Fst, θ, Tajima's D and more on VCF/BCF files, and does filtering, LD
pruning, indexing and conversion.

## Installation

```bash
cargo install --path crates/reed-cli     # requires Rust 1.91+
reed --version
```

## Usage

```bash
# Windowed pi, dxy, Fst, theta and Tajima's D
reed popg allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out

# Variant sites only
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --fst --window 100k --step 50k -o out

# QC reports and summary
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# Histograms
reed stats clean.vcf.gz --hist qual,dp-site,gq,INFO/DP -o qc

# Rename samples, contigs and missing variant IDs
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

# A hapmap back to VCF (-Oz writes BGZF; without -o the output goes to stdout)
reed convert panels.hmp.txt --vcf -Oz -o panels.vcf.gz

# The projected SFS
reed convert clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out

# Principal components (plink-compatible output)
reed pca clean.vcf.gz --pcs 20 --missing pairwise -o pcs

# Individual-by-individual matrices
reed distance clean.vcf.gz --type grm -o rel      # .rel + .rel.id
reed distance clean.vcf.gz --type phi -o kin      # vcftools-compatible .relatedness2

# A tree from a distance matrix (ape/poppr-compatible newick)
reed distance clean.vcf.gz --type 1-ibs --tree nj -o tree    # tree.nwk
# --tree needs a distance matrix: the similarity types (grm, king, phi,
# ibs, pihat) are rejected because their diagonal is not zero.

# Bootstrap the tree: resample loci, write branch support as node labels
reed distance clean.vcf.gz --type 1-ibs --tree nj --bootstrap 100 -o tree

# Discriminant analysis of principal components (adegenet::dapc semantics)
reed pca clean.vcf.gz --dapc --pop pops.txt --pcs 20 -o dapc

# Cross-validate it: picks n_pca, writes PREFIX.dapc.xval
reed pca clean.vcf.gz --dapc --pop pops.txt --xval --xval-rep 100 -o dapc
```

## Commands

- `filter` .. selection, genotype masking, QC thresholds, LD pruning
- `stats` .. per-site and per-sample QC: freq, hardy, het, missing, depth, summary
- `popg` .. windowed statistics: pi, dxy, Fst, PBS, SFS, f-stats, LD decay
- `pca` .. principal components of the sample relationship matrix (`--dapc` adds a discriminant analysis; `--xval` cross-validates it)
- `distance` .. individual-by-individual matrices: distance, IBS, GRM, KING, relatedness, discordance; `--tree upgma|nj` builds a newick tree, `--bootstrap N` adds branch support
- `convert` .. phylogenetic matrices, genotype matrices, hapmap, MSA, SFS
- `rename` .. identifier edits: samples, contigs, variant IDs, header lines
- `index` .. build a CSI/TBI index

Run `reed` to list the commands; `reed <command> --help` lists every option.

## reed popg [OPTIONS] FILE

Population-genetics statistics per window, contig or genome.

### Statistics options

| Option | Description |
|---|---|
| `--pi` | π per population (`PREFIX_pi.txt`) |
| `--dxy` | dxy per population pair (`PREFIX_dxy.txt`) |
| `--fst [wc84\|hudson]` | Fst for all population pairs, wc84 by default (`PREFIX_fst.txt`) |
| `--fst-global` | Global WC84 Fst over all populations, with `--fst` (`PREFIX_fst_global.txt`) |
| `--theta` | Watterson's θ per population (`PREFIX_watterson_theta.txt`) |
| `--tajima-d` | Tajima's D per population (`PREFIX_tajima_d.txt`) |
| `--pbs A,B,C` | PBS for the named triple, repeatable; `--pbs-normed` writes PBSn1 (`PREFIX_pbs.txt`) |
| `--sfs [folded]` | Projected SFS, unfolded by default; needs `--proj N` (`PREFIX_sfs/dadi/`, `PREFIX_sfs/fastsimcoal2/`) |
| `--f2 A,B`, `--f3 C;A,B`, `--f4 A,B;C,D`, `--dstat A,B;C,D`, `--f4-ratio A,B;C,D;X` | Patterson f-statistics, block jackknife with 5 Mb blocks (`PREFIX.fstats`) |
| `--ld-decay [MAX_BP]` | Pairwise r² in 1 kb bins within MAX_BP, default 100k (`PREFIX.ld.decay`); `--ld-pairs` adds the per-pair table (`PREFIX.ld`) |

### Window options

| Option | Description |
|---|---|
| `--window SIZE` | Window size in bp; accepts `100k`, `1M`. Default: 100000 |
| `--step STEP` | Window step. Default: the window size |
| `--levels LIST` | Rows to report: `window`, `contig`, `genome`. Default: `window` |
| `--skip-empty` | Drop windows that no site fell into |
| `--pop FILE` | `sample<TAB>population` rows, no header; required by the windowed modes |
| `--variant-only` | Declare an all-variant input: window positions the file does not contain count as monomorphic with all samples called |
| `--bypass-invariant-check` | Compute over the sites the file contains (overestimates π/dxy) |
| `--contig-lengths FILE` | `.fai`-style `name<TAB>length`, for headers without `##contig` lines |

## reed pca [OPTIONS] FILE

Principal components of the sample relationship matrix (the GRM), matching
the plink2/GCTA/flashpca standardization.

| Option | Description |
|---|---|
| `--pcs INT` | Number of principal components to write. Default: 10 |
| `--stand binom2\|binom` | Standardization: `binom2` = `sqrt(2p(1-p))` (GCTA/plink2/flashpca, default), `binom` = `sqrt(p(1-p))` (older EIGENSOFT). The two give identical components with eigenvalues exactly 2x apart |
| `--missing pairwise\|meanimpute` | Missing-genotype policy: `pairwise` (per-pair denominator, plink2/GCTA default) or `meanimpute` (flashpca, plink2 `meanimpute`) |
| `--solver auto\|dense\|stream` | `auto` (default) picks `dense` when the N×N matrix fits the memory budget, `stream` otherwise. `dense` extracts the full spectrum, `stream` is matrix-free block Lanczos |
| `-o PREFIX` | Required: writes `PREFIX.eigenvec`, `PREFIX.eigenval`, `PREFIX.pve` |
| `-s/-S`, `--keep`, `--remove` | Sample subsets (same spelling as `reed filter`) |
| `--threads N`, `--memory MB` | Threads and memory budget |

`PREFIX.eigenvec` uses plink's header (`#FID IID PC1..PCk`) so existing
scripts read it unchanged. Only biallelic sites with a diploid GT contribute;
LD pruning is an upstream step, so run `reed filter --indep-pairwise` first
rather than passing it here.

## reed stats [OPTIONS] FILE

Per-site and per-sample QC reports.

| Option | Description |
|---|---|
| `--freq` | Per-site allele frequencies (`PREFIX.afreq`) |
| `--hardy` | HWE exact test, observed and expected heterozygosity (`PREFIX.hardy`) |
| `--het` | Per-sample inbreeding coefficient F (`PREFIX.het`) |
| `--missing` | Per-sample and per-site missingness (`PREFIX.smiss`, `PREFIX.vmiss`) |
| `--depth` | Per-sample and per-site mean depth from FMT/DP (`PREFIX.idepth`, `PREFIX.ldepth.mean`) |
| `--summary` | QC summary (`PREFIX.stats`) |
| `--hist LIST` | One `PREFIX.<m>.hist` histogram per metric: built-ins `qual`, `dp`, `dp-site`, `gq`, `idl`, `af`, or any numeric `INFO/<TAG>`/`FMT/<TAG>`; integer metrics carry `#mean`/`#median`/`#p95` |
| `--singletons` | Per-allele singletons and private doubletons (`PREFIX.singletons`) |
| `--pop FILE` | `sample<TAB>population` rows, no header; with `--freq`/`--missing` writes the stratified `PREFIX.afreq.pop` and `PREFIX.vmiss.pop` (one row per population), with `--singletons` adds a POP column |

`PREFIX.dp.hist`:

```
bin_lower	bin_upper	count	percent
5	5	12	1.234568
6	6	48	4.938272
#mean	8.104938
#median	7
#p95	21
```

## reed rename [OPTIONS] FILE

Edit sample names, contig names, variant IDs and header lines.

| Option | Description |
|---|---|
| `-s "old new,…"`, `-N FILE` | Rename samples (`old<TAB>new` rows in the file) |
| `--replace-sample-name PATTERN` | Literal rewrite (`pre*` strips a prefix, `*suf` a suffix, else removes the text) |
| `--rename-chrs FILE` | `old<TAB>new` contig names; the `##contig` lines and the CHROM column change together |
| `--chr-style plain\|ucsc\|mt-is-m` | Normalize contig names (no `chr`, `chr`, `chrM`) |
| `--set-id PATTERN` | Rewrite every ID from `%CHROM %POS %REF %ALT` |
| `--set-missing-var-ids PATTERN` | Fill in IDs that are `.` only |
| `--update-name FILE`, `--idspace-to SEP`, `--rm-dup-ids`, `--clear-ids` | ID mapping, separator swap, unique suffixes, all-dot |
| `-h FILE`, `-f FILE`, `--add-meta LINE`, `--remove-meta KEY` | Replace the header, set `##contig` lengths from a `.fai`, add or drop header lines |

Output: `-o FILE` (stdout by default), `-O v|z[level]`.

## reed filter [OPTIONS] FILE

Subset, filter and mask VCF/BCF; prune by LD; write VCF or BCF.

| Option | Description |
|---|---|
| `-r LIST` / `-R FILE` | Regions; needs a `.csi`/`.tbi` next to the input, built with `reed index`. Files may be a plain region list or a **BED** (0-based half-open) |
| `-t LIST` / `-T FILE` | Target positions; streamed, no index needed. Same region spellings; `CHROM:POS` selects one site, files may be BED or `CHROM<TAB>POS` |
| `-s/-S`, `--keep`, `--remove` | Sample subsets; `--force-samples` downgrades unknown names to warnings |
| `-v/-V`, `--snps-only` | Select by variant type; `--snps-only` keeps SNPs only |
| `-m/-M INT` | Min/max alleles per record; default `-M 2`, multiallelic records rejected |
| `-f LIST` | Keep sites whose FILTER column contains a value from LIST |
| `-i/-e EXPR` | Include/exclude expressions |
| `-q/-Q/--min-af/--max-af FLOAT[:type]`, `--geno`, `--maf`, `--max-maf`, `--mac`, `--max-mac`, `--hwe <P>` (the mid-p variant is not implemented in this version) | Site filters, applied after masking: allele-frequency bounds (`:type` is `nref\|alt1\|minor\|major\|nonmajor`), missing rate, MAC (`--max-mac 0` is refused — `0` is the internal "off" value), HWE |
| `--min-dp`, `--max-dp`, `--min-gq` | Mask genotypes from FMT/DP and FMT/GQ; masked genotypes count as missing for the site filters |
| `--set-gts .` | With `-i/-e`, set matching FORMAT genotypes to `.` |
| `--mind FLOAT` | Drop samples whose missing rate exceeds FLOAT |
| `--indep-pairwise W[kb] [S] R2` | LD pruning, writes the kept sites: `W` is a variant count as a bare number, an interval with a `kb` suffix (step fixed at 1); `S` defaults to 1; `R2` is the r² ceiling and must be a finite number in 0–1 |
| `-I/--no-update` | Do not update `INFO/AC` and `INFO/AN` (updated with the sample set by default; records without a GT keep their input values) |
| `-x/--remove-info LIST` | Drop INFO tags; default `AC,AN,AF,MLEAC,MLEAF,ExcessHet,InbreedingCoeff,DP`, `-x ''` keeps every tag |
| `--fill-tags LIST` | Recompute AC, AN, NS, MAF or HWE (case-insensitive) |
| `-O v\|z\|b\|u[0-9]`, `-o FILE`, `-W[=csi\|tbi]` | Output type (optional compression level), output file (stdout by default), write an index (`-W` needs a BGZF output and `-o FILE`; anything else is a usage error) |
| `--report`, `--timing`, `--no-version`, `--no-command` | JSON run report, per-phase timings, suppress the `##reed_*` lines |

## reed index [OPTIONS] FILE

| Option | Description |
|---|---|
| `-c/--csi` | Write a CSI index (default) |
| `-t/--tbi` | Write a TBI index |
| `-o FILE` | Write the index elsewhere |
| `-f` | Overwrite an existing index |

The input must be a file.

## reed convert [OPTIONS] FILE

Convert VCF/BCF to other formats; a hapmap input converts back to VCF.

| Option | Output |
|---|---|
| `--phylip` / `--fasta` / `--nexus` | `PREFIX.phy` / `.fasta` / `.nexus` nucleotide matrices |
| `--nexus-binary` | `PREFIX.bin.nexus` |
| `--012` | `PREFIX.012` + `.012.pos` + `.012.indv` |
| `--hapmap` | `PREFIX.hmp.txt` |
| `--vcf [-O v\|z[level]]` | a hapmap input back to VCF |
| `--msa --reference ref.fa [--region CHR:B-E]` | `PREFIX_<contig>.fasta` |
| `--sfs --proj N[,N…] [--pop FILE]` | `PREFIX_sfs/dadi`, `PREFIX_sfs/fastsimcoal2` |
| `--outgroup NAME`, `--min-samples-locus N`, `--resolve-iupac`, `--write-used-sites` | controls, see `--help` |

Differences are recorded in `docs/convert-parity.md`.

## Expressions (`-i/-e`)

`-i` keeps sites where the expression is true, `-e` drops them.

| Category | Values |
|---|---|
| Columns | `CHROM`, `POS`, `ID`, `QUAL`, `FILTER`, `TYPE` |
| Tags | `INFO/TAG`, `FMT/TAG`; a bare `TAG` looks in INFO first |
| Computed | `N_ALT`, `N_SAMPLES`, `N_MISSING`, `F_MISSING`; `AC`, `AN`, `AF`, `MAC`, `MAF` (the INFO value when present, else from the genotypes) |
| Operators | `==` `=` `!=` `<` `<=` `>` `>=`; regex `~` `!~`; `&&` `\|\|` `!` and parentheses |
| Missing | `"."`; comparisons against missing are false |

```
reed filter in.vcf.gz -i 'QUAL>30 && MAF>=0.05' -Oz -o out.vcf.gz
reed filter in.vcf.gz -e 'FMT/DP<5' --set-gts . -Oz -o out.vcf.gz

# GATK SNP hard filter
reed filter in.vcf.gz -e 'QD<2 || QUAL<30 || SOR>3 || FS>60 || MQ<40 || MQRankSum<-12.5 || ReadPosRankSum<-8' -Oz -o out.vcf.gz
# GATK indel hard filter
reed filter in.vcf.gz -e 'QD<2 || QUAL<30 || FS>200 || ReadPosRankSum<-20' -Oz -o out.vcf.gz
```

## Output

`filter`, `rename`, `convert` and `index` write the file named by `-o`;
`stats` and `popg` take `-o` as a prefix (`stats` has no stdout mode). `popg`
also writes a `PREFIX.log`.

`PREFIX_pi.txt`:

```
pop	chromosome	window_pos_1	window_pos_2	avg_pi	no_sites	count_diffs	count_comparisons	count_missing
popA	chr1	1	10000	0.001234	95	12	812	3
```

Estimates with a zero denominator print `NA`.

## Library use

The `reed` crate is the public API.

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

A `sample -> population` mapping can be passed directly, without a
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

The crate root re-exports: `run_filter`, `run_popgen`,
`run_stats`, `run_rename`, `run_index`, `run_convert`, the request types
(`PopgenRequest`, `FilterSpec`, `StatsReportRequest`, `RenameSpec`,
`ConvertRequest`), the reports (`PopgenReport`, `FilterReport`,
`StatsOutcome`, `RenameReport`, `IndexOutcome`, `ConvertReport`), the writers
(`reed::windows`, `reed::reports`, `reed::histogram`, `reed::sfs`,
`reed::fstats`, `reed::ld`, `reed::summary`), and a single `reed::Error`.
`PopgenRequest::to_args(&group_names)` and `FilterSpec::to_args()` render a
request back to the CLI spelling; requests the CLI cannot spell (per-bp
kinds, a `Regions` grid) return an error. Only the crate root carries a
stability promise; the lower crates under `crates/` are internal. Classify
errors with `Error::kind()`, read messages with `Display`/`source()`.

## Citation

reed implements published estimators. Cite the methods it uses:

- Nei & Li 1979 for π; Nei 1987 for dxy; Watterson 1975 for θ; Tajima
  1989 for D.
- Hudson 1992 and Bhatia et al. 2013 for Hudson Fst; Weir & Cockerham 1984
  for WC84 Fst.
- Patterson et al. 2012 for f-statistics; Yi et al. 2010 for PBS.
- Korunes & Samuk 2021 for the missing-data treatment of π, dxy, theta and
  D.

## License

MIT OR Apache-2.0. See [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).
