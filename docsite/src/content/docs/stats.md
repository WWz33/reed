---
title: stats
---

## Usage

```
reed stats [OPTIONS] --freq --hardy ... -o PREFIX <INPUT>
```

Per-site and per-sample QC reports. Single pass; `-o` is the output prefix and reports are written to files. Missing `-o`, no report selected, or `--pop` without a report it can stratify are usage errors.

## Options

### Reports

```
    --freq                  Per-site allele frequency → PREFIX.afreq
    --hardy                 HWE exact test and observed/expected heterozygosity
                            → PREFIX.hardy
    --het                   Per-sample inbreeding coefficient F → PREFIX.het
    --missing               Per-sample and per-site missingness →
                            PREFIX.smiss, PREFIX.vmiss
    --depth                 Mean depth from FMT/DP → PREFIX.idepth,
                            PREFIX.ldepth.mean
    --summary               QC summary → PREFIX.stats
    --hist LIST             Histograms, one file per metric →
                            PREFIX.<metric>.hist. Built-ins: qual, dp (FMT/DP),
                            dp-site (INFO/DP), gq, idl (signed indel length),
                            af; also accepts INFO/<TAG> and FMT/<TAG>
    --singletons            Per-allele singletons and private doubletons →
                            PREFIX.singletons
    --sites-file FILE       CHROM<TAB>POS (1-based); report only these sites
-o, --output PREFIX         Output prefix; reports go to files, not stdout.
                            Alias --out
```

### Population stratification

```
    --pop FILE              sample<TAB>population, two columns, no header. Adds
                            stratified files to --freq/--missing and a POP
                            column to --singletons. Alias --populations
    --force-samples         Unknown samples in the population file only warn
                            and are skipped; default is an error
```

### Input and regions

```
<INPUT>                     VCF/VCF.gz/BCF; `-` reads stdin
-r, --regions LIST          Keep these regions, comma-separated; needs an index
                            next to the input
-R, --regions-file FILE     Region list file, same format as filter; needs an index
-t, --targets LIST          Same region syntax, streamed, no index needed
-T, --targets-file FILE     Region list file; a leading ^ excludes. No index
```

### Global options

```
    --threads INT           Worker threads, default all available cores
    --memory MB             In-flight memory budget in MiB, default half of
                            physical memory
    --tmp-dir DIR           Accepted for compatibility; single-pass report runs
                            never spill
    --verbosity INT         0 silent, 1 summary (default), 2 stage timing, 3 debug
```

### Rejected options

```
-s, --samples / -S, --samples-file / --keep / --remove
                            Not accepted by report mode; the error asks for
                            `reed filter` first
    --mind                  Separate error: the sample rule needs a full pass,
                            the message asks for reed filter --mind X
    --indep-pairwise        Separate error: pruning needs a second pass, the
                            message asks for reed filter --indep-pairwise
-v, --types / -V, --exclude-types / --snps-only
-m, --min-alleles / -M, --max-alleles / -f, --apply-filters
-i, --include / -e, --exclude / -g, --genotype
    --min-dp / --max-dp / --min-gq / --set-gts
                            (masking aliases: --vcf-min-dp / --vcf-max-dp /
                            --vcf-min-gq)
-q, --min-af / -Q, --max-af / --geno / --maf / --max-maf / --mac /
--max-mac / --hwe
                            Site rules and masking are rejected with the same
                            report-mode error
```

The report-mode error:

```
reed: report mode (--freq/--hardy/--het/--missing/--depth/--summary) does not apply
-s/-S/--keep/--remove or site/threshold filters yet; run `reed filter` first and feed its output to reed stats
```

## Output format

`qc/raw.singletons` (the file carries a header, one allele per line; S means a single copy, D means two copies with a homozygous carrier):

```
CHROM	POS	SINGLETON/DOUBLETON	ALLELE	INDV
chr1	501	D	G	s3
```

`qc/raw.vmiss`:

```
#CHROM	ID	MISSING_CT	OBS_CT	F_MISS
chr1	rs2	1	10	0.1
chr1	rs4	1	10	0.1
chr1	rs8	1	10	0.1
```

`qc/raw.dp.hist` (buckets are the DP values themselves; empty buckets are omitted):

```
bin_lower	bin_upper	count	percent
2	2	126	2.520000
3	3	47	0.940000
4	4	58	1.160000
5	5	70	1.400000
```

Multiallelic sites are not split: ALT_FREQS in `.afreq` is the combined frequency of all ALT alleles (`CollapseAlt` semantics), consistent across the family.

## HWE exact test

`--hardy` runs the same test as `filter --hwe`: given the number of called diploids `n` and the minor allele count `mac`, the conditional distribution of heterozygotes `nh` is `P(nh) ∝ n!/(hm!·nh!·hM!)·2^nh` (`hm = (mac−nh)/2`), recursed in log space as `P(nh+2)/P(nh) = 4·hm·hM/((nh+1)(nh+2))`; after normalization the two-sided p is the sum of all probabilities not above the observed one. Monomorphic sites get p = 1. The mid-p variant (`p_mid = p_exact − 0.5·p_obs`) exists in the library but is not exposed as a command-line option.

## Histograms

Built-in metrics:

| Metric | Value |
|---|---|
| `qual` | QUAL |
| `dp` | `FMT/DP` (genotype level) |
| `dp-site` | `INFO/DP` (site level) |
| `gq` | `FMT/GQ` |
| `idl` | Signed indel length (deletions negative, insertions positive) |
| `af` | Allele frequency |

Any `INFO/<TAG>` or `FMT/<TAG>` declared as Integer or Float in the header works too. Bucket boundaries use the quantity's own unit: QUAL by quality value, `af` by frequency. Float metrics get float boundaries; integer metrics carry three exact quantile lines at the end:

```
#mean	31.012800
#median	30
#p95	58
```

A `#missing` line appears when observations are missing; integer observations beyond the accumulator range are written as a `>max` line and stay out of the quantiles. The tag must be declared in the header with a numeric type, otherwise it is a usage error.

## Population stratification

`--pop FILE` (`sample<TAB>population`, no header) adds stratified data to three reports:

- `--freq` also writes `PREFIX.afreq.pop` (columns `#CHROM ID POP REF ALT ALT_FREQS OBS_CT`), one row per population;
- `--missing` also writes `PREFIX.vmiss.pop` (columns `#CHROM ID POP MISSING_CT OBS_CT F_MISS`);
- `--singletons` gains a sixth POP column (the carrier's population; samples outside the population file get `.`, counts still cover all samples).

```
qc.afreq.pop
#CHROM	ID	POP	REF	ALT	ALT_FREQS	OBS_CT
chr1	rs2	popA	A	C	0.5	10
chr1	rs2	popB	A	C	0.5	8

qc.vmiss.pop
#CHROM	ID	POP	MISSING_CT	OBS_CT	F_MISS
chr1	rs2	popA	0	5	0
chr1	rs2	popB	1	5	0.2
```

A stratified value equals the corresponding report run on that population's samples alone. Unknown samples in the population file are an error by default; `--force-samples` warns and skips.

## Examples

```bash
mkdir -p qc

# full QC
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc/raw

# histograms for threshold selection
reed stats raw.vcf.gz --hist qual,dp-site,gq,INFO/DP -o qc/raw

# per-allele singletons
reed stats extra.vcf.gz --singletons -o qc/raw

# stratified allele frequency and missingness
reed stats clean.vcf.gz --pop pops.txt --freq --missing -o qc/by_pop

# depth over one region
reed stats raw.vcf.gz -r chr1:1-1000000 --depth -o qc/chr1

# only the sites on a list
reed stats raw.vcf.gz --sites-file candidate_sites.txt --freq -o qc/cand

# filter first, then report
reed filter raw.vcf.gz --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 -o clean.vcf.gz -Oz
reed stats clean.vcf.gz --freq --missing -o qc/clean
```
