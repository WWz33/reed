---
title: filter
---

## Usage

```
reed filter [OPTIONS] <INPUT>
```

Hard filtering: site and sample selection, genotype masking, QC thresholds, LD pruning. Biallelic sites are kept by default (`-M 2`); output goes to stdout unless `-o` names a file. The command writes VCF only, no statistics.

## Options

### Input and regions

```
<INPUT>                     VCF/VCF.gz/BCF; `-` reads stdin. Multi-pass runs
                            (--mind, LD pruning) spill stdin to --tmp-dir first
-r, --regions LIST          Keep these regions, comma-separated: chr1, chr1:100,
                            chr1:100-200, chr1:100-, chr1:-200. Needs a
                            .csi/.tbi index next to the input
-R, --regions-file FILE     Region list, one per line. Three or more columns are
                            read as BED (0-based half-open), two columns as
                            CHROM POS (1-based), one column as a region string.
                            Needs an index
-t, --targets LIST          Same region syntax, streamed without an index
-T, --targets-file FILE     Region list file, same format as -R; a leading ^
                            excludes the region. No index needed
```

`#` comment lines and blank lines are skipped.

### Sample selection

```
-s, --samples [^]LIST       Keep these samples, comma-separated; a ^ prefix
                            complements (^A,B = everything except A and B)
-S, --samples-file FILE     Sample list in a file, one name per line
    --keep FILE             Same as -S, alternative spelling
    --remove FILE           Samples to drop, one name per line
    --force-samples         Names absent from the header only warn, no error
    --mind FLOAT            Drop samples whose missing rate exceeds FLOAT,
                            counted as missing alleles / 2·sites over the sites
                            that pass the pre-mask rules and after masking;
                            runs before --geno/--maf/--hwe; boundary kept
```

`-S` with `--keep`, `-s` with `-S/--keep`, and `-s ^…` with `--remove` are usage errors.

### Site selection

```
-v, --types LIST            Keep these variant types: snps, mnps, indels, other, ref
-V, --exclude-types LIST    Drop these variant types
    --snps-only[=just-acgt] Keep SNP sites. Bare means any single-base
                            substitution; =just-acgt also requires every
                            REF/ALT to be in A/C/G/T
-m, --min-alleles INT       Minimum allele count, REF+ALT
-M, --max-alleles INT       Maximum allele count, REF+ALT. Default 2; -M 1
                            keeps only invariant sites
-f, --apply-filters LIST    Keep records whose FILTER is in LIST, comma-separated
-i, --include EXPR          Keep records that satisfy the expression, repeatable
                            (see the filter expressions page)
-e, --exclude EXPR          Drop records that satisfy the expression, repeatable
-g, --genotype [^](hom|het|miss)
                            Select sites by genotype composition: hom a
                            homozygous-alt call, het a heterozygous call, miss a
                            missing call; repeatable, combined with AND; ^
                            negates. A half-call counts as missing
```

### Genotype masking

```
    --min-dp INT            Genotypes with FMT/DP below INT count as uncalled;
                            a missing DP never masks. Alias --vcf-min-dp
    --max-dp INT            Genotypes with FMT/DP above INT count as uncalled.
                            Alias --vcf-max-dp
    --min-gq INT            Genotypes with FMT/GQ below INT count as uncalled;
                            a missing GQ never masks. Alias --vcf-min-gq
    --set-gts .             With -i/-e: the FORMAT expression is evaluated per
                            sample, matching genotypes become missing and the
                            site stays. Only `.` is accepted
```

Masked genotypes count as missing in every site threshold below.

### Frequency and missingness thresholds

```
    --geno FLOAT            Drop sites whose missing-allele fraction exceeds
                            FLOAT, counted after masking. Per allele: a
                            half-call counts as 1 called allele, a fully missing
                            genotype as 2 missing; boundary kept
    --maf FLOAT             Drop sites whose minor allele frequency is below
                            FLOAT (after masking)
    --max-maf FLOAT         Upper bound on the minor allele frequency
    --mac INT               Lower bound on the minor allele count
    --max-mac INT           Upper bound on the minor allele count. Do not write
                            --max-mac 0: 0 is the internal "off" value and is
                            rejected; for MAC ≤ 1 write 1
-q, --min-af FLOAT[:type]   Lower bound on allele frequency. Reads INFO/AC and
                            INFO/AN when the sample set is unchanged, otherwise
                            recomputes from genotypes. :type picks the allele:
                            nref non-reference (default), alt1 first ALT, minor,
                            major, nonmajor
-Q, --max-af FLOAT[:type]   Upper bound on allele frequency, same :type values
    --hwe P [midp]          Drop sites whose HWE exact-test p is below P;
                            optional midp uses the mid-p value
```

`-q/-Q` prefer AC/AN from INFO; `--maf` always recomputes from genotypes. When a sample subset or masking changes the sample set, `-q/-Q` switch to genotypes and report stale INFO.

### LD pruning

```
    --indep-pairwise W[kb] [S] R2
                            Window, step, r² ceiling. The window counts variant
                            sites (write 50kb for base-pair distance, in which
                            case the step must be 1). For pairs inside the window
                            with r² above R2, the site with the lower minor
                            allele frequency is dropped. The output is the
                            pruned VCF
```

### Output

```
-O, --output-type v|z|b|u[0-9]
                            v plain VCF, z BGZF VCF, b BCF, u uncompressed BCF;
                            an optional compression level (-Oz6, -Ob1)
-o, --output FILE           Output file; without it, stdout. Pointing at the
                            input itself is rejected. Without -O the extension
                            decides (.gz compressed, .bcf for BCF). Alias --out
-W, --write-index[=csi|tbi] Build the index while writing: -W or -W=csi
                            (default), -W=tbi; omit for no index. BGZF output
                            only, and requires -o
    --report FILE           JSON run report: per-rule cumulative passes, mask
                            attribution, --mind drops, LD pruning counts, time
    --timing FILE           Milliseconds per stage
    --no-version            No ##reed_* provenance header lines
    --no-command            Omit only the ##reed_command line (the raw argv with
                            machine paths); other provenance lines stay
    --tmp-dir DIR           Spill directory for multi-pass stdin runs
    --threads INT           Worker threads, default all available cores
    --memory MB             In-flight memory budget in MiB, default half of
                            physical memory
    --verbosity INT         0 silent, 1 summary (default), 2 stage timing, 3 debug
```

### INFO tags

```
-I, --no-update             Leave INFO/AC and INFO/AN alone. By default they
                            are recomputed whenever the sample set changes
-x, --remove-info LIST      Drop these INFO tags. Default: AC,AN,AF,MLEAC,MLEAF,
                            ExcessHet,InbreedingCoeff,DP; -x '' drops none
    --fill-tags LIST        Recompute and write these tags from the current
                            sample set: AC, AN, NS, MAF, HWE, case-insensitive;
                            unsupported tags are usage errors
```

## Evaluation order

Every record passes eight stages in a fixed order. The order decides the funnel counts, not the final result; the one rule that changes the result is `--mind` running after masking:

1. Sample subset and region/position selection
2. Load-time site filters: `-v/-V/--snps-only/-m/-M/-f`
3. Genotype-independent `-i/-e`
4. Genotype masking: `--min-dp/--max-dp/--min-gq`, `--set-gts .` with FORMAT expressions
5. `--mind`
6. `--geno`, `--hwe`, `--maf/--max-maf/--mac/--max-mac`, genotype-dependent `-i/-e`
7. LD pruning (`--indep-pairwise`, needs the global site order)
8. INFO updates (`-x` → automatic AC/AN update → `--fill-tags`) and output

Rules that need a first pass over the whole input (`--mind`, LD pruning) scan once, then emit on a second pass.

Removed spellings have no silent alias: `--max-missing 0.9` is rejected with an error that asks for `--geno 0.1` (opposite meaning: 0.9 is the fraction that must be called).

## Examples

```bash
# region query, needs a reed index first
reed filter raw.vcf.gz -r chr1,chr2:1000-200000 -R targets.bed -o region.vcf.gz -Ov

# variant types and expression filters combined
reed filter raw.vcf.gz --snps-only=just-acgt -M 2 -f PASS -i 'QUAL>30' -e 'F_MISSING>0.1' -o snp.vcf.gz -Oz

# sample subset, dropping samples by missing rate on the way
reed filter raw.vcf.gz -s sampleA,sampleB --mind 0.2 -o subset.vcf.gz -Oz

# genotype masking: low depth, high depth and low GQ all count as uncalled
reed filter raw.vcf.gz --min-dp 5 --min-gq 20 --max-dp 100 -o masked.vcf.gz -Oz

# mask by expression, --set-gts turns matching genotypes into missing calls
reed filter raw.vcf.gz -e 'FMT/DP<5 || FMT/GQ<20' --set-gts . -o masked.vcf.gz -Oz

# site thresholds
reed filter clean.vcf.gz --geno 0.1 --maf 0.05 --hwe 1e-6 -o qc.vcf.gz -Oz

# LD pruning; the output is the pruned VCF
reed filter clean.vcf.gz --indep-pairwise 50 5 0.2 -Oz6 -o pruned.vcf.gz

# recompute INFO tags, dropping none of the old ones
reed filter raw.vcf.gz --fill-tags AC,AN,NS,MAF,HWE -x '' -o annotated.vcf.gz -Oz

# JSON report, per-stage timing and an index
reed filter raw.vcf.gz -o clean.vcf.gz -Oz6 -W --report clean.json --timing clean.timing
```

On completion stderr carries one summary line plus the filter funnel:

```
reed filter: 500 -> 131 records (10 of 10 samples), strategy single_pass, 0.0s
funnel:
  site.biallelic_gate          500
  after_mask.missing           177
  after_mask.maf               131
  after_mask.hwe               131
```

Each `funnel` line is the number of sites left after one rule. Unimplemented parameters such as `--mqf` or `--min-mean-dp` fail at parse time, not halfway through a run.
