# Quick start

One analysis line that runs the main reed commands end to end: index, QC, filter, LD pruning, population statistics, PCA, tree building, format conversion.

## Data

The input is a VCF or BCF. This page assumes:

- `raw.vcf.gz`: raw calls, invariant sites included;
- `pops.txt`: population file, `sample<TAB>population` in two columns, no header:

```
sampleA	popA
sampleB	popA
sampleC	popB
sampleD	popB
```

A sample name in the population file that is absent from the header is an error. Add `--force-samples` to warn and skip it instead.

## Step 1: index

```bash
reed index raw.vcf.gz
```

This writes `raw.vcf.gz.csi`, needed later for `filter -r/-R` region queries. Add `-t` for a TBI index.

## Step 2: QC, see what the data looks like

```bash
mkdir -p qc
reed stats raw.vcf.gz --freq --hardy --het --missing --depth --summary -o qc/raw
reed stats raw.vcf.gz --hist qual,dp-site,gq -o qc/raw
```

The first line reports per-site allele frequency, HWE exact test, per-sample inbreeding coefficient, missingness, mean depth and a summary. The second writes histograms of QUAL, site depth and GQ for choosing filter thresholds. Outputs use `qc/raw` as the prefix (`qc/raw.afreq`, `qc/raw.hardy`, `qc/raw.dp-site.hist`, ...).

A histogram file looks like this; integer metrics carry mean, median and 95th percentile at the end:

```
bin_lower	bin_upper	count	percent
2	2	126	2.520000
3	3	47	0.940000
4	4	58	1.160000
5	5	70	1.400000
…
#mean	31.012800
#median	30
#p95	58
```

## Step 3: hard filtering

Site thresholds and genotype masks go in one command. The evaluation order is: sample subset → variant types → expressions → genotype mask → sample missingness → site missingness, HWE and frequency thresholds. Genotypes are masked first and then take part in site filtering:

```bash
reed filter raw.vcf.gz -o clean.vcf.gz -Oz6 \
    --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 --hwe 1e-6 -W
```

On completion stderr carries one summary line and the filter funnel, the sites left after each rule:

```
reed filter: 500 -> 131 records (10 of 10 samples), strategy single_pass, 0.0s
funnel:
  site.biallelic_gate          500
  after_mask.missing           177
  after_mask.maf               131
  after_mask.hwe               131
```

- `--min-dp 5 --min-gq 20`: genotypes with DP below 5 or GQ below 20 become `./.`. Genotypes whose DP or GQ is itself missing are untouched.
- `--geno 0.1`: drop sites whose missing-allele fraction after masking exceeds 10%.
- `--maf 0.05`: drop sites whose minor allele frequency is below 5%. Invariant sites (MAF=0) go too; pass `--variant-only` when feeding the filtered file to `popg`.
- `--hwe 1e-6`: drop sites whose HWE exact-test p is below the threshold.
- `-Oz6`: BGZF compression, level 6.
- `-W`: build the index while writing (needs BGZF output and `-o`).

Without `-W` no index is written. `-O` also accepts `v` (plain VCF), `b` (BCF) and `u` (uncompressed BCF).

Expression filters are the other way, for rules such as "pick sites by an INFO value":

```bash
# drop sites with QUAL below 30 or MAF below 0.05
reed filter raw.vcf.gz -i 'QUAL>30 && MAF>=0.05' -Oz -o clean2.vcf.gz

# set shallow genotypes to missing, keep the site
reed filter raw.vcf.gz -e 'FMT/DP<5' --set-gts . -Oz -o masked.vcf.gz
```

Expression syntax is on the [Filter expressions](expressions.md) page.

## Step 4: LD pruning

Prune LD before PCA. `--indep-pairwise` writes the pruned VCF directly; one command selects the sites and extracts them:

```bash
reed filter clean.vcf.gz --indep-pairwise 50 5 0.2 -Oz6 -o pruned.vcf.gz
```

The three numbers are window, step and r² ceiling: a window of 50 variant sites, advancing 5 at a time; for site pairs with r² strictly above 0.2 inside the window, the one with the lower minor allele frequency is dropped. Write the window as `50kb` to count by base-pair distance instead (the step must then be 1).

## Step 5: windowed population statistics

```bash
reed popg clean.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out
```

Several statistics in one run, each in its own file (`out_pi.txt`, `out_dxy.txt`, `out_fst.txt`, `out_watterson_theta.txt`, `out_tajima_d.txt`).

The π and dxy denominators rest on invariant sites, and `popg` assumes the file has them (whole-genome calls); it checks at startup. Frequency filters such as `--maf` remove invariant sites too, so running on a filtered file fails with this error:

```
reed: job failed: input has no invariant sites (122 records, none with ALT = "."):
declare --variant-only for an all-variant file, or --bypass-invariant-check to
compute over present sites only (both overestimate π/dxy denominators)
```

Two typical cases:

- Whole-genome calls: run as is. Empty windows are still written with `NA`; add `--skip-empty` to drop them.
- A file reduced to variant sites by filtering (such as `clean.vcf.gz`): add `--variant-only`; positions the window does not contain count as monomorphic sites with all samples called.

To take neither route and compute only over the sites the file holds, use `--bypass-invariant-check`; reed warns that π and dxy are overestimated.

Rows of `out_pi.txt`:

```
pop	chromosome	window_pos_1	window_pos_2	avg_pi	no_sites	count_diffs	count_comparisons	count_missing
popA	chr1	1	10000	0.0003778491492837536	9	170	320	85
popB	chr1	1	10000	0.00040893523488217553	9	184	354	51
```

PBS, f-statistics and LD decay:

```bash
# PBS: A is the focal population; --pbs repeats for several triples
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --pbs popA,popB,popC --pbs popA,popB,popD --window 50k -o out

# f-statistics: block-jackknife standard errors, 5 Mb blocks by default, --block-size changes that
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --f4 popA,popB;popC,popD --dstat popA,popB;popC,popD -o out

# LD decay: r² binned at 1 kb, 100 kb span by default; --ld-pairs adds the per-pair table
reed popg clean.vcf.gz --pop pops.txt --variant-only --ld-decay 200k -o out

# projected SFS: dadi and fastsimcoal2 layouts; two populations get two projection sizes
reed popg clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out
```

Population names must match `pops.txt`; f-statistic population combinations are written `A,B`, `C;A,B`, `A,B;C,D` and `A,B;C,D;X`.

## Step 6: PCA and discriminant analysis

```bash
reed pca pruned.vcf.gz --pcs 20 -o pcs
```

This writes `pcs.eigenvec` (header `#FID IID PC1..PCk`), `pcs.eigenval` and `pcs.pve`:

```
#FID	IID	PC1	PC2	PC3
s1	s1	0.253837528	-0.199416559	0.047029802
s2	s2	0.498821822	-0.283700482	-0.197574393
```

Discriminant analysis with population labels, and cross-validation to choose the number of PCs:

```bash
reed pca pruned.vcf.gz --dapc --pop pops.txt --pcs 5 -o dapc

# cross-validation: 30 repetitions by default, 90% training, fixed seed
reed pca pruned.vcf.gz --dapc --pop pops.txt --xval -o dapc
```

Cross-validation prints the `n_pca` it selected and refits with it in place of `--pcs`:

```
reed pca --dapc --xval: best n_pca=5 (mean success 0.8167), highest mean success at n_pca=5, random chance 0.5000
  note: --pcs 10 was overridden by the cross-validated choice (5)
```

The candidate curve and summary go to `dapc.dapc.xval`.

## Step 7: between-individual matrices and trees

```bash
reed distance clean.vcf.gz --type grm -o rel        # rel.rel + rel.rel.id
reed distance clean.vcf.gz --type king -o kin       # kin.king + kin.king.id
reed distance clean.vcf.gz --type phi -o rela       # rela.relatedness2
```

Tree building from a distance matrix, with bootstrap support:

```bash
reed distance clean.vcf.gz --type 1-ibs --tree nj --bootstrap 100 -o tree
```

In `tree.nwk` the numbers at internal nodes are bootstrap support:

```
(((s1:0.2010,s4:0.2173)35:0.0187,((s2:0.2185,s6:0.2104)25:0.0134,
  (s5:0.1824,s9:0.1875)60:0.0347)20:0.0148)20:0.0087,
 ((s3:0.1950,s10:0.1764)40:0.0181,s8:0.2309)20:0.0090,s7:0.2161)100;
```

`--tree upgma|nj` accepts only distance matrices with a zero diagonal; similarity kinds such as `grm`, `king` and `ibs` are rejected with a hint to use `1-ibs`, `allele-ct`, `ajk` or `discordance`. The random seed is 42 by default, so runs are reproducible.

## Step 8: format conversion

```bash
# phylogenetic and genotype matrices, one read, several targets
reed convert clean.vcf.gz --phylip --fasta --nexus --nexus-binary --012 --hapmap -o matrices

# reference-driven MSA
reed convert clean.vcf.gz --msa --reference ref.fa --region chr1:1-100000 -o aln

# projected SFS (same output as popg --sfs)
reed convert clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out

# hapmap back to VCF
reed convert panels.hmp.txt --vcf -Oz -o panels.vcf.gz
```
