# pca

## Usage

```
reed pca [OPTIONS] -o PREFIX <INPUT>
```

Principal components of the sample relationship matrix (GRM). The first 10 components are written by default, with `binom2` standardization and `pairwise` missing-data handling. Only biallelic sites with diploid GT are counted: multiallelic sites, haploid and half-calls, and sites with no called ALT are skipped; with no usable site at all the command reports `no biallelic sites with a diploid GT were found`.

## Options

### Principal components

```
<INPUT>                     VCF/VCF.gz/BCF; `-` reads stdin
    --pcs INT               Number of components to write, default 10, minimum 1
    --stand binom2|binom    Standardization: binom2 uses √(2p(1−p)) (default),
                            binom uses √(p(1−p)). Both give the same components
                            with eigenvalues differing by 2x
    --missing pairwise|meanimpute
                            Missing-genotype policy: pairwise (default) uses
                            each sample pair's own jointly-called count as the
                            denominator; meanimpute fills missing dosages with
                            the site mean and divides by all sites
    --solver auto|dense|stream
                            Solver. auto (default) uses dense when the N×N
                            matrix fits the budget (full spectrum), otherwise
                            stream (matrix-free block Lanczos)
```

### DAPC and cross-validation

```
    --dapc                  Discriminant analysis on the components; needs
                            --pop; writes PREFIX.dapc.ind, PREFIX.dapc.eig and
                            PREFIX.dapc.var
    --pop FILE              Population file for DAPC, sample<TAB>population,
                            no header. Alias --populations
    --n-da INT              Discriminant axes to retain, default every
                            significant axis (svd² > 1e-10)
    --force-samples         Unknown samples in the population file only warn and
                            are skipped; default is an error
    --xval                  Cross-validate to choose the PC count; writes
                            PREFIX.dapc.xval. Needs --dapc
    --xval-rep INT          Repetitions per candidate, default 30
    --xval-training FLOAT   Training fraction, default 0.9; lowered
                            automatically when the smallest population cannot
                            leave a validation sample
    --xval-result overall|group-mean
                            Accuracy: overall (default) or group-mean, each
                            population weighted equally
    --xval-seed INT         Random seed, default 42
    --xval-pca LIST         Candidate PC counts, comma-separated; by default
                            about 10 candidates on a fixed scale, deduplicated
```

### Samples and output

```
-s, --samples [^]LIST       Sample subset, comma-separated, ^ prefix complements
-S, --samples-file FILE     Sample list file, one name per line
    --keep FILE             Same as -S, alternative spelling
    --remove FILE           Samples to drop
-o, --out PREFIX            Output prefix
    --threads INT           Worker threads, default all available cores
    --memory MB             Memory budget in MiB, default half of physical memory
    --verbosity INT         0 silent, 1 summary (default), 2 stage timing, 3 debug
```

`--missing pairwise` and `--solver stream` are incompatible; together they raise an error listing three ways out: use `dense`, use `meanimpute`, or raise `--memory`. On low memory the error suggests running `reed filter --indep-pairwise` first to cut sites.

## Algorithm

The data flow has three steps.

1. Standardization. Each site's allele frequency `p` is estimated with called dosages only; dosage `x ∈ {0,1,2}` is transformed to

```
--stand binom2 (default):  W = (x − 2p) / √(2p(1−p))
--stand binom:             W = (x − 2p) / √(p(1−p))
```

A missing dosage becomes 0 in standardized space (the site mean); monomorphic or fixed-ALT degenerate columns are zeroed whole and enter no numerator.

2. Relationship matrix.

```
--missing pairwise (default):  A_ij = Σ_{s: both i and j called} W_is·W_js / (jointly called sites)
--missing meanimpute:          A_ij = Σ_s W_is·W_js / m            (m = total sites)
```

`pairwise` counts only sites called in both samples, both in numerator and denominator, each pair with its own denominator; `meanimpute` has missing filled with the mean, so the denominator is all sites. The two denominators do not produce the same `A_ij`: `meanimpute` pulls relationships toward the mean, `pairwise` is noisier for pairs with low coverage. Coordinates are only comparable across runs that use the same missing-data policy.

3. Eigen-solve. `dense` materializes the N×N matrix and runs a dense symmetric eigendecomposition over the full spectrum; `stream` scans genotype blocks and iterates, with a fixed seed, so results are reproducible.

## Output

| File | Content |
|---|---|
| `PREFIX.eigenvec` | Header `#FID IID PC1..PCk`, tab-separated; FID and IID both hold the sample name; eigenvectors have unit norm, 9 decimals. The overall sign is arbitrary; compare runs by direction |
| `PREFIX.eigenval` | Eigenvalues, descending, one per line, 9 decimals (= the GRM eigenvalues) |
| `PREFIX.pve` | Each component's share of the full spectrum: λᵢ over the spectrum sum (`stream` uses `Σw²/m` as the sum, equivalent to the trace) |

```
PREFIX.eigenvec
#FID	IID	PC1	PC2	PC3
s1	s1	0.253837528	-0.199416559	0.047029802
s2	s2	0.498821822	-0.283700482	-0.197574393
s3	s3	-0.238191065	0.372053836	-0.082970657
…

PREFIX.eigenval
2.912527165
2.560773428
2.291419510

PREFIX.pve
0.218853435
0.192421917
0.172182096
```

## DAPC discriminant analysis

`--dapc` runs a discriminant analysis on top of the principal components: the first `n_pca` components are taken (default `min(samples−1, 50)`, components with eigenvalues below 1e-14 excluded from the rank), sample scores are `U·√λ`, and a linear discriminant follows (within-group whitening plus between-group mean SVD).

Output:

| File | Content |
|---|---|
| `PREFIX.dapc.ind` | Header `#IID POP ASSIGN LD1..LDn_da`: sample, true population, assigned population, discriminant coordinates |
| `PREFIX.dapc.eig` | Discriminant eigenvalues, one per line |
| `PREFIX.dapc.var` | One row `#n_pca n_da var`: retained PC count, discriminant axis count, variance share carried by those PCs |

```
PREFIX.dapc.ind
#IID	POP	ASSIGN	LD1
s1	popA	popA	3.802819344
s2	popA	popA	1.383581045
s3	popA	popA	3.283469062
…

PREFIX.dapc.var            PREFIX.dapc.eig
#n_pca	n_da	var          75.856959395
5	1	0.813187618
```

## Cross-validation

`--xval` selects the PC count automatically. PCA runs once on the full data; each repetition derives a random stream from the fixed seed, samples training individuals per population by the training fraction, refits the discriminant on the training set and scores accuracy on the held-out set. Repetitions are summarized per candidate `n_pca`.

The criterion is `√(mean((success rate − 1)²))`, minimized over candidates; ties go to the largest candidate. The random baseline `1/population count` is reported alongside. Results go to `PREFIX.dapc.xval`: one row per candidate under `#n_pca mean_success`, followed by the summary lines `#best_n_pca`, `#best_mean_success_n_pca` and `#random_chance`.

```
PREFIX.dapc.xval
#n_pca	mean_success
1	0.710000000
2	0.600000000
3	0.400000000
4	0.615000000
5	0.840000000
6	0.735000000
7	0.495000000
8	0.355000000

#best_n_pca	5
#best_mean_success_n_pca	5
#random_chance	0.500000000
```

The chosen PC count overrides `--pcs`, the final DAPC is refit, and one note line is printed:

```
reed pca --dapc --xval: best n_pca=5 (mean success 0.8400), highest mean success at n_pca=5, random chance 0.5000
  note: --pcs 10 was overridden by the cross-validated choice (5)
```

## Notes

- LD pruning is an upstream step: run `reed filter --indep-pairwise` first and feed the pruned file to `pca`. This command does not accept `--indep-pairwise`.
- Large-effect sites associated with the trait contaminate PCs; for association analysis the usual practice is to use PCs computed after removing the target sites.
- Eigenvector signs are arbitrary: compare two sets of PCs by direction (absolute inner product), not element-wise.

## Examples

```bash
# write the first 20 components
reed pca pruned.vcf.gz --pcs 20 -o pcs

# switch standardization and missing-data policy
reed pca pruned.vcf.gz --stand binom --missing meanimpute -o pcs

# force the streaming solver with a memory cap
reed pca pruned.vcf.gz --solver stream --memory 2000 -o pcs

# only part of the samples
reed pca pruned.vcf.gz -S samples.txt -o pcs

# DAPC on 5 components
reed pca pruned.vcf.gz --dapc --pop pops.txt --pcs 5 -o dapc

# DAPC with cross-validation, 100 repetitions
reed pca pruned.vcf.gz --dapc --pop pops.txt --xval --xval-rep 100 -o dapc
```
