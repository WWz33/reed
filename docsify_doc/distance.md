# distance

## Usage

```
reed distance [OPTIONS] --type KIND -o PREFIX <INPUT>
```

Individual-by-individual matrices: allele-difference distance, IBS, GRM, covariance, KING, PI_HAT, relatedness, discordance; `--tree upgma|nj` builds a Newick tree from a distance matrix and `--bootstrap` resamples sites for branch support. Defaults are `--type grm` and `--format triangle`.

The input gate matches `pca`: only biallelic sites with diploid GT are counted, and sites with no called ALT are skipped.

## Options

### Matrices

```
<INPUT>                     VCF/VCF.gz/BCF; `-` reads stdin
    --type KIND             Matrix type, default grm. Values: allele-ct (alias
                            dist), ibs, 1-ibs, grm (alias rel), cov, king,
                            pihat (alias genome), phi (alias relatedness2),
                            ajk (alias relatedness), discordance (alias
                            gtcheck). Formulas below
    --format FMT            Layout for matrix outputs: triangle (default, lower
                            triangle with diagonal), square, square0 (square
                            with the upper right zeroed), bin (little-endian f64)
    --stand binom2|binom    Standardization, as in pca, default binom2
    --mean-impute           grm/cov only: divide by all sites instead of each
                            pair's jointly-called count
    --unbounded             pihat only: skip the [0,1] clipping of PI_HAT
    --nudge                 pihat only: when PI_HAT² < Z2, snap to the nearest
                            feasible configuration
    --mismatch-count        discordance only: report the mismatch count instead
                            of the probability score
    --error-probability INT discordance only: Phred-scaled genotyping error
                            rate, default 40
```

### Trees and bootstrap

```
    --tree upgma|nj         Build a tree instead of writing the matrix, to
                            PREFIX.nwk. upgma is average linkage (rooted), nj
                            is neighbor joining (unrooted)
    --digits INT            Decimal places in Newick branch lengths, default 10;
                            0 is an error (every branch would round to zero)
    --bootstrap INT         Bootstrap: resample sites with replacement this many
                            times, support becomes internal node labels
    --bootstrap-seed INT    Bootstrap random seed, default 42
    --bootstrap-cutoff FLOAT
                            Support below this percentage gets an empty label,
                            default 0 (keep all)
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

## Matrix types and formulas

Standardization uses `--stand binom2` (`√(2p(1−p))`) or `--stand binom` (`√(p(1−p))`), with `p` estimated from called dosages. Some types need the raw dosage, recovered from the standardized value as `dose = round(W·sd + mean)`, taking 0/1/2 with -1 for missing.

| `--type` | Formula | Output file |
|---|---|---|
| `allele-ct` (alias `dist`) | Sum of allele differences over jointly called sites `Σ|x_i − x_j|`, with HWE weight correction for missingness: weight `w = p(1−p)(p²−p+1)`, the denominator is the weight sum minus the weights of either-side and both-side missing, multiplied back by the total weight | `PREFIX.dist` |
| `ibs` | `1 − d/2m`, d the allele-ct distance, m the site count | `PREFIX.mibs` |
| `1-ibs` | `d/2m` | `PREFIX.mdist` |
| `grm` (alias `rel`) | `A_ij = Σ W_is·W_js / denominator`; the denominator is the jointly-called site count by default, `--mean-impute` switches to all sites m. The diagonal is `ΣW²/m`, not forced to 1 | `PREFIX.rel` |
| `cov` | Same as grm but without dividing by sd (centering kept), giving covariance scale | `PREFIX.rel` |
| `king` | Counts four genotype combinations at jointly called sites (ibs0, double het, het1-hom2, het2-hom1), `φ = 0.5 − (4·ibs0 + h1h2 + h2h1) / (4·(hethet + min(h1h2, h2h1)))`; a zero denominator writes NaN | `PREFIX.king` |
| `pihat` (alias `genome`) | Solves the IBS coefficients per pair: `z0 = ibs0/(e00·n)`, `z1 = ibs1/(e11·n)`, `z2 = (oo − n(z0·e02 + z1·e12))/n`, expected coefficients accumulated from site allele frequencies; clipped to [0,1] (`--unbounded` turns clipping off; `--nudge` snaps to `z0=(1−π)², z1=2π(1−π), z2=π²` when `π² < z2`); `PI_HAT = z1/2 + z2` | `PREFIX.genome` |
| `phi` (alias `relatedness2`) | `φ = (N_AaAa − 2·N_AAaa)/(N_Aa[i] + N_Aa[j])`, the denominator being the sum of both samples' heterozygous site counts (marginal counts); diagonal 0.5 | `PREFIX.relatedness2` |
| `ajk` (alias `relatedness`) | Only sites with frequencies strictly inside (0,1): `A_ij = mean (a−2p)(b−2p)/(2p(1−p))`; the diagonal adds a `x² − (1+2p)x + 2p²` term to the same numerator | `PREFIX.relatedness` |
| `discordance` (alias `gtcheck`) | Genotype error model `e = 10^(−E/10)` (`--error-probability` default 40); per site, `cost(a,G)+cost(b,G)` is minimized over every possible true genotype and summed; `--mismatch-count` counts mismatches directly | `PREFIX.gtcheck` |

Matrix outputs (dist, mibs, mdist, rel, king) carry a `PREFIX.<suffix>.id` sidecar (two sample-name columns). Table outputs (genome, relatedness2, relatedness, gtcheck) follow their own table layouts; the DST/PPC/RATIO columns of `.genome` are placeholder 0.

## Output layout

`--format` applies to matrix outputs only:

| Format | Layout |
|---|---|
| `triangle` (default) | Lower triangle with the diagonal: row i holds `d(i,0)..d(i,i)`, 9 decimals, tab-separated |
| `square` | Full N×N matrix |
| `square0` | Matrix with the upper right zeroed |
| `bin` | n² little-endian f64 values, no header, directly readable from binary |

## Trees

`--tree` accepts only distance matrices with a zero diagonal. Similarity types (grm, cov, king, ibs, pihat, phi) are rejected and the error lists usable distance types: `1-ibs`, `allele-ct`, `ajk`, `discordance`.

| Method | Algorithm |
|---|---|
| `nj` | Neighbor joining: merge the pair minimizing the Q matrix `Q(i,j) = (n−2)d(i,j) − S_i − S_j`, new branch `d(i,u) = d(i,j)/2 + (S_i−S_j)/(2(n−2))`, reduction `d(u,k) = (d_ik + d_jk − d_ij)/2`; negative branches are kept, not truncated; the last three taxa merge into a central node; the result is unrooted with `2n−3` edges |
| `upgma` | Average linkage (Lance-Williams update), the merged distance split in half for branch lengths; ultrametric input is recovered exactly; the result is rooted |

Before building, the distance matrix is checked: square, finite, zero diagonal, symmetric (tolerance 1e-9). Trailing zeros are trimmed in Newick output and reserved characters in labels are single-quoted.

## Bootstrap support

- Support is counted per branch: the fraction of replicates in which an internal node's descendant tip set appears (for unrooted trees, the split containing tip 0), written as a percentage; the root label is 100.
- Each replicate re-standardizes, reruns the same `--type` and rebuilds the same tree.
- Failed replicates (for example no usable site after resampling) are skipped and counted; the support denominator shrinks accordingly.

## Examples

```bash
# genetic relationship matrix
reed distance clean.vcf.gz --type grm -o rel

# covariance scale (no division by sd)
reed distance clean.vcf.gz --type cov -o cov

# kinship matrix (autosomal sites expected)
reed distance clean.vcf.gz --type king -o kin

# relatedness tables
reed distance clean.vcf.gz --type phi -o rela
reed distance clean.vcf.gz --type ajk -o yang

# distance matrix and IBS
reed distance clean.vcf.gz --type allele-ct -o d
reed distance clean.vcf.gz --type 1-ibs -o ibsd

# discordance with Phred error 30
reed distance clean.vcf.gz --type discordance --error-probability 30 -o disc

# full square and binary output
reed distance clean.vcf.gz --type grm --format square -o rel
reed distance clean.vcf.gz --type 1-ibs --format bin -o ibsd

# NJ tree, 100 bootstraps, labels kept only at support >= 50
reed distance clean.vcf.gz --type 1-ibs --tree nj \
    --bootstrap 100 --bootstrap-cutoff 50 -o tree
```

`rel.rel` (`--format triangle` default layout, row i up to column i) and `rel.rel.id`:

```
1.310800633
-0.206295975	1.437578015
-0.296796369	-0.411552149	1.284675298

s1	s1
s2	s2
s3	s3
```

`rela.relatedness2` (`--type phi`, both directions written):

```
INDV1	INDV2	N_AaAa	N_AAaa	N1_Aa	N2_Aa	RELATEDNESS_PHI
s1	s1	37	0	37	37	0.500000000
s1	s2	12	26	37	31	-0.588235294
s1	s3	6	29	37	30	-0.776119403
```

`kin.genome` (`--type pihat`, strictly lower triangle; DST/PPC/RATIO are placeholder columns):

```
FID1	IID1	FID2	IID2	RT	EZ	Z0	Z1	Z2	PI_HAT	PHE	DST	PPC	RATIO
s2	s2	s1	s1	OT	NA	1.0000	0.0000	0.0000	0.0000	-9	0	0	0
s3	s3	s1	s1	OT	NA	1.0000	0.0000	0.0000	0.0000	-9	0	0	0
```

`tree.nwk`, with bootstrap support at internal nodes:

```
(((s1:0.2010492714,s4:0.2172579309)35:0.01874847015,
  ((s2:0.2185298336,s6:0.2104283332)25:0.01339241676,
    (s5:0.182356414,s9:0.1875419981)60:0.03474948761)20:0.01475450946)20:0.008743806348,
 ((s3:0.1950385303,s10:0.1764491831)40:0.01806552434,s8:0.2309362055)20:0.009008248687,
 s7:0.2161362545)100;
```
