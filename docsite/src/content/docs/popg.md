---
title: popg
---

## Usage

```
reed popg [OPTIONS] --pop FILE --pi --fst ... -o PREFIX <INPUT>
```

Windowed population-genetics statistics. One run can compute several statistics at once, each writing its own file; a `PREFIX.log` records the command, version, method, literature and filter counts for every statistic. Alias `popgen`. Invariant sites are assumed present by default; without them the run is refused.

## Options

### Populations and statistics

```
    --pop FILE              sample<TAB>population, two columns, no header.
                            Optional when only --sfs or LD statistics are
                            requested: all samples then form one group named all.
                            Alias --populations
    --pi                    π per population (Nei & Li 1979)
    --dxy                   dxy per population pair (Nei 1987)
    --theta                 Watterson θ per population (Watterson 1975)
    --tajima-d              Tajima's D per population (Tajima 1989)
    --fst [wc84|hudson]     Fst for all population pairs; wc84 by default
                            (Weir & Cockerham 1984), hudson is the ratio-of-sums
                            form (Hudson 1992 / Bhatia 2013)
    --fst-global            With --fst: report the merged r-population WC84
                            global Fst to PREFIX_fst_global.txt
    --pbs A,B,C             PBS for the triple, A being the focal population
                            (Yi et al. 2010); repeatable
    --pbs-normed            Report PBSn1 (Malaspinas et al. 2016) instead of
                            the raw PBS
    --f2 A,B                Patterson f2(A, B); repeatable
    --f3 C;A,B              Patterson f3(C; A, B); repeatable
    --f4 A,B;C,D            Patterson f4(A, B; C, D); repeatable
    --dstat A,B;C,D         Patterson's D(A, B; C, D); repeatable
    --f4-ratio A,B;C,D;X    f4-ratio, the alpha estimate; repeatable
    --block-size SIZE       Block-jackknife block length for f-statistics,
                            default 5M
    --ld-decay [MAX_BP]     LD decay: pairwise r² within MAX_BP, binned at 1 kb,
                            default 100k
    --ld-pairs              Also write the per-pair LD table PREFIX.ld; requires
                            --ld-decay
    --sfs [folded]          Projected SFS, genome-wide; writes PREFIX_sfs/dadi
                            and PREFIX_sfs/fastsimcoal2. Unfolded by default
                            (REF as ancestral), --sfs folded folds it.
                            Requires --proj
    --proj N[,N...]         Projection sizes in chromosomes. One value
                            broadcasts to every population; several values map
                            to populations in first-appearance order
    --sfs-int               Round projected values to integers
```

### Windows and sites

```
    --window SIZE           Window size, default 100,000 bp
    --step SIZE             Step, default equal to the window size (0 is also
                            treated as the window size)
    --bed-file FILE         BED3 (0-based half-open) custom windows, one window
                            per interval; exclusive with --window/--step
    --sites-file FILE       CHROM<TAB>POS (1-based); count only the listed sites
    --levels LIST           Which rows to report: window,contig,genome, default
                            window; a contig row folds the whole contig, a
                            genome row folds everything
    --skip-empty            Drop windows no site falls into (by default they are
                            written with NA)
    --contig-lengths FILE   When the header lacks ##contig lines, take lengths
                            from a .fai-style name<TAB>length file; still
                            missing is an error, nothing is invented
```

### Invariant sites

```
    --variant-only          Declare an all-variant input: positions the window
                            does not contain count as monomorphic sites with all
                            samples called, entering the π/dxy/θ denominators
    --bypass-invariant-check
                            Compute only over the sites the file holds, warning
                            that π and dxy are overestimated
```

With neither flag, an input without `ALT = "."` records fails at run time and the error gives both routes. Requesting only SFS or f-statistics switches to `--variant-only` automatically.

### Input and regions

```
<INPUT>                     VCF/VCF.gz/BCF; `-` reads stdin
-r, --regions LIST          Keep these regions, comma-separated; needs an index
                            next to the input
-R, --regions-file FILE     Region list file; needs an index
-t, --targets LIST          Same region syntax, streamed, no index needed
-T, --targets-file FILE     Region list file; a leading ^ excludes. No index
```

### Samples and sites

```
-s, --samples [^]LIST       Sample subset, comma-separated, ^ prefix complements
-S, --samples-file FILE     Sample list file, one name per line
    --keep FILE             Same as -S, alternative spelling
    --remove FILE           Samples to drop
    --force-samples         Unknown samples in the population file only warn and
                            are skipped; default is an error
-v, --types LIST            Keep these variant types: snps, mnps, indels,
                            other, ref
-V, --exclude-types LIST    Drop these variant types
    --snps-only[=just-acgt] Keep SNP sites; =just-acgt requires every REF/ALT
                            to be in A/C/G/T
-m, --min-alleles INT       Minimum allele count, REF+ALT
-M, --max-alleles INT       Maximum allele count, REF+ALT
-f, --apply-filters LIST    Keep records whose FILTER is in LIST
-i, --include EXPR          Keep records that satisfy the expression, repeatable
-e, --exclude EXPR          Drop records that satisfy the expression, repeatable
-q, --min-af FLOAT[:type]   Lower bound on allele frequency, :type as in filter
-Q, --max-af FLOAT[:type]   Upper bound on allele frequency
    --geno FLOAT            Drop sites whose missing-allele fraction exceeds
                            FLOAT (after masking)
    --maf FLOAT             Lower bound on the minor allele frequency
    --max-maf FLOAT         Upper bound on the minor allele frequency
    --mac INT               Lower bound on the minor allele count
    --max-mac INT           Upper bound on the minor allele count; 0 is the
                            internal "off" value and is rejected
    --hwe P [midp]          Lower bound on the HWE exact-test p
-g, --genotype [^](hom|het|miss)
                            Select sites by genotype composition, repeatable,
                            combined with AND
    --min-dp INT            Genotype mask, FMT/DP lower bound. Alias --vcf-min-dp
    --max-dp INT            Genotype mask, FMT/DP upper bound. Alias --vcf-max-dp
    --min-gq INT            Genotype mask, FMT/GQ lower bound. Alias --vcf-min-gq
    --set-gts .             Per-sample masking with -i/-e; only `.` is accepted
    --indep-pairwise W[kb] [S] R2
                            Not accepted here: the error asks for a prior
                            reed filter --indep-pairwise, then feed the pruned
                            file in
```

### Output and global options

```
-o, --output PREFIX         Output prefix. Required in window mode; -o - is an
                            error; a bare command (no statistic) writes the
                            summary TSV to stdout. Alias --out
    --threads INT           Worker threads, default all available cores
    --memory MB             In-flight memory budget in MiB, default half of
                            physical memory
    --verbosity INT         0 silent, 1 summary (default), 2 stage timing, 3 debug
```

Population ids are numbered in first-appearance order in the population file; duplicate names, empty populations and a sample assigned to two populations are usage errors.

## Statistics and algorithms

### π (`--pi`, Nei & Li 1979)

For every pair of chromosomes within a population, the number of differences at a site is `n_ref·n_alt` and the number of comparisons is `C(n,2)` (`n` being the called chromosome count of that population at that site). The window value is

```
π = Σ site differences / Σ site comparisons
```

The denominator counts only pairs where both chromosomes are called at that site (pairwise complete); the number of missing pairs goes to the output column `count_missing`. A window with no comparable site writes `NA`.

### dxy (`--dxy`, Nei 1987)

For population A and B, the number of differences per cross-population chromosome pair is `a_ref·b_alt + a_alt·b_ref` and the number of comparisons is `(a_ref+a_alt)·(b_ref+b_alt)`. The window value is again a ratio of sums across sites. With `--variant-only`, the comparison count gains `(window length − file site count)·2N_A·2N_B`.

### Fst (`--fst [wc84|hudson]`, wc84 by default)

Both estimators share one site traversal and all output goes to `PREFIX_fst.txt`.

WC84 (Weir & Cockerham 1984, r = 2 pairwise): variance components from allele copy counts. Let `p̄` be the pooled allele frequency, `h̄` the mean heterozygosity and `n̄` the mean gene copy number (each population's called copies divided by 2, averaged):

```
s² = Σ nᵢ(pᵢ−p̄)² / (n̄(r−1))
a = (n̄/n_C)·(s² − (p̄(1−p̄) − (r−1)s²/r − h̄/4)/(n̄−1))
b = (n̄/(n̄−1))·(p̄(1−p̄) − (r−1)s²/r − (2n̄−1)h̄/(4n̄))
c = h̄/2
window value = Σa / (Σa + Σb + Σc)
```

`n_C = (Σnᵢ − Σnᵢ²/Σnᵢ)/(r−1)`. A population with fewer than 2 called chromosomes contributes only the `c` component at that site.

Hudson (Hudson 1992; the ratio-of-sums form of Bhatia et al. 2013): the numerator is the between-group mean difference minus the within-group mean difference, the denominator is the between-group mean difference:

```
within   mpd_a = n_ref_a·n_alt_a / C(n_a,2)     (requires n_a ≥ 2)
between  mpd_ab = (a_ref·b_alt + a_alt·b_ref)/(n_a·n_b)
Fst  = Σ(mpd_ab − (mpd_a+mpd_b)/2) / Σ mpd_ab
```

A site missing a numerator term contributes 0 to the numerator sum; the denominator sum accumulates independently.

Global Fst (`--fst-global`) with three or more populations computes the merged r-population WC84 (the formulas above with general r). Hudson has no r-population form, so `--fst-global` with `--fst hudson` is an error; with two populations the global value equals the pairwise value and only a note is printed.

### Watterson θ (`--theta`, Watterson 1975)

Sites are binned by observed allele count `n` (the population's called chromosome count); segregating sites with `0 < alt < n` count as one mutation. The window value is

```
θ_W = (Σ site mutations/a₁(n)) / n_sites,   a₁(n) = Σ_{k=1}^{n−1} 1/k
```

The harmonic table is precomputed up to 2¹⁶ and summed directly with compensated summation beyond that. Missing data enters through each site's own `n`; sites are never dropped wholesale. With `--variant-only` the denominator is the window length. Fixed-ALT sites contribute 0.

### Tajima's D (`--tajima-d`, Tajima 1989)

`raw_pi` is the raw π sum, `raw_theta` the raw θ_W of the previous section, `S` the window's total mutation count. The variance constants are evaluated at the window's mean allele count (rounded half to even):

```
a₁ = Σ1/k,  a₂ = Σ1/k²   (k = 1..n−1)
b₁ = (n+1)/(3(n−1)),      b₂ = 2(n²+n+3)/(9n(n−1))
c₁ = b₁ − 1/a₁,           c₂ = b₂ − (n+2)/(a₁n) + a₂/a₁²
e₁ = c₁/a₁,               e₂ = c₂/(a₁²+a₂)
variance = e₁·S + e₂·S(S−1),   D = (raw_pi − raw_theta)/√variance
```

A window with no sites, no mutations or a mean allele count < 2 writes `NA`. The output carries `tajima_d_stdev`.

### PBS (`--pbs A,B,C`, Yi et al. 2010)

The input is three pairwise windowed Hudson Fst values (each pair decides its own sample coverage; no three-way intersection). Fst is truncated to `[0, 0.99999]` and transformed as `t(f) = −ln(1−f)`:

```
PBS = (t_AB + t_AC − t_BC)/2      (A is the focal population)
```

`--pbs` is repeatable; each triple writes `PREFIX_pbs.<A>-<B>-<C>.txt` and a single triple writes `PREFIX_pbs.txt`. The output carries the three untruncated pairwise Fst columns.

`--pbs-normed` replaces PBS with PBSn1 (Malaspinas 2016): `PBS/(1 + (t_AB+t_AC+t_BC)/2)`.

### f-statistics (unbiased-corrected form of Patterson et al. 2012)

Allele frequency `p = n_alt/an` (by ALT), with `ĥ = n_ref·n_alt/(an(an−1))` the finite-sample heterozygosity correction:

| Statistic | Per-site contribution |
|---|---|
| `--f2 A,B` | `(pA−pB)² − ĥA/anA − ĥB/anB` |
| `--f3 C;A,B` | `T = (pC−pA)(pC−pB) − ĥC/anC`; the normalized variant also reports `ΣT/ΣB`, `B = 2ĥC` |
| `--f4 A,B;C,D` | `(pA−pB)(pC−pD)` |
| `--dstat A,B;C,D` | `Σ(pA−pB)(pC−pD) / Σ(pA+pB−2pA·pB)(pC+pD−2pC·pD)` |
| `--f4-ratio A,B;C,D;X` | f4 numerator sum including X ÷ f4 numerator sum without X |

All f-statistics take block-jackknife standard errors: the genome is cut into fixed-length blocks (5,000,000 bp by default, changeable with `--block-size`), each block is deleted in turn and pseudovalues computed as

```
v_k = (Σnum − num_k)/(Σden − den_k),   SE = √((n−1)/n · Σ(v_k − v̄)²),   z = estimate/SE
```

Blocks are keyed by `(contig, pos/block length)`, so multi-threaded block merging is safe. Output is `PREFIX.fstats` with columns `stat pops estimate se z n_blocks n_snps`; `--f2/--f3/...` may all be repeated and are summarized in one file. Zero denominators write `NA`.

### Projected SFS (`--sfs [folded] --proj N[,N…]`)

Genome-wide (no windows). One-dimensional per population; from two populations on, each pair also gets a two-dimensional joint spectrum. Projection uses the hypergeometric distribution

```
P(k | n, m, alt) = C(alt,k)·C(m−alt,n−k)/C(m,n)
```

(with m the site's called count and n the projection target), computed by a mode-anchored ratio recursion, truncating relative weights below 1e-15 of the peak. Output comes in two layouts:

- `PREFIX_sfs/dadi/<pop>-<n>.sfs` (1D), `PREFIX_sfs/dadi/<A>-<B>.sfs` (2D), three lines: dimensions, values, mask;
- `PREFIX_sfs/fastsimcoal2/<pop>_MAFpop0.obs` (folded) / `<pop>_DAFpop0.obs` (unfolded), `<base>_jointMAFpop<idxB>_<idxA>.obs` (2D, transposed layout).

The mask line marks the structural zeros (k=0 and the unfolded top, the folded upper triangle).

### LD decay (`--ld-decay [MAX_BP]`, `--ld-pairs`)

r² is computed on unphased 0/1/2 dosages (an integer implementation of the unphased Rogers-Huff r²), using samples called at both sites (pairwise removal). `MAX_BP` defaults to 100,000; distances are binned at a fixed 1 kb, each bin recording the mean r² and pair count. Output is `PREFIX.ld.decay` with columns `pop bin_start bin_end n_pairs mean_r2`: all samples first (pop = `all`), then one group per population when several exist. `--ld-pairs` also writes the per-pair table `PREFIX.ld` (columns `#CHROM_A POS_A ID_A CHROM_B POS_B ID_B UNPHASED_R2`).

## Window mechanics

The window grid is anchored at integer multiples of the step: the first window starts at the nearest anchor before the first site and the window end is clipped to the contig length. Output coordinates are 1-based closed intervals in `window_pos_1` and `window_pos_2`.

## Output files and columns

| File | Columns |
|---|---|
| `PREFIX_pi.txt` | `pop chromosome window_pos_1 window_pos_2 avg_pi no_sites count_diffs count_comparisons count_missing` |
| `PREFIX_dxy.txt` | `pop1 pop2 chromosome window_pos_1 window_pos_2 avg_dxy no_sites count_diffs count_comparisons count_missing` |
| `PREFIX_fst.txt` | `pop1 pop2 chromosome window_pos_1 window_pos_2 avg_wc_fst \| avg_hudson_fst no_snps` |
| `PREFIX_fst_global.txt` | `pops chromosome window_pos_1 window_pos_2 avg_wc_fst no_snps` |
| `PREFIX_watterson_theta.txt` | `pop chromosome window_pos_1 window_pos_2 avg_watterson_theta no_sites raw_watterson_theta no_var_sites weighted_no_sites` |
| `PREFIX_tajima_d.txt` | `pop chromosome window_pos_1 window_pos_2 tajima_d no_sites raw_pi raw_watterson_theta tajima_d_stdev` |
| `PREFIX_pbs.txt` (multiple triples: `PREFIX_pbs.<A-B-C>.txt`) | `popA popB popC chromosome window_pos_1 window_pos_2 pbs no_snps fst_12 fst_13 fst_23` |
| `PREFIX.fstats` | `stat pops estimate se z n_blocks n_snps` |
| `PREFIX.ld.decay` | `pop bin_start bin_end n_pairs mean_r2` |
| `PREFIX.ld` | `#CHROM_A POS_A ID_A CHROM_B POS_B ID_B UNPHASED_R2` |

`no_sites` is the number of file sites in the window (multiallelic records that were skipped included); `no_snps` is every site fed to the statistic. Missing values write `NA`. Negative Fst is normal: variance-component estimates in small windows can fall below zero, and are usually read as 0.

## Examples

```bash
# several statistics in one run (whole-genome calls input)
reed popg allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k --threads 56 -o out

# variant-only input with half-overlapping windows
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --fst hudson --window 100k --step 50k -o out

# global WC84 Fst
reed popg clean.vcf.gz --pop pops.txt --variant-only --fst --fst-global --window 1M -o out

# PBS: focal population popA against two references
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --pbs popA,popB,popC --pbs-normed --window 50k -o out

# f4-ratio test
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --f4-ratio popA,popB;popC,popD;popX --block-size 2M -o out

# custom windows (BED) and a site list
reed popg clean.vcf.gz --pop pops.txt --variant-only --pi --bed-file windows.bed -o out
reed popg clean.vcf.gz --pop pops.txt --variant-only --pi --sites-file sites.txt -o out

# filter and compute in one pass
reed popg raw.vcf.gz --pop pops.txt --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 \
    --pi --window 100k -o out
```

`out_pi.txt`:

```
pop	chromosome	window_pos_1	window_pos_2	avg_pi	no_sites	count_diffs	count_comparisons	count_missing
popA	chr1	1	10000	0.0003778491492837536	9	170	320	85
popB	chr1	1	10000	0.00040893523488217553	9	184	354	51
popA	chr1	10001	20000	0.0005067049510407453	11	228	461	34
popB	chr1	10001	20000	0.00046229207969204235	11	208	427	68
```

`out_tajima_d.txt`:

```
pop	chromosome	window_pos_1	window_pos_2	tajima_d	no_sites	raw_pi	raw_watterson_theta	tajima_d_stdev
popA	chr1	1	10000	2.0866267869882233	9	4.843650793650793	3.3423177051634436	0.7195024514442903
popB	chr1	1	10000	1.9896866772136814	9	4.70952380952381	3.2779393676625217	0.7195024514442903
```

`out_fst.txt`:

```
pop1	pop2	chromosome	window_pos_1	window_pos_2	avg_wc_fst	no_snps
popA	popB	chr1	1	50000	-0.03926597324664331	63
popA	popB	chr1	50001	100000	-0.026026648105796968	59
```

`out.fstats` (one line per f-statistic or D test):

```
stat	pops	estimate	se	z	n_blocks	n_snps
f4	popA,popB,popC,popD	0.4027777777777775	NA	NA	0	122
dstat	popA,popB,popC,popD	0.013145965548504073	NA	NA	0	122
```

This is the 122-site demo data: the whole contig falls into one block, the block count is 0 and the standard error writes `NA`. Real data spans dozens of blocks before `se` and `z` carry content.

`out.ld.decay` (all samples first, then one group per population with several populations):

```
pop	bin_start	bin_end	n_pairs	mean_r2
all	0	1000	118	0.15795037970268636
all	1000	2000	144	0.15256733855146404
all	2000	3000	155	0.10808903029910441
all	3000	4000	144	0.13556198515852427
```

## Checks and common errors

- `--proj` needs `--sfs`; `--sfs` needs `--proj N`; the number of projection values must equal the number of population groups.
- `--block-size` needs f-statistics; `--pbs-normed` needs `--pbs`; `--ld-pairs` needs `--ld-decay`.
- A `--pbs` triple must be three distinct known population names.
- f-statistic combinations: `A,B` / `C;A,B` (also `C;A;B`) / `A,B;C,D` / `A,B;C,D;X`; every population name must appear in the `--pop` file.
- Window mode without `-o`, or with `-o -`, is an error (window mode writes files, never stdout); a bare command (no statistic) writes the summary TSV to stdout.
