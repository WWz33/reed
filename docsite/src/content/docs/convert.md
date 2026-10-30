---
title: convert
---

## Usage

```
reed convert [OPTIONS] --phylip --012 ... -o PREFIX <INPUT>
reed convert --vcf [OPTIONS] -o FILE <HAPMAP>
```

Format conversion. One read, several targets; without `-o` the output prefix is the input file name minus its extension. `--vcf` is the only reverse conversion (hapmap to VCF) and is exclusive with all VCF-input targets.

## Options

### Targets

Give the target and it is written; several can be combined.

```
    --phylip                PHYLIP matrix (relaxed layout) → PREFIX.phy
    --fasta                 FASTA matrix → PREFIX.fasta
    --nexus                 NEXUS (DNA) matrix → PREFIX.nexus
    --nexus-binary          Binary NEXUS (SNP datatype), biallelic diploid
                            sites only → PREFIX.bin.nexus
    --012                   0/1/2/-1 genotype dosage matrix → PREFIX.012,
                            PREFIX.012.pos, PREFIX.012.indv
    --hapmap                Hapmap matrix (11 fixed columns plus one per
                            sample) → PREFIX.hmp.txt
    --msa                   Reference-driven whole-site MSA, one file per contig
                            or region. Needs --reference →
                            PREFIX_<contig>[_<start>-<end>].fasta
    --sfs                   Projected SFS, sharing the implementation of
                            reed popg --sfs. Needs --proj → PREFIX_sfs/dadi,
                            PREFIX_sfs/fastsimcoal2
    --vcf                   Reverse: convert a hapmap input back to VCF;
                            exclusive with every other target → -o FILE
```

### Matrix content

```
    --outgroup NAME         Write this sample first in the alignment
    --min-samples-locus N   Keep a site only if at least N samples have a
                            genotype. Default 4 for phylogenetic targets;
                            --012 and --hapmap apply the gate only when given
    --resolve-iupac         Resolve heterozygous sites to one base; the draw is
                            deterministic (seed 0), so runs repeat
    --write-used-sites      Also write PREFIX.used_sites.tsv listing the sites
                            that entered the matrices and their sample counts
```

### MSA

```
    --reference FASTA       Reference genome; --msa errors without it. Must match
                            the VCF
    --region CHR:B-E        Restrict to this segment (1-based closed interval),
                            repeatable; default the whole contig
    --flank INT             Extend each region by INT bp on both sides, default 0
    --collapse-suffix STR   Truncate sample names at the first STR; samples that
                            collide after truncation are rejected
```

### SFS

```
    --sfs[=folded]          Unfolded spectrum (default, REF as ancestral);
                            --sfs folded folds it by the minor allele
    --proj N[,N...]         Projection size per population in chromosomes; one
                            value broadcasts, several map in first-appearance order
    --sfs-int               Round projected values to integers
    --pop FILE              Population file, sample<TAB>population, no header;
                            without it all samples form one group.
                            Alias --populations
    --force-samples         Unknown samples in the population file only warn and
                            are skipped
```

### Hapmap to VCF

```
    --vcf                   Input is .hmp.txt (gzip/BGZF auto-detected, stdin
                            works too)
-O, --output-type v|z[0-9]  Plain VCF or BGZF; no BCF output
    --no-version            No ##reed_* provenance header lines
    --no-command            No ##reed_command line
```

### Input and samples

```
<INPUT>                     VCF/VCF.gz/BCF (hapmap with --vcf); `-` reads stdin
-s, --samples [^]LIST       Sample subset, comma-separated, ^ prefix complements
-S, --samples-file FILE     Sample list file, one name per line
    --keep FILE             Same as -S, alternative spelling
    --remove FILE           Samples to drop
```

### Output and global options

```
-o, --output PREFIX         Output prefix; default is the input file name
                            without its extension. With --vcf it is the output
                            file. Alias --out
    --threads INT           Worker threads, default all available cores
    --memory MB             In-flight memory budget in MiB, default half of
                            physical memory
    --verbosity INT         0 silent, 1 summary (default), 2 stage timing, 3 debug
```

## General rules

- Phylogenetic matrices and hapmap accept only sites whose REF and every ALT are single characters; multi-character sites are skipped and counted. `--012` skips only multiallelic sites.
- A target's gate affects that target only: sites skipped by the phylogenetic stream still enter `--012`, hapmap and MSA.
- Malformed records (wrong column count, bad POS, unrecognized allele letters) are skipped and counted, with one warning on stderr; totals go to `PREFIX.log`.

## Phylogenetic matrices

Homozygotes write the base itself, heterozygotes the IUPAC ambiguity code (`A/G→R`, `C/T→Y`, `A/C→M`, `G/T→K`, `A/T→W`, `C/G→S`; three bases `V/H/D/B`; four `N`), missingness `N`.

`PREFIX.phy` (relaxed layout, the first line is `samples sites`):

```
3 5
sampleA     ACGTA
sampleB     ACGTC
sampleC     NGTTA
```

A real run (10 samples, 122 sites, sequence rows truncated):

```
10 122
s1    MAGGRCAAMAMTARGNRNCTTCNRWRCGWACANGTTAGACRWGCGCMMCTWTMAAAAT…
s2    MRRAACACAGCWTGRARMAAWNTATGCNNTNGANWATAACRAGCGCCMAAWWCNANGT…
s3    MGAGNCCACGCTAAAMRMATTMTAARMGAWAARGWAAAGAGTRANCANMWNTACTTRT…
…
```

`PREFIX.fasta`:

```
>s1
MAGGRCAAMAMTARGNRNCTTCNRWRCGWACANGTTAGACRWGCGCMMCTWTMAAAATNANNTAATWTMWMAACAMMCANARATMMACTRRAAWWGAACC…
```

`PREFIX.nexus`:

```
#NEXUS

BEGIN DATA;
	DIMENSIONS NTAX=10 NCHAR=122;
	FORMAT DATATYPE=DNA MISSING=N GAP=- ;
MATRIX
s1    MAGGRCAAMAMTARGNRNCTTCNRWRCGWACANGTTAGACRWGCGCMMCTWTMAAAAT…
s2    MRRAACACAGCWTGRARMAAWNTATGCNNTNGANWATAACRAGCGCCMAAWWCNANGT…
;
END;
```

`PREFIX.bin.nexus` (`--nexus-binary`, `0/0→0`, `0/1→1`, `1/1→2`, everything else including `./.` writes `?`):

```
#NEXUS

BEGIN DATA;
	DIMENSIONS NTAX=10 NCHAR=122;
	FORMAT DATATYPE=SNP MISSING=? GAP=- ;
MATRIX
s1    102212001012012?1?2222?111221020?2220202112222112212100002?0??200212111002…
s2    111002020221221011001?20222??2?20?1020021022222100112?0?222022211?00?2100…
;
END;
```

## Genotype matrix (`--012`)

Three files. The leading index is the sample's position in the original VCF header; missing, half-calls and haploids all write `-1`:

```
PREFIX.012                                    PREFIX.012.pos    PREFIX.012.indv
0  1 0 2 2 1 2 0 0 1 … 2 -1 0 -1 -1 2 0 0 0 1  chr1  400         s1
1  1 1 1 0 0 2 0 2 0 … 0 0 1 0 0 1 0 1 2 2 0   chr1  800         s2
2  1 2 0 2 -1 2 2 0 2 … 0 -1 0 0 0 1 0 1 0 1   …                 …
```

## Hapmap (`--hapmap`)

`PREFIX.hmp.txt`, 11 fixed columns plus one per sample:

```
rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10
rs2	A/C	chr1	400	+	NA	NA	NA	NA	NA	NA	M	M	M	C	A	M	N	A	C	M
rs4	A/G	chr1	800	+	NA	NA	NA	NA	NA	NA	…	…	…	…	…	…	…	…	…	…
```

How to read it: in the first row `s4` is `C`, i.e. homozygous ALT for A/C; `s1` is `M`, het A/C; `s7` is `N`, missing. A few details:

- Cells use single-letter IUPAC codes, then an alphabet contraction: only `A C G T N R Y S W K M -` survive; anything else (lowercase letters, three-base codes) writes `N`.
- `rs#` comes from the VCF ID; an ID of `.` writes `chrom:pos`.
- The `alleles` column holds the declared REF/ALT set (slash-joined, deduplicated, `ALT='.'` skipped, `*` written `-`).
- strand writes `+`, the six metadata columns in between write `NA`.

With `--write-used-sites`, `PREFIX.used_sites.tsv` looks like this:

```
#CHROM	POS	NUM_SAMPLES
chr1	400	9
chr1	800	9
```

## Reference MSA (`--msa`)

Each contig (or `--region`) writes `PREFIX_<contig>[_<start>-<end>].fasta`, with records `>contig:start-end_sample`:

```
>chr1:1-1000_s1
AAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA…
>chr1:1-1000_s2
AAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA…
```

Per-site sample values: single-character alleles merge into the IUPAC code; indel alleles keep their literal text, and when all alleles are indels the shortest wins (ties by lexicographic order); a fully missing genotype writes the reference base. Each sample's value is right-padded with `-` to the site's maximum length; no external aligner is called. The input must be sorted by contig; several records at one position take the last one.

## Projected SFS (`--sfs`)

Shares the implementation of `reed popg --sfs`; same prefix and input give identical output.

dadi layout (`PREFIX_sfs/dadi/<pop>-<n>.sfs`, three lines: dimensions, values, mask):

```
7 unfolded
2.7357142857142858 10.776190476190477 26.226190476190478 35.642857142857146 30.773809523809526 14.073809523809524 1.7714285714285716
1 0 0 0 0 0 1
```

fastsimcoal2 layout (`PREFIX_sfs/fastsimcoal2/<pop>_DAFpop0.obs`; DAF for unfolded, MAF for folded):

```
1 observation
d0_0	d0_1	d0_2	d0_3	d0_4	d0_5	d0_6
2.7357142857142858 10.776190476190477 26.226190476190478 35.642857142857146 30.773809523809526 14.073809523809524 1.7714285714285716
```

From two populations on, each pair also gets a two-dimensional joint spectrum (`dadi/<A>-<B>.sfs`, `fastsimcoal2/<base>_jointMAFpop<idxB>_<idxA>.obs`). The directory structure for the output above:

```
out_sfs/dadi/            popA-6.sfs  popB-6.sfs  popA-popB.sfs
out_sfs/fastsimcoal2/    popA_DAFpop0.obs  popB_DAFpop0.obs  out_jointDAFpop1_0.obs
```

## Hapmap to VCF (`--vcf`)

Translation rules:

- `CHROM/POS/ID` are copied from the chrom, pos and rs# columns; REF/ALT follow the declared order in `alleles`, no flipping, no reordering. A row with `alleles=NA` writes `./.` for every genotype.
- Cells expand to GT: ambiguity codes such as `R → 0/1` are looked up and expanded; only `N` and `?` are missing (`./.`); expanded bases must be declared alleles (case-insensitive); `-`, `+` and `0` index into the `alleles` column.
- The output header carries `##fileformat=VCFv4.3` and `##FORMAT=<ID=GT>`, with no `##contig` lines.
- Sample subsets and VCF-input target options are not accepted.

## Examples

```bash
# every nucleotide and dosage matrix in one read
reed convert clean.vcf.gz --phylip --fasta --nexus --nexus-binary --012 --hapmap -o matrices

# phylogenetic matrices: outgroup, per-site sample gate, IUPAC resolution,
# record the sites used
reed convert clean.vcf.gz --phylip --outgroup sampleZ \
    --min-samples-locus 10 --resolve-iupac --write-used-sites -o matrices

# reference-driven MSA over one segment with 500 bp flanks
reed convert clean.vcf.gz --msa --reference ref.fa --region chr1:100000-200000 --flank 500 -o aln

# projected SFS
reed convert clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out

# hapmap back to VCF
reed convert panels.hmp.txt --vcf -Oz6 -o panels.vcf.gz
```
