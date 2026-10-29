# Commands and global conventions

`reed --help` lists commands by group: Filtering (`filter`), Statistics (`stats`, `popg`, `pca`), Relatedness (`distance`), Transform (`convert`, `rename`), File Operations (`index`). `reed <command> --help` lists every option of one command.

## Global options

Every subcommand accepts:

```
<INPUT>                 VCF/VCF.gz/BCF; `-` reads stdin (index rejects stdin,
                        an index needs a seekable file)
    --threads INT       Worker threads, default all available cores
    --memory MB         In-flight memory budget in MiB, default half of
                        physical memory
    --verbosity INT     0 silent; 1 summary (default); 2 stage timing; 3 debug
                        logs (via RUST_LOG, default debug); above 3 is an error
```

`filter` additionally has:

```
    --tmp-dir DIR       Spill directory for multi-pass stdin runs
    --report FILE       JSON run report: filter funnel, sample rules, counters, strategy
    --timing FILE       Per-stage timing (key-value text)
```

`filter` and `rename` also accept `--no-version` (drop every `##reed_*` provenance line) and `--no-command` (drop only `##reed_command`, the raw argv with machine paths).

## Length and number syntax

Length parameters (`--window`, `--step`, `--block-size`, `--flank`, `--ld-decay`, the window of `--indep-pairwise`) accept `100000`, `100k`, `100kb`, `1M`, `1Mb`; suffixes are case-insensitive, k is 10³ and M is 10⁶. Missing-rate thresholds (`--geno`, `--mind`, `--maf`, `--hwe`) are decimals in [0,1].

## Shared option groups

Three groups keep the same spelling and semantics across commands:

- Sample subsets: `-s/--samples [^]LIST` (comma-separated, a `^` prefix complements), `-S/--samples-file FILE`, `--keep FILE`, `--remove FILE`, `--force-samples`. Used by `filter`, `stats`, `popg`, `pca`, `distance` and `convert`; `-S` with `--keep`, `-s` with `-S/--keep`, and `-s ^…` with `--remove` are usage errors.
- Region selection: `-r/--regions LIST`, `-R/--regions-file FILE` (need an index); `-t/--targets LIST`, `-T/--targets-file FILE` (streamed, no index). File formats are auto-detected: three or more columns read as BED (0-based half-open), two as `CHROM POS` (1-based), one as a region string (`Chr01`, `Chr01:100-200`). `#` comment lines and blank lines are skipped; a `^` prefix in a `-T` file means exclusion.
- Site filters and thresholds: `-v/-V`, `--snps-only`, `-m/-M`, `-f`, `-i/-e`, `-q/-Q`, `--geno`, `--maf/--max-maf`, `--mac/--max-mac`, `--hwe`. `filter` implements all of them; `popg` has the same site rules with "filter while computing" semantics; `stats` report mode rejects them and asks for `filter` first.

## Exit codes and error format

| Code | Meaning |
|---|---|
| 0 | success |
| 1 | runtime failure (the message carries a `contig:pos` location) |
| 2 | usage error: conflicting options, undeclared tags, unknown populations, expression compile failure, failed invariant-site check, withdrawn options |

Errors go to stderr as `reed: <reason>`, warnings as `reed: warning: <reason>`. Parse-time errors suggest a rewrite, for example `--max-missing 0.9` is rejected with `use --geno 0.1 (the fraction that may be missing)`; there is no silent alias.

## Evaluation order (filter pipeline)

`filter` runs every record through eight stages in a fixed order. The order decides the funnel counts, not the final result; the one rule that changes the result is `--mind` running after masking:

1. Sample subset and region/position selection
2. Load-time site filters: `-v/-V/--snps-only/-m/-M/-f`
3. Genotype-independent `-i/-e`
4. Genotype masking: `--min-dp/--max-dp/--min-gq`, `--set-gts .` with FORMAT expressions
5. `--mind`
6. `--geno`, `--hwe`, `--maf/--max-maf/--mac/--max-mac`, genotype-dependent `-i/-e`
7. LD pruning (`--indep-pairwise`, needs the global site order)
8. INFO updates (`-x` → automatic AC/AN update → `--fill-tags`) and output

Rules needing first-pass information (`--mind`, LD pruning) scan the input once, then emit on a second pass.
