---
title: Filter expressions
---

`filter`'s `-i/--include` and `-e/--exclude` take expressions: `-i` keeps sites where the expression is true, `-e` drops them. Both are repeatable and apply in command-line order.

## Values

| Category | Syntax | Evaluation time |
|---|---|---|
| Record columns | `CHROM` `POS` `ID` `QUAL` `FILTER` `TYPE` | before masking |
| INFO tags | `INFO/TAG` or the bare name `TAG` (INFO is checked first) | before masking |
| FORMAT tags | `FMT/TAG` | only meaningful in per-sample masking with `--set-gts .` |
| Computed values | `N_ALT` `N_SAMPLES` `N_MISSING` `F_MISSING` | `N_MISSING` and `F_MISSING` after masking |
| Allele statistics | `AC` `AN` `AF` `MAC` `MAF` | INFO when present; recomputed from genotypes when the sample set changed |
| Literals | numbers (including negative), strings, `"."` (missing) | |

The recomputation formulas for `AC/AN/AF/MAC/MAF`: `AF = ac/an`, `MAC = an − max(ac, an−ac)`, `MAF = min(af, 1−af)`. An expression that reads AC/AN already present in INFO triggers a stale-INFO note suggesting `--maf/--mac` or `-x INFO/AC,INFO/AN` first.

## Operators

| Category | Operators |
|---|---|
| Comparison | `==` `=` `!=` `<` `<=` `>` `>=` |
| Pattern | `~` matches, `!~` does not (literal patterns with `^`/`$` anchors only; regex metacharacters `.*+?[](){}\|` are errors) |
| Logic | `&&` `\|\|` `!`, parentheses group; precedence `!` > comparison > `&&` > `\|\|` |

Semantics:

- Comparisons with missing values are always false. `QUAL=="."` matches sites with a missing QUAL; `QUAL>30` is false for a missing QUAL.
- Flag INFO tags read and write by presence: `DB=1` tests presence, `DB=0` tests absence. A bare flag (a lone `DB`) is a compile error.
- Multi-value INFO (`Number=A/R/G/.` or a count >1) compares as "any element satisfies".
- A `Number=1` tag given several values counts as missing; a value that fails numeric parsing counts as missing, never string-compared.

## Not supported

- FORMAT subscripts `FMT/DP[0]`, `FMT/DP[*]` (parse-time errors);
- arbitrary arithmetic (no `+ - * /`);
- general regular expressions;
- expressions under `--set-gts .` mixing FORMAT values with site-level values derived from genotypes; FORMAT tags must be `Number=1` and numeric.

## Examples

```bash
# QUAL and minor allele frequency thresholds together
reed filter in.vcf.gz -i 'QUAL>30 && MAF>=0.05' -Oz -o out.vcf.gz

# shallow genotypes become missing, the site stays
reed filter in.vcf.gz -e 'FMT/DP<5' --set-gts . -Oz -o masked.vcf.gz

# same, combining DP and GQ
reed filter in.vcf.gz -e 'FMT/DP<5 || FMT/GQ<20' --set-gts . -Oz -o masked.vcf.gz

# variant type and ID
reed filter in.vcf.gz -i 'TYPE="snp" && ID!="."' -Oz -o snps_named.vcf.gz

# FILTER flag presence
reed filter in.vcf.gz -e 'LowQual=1' -Oz -o keep.vcf.gz

# sites with a missing QUAL handled separately
reed filter in.vcf.gz -e 'QUAL=="."' -Oz -o defined_qual.vcf.gz

# two hard-filter sets (different thresholds for indels and SNPs)
reed filter in.vcf.gz -e 'QD<2 || QUAL<30 || SOR>3 || FS>60 || MQ<40 || MQRankSum<-12.5 || ReadPosRankSum<-8' -Oz -o snp.vcf.gz
reed filter in.vcf.gz -e 'QD<2 || QUAL<30 || FS>200 || ReadPosRankSum<-20' -Oz -o indel.vcf.gz
```
