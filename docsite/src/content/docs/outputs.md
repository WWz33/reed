---
title: Output and logging
---

Output file names, numeric formats, logs and reports follow common conventions collected on this page; the output files and columns specific to each command live on its own page.

## What `-o` means

| Command | `-o` means | Default |
|---|---|---|
| `filter`, `rename` | output file path | stdout |
| `index` | index file path | next to the input |
| `convert` | output prefix (a file path with `--vcf`) | input file name minus extension; stdin without `-o` is an error |
| `stats`, `popg`, `pca`, `distance` | output prefix | no default, an empty value is an error |

`stats` and `popg` have no stdout mode; an empty `-o` is an error (it would write hidden files in the working directory).

## Numeric format

- Estimates with a zero denominator, no coverage or an empty window write `NA`, not 0.
- Coordinates are 1-based closed intervals; BED input is read as 0-based half-open and converted to 1-based.
- Floats carry many decimals (9 for matrices); quantiles and proportions keep full precision.

## Provenance header

`filter` and `rename` insert into the header:

```
##reed_version=<version>
##reed_command=<raw argv>          # --no-command omits it
##reed_filterCommand=<normalized command>  # filter only
##reed_filter_spec=<normalized JSON>       # filter only
##reed_filter_sha256=<hash of the normalized JSON>
##reed_rename=<planes touched>     # rename only
```

`--no-version` removes every `##reed_*` line; `--no-command` removes only the raw-argv line (the one with machine paths).

## PREFIX.log

`popg` in window mode and `convert` both write `PREFIX.log`: reed version, the full command line, method and literature per statistic (for π, the estimator and its missing-data treatment), invariant-site mode, window and record counts, and for `convert` the targets and counters (total records, samples skipped per site, multi-character sites, malformed rows).

## filter's JSON report

`--report FILE` writes the full accounting of one filtering run:

```json
{
  "n_input_records": 500,
  "n_output_records": 373,
  "n_samples_in": 10,
  "n_samples_out": 10,
  "funnel": [
    { "rule": "site.biallelic_gate", "n_passed": 500 },
    { "rule": "after_mask.missing", "n_passed": 373 }
  ],
  "n_masked_gt": 610,
  "masked_by_rule": [
    { "rule": "genotype.mask.DP", "n_masked": 639 }
  ],
  "strategy": "single_pass",
  "elapsed_secs": 0.0
}
```

The funnel lists cumulative passes per rule; rule ids correspond to the stages in [evaluation order](../commands/) (`site.regions`, `after_mask.missing`, `after_mask.ld_prune`, ...). `--mind` sample rules appear separately under `sample_rule`; LD pruning counts under `ld_prune`. `--timing FILE` adds milliseconds per stage.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | success |
| 1 | runtime failure (the message carries a `contig:pos` location) |
| 2 | usage error (conflicting options, undeclared tags, expression compile failure, missing populations, ...) |
