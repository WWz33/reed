# rename

## Usage

```
reed rename [OPTIONS] <INPUT>
```

Identifier edits: sample names, contig names, variant IDs, header lines. Alias `reheader`. Output goes to stdout by default; genotype columns take part in no edit, so the data never changes.

## Options

### header

```
-h, --header FILE           Replace the whole header with the file's contents;
                            its #CHROM line must match the input's samples
                            exactly (same names, same order) or the command
                            errors. -h is --header here, use --help for help
-f, --fai FILE              Set ##contig lengths from a .fai (name<TAB>length):
                            existing lines get length=, contigs absent from the
                            header are inserted before #CHROM
    --add-meta LINE         Append a ##key=value header line, repeatable
    --remove-meta KEY       Remove every ##KEY= line, repeatable
```

### contig

```
    --rename-chrs FILE      old<TAB>new mapping (whitespace-separated and #
                            comments accepted); ##contig lines and the CHROM
                            column change together. Old names must exist in the
                            header
    --chr-style plain|ucsc|mt-is-m
                            Normalize contig names: plain strips the chr prefix
                            (1, MT); ucsc adds chr (chr1, chrMT); mt-is-m adds
                            chr and writes the mitochondrion as chrM. The
                            --rename-chrs mapping applies first, then the style
```

Duplicate contig names after renaming are rejected (they would produce duplicate `##contig` lines and merge two contigs' records). A stale index next to the output gets a warning suggesting `reed index -f`.

### ID

```
    --set-id PATTERN        Rewrite every ID. The template takes %CHROM %POS
                            %REF %ALT and %% for a literal %, e.g.
                            %CHROM:%POS:%REF:%ALT
    --set-missing-var-ids PATTERN
                            Same, but only fills IDs that are `.`
    --update-name FILE      Rename IDs by an old<TAB>new mapping; any key that
                            never matched at the end of the run is an error (a
                            mistyped key cannot pass silently)
    --idspace-to SEP        Replace : in IDs with SEP
    --rm-dup-ids            Suffix repeated IDs with _1, _2…; . is not
                            deduplicated
    --clear-ids             Set every ID to .
```

### Samples

```
-s, --samples PAIRS         Inline rename pairs, comma-separated, each pair two
                            whitespace-separated names:
                            -s "old new,old2 new2"
-N, --samples-file FILE     old<TAB>new, one pair per line
    --replace-sample-name PATTERN
                            Literal rewrite (no regex): pre* strips a prefix,
                            *suf strips a suffix, otherwise every occurrence of
                            the text is removed. Regex metacharacters are usage
                            errors
```

Rename pair checks: the old name must exist; the new name cannot be empty; a collision after renaming is an error.

### Output

```
-o, --output FILE           Output file; stdout by default. Without -O the
                            extension decides compression. Do not point it at
                            the input itself. Alias --out
-O, --output-type v|z[0-9]  Plain VCF or BGZF, optionally compressed. BCF input
                            is converted to VCF first; -O b/-O u are rejected
-W, --write-index[=csi|tbi] Build the index while writing, from the renamed
                            header contig table (rejected when the header was
                            replaced and has no ##contig lines); needs -o FILE
                            and BGZF output
    --no-version            No ##reed_* provenance header lines
    --no-command            No ##reed_command line
    --threads INT           Worker threads, default all available cores
    --memory MB             In-flight memory budget in MiB, default half of
                            physical memory
    --verbosity INT         0 silent, 1 summary (default), 2 stage timing, 3 debug
```

## Evaluation order

Edits apply in a fixed order: header replacement → contig rename → ID operations → sample renames → header metadata adds/removes. `%CHROM` in ID operations uses the already-renamed contig name.

## Conflicts

These combinations are usage errors:

- `--set-id` with `--update-name` (the template rewrites every ID, so the mapping would never match);
- `--clear-ids` with any other ID option;
- `--set-id`/`--set-missing-var-ids` with `--idspace-to` (put the separator in the template);
- no edit option at all.

## Output details

The header gains `##reed_version`, `##reed_command` and `##reed_rename=<planes touched>` at the end (`--no-version`/`--no-command` turn them off). The closing summary line gives record count, contig count, sample count, changed ID count and elapsed time.

## Examples

```bash
# sample rename + contig rename + fill missing IDs in one command
reed rename in.vcf.gz -s "s0 sampleA s1 sampleB" \
    --rename-chrs chrs.txt --set-missing-var-ids '%CHROM:%POS' \
    -Oz -o renamed.vcf.gz -W

# batch sample rename from a file
reed rename in.vcf.gz -N sample_map.tsv -Oz -o renamed.vcf.gz

# strip a sequencing-batch prefix
reed rename in.vcf.gz --replace-sample-name 'Batch3_' -Oz -o clean_names.vcf.gz

# remove colons from IDs (some downstream tools reject them) + deduplicate
reed rename in.vcf.gz --idspace-to '_' --rm-dup-ids -Oz -o ids.vcf.gz

# rewrite every ID as chrom:pos:ref:alt
reed rename in.vcf.gz --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ids.vcf.gz

# swap the header, add contig lengths, drop old reference annotations
reed rename in.vcf.gz -h new_header.txt --remove-meta INFO -Oz -o newh.vcf.gz
reed rename in.vcf.gz -f ref.fa.fai -Oz -o contigs.vcf.gz

# normalize contig names to UCSC style
reed rename in.vcf.gz --chr-style ucsc -Oz -o ucsc.vcf.gz
```
