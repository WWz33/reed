---
title: index
---

## Usage

```
reed index [OPTIONS] <FILE>
```

Build a CSI or TBI index for a VCF/BCF, CSI by default. `filter -r/-R` region queries and any downstream tool doing random access by region rely on an index next to the input.

## Options

```
<FILE>                      Input file, must be BGZF-compressed VCF or BCF.
                            stdin is rejected: the index records byte offsets
                            and reads the file back
-c, --csi                   Write a CSI index (default; exclusive with -t)
-t, --tbi                   Write a TBI index
-o, --output FILE           Write the index elsewhere (default: next to the
                            input, with a .csi/.tbi suffix)
-f, --force                 Overwrite an existing index; without it an existing
                            index is an error with a hint
    --threads INT           Worker threads, default all available cores
    --verbosity INT         0 silent, 1 summary (default), 2 stage timing, 3 debug
```

## Behavior

- The input must be a file. Plain gzip or uncompressed VCF is an error: indexing them yields all-zero virtual offsets.
- BCF cannot take TBI; BCF is CSI only.
- TBI coordinates top out at 2²⁹−1; larger positions need CSI.
- Sorted input is required: a position going backwards within a contig raises `positions are not ascending (indexing needs sorted input)`.
- The CSI bin structure and chunk depth use the standard parameters, and the product can be read by any program that reads CSI.

## Examples

```bash
reed index raw.vcf.gz                  # raw.vcf.gz.csi
reed index -t raw.vcf.gz               # raw.vcf.gz.tbi
reed index -f raw.vcf.gz               # overwrite the old index
reed index -o /fast/disk/raw.csi raw.vcf.gz
```
