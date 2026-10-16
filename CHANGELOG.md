# Changelog

All notable changes to reed are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the
lockstep workspace versions, and the CLI flags plus the output file formats
are the stability promise (`docs/spec/cli.md` §0).

## [Unreleased]

### Fixed

- INFO values survive the BCF round trip: an Integer outside the 32-bit
  range becomes missing instead of wrapping (`3000000000` used to be
  written as `-1294967296`; htslib and VCFv4.5 §1.3 both forbid the value),
  a Flag carrying a value keeps it (`DB=1` lost it, and `DB=0` — bcftools'
  "clear the flag" — came back as a set flag), and a String value of `.` is
  no longer turned into an empty string (doc_v12 D12-3/4/5).
- `-q/-Q` follow the genotypes once the sample set changed: after a static
  subset, a mask or `--mind`, the input's `INFO/AC,AN` describe a different
  cohort, so the decision is recomputed from the kept genotypes like
  bcftools does (vcfview.c:341-350); the promised stale-INFO warning is
  emitted on the INFO path (doc_v12 D12-6).
- `--mind` no longer evaluates a tail `-i/-e` expression blind: the sketch
  path has no record bytes, so `--mind` with such an expression takes the
  two-pass path, and `--mind --indep-pairwise` with one is a usage error
  instead of silently dropping every record (doc_v12 D12-1).
- Expressions see the header's tag types: `DB=1`/`DB=0` test a flag's
  presence (they used to match nothing), a multi-valued INFO tag compares
  any element (`DP4=0,0,30,40` at `DP4>10` used to miss), and a lone tag
  with no comparison is a compile-time error instead of a silent
  constant (bcftools evaluates it to false; doc_v12 D12-7).
- One INFO number grammar for the filter and the writer: `XI=1e3` no longer
  filters as 1000 while writing `.` to BCF, and an unparsable numeric tag
  is missing rather than compared byte-wise (`XI=abc` used to pass
  `XI>4`; doc_v12 D12-14).
- A GT-less record keeps its input `AC/AN` under a sample subset (it used
  to be recomputed from nothing as `AC=0;AN=0`; bcftools' `bcf_get_fmt`
  null check, vcfview.c:343; doc_v12 D12-2).
- `ALT='.'` counts zero ALT alleles, so `-m 2`/`-M 1`/`N_ALT=0` match
  bcftools on invariant sites; `-M 1` means "invariant only" instead of
  being dropped as "off" (doc_v12 D12-10).
- Records whose INFO/FORMAT counts do not fit the BCF header fields are
  refused instead of packed with the overflow leaking into the next field
  (65536 INFO keys or 300 FORMAT keys produced an unreadable file;
  doc_v12 D12-9).
- `--fill-tags` validates its list: an unsupported tag (AF, F_MISSING, a
  typo) is a usage error naming the supported set, like `bcftools
  +fill-tags`; it used to be a silent no-op that the default `-x` then
  turned into a deletion (doc_v12 D12-11).
- `[[site.info]]`'s `quantifier = "any"` requires an element to satisfy the
  comparison (it was true whenever an element parsed; doc_v12 D12-8).
- An INFO `AC` vector whose length does not match the ALT count falls back
  to the genotype counts (htslib's `bcf_calc_ac` check, vcfutils.c:56-64;
  doc_v12 D12-13).

### Added

- `reed convert` (docs/rounds/v13): one parallel pass over the input, any
  combination of targets. Phylogenetic matrices with vcf2phylip 2.9
  semantics (Ortiz 2019, DOI 10.5281/zenodo.2540861) — `--phylip`,
  `--fasta`, `--nexus`, `--nexus-binary` (SNAPP), plus `--outgroup`,
  `--min-samples-locus`, `--resolve-iupac` (deterministic, seed 0),
  `--write-used-sites`; the vcftools `--012` triple; the TASSEL hapmap
  (`.hmp.txt`); a reference-driven whole-site MSA with vcf2msa semantics
  (`--msa --reference`, `--region`, `--flank`, `--collapse-suffix`); and
  `--sfs`, byte-identical to `reed stats --sfs` (one shared request builder,
  one SFS fold). Outputs are independent of `--threads`; parity is gated by
  `gold/popgen/run_convert_parity.sh` (vcf2phylip five files, vcftools
  `--012`, TASSEL read-back, vcf2msa sequences). Deliberate divergences from
  the reference tools are registered in `docs/spec/cli.md` §8b: malformed
  records are skipped and counted, `--012` writes `-1` for half calls and
  haploids (vcftools writes 1), the hapmap `alleles` column lists the
  declared REF/ALT alleles, the MSA pads indels deterministically instead of
  calling an external aligner, and an allele set outside the IUPAC table
  skips the record where vcf2phylip crashes.
- `reed-core`'s `iupac` module: the IUPAC ambiguity table (62 entries), the
  binary-NEXUS genotype table, and the inverse lookup (`ambig_expand`: `R` →
  `AG`), the single implementation shared by every convert target.
- `reed convert --vcf`: the hapmap → VCF reverse direction (v13 C4), so the
  two formats convert both ways. The input is a hapmap (`.hmp.txt`, one- or
  two-letter cells auto-detected); `chrom`/`pos`/`rs#`/`alleles` are copied
  verbatim (`ID` = `rs#`, REF/ALT in declared order), each IUPAC cell expands
  to a `GT` (`R` → `0/1`), only `N`/`?` read as missing, and `alleles=NA`
  becomes an all-missing record. `-O v` (default) or `-Oz[level]` writes
  plain text or BGZF, `-o` omitted streams to stdout; the header carries no
  `##contig` (TASSEL writes none) and its `##reed_*` provenance lines are
  controlled by `--no-version`/`--no-command`. Exclusive with every
  VCF-input target and with sample subsetting; malformed rows are skipped and
  counted. On the lossless domain (biallelic diploid SNPs, no missing) a VCF
  → hapmap → VCF round trip reproduces every record, and TASSEL 5 reads
  reed's VCF back (`gold/popgen/run_convert_parity.sh` gate 6).
- `-g/--genotype [^](hom|het|miss)` (bcftools `-g`) for `reed filter` and
  `reed stats`: a post-mask site rule over the genotype classes
  (`TailRules::gt_class` → `TailRule::GtClass`). A half-call (`0/.`) reads
  as missing (bcftools `bcf_gt_type` UNKN), so `SiteCounts` gains
  `n_half_called`, counted by `count_site`, the plane sketch and the
  single-pass pipeline. Repeated flags AND (bcftools keeps only the last
  one); exact duplicates collapse and opposite polarities (`hom` + `^hom`)
  are a compile-time error. Funnel ids are `after_mask.gt_class.(hom|het|
  miss|^hom|^het|^miss)`. Diploid parity against bcftools 1.24 PASS over
  every call shape; haploid sites stay policy-gated (`require_diploid`
  rejects, `allow_haploid` counts them as half calls, so `-g hom`/`^miss`
  diverge from bcftools there by the same convention as `--maf`/`--mac`).
  A record whose FORMAT lacks the GT field skips every `-g` rule and
  survives (`SiteCounts::has_gt`; bcftools' `bcf_get_fmt` null check does
  the same) — distinct from a GT whose values are all `./.`. Class names
  are matched exactly, case-insensitively (bcftools `strcasecmp`; no
  surrounding whitespace is trimmed). Like bcftools, the rules evaluate
  after subsetting, over the kept samples only.

### Changed

- WC84 Fst (`--fst wc84` pairwise and `--fst-global`): a site where a
  population has no called alleles now contributes only the finite
  Weir-Cockerham `c` component (h̄/2 over the typed populations) instead of
  being skipped — pixy / scikit-allel `nansum` semantics, for both the
  pairwise and the r-population forms (doc_v7 D7-20). vcftools and plink2
  drop such sites (previously reed's behavior); hierfstat recomputes them
  with the typed populations only. Parity: pairwise == pixy to 13
  significant digits on a masked Chr01 with 1–3% all-missing sites
  (0.0086203 / 0.0020595 / −0.0005171), global == scikit-allel r=3
  (0.003237); the vcftools gold gate still passes (the benchmark data has
  no such sites, and values there are bit-identical).

- The engine's abort signal is an atomic flag (doc_v6 A6-1, Y5): `is_aborted`
  no longer flips back when a `select!` waiter consumes a channel token, the
  wake-up channel is sized by `--threads` (which may exceed the core count),
  and registered byte-permit pools are notified on trigger, so
  `BytePermits::acquire` waits on a condvar instead of polling every 10 ms
  (`docs/spec/architecture.md` §7.5). Every engine entry point now builds its
  abort channel from the run's thread budget (`abort_for`).
- In-flight chunk permits live in a drop guard (A6-1, Y5): a worker failure
  or panic, an aborted send, or a dropped pending message all return the
  permit and the in-flight accounting, so the reader can no longer block
  forever on a leaked permit. `StageMetrics::in_flight_final_bytes` makes the
  residue assertable.
- `threads <= 2` now takes a real sequential path (A6-3, Y5): one thread
  reads, inflates, computes and writes, with the same seam stitcher and sink
  as the parallel path; output is byte-identical (locked by a test). The
  `is_sequential` selector documented in §7.3 finally has a caller.
- The seam stitcher is bounded (A6-2, Y5): a record reassembled across chunk
  seams may not exceed `reed_engine::pipeline::MAX_RECORD_BYTES` (1 GiB,
  shared with the BCF frame chain), so a newline-free input fails the run
  instead of growing the reducer's fragment to the file size.
- `StageMetrics` gains `reader_wait`, `reducer_wait`, `inflate_time`,
  `job_time`, `deflate_time` and `in_flight_final_bytes` (A6-4, Y5). Output
  deflate time rides `ChunkOutput` to the sink (`OrderedSink::deflate_time`),
  so the three worker stages are measured instead of inferred; `--timing`
  and `--verbosity 2` print `reader_wait_ms`, `reduce_wait_ms`, `inflate_ms`,
  `job_ms` and `deflate_ms`.

- One pass0 serves every filter-side global analysis (doc_v6 03 §7, Y1e):
  the `MatrixJob`/`MatrixSink` pair fills per-contig site-major plane
  matrices that `--mind`, `--indep-pairwise` and the sketch path's mid
  stage all read back. The three previous normalizations
  (`MissingJob`, `SketchJob`, `LdPruneJob`) and `reed_filter::sketch` are
  deleted. `--mind` counts come from the miss-plane column popcounts
  (`called(s) = 2·sites − pop(slot0_miss ∨ mask) − pop(slot1_miss ∨ mask)`);
  `--mind` + `--indep-pairwise` now combine in one run (they previously
  could not share the first open).
- The pruning runs after every site rule: without `--mind` the pass0
  evaluates the tail rules inline and the pruner sees only post-tail
  sites (plink2's load-then-prune order, 02-cli §6 step 6b); with `--mind`
  the mid stage's survivor bitmap filters the pruner. The report's
  `ld_prune` counts are post-tail (`n_seen` shrinks when tail rules drop
  sites).
- `PipelineOptions::scan_head_failures` and `RecordView::head_passed` are
  gone (their only user was the sketch pass).
- The window-layout output writers are `reed::windows` (`WindowLayoutOptions`,
  `window_kinds`) instead of the previous vendor-named module (doc_v6 05 §7.2,
  E6-13); the output bytes are unchanged.
- The unsound AVX-512 counting kernel and `CountTable` were deleted
  (doc_v6 03 §5, A6-22); the mmap is now the workspace's only audited
  `unsafe` site.
- `ContigId` is a single newtype shared by every layer (doc_v6 A6-23).
- The whole tree is `cargo fmt`-clean under the pinned 1.99 toolchain, and the
  CI clippy command (`--all-targets -- -D warnings`) passes (doc_v6 05 §1–§2).

- I/O hot paths reuse their buffers (doc_v6 03 §3, A6-5/6/8, Y2): the
  deflater keeps a resident output scratch instead of zero-filling the
  worst-case bound per block, `inflate_chunk` walks the chunk's block table
  into exact-size windows (no second header parse, no per-block 64 KiB
  memset), and the chunk tail scan uses `memrchr`. The filter gold run
  drops from 460 to 400 CPU-s.
- The streaming BGZF source keeps empty members in the chunk (A6-7, Y2):
  dropping them used to shift every later chunk's `compressed_start`, so
  virtual offsets computed from stdin/pipe input could be wrong. The
  streaming and mmap sources now produce identical chunk metadata, locked
  by a test with an empty-member fixture.
- Per-record allocations left the hot paths (A6-11/12, Y2): the error
  `Locus` is built from borrowed bytes only when an error carries it
  (`LocusSource`), INFO edits write into the output buffer directly, and
  the index writer no longer builds a contig-name `String` per record on
  the (default) CSI path. The BCF GT decode table is one `Arc` per run
  instead of one per worker (A6-10).
- The filter's survivor emit is one method (`emit_survivor`) for the text
  and native BCF paths, the stale-tag list builds once per chunk, and the
  mask rules run columnar — one contiguous pass per rule over the scanned
  FORMAT columns, integer thresholds compared in the integer domain
  (A6-15/17/20, Y3). The `--set-gts .` expression path reuses the
  pipeline's scanned FORMAT table in place (no per-sample allocation) and
  expressions compile to operator enums at parse time (A6-18). INFO reads
  go through a per-record fragment index shared by the typed head rules
  and the expressions (`reed_io::vcf::info::InfoIndex`, A6-18). The
  exact-test HWE (`--fill-tags HWE`, `--hwe`, `--hardy`) is memoized per
  worker through `HweCache` (A6-16); the emitted p-values are unchanged.
  Full-gold CPU moved 398.2 → 395.6 s; the measurements behind the
  roadmap's larger estimate are recorded in `docs/rounds/v6/STATUS.md`
  (Y3 deviation).

### Added

- Streaming per-contig pruning (Y1f): each contig's matrix feeds its own
  prune thread as soon as the stream closes it, so one matrix is resident
  at a time; a single contig over the `--memory` budget fails with
  `MemBudgetExceeded`.

- `LICENSE-MIT` and `LICENSE-APACHE` (the manifests declared
  `MIT OR Apache-2.0` with no license text until now).
- `reed --version` now prints the source state the binary was built from
  (`reed 0.3.0 (v0.3.0-12-g8709f1f 2026-10-10, rustc 1.99.0)`), injected by
  `build.rs` (E6-11).
- README sections: Installation, Citation, Input requirements, Exit codes
  (E6-12).

## [0.3.0] - 2026-10-10

The current development line (the W0–W6 / X0–X7 rounds in `docs/rounds/`).

### Added

- `reed stats`: π, dxy, Fst (`wc84`/`hudson`, pairwise and global), Watterson
  θ, Tajima's D, PBS, projected SFS, f-statistics (`--f2/--f3/--f4/--dstat/--f4-ratio`),
  LD (`--ld-decay`, `--ld-pairs`), per-site/per-sample reports
  (`--freq/--hardy/--het/--missing/--depth`) and the bcftools-stats
  compatible `--summary`, all in one pass.
- `reed filter`: bcftools-view-style hard filtering, genotype masking,
  INFO recomputation (`--fill-tags`), `-q/-Q` with the `:type` suffixes,
  LD pruning (`--indep-pairwise`) as a site rule, index writing (`-W`/`--write-index`).
- `reed index`: CSI (default) and TBI, parallel.
- Window layout output (`PREFIX_pi.txt`, `PREFIX_dxy.txt`, …) verified
  column-for-column against pixy 2.2.4; `docs/rounds/v5/05-algorithms.md`.
- BCF read and write (native, no text round-trip for the filter decision),
  BGZF chunk streaming with mmap and streaming sources.
- The deterministic chunk engine: ticket window, ordered reducer, byte-budget
  backpressure; output byte-identical across thread counts and strategies.

### Verified against

bcftools, plink2, vcftools, pixy, scikit-allel, easySFS, admixtools-style
f-statistics; see `gold/` and the `Verification` section of the README.

## [0.2.0]

The pre-migration architecture (`docs/rounds/v1/design.md`), frozen as the
differential-testing oracle (`crates/oracle`). Not a published release; the
history is recorded in `docs/rounds/v2/STATUS.md`.
