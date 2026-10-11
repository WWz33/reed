# Changelog

All notable changes to reed are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the
lockstep workspace versions, and the CLI flags plus the output file formats
are the stability promise (`docs/spec/cli.md` §0).

## [Unreleased]

### Changed

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
