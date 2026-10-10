# Changelog

All notable changes to reed are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the
lockstep workspace versions, and the CLI flags plus the output file formats
are the stability promise (`docs/spec/cli.md` §0).

## [Unreleased]

### Changed

- The window-layout output writers are `reed::windows` (`WindowLayoutOptions`,
  `window_kinds`) instead of the previous vendor-named module (doc_v6 05 §7.2,
  E6-13); the output bytes are unchanged.
- The unsound AVX-512 counting kernel and `CountTable` were deleted
  (doc_v6 03 §5, A6-22); the mmap is now the workspace's only audited
  `unsafe` site.
- `ContigId` is a single newtype shared by every layer (doc_v6 A6-23).
- The whole tree is `cargo fmt`-clean under the pinned 1.99 toolchain, and the
  CI clippy command (`--all-targets -- -D warnings`) passes (doc_v6 05 §1–§2).

### Added

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
