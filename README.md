# reed

[English](README.md) | [简体中文](README.zh-CN.md)

**reed** = population-genetics statistics library + hard-filter engine, in pure Rust.

Two independent crates under the same umbrella:

1. **reed-core**: windowed population-genetics statistics (π, dxy, Fst, PBS, SFS, Tajima's D, LD r², HWE, f-statistics). See `reed-rs/crates/reed-core`.
2. **reed-filter**: streaming hard-filter engine for VCF/BCF (sample missing-rate reduction → filter funnel → genotype mask → recomputation → emit). See `reed-rs/crates/reed-filter`.

Both share `reed-io` as their transport layer (VCF/BCF record streamer based on noodles).

---

## 1) reed-core: windowed statistics (v0.2)

Compute site/window statistics from VCF/BCF with correct ascertainment semantics and exact integer tallies.

**Supported statistics**

| Category | Statistics | Estimators / notes |
|---|---|---|
| Diversity | π | PerPairComparable / AccessibleBp; empty/all-missing windows return None |
| | Watterson θ | per-site θ_W + total; empty windows → None; variance shared with Tajima |
| | Tajima's D | min_sites guard; uses Watterson θ total & per-bp variants |
| Differentiation | dxy | two estimators matching allele-level deletion vs accessible bp |
| | Fst | Hudson (2013), Weir-Cockerham (1984), Nei (1986); corrected implementations (no pooled ratios) |
| | PBS | Yi (2010) / PBSn1; pair-audit trail; complete intersection semantics documented |
| Allele frequency | AC, AF, SFS | folded/unfolded spectra; sample-size corrected allele counts |
| LD & HWE | LD r² (Rogers-Huff) | pairwise; optional dosage integerization (+1.2~2x speed, deterministic accuracy) |
| | HWE exact test | Wigginton (2005) Levene–Haldane recursion; double vs rational arithmetic verified |
| | Heterozygosity | Hobs, Hexp, Fis; p-denominator fixed (alt/(ref+alt)) |
| Admixture | f₂/f₃/f₄, D | sample-mean definitions; Patterson correction reserved (v0.3) |

**Design guarantees**

- Count-first semantics: every float result has integer numerator/denominator available for audit.
- Ascertainment awareness: SiteCompleteness model (WithInvariants / VariantOnly(Filtered/Complete/Unknown)); Watterson θ rejects VO(Unknown), others bias-marked or exempted where justified.
- Correct error paths: missing contig lengths, non-diploid calls, multiallelic sites, out-of-order records, positions beyond contig length — all explicit errors (no silent drops).
- No per-sample random access: series iteration (column-major, linear) over genotypes; O(n²) indexing is banned by API convention.

**Running examples** (concepts only; full code in `reed-rs/crates/reed-core/tests/*`)

```rust
use reed_core::{ReadOptions, WindowBuffer, PopulationView};

let opts = ReadOptions::builder()
    .window_size(1_000)
    .step(1_000)
    // contig_len required; Beagle outputs stripped ##contig are rejected
    .build()?; 

for window in WindowBuffer::from_file("clean.bcf", opts, &header)? {
    let fv = PopulationView::from_groups(...).validate_against(&window)?;
    
    let pi = reed_core::estimators::pi::compute(window.as_ref(), &fv); 
    let fst_ab = reed_core::estimators::fst::hudson(window.as_ref(), &fv, group_a, group_b)?;
}
```

See `reed-rs/docs/filter-architecture.md` §12 & `DESIGN.md` for implementation status & roadmap.

---

## 2) reed-filter: hard-filter engine (v0.2)

Parallel, field-preserving filtering with gold-parity guarantees (bit-for-bit aligned with bcftools on real data).

**Key features**

- Two-pass pipeline: pass0 (genome-wide sample missing-rate reduction) → pass1 (filter funnel + edit + emit).
- Exact per-rule funnel counters (RegionRule → biallelic SNP → PASS → QUAL → INFO QD/FS/MQ/SOR → GT mask → site missing → MAF → HWE).
- Genotype mask by DP/GQ thresholds; keeps DP/GQ traceable (masked samples become `./.` while original fields preserved).
- Recompute AC/AN/NS/MAF/HWE post-mask from GT (never use stale INFO).
- Deterministic output across thread counts; cross-thread reproducibility tests included.
- Gold parity: L_site / L_gt / L_tags verification against bcftools/vcftools/plink2 fixtures (gold/).

**FilterSpec** (type-safe, no strings)

```rust
let spec = FilterSpec::builder()
    .cohort(CohortRule::max_missing(0.20))         // sample missing-rate threshold
    .regions(RegionRule { include: ..., exclude: ... }) // targets-style predicates (streaming, index future work)
    .sites(SiteRules {
        biallelic_snps: true,
        require_pass: true,
        min_qual: Some(30.0),
        info: InfoRules { min_qd: 4.0, max_fs: 50.0, min_mq: 40.0, max_sor: 3.0 },
        max_missing: Some(0.10),
        min_maf: Some(0.05),
        min_hwe_p: Some(1e-6),
    })
    .genotypes(GtMask { dp_min: Some(5), dp_max: Some(100), gq_min: Some(20) })
    .recompute(Recompute::AC_AN_NS_MAF_HWE)
    .build()?;
```

**CLI usage**

```bash
reed-filter [-c conf.conf] [-o out.vcf.gz|out.bcf] [--threads N] <input.vcf.gz|input.bcf>
```

- `-c`: filter spec in KEY=VALUE format (matches gold/config/filters.conf).
- `-o`: omit for stdout plain VCF; `.vcf.gz`/`.bcf` auto-dispatch to MT BGZF pool.
- `--threads`: compute pool size + compression workers; no strategy flags.

Run report (funnel) to stderr; exit codes: 0 success, 1 runtime, 2 usage.

**Performance snapshot (v0.2)**

On clean.vcf.gz (6.6GB compressed, 348 samples, 3.3M variants, 56-core server):
- Current batched RecordBuf transport: ≈9.6 minutes (two-pass parse materialized for all records).
- Primary bottleneck: double-parse cost + remote deallocation contention (diagnosis in architecture doc §7).
- Planned lazy-batch transport + survivor-only materialization + MT inflate → target 1–2 minutes.

See `reed-rs/docs/filter-architecture.md` §7 for detailed performance diagnosis & roadmap.

---

## 3) Common ground: reed-io (transport layer)

- Single backend: noodles 0.117 (no htslib).
- Supports VCF/VCF.gz/BCF (`for_each_record`, `record_header`).
- Record view decoding (GT/DP/GQ, INFO QD/FS/MQ/SOR), surgical edits (apply), header manipulation (filter_header, ensure_info_defs).
- Output sinks: stdout/plain-file (text), VcfGz/Bcf (multithreaded BGZF compression).

---

## Quick start

### Library (reed-core)

```rust
use reed_core::{ReadOptions, WindowBuffer, PopulationView, estimators};

let opts = ReadOptions::builder().window_size(1_000).step(1_000).build()?;
for w in WindowBuffer::from_file("clean.bcf", opts, &header)? {
    let stats = estimators::pi::compute(w.as_ref(), &population_view);
    println!("π={:?}", stats.value());
}
```

### Library (reed-filter)

```rust
use reed_filter::{exec, emitter, exec::ExecOptions, FilterSpec, spec};

let spec = FilterSpec::builder()
    .cohort(/* CohortRule::... */ None)
    .sites(spec::SiteRules { max_missing: Some(0.10), min_maf: Some(0.05), ..Default::default() })
    .genotypes(Some(spec::GtMask { dp_min: Some(5), dp_max: Some(100), gq_min: Some(20) }))
    .build()?;

let mut sink = emitter::open_sink(Some(&path), 8)?;
let report = exec::execute(&input, &spec, sink.as_mut(), &ExecOptions::default())?;
sink.finish()?;
// use report.funnel, report.elapsed, report.n_output_records, etc.
```

### CLI (reed stats) — stats v2 (v0.3, M4)

Windowed π/dxy/Fst/PBS/SFS/θ_W/Tajima's D over a BGZF VCF, computed as
additive accumulators over per-site group counts (oracle-verified:
integers exact, floats ≤ 1e-12; byte-identical output across thread
counts):

```bash
reed stats input.vcf.gz --groups samples.3pop.txt --size 100000 --step 50000 \
    --levels window,contig,genome [--skip-empty] [--threads 0] [--out out.tsv]
```

- `--groups`: `sample<TAB>population` rows, first line skipped; groups are
  indexed in first-appearance order.
- `--levels`: window rows (default), per-contig rows, and a final `genome`
  row — all plain accumulator merges, so PBS stays re-aggregable.
- The full per-group/per-pair statistic set is emitted (38 columns for 3
  groups); multiallelic records are skipped and counted.
- Reference throughput: 3.3M records × 348 samples → 19k windows in
  ≈0.9 s (56 threads), ≈23 s single-threaded.

### CLI (reed filter) — filter v2 (v0.3, M5)

Streaming hard filter on the chunk pipeline: head → GT-mask → post-mask
tail → byte-level line edit → worker-local BGZF compression, written in
chunk order with a synchronized tabix index (`--write-index`). Gold-parity
gates (docs/filter-architecture.md §10): site set, genotypes and
recomputed tags (AC/AN/NS integer-exact; MAF/HWE ≤ 1e-9) all match the
bcftools/vcftools/plink2 three-tool gold standard; record output is
byte-identical across thread counts.

```bash
reed filter input.vcf.gz --conf filters.conf --out out.vcf.gz \
    --threads 48 --write-index
reed filter input.vcf.gz --spec filter.toml --out out.vcf.gz   # TOML spec (§11.1)
```

- `--conf`: the gold `KEY=VALUE` compatibility reader (`compat::from_gold_conf`).
- `--spec`: the typed TOML spec (cohort / targets / site / genotype masks /
  site_after_mask / recompute / policy), validated at plan-compile time
  against the input header.
- The output header carries `##reed_version/command/filter_spec/
  filter_sha256/ascertainment` provenance (§11.7–11.8); `reed stats` picks
  the ascertainment up automatically.
- Reference throughput: 3.3M records × 348 samples → 1.77M output records
  in **12.6 s** (TwoPassReparse) / **14.7 s** (SketchThenEmit), 48 threads;
  gates ≤ 120 s (M5) and ≤ 60 s (M5b); cohort reduction + filter +
  recompute + compress + index.
- Strategies (§11.4): `SinglePass` (no cohort rule or keep/remove lists),
  `TwoPassReparse` (max_missing cohort), and `SketchThenEmit` (max_missing
  cohort when the pass0 sketch fits the memory budget: pass0 captures
  per-site genotype planes + per-rule mask bitmaps, the mid stage derives
  tail decisions from `keep ∧ planes`, pass1 re-parses survivors only) —
  all three byte-identical, test-locked; the chosen strategy is auto-picked
  from a first-chunk records/byte extrapolation and reported.

### CLI (reed qc) — QC report (v0.3, M6)

The full §10.1 QC report on the chunk pipeline, with a `bcftools stats`-
compatible text output (feedable to `plot-vcfstats`/MultiQC) plus a
structured JSON report:

```bash
reed qc input.vcf.gz --out report.stats.txt --json report.qc.json
```

Sections: `summary` (SN), `tstv` (TSTV), singletons (SiS), `af_spectrum`
(AF bins, bcftools binning formula), `qual` (QUAL), `substitution` (ST),
`indel_length` (IDD), depth distributions (DP, per-site and per-genotype),
and per-sample counters (PSC) — all additive accumulators over chunk
partials, byte-identical across thread counts.

**Parity gate (M6)**: every bcftools-overlapping section matches
`bcftools stats -s -` line-for-line on the full 3.3M-record gold input
(SN/TSTV/SiS/ST/AF/QUAL/DP/IDD), including the exact `dist_t`
log-binning for high QUAL, the AF singleton-bin fold, and htslib's
INFO/AC+AN preference over GT. Runtime: **3.2 s** (48 threads; gate ≤ 30 s).

### CLI (reed-filter, v0.2 legacy binary — retained until 0.5.0, D14)

```bash
reed-filter -c filters.conf -o out.vcf.gz --threads 8 input.vcf.gz
```

---

## Testing & validation

**reed-core**

- Contracts + families + unequal-groups: bit-for-bit alignment with scikit-allel.
- Unit tests for edge cases: haploid/polyploid rejection, non-biallelic, `.`, out-of-order records, contig revisit, position beyond length.
- clippy 0 warnings.

**reed-filter**

- Determinism across threads (synthetic 10k-record fixture; bytes identical for threads=0/1/3/7).
- Funnel stability tests.
- Sink roundtrip tests (plain VCF file, VCF.gz, BCF).
- Integration with gold/ parity harness (future gate; scripts ready).

Run tests:

```bash
cargo test --release -p reed-core
cargo test --release -p reed-filter
cargo test --release -p reed-io
```

Performance probe: set `REED_PROBE_INPUT=<vcf.gz>` and run `cargo test --release -- --ignored probe_reader_variants` (reed-filter).

---

## Architecture & docs

- `reed-rs/docs/filter-architecture.md`: complete narrative for reed-filter (motivation, model shifts, spec/plan/exec details, v0.2 instrumentation & redesign, roadmap).
- `reed-rs/DESIGN.md`: seed chapter for popgen statistics (ascertainment, availability gates, correctness bugs & fixes, roadmap).
- `reed-rs/docs/naming-conventions.md`: canonical naming across crates (public API surface).

---

## Roadmap & priorities

**reed-core (popgen statistics)**

- P0/P1 completed (v0.2): audits for PBS/Fst/Watterson (missing accounting, boundary guards), metadata fields on results, value() methods, availability doors coverage.
- P2 (v0.2 finished in place): validate population against windows, Result<T,Unavailable> everywhere.
- Future: block jackknife, Patterson correction for f-stats, joint SFS, f4-ratio, individual inbreeding F, LD D/D'.

**reed-filter (hard filter)**

- P0: mechanism confirmation (probe_reader_variants baseline measurement).
- P1: lazy-batch transport (reed-io).
- P2: migrate exec to lazy batches (survivor-only materialization).
- P3: MT BGZF decompression input.
- P4: allocator tuning, buffer pooling (only if profiling shows gain).
- Future: region-index optimizations, ordered-carry (LD pruning), ascertainment v2 manifest integration, WindowFeeder zero-intermediate mode.

---

## License

MIT OR Apache-2.0.

---

## Notes

- Text outputs clamp VCF version ≤4.3 to avoid noodle/GT rendering incompatibilities.
- All reed components use noodles 0.117; no htslib dependency anywhere.
- Pass0 is genome-wide sample reduction and runs strictly before pass1; they are not parallelizable across passes (but parallel inside each pass: decode/compress pools).
- For large-scale runs, keep an eye on memory footprint during in-flight batching; current limits are safe for hundreds of GB input under typical RAM budgets.
