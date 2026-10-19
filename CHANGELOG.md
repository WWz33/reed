# Changelog

All notable changes to reed are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the
lockstep workspace versions, and the CLI flags plus the output file formats
are the stability promise (`docs/spec/cli.md` §0).

## [Unreleased] — doc_v15 缺陷修复（2026-10-18）

对一份全库评审（6 子代理静态评审 + 逐条实证复核，40 条发现）的修复；
规范与证据见 `docs/rounds/v15/`。

### Fixed

- **filter**：Float 类型 FORMAT 字段的基因型掩码不再整体失效（D15-1）；
  掩码与 `--set-gts .` 改变 cohort 后 `-q/-Q` 一律按掩码后基因型重算，
  不再读陈旧 INFO/AC,AN（D15-2，cli.md §3.2）；`--proj` 超过群体染色体
  数、pop 文件把同一样本分配进两个群体、`--levels` 空列表均报用法错并
  点名（D15-9/10/21）；`--bed-file` 与 `--step` 解析期互斥、`-o -` 写
  stdout 而非字面文件 `-`（D15-22/39）。
- **popg**：缺少 `-o` 在打开输入前即报错（D15-13）；stdout 管道关闭由
  panic 改为普通运行错误（D15-14）；裸 popg 不再忽略
  `--variant-only/--bypass-invariant-check`，不变位点检查按窗口路由同款
  执行（D15-15）；`--f2/--ld-decay` 不再静默把 π/dxy/θ 切到 variant-only
  口径（D15-18）；窗口发射改为"严格更后位点才 emit"，跨 BGZF 块接缝的**合法同 POS 重复**（VCFv4.3："It is permitted to have multiple records with the same POS"）落入仍开窗口，窗口和与 contig 行逐项一致（D15-19）；
  尾随三态记录之后的窗口正常输出（D15-20）；`--pbs` **可重复**：单三元组写 `PREFIX_pbs.txt`（逐字节不变），多三元组每份
  `PREFIX_pbs.<A>-<B>-<C>.txt` 各持自身表头（D15-8，命名规则见 cli.md §7）。
- **convert**：`--phylip` 组合运行不再吞掉 `--012/--hapmap/--msa` 的位点
  （D15-4）；`--012` 双等位门修正，三态位点按 vcftools 语义跳过（D15-5）；
  仅 `--nexus-binary` 时 `--outgroup` 生效且未知名报错（D15-6）；MSA 长
  等位不再以 u8 截断（D15-26）；`POS=0` 按 malformed 跳过（D15-27）；
  gzip hapmap 可走 stdin（D15-29）。
- **rename**：`--rename-chrs` 撞名报错而非写出重复 `##contig` 的非法
  VCF（D15-25）；`--idspace-to` 与 `--set-missing-var-ids` 互斥（D15-30）；
  `--rm-dup-ids` 改为两阶段（串行定序 + 批次并行 deflate：worker 不再做
  被丢弃的压缩、索引键随记录并行算好，有序 sink 按字节封口批次并行压缩）：
  3M 记录 `-Oz` 1.39s → 0.23s（对照 15.4× → 2.5×），压缩流与 sidecar 对
  任意 `--threads` 逐字节相同，decompressed 输出与改前逐字节一致（D15-28）。
- **engine**：输入源线程的 panic 与 worker/reducer 一样收敛到单一错误
  出口，报 `reader panicked: …`（`EngineError::ReaderPanicked`；此前该闭包
  不在 `catch_unwind` 下，panic 会绕过引擎的错误通道，D15-40）。
- **io**：BCF 头恰落 BGZF 块边界（头文本长度 ≡ 0 mod 65280）不再整体
  错位一字节（D15-3）；`##contig` 行的显式 `IDX=` 按 htslib 规则生效
  （D15-24）。
- **stats**：Hudson Fst/PBS 在窗口内无可用分子时输出 0.0 而非 NA
  （scikit-allel `nansum` 口径，D15-11）；oracle 的 PBS 与 f2/f3/f4 改为
  与 shipped/allel 同口径（D15-12/31）。
- **性能**：`--indep-pairwise` 的保留位点查询由逐 contig 线性扫描改为
  二分（1 万 contig、30 万位点单线程 ≈11×，D15-17）；`--rm-dup-ids` 的
  去重输出由"sink 串行重压"改为两阶段并行（3M 记录 `-Oz` 1.39s →
  0.23s，D15-28）；GT 浮点渲染 `format_g6_into` 改为栈缓冲零分配
  （242.9 → 182.9 ns/次，D15-35）；`sketch_size_bytes`
  补上 mask 位平面行，`--mind`+`--indep-pairwise` 的 sketch 分支加预算
  检查（D15-16，预防性）。

### Changed

- `reed popg` 不再接受 `--tmp-dir`（该命令不落临时文件）；`reed stats`
  保留该参数并在帮助中写明"单趟报告从不落盘，此参数仅为兼容接受"
  （doc_v15 D15-23，cli.md §2 同步）。
- 内部收敛：Patterson 只保留校正式（`PattersonKind::Naive` 与相关死结构删除，
  `f2_site`/`f3_site` 去掉 kind 参数）；QC summary 只剩一个实现
  （`summary::SummaryAcc` 再导出 `report` 的活结构，`QcReportAcc::merge` 委托）；
  fst/PBS 审计改用 `PairAcc::add_partial`，不再假装跟踪 missing（D15-38）；
  `notify_abort` 持锁通知关闭丢失唤醒窗口（D15-32）；`--indep-pairwise` 的
  保留位点查询与 convert/popg 的 `--proj` 校验收口到共享层。

### Changed

- `convert` 的输出缓冲改为两级落盘（doc_v15 D15-7）：单缓冲超 4 MiB 或
  sink 缓冲总量超 8 MiB 即落临时文件、渲染时流回；覆盖 hapmap 行、
  `--012` 矩阵与 `.pos`、phylo 列与 `used_sites`、MSA 帧流（渲染改为逐样本遍历，
  并按样本 seek 只读该样本切片——不做整帧重读，延续 v13 D13-2 的禁止重复 IO）。
  800k 位点/50 样本夹具
  单线程峰值内存：hapmap 123.4 → 15.6 MB、`--012` 70.8 → 26.0 MB、MSA 151.4 → 15.8 MB，
  输出逐字节不变；临时文件在成功、失败与中止路径均清理。MSA 的
  per-contig 缓冲仍驻留内存（留待单独设计）。

### Removed

- 死代码：`JointSfsAcc`（连同 `tests/sfs_joint.rs`）、`PattersonKind::Naive`、
  `Cursor::skip`、`SiteView::info_int/info_float`、`HeaderScan::n_samples`、
  `WindowGrid::bin_of_pos0`、`keep_mode` 死参、`header_level` 弃值、stats sink 的
  死 `or(take)`、convert 的空臂（doc_v15 D15-36）。

## [0.4.0] - 2026-10-17

The doc_v14 command-surface round: `reed stats` is split into the QC command
(`stats`) and the windowed statistics command (`popg`); `--hist` and
`reed rename` are new; `reed --help` groups its commands seqkit-style. The
split is byte-identical for every existing output
(`gold/popgen/run_split_parity.sh`, 22 cases).

### Changed

- **`reed stats` is now the QC-report command** (`--freq/--hardy/--het/
  --missing/--depth/--summary`); the windowed population-genetics statistics
  moved to **`reed popg`** (alias `popgen`; doc_v14 §2). The moved flags
  (`--pi --dxy --theta --tajima-d --fst --fst-global --pbs --pbs-normed
  --f2 --f3 --f4 --dstat --f4-ratio --block-size --ld-decay --ld-pairs
  --sfs --proj --sfs-int --variant-only --bypass-invariant-check
  --contig-lengths --window --step --levels --skip-empty --bed-file`) on
  `reed stats` are parse-time errors naming the exact rewrite (doc_v14
  D14-2) — never silent aliases.
- A bare `reed stats` (no report flag) is a usage error naming the flags.
- The library facade renamed its drivers to follow the commands:
  `run_stats` → `run_popgen` (`StatsRequest` → `PopgenRequest`,
  `StatsReport` → `PopgenReport`), `run_reports` → `run_stats`
  (`ReportRequest` → `StatsReportRequest`, `ReportsOutcome` →
  `StatsOutcome`). The library has no stability promise (cli.md §0.9).
- `reed --help` lists the commands in groups (Filtering, Statistics,
  Transform, File Operations, Miscellaneous); command paths are unchanged
  (doc_v14 §1).

### Added

- Expressions (`-i/-e`) accept negative number literals (`MQRankSum<-12.5`,
  the GATK hard-filter shape; bcftools' unary minus, `filter.c:3902`). A
  `-` before a digit starts the number; a leading `+` stays an error, as in
  bcftools, which lexes it as the addition operator (`filter.c:332`); `-`
  inside identifiers stays part of the name, and there is no arithmetic, so
  `A - B` remains an error.
- `reed popg`: the windowed statistics (π, dxy, Fst, θ, Tajima's D, PBS,
  SFS, Patterson f-statistics, LD decay), with the 0.3.x flags and
  byte-identical outputs.
- `reed stats --hist LIST`: histograms from the same pass (doc_v14 §3) —
  built-ins `qual`, `dp`, `dp-site`, `gq`, `idl`, `af`, plus generic
  `INFO/<TAG>` and `FMT/<TAG>` requests; one `PREFIX.<m>.hist` TSV per
  metric with exact mean/median/p95 for integer metrics.
- `reed stats --singletons`: the per-allele singleton/private-doubleton
  table in vcftools' layout (`CHROM POS SINGLETON/DOUBLETON ALLELE INDV`,
  `PREFIX.singletons`), with a POP column under `--pop`.
- `reed stats --pop FILE` (with `--force-samples`): stratified QC
  (doc_v14 §2.2) — `--freq`/`--missing` write `PREFIX.afreq.pop` and
  `PREFIX.vmiss.pop` (one row per population, first-appearance order),
  `--singletons` gains a POP column (`.` for a carrier outside every
  group; counting still runs over the whole cohort). The stratified values
  equal a per-population plink2 `--freq/--missing` run
  (`gold/popgen/run_popqc.sh`); the cohort files are unchanged. `--pop`
  without one of `--freq/--missing/--singletons` is a usage error.
  Multiallelic ALT frequencies stay collapsed (the `.afreq` convention).
- `reed rename` (alias `reheader`): identifier edits without decoding
  genotypes (doc_v14 §4) — samples (`-s/-N/--replace-sample-name`), contigs
  (`--rename-chrs`, `--chr-style plain|ucsc|mt-is-m`), variant IDs
  (`--set-id`, `--set-missing-var-ids`, `--update-name`, `--idspace-to`,
  `--rm-dup-ids`, `--clear-ids`), header lines (`-h -f --add-meta
  --remove-meta`); VCF byte editing with `-Ov/-Oz`.

### Fixed

- The `-q/-Q` stale-INFO warning fired on every run whose header declares
  `INFO/AC` and `INFO/AN` with an unchanged cohort (every GATK header does),
  asserting "-q/-Q reads INFO/AC,INFO/AN" of commands that never passed
  `-q/-Q`. It now fires only when a `nref`-form threshold actually takes
  the INFO path (spec cli.md §3.2: 取到时).
- `reed rename` wrote a blank line at every chunk boundary (the seam record
  arrives with its terminator and got a second one), so every output over
  one chunk carried stray empty lines; with `--rm-dup-ids` the same doubled
  newline made the dedup pass fail outright (`job failed: dedup: unparsable
  record`). The suffix probe no longer restarts at `_2` for every duplicate,
  so a run of identical IDs is linear instead of quadratic.
- `reed rename` ran its record planes on inputs whose CHROM is not declared
  in `##contig`: a contig-less VCF (legal VCF) and a sites-only 8-column VCF
  both failed mid-run. Plain output now passes them through, `-W` says what
  the header is missing, and the record planes are refused up front on
  sites-only input.
- `reed rename -h FILE -W` built the index from the *input* header's contig
  order while writing the replacement header, so region queries could
  resolve to the wrong contig; the index now follows the written header.
- `reed rename -f FILE` appended the missing `##contig` lines after the
  `#CHROM` line (not a valid VCF header); they are inserted before it.
- `--rename-chrs` and `--chr-style` overwrote each other instead of
  stacking (the style was applied to the pre-rename name and discarded the
  explicit rename).
- A CRLF header's last sample kept its `\r`, so sample renames counted the
  edit but silently did not apply it.
- `--replace-sample-name` used to accept regex syntax and silently do
  nothing (or strip the wrong text); the flag is documented as literal, so
  metacharacters are now a usage error.
- `--update-name` did not check that a mapped name exists (the one
  validation-table row that was missing): a typo renamed nothing and exited
  0; the run now fails naming the key.
- `--no-command` was parsed and ignored by `reed rename` (the
  `##reed_command` line was written anyway).
- `--hist` float tags were labelled one 0.1 step high (`iqual` divided by
  10 instead of `(iqual−1)/10`), so `SC=1.0` printed as the bin `1.1 1.2`;
  a missing QUAL became a real `0.0` bin instead of `#missing`; the `af`
  bins used the wrong denominator and dropped indels; integer metrics
  passed the overflow bucket off as an exact value (`DP 600` → row `501`,
  mean 501). Bins are now labelled in the quantity's own units, capped
  observations print `>max` and stay out of the exact quantiles, and only
  occupied bins are written (a single huge value used to materialize a
  half-million-row file).
- `--hist gq`/`dp`/`dp-site` did not check that the header declares the
  field and wrote an all-zero histogram with exit 0; an all-empty
  `--hist ""` list did the same. Both are usage errors now.
- `--hist FMT/<TAG>` counted one missing observation per record when the
  FORMAT key was absent (it is one per sample) and never split
  Number=A/G/R subfields, so `FMT/AD` looked empty.
- `StatsOutcome::hist` (the library's histogram output) was always empty.
- `reed popg --ld-pairs` without `--ld-decay` exited 0 and produced no
  output; it is a usage error now.
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
  `--sfs`, byte-identical to `reed popg --sfs` (one shared request builder,
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
  `reed popg`: a post-mask site rule over the genotype classes
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
