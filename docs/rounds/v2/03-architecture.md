# 03 架构重构（允许"爆破"）

> 目标：一套实现、一个记录解码器、诚实的 crate 边界、真正的背压。
> 性能已经达标，本次重构**不以提速为目的**，以消除重复实现和正确性隐患为目的；
> 每一步都要求 gold 对拍与性能不回退（06-roadmap.md）。

---

## 1. 现状结构问题

| # | 问题 | 证据 |
|---|---|---|
| S1 | 新旧两套实现并存：legacy reader/window/estimators（reed-core、reed-io）、legacy 过滤（reed-filter `exec.rs` + `spec/legacy.rs` + `plan/legacy.rs` + `main.rs` 二进制）与 v2 并存 | `reed-core/src/{window,population,result,bits}.rs` + `estimators/*`；`reed-io/src/{lib.rs,record/*,sink.rs}`；`reed-filter/src/{exec.rs,main.rs,emitter.rs}` |
| S2 | v2 驱动仍然调用 legacy：读 header 用 `reed_io::record::record_header`（noodles 重新打开文件）、剔除样本列的 header 用 `record::filter_header` | `reed/src/filter.rs:87,158`；`reed/src/main.rs` |
| S3 | GT 文本解析有 5 份，语义各不相同，只有一份执行 ploidy/multiallelic 策略 | `reed-io/src/vcf/sample_scan.rs`（`fast_gt`+`write_gt`，有策略）、`reed-io/src/vcf/gt.rs::tally_gt`（stats 用）、`reed/src/filter.rs::parse_gt` 与 `is_fully_missing`（filter 用）、`reed/src/qc_report.rs::parse_gt_pair`（qc 用） |
| S4 | crate 职责不符：`reed-filter` 依赖 noodles/clap/toml/tracing-subscriber/crossbeam，带二进制；过滤驱动的主体在 `reed/src/filter.rs`（1,699 行）；`reed` 既是门面库又带 clap 二进制 | `reed-filter/Cargo.toml`；`reed/Cargo.toml` |
| S5 | 输出层逻辑（切 BGZF 块、记录虚拟偏移、tid 查找、索引）在 FilterJob 与 EmitJob 中各写一遍；SketchJob 与 FilterJob 的逐样本循环各写一遍 | `filter.rs` `cut_block`×2、`tid_of`×2 |
| S6 | 引擎的顺序模式是另一份实现 | `pipeline.rs::run_sequential` |
| S7 | 错误跨层时变成字符串：job 内的类型化错误（带 `Locus`）被 `EngineError::Job(String)` 吞掉 | `filter.rs`、`stats.rs` 中大量 `.map_err(\|e\| EngineError::Job(e.to_string()))` |
| S8 | M1 的模型层（`GroupTable`、清单 v2、`PlaneBlock`、计数内核、投影、可及性掩码）在生产路径中基本没有被使用 | 01-review §2 M1 |

---

## 2. 删除清单（一次性完成，不保留兼容层）【规】

理由：reed 还没有外部用户，0.x 版本可以直接破坏 API。兼容层的维护成本（S1、S2）已经超过它带来的价值。
oracle crate 保留，它是差分测试的基准。

| 删除 | 替代 |
|---|---|
| `reed-core/src/{bits,window,population,result}.rs`、`reed-core/src/estimators/`（全部） | `reed-core` 只保留模型、内核与数学；估计器只在 `reed-stats` |
| `reed-core/tests/{core,contracts,families,unequal_groups}.rs`（针对 legacy API） | 等价用例迁移到 `reed-stats/tests/`，对 oracle 和 allel 正典对拍 |
| `reed-io/src/lib.rs` 中的 `for_each_window`/`read_windows`/`header_info`/`Assembler` | `reed::stats` 驱动 |
| `reed-io/src/record/`（`for_each_record`、`RecordView`、`decode/apply/count_missing/filter_header`、`GtMask`、`Recompute`、`STALE_INFO`） | `reed-io/src/vcf/*`、`reed-io/src/header.rs`、`reed-filter` |
| `reed-io/src/sink.rs`、`reed-io/examples/scan.rs` | `reed-io/src/write/*` |
| `reed-filter/src/{exec.rs,emitter.rs,main.rs}`、`spec/{legacy,compat}.rs`、`plan/legacy.rs` | v2 + CLI 参数（02-cli.md） |
| `reed-io/src/vcf/gt.rs`、`filter.rs::{parse_gt,is_fully_missing}`、`qc_report.rs::parse_gt_pair` | 统一的 `GenotypeScanner`（§3） |
| `toml` 依赖、`FilterSpecV2` 的 `Deserialize` | 02-cli.md §7 |
| `REED_THREADS` 环境变量 | `--threads` |

保留：`reed-core-v02-oracle`（`publish = false`，只作为 dev-dependency）；`bench/suites/vs-htslib`、`bench/suites/vs-noodles`（冻结的历史对照）。

---

## 3. 统一记录解码（核心重构）【规】

### 3.1 目标

**每种输入格式只有一条解码路径，每个基因型只解析一次，策略只在一处执行。**
stats、filter、qc、index、sketch 全部消费同一种中间表示。

### 3.2 结构

```
reed-io/src/vcf/
  line.rs          LineRef：列偏移（现有）。扩展：样本列数在扫描时顺带得到，不再单独扫一遍
  site.rs          SiteView：CHROM→ContigId（带上次命中缓存）、POS、REF/ALT 分类、N_ALT、QUAL、FILTER、
                   INFO 按 TagId 取类型化值（零分配）
  format_layout.rs FormatLayoutCache：FORMAT 串 → Layout（GT 下标 + 计划内字段下标 + 最后一个需要的下标）。
                   返回 &Layout，不 clone
  genotype_scan.rs GenotypeScanner：单遍扫描样本区，输出 GtCode 序列 + 计划内 FORMAT 数值列；
                   策略（ploidy/multiallelic/半调用）在这里执行并计数
  line_edit.rs     LineEditor（现有，按 03 §5 调整）
reed-io/src/bcf/
  record.rs  decode.rs   BcfRecordRef / 解码到同样的 SiteView + GtCode
```

### 3.3 `GtCode`：每个样本一个字节

```rust
/// 一个样本的基因型编码。低 4 位 = slot0 状态 | slot1 状态 << 2，
/// 状态：0 = REF，1 = ALT，2 = 缺失；bit4 = 定相；bit5 = 原始为单倍体。
#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(transparent)]
pub struct GtCode(u8);

impl GtCode {
    pub const MISSING: Self;                  // ./.
    pub fn state(self) -> u8;                 // 低 4 位，0..=10，共 9 种有效组合
    pub fn masked(self) -> Self;              // 掩码后：两个槽都置为缺失，保留 haploid 位
    pub fn is_fully_missing(self) -> bool;
    pub fn is_called_any(self) -> bool;
}
```

- **为什么不是直接写位平面**：位平面适合 LD、sketch 这类需要整列位运算的场景；但 stats、filter、qc 的主要操作是
  "按群体累加计数"。用字节编码时可以做到**无分支的直方图累加**：
  `hist[group][code.state()] += 1`，每个位点结束时用一张 9 项的常量表把直方图换算成 `SiteGroupCounts`。
  需要位平面的消费者调用 `PlaneBlock::push_codes(&codes)` 转换（一次线性扫描）。
- 多等位策略：`CollapseAlt` 时下标 ≥ 1 记为 ALT；`Reject` 时报错；`SplitPerAlt`（M7）时由扫描器为每个 ALT 产出一组编码。
- 单倍体：`HaploidAsHalfCall` 时 slot1 记缺失并设 haploid 位，掩码时写回 `.` 而不是 `./.`（修复 01-review D17）。
- 这些策略不是 reed 的统计语义，只是解码层如何如实表示输入。统计层按 04-algorithms.md §3 逐等位计数（pixy 的多等位推广），
  因此解码器必须把多等位信息完整保留下来，而不是在解码时折叠。

### 3.4 `GenotypeScanner` 契约

```rust
pub struct GenotypeScanner { layouts: FormatLayoutCache, plan: FieldPlan, policies: GenotypePolicies }

pub struct ScanOutput<'s> {
    pub codes: &'s [GtCode],               // len = n_samples
    pub int_columns: &'s [Vec<u32>],       // 计划内整数 FORMAT 列，缺失为 INT_MISSING
    pub float_columns: &'s [Vec<f32>],     // 计划内浮点 FORMAT 列，缺失为 NaN
}

impl GenotypeScanner {
    pub fn scan<'s>(&'s mut self, line: &LineRef<'_>, site: &SiteView<'_>,
                    counters: &mut PolicyCounters) -> Result<ScanOutput<'s>, IoError>;
}
```

- 每个样本：快速路径（3 字节 GT）→ 否则慢速路径；扫描到 `Layout::last_needed` 后用 memchr 直接跳到下一个 `\t`。
- 在同一遍中统计样本列数，不一致即报 `MalformedRecord`（替代 `LineRef::n_sample_cols` 的额外扫描）。
- 零分配：所有输出缓冲由扫描器持有并复用（修复 01-review D10）。
- 下标解析使用 checked 算术（修复 D14）。
- 差分测试：`GenotypeScanner` 对 noodles 参考解码器（现有 `sample_scan_vs_noodles.rs` 扩展）；fuzz 目标改为它。

### 3.5 消费者

| 消费者 | 用法 |
|---|---|
| stats（窗口统计） | `codes` → 每群体直方图 → `SiteGroupCounts` → 各统计量的累加器 |
| filter（单遍） | 掩码规则读 `int_columns`，命中则 `code.masked()`；随后与 stats 共用同一个直方图计数函数（cohort keep 视为一个群体） |
| filter（sketch） | `PlaneBlock::push_codes` + 掩码位图 |
| filter（reduce pass） | 只要求 GT 的 FieldPlan，`is_fully_missing` 计数 |
| stats（`--summary/--het/--missing/--depth`） | `codes` 直接驱动逐样本计数；`int_columns` 驱动 DP/GQ 直方图 |
| LD（M7） | `PlaneBlock` |

---

## 4. 引擎（reed-engine）修复与简化

### 4.1 真正的票据【规】（修复 01-review D1）

```rust
/// 在途 chunk 的许可。reader 在切出 chunk 前获取，reducer 在按序交付后归还。
struct Tickets { tx: Sender<()>, rx: Receiver<()> }   // 容量 T 的令牌通道，初始放满 T 个令牌
```

- reader：`select!{ recv(tickets.rx) => 取到令牌, recv(abort) => 退出 }`，然后再读块、切 chunk。
- reducer：`take_ready()` 每交付一个 chunk，`tickets.tx.send(())` 归还一个令牌。
- 由此保证在途 chunk id ∈ `[expected, expected + T)`，重排环容量 = T，`insert` 的断言永远成立。
- chunk 通道与结果通道可以设为无界（在途数量已由令牌限制），也可以保持有界，二者都正确。
- 测试：在 `output_is_invariant_across_thread_counts` 之外，增加"chunk 0 人为延迟 500 ms"的测试，验证不 panic、内存不增长。

### 4.2 按字节的内存预算（修复 D5）

令牌改为按字节计：每个 chunk 获取 `compressed_len × 放大系数` 个字节许可（系数默认 8：解压后约 5 倍，加上输出块与部分结果）。
超出预算时 reader 阻塞，而不是报错。`MemBudgetExceeded` 只保留给"单个 chunk 本身就超过预算"这种无法继续的情况。

### 4.3 单一实现

- `run_sequential` 删除。threads ≤ 2 时，用同一个 `run_parallel`，配置 1 个 worker，reader 与 reducer 在同一个线程上交替执行
  （或者直接使用 1 个 worker 线程）。确定性测试因此只需覆盖一条代码路径。
- 去掉 watcher 线程：外部 `AbortSignal` 与内部中止共用同一个关闭通道（调用方传入 handle，内部 clone）。

### 4.4 类型化错误（修复 S7）

```rust
pub enum EngineError {
    Io(std::io::Error),
    Input(reed_io::IoError),                     // 带 Locus
    Job(Box<dyn std::error::Error + Send + Sync>), // 保留原始错误类型与 source 链
    WorkerPanicked { chunk: ChunkId, message: String },
    MemBudgetExceeded { needed_bytes: u64, budget_bytes: u64 },
    Aborted,
}
```

---

## 5. 输出层统一（reed-io/src/write）【规】

```rust
/// worker 本地的 BGZF 写出缓冲。块边界与记录边界无关（与 htslib 一致）。
pub struct BlockWriter { level: u8, deflater: Deflater, pending: Vec<u8>, blocks: Vec<OutBlock>, records: Vec<IndexRecord> }

impl BlockWriter {
    /// 追加一条记录；每当 pending 满 65_280 字节就切出一个块，记录可以跨块。
    pub fn push_record(&mut self, bytes: &[u8], key: IndexKey) -> Result<(), IoError>;
    pub fn finish_chunk(&mut self) -> Result<ChunkOutputBlocks, IoError>;
}
pub struct IndexKey { pub contig: ContigId /* u32 */, pub start1: u64, pub end1: u64 }
pub struct IndexRecord { key: IndexKey, block_idx: u32, offset_in_block: u16 }
```

- 修复 01-review D11：记录跨块时，虚拟偏移指向记录起始所在的块与块内偏移；索引的 chunk 终点为下一条记录的起始偏移。
- 修复 D4：`ContigId` 为 u32，来自 header 的 contig 表；不在 header 中的 contig 报错（或在 `--contig-lengths` 中声明）。
- `OrderedOutput`（reducer 端）：按序 append 压缩块、累计压缩偏移、把 `IndexRecord` 换算成虚拟偏移并交给索引器。
  索引器按 `-W` 选择 CSI（默认，min_shift = 14）或 TBI（位置 ≥ 2^29 时报错）。
- `-O v/u`（未压缩）时 `BlockWriter` 退化为普通字节缓冲，接口不变。
- `-O b`（BCF）：在原生 BCF 编码器完成之前（R6），**直接报错**"BCF output not yet supported, use -Oz"，不得静默写文本（修复 D12）。
- `deflate_block` 对超过 65,280 字节的输入返回错误，`BSIZE` 用 `u16::try_from` 校验。

filter（单遍与 emit）与 stats 的逐位点输出（`--freq`、`--hardy`、`--missing` 的位点表）都通过同一个 `BlockWriter`/`OrderedOutput` 写出。

---

## 6. 输入层补全（reed-io/src/chunk）

| 输入 | 现状 | 目标 |
|---|---|---|
| BGZF VCF（文件） | ✅ | 修复 D2（空块）、D3（XLEN 遍历） |
| BGZF VCF（stdin） | 代码可用但驱动先用 noodles 重开文件读 header，不支持 `-` | header 从 `scan_header` 的文本解析（`HeaderFacts` + noodles `Header` 由同一份文本构造），删除对 `record_header` 的调用 |
| 未压缩 VCF | ❌ | `PlainChunkSource`：文件用 mmap 按 1 MiB 切；stdin 按 1 MiB 读入切，接缝规则相同 |
| 普通 gzip | ❌ | `StreamChunkSource`：一个线程串行 inflate，再按字节切；警告"非 BGZF，解压是串行的" |
| BCF | ❌ | 有索引：按 CSI bin 的虚拟偏移切；无索引：并行 inflate + 串行按长度前缀定界 |
| 区域 `-r` | ❌ | 用索引得到每个区域的虚拟偏移区间，作为 chunk 源的起止 |

所有输入路径统一由 `reed::input::open(path_or_dash, &RegionOpts) -> (Box<dyn ChunkSource>, HeaderBundle)` 构造，
`Format::detect` 在这里接入。

---

## 7. crate 边界（v2 终版）【规】

```
reed-core     模型 + 内核 + 数学                         依赖：bitflags、serde(仅 Serialize)、thiserror
reed-io       格式、BGZF、chunk 源、解码、编辑、写出、索引   依赖：reed-core、noodles、libdeflater、memchr、memmap2
reed-engine   chunk 流水线、票据、重排、中止               依赖：reed-core、reed-io、crossbeam-channel
reed-stats    所有统计量（纯），每个统计量一个文件，按科学方法命名   依赖：reed-core
              （reed-qc 并入，成为 summary/；crate 数从 8 减为 7）
reed-filter   过滤 Spec、表达式、plan、求值、sketch（纯）   依赖：reed-core（不依赖 noodles/clap/toml/crossbeam）
reed          驱动：filter/stats/index 的执行编排            依赖：以上全部、rayon
reed-cli      二进制 `reed`：clap 参数、单位解析、参数→Spec 翻译、错误渲染、mimalloc、tracing-subscriber
```

- reed-filter 的 plan 编译需要 header 中的 INFO/FORMAT 类型：在 reed-core 定义 `TagTable`（`TagId`、名字、`Number`、`Type`），
  由 reed-io 从 header 构建后传入。reed-filter 因此不再需要 noodles。
- `reed` 库不依赖 clap。CLI 相关的一切都在 reed-cli。
- **reed-stats 按统计量组织，文件头写明方法与文献**：每个文件只实现一个统计量的唯一定义（04-algorithms.md §3），
  文件头列出文献、与各软件的已知差异（04 §2）以及测试用的参照物（04 §5）。不存在"同一统计量的多种软件实现"。
- 输出写出器按文件格式分文件（`reed/src/stats/writers/{window_table.rs, plink2_table.rs, vcftools_depth.rs, sfs_files.rs, bcftools_stats.rs, log.rs}`）。
  格式只是文件约定，与统计方法无关。

---

## 8. 目标文件树（只列有变化的部分）

```
crates/
├── reed-core/src/
│   ├── lib.rs  id.rs  coord.rs  contig.rs  sample.rs  group.rs  bitset.rs  site.rs  grid.rs  error.rs
│   ├── tag.rs                      # TagTable / TagId / TagNumber / TagType（新增）
│   ├── genotype/{mod.rs, code.rs, planes.rs, policy.rs}   # code.rs = GtCode（新增）
│   ├── count/{mod.rs, histogram.rs, kernel.rs}            # histogram.rs：GtCode 直方图 → SiteGroupCounts（新增）
│   └── math/{harmonic.rs, hypergeom.rs, hwe.rs, sum.rs}
├── reed-io/src/
│   ├── lib.rs  error.rs  format.rs  header.rs  field_plan.rs
│   ├── bgzf/{block.rs, inflate.rs, deflate.rs}
│   ├── chunk/{mod.rs, bgzf.rs, plain.rs, stream.rs, bcf.rs, region.rs}
│   ├── vcf/{line.rs, site.rs, info.rs, format_layout.rs, genotype_scan.rs, line_edit.rs}
│   ├── bcf/{record.rs, decode.rs}
│   └── write/{block_writer.rs, ordered_output.rs, index.rs}
├── reed-engine/src/{lib.rs, pipeline.rs, tickets.rs, reorder.rs, abort.rs, budget.rs, metrics.rs, error.rs}
├── reed-stats/src/
│   ├── lib.rs  request.rs  accumulator.rs  window.rs     # window.rs：分箱 + 唯一的窗口网格（04 §4）
│   ├── pairs.rs                                          # 逐位点 diffs/comps（π、dxy、Hudson 共用的计数）
│   ├── pi.rs  dxy.rs                                     # Nei & Li 1979 / Nei 1987；缺失数据：Korunes & Samuk 2021
│   ├── fst_wc84.rs  fst_hudson.rs                        # Weir & Cockerham 1984 / Hudson 1992 + Bhatia 2013
│   ├── watterson_theta.rs  tajima_d.rs                   # Watterson 1975 / Tajima 1989
│   ├── pbs.rs                                            # Yi et al. 2010；Malaspinas et al. 2016
│   ├── sfs_projection.rs                                 # Marth et al. 2004；Gutenkunst et al. 2009
│   ├── freq.rs  hardy.rs  het.rs  missing.rs  depth.rs   # 逐位点/逐样本报告；hardy：Wigginton 2005，mid-p：Graffelman & Moreno 2013
│   ├── invariant_check.rs                                # pixy 的不变位点检查
│   └── summary/{mod.rs, sections.rs, histogram.rs}       # 原 reed-qc：bcftools stats 格式的 QC 汇总
├── reed-filter/src/
│   ├── lib.rs  error.rs  funnel.rs  sketch.rs
│   ├── spec/{mod.rs, sample.rs, site.rs, genotype.rs, after_mask.rs, info_update.rs, policy.rs, render.rs}  # render.rs = to_args
│   ├── expr/{mod.rs, lexer.rs, parser.rs, ast.rs, check.rs, eval.rs}   # bcftools 表达式子集，语义照搬
│   ├── plan/{mod.rs, strategy.rs, field_plan.rs}
│   └── eval/{head.rs, mask.rs, tail.rs, fill_tags.rs}               # fill_tags.rs：复现 bcftools +fill-tags（plugins/fill-tags.c）
├── reed/src/
│   ├── lib.rs  input.rs  output.rs  provenance.rs
│   ├── filter/{mod.rs, record_pipeline.rs, reduce.rs, single_pass.rs, sketch.rs, mid_stage.rs, emit.rs, report.rs}
│   ├── stats/{mod.rs, job.rs, sink.rs}
│   ├── stats/writers/{window_table.rs, plink2_table.rs, vcftools_depth.rs, sfs_files.rs, bcftools_stats.rs, log.rs}
│   └── index.rs
└── reed-cli/src/
    ├── main.rs  cli.rs  units.rs  error_render.rs
    ├── args/{common.rs, filter_flags.rs, filter.rs, stats.rs, index.rs}   # filter_flags.rs：filter 与 stats 共用的过滤参数组
    └── translate/{filter.rs, stats.rs}
```

规模约束：单个文件不超过约 600 行（不含测试）；`reed/src/filter.rs` 按上表拆成 8 个文件。
`record_pipeline.rs` 是"头部 → 扫描 → 掩码 → 计数 → 尾部 → 编辑"的唯一实现，single_pass 和 sketch 都调用它（修复 S5）。

---

## 9. 其他结构性建议

- **中段并行**：`mid_stage` 用 rayon 按 chunk 并行（每个 chunk 的判定相互独立），结果按 chunk id 存入 `Vec`（下标即 id），
  不用 `HashMap`。逐位点的掩码字与 INFO 值存放在每个 chunk 一块的扁平数组里，不再为每个存活位点分配 `Vec`。
- **sketch 溢写**：超过预算时按 chunk 写入 `--tmp-dir`（每个 chunk 一个文件，二进制定长格式），emit 阶段按需读回。
- **stats 输出流式化**：`StatsSink` 不再把全部输出拼成一个 `String`，改为按文件流式写出（逐位点输出时是必需的）。
- **`ContigBins`**：现用 `BTreeMap<u64, Vec<StatAcc>>`。在水位线推进的前提下，活跃分箱是连续的一段，可以改用 `VecDeque`；
  `seed.clone()` 只在新分箱首次被触及时发生（现已如此）。属于可选优化。
- **header 的唯一来源**：`HeaderBundle { text: Vec<u8>, noodles: vcf::Header, facts: HeaderFacts, tags: TagTable, contigs: ContigTable, samples: SampleTable }`，
  由 `scan_header` 的文本一次性构造，所有驱动只用它。
