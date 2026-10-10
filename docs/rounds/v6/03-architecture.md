# 03 架构：建成形态、缺陷清单、最优重构方案

> 结论先行：**并行骨架是对的、已在 Amdahl 极限上，不要动它**；要动的是三件事——① 全局分析被塞进串行归约点（02 §2.2）；② `-Ob`/BCF 输入的文本往返与一批热路径浪费；③ 三套 pass0 与两套按记录体的重复代码。每条带 `file:line`；来自子审计（engine/core、reed-io）的结论都经我复核。

---

## 1. 建成形态

```
reed-core 4.0k  ← reed-io 7.0k ← reed-engine 1.5k
reed-core      ← reed-stats 4.0k
reed-core, io  ← reed-filter 4.2k
全部           ← reed 7.5k（门面，lib.rs 只 pub use）← reed-cli 1.9k（只依赖 reed）
oracle 2.9k    仅 dev-dependency
```

前后端边界（v5 08 §6）：`reed-cli/Cargo.toml` 只依赖 `reed`（守卫 ok）；后端 `eprintln!` 归零；门面 0 条 `pub mod`。未完成：`translate/stats.rs` 仍调 `reed::plan::compile` 并自行 `assign_groups`；没有统一的 `reed::Error`+`kind()`；stats/index 的往返测试缺（v5 STATUS "未完成/偏差"）。

关键机制（已核实）：BAPP 块对齐 chunk（mmap 源）、票据 + 字节许可、重排环、BCF 长度前缀链（先发布后解码、毒化、全局静默看门狗）、统一记录管线 `RecordPipeline`（文本与 BCF 原生两个入口共用 `run_stages`）、worker 本地 `BlockWriter` + 有序追加、原生 BCF 编码器。

---

## 2. 引擎（reed-engine）：形态正确，补四处

| # | 项 | 证据 | 处理 |
|---|---|---|---|
| A6-1 | 中止通道容量按 `available_parallelism()+4`（`abort.rs:63-67`），不按 `--threads`；worker 失败路径不释放字节许可（`pipeline.rs:355-368` 直接 `break`） | `--threads` 大于核数 + 小内存预算时 reader 可能卡在 `BytePermits::acquire` 轮询（`budget.rs:285-302`）而永不退出 → 挂起而非 `Aborted` | 通道容量取 `budget.threads`；`AbortSignal` 加 `AtomicBool`；在途 chunk 用 drop guard 释放许可 |
| A6-2 | 文本路径接缝片段无上限（`pipeline.rs:533-574`），BCF 路径有 1 GiB 上限（`bcf_input.rs:26-28`） | 无换行的病态输入把 reducer 的片段撑到文件大小 | 共享 `MAX_RECORD_BYTES`；抽出 `SeamStitcher`，换行与 BCF 两套接缝实现合一 |
| A6-3 | 顺序模式：`budget.rs:46-50` 与 `docs/spec/architecture.md:804` 承诺 `threads ≤ 2` 走顺序路径，`run()` 无条件 `run_parallel`（`pipeline.rs:157`） | 规范偏差 | 实现 `run_sequential`（复用 `stitch_seam`）或改规范 |
| A6-4 | 计量缺 reader 等票据/许可、reducer 等 recv 的时间；inflate / 计算 / deflate 融合在 `worker_busy` 里 | "reader 受限"仍靠推断；deflate 占比只能用级别对照实验反推（04 §4） | `StageMetrics` 加 `reader_wait`/`reducer_wait`；`ChunkJob::stage_nanos()` 回报 inflate/job/deflate 三段 |

不改：票据 + 重排环（容量 = 票据数使窗口不变量可证，`reorder.rs:177-191`）、长度前缀链、单 reducer 线程、融合 worker。工作窃取 / 两段流水 / NUMA 都没有收益：chunk 均匀（1 MiB），filter 四个线程点 user 恒定、效率 96–100%（04 §3）。

---

## 3. I/O 层（reed-io）：热路径浪费与一处潜伏错误

| # | 项 | 证据 | CPU / 风险 | 处理 |
|---|---|---|---|---|
| A6-5 | deflate：每块 `Vec::with_capacity` + `dst.resize(…, 0)` 清零约 65 KiB（`block_writer.rs:183`、`deflate.rs:134-135`） | ≈ 350K 块 × 65 KiB ≈ 23 GB memset | 5–8 CPU-s | worker 常驻压缩缓冲；级别策略见 04 §8 P6-2 |
| A6-6 | inflate：每块 `dst.resize(start + 64 KiB, 0)`（`inflate.rs:49`）并重复解析块头（`inflate.rs:36` vs `chunk/mod.rs:182-188`） | ≈ 230K 块 × 64 KiB ≈ 15 GB memset | 3–5 CPU-s | 常驻缓冲、传入已解析的 `total` |
| A6-7 | **流式 chunk 源的偏移漂移**：`fill_to_target` 把空成员从缓冲里 `drain` 掉（`chunk/bgzf.rs:131-135`），随后 `carry_offset = compressed_start + cut`（`:310`）按被删短的缓冲计算，之后所有 chunk 的 `compressed_start` 比真实文件偏移少了被删字节；`VirtualMap`（`bcf_input.rs:141`、`index.rs:156,230`）用它算虚拟偏移。mmap 源保留空成员（`:425-432`），两种源的 chunk 元数据不同，与 `:326-330` 的"相同"注释矛盾 | 普通文件走 mmap，所以 `reed index` 与区域查询今天没错；stdin/管道走流式源，含空成员时虚拟偏移错。**潜伏、无测试** | 记录被丢弃的字节数并计入 `carry_offset`；加"流式 vs mmap 元数据相等"测试（含空成员夹具） |
| A6-8 | 1 MiB chunk 的尾部定位用标量 `rposition`（`chunk/mod.rs:231-246`） | 7 GB 标量扫描 | 1–3 CPU-s | `memchr::memrchr` |
| A6-9 | `-Ob`：编辑后的文本被原生编码器重新解析（`native_encode.rs:131-139`），每样本一个 `Vec`（`:377,615-654`） | 实测 +61 CPU-s | 13% | §8 的记录模型 |
| A6-10 | BCF GT 查表（128 KiB，`bcf/gt.rs:195-196`）每 worker 建一次（`record_pipeline.rs:330`） | 56 × 65K 次解码 | 1–3 CPU-s | `Arc` 全进程一份 |
| A6-11 | 渲染/编辑中的分配：`format_g6` 返回 `String`（`line_edit.rs:250,276-302`）、`edit_info` 每记录 `edited`/`seen` 两个 Vec（`:123-124`）、`Locus` 每记录 `String`（`vcf/site.rs:373-376`、`bcf/render.rs:73-80`） | 每存活记录数次分配 | 3–6 CPU-s | 栈缓冲 / 复用 |
| A6-12 | 索引写出：每记录 contig 名转 `String`、构造 noodles `Position`/`Chunk`（`ordered_output.rs:436-440`、`write/index.rs:591-618`）；htslib 的 `hts_idx_push` 是纯整数 | 只影响 `-W` | 按 tid 缓存 |
| A6-13 | `field_plan.rs` 无任何消费者；`error.rs:61-63` 引用不存在的 `ReadError`；`line_edit.rs:8-10` 的"无编辑不走编辑器"没有调用点兑现（`single_pass.rs:326,388` 无条件调用） | 死代码 / 漂移 | 删 / 修 / 空计划时直接拷贝 |
| A6-14 | 渲染两处与 htslib 的差：VCF ≥ 4.4 前缀相位（`vcf.c:6500-6518`）未实现；INFO 浮点 `-0.0` 打成 `0`（`line_edit.rs:260-262`，`kputd` 保留 `-0`） | 手工构造的相位位 / 负零会与 bcftools 不同字节 | 对齐并加夹具 |

noodles 仍在 header 解析、索引读写、BCF header 文本（一次性工作）与测试 oracle 里；热路径无 noodles。合理，但 `HeaderFacts.header`/`FileIndex`/`Indexer` 把 noodles 类型放进了内部公开签名——升级 noodles 即破坏 API，接受但记录。

---

## 4. 过滤层（reed-filter + reed/src/filter）：每记录/每样本的浪费

| # | 项 | 证据 | 处理 |
|---|---|---|---|
| A6-15 | 每存活记录 `stale: Vec<&[u8]>`（`single_pass.rs:290,332`）、`fill_info_set` 把 INFO 值格式化成文本（`edit_plan.rs:18-37`）、`edit_info` 的两个 Vec | 1.77M 次 × 数次分配 | `stale` 预先算在 `FilterPlan` 里；输出用栏缓冲 |
| A6-16 | `--fill-tags HWE`：`SiteCounts::hwe_p` 调 **未缓存** 的 `hwe_exact_p`（`site_counts.rs:71-81`），而 `HweCache` 存在且被 `--hardy` 用（`math/hwe.rs:130-172`） | 每存活位点一次完整分布计算 | 走缓存（`(n_called, mac)` 大量重复） |
| A6-17 | 掩码阶段逐样本把整型 DP/GQ 转 f64 再比较（`record_pipeline.rs:357-372`）；`--set-gts .` 表达式路径**每样本分配一个 `Vec<FmtValue>`**（`:395-399`） | 348 × 规则数 × 3.3M 次转换；表达式路径每记录 348 次分配 | 列式掩码：对 `Vec<u32>` 列整列比较出位图（整数比较，可 SIMD），再 OR 进 `masked`；表达式路径复用一个 `Vec` |
| A6-18 | 表达式求值是 AST 递归 + 每个 INFO tag 重扫 INFO 串（`expr/mod.rs` `eval`/`value`/`info_value`），而 head 规则路径有一次性片段索引（`eval/head.rs:190-214`） | gold 表达式 4 个 tag → 每记录最多 4 次 INFO 扫描 + 4 次 `from_utf8().parse::<f64>()` | 一记录一次片段索引供两路共用；预编译为闭包/字节码（AST 节点不再在热路径匹配字符串算符 `"&&"`） |
| A6-19 | `masked` 位图逐位清零（`record_pipeline.rs:345-348`） | 348 次 `set(false)`/记录 | `fill(0)` |
| A6-20 | 两个几乎相同的"编辑 + 写出"闭包（文本路径 `single_pass.rs:353-412`、BCF 路径 `:278-346`） | 约 40 行重复 | 抽成 `emit_survivor(view, …)` |
| A6-21 | 三套 pass0（`MissingJob`、`SketchJob`、`LdPruneJob`）各自实现 `run_chunk`/`run_seam`/漏斗重置与一套 sink | 维护成本；BCF 输入强制 two-pass（`filter/mod.rs:209-211`） | 统一为 §7 的矩阵 pass0 |

策略选择（`plan/strategy.rs:48-78`）与计划编译的分层是清楚的；`FilterSpec`/`MissingRate`/`to_args()` 往返的设计保留。

---

## 5. reed-core：内核的 soundness 与死面

| # | 项 | 证据 | 处理 |
|---|---|---|---|
| A6-22 | `count_site_group` 只用 `debug_assert_eq!` 校验掩码/平面长度（`kernel.rs:46`），AVX-512 循环按平面长度做裸加载（`:121-133`）：release 下更短的 `group_words` = 越界读 = **安全代码可达的 UB**；SAFETY 注释说"dispatcher 断言"（`:99-104`）在 release 不成立。且 348 样本下向量循环不执行（W=6<8），工作区无生产调用者 | 受批准的 unsafe 例外实际不健全 | 二选一：**删除** SIMD 路径与 `CountTable`（计数 <5% CPU，不值得维护），或改为 release 断言 + 带检查的 `SitePlanes` 构造器 + 填充到 8 字 |
| A6-23 | 两个 `ContigId`（02 §1.5） | 根导出裸 u32 | 统一 newtype |
| A6-24 | `math::hypergeom::project` 与 `reed_stats::sfs::projection_weights` 两套投影权重、截断规则不同；前者无调用者 | 一问两答 | 删一个 |
| A6-25 | 死面：`ChunkBytesExport`、`MemBudget::from_half_memory_for_tests`、`is_sequential`（无调用）、`PlaneBlock` 复用块不清零依赖 `codes.len()==n_samples`（`planes.rs:175-179`） | — | 删 / 加 debug 断言 |

---

## 6. 统计层与门面

- `StatKind`/`StatAcc` 封闭枚举 + 并列的 genome-wide partial 是对的（02 §1.3）；`job.rs` 1,167 行把请求构造、运行、写出放在一起，v5 C-7 未做。
- `StatsSink` 既产出结构化 `WindowRow` 又自己渲染 `rows: String`（`stats/job.rs:1296-1335`）：数据与文本混在 sink 里；按 v5 08 §6.2 应是 sink 产数据、`reed::writers` 产文本。
- LD 的串行 sink（02 §2.1）。
- CLI：`translate/stats.rs` 编译计划、分配群体——后端逻辑在前端（v5 STATUS 已登记）；`--ld-decay` 被归为"窗口统计"而要求 `--pop`（02 §3.2）。

---

## 7. 最优方案一：全局分析的矩阵化 pass0（替代三套 pass0）

见 02 §2.2 的模型与代码示意。架构上的落点：

```
reed/src/filter/
  pass0/
    mod.rs        // 统一入口：需要哪些全局量（样本已调用计数、位平面矩阵、主等位频率）由 FilterSpec 决定
    matrix.rs     // MatrixJob / MatrixSink：worker 填 PlaneBlock，reducer 拼接
    analyses.rs   // mind()（样本保留位图）、prune()（contig 并行）、decay()（区块并行）
  emit/           // pass1：现有 single_pass 的发射，输入 = 保留位图 + 样本位图
```

- `--mind`：从矿阵的 `slot*_miss` 平面按列 popcount 得每样本已调用数（不再需要 `MissingJob`/`SketchJob` 的单独逻辑）。
- `--indep-pairwise`：`prune_all`（contig 并行）。
- stats 的 `--ld-decay`：同一 `ContigMatrix`，`ld_decay()` 区块并行；stats 的 LD 不再走窗口 partial 的 head/tail 携带。
- 内存：`sites × 192 B` 进预算；超预算按 contig 分段。
- 门：剪枝集合与 plink2 逐 ID 相同（夹具 + 区域 + 全量 235,403）；衰减 r² 与现有实现逐对相同；`--mind` 结果与现有 sketch/two-pass 逐字节相同；1/8/56 线程输出相同；剪枝 ≤ 3 s、衰减 ≤ 10 s。

## 8. 最优方案二：`-Ob` 与 BCF 输入的记录模型

```rust
// reed-io/src/record.rs —— 格式中立、可修改的一条记录（只在需要重编码的路径上使用）
pub struct Record<'a> {
    pub site: SiteFields<'a>,                // CHROM/POS/ID/REF/ALT/QUAL/FILTER（借用原字节或 BCF 字典下标）
    pub info: Vec<InfoItem<'a>>,             // (tag id, 类型化值或原字节)
    pub format: FormatBlock<'a>,             // 键序 + 每键的类型化列（GT 为 GtCode 列）
    pub n_samples: usize,
}
pub trait RecordSource { fn decode(&mut self, bytes: &[u8], out: &mut Record<'_>) -> Result<(), IoError>; }   // VcfText / Bcf
pub trait RecordSink   { fn encode(&mut self, rec: &Record<'_>, out: &mut Vec<u8>) -> Result<(), IoError>; }   // VcfText / Bcf
// 编辑 = 改模型：rec.info.retain(not stale) + push(recomputed)；rec.format.gt[masked] = GtCode::MISSING
```

- `-Oz` 文本输出**保持**现在的字节级编辑（最省）；`-Ob` 与 BCF→任何输出走记录模型：解码一次、改模型、直接编码，去掉"渲染→解析→编码"。
- BCF→BCF 时未修改的 FORMAT 列按 htslib `bcf_write` 的思路**原字节拷贝**（只重写 INFO 与 GT 列）。
- 门：`-Ob` 读回正文 ≡ `-Oz` 正文（现有测试）；noodles 编码器 oracle 字节相等；全量 `-Ob` CPU −50 s 以上。

## 9. 最优方案三：列式掩码 + 一次性 INFO 片段索引

```rust
// 掩码：对每条规则整列比较，产出位图，再 OR（整数比较；编译器可自动向量化）
fn mask_column_lt(col: &[u32], thresh: u32, missing: u32, out: &mut Bitset) { for (i, &v) in col.iter().enumerate() { if v != missing && v < thresh { out.set(i, true); } } }
// INFO：一记录一次片段索引（HeadScratch.fragments 已有），表达式与 head 规则共用
struct InfoIndex<'a> { frags: &'a [FragRange] } impl InfoIndex<'_> { fn get(&self, tag: &[u8]) -> Option<&[u8]> { /* 线性/小表查找，无重扫 */ } }
```

门：子集 263,420 / 381,672 不变；gold md5 不变；`stats_site_filters_match_the_filter_then_stats_composition` 通过。

## 10. 不做的事

- 不改并行骨架（票据/环/链/单 reducer/融合 worker）。
- 不把 2 位编码定为规范表示（02 §1.1）。
- 不引入 `target-cpu=native`（可移植性规范）。
- 不做 Python 绑定（裁决）。
