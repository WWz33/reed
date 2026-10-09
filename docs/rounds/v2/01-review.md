# 01 复审：M0–M5b

> 复审对象：`6e41bfd`（M0+M1）→ `2df6eb0`（M2）→ `04c1f1a`（M3）→ `e7de0c4`（M4）→ `b6d1c4c`（M5）→ `771c71c`（M5b），
> 以及工作区中未提交的 reed-qc 进行中改动（`reed-qc/src/{histogram,report}.rs`、`reed/src/qc_report.rs`）。
> 方法：逐文件阅读 + `cargo build/test/clippy --release --workspace` + 全量数据实测 + 与 bcftools gold 对拍。
> 所有数字均为 2026-10-06 本机实测（热页缓存），命令见各节。

---

## 0. 总评

实现者把 v1 spec 的骨架**落地得很好**：BAPP、融合式 chunk 流水、worker 内 deflate、边写边建索引、
可加累加器 + 分箱、sketch-then-emit 全部建成，并且性能比 spec 的预算好一个数量级。
过滤输出与 bcftools gold 完全一致。

问题集中在四类：

1. **几个会在生产数据上出事的缺陷**（P0）：重排窗口 panic（测试已红）、BGZF 中途空块截断、
   超长记录（> 65 280 字节）写出损坏的 BGZF、tabix 索引的 contig 编号截断与 512 Mb 上限、`.bcf` 输出实际写成文本。
2. **统计层的输入契约在重写时丢失了**：contig 长度必需、位点升序校验。另外，统计定义本身（complete-case、清单门、reed 自定义窗口网格）
   不对应任何主流实现；按 2026-10-06 的对齐裁决（不做算法创新），这些需要换成经典方法 + pixy 缺失数据处理的唯一定义（04-algorithms.md）。
3. **新旧两套实现并存**，加上 5 份 GT 解析器，策略（ploidy/multiallelic）只在其中一份里生效。
4. **CLI 形态被否决**（见 02-cli.md）。

---

## 1. 实测结果【据】

### 1.1 性能（本机全量，56 线程，热缓存）

| 命令 | 实测 wall | v1 目标 | 对照 |
|---|---|---|---|
| `gt_scan clean.vcf.gz`（M2 门，bench 日志） | 3.14 s | ≤ 20 s | v0.2 reader 283 s/遍 |
| `qc_summary`（M3 门，bench 日志） | 2.38 s | — | bcftools stats 单线程 |
| `reed stats clean.vcf.gz -g clean.3pop.txt`（本次复测） | **2.65 s** | ≤ 30 s | 18,992 窗 × 51 列 |
| `reed filter … --conf filters.conf -o x.vcf.gz --write-index`（本次复测） | **14.2 s** | ≤ 60 s | bcftools 约 40 min，plink2 约 2.5 min，v0.2 9.6–40 min |

**问题**：M4、M5、M5b 的验收数字**没有写入** `bench/suites/e2e/results.tsv`（只有 M0/M2/M3 三行）。
v1 §17.3 规定"文档中引用的任何性能数字都必须能在日志里找到"，这条没有执行。
上表后两行应补登（R0，06-roadmap.md）。

### 1.2 过滤正确性（全量对拍）

| 检查 | 命令 | 结果 |
|---|---|---|
| 存活位点数 | `reed filter` 报告 | 3,315,604 → 1,768,759（与 gold 一致） |
| 位点 + 全部 GT | `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'` md5 | **一致**（`e9efcc63…`） |
| AC/AN/NS | `bcftools query -f '%CHROM\t%POS\t%AC\t%AN\t%NS\n'` md5 | **一致**（`09459e65…`） |
| MAF/HWE | Chr01:1–3,000,000（5,299 位点）逐行 diff | **一致** |
| BGZF 完整性 | `bgzip -t` | 通过 |
| 索引可用性 | `bcftools view -r Chr05:1000000-1100000` | 119 行，正常 |
| 与 v0.2 reed 输出逐字节对比 | 正文 md5 | 不同（预期内）：INFO 顺序和浮点渲染改为 bcftools 风格 |

注：新输出的 INFO 追加顺序为 `AC;AN;NS;MAF;HWE`，gold（bcftools fill-tags）为 `NS;AN;MAF;AC;HWE`。
值完全一致，L_tags 通过；如果要做到正文逐字节一致，需要在 `fill_info_set` 中采用 fill-tags 的追加顺序（P2）。

### 1.3 构建与测试

- `cargo build --release --workspace`：通过，有 13 条 rustc 警告（未使用变量、多余 `mut`）。
- `cargo test --release --workspace`：**失败 1 项**——
  `reed-engine/tests/pipeline.rs::output_is_invariant_across_thread_counts`，
  panic 于 `reed-engine/src/reorder.rs:58`（`chunk … outside the reorder window`）。连续运行 5 次失败 4 次（偶发但高频，取决于线程调度），根因见 §3 D1。
- `cargo clippy`（pedantic=warn）：约 700 条警告，其中 277 条 `expect()`、18 条 `u64→u32` 截断、
  17 条 `usize→u32` 截断、15 条"大小写敏感的扩展名比较"。截断类警告对应真实缺陷（D4、D7）。

---

## 2. 逐里程碑评价

### M0 卫生 — 基本完成

- ✅ workspace 统一 package/deps/lints/profile、edition 2024、clap、tracing 依赖、`rust-toolchain.toml`、`deny.toml`。
- ⚠ `reed-filter` 的 legacy 二进制仍然存在（`reed-filter/src/main.rs`），与 `reed` 二进制功能重叠。
- ⚠ 基线测量（585.7 s，subset）入了日志，但 full 基线没有重测。

### M1 core 模型 — 建成但**大部分没被使用**

- ✅ `id/coord/sample/group/bitset/genotype(planes,policy)/count(kernel)/manifest v2/grid/math(harmonic,hypergeom,hwe,sum)` 均已实现，有 proptest。
- ✅ `reed-core-v02-oracle` 冻结并在 `reed-stats/tests/oracle_parity.rs` 中使用，这是正确做法。
- ❌ **生产路径几乎不用 M1 模型**：
  - stats 驱动使用 `reed_io::vcf::gt::tally_gt` 逐样本累加（`reed/src/stats.rs:299`），不使用 `PlaneBlock`/计数内核；
  - stats 驱动用自写的 `assign_groups`（`stats.rs:155-178`），不用 `GroupTable::from_tsv`；
  - 清单 v2（`Ascertainment`/`Availability`/`BiasCode`）只在 filter 的 provenance 中出现，统计层完全不读（按对齐裁决，这一层整体删除，04 §1）；
  - `math::hypergeom`（投影）、`AccessibilityMask`、`WindowGrid::Sites/Regions` 未被任何驱动使用。
- ❌ legacy 模块（`window.rs`、`population.rs`、`result.rs`、`bits.rs`、`estimators/*`）仍在 `reed-core` 中且对外公开，
  与 `reed-stats` 中的新实现并存。

### M2 io 传输（BAPP）— 性能优秀，有两个正确性缺陷

- ✅ 只读块头的 reader、libdeflate、chunk 头/体/尾划分、接缝链式拼接、差分测试（`sample_scan_vs_noodles.rs`）。
- ❌ D2：流中间的 EOF 空块如果恰好落在 chunk 起点，读取就此结束（`chunk/bgzf.rs:90`），拼接文件会被概率性截断。
- ❌ D3：`parse_block_size` 假设 `BC` 子字段固定在偏移 12（`bgzf/block.rs:36-44`），不遍历 XLEN 子字段。
- ⚠ `SampleScanner` 有分配和多余扫描（D10）。
- ⚠ 只支持 BGZF VCF；plain VCF、stdin、BCF、普通 gzip 都没有 v2 输入路径，`Format::detect` 已写但未接入任何驱动。

### M3 engine — 骨架正确，背压模型错误

- ✅ 首错即停、`select!` 中止、worker panic 捕获、按序归约、接缝在 reducer 处理。
- ❌ D1：票据窗口没有真正实现，重排环可以溢出并 panic（测试已红）。
- ❌ D5：内存预算是"超出就报错"，不是背压；而且三个驱动都传入 `u64::MAX`，预算形同虚设。
- ⚠ 顺序模式 `run_sequential` 是另一份独立实现，与并行模式靠测试保持一致；应改为同一实现的退化形态（03 §4）。
- ⚠ 每次 `run` 额外起一个 watcher 线程转发外部中止。

### M4 stats v2 — 快且与 oracle 一致，但丢了护栏

- ✅ 七族累加器 + 分箱 + 水位线输出 + contig/genome 层级；与 oracle 逐位对拍；累加器定律测试。
- ⚠ D6：清单门丢失。按对齐裁决不再恢复，改为 pixy 的不变位点检查 + `--variant-only` 声明。
- ❌ D7：缺少 contig 长度时用 `last_pos0 + 1` 编造（`stats.rs:455,495`），违反 DESIGN §3 的"contig_len 一等公民"。
  per-bp 统计量的分母因此静默出错。
- ❌ D8：不再校验 contig 内位置升序和 PositionBeyondContigLength；乱序输入会把位点推入已经输出过的分箱，被静默丢弃。
- ❌ D9：群体文件首行被无条件当表头跳过（`stats.rs:162`）。pixy 和 vcftools 的群体文件没有表头，第一个样本会被静默丢掉。
- ⚠ 输出列名用组下标（`pi.pair.0`），群体名只打印在 stderr；宽表里 SFS 用 `;` 拼接。
- ⚠ 默认计算全部 51 个统计量，没有选择参数。
- ⚠ 整个输出先在内存里拼成一个 `String`，最后才写出（窗口数量少时可以接受，逐位点输出时不可接受）。
- ❌ SFS/θ/Tajima 是 complete-case，这不对应任何主流实现（pixy、vcftools、allel 都不丢弃有缺失的位点），按 04-algorithms.md §1 撤下。
  这里修正 v1 的一个表述：**原始数据**上 complete-case 还能保留大部分位点（Chr01 首窗 pop1：483 中 375，78%）；
  但在**过滤掩码后的数据**上只保留 1.2%（gold 输出 Chr01:1–2 Mb，pop1：4,051 中 49）。
  v1 引用的 profile 缺失率是掩码后口径。用户实际拿去做统计的正是过滤后的数据，所以 complete-case 在实际使用中不可用。

### M5 / M5b filter v2 — 正确且快，但输出层有致命缺陷，代码结构需要拆

- ✅ plan 编译期按 header 校验 INFO 类型/Number；漏斗按规则归因；掩码归因；provenance 与 ascertainment 写入 header；
  SketchThenEmit 在本数据集自动选中，14.2 s；GT/标签与 gold 一致。
- ❌ D4：tabix 索引的 tid 用 `as u16` 截断，未知 contig 默认为 0（`filter.rs:863-864, 1647-1648`）；
  `pos1` 超过 u32 时变成 `u32::MAX`；并且总是写 TBI，TBI 不能索引超过 2^29（512 Mb）的位置，小麦、大麦等植物染色体会超出。
- ❌ D11：单条记录超过 65 280 字节时，块缓冲直接超过 BGZF 上限，`deflate.rs:55` 的 `(total - 1) as u16` 截断，写出损坏的块。
  宽队列（例如 2,000 样本 × PL）很容易出现这种行。
- ❌ D12：输出路径不是 `.gz` 就写纯文本（`filter.rs:348` `is_gz`），`-o out.bcf` 得到的是一个名为 `.bcf` 的文本 VCF。
- ❌ D13：`[policy] ploidy/multiallelic` 在 v2 过滤路径里没有生效：`parse_gt`（`filter.rs:906`）把单倍体一律当半调用，
  多于两个等位静默截断，非数字记为缺失；`policy_counters` 永远是默认值。
- ⚠ 掩码路径仍然用 `keep.get(i).copied().unwrap_or(true)`（`filter.rs:756`），v1 的 A13 没修。
- ⚠ 中段（mid stage）是串行循环，每位点分配 `rule_hits_of` 和 `mask_words` 两个 `Vec`。
- ⚠ `reed/src/filter.rs` 1,699 行，包含三个 job（Missing/Filter/Sketch+Emit）、两份 `tid_of`/`cut_block`/块切分逻辑。
- ⚠ `reed-filter` crate 依赖 `noodles`、`clap`、`toml`、`tracing-subscriber`、`crossbeam-channel`，并且带 legacy `exec.rs`（988 行）。
  v1 要求的"L3 纯策略层"没有实现。

### 进行中：reed-qc（未提交）

- ✅ bcftools stats 各节语义按 `vcfstats.c` 钉住（SN/TSTV/SiS/AF/QUAL/ST/DP/IDD/PSC）。
- ⚠ `reed/src/qc_report.rs:406` 又写了一个 `parse_gt_pair`（第 5 份 GT 解析器）。
- ⚠ CLI 没有接入。按 02-cli.md §1，它并入 `reed stats --summary`（bcftools stats 格式），代码并入 reed-stats 的 `summary/`。

---

## 3. 缺陷清单

| # | 级 | 位置 | 缺陷 | 触发条件 / 后果 | 修法 |
|---|---|---|---|---|---|
| D1 | **P0** | `reed-engine/src/pipeline.rs:285`、`reorder.rs:58` | "在途 ≤ tickets + workers"的推理不成立：reducer 把结果从通道搬进环，腾出通道空位，reader 可以继续发 chunk。一个慢 chunk 就会使后续 id 超出环窗口 | 主线程 `assert!` panic，进程崩溃。测试已经复现 | 真正的票据：reader 在发 chunk 前获取许可（`Semaphore` 或容量为 T 的令牌通道），**只在 reducer 按序交付后归还**。环容量 = T。见 03 §4.1 |
| D2 | **P0** | `reed-io/src/chunk/bgzf.rs:90` + `next_chunk` | `read_block` 把 EOF 空块当作"流结束"返回 `None`。中途的 EOF 块落在 chunk 内部时只会让当前 chunk 提前结束（无害）；但如果恰好是新 chunk 读到的第一个块，`next_chunk` 返回 `None`，整个流就此结束 | `cat a.vcf.gz b.vcf.gz`、`bcftools concat --naive` 产生的拼接文件，每个接合处约有 1/16 的概率（一个 chunk 约 16 个块）**静默丢弃其后全部数据**，且结果取决于块对齐，难以复现 | 空块（ISIZE=0）一律跳过，只有底层读到 EOF 才结束；最后一块不是 EOF 标记时发出警告（htslib 行为）。增加"拼接文件"fixture，并在每个块偏移处放置 EOF 块做穷举测试 |
| D3 | P1 | `reed-io/src/bgzf/block.rs:36-44` | 假设 `BC` 子字段固定在偏移 12 | 合法但带其他 extra 子字段的 BGZF（非 bgzip 产生）→ 被判为"非 BGZF" | 按 XLEN 遍历子字段查找 `BC` |
| D4 | **P0** | `reed/src/filter.rs:863-864,1647-1648`；`OutputSink::accept_block` | tid 用 `as u16` 截断、未知 contig 默认为 0；pos 转 u32；只写 TBI | 超过 65,535 个 contig 的草图基因组 → 索引错位；header 未声明的 contig → 记到第 0 号 contig；位置 > 512 Mb → TBI 失效 | tid 用 `usize`/`u32` 并对未知 contig 报错；`-W` 默认写 CSI（与 bcftools 一致），`-W=tbi` 时若有位置 ≥ 2^29 则报错 |
| D5 | P1 | `reed-engine/src/pipeline.rs:183`；三个驱动传 `MemBudget::from_bytes(u64::MAX)` | 内存预算超出即失败，不阻塞；而且实际没有设置预算 | 预算功能不存在 | 把字节预算并入票据（按字节计的许可），超出时阻塞 reader；`--max-mem` 由 CLI 传入 |
| D6 | P1（对齐） | `reed-stats/*`、`reed/src/stats.rs` | 统计定义不是任何经典方法 + 主流缺失数据处理的组合：π/dxy 是 reed 的两种估计器，SFS/θ/Tajima 是 complete-case，窗口网格是 reed 自定义的；清单门丢失 | 输出数字无法与 pixy/vcftools/allel 中任何一家逐数对上 | 按 04-algorithms.md §3 的唯一定义重写（经典方法 + pixy 缺失数据处理）；不恢复清单门，改为 pixy 的不变位点检查与 `--variant-only` |
| D7 | **P0（语义）** | `reed/src/stats.rs:455,495` | 缺少 contig 长度时编造 `last_pos0 + 1` | Beagle 输出等缺 `##contig` 的文件 → 末窗长度和 per-bp 分母静默出错 | 缺长度即报错；可用 `--contig-lengths FILE`（.fai 格式）补充 |
| D8 | **P0（语义）** | stats 驱动整体 | 不再校验位置升序与是否超出 contig 长度 | 未排序输入 → 位点推入已输出并释放的分箱，静默丢失 | worker 检查 chunk 内部顺序，reducer 检查跨 chunk 的衔接（v1 §8.3 第 4 条），错误类型复用 `PositionNotAscending`/`PositionBeyondContigLength` |
| D9 | P1 | `reed/src/stats.rs:162` | 群体文件首行无条件跳过 | 无表头文件（pixy 格式，正是 `--pop` 采用的格式）丢第一个样本 | 按 pixy 的读法：两列、无表头；统一使用 `GroupTable::from_tsv`；群体文件中出现 VCF 没有的样本时，行为与 pixy 一致（报错），`--force-samples` 时警告 |
| D10 | P1（性能） | `reed-io/src/vcf/sample_scan.rs:149,167,179`；`line.rs:107` | `scan_fields` 每个**样本**分配一个 `Vec`；`layout_for` 每条记录 clone 一次 `Layout`（含 `Vec`）；`n_sample_cols` 额外完整扫描一遍样本区；`split_subfields` 切分全部子字段而不是在最后一个需要的子字段处停止 | 过滤热路径上每条记录 348 次分配；本数据集上仍然只要 14 s，但宽队列会被放大 | 统一记录解码器（03 §3）：零分配回调、FORMAT 布局按 `&Layout` 借用、样本数在扫描时顺带计数、早停 |
| D11 | **P0** | `reed/src/filter.rs` `cut_block`/`emit_line`；`reed-io/src/bgzf/deflate.rs:55` | 块按"记录边界"切，单条记录 > 65,280 字节时块超限；`bsize as u16` 截断 | 宽队列输出损坏的 BGZF，下游 htslib 报错甚至误读 | 写出层与记录边界解耦：按 65,280 字节切块，记录可以跨块（htslib 行为）；`deflate_block` 对 `src.len() > 65_280` 返回错误而不是截断 |
| D12 | **P0** | `reed/src/filter.rs:348` | 只根据 `.gz` 后缀决定是否压缩 | `-o out.bcf` 写出文本 VCF；`out.vcf.bgz` 写出未压缩文本 | 用 `-O {v,z,b,u}` 显式指定，未指定时按扩展名推断（`.vcf.gz/.vcf.bgz` → z，`.bcf` → b，其他 → v）；不支持的格式报错 |
| D13 | **P0（语义）** | `reed/src/filter.rs:906` `parse_gt`；`policy_counters` 未赋值 | ploidy/multiallelic 策略在 v2 过滤路径中没有执行 | 多倍体被静默截断成前两个等位；非数字等位被当作缺失 | 统一 GT 解码器 + 策略在一处执行（03 §3） |
| D14 | P1 | `reed-io/src/vcf/sample_scan.rs:323` | `parse_allele` 允许 5 位数字，`u16` 乘法溢出 | `70000` 这样的下标在 release 下回绕成错误值 | `checked_mul/checked_add`，溢出报 `MalformedRecord` |
| D15 | P1 | `reed/src/filter.rs:756` | 掩码与 keep 的短掩码默认 keep（v1 A13 未修） | 目前不会触发（长度一致），但属于静默默认 | 用 `SampleMask`（长度强制相等） |
| D16 | P1 | `reed/src/filter.rs` `estimate_records` | 用第一个 chunk 的密度外推记录数，用于策略选择 | 文件头部记录宽度与平均值差异大时估计偏差很大；默认内存预算取物理内存一半而不是用户值 | 优先用索引（CSI/TBI 的 n_mapped）；没有索引时用前 K 个 chunk；运行时超出预算则溢写，而不是依赖估计正确 |
| D17 | P2 | `reed-io/src/vcf/line_edit.rs` | 掩码 GT 一律写 `./.`，单倍体应为 `.`；每条记录对 `mask_gt` 做一次 `count_ones`；`edit_info` 每条记录分配两个 `Vec` | 少见；小幅性能损失 | 按原 GT 的倍性写缺失；编辑器复用 scratch |
| D18 | P2 | `reed/src/filter.rs` `fill_info_set` | 新增 INFO 的追加顺序与 bcftools fill-tags 不同 | 正文不能与 gold 逐字节一致（值一致） | 采用 fill-tags 顺序，或在文档中声明顺序约定 |
| D19 | P1（流程） | `bench/suites/e2e/results.tsv` | M4/M5/M5b 验收数字未入日志；测试红着就提交 | 无法追溯 | 06-roadmap.md §3 流程纪律 |

---

## 4. 值得保留并推广的做法

- **oracle crate**：冻结 v0.2 作为差分基准，是这次重写能保证 M4 正确的关键。后续 legacy 删除时，oracle 继续保留（03 §2）。
- **worker 内 deflate + reducer 只 append**：性能的主要来源，结构正确。
- **按规则的漏斗与掩码归因**：报告质量超过三家参考工具。
- **plan 编译期按 header 校验类型**：应推广到 CLI 表达式（02-cli.md §4）。
- **bench 日志的 note 字段写明"门"和"对照"**：格式正确，只是需要坚持每个里程碑都写。
