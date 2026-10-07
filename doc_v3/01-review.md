# 01 第三轮复审（R0–R8 落地后，2026-10-08）

> 复审对象：`28bc26d`（R0）… `77a5cdd`（R5 FMT 表达式）共 26 个提交，以及 `doc_v2/STATUS.md` 的自述。
> 方法：逐模块阅读、`cargo build/test/clippy`（测试连续运行 3 次）、全量数据对照实验（本机 56 核，热页缓存，bash `time`）。
> 本机没有 `perf`，只有 valgrind；CPU 去向靠"同一命令只改一个变量"的对照实验分离。

---

## 0. 总评

实现质量高：doc_v2 的 P0 缺陷全部修复并有回归测试；R2（统一解码）、R3（真正的票据、字节背压、单一代码路径）、
R4（输出层）、R5（flag CLI、bcftools 表达式子集含 `FMT/` 与 `--set-gts .`）、R7（π/dxy/θ/Tajima/Fst 与 pixy 逐位一致、
`--variant-only` 与 vcftools 一致、PBS 与 allel 一致、plink2 格式报告）都已落地，`STATUS.md` 的记录方式值得保持。

剩下的问题分三类：

1. **R6 的 BCF 部分没做**，方案见 `02-bcf.md`。
2. **几个串行瓶颈和资源问题**：报告类统计在 reducer 里串行计算，未压缩输出的分配抖动，`reed index` 单线程，stdin 整体落盘且不删除。
3. **两个语义缺陷**：`--mind` 的实现与声明的 plink2 语义不符；样本规则互斥导致参数被静默忽略。
4. **对齐裁决的残留**：`##reed_ascertainment` 仍在写出、reed-stats 仍导出已撤下的方法、reed-filter 仍依赖 noodles、`reed qc` 与 `stats --summary` 两个入口，以及 CLI 与 doc_v2 规范的若干出入。

---

## 1. 状态核实

| 项 | 结果 |
|---|---|
| 构建 | `cargo build --release --workspace` 通过；reed-cli 有 17 条 rustc 警告 |
| 测试 | `cargo test --release --workspace` 连续 3 次全绿（doc_v2 时期的偶发失败已消失） |
| clippy | 无 `clippy::all` 错误；pedantic 约 700 条警告（`expect()` 404 条） |
| 正确性 | 沿用 STATUS 的对拍记录；本轮复测 gold 等价命令：3,315,604 → 1,768,759，与 gold 一致 |
| 跨格式数据 | `clean.bcf` 与 `clean.vcf.gz` 记录数相同（3,315,604），抽查的首 2,000 条与 Chr20 区间经 `bcftools view` 渲染后逐字节一致，可直接作为 BCF 的全量验收数据 |

---

## 2. 性能实测（2026-10-08，56 线程，热缓存）

| 命令 | wall | user CPU | 说明 |
|---|---|---|---|
| `cat clean.vcf.gz > /dev/null` | 0.65 s | — | 读取不是瓶颈（约 10 GB/s） |
| filter（gold 规则，`-S keep`，单遍）`-Oz6` | **8.8 s** | 459 s | 默认级别 |
| 同上 `-Oz1` | 5.6 s | 288 s | 输出 4.88 GB（级别 6 为 3.90 GB） |
| 同上 `-Ov -o /dev/null` | 5.8 s | 202 s | 不压缩 |
| filter `--mind 0.2`（sketch-then-emit） | 16.2 s | 492 s | 两遍 |
| filter 只有 `--geno/--maf`（全部 3,315,604 条存活）`-Oz1` | 5.5 s | 270 s | |
| 同上 `-Ov -o /dev/null` | **13.2 s** | 116 s | 8 线程时 18.6 s：未压缩路径的串行瓶颈（F3） |
| stats `--pi` | 2.6 s | 72 s | GT-only |
| stats `--summary` / `reed qc` | 2.8 / 2.9 s | 104 s | |
| stats `--freq` / `--het` / `--missing` | 2.8 / 2.8 / 2.9 s | 约 105 s | 每项的 CPU 都相同：所有报告总是全部计算（F2） |
| stats `--hardy` | **13.5 s** | 116 s | 串行瓶颈（F1） |
| `reed index`（3.9 GB 输出） | **31.0 s** | 30 s | 单线程（F4）；bcftools 79.7 s |

由此得到：

- 并行路径是 **CPU 受限**的：wall ≈ 总 CPU ÷ 54 个 worker（-Oz6：459 / 54 = 8.5 s，实测 8.8 s）。所以后续优化看的是每条记录的 CPU。
- **deflate 6 级占 filter 总 CPU 的 56%**（459 − 202 = 257 CPU-s）。级别 6 与 htslib/bgzip 的默认一致，不应改默认值（03-performance.md O8）。
- 不含压缩时，filter 的 CPU 约为 175 ns/样本·记录（含 inflate、GT+DP+GQ 扫描、掩码、编辑）；GT-only 的 stats 约为 63 ns（含 inflate）。

---

## 3. 发现

| # | 级 | 位置 | 问题 | 证据 | 修法（详见） |
|---|---|---|---|---|---|
| F0 | **P1** | `reed/src/filter/{reduce,sketch}.rs`；`reed-cli/src/cli.rs` 的 `--mind` 帮助 | `--mind` 宣称 plink2 语义（掩码后、载入期过滤后的位点），实际统计的是**全部记录上的原始缺失**（gold 口径）。同时使用掩码参数时，保留的样本与 plink2 不同 | `reduce.rs:60-70` 对每条记录、掩码前计数；plink2 `MindFilter`（`plink2.cc:2061`）在导入掩码之后、`--geno`（`:2707`）之前 | 按 plink2 实现，gold 继续用 `-S`（04 C0） |
| F0b | **P1** | `reed-cli/src/translate/filter.rs:23-36` | 样本规则是互斥的 `if/else if`：`-S keep.txt --mind 0.2` 时 `--mind` 被**静默忽略**，`-S` 与 `--keep` 同时给出时 `--keep` 被忽略 | 代码；`spec.cohort` 只能容纳一条规则 | 静态子集与 `--mind` 分开建模，冲突参数报错（04 C0b） |
| F1 | **P1** | `reed/src/reports.rs` `ReportSink::write_site_rows`（在 `accept` 中调用） | `.hardy` 的 HWE 精确检验与所有逐位点报告的格式化都在 reducer 上串行执行；每行还有多次 `String` 分配（`from_utf8_lossy(...).into_owned()`、`format!`） | `--hardy` 13.5 s，其他报告约 2.8 s，CPU 只多 11 s | worker 计算并格式化成字节，reducer 只做按序 append（03 O1） |
| F2 | P1 | `reed/src/reports.rs` `ReportJob::push_line` | 不论请求了哪项报告，worker 都计算全部（深度、个体 F、位平面、逐位点行，以及每行 4 次 `to_vec()`） | 任意单项报告都是约 105 CPU-s | 按请求装配（03 O2） |
| F3 | **P1** | `reed-io/src/write/block_writer.rs`（`raw` 每 chunk `mem::take`）、reed-cli 未启用全局分配器 | 未压缩输出时，每个 chunk 约 5 MB 的输出缓冲从容量 0 反复扩容，并在 reducer 线程跨线程释放（大块走 mmap/munmap）；workspace 声明了 mimalloc，但二进制没有启用（v1 决策 D7） | 全部存活时 `-Ov` 13.2 s，`-Oz1` 只要 5.5 s；8 线程时更慢（18.6 s） | 缓冲池（reducer 写完后把缓冲还给 worker）+ 预分配 + 启用 mimalloc（03 O3） |
| F4 | P1 | `reed/src/index.rs::run_index` | 单线程：noodles BGZF reader 加逐行 `read_until`；每条记录线性查找 contig、分配 `String` | 3.9 GB 用 31 s，6.6 GB 用 48.7 s | 用 BAPP 并行：worker 从块布局算出每条记录的虚拟偏移，reducer 按序调用 `Indexer::add`（03 O4） |
| F5 | **P1** | `reed/src/input.rs::stdin_spill` | stdin 一律整体落盘到 `/tmp/reed-stdin-PID.vcf`，**从不删除**（本机 /tmp 已残留 5 个）；必须等上游写完才能开始（管道失去流式意义）；根分区只剩 113 GB，而本数据解压后有 35 GB | 代码 + `/tmp` 现场 | 单遍任务直接流式读取 stdin；只有多遍策略才落盘，落到 `--tmp-dir`，用 RAII 保证退出（包括出错）时删除（03 O5） |
| F6 | P2 | `reed-filter/src/plan/{ascertainment,provenance}.rs`、`reed-core/src/manifest.rs`；`reed/tests/e2e_filter.rs:325` | `##reed_ascertainment` 仍在写出，清单 v2 类型仍在 core，测试还在断言它；doc_v2 04 §1 已裁决删除 | grep | 删除并改测试（04 C1） |
| F7 | P2 | `reed-stats/src/lib.rs` | 仍然导出已撤下的方法（`NeiMultiAcc`、`HudsonFiniteAcc`），以及没有 CLI 出口、没有参照物对拍的累加器（LD、Patterson f 统计、jackknife、个体 F/联合 SFS、投影 SFS、`FstMatrixAcc`） | `pub use` 列表 | 删除已撤下的；其余在接入 CLI 并通过对拍之前设为 `pub(crate)` 或移入 feature（04 C2） |
| F8 | P2 | `reed-filter/src/plan/mod.rs` | 仍依赖 noodles 的 header 类型，L3 纯策略层的目标没有达成（03 §7） | grep | 在 reed-core 定义 `TagTable`（04 C3） |
| F9 | P2 | `reed qc` 与 `reed stats --summary` | 两个入口做同一件事；reed-qc crate 没有并入 reed-stats（doc_v2 02 §1、03 §7） | CLI help | 只保留 `stats --summary`（`qc` 作隐藏别名一个版本后删除）（04 C4） |
| F10 | P2 | `reed-cli/src/cli.rs` | CLI 与 doc_v2 规范的出入（完整列表见 04-cleanup.md §2）：`-g/--groups` 为主名而 `--pop` 只是隐藏别名；stats 没有过滤参数与 `-r/-t/-s`；`--pbs` 不接受群体参数；`--fst-pops/--sfs/--proj/--bed-file/--sites-file` 缺失；`index -t` 的含义与 bcftools 相反；`--compress-level` 与 `-Oz[0-9]` 重复；多处 help 文本过时 | `--help` 输出 | 04 §2 |
| F11 | P2 | `reed-engine/src/error.rs` | 类型化 `EngineError` 已经两次推迟（doc_v2 03 §4.4）；BCF 会引入大量新的错误来源 | STATUS 偏差 1 | 在 BCF 之前完成（02-bcf.md B0） |
| F12 | P2 | `bench/e2e/results.tsv` | 多数基准行在有未提交改动的工作区上测得（`+dirty`），无法复现；同一命令不同行的口径不清（`window-scan` 21.2 s 是哪个统计？） | 日志 | 只在干净提交上测；`dirty` 时记录 `git diff` 的哈希；`command` 列写完整命令（03 §3） |
| F13 | P3 | `reed-io/src/chunk/bgzf.rs::read_block`、`chunk/plain.rs` | 每块分配一个 `Vec` 并零填充，再拷贝进 chunk；plain 源再 `Arc::from(slice)` 拷贝一次 | 代码 | 文件输入用 mmap 零拷贝切 chunk；不紧急（读取不是瓶颈） |
| F14 | P3 | `reed-io/src/vcf/genotype_scan.rs::scan` | 每个样本调用多次 `memchr`（找 `\t` 一次，切子字段每个 `:` 一次）；对 30 字节左右的短字段，函数调用开销占主导 | 175 ns/样本 | 单遍字节循环同时跟踪 `:` 与 `\t`，以微基准为门（03 O6） |
| F15 | P3 | `reed-io/src/vcf/line_edit.rs::edit_samples` | 逐样本 `extend_from_slice`；未被掩码的连续样本可以一次拷贝 | 代码 | 按"连续未编辑样本段"批量拷贝（03 O7） |

---

## 4. 值得保留的做法

- `doc_v2/STATUS.md`：每项都写明 commit、验证方法与数值，偏差单独登记。doc_v3 沿用这一格式。
- 对拍覆盖面：pixy（最大相对误差 1e-14 级）、vcftools（`--window-pi`）、allel（PBS）、plink2（报告）、bcftools（表达式与 `FMT/` 掩码），每项都有数字。
- 所有 P0 缺陷都先写失败用例再修复。
- 对参照物缺陷的处理：bcftools 1.24 的 `!~` 与 `~` 行为相同，与手册相反，reed 按手册实现并在 STATUS 中记录。这是对齐原则的正确用法。
