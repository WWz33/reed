# 03 性能优化计划（基于 2026-10-08 实测）

> 前提（01-review.md §2）：并行路径 CPU 受限，wall ≈ 总 CPU ÷ worker 数；读取不是瓶颈。
> 因此优化分两类：①消除串行段（直接降 wall）；②降低每条记录的 CPU（按比例降 wall）。
> 每项都要求前后对照入 `bench/e2e/results.tsv`，没有收益就撤回（v1 §17.4）。

---

## 1. 消除串行段（P1，先做）

### O1 报告类统计在 worker 中计算与格式化（F1）

- **现状**：`ReportSink::write_site_rows` 在 reducer 中逐行计算 HWE、渲染浮点、分配 `String`。`--hardy` 13.5 s。
- **改法**：worker 在 `push_line` 中直接把 `.afreq/.hardy/.vmiss/.ldepth.mean` 的行渲染进各自的字节缓冲（每个输出文件一个），
  HWE p 值在 worker 中计算；partial 只携带这些字节缓冲与逐样本累加量。reducer 按 chunk 序 `write_all`，与 filter 的写出方式相同。
  格式化使用零分配的整数/浮点渲染（`itoa` 式手写、`g6` 改为写入 `&mut Vec<u8>`）。
- **预期**：`--hardy` 13.5 s → 约 3 s（与 `--freq` 相同量级）。
- **门**：四个报告文件与改动前逐字节相同；wall 入日志。

### O2 报告按请求装配（F2）

- **现状**：任意单项报告都计算深度、个体 F、位平面与逐位点行，约 105 CPU-s。
- **改法**：`ReportRequest` 决定 `FieldPlan`（不要 DP 就不解 DP）与计算项（不要 `--het` 就不建位平面、不跑 `IndividualFAcc`）；
  逐位点行只在需要时生成，且不再 `to_vec()` CHROM/ID/REF/ALT（O1 后行在 worker 内直接渲染，不需要保留）。
- **预期**：`--freq` 单项约 75 CPU-s（接近 `--pi` 的 72）。
- **门**：输出逐字节不变。

### O3 输出缓冲池 + mimalloc（F3）

- **现状**：未压缩输出时每个 chunk 新建一个 `Vec` 并在 reducer 释放；全部存活时 `-Ov` 13.2 s（`-Oz1` 只要 5.5 s），8 线程 18.6 s。
- **改法**：
  1. 缓冲池：引擎提供一个有界的空缓冲回收通道，reducer 写完一个 partial 后把其中的大缓冲（plain 输出字节、压缩块、报告行缓冲）还回去，
     worker 从池中取缓冲（`clear()` 后复用，容量保留）。稳态下零分配。
  2. 预分配：`BlockWriter` 的 plain 缓冲按"上一个 chunk 的输出量"预留容量。
  3. 在 reed-cli 启用 `mimalloc` 作为全局分配器（v1 决策 D7，workspace 已声明依赖但未启用），以吸收其余的跨线程释放。
- **预期**：全部存活的 `-Ov` 13.2 s → ≤ 6 s；其他路径小幅收益。
- **门**：输出逐字节不变；分别测"只开池"与"池 + mimalloc"，各自无收益则撤回。

### O4 并行 `reed index`（F4）

- **现状**：单线程，3.9 GB 输出 31 s，6.6 GB 输入 48.7 s。
- **改法**：
  - reader 切 chunk 时已经解析每个块的 BSIZE；再从块尾读出 ISIZE（不需解压），把每块的
    `(压缩偏移, 解压长度)` 记入 `ChunkBytes::blocks`，并记录 chunk 在文件中的压缩起点。
  - worker：对 chunk 内每条完整记录，用块的解压长度前缀和把"记录起点在 chunk 文本中的偏移"换算成虚拟偏移
    `(chunk 起点 + 块压缩偏移) << 16 | 块内偏移`，与 `(tid, pos, end)` 一起放进 partial。
  - 接缝记录：起点在前一个 chunk 的尾片段中，其虚拟偏移由前一个 chunk 的 worker 计算并随 partial 带出。
  - reducer：按序调用 `Indexer::add`（330 万次，廉价），计算每条记录 chunk 的终点 = 下一条记录的起点。
  - contig 查找用"上次命中"缓存，不再逐条分配 `String`。
- **预期**：全量 ≤ 3 s。BCF 输入（02-bcf.md B1）用同一机制。
- **门**：CSI/TBI 与现有实现逐字节相同；bcftools 用该索引查询的结果与流式 `-t` 一致。

### O5 stdin 流式读取（F5）

- **现状**：整体落盘到 `/tmp/reed-stdin-PID.vcf`，从不删除，并且要等上游写完才能开始。
- **改法**：
  - 单遍任务（`stats`、单遍 `filter`、`index` 以外的所有命令）：`BgzfChunkSource`/`PlainChunkSource` 直接包装 stdin（它们本来就只需要 `Read`），
    格式探测读入的前缀用 `Read::chain` 接回去。
  - 多遍任务（`filter --mind`、两遍策略）：才需要回放。落盘到 `--tmp-dir`（默认 `$TMPDIR`），文件名带进程号与随机后缀，
    用 RAII 守卫在正常结束、出错和 Ctrl-C（`ctrlc` 处理）时删除；只落**原始压缩字节**（与现在相同）。
  - `-r/-R` 与 stdin 同时出现时报错（已有）。
- **门**：`bcftools view x.vcf.gz | reed stats - …` 的结果与直接读文件一致；进程结束后 `/tmp` 中没有残留；单遍任务不产生任何临时文件。

---

## 2. 降低每条记录的 CPU（P2，以微基准为门）

### O6 单遍样本列扫描

- **现状**：每个样本调用一次 `memchr` 找 `\t`，再为每个需要的子字段调用 `memchr` 找 `:`（`split_needed`）。
  30 字节左右的短字段上，调用开销占主导；filter 约 175 ns/样本，stats（GT-only）约 63 ns/样本（均含 inflate）。
- **改法**（逐级尝试，每级都要微基准证明收益）：
  1. 单个字节循环同时跟踪 `:` 计数与 `\t`，读到 `last_needed` 之后的子字段时转用 `memchr` 跳到 `\t`；
  2. 一次性用 SIMD（AVX2 / AVX-512，运行时分派）为 64 字节块生成 `\t` 与 `:` 的位图，按位遍历字段边界（simdjson 的结构索引思路）；
  3. GT 快路径保留 3 字节查表。
- **基准**：新增 `benches/genotype_scan.rs`，输入取 `clean.vcf.gz` 前 1 万行，测 GT-only 与 GT+DP+GQ 两种计划的 ns/样本。
- **预期**：扫描部分 1.5–3 倍。filter 的 wall 主要由 deflate 决定（O8），所以这一项主要惠及 stats、`-Oz1/-Ov`、以及核数较少的机器。
- **门**：差分测试与 fuzz 不变；全量输出逐字节不变。

### O7 编辑时按连续样本段拷贝

- **现状**：`edit_samples` 逐样本 `extend_from_slice`。
- **改法**：遍历样本时只记录"需要改动的样本"（被掩码或被剔除）的位置，其间的连续未改动区段一次拷贝。
  本数据约 14% 的基因型被掩码，平均每段约 7 个样本，`extend` 调用约减为 1/7。
- **门**：输出逐字节不变；`benches/line_edit.rs`。

### O8 压缩级别（不改默认值）

- 实测：deflate 6 级占 filter 总 CPU 的 56%；`-Oz1` 的 CPU 少 37%，输出大 25%。
- 默认级别保持 6（与 bgzip/htslib 默认一致，属于对齐范围）。在文档与 `--help` 中说明：
  中间文件用 `-Oz1` 或 `-Ou`（BCF 落地后），最终交付用默认级别。
- 不引入其他压缩库：libdeflate 在 6 级已是同类最快之一，换库收益不确定，不做。

### O9 sketch-then-emit 的第二遍

- **现状**：`--mind` 时 16.2 s，单遍 8.8 s。
- **先测量**：实现 `--timing`（O10）后拆出 pass0、中段、pass1 各自的耗时，再决定以下两项是否值得做：
  1. 中段（`mid_stage`）用 rayon 按 chunk 并行（doc_v2 03 §9，未做）；
  2. 投机单遍：pass0 同时按"全部样本保留"产出输出；若 `--mind` 最终没有剔除任何样本（常见情形），直接采用，省去 pass1；
     否则丢弃并运行 pass1。按"剔除样本的概率 × pass1 成本"决定是否默认开启。

### O10 计时与剖析工具

- 实现 `--timing FILE`（doc_v2 02 §2 已列出，未实现）：输出各阶段的忙碌/空闲时间、在途字节峰值、各 pass 耗时（引擎的 `StageMetrics` 已有基础）。
- 本机 `perf_event_paranoid = 2`，允许剖析自己的进程：建议 `cargo install samply`，生成火焰图作为每个 P2 优化的前后证据。

### O11 读入零拷贝（F13，低优先级）

文件输入用 `memmap2` 映射，`ChunkBytes` 持有 `Arc<Mmap>` 加区间，reader 只解析块头。读取目前不是瓶颈，只在其他项完成后再做。

### O12 HWE 记忆化（低优先级）

O1 之后若 `--hardy` 或 `--hwe` 仍可见，按 v1 §9.6 以 `(n_called, mac)` 为键缓存整条概率分布。数值必须与现实现逐位一致。

---

## 3. 基准日志协议补充（F12）

- 只在干净提交上测；必须在有改动的工作区测时，`dirty` 列写 `git diff | sha256sum` 的前 12 位。
- `command` 列写完整命令行，`note` 写"测的是哪个统计/哪个阶段"。
- 每次至少两次运行取中位数；页缓存状态写明（本机 /sri 是 dpc 并行文件系统，冷读与热读差异很大）。

---

## 4. 优化后的目标（本机全量，56 线程，热缓存）

| 命令 | 现在 | 目标 | 依赖 |
|---|---|---|---|
| `stats --hardy` | 13.5 s | ≤ 3.5 s | O1 |
| `stats --freq --hardy --het --missing` | 15.0 s | ≤ 4 s | O1、O2 |
| filter 全部存活 `-Ov` | 13.2 s | ≤ 6 s | O3 |
| `reed index`（6.6 GB） | 48.7 s | ≤ 3 s | O4 |
| `bcftools view … \| reed stats -` | 先整体落盘再计算 | 与读文件相同 | O5 |
| filter gold 规则 `-Oz6` | 8.8 s | ≤ 8.8 s（不回退） | — |
| stats `--pi` | 2.6 s | ≤ 2.6 s；O6 后 ≤ 2.0 s | O6 |
