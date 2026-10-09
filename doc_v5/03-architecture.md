# 03 架构：W0–W4 之后的建成形态与残留

> doc_v2/03 的目标架构已建成，doc_v4 的 W3（mmap reader）与 W4（BCF 原生过滤输入）把最后两处"多付的工作"去掉了。本文记录建成形态、本轮核实的偏差、crate 内部的历史残留（用户追问"crates 现在不乱么"的答案）与结构性修法。

---

## 1. 建成形态

### 1.1 crate 图（8 个，分层单向，无环）

```
reed-core     模型 + 计数内核（标量 / AVX-512 VPOPCNTDQ）+ 数学 + TagTable     依赖：bitflags、serde、thiserror、tracing
reed-io       BGZF、chunk 源（流式 / mmap）、VCF/BCF 解码与编辑、写出、索引   依赖：reed-core、noodles、libdeflater、memchr、memmap2、flate2、bstr
reed-engine   chunk 流水线、票据、字节许可、重排、中止、长度前缀链、StageMetrics 依赖：reed-core、reed-io、crossbeam-channel
reed-stats    统计量累加器 + summary（bcftools stats）                         依赖：reed-core、serde
reed-filter   Spec、表达式、plan、求值、sketch、funnel、provenance（无 noodles） 依赖：reed-core、reed-io、memchr、serde、sha2
reed          驱动：filter / stats / index、统一输入入口、BCF 原生输入        依赖：以上 + rayon、noodles、bstr、serde_json
reed-cli      二进制 `reed`：clap、单位、翻译、verbosity                       依赖：reed、reed-engine、reed-filter（Cargo.toml 多声明了 core/io/stats，见 §4 C-8）
reed-core-v02-oracle  冻结的 v0.2（仅 dev-dependency）                       → 改名 crates/oracle（§4 C-9）
```

体量（src 行）：reed 7.5k、reed-io 7.0k、reed-filter 4.2k、reed-stats 4.0k、reed-core 4.0k、oracle 2.9k、reed-cli 1.9k、reed-engine 1.5k。

### 1.2 关键机制（已核实）

| 机制 | 位置 | 状态 |
|---|---|---|
| BAPP 块对齐并行解析 | `reed-io/src/chunk/bgzf.rs` | reader 只读块头切 ~1 MiB 压缩 chunk；W3 后文件输入走 **mmap**（`MmapChunkSource`，chunk = 映射区间，`ChunkPayload::{Owned, Mapped}`），非普通文件回退流式 |
| 票据 + 字节许可 | `reed-engine/src/{pipeline,budget}.rs` | 在途 ≤ tickets；在途峰值随线程线性（56 线程 912 MB） |
| 长度前缀链（BCF 定界） | `reed-engine/src/frame_chain.rs` | 先发布后解码、定向唤醒、毒化、**全局静默 120 s** 看门狗（W2/F4-14） |
| 统一记录管线 | `reed/src/filter/record_pipeline.rs` | `process`（文本）与 `process_scanned`（原生 codes/列）共用 `run_stages`；single / two-pass / sketch 与 `FilteredJob` 共用 |
| BCF 原生过滤输入（W4） | `reed-io`：`render_site_line` + `decode_format_columns`；`reed/src/filter/mod.rs:191-285,506-541` | 判定不渲染；只有存活记录才 `render_record` → `edit_vcf_line` → 写块；适用条件：BCF 输入 ∧ 掩码无 `--set-gts .` FMT 表达式 ∧ 计划内 FORMAT 列均 Number=1，否则回退渲染路径 |
| 边过滤边统计 | `reed/src/filtered.rs` | `FilteredJob::with_records` 把掩码后的 `GtCode` 直接交内层（W2/O-2） |
| 报告一次读取 | `reed/src/stats/job.rs` | `ReportRequest.summary`：summary 累加器并入 report partial（W2/O-1） |
| 输出层 | `reed-io/src/write/{block_writer,ordered_output,index}.rs` | 记录可跨块、按序 append、`-W` 才建索引（W1 修正） |
| BCF 编码 | `reed-io/src/bcf/native_encode.rs` | 原生编码器；noodles 编码器仍在生产路径作可选项（§4 C-4） |
| 计时 | `reed-engine/src/metrics.rs` | `StageMetrics`：chunks_read、reader、reduce、chain_wait、worker_busy、在途峰值（W2/O-5） |

### 1.3 reader 之后的瓶颈位置

W3 把 reader 从 2.2 s 压到 0.35 s 之后，stats 的 wall 由 CPU 决定（五统计 1.7 s / 1.5 s）；filter 的 wall 由 deflate 决定（6 级 ≈ 56% CPU）。剩下的"多付的工作"只有一处：`--summary` 用自己的扫描器再解一遍 GT（§5）。

---

## 2. `unsafe` 策略（本轮核实出的偏差）

工作区 `unsafe_code = deny`。实测有 **两处** 允许，不是 STATUS 说的一处：

| 位置 | 用途 | 评审 |
|---|---|---|
| `reed-core/src/count/kernel.rs:16` `#![allow(unsafe_code)]`（模块级，6 个 `unsafe` 块） | AVX-512 `VPOPCNTDQ` 计数内核，运行时分派，标量路径作 parity oracle | 自 doc_v1 §14.4 起就是"受批准的例外"，每块有 `// SAFETY:` |
| `reed-io/src/chunk/bgzf.rs:361` `#[allow(unsafe_code)]`（单语句） | `memmap2::Mmap::map` | W3 新增；契约：只读映射、生命周期覆盖全部 chunk、输入不被写入、非普通文件回退流式 |

处理：`bgzf.rs:356` 的注释 "the single audited unsafe in the workspace" 改为 "one of the two audited `unsafe` sites（另一处 `reed-core/src/count/kernel.rs`）"；`doc_v4/STATUS.md` 的偏差行同步；`scripts/check-doc-drift.sh` 加一条守卫：`grep -rl 'allow(unsafe_code)' crates --include=*.rs | wc -l` 必须等于 2（新增 unsafe 必须先改文档再改守卫）。

---

## 3. 残余差距

| # | 目标 | 现状 | 处理 |
|---|---|---|---|
| A5-1 | reed-filter 只依赖 reed-core（doc_v2 03 §7 L3） | 仍依赖 reed-io（`LineRef`、`edit_vcf_line` 在 eval 中使用） | 维持（方向正确、无环），不作为阻塞项 |
| A5-2 | STATUS 的"验证"列引用测试名或命令 | 做到；但"唯一 unsafe"的陈述与代码不符（§2） | X1 改文档与守卫 |
| A5-3 | 代码注释与规范同形 | `reed-filter/src/funnel.rs:67,72` 注释仍写 `sample.max_missing`（常量已是 `sample.mind`）；`reed-stats/src/lib.rs:47` 注释写 `--fst-pops`；stats `--tmp-dir` 帮助写 "read twice" | X1 一次 grep 清理，并把 `max_missing`/`fst-pops` 加入守卫的禁词（排除撤回参数的报错函数） |
| A5-4 | 工作树可编译 | 2026-10-09 18:08 实测：`cargo test --workspace` **编译失败**（`reed-stats/tests/accumulator_laws.rs:60-118` 的 `match` 未覆盖未提交的 `StatAcc::ProjectedSfs`；`git status` 显示 `reed-stats/src/{request,sfs}.rs` 有 W5 在途改动）。HEAD（c5b6bfe）单独验证见 01 §1 | X0：在途改动要么补齐测试臂并提交，要么 stash；主干不得留在"测试不编译"的状态 |
| A5-5 | reports 两族 ≤ 3.2 s（W2 门） | 3.58 s，CPU 181 s；`--summary` 的 `qc_report.rs:153` 自建 `GenotypeScanner::tracking_alleles()` 再扫一遍 | §5 的结构修法（X3） |

---

## 4. `crates/` 内部的历史残留（C-1…C-9）

分层不乱；乱在两代并存、名字漂移、散件三类，全部核实：

| # | 位置 | 问题 | 处理 |
|---|---|---|---|
| C-1 | `reed/src/stats/writers/qc.rs`（146 行，"M3 first consumer"） | 与 `qc_report.rs`（759 行）是两代；只在 `writers/mod.rs:5` 被声明，无调用者 | 删除 |
| C-2 | `reed-cli`：`cli.rs:49` `QcCli`、`translate/qc.rs`、`main.rs:25` `Command::Qc` | 隐藏的 deprecated `reed qc`，注释指向已合并的 reed-qc | 删除；`--summary` 已覆盖 |
| C-3 | "qc" 一词 | CLI `--summary`、库 `reed_stats::summary`，驱动与产物却叫 `qc_report.rs`、`.qc.json` | `qc_report.rs` → `stats/writers/summary.rs`；`.qc.json` → `.summary.json` |
| C-4 | `reed-io/src/bcf/encode.rs`（noodles 编码器） | B4 已被 `native_encode.rs` 取代，但 `reed/src/filter/single_pass.rs:31` 仍以 `Option<BcfEncoder>` 带在生产路径 | 生产路径只留 `NativeEncoder`；noodles 版只作 `reed-io/tests/bcf_encode.rs` 的 oracle |
| C-5 | `reed-stats/src/individual.rs` | 一个文件两个统计量（individual F、joint SFS） | 拆 `het.rs` 与 `sfs/{projected,joint}.rs`（与 W5 的 SFS 接入同一 PR） |
| C-6 | `reed/src/{filtered.rs, bcf_input.rs, input.rs}` 散在顶层 | 顶层应只有三条命令驱动 + 输入层 | `input/{mod.rs, bcf.rs}`；`filtered.rs` → `stats/filtered.rs` |
| C-7 | `reed/src/stats/job.rs` 1,167 行 | 建请求、跑、写三件事一个文件 | 拆 `stats/{request.rs, run.rs, sink.rs}`，与 `filter/{single_pass,reduce,emit}` 对称 |
| C-8 | `reed-cli/Cargo.toml` | 声明 reed-core/io/stats，源码只 `use reed_engine`、`reed_filter` | 删未用依赖 |
| C-9 | `reed-core-v02-oracle` | 版本后缀 | `crates/oracle`，包名 `reed-oracle`，`publish = false` |

不改：`reed-core/src/grid.rs` 与 `reed-stats/src/window.rs` 是分层不是重复；`reed-filter/src/{sketch,funnel}.rs` 放顶层可接受；`reed-io` 的 `bgzf/ chunk/ vcf/ bcf/ write/` 分法清楚。

顺序与门：C-1/C-2/C-8/C-9 零风险，随目录搬迁（06）同一 PR；C-3/C-5/C-6 只改名归属，门是测试全绿 + 全部输出 md5 不变；C-4/C-7 触及生产路径，各自单独 PR，门是 gold 正文 md5 与 `-Ob` 正文 sha256 不变。

---

## 5. reports 门的结构修法（X3）

现状：`--freq/--hardy/--het/--missing/--depth` 共用一次 `GenotypeScanner` 扫描产出 `ScanOutput`；`--summary` 在同一 job 里**再建一个** `tracking_alleles()` 扫描器重扫（`qc_report.rs:153`），为的是多等位下标。全量 CPU 181 s 中约 25–30 s 是这次重扫。

修法：
1. reports job 在请求含 `--summary` 时一开始就用 `tracking_alleles()` 的扫描器（多一个下标数组，成本远小于重扫）；
2. summary 累加器改为消费 `ScanOutput`（`QcAccumulator::push_scan(&ScanOutput, &SiteView)`），删掉它自己的 GT 解析；
3. 门：8 个报告文件 md5 不变；全量 ≤ 3.2 s、CPU ≤ 155 s（负载 < 1，3 次中位数）。

---

## 6. 守卫扩展（`scripts/check-doc-drift.sh` → 搬迁后 `ci/`）

新增：
- `allow(unsafe_code)` 恰好 2 处（§2）；
- 代码里不再出现 `sample.max_missing`、`--fst-pops`（排除撤回参数的报错函数与测试里对 vcftools 的调用）；
- `reed qc` 不存在（`! grep -q 'Command::Qc' crates/reed-cli/src/main.rs`，C-2 后生效）；
- `cargo test --workspace` 编译通过是 CI 的第一条（A5-4 的教训）。

---

## 7. 前后端边界（08 §6 的架构落点）

```
前端  reed-cli ──只依赖──► 门面 reed ──► reed-filter / reed-stats / reed-engine / reed-io / reed-core   ← 后端
      （唯一前端；reed-py 已裁决不做）
```

三条硬规则：**前端只依赖门面**；**后端不打 stderr、不认识 clap、不定退出码**（警告走 `tracing::warn!`，错误带 `kind()`）；**门面只 `pub use`**（§6.2 的清单，不 `pub mod` 实现模块）。现状违规（前端 20 处直接拿下层 crate、后端 5 处 `eprintln!`、`plan::compile`/`assign_groups` 落在前端）与守卫见 08 §6.3–6.4，动作并入 X1。
