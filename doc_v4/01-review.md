# 01 第四轮复审（V0–V6 落地后，2026-10-09）

> 复审对象：`b1b39e2`（V0）… `61f24a2`（BCF 原生记录路径）共 21 个提交，以及 `doc_v3/STATUS.md` 的自述。
> 方法：逐模块阅读、`cargo build/test/clippy`（测试连续 3 次）、全量与子集实测（56 线程；每次记录 1 分钟负载）、
> 三工具边界对拍（vcftools、plink2、bcftools）。本机无 `perf`；CPU 去向靠 `--timing`（O10 已落地）与对照实验。

---

## 0. 总评

doc_v3 的路线基本全部落地，质量与 v3 时一致：每项都有 commit、数字和回归测试；BCF 输入/输出与 VCF 路径的结果逐字节一致；`--timing` 让性能分析有了依据。

本轮的问题分三类：

1. **参数语义方向不一致**（用户本轮的首要问题）：`--geno`、`--mind` 用"丢弃条件"，其他阈值用"保留条件"。证据与修法见 02-cli-semantics.md。
2. **文档与代码的漂移**：STATUS 说 C3 未做，代码已做；`--verbosity` 被接受但在 stats 中无效；`reed index -` 的帮助说支持 stdin，实际报 "No such file"。
3. **若干可优化的实现点**：BCF 输入的 filter 走"渲染→文本→重扫"（CPU +38%）；`stats` 两族报告同时请求时读两遍；`--timing` 没有 reader 阶段的分解。

---

## 1. 状态核实

| 项 | 结果 |
|---|---|
| 构建 | `cargo build --release --workspace` 0 警告 |
| 测试 | `cargo test --release --workspace`：**344 通过 / 0 失败**，连续 3 次 |
| clippy | 0 error；剩余 pedantic 提示（`format!` 追加、naive bytecount、`u64→u32` 截断等）按 doc_v3 C7 的豁免政策保留 |
| crate | 8 个：reed-core / reed-io / reed-engine / reed-stats / reed-filter / reed / reed-cli / reed-core-v02-oracle；reed-qc 已并入 reed-stats |
| 规模 | 约 38,300 行（含 oracle 约 3,000 行、测试约 5,000 行） |
| 测试文件 | `reed/tests/{determinism,e2e_bcf,e2e_filter,e2e_index,e2e_pixy,e2e_qc,e2e_reports,e2e_stats}.rs`、`reed-cli/tests/stdin_streaming.rs`、reed-io/reed-stats/reed-engine/reed-core 各自的集成测试 |

---

## 2. 正确性核实（2026-10-09）

| 检查 | 结果 |
|---|---|
| gold 等价命令（VCF.gz 输入，`-Oz6`） | 3,315,604 → 1,768,759（与 gold 一致） |
| `-Ob` 正文 vs `-Oz` 正文 | `bcftools view -H` 正文 md5 相同 |
| BCF 输入 → `-Oz` 正文 vs VCF 输入 | md5 相同 |
| 子集 BCF→BCF（`-Ob`） vs VCF→VCF | md5 相同 |
| stats 五统计（π/dxy/θ/D/Fst）BCF vs VCF 输入 | `_pi.txt` md5 `45c2f2e1…` 相同（其余四文件同） |
| `-Ou` | 0 级 BGZF（首块 22,153 B 压缩 / 22,122 B 解压），`bgzip -t` 通过，bcftools 可读 |
| 掩码后缺失阈值边界（subset，`--min-dp 5 --min-gq 20`） | reed `--geno 0.1` = vcftools `--max-missing 0.9` = plink2 `--geno 0.1` = **263,420**；`--geno 0.2` = vcftools `0.8` = **381,672** |
| 半调用 | Chr01:1–300 kb：0 个半调用（351,980 已调用 / 2,604 缺失），VCF 与 BCF 一致 |
| stdin | 单遍 stats 从管道读 4.45 s，输出与文件一致，`/tmp` 无残留 |

---

## 3. 性能实测（2026-10-09，56 线程，热缓存，1 分钟负载 0.3）

### 3.1 全量（clean.vcf.gz 7.0 GB / clean.bcf 8.1 GB，3,315,604 条）

| 命令 | wall | user CPU | 说明 |
|---|---|---|---|
| filter gold 规则 `-Oz6`（VCF 输入） | **8.8 s** | 455 s | 单遍 |
| 同上 `-Ob` | 9.9 s | 519 s | 原生编码器；+13% |
| 同上，BCF 输入 `-Oz6` | **12.1 s** | 630 s | 渲染为文本后重扫（F4-5） |
| filter `--mind 0.2 --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 -Oz1` | 6.2 s | 272 s | sketch：pass0 2.80 s + mid 0.08 s + pass1 3.04 s；剔除 25 个样本 |
| filter 无位点规则 `-Oz1`（全部存活） | 5.5 s | 269 s | CPU 以 deflate 为主 |
| stats π/dxy/Fst/θ/D（VCF） | **2.7 s** | 71 s | |
| 同上（BCF） | 3.1 s | 56 s | CPU 少 21%，wall 多 12%：reader 受限（doc_v3 STATUS 已分析） |
| stats `--min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 --pi`（边过滤边统计） | 3.4 s | 170 s | 过滤包装层对每条存活记录做 `edit_vcf_line`（F4-7） |
| stats `--freq --hardy --het --missing --depth --summary` | 5.3 s | 232 s | 两族报告 → 读两遍（F4-6） |
| `reed index`（VCF.gz） | 2.4 s | 77 s | bcftools 142 s |
| `reed index`（BCF） | 2.2 s | 46 s | |

### 3.2 子集（Chr01–03，480,952 条，1.0 GB VCF.gz / 1.18 GB BCF）

| 命令 | wall | user | worker_busy |
|---|---|---|---|
| filter VCF 输入 `-Oz6`（`-M 2 --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05`） | 1.24 s | 59 s | 61.8 s |
| filter BCF 输入 `-Oz6` | 1.69 s | 85 s | 86.4 s（+40%） |
| filter BCF 输入 `-Ob` | 1.86 s | 94 s | 95.2 s |

### 3.3 线程扩展（全量，`--geno 0.1 --maf 0.05 -Oz1`，全部存活）

| 线程 | wall | user | worker_busy ÷ workers | 判定 |
|---|---|---|---|---|
| 4 | 131 s | 263 s | 262 / 2 = 131 s | CPU 受限 |
| 8 | 44.7 s | 263 s | 266 / 6 = 44 s | CPU 受限 |
| 16 | 46.1 s | 264 s | 267 / 14 = 19 s | **污染**（wall 2.4× 预期） |
| 32 | 35.0 s | 266 s | 268 / 30 = 9 s | 污染 |
| 56（同一时段） | 29.3 s | 269 s | 280 / 54 = 5.2 s | 污染 |
| 56（负载 0.3 时重测两次） | **5.5 / 5.6 s** | 269 s | 5.0 s | CPU 受限 |

结论：引擎在 4–8 线程下精确表现为 CPU 受限（wall = worker_busy ÷ workers），没有串行段；16/32/56 线程那一组的偏差来自共享机器的其他负载（当时 1 分钟负载 2.97，13 个登录用户），**不是代码缺陷**。
由此得到的测量纪律写入 04-performance.md §4：记录运行时负载，3 次取中位数，污染样本不入日志。

### 3.4 `--mind` 的两遍分解（`--timing`）

`open 3 ms、pass0 2,804 ms、mid 75 ms、pass1 3,036 ms，worker_busy 145 s，在途峰值 912 MB`。pass0 与 pass1 各约 3 s，mid 并行化后可忽略（doc_v3 O9 的结论得到确认）。

---

## 4. 发现

| # | 级 | 位置 | 问题 | 证据 | 修法 |
|---|---|---|---|---|---|
| F4-1 | **P1（语义）** | `reed-cli/src/translate/filter.rs:116-117`、`reed-filter/src/eval/tail.rs:33-36`、`reed/src/filter/mod.rs:256-268` | `--geno`/`--mind` 是"丢弃条件"，其余阈值是"保留条件"；`TailRules.max_missing` 字段存的是丢弃比例 | §2；用户裁决 | 02-cli-semantics.md |
| F4-2 | P2（文档） | `doc_v3/STATUS.md` "未完成" 表 | 写 "C3 `TagTable` 未做"，但 `reed-core/src/tag_table.rs` 存在、`reed_io::header::tag_table()` 已接入 `plan::compile`、`reed-filter/Cargo.toml` 已无 noodles | `grep noodles crates/reed-filter/src` 为空 | STATUS 改为"已做"，并补 C3 的验证行 |
| F4-3 | P2（CLI） | `reed-cli/src/translate/stats.rs:341` `let _ = verbosity;` | `reed stats --verbosity` 被接受但无任何效果；filter 中只映射为 `quiet` | 代码 | 两个命令统一：0 静默、1 摘要、2 阶段计时（= `--timing` 到 stderr）、3 调试（`RUST_LOG=debug`）；没实现的级别不接受 |
| F4-4 | P2（CLI） | `reed index` 帮助 "or `-` for stdin"；`reed/src/index.rs:457` 用 `std::fs::metadata(input)` 取文件尾虚拟偏移 | `cat x.vcf.gz \| reed index -` 报 `No such file or directory`，而非"stdin 不支持"的明确错误 | 实测 | 索引需要文件尾偏移与可重读，stdin 不可能；帮助文本改为 "file only"，参数校验报清楚的错误 |
| F4-5 | P2（性能） | `reed/src/bcf_input.rs`（`native` 只在无过滤计划时启用）、`reed/src/filter/mod.rs:191` | BCF 输入的 filter 把每条记录渲染成 VCF 文本再走文本扫描器：全量 +175 CPU-s（+38%），子集 +40% worker 时间 | §3.1、§3.2 | 见 04 §2 O-4：为 `RecordPipeline` 提供 BCF 原生输入（GT 表 + 计划内 FORMAT 列的类型化解码，`SiteView` 的 BCF 实现），跳过渲染；输出侧仍渲染→编辑→（编码） |
| F4-6 | P2（性能） | `reed-cli/src/translate/stats.rs:264-270` | `--freq/--hardy/--het/--missing/--depth` 与 `--summary` 同时请求时读两遍（stdin 还要落盘） | 代码；全量 5.3 s vs 单族 3.0/2.8 s | 把 summary 累加器并入 report job 的 partial（都是可加累加器），一次读取 |
| F4-7 | P3（性能） | `reed/src/filtered.rs` `filter_line`：每条存活记录 `edit_vcf_line` 重写 | 边过滤边统计时 CPU 从 71 s 升到 170 s；重写只为把掩码基因型改成 `./.`，而内层 stats 只消费 `GtCode` | 代码；§3.1 | 让 `FilteredJob` 直接把掩码后的 `codes`（和计划列）交给内层（内层接口改为消费 `ScanOutput` 而不是文本），不重写文本；文本重写只保留给需要字节的消费者 |
| F4-8 | P3（可复现） | `reed-filter/src/plan/provenance.rs:58` `##reed_command` 写原样 argv（含二进制绝对路径、`-S` 文件路径、`--threads`） | 同一语义在不同机器/路径下 header 不同；确定性测试只能比正文 | 实测 header | 与 bcftools 一致地保留 argv，但把路径参数按给定形式而非绝对化；或在 `--no-version` 之外提供 `--no-command`；文档注明"正文确定、header 含环境信息" |
| F4-9 | P3（CLI） | `reed stats` 帮助 `--pop` 行 "`-g, --pop`" | 短参数 `-g` 仍是主短名（doc_v3 04 §2.1 要求 `-g/--groups` 隐藏） | `--help` | 去掉 `-g` 短名（保留隐藏长别名 `--groups` 一个版本） |
| F4-10 | P3（CLI） | `reed stats` 报告模式对 `-t/-T` 报 "not supported yet" | 文档已记录；但 `-t` 对窗口统计可用、对报告不可用，用户难以预期 | 实测 | 报告模式实现 `-t/-T`（谓词已在 `FilteredJob` 路径里，只需把 `targets` 作为 head 规则传入报告 job） |
| F4-11 | P3（测试） | `bench/e2e/results.tsv` | V0–V6 的数字只记在 `doc_v3/STATUS.md`，未入日志（doc_v3 03 §3 协议未执行） | 文件 | 补登；每行带负载与中位数（04 §4） |
| F4-12 | P3（代码） | `crates/reed-filter/src/spec/{mod,render}.rs`、`plan/{mod,provenance}.rs` | 类型仍叫 `FilterSpecV2`（doc_v2 05 §4 禁止版本后缀；legacy 已删，后缀失去意义） | grep 15 处 | 改名 `FilterSpec` |
| F4-13 | P3（文档） | `reed filter --help` `<INPUT>` 行 | 仍写 "Input BGZF VCF (.vcf.gz)"，而 BCF/plain/gzip/stdin 都已支持（stats/index 的帮助已更新） | `--help` | 改为 "Input VCF/VCF.gz/BCF, or `-` for stdin" |
| F4-14 | P3（引擎） | `reed-engine/src/frame_chain.rs:120` `STALL_LIMIT = 120 s` | 链等待超时是固定 120 s；在极端慢存储上一个 1 MiB chunk 的解压+走链不可能超过它，但 NFS 卡顿时会误报 | 代码 | 超时改为"自上次任何 chunk 发布以来"的全局静默时间而非单个 cell 的等待时间；并在错误里写明 chunk id 与等待时长 |

doc_v3 的发现 F0–F15 处理状态：F0/F0b/F1–F5/F6/F7/F9/F11/F13（部分）/F15 已修；F8（reed-filter 去 noodles）已修但 STATUS 未更新（F4-2）；F12（基准日志）未做（F4-11）；F14（扫描微优化）做了第 1 级。

---

## 5. 值得保留的做法

- `--timing` 的阶段报告：本轮所有性能判断都基于它，应成为每个基准行的附件。
- 三策略（single/two-pass/sketch）输出逐字节相同的门，与 BCF 输入/输出的跨格式 md5 门。
- 对参照物行为的实测记录（vcftools `--max-missing` 的方向、bcftools `!~` 的缺陷、plink2 的 ε 规则）。
