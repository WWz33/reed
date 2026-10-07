# doc_v3 实施状态（2026-10-08）

> 对照 `05-roadmap.md` 的 V0–V6。每项注明 commit 与验证方法；偏差单独登记（格式沿用 doc_v2/STATUS.md）。
> 全量门：`cargo test --release --workspace` 全绿；`cargo clippy --release --workspace` 无 deny 级问题。

## 已完成

### V0 语义修正（C0、C0b）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| C0 | `--mind` 改为 plink2 语义：统计**通过头部规则**的记录上、**掩码之后**的样本缺失；分母 = 头部存活的记录数。`MissingJob`（SinglePass/TwoPass 的 pass0）改走与单遍相同的 `RecordPipeline`（含掩码字段扫描），sketch pass0 同步改为 head 存活 + post-mask | 见下 | 与 plink2 对拍：fixture 上 `--vcf-min-gq 20 --mind 0.1` 两边都保留 {A,C}（B 掩码 2/10、D 缺失 3/10 被剔除，C 恰为 0.1 保留）；另有头部规则收缩分母的回归用例（s0 旧口径 91/100→剔除、新口径 1/10→保留），三策略输出一致 |
| C0b | `cohort: Option<CohortRule>` 拆为 `sample_select: Option<SampleSelect>`（`-s/-S/--keep/--remove`，keep ∩ ¬remove）与 `mind: Option<f64>`；先静态子集、后 `--mind`（plink2 同序）。`-s ^A,B` 排除语法；同极性参数（`-s`+`-S`、`-s`+`--keep`、`-s ^`+`--remove`、重复 `-s`）报错而非静默择一 | 见下 | 组合用例：`-S keep.txt --mind 0.2` 保留集合正确；`to_args()` 往返测试覆盖组合与冲突报错 |
| 防御 | `eval_head`/`eval_tail` 的计数器长度与规则数不一致时断言失败（此前 zip 静默跳过规则，导致 pass0 漏评估头部规则） | 见下 | 该断言在开发中直接定位了 pass0 funnel 未初始化的缺陷 |

开发中发现并修复的缺陷（非文档列出）：pass0 的 `FunnelPartial` 未按规则数初始化，`eval_head` 的 `zip` 空转使头部规则被跳过（SinglePass/TwoPass 的 `--mind` 因此按旧口径统计）；pass0 的扫描器未带计划 FORMAT 列，掩码未参与缺失计数。

### V1 串行瓶颈与资源

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| O5 | stdin 单遍任务直接流式（探测前缀用 `Read::chain` 接回）；多遍任务（`--mind`）落盘到 `--tmp-dir`，`StdinSpill` 的 Drop 在正常/错误/提前返回时删除；`reed stats`/`qc`/reports 的 CLI 改为一次 open 并把 source 交给 driver（`run_*_opened`）；报告模式两族同用时 stdin 落盘一次 | 见下 | 新 CLI 测试 3 条：单遍 stdin 与文件输出一致且无临时文件；多遍 stdin 与文件输出正文一致且 spill 删除；错误路径同样删除。`/tmp` 历史残留 5 个已清空 |
| O4 | `reed index` 并行：`ChunkBytes` 携带 `compressed_start` + 每块 `(压缩长, ISIZE)`；worker 从 chunk 文本偏移换算虚拟偏移，reducer 按序配对相邻记录起点（末条到文件尾）并调用 `Indexer::add`；`eval` 无；`-c/-t/-o/-f/--threads` 对齐 bcftools | 见下 | 1 GB 子集 480,952 条：56 线程 0.47 s、单线程 10.4 s，1/3/56 线程 CSI sha256 相同；全量 `clean.vcf.gz`（7.0 GB，3,315,604 条）2.7 s（bcftools index 142 s）；bcftools 用 reed 的索引查询与流式扫描逐字节一致（子集与全量各 2 个区域） |

| O3 | **根因不是缓冲池**：`-O v|z` 忽略 `-o`，`-Ov -o /dev/null` 实际把 35 GB 写到 stdout，基准测的是管道消费者（复现 13.3 s，reducer 的 `accept` 占 13.1 s）。修法：`-O` 只选格式、`-o` 命名文件；扩展名推断仅在无 `-O` 时生效。缓冲池与 mimalloc 都做了 A/B（`-Ov` 全量：系统分配器 2.9/3.0 s vs mimalloc 3.2/3.2 s；缓冲复用与 `mem::take` 均 2.9 s），按 03 §1 的"无收益则撤回"全部撤回，不引入池 | 见下 | `-Ov` 全量 56 线程 2.9–3.2 s（目标 ≤6 s）；`-Oz1` 5.4 s；gold 等价命令 9.0/9.3 s、1,768,759 条，`gold/scripts/check_parity.py` 对 `bcftools.full` **PASS**（位点 jaccard 1.0，基因型不一致 0）。新增 CLI 回归测试固定 `-O`/`-o` 组合与扩展名推断 |

| O1+O2 | 报告改为 worker 内计算与渲染：`ReportPartial` 携带每个报告文件的字节缓冲（afreq/hardy/vmiss/ldepth），HWE 与 `%g`（`g6_into`，指数形式才用临时 String）都在 worker 完成，reducer 只按序 append。按请求装配：`SiteRow` 的四次 `to_vec()` 删除；计数改走 `SiteGroupCounts::add_code`（不再构平面），DP 只在 `--depth` 时解码，平面只在 `--het` 时构建 | 见下 | 全量单报告：`--hardy` 13.3 → **2.7 s**（目标 ≤3.5），`--freq` 2.86→2.61、`--missing` 2.82→2.58、`--het` 2.80→2.74、`--depth` 3.71→2.64；组合 `--freq --hardy --het --missing --depth` **3.0 s**（目标 ≤4）；7 个输出文件与改动前 **md5 全部相同**（全量 3,315,604 条） |

### 待办

- V2 B0、V3 B1、V4 B3a、V5 B2、V6 B3b；清理 C1–C4、C7 与 CLI 一致性（04 §2）。V1（O1–O5）完成。

### 与文档的偏差

1. **F3/O3 的诊断与实测不符**：01-review 的"F3 未压缩输出 13.2 s"与 O3 的缓冲池方案建立在错误的测量上——`-O v` 当时忽略 `-o`，全量 35 GB 文本进了基准脚本的管道。修正路由后未压缩路径 2.9–3.2 s，无需缓冲池。O3 的验收目标（≤6 s）以修正后的同口径达成。
