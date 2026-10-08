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

### V2–V4 BCF 输入（B0+B1+B3a）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 解码/渲染 | `reed_io::bcf`：`BcfDict`（header `IDX=`，下标不重排）、`BcfRecordRef`/`Cursor`（typed 描述符、int/float 哨兵、扩展计数）、`render_record`（htslib `vcf_format1` 规则：GT 按各等位定相位、`kputd` 六位有效数字、缺失 `.`、向量结束截断、flag 只出键） | 见下 | **与 bcftools 逐字节一致**：`clean.bcf` 全量 3,315,604 条渲染正文 sha256 = `548de25b…`（与 `bcftools view -H` 相同）；边界 fixture（单倍体、多等位、flag、字符串、半缺失）逐字节对齐（`reed-io/tests/bcf_render.rs`） |
| 修复 format_g6 | 既有 `%g` 渲染在"指数 ≥ 5 的整数"（如 QUAL=128440）会误删小数点前的尾零；先按 6 位有效数字取整再选风格（C %g 语义） | 见下 | `128440→"128440"`、`999999.5→"1e+06"` 等用例；gold VCF 的 80,228 个渲染浮点循环测试继续通过 |
| 长度前缀链 | `RecordFraming::{Newline,Raw}` 进引擎；`FrameChain` 沿 chunk 序传递"跨界记录"（字节 + 已有字节数 + 虚拟偏移），worker **先定界发布、后渲染**；`PublishGuard` 在 worker panic/早退时毒化链（abort token 可能被阻塞的 select! 消费，不能只靠 abort）；120 s 看门狗兜底 | 见下 | 全量 stats BCF 输入与 VCF.gz 输入**逐字节一致**（π/dxy/θ/Tajima/Fst 五文件 md5 相同）；1/56 线程一致 |
| 接入 | `InputJob<J>` 包装器：stats/reports/qc/filter（含 `--mind` pass0）全部走同一文本消费实现；BCF 的 filter 强制 TwoPassReparse；`-r/-R` 区域路径也识别 BCF（探测格式 + 字典 + 跨界截断合法）；`reed index` BCF：链携带跨界记录起点 → worker 算虚拟偏移 → sink 配对相邻记录 | 见下 | fixture e2e（stats/reports/filter 三组，输出 md5 相同）；全量 `reed index clean.bcf` 10.9 s（3.3M 条），bcftools 用该索引查询与流式 `-t` 一致（2 区域 MATCH） |

### V5 BCF 输出（B2 + B4）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| B2（noodles） | `BcfEncoder`（noodles `write_variant_record`，lazy parse + 全缺失字符串的 RecordBuf 回退）：每条记录 → BCF 字节 → 现有 `BlockWriter`（切块规则不变） | 见下 | fixture 与 50k 条：bcftools 读回与 VCF 输出逐字节一致 |
| B4（原生） | `NativeEncoder`：一次文本遍历直接写 BCF（镜像 `render.rs`），整数宽度按 `bcf_enc_vint`（−120…127 / −32760…32767），GT 按 htslib 前缀定相规则（首等位在后继全定相时置位；缺失等位写 0），字典下标来自输出的 `StringMaps`；`REED_BCF_ENCODER=noodles` 可切回 B2 做 A/B | 见下 | 与 noodles 编码器**逐字节一致**（`bcf_encode.rs`，fixture 全部记录）；往返渲染 == 原始记录文本 |
| header IDX | `bcf_header_text`：BCF 输出的 header 文本给每个 FILTER/INFO/FORMAT 行补显式 `IDX=`（htslib 约定），把 htslib 的字典钉到 `StringMaps` 的编号上；`-Ob` + `-W=tbi` 报错（BCF 只有 CSI） | 见下 | 字典删除/新增用例：`-x INFO/QD,INFO/FS --fill-tags AC,AN -Ob` 后 bcftools 读回与 `-Oz` 输出逐记录一致；`bcftools index -s`/`-r` 用 reed 写的 CSI 查询与流式一致 |
| 性能 | `-Ob` 全量（3,315,604 条，8.1 GB 输出）**13.3 s** vs `-Oz` 11.1 s = 额外 **20%**（门 ≤30%）；noodles 路径为 41.8 s（+277%，未达门，故按 B4 的触发条件落地原生编码器） | 见下 | `bcftools view -H` 逐条与 VCF 输出 sha256 相同（`b0b22079…`） |

### 清理与 CLI 一致性（04 §1/§2，部分）

| 项 | 内容 | 验证 |
|---|---|---|
| C1 | 删除清单 v2 残留：`reed_ascertainment` 行、`reed-core/src/manifest.rs`、`reed-filter/src/plan/ascertainment.rs`、`Plan::ascertainment`/`with_cohort_subset` 及测试断言（改为断言**不出现**） | `cargo test --workspace` 全绿；全量 filter 输出与 gold 对拍不变（`-S` 路径不受影响） |
| C2（部分） | 删除已撤下的 `NeiMultiAcc`、`HudsonFiniteAcc`（含其测试）；未接入 CLI 的累加器（`FstMatrixAcc`/`WcMultiAcc`/`JointSfsAcc`/`BlockJackknife`/LD/Patterson/projected-SFS）改为 `#[doc(hidden)]`，不再出现在文档化 API 中 | `cargo test --workspace` 全绿（测试数 336 → 331，删除的 5 条是被撤下累加器的测试） |
| C6 | 进程结束无临时文件：O5 的 3 条 CLI 测试覆盖（单遍不产生、多遍删除、错误路径删除） | `cargo test -p reed-cli` |
| §2.3 | 表达式读 `AC/AN/AF/MAC/MAF` 且 header 声明了同名 INFO 时打印"陈旧 INFO"警告（提示 `--maf/--mac` 或 `-x INFO/...`）；`stale_info_hits` 纯函数 + 单元测试 | `cargo test -p reed-filter`；`reed filter x.vcf.gz -e 'AC<5'` 打印警告 |
| C4（可见部分） | `reed qc` 变为 `stats --summary` 的**隐藏**弃用别名（stderr 打印弃用警告，主 help 不再列出）；crate 合并且未做 | `reed qc x.vcf.gz` 打印警告；`reed --help` 无 qc |
| CLI `--pbs A,B,C` | `--pbs` 改为接受群体名三元组（A 为焦点，可重复）；列名按请求的三元组命名（此前行数据用 0/1/2 索引、列名硬编码 pop1/pop2/pop3，是错的）；`parse_pbs_triples` + 单元测试 | `--pbs pop3,pop2,pop1` 的 `pbs` 值 = `--pbs pop1,pop2,pop3` 的相反数、fst 列随 A/B/C 重排；错误用例（两个名字、未知群体）报错 |
| BCF 确定性 | BCF 输出与线程数无关（1 vs 4 线程字节相同） | `e2e_bcf.rs::bcf_output_is_thread_count_independent` |
| CLI §2 | `-o/--output`（filter/stats/qc/index）；stats 主名 `--pop`（`--populations` 可见别名，`-g/--groups` 隐藏）；`--contig-lengths`（旧名别名）；filter 的 plink2 别名 `--vcf-min-dp/--vcf-max-dp/--vcf-min-gq`；`--compress-level` 隐藏；`--memory`、`--verbosity`、`--no-version` 实现；输入/`--threads` 帮助文本更新 | 手动验证：`--no-version` 后无 `##reed_*` 行；`--pop`/`-g` 两种拼写输出一致；`--help` 文本 |

### 未完成（明确记录）

| 项 | 状态 | 原因 |
|---|---|---|
| C3 `TagTable`（reed-filter 去 noodles） | 未做 | 架构纯净性改造，不影响功能；plan::compile 仍直接接收 noodles 头 |
| C4 合并 QC 入口 / reed-qc 并入 reed-stats | 未做 | 结构重构（crate 合并 + 文件迁移），与功能无关 |
| C2 的 `pub(crate)` 化 | 改为 `#[doc(hidden)]` | 对拍套件是 integration test，需要跨 crate 访问这些累加器；移入单元测试是纯机械工作 |
| C5 类型化 `EngineError` | 已做（doc_v2 03 §4.4 的形状，doc_v3 C5 引用该节）：`Io(std::io::Error)`、`Input(reed_io::IoError)`（保留 Locus）、`Job(Box<dyn Error + Send + Sync>)`（保留类型与 source 链）、`WorkerPanicked`/`MemBudgetExceeded`/`Aborted`；`EngineError::job(msg)` 造消息型错误，`From<CoreError>`/`From<IoError>`/`From<io::Error>` 三处转换 | 构造点全部改走 typed 转换；`source_errors_propagate` 断言传输错误落到 `Input(IoError::Io(_))`；335 测试全绿 |
| C7 clippy pedantic | 部分 | `clippy::all` 无错误；生产代码的 `expect()` 已全部带 `invariant:` 前缀（剩余 18 处缺口都在 `#[cfg(test)]`/examples）；`cast_possible_truncation` 约 118 处尚未逐条处理 |
| CLI §2 的 stats 过滤参数复用、`-q/-Q` | 未做 | 需在 stats 驱动接入过滤计划（head/mask/tail）与新的 AF 尾规则（bcftools 的 any/all 语义需要单独对拍） |
| B 门补充 | 已补：**noodles 差分**（`bcf_vs_noodles.rs`：fixture 全部记录 + `clean.bcf` 前 10 万条，逐 token 值等价；noodles 的浮点文本与 htslib 不同，已注明）、**链切分穷举**（`e2e_bcf.rs::chain_framing_is_boundary_independent`：同一 BCF 按 13–512 字节重切块，输出记录不变）、**BCF 输出确定性**（1/2/3/7/56 线程字节相同） | 未做：typed-value proptest（该设计里 typed 值不单独成层） |
| V6：B3b 原生 BCF 手术 | 未做（文档标注"可选/性能阶段"） | 渲染+编码路径已正确；原生手术是加速项 |
| O6/O7 扫描与编辑微优化 | 未做 | 需微基准门（`benches/genotype_scan.rs`、`benches/line_edit.rs`） |
| 基准日志（03 §3 协议） | 部分：本轮的实测数字都记在本文档；未写入 `bench/e2e/results.tsv` | — |
| O7 按连续段拷贝 | 已做（+1%） | `benches/line_edit.rs`：3.06 µs/记录（348 样本，~14% 掩码）＝ 8.8 ns/样本；撤下改动后 3.09 µs（criterion `change +1.05%`，p<0.05）——收益真实但在噪声边缘 |
| O6 单遍样本列扫描（第 1 级） | 已做 | `benches/genotype_scan.rs`（criterion A/B）：`split_needed` 由"每个子字段一次 memchr"改为单遍字节循环，GT-only 4.05 → **3.84 µs**（−4.7%）、GT+DP+GQ 9.56 → **8.65 µs**（−9.5%），均 p<0.05；第 2 级（AVX2 位图）未做 |

## 基准（2026-10-08，本机 56 线程，热缓存，全量 clean.vcf.gz 3,315,604 条）

| 命令 | 改前 | 改后 | 依据 |
|---|---|---|---|
| `filter` gold 等价 `-Oz6` | 9.7 s / 9.0–9.4 s | 9.0–9.4 s（无回退；对拍 gold PASS，body md5 `7aedc59a…`） | V0/V1 |
| `filter` 全部存活 `-Ov` | 15.0 s（实为写入 stdout 管道，测量口径错误） | **2.9–3.2 s**（`-o /dev/null`） | O3 根因 |
| `filter -Oz1` | 5.5 s | 5.4 s | — |
| `filter -Ob`（原生编码器） | 41.8 s（noodles） | **13.3 s**（+20% vs `-Oz`，门 ≤30%） | B4 |
| `stats --hardy` | 13.3 s | **2.7 s** | O1 |
| `stats --freq` / `--missing` / `--het` / `--depth` | 2.86 / 2.82 / 2.80 / 3.71 s | 2.61 / 2.58 / 2.74 / 2.64 s | O1+O2 |
| `stats --freq --hardy --het --missing --depth` | ~15 s | **3.0 s** | O1+O2 |
| `reed index`（7.0 GB VCF.gz） | 48.7 s | **2.7 s** | O4 |
| `reed index`（8.1 GB BCF） | — | 10.9 s | B1 |
| `reed stats` on BCF（π/dxy/θ/D/Fst） | — | **5.1 s**（VCF.gz 3.0 s；逐字节一致） | B1（渲染路径，未达"不慢于"门） |
| stdin 单遍 | 整体落盘 | 直接流式，与文件一致；无临时文件残留 | O5 |

### 与文档的偏差

1. **F3/O3 的诊断与实测不符**：01-review 的"F3 未压缩输出 13.2 s"与 O3 的缓冲池方案建立在错误的测量上——`-O v` 当时忽略 `-o`，全量 35 GB 文本进了基准脚本的管道。修正路由后未压缩路径 2.9–3.2 s，无需缓冲池。O3 的验收目标（≤6 s）以修正后的同口径达成。
