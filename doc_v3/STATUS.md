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
| C0c（补，非文档列出） | 默认双等位门与 R7 之前的 AC 一致性：`plan::compile` 增加 `biallelic_gate` 参数（filter=true、stats=false，stats 的多等位处理按各报告器定义）；spec 不允许多等位输入（`multiallelic != "allow"` 且 `max_alt_alleles <= 1`）时压入 `HeadRule::BiallelicGate`（funnel id `site.biallelic_gate`），语义与 §11.6"允许多等位输入（`multiallelic != reject` 或 `max_alt_alleles > 1`）"及 doc_v2 STATUS"位点级双等位门"一致；spec 允许多等位而 AC/MAF/HWE 需要重算时（`-s/--mind` 会自动打开 AC/AN 更新）在写输出前报错——A7 的"在 plan 期拒绝"分支（逐 ALT 计数见 R7） | 缺陷复现：`A>G,T`（GT `1/2`）+ `-s s1`，改前写 `AC=2` 而 bcftools 写 `AC=1,1`（AC 头声明 Number=A，单值非法）；改后默认门先剔除该记录，显式允许时给出带 `-M 2`/`-I` 指引的错误 | 新增 e2e `default_biallelic_gate_drops_multiallelic_records`、`ac_update_with_allowed_multiallelics_is_rejected`；全量 `clean.bcf` 扫描 0 条多等位（`bcftools view -m3`，8.1 GB 全文件）；区域 A/B：隐式门 ≡ 显式 `-M 2`（正文 md5 相同），条数 == bcftools（8364 条）。注意：spec 用 `max_alt_alleles == 0` 表示"未设"，因此 CLI `-M 1`（只留无 ALT 记录）本就无法表达（旧行为=不约束，新行为=≤1 ALT），属既有编码限制 |
| rustc 警告清零 | `unreachable_pub`（bin crate 的 `pub` → `pub(crate)`，cargo fix 机械修改）、未用 import/`mut`/变量、死代码（`Shared.first`、`QcScratch.n_samples`、`IndexScratch.starts`、`ChunkDecision.n_records`、`ProjectedSfsAcc.group`、`FromI64`、`parse_gt_pair`、测试 helper）、重复 `#[test]`、`items after a test module`、bench 的 `missing_docs` | `cargo build --workspace --all-targets` **0 warning** |

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
| CLI stats 复用过滤参数（区域/样本） | `reed filter` 与 `reed stats` 共用抽出的 clap 参数组（`cli/filter_flags.rs`：`-r/-R/-t/-T` 与 `-s/-S/--keep/--remove`）；stats 的 `-r/-R` 走索引 I/O 且同样的区域同时作为记录谓词（索引区间是整块，必须再判），`-t/-T` 走 `region_allows`（`^` 排除），`-s/-S/--keep/--remove` 通过把被排除样本的 group 置空实现（分组分母随之收缩）；报告模式暂不支持 `-t/-T`（明确报错） | `-r` 与 `-t` 在子集上逐行一致；`-s s0` 与"先用 bcftools 物理子集再跑 stats"的 popA 行**逐字节相同**；未知样本报错；CLI 测试 `stats_accepts_filter_region_and_sample_flags` |
| CLI `-q/-Q`（min-af/max-af） | `-q/--min-af`、`-Q/--max-af` 落地：header 声明了 INFO/AC+AN 时在位点阶段用 INFO 值判定（不要求解码基因型，与 bcftools 的 `view` 阶段一致），缺 INFO 的记录回退到尾部按基因型计数；阈值按 bcftools 的 `sscanf("%f")` 量化为 f32（否则 AC=34/AN=680 对 0.05 这类恰好相等的情形会与 bcftools 反向） | 与 `bcftools view -q/-Q` 在 50,000 条上逐记录集合一致（7 组阈值，含 0.05、0.1234567 等边界）；fixture 三组断言入 `e2e_bcf.rs` |
| CLI `--pbs A,B,C` | `--pbs` 改为接受群体名三元组（A 为焦点，可重复）；列名按请求的三元组命名（此前行数据用 0/1/2 索引、列名硬编码 pop1/pop2/pop3，是错的）；`parse_pbs_triples` + 单元测试 | `--pbs pop3,pop2,pop1` 的 `pbs` 值 = `--pbs pop1,pop2,pop3` 的相反数、fst 列随 A/B/C 重排；错误用例（两个名字、未知群体）报错 |
| BCF 确定性 | BCF 输出与线程数无关（1 vs 4 线程字节相同） | `e2e_bcf.rs::bcf_output_is_thread_count_independent` |
| CLI §2 | `-o/--output`（filter/stats/qc/index）；stats 主名 `--pop`（`--populations` 可见别名，`-g/--groups` 隐藏）；`--contig-lengths`（旧名别名）；filter 的 plink2 别名 `--vcf-min-dp/--vcf-max-dp/--vcf-min-gq`；`--compress-level` 隐藏；`--memory`、`--verbosity`、`--no-version` 实现；输入/`--threads` 帮助文本更新 | 手动验证：`--no-version` 后无 `##reed_*` 行；`--pop`/`-g` 两种拼写输出一致；`--help` 文本 |

### 未完成（明确记录）

| 项 | 状态 | 原因 |
|---|---|---|
| C3 `TagTable`（reed-filter 去 noodles） | **已做**（doc_v4 01 F4-2 更正：本行此前记为未做，与代码不符） | `reed-core/src/tag_table.rs` 存在；`reed_io::header::tag_table()` 已接入 `reed_filter::plan::compile` 的入口（`crates/reed/src/filter/mod.rs`、`crates/reed-cli/src/translate/stats.rs` 均从 header 构建）；`crates/reed-filter/Cargo.toml` 无 noodles（`grep noodles` 为空）。守卫：`bash scripts/check-doc-drift.sh` |
| C4 合并 QC 入口 / reed-qc 并入 reed-stats | 已做 | `reed qc` 保留为隐藏弃用别名（stderr 警告）；`reed-qc` crate 删除，代码并入 `reed-stats::summary`（`summary/{mod.rs,histogram.rs,report.rs}`，`serde` 依赖随之移入），其测试改引用新路径；`reed/src/{stats,reports,pixy,qc,qc_report}.rs` 按 doc_v2 03 §8 移到 `reed/src/stats/{job.rs,writers/*.rs}`，旧公共路径（`reed::reports`/`reed::pixy`/`reed::qc_report`）由 crate 根 re-export 保持 | 342 测试全绿、clippy 干净、工作区从 9 个 crate 减到 8；`writers/*` 内部尚未按格式细分（`plink2_table`/`vcftools_depth`/`bcftools_stats`/`log`），文件保持现有单文件结构 |
| C2 的 `pub(crate)` 化 | 改为 `#[doc(hidden)]` | 对拍套件是 integration test，需要跨 crate 访问这些累加器；移入单元测试是纯机械工作 |
| C5 类型化 `EngineError` | 已做（doc_v2 03 §4.4 的形状，doc_v3 C5 引用该节）：`Io(std::io::Error)`、`Input(reed_io::IoError)`（保留 Locus）、`Job(Box<dyn Error + Send + Sync>)`（保留类型与 source 链）、`WorkerPanicked`/`MemBudgetExceeded`/`Aborted`；`EngineError::job(msg)` 造消息型错误，`From<CoreError>`/`From<IoError>`/`From<io::Error>` 三处转换 | 构造点全部改走 typed 转换；`source_errors_propagate` 断言传输错误落到 `Input(IoError::Io(_))`；335 测试全绿 |
| C7 clippy pedantic | 已做 | pedantic 警告 ~700 → **44**：`expect()` 逐条审查（生产代码全部带 `invariant:` 前缀，输入相关条件返回错误；测试侧由 clippy.toml 的 allow-*-in-tests 覆盖）→ 豁免该 lint 并在 Cargo.toml 注明审查结论；`cast_possible_truncation` 逐条处理（修复 BCF 记录长度前缀与 GT 等位下标的真实截断风险，其余 22 个文件按"构造上有界"加文件级 `#![allow]` + 不变量说明）；`cast_lossless/possible_wrap/sign_loss/precision_loss` 四类按 04 §C7"风格类按项豁免"处理；保留 `format_push_string`/`naive_bytecount`/`float_cmp` 等 44 条有信息量的提示；测试/基准目标另修 3 处 deny 级（`items after a test module`、多余 `.into()`、无效果的 struct update） | `cargo clippy --workspace` 0 error、44 warning（全部为性能/正确性提示）；`--all-targets` 0 error、164 warning（含测试/基准） |
| stats 的位点/阈值过滤参数 | 已做 | CLI 层共用 `SiteFilterArgs`（`-v/-V/--snps-only/-m/-M/-f/-i/-e/--min-dp/--max-dp/--min-gq/--set-gts/-q/-Q/--geno/--maf/--max-maf/--mac/--max-mac/--hwe`，INFO 更新类仍只属 filter）；执行用 `reed::filtered::FilteredJob` 包住 stats job——单遍内跑 filter 的 `RecordPipeline`（head→mask→tail）并用 `edit_vcf_line` 把掩码基因型改写成 `./.`，内层消费者看到的就是"物理 filter 过"的字节 | 五种参数组合下 `stats --X` 与"`filter --X` 落盘 + `stats`"的 π 表**逐字节相同**（含 BCF 输入；CLI 测试 `stats_site_filters_match_the_filter_then_stats_composition`） |
| stats 的 `--mind` | 未做 | 需要先跑一遍全局缺失还原（filter 的 pass0），属多遍策略；未暴露该参数以免静默忽略 |
| B 门补充 | 已补：**noodles 差分**（`bcf_vs_noodles.rs`：fixture 全部记录 + `clean.bcf` 前 10 万条，逐 token 值等价；noodles 的浮点文本与 htslib 不同，已注明）、**链切分穷举**（`e2e_bcf.rs::chain_framing_is_boundary_independent`：同一 BCF 按 13–512 字节重切块，输出记录不变）、**BCF 输出确定性**（1/2/3/7/56 线程字节相同） | 未做：typed-value proptest（该设计里 typed 值不单独成层） |
| 多等位逐等位重算（R7 的一部分） | 未做（有守卫，不再静默写错） | doc_v3 05-roadmap R7："V1 之后任意时间"；现状：默认门剔除多等位记录；显式允许（`-M 3+` 或 `multiallelic = "allow"`）且需 AC/MAF/HWE 重算时在写输出前报错，不会写出 Number=A 单值 |
| V6：B3b 原生 BCF 手术 | 不做（文档标注"可选/性能阶段"） | 渲染+编码路径已正确；原生手术（原样拷贝/原位掩码/INFO 重写）需要在 BCF 字节上重写一遍过滤循环——正是渲染路线刻意避免的"两份消费者"；现状代价已实测：BCF→BCF 17.0 s（VCF.gz→VCF.gz 11.1 s），都在可接受范围 |
| O6/O7 扫描与编辑微优化 | 未做 | 需微基准门（`benches/genotype_scan.rs`、`benches/line_edit.rs`） |
| 基准日志（03 §3 协议） | 部分：本轮的实测数字都记在本文档；未写入 `bench/e2e/results.tsv` | — |
| O7 按连续段拷贝 | 已做（+1%） | `benches/line_edit.rs`：3.06 µs/记录（348 样本，~14% 掩码）＝ 8.8 ns/样本；撤下改动后 3.09 µs（criterion `change +1.05%`，p<0.05）——收益真实但在噪声边缘 |
| O10 `--timing FILE` | 已做 | 输出各阶段 wall（open/pass0/mid/pass1/total）、records、chunks、worker_busy 与在途峰值；全量 `--mind 0.2` 的分解：pass0 2.7 s、mid 0.94 s（并行后 0.03–0.05 s）、pass1 9.8 s |
| O9 sketch 中段并行 | 已做（第 1 项） | `mid_stage` 按 chunk `par_iter` 计算决策、再按 chunk 序折叠计数（整数和与序无关，确定性不变）；全量 `--mind 0.2` 的 mid **938 → 32–47 ms**，wall 13.7 → 12.7–12.9 s |
| O9 投机单遍（第 2 项） | 不做（已测量后决定） | 只有"最终没有剔除任何样本"才省下 pass0（2.7 s）；一旦剔除了样本，投机输出必须丢弃并**重读输入**（35 GB 无法缓冲），代价是整整一遍 pass1。文档要求按概率权衡，本数据 `--mind` 的用户通常正是要剔除样本，因此不默认开启 |
| O6 单遍样本列扫描（第 1 级） | 已做 | `benches/genotype_scan.rs`（criterion A/B）：`split_needed` 由"每个子字段一次 memchr"改为单遍字节循环，GT-only 4.05 → **3.84 µs**（−4.7%）、GT+DP+GQ 9.56 → **8.65 µs**（−9.5%），均 p<0.05；第 2 级（AVX2 位图）未做 |

## 最终验证（2026-10-08，本轮所有改动之后）

| 检查 | 结果 |
|---|---|
| `cargo test --release --workspace` | **344 通过 / 0 失败**（含 5 个 v02-oracle 回归、BCF 输入/输出、三条链门、双等位门与 R7 守卫用例） |
| `cargo build --release --workspace --all-targets` | **0 warning**（此前 25 条 `unreachable_pub` + 未用 import/死代码等已清零） |
| `cargo clippy --release --workspace` | 0 error、47 warning（全部为性能/正确性提示）；`--all-targets` 0 error、170 warning |
| 双等位门 | 全量 `clean.bcf` 0 条多等位（`bcftools view -m3`）；`-r Chr01:1-2000000` 隐式门 ≡ 显式 `-M 2`（正文 md5 `18676382…` 相同），条数 8364 == bcftools |
| gold 对拍（双等位门之后重跑） | `gold/scripts/run_reed.sh` 全量 9 s（56 线程），正文 md5 `7aedc59ac798…` 与 V0/V1 记录相同；`check_parity.py` 对 `bcftools.full`：**PASS**（1,768,759 位点，gold-only 0 / cand-only 0，jaccard 1.0，基因型不一致 0）——数据本身 0 条多等位且 gold 命令显式 `-M 2`，故输出不动 |
| 跨格式一致 | 全量 `clean.bcf` ≡ `clean.vcf.gz`：stats 五个输出文件 md5 两两相同（π/dxy/θ/Tajima/Fst，9501 窗） |
| 渲染一致 | `clean.bcf` 全量 3,315,604 条渲染 sha256 == `bcftools view -H`（`548de25b…`） |
| BCF 输出一致 | 全量 `-Ob` 输出经 bcftools 读回的正文 sha256 == `-Oz` 输出（`b0b22079…`） |
| 微基准 | `benches/line_edit.rs`、`benches/genotype_scan.rs`（O6/O7 的 criterion 门） |

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
| `reed stats` on BCF（π/dxy/θ/D/Fst） | 5.1 s（渲染路径） | **2.30–2.48 s**（同数据 VCF.gz 2.17 s；输出与渲染路径逐字节相同） | V3/B1 原生记录路径；**门仍未达**（BCF 慢 ~6%），根因见下节 |
| stdin 单遍 | 整体落盘 | 直接流式，与文件一致；无临时文件残留 | O5 |

### V3 性能门复盘：`stats --pi` 在 BCF 上仍慢 ~6%（2026-10-08 第二轮）

02-bcf.md §4 的门是 `clean.bcf` 上 `reed stats --pi` **≤** VCF.gz 上的同一命令（预期更快："只解 GT 字节"）。
现状：`--pi` BCF **2.30–2.48 s** vs VCF.gz **2.17 s**（56 线程，热缓存，各 3 次）；五统计（π/dxy/θ/D/Fst）BCF **2.45 s** vs VCF **2.44 s**（持平）。输出 md5 `825cd749…` 两格式相同、且等于渲染路径旧输出。门未达，但差距从 +70% 收到 ~+6%（五统计已持平）。

本轮已做（全部保留）：

| 改动 | 内容 | 效果 |
|---|---|---|
| BCF 原生记录路径（02-bcf.md §2.3–§2.5） | 新增 `bcf::gt`：`site_header`（只读记录头，shared 块其余部分与未计划的 FORMAT 字段按描述符跳过）+ `GtTable`（64 KiB 字节对查表，条目由渲染器的 `fmt_gt` + 文本扫描器的 `decode_gt` 生成，与文本路径逐位一致，65,536 项单测锁定）+ `decode_record_gt`；`InputJob` 增加 `BcfRecordFn` 原生入口，stats 在**无过滤计划**时直接消费定界后的记录字节（不再渲染整行文本）；`FilteredJob` 在 plan=None 时不再包一层 | BCF stats 的 worker CPU 72.5 s → **57.6 s**（user）；wall 3.08 → 2.3 s |
| `FrameChain` 定向唤醒 | `publish` 原为每次 `notify_all`；在 112 个在途 chunk 的窗口下每个 chunk 唤醒全部等待者（惊群+单锁竞争）。改为按 cell 的等待位图，仅当该 cell 有等待者时通知 | BCF wall **−0.3 s**（2.8 → 2.5 s 档） |
| reader 直读 | `next_chunk` 预 reserve 整个 chunk 缓冲；块头解析后把 payload **直接读进 chunk 缓冲**（此前每块先读进临时 `Vec` 再 `extend_from_slice`：每字节两次拷贝 + 每块一次分配） | reader cut（探针）：BCF **3.0 → 2.34 s**、VCF **2.55 → 2.15 s** |

实测的瓶颈定位（根因）：

1. **墙是单线程 reader（BAPP 切块），不是 worker**。`--timing` 探针显示 reader 的 `next_chunk` 累计 2.2–2.3 s（BCF）/2.1 s（VCF），与总 wall 2.3/2.17 s 基本相等；worker busy 56.3 s（BCF）/72.6 s（VCF）÷56 线程 = 1.0/1.3 s，低于 reader。
2. **单位字节的 reader 成本相同**（~3.4–3.7 GB/s），而 BCF 文件大 15%（8.06 vs 7.0 GB），所以 BCF 必然慢 ~0.15–0.2 s；worker 侧 BCF 省下的 CPU（−21%）在 56 线程下不构成墙。
3. reader 并非设备上限：`dd bs=1M` 11.2 GB/s、8 路并行合计 43 GB/s；64 KiB 粒度顺序读的实测上限约 5.2 GB/s（python 对照），即每字节的 syscall+拷贝开销是主项。

试过并回退（留档，避免下次重复）：

| 尝试 | 结果 | 处置 |
|---|---|---|
| 池化 chunk 缓冲（`ChunkBuf` + 全局池，reader 直读进池化页、尾部带走） | 逻辑正确（346 测试全绿、输出逐字节相同），但 BCF 10.2 s / VCF 9.6 s，**比基线慢 4 倍**；原因未及定位（候选：每 chunk 新缓冲 `resize` 清零、尾部搬移、读粒度），停止实验 | 回退 |
| `TICKETS_PER_WORKER` 2 → 8 | 两种格式都略慢（BCF 3.29 / VCF 2.91），在途内存压力 | 回退（保持 2） |

要达门还需要（下一轮的方向）：把 reader 的每字节成本压到 worker 之下——把 64 KiB 块读合并成 MB 级读且不引入第二次拷贝（池化缓冲做对版：缓冲以 `len == capacity` 停放、避免清零；或 reader 拆成两线程；或用 `posix_fadvise` 预取）。只要 reader 降到 ~1.2 s，BCF（worker 1.0 s + 文件更小优势）就会快于 VCF。

### 与文档的偏差

1. **F3/O3 的诊断与实测不符**：01-review 的"F3 未压缩输出 13.2 s"与 O3 的缓冲池方案建立在错误的测量上——`-O v` 当时忽略 `-o`，全量 35 GB 文本进了基准脚本的管道。修正路由后未压缩路径 2.9–3.2 s，无需缓冲池。O3 的验收目标（≤6 s）以修正后的同口径达成。
