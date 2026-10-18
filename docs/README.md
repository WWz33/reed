# 文档索引（当前规范 / 历史）

> docs/rounds/v4/03-architecture.md §2 A4-7 定的规则：每轮文档并存且互相引用，
> 这一页说明**哪一份是当前规范、哪些是历史记录**，以及每轮被谁取代。
> 规则的默认方向：后一轮只修订它明确点到的章节，其余条款继续有效。
>
> 目录搬迁（docs/rounds/v5/06-layout.md）后：`docs/spec/` 是当前规范，
> `docs/rounds/vN/` 是历史轮次（只增不改）。

## 当前规范（按主题）

| 主题 | 当前文档 | 说明 |
|---|---|---|
| CLI 参数语义 | `docs/spec/cli.md` | **合并版 CLI 规范**（= docs/rounds/v5/02），当前唯一有效；v4/v2 的 02 转为历史 |
| 统计出口设计 | `docs/rounds/v5/05-algorithms.md` | SFS / `--fst-global` / f 统计 / LD 的出口形状；对齐台账在 v2/04 + v4/05 |
| 统计与过滤算法 | `docs/rounds/v2/04-algorithms.md` + `docs/rounds/v4/05-algorithms.md` | 方法与文献定义在 v2；对齐状态与缺口台账在 v4/05 |
| 架构 | `docs/spec/architecture.md`（v1 基础）+ `docs/rounds/v2/03-architecture.md` + `docs/rounds/v3/02-bcf.md` + `docs/rounds/v4/03-architecture.md` + `docs/rounds/v5/03-architecture.md` + **`docs/rounds/v6/03-architecture.md`** | v1 长期条款仍有效；BCF 在 v3/02；`unsafe` 策略与前后端边界在 v5/03；**缺陷清单 A6-1…A6-25 与三个最优方案（矩阵化 pass0、`-Ob` 记录模型、列式掩码）在 v6/03** |
| 模型（数据 / 执行 / 统计） | **`docs/rounds/v6/02-models.md`** | 基因型三种表示与 2 位编码的取舍、位平面矩阵作为全局分析的规范输入（含代码示意）、统计模型对拍与偏差台账、M6-1…M6-7 |
| 工程（CI / 测试 / 依赖 / 文档漂移） | **`docs/rounds/v6/05-engineering.md`** | 健康状态、CI 缺项、测试缺口、漂移清单、E6-1…E6-9 |
| 性能 | **`docs/rounds/v6/04-performance.md`**（当前基线）+ `docs/rounds/v5/04-performance.md` + `docs/rounds/v4/04-performance.md` | v6/04：3 次中位数基线、线程效率、deflate 级别研究、与 plink2/bcftools 对照、LD 执行模型问题、P6-1…P6-10；测量纪律在 v5/04 §4 + v6/04 §9 |
| 命名 | `docs/spec/naming.md` + `docs/rounds/v2/05-naming.md` | 版本后缀禁令（v2 05 §4）由 `ci/check-doc-drift.sh` 守卫 |
| 项目定位 | `docs/rounds/v5/08-api-surface.md` | 前后端形态（`reed-cli` 是唯一前端）与门面边界（§6） |
| 参照工具差异与对拍 | **`docs/convert-parity.md`** | README 只留参数与输出；convert 的差异清单与布局细节、f 统计 SE / Tajima's D / LD 剪枝的参照约定、全部对拍脚本总表集中在此 |
| 路线与验收门 | **`docs/rounds/v6/06-roadmap.md`** | Y0–Y7；v5/07 的 X0–X7 已落地（未完成项并入 Y6），v4/06 的 W0–W6 已完成 |
| 缺陷修复（外部审查核实） | **`docs/rounds/v7/01-review.md`** + `docs/rounds/v7/06-roadmap.md` | 对一份外部审查报告的独立核实与修复规范（D7-1…D7-29）；路线 Z0–Z4。v7 是缺陷修复轮，不重开 v6 的 Y 系列结构项 |
| 线程模型与并行预算 | **`docs/rounds/v8/01-review.md`** + `docs/rounds/v8/06-roadmap.md`** | 全部线程创建点台账与 spec §7.3 的逐条对照；rayon 三阶段与 stats LD sink 脱离 `ThreadBudget` 的修复规范（D8-1…D8-5）；路线 P0–P2。v8 是缺陷修复轮，不重开 v6/v7 的结构项 |
| 引擎错误路径与终止语义 | **`docs/rounds/v9/01-review.md`** + `docs/rounds/v9/02-architecture.md` + `docs/rounds/v9/06-roadmap.md` | 阻塞点—唤醒来源台账；reducer 侧错误/panic **永久挂死**的实测现场（`eu-stack` 5 线程全 `futex_wait`）与规约 W1/W2（唯一错误出口）；通道发送结果语义 C1；输入截断契约（文本丢弃 / index 收录 / BCF 报错）；seam 与 chunk 体的错误一致性 E1；指标口径（peak = 压缩 × `MEM_FACTOR`、reducer 上的 seam 时间漏计）；缺陷 D9-1…D9-8，路线 Q0–Q3。v9 是缺陷修复轮，不重开 v6/v7/v8 的结构项 |
| BGZF/BCF/VCF 格式兼容性 | **`docs/rounds/v10/01-review.md`** + `docs/rounds/v10/02-architecture.md` + `docs/rounds/v10/06-roadmap.md` | 以工作区 htslib/bcftools 源码为基准的格式处理对照：GT 选宽按编码值（N1，修大等位编号被 `as i16` 静默截断）；NULL 描述符非法（N2）；IDX 限界（N3）；BGZF 最小块长（N4）；文本 index 残片承载（N5，对齐 `BcfIndexJob` 的 FrameChain 模式）；多成员 gzip（N6）；PASS 注入（N7）；计数乘法 checked_mul（N8）。缺陷 D10-1…D10-8，路线 R0–R3。v10 是缺陷修复轮，不重开 v6–v9 的结构项 |
| 门面 API 公开面 | **`docs/rounds/v11/01-review.md`** + `docs/rounds/v11/02-architecture.md` + `docs/rounds/v11/06-roadmap.md` | `reed` 门面（`crates/reed/src/lib.rs` 公开面）的 API 复审，量尺是 v5/08 §6.2 与 README「Library use」：默认 `Levels` 产出 0 行（A11-1/N1）；`run_stats` 按值取请求使 README 示例不编译（A11-2/N2）；`to_args` 对本仓 `default_kinds` panic（A11-3/N3）；门面横跨六种错误类型、`Error` payload 不可命名（A11-4/5/N4/N5；**N5 已定方案②b：`Error` 变体封口**）；`run_reports`/报告类型不在根（A11-6/N6）；`is_empty` 漏 `summary`（A11-7/N7）；死再导出、`OpenedInput` 内部量、测试钩子、缺内存态 `Populations`、命名与重复 helper（A11-8…13/N8…N13）。证据通道：实跑 / 下游探针编译 / 源码阅读。缺陷 A11-1…A11-13，路线 R0–R3。v11 是 API 复审轮，不重开 v5/08 §6 的分工与 v6–v10 的结构项 |
| INFO 字段兼容性 | **`docs/rounds/v12/01-review.md`** + `docs/rounds/v12/02-architecture.md` + `docs/rounds/v12/06-roadmap.md` | v10 结构层对照向 INFO **语义层**的延伸，真值门 bcftools 1.24 / 工作区 htslib：**读侧已逐字节一致**（差分夹具 + `-q :type` 五种 + 有 GT 的子集更新）；写侧值编码偏差——超 int32 整数 `-Ob` 回绕（D12-3/N3，htslib 置 missing）、Flag 带值丢值且 `DB=0` 语义反转（D12-4/N4）、String 值 `.` 变空串（D12-5/N4）、打包计数溢出写坏 BCF（D12-9/N8）；过滤语义偏差——sketch 策略把尾部表达式放在无记录字节的上下文里求值致**静默丢记录**（D12-1/N1）、无 GT 记录的 AC/AN 被子集更新清零（D12-2/N2）、`-q/-Q` 在 cohort 变化后仍用陈旧 INFO 且缺 cli.md:78 承诺的警告（D12-6/N5）、表达式层拿不到 header 类型致 Flag 比较恒假 / 向量只比首元素 / 裸标签真值相反（D12-7/N6）、`quantifier = "any"` 不看比较结果（D12-8/N7）、`ALT='.'` 计 1 个 ALT（D12-10/N9）；契约可见性——`--fill-tags` 不校验标签（D12-11/N10）、默认 `-x` 陈旧清单在 help/README/cli.md 均无提示（D12-12/N11）；共性代码收口——AC/AN 读取三份实现（D12-13/N12）、INFO 数值解析三份致 **reed 内部自相矛盾**（同记录过滤判定值≠`-Ob`写出值≠`-Ov`写出值，与 htslib 23 处差异；D12-14/N13）、`edit_vcf_line` 死代码（D12-15/N14）。**§0.5 记 10 个共性代码根因簇（12/15 条缺陷可归因；两个元模式：多份实现各自漂移、元数据存在但消费者没查）**。缺陷 D12-1…D12-15，路线 R0–R3。v12 是缺陷修复轮，与 v11 正交（无共享文件）可并行 |
| `convert` 子命令（Phylo 多格式 + SFS 出口迁移 + hapmap↔VCF 双向，已落地） | **`docs/rounds/v13/README.md`** + `docs/rounds/v13/01-review.md` + `docs/rounds/v13/02-architecture.md` + `docs/rounds/v13/06-roadmap.md` | **新功能轮**（非缺陷修复；**已落地**，对拍门全绿）：`reed convert` 子命令，六组转换目标——SFS（从 `stats --sfs` 等价迁移，逐字节等价契约、stats 入口不删，翻译层共享单一 `build_sfs_request`、折叠单一实现）+ Phylo 四格式（vcf2phylip 2.9 的 PHYLIP/FASTA/NEXUS/二元 NEXUS，三格式共享单一列模型）+ **`--012`（vcftools 0/1/2 矩阵）+ `--hapmap`（TASSEL .hmp.txt）+ `--msa`（vcf2msa 参考序列全位点 MSA，第二轮增补）+ `--vcf`（hapmap→VCF 反向，第三轮增补：用户点名"hapmap 应该可以相互转换"，record→record 流式翻译复用 `-O` 写出）**；全部**并行化**（复用融合式 chunk 流水线：worker 内一次行解析多目标消费、partial 加法合并按 chunk 序、天然确定性；vcftools/vcf2phylip 的每样本临时文件转置不复现）；门面 `ConvertRequest`/`run_convert`（七种 `ConvertTarget`，请求按引用、`to_args` 可失败，v11 纪律沿用）；参照语义逐行号钉死（AMBIG 62 项/GEN_BIN、`<NON_REF>` 与 `-`→`*` 处理序、min-samples-locus 数样本不数等位、012 的半调用分支、hapmap 11 固定列、msa 的无记录=REF 与 indel 填充回退路径、hmp→VCF 的 `ExportUtils.writeToVCF` 真值）；vcftools 其余格式（plink/IMPUTE/BEAGLE/LDhat/LDhelmet）台账钉死在 01 §1b.4 但不进首发轮（C13-5/6）；有意分歧 12 条（B13-1…12）；待拍板 6 条（C13-1…6，默认均已落地）；vcf2phylip 与 vcf2msa 为 GPLv3——**只实现行为不复制代码**，Credits 引 Ortiz 2019/vcftools/TASSEL/Chafin。约束 D13-1…8（含 012/hapmap/msa 共享行视图单一实现、hmp→VCF 无新压缩路径），路线 C0–C4（真值门 fixture → phylip 垂直切片 → 多格式渲染 → 012/hapmap/msa 三连 → SFS 迁移 + 并行验证 + 文档收口 → hapmap↔VCF 双向）。v13 与 v12 正交（不触碰 INFO 路径） |
| Y5 引擎落地规格 | **`docs/spec/v6_y5/README.md`** + `docs/spec/v6_y5/TODO.md` | A6-1…A6-4 的落地形态（中止原子位 / 许可 drop guard / 接缝上限 / 顺序模式 / 计量）、与 roadmap 附录的差异表、门与执行清单 |
| Y6 前后端落地规格 | **`docs/spec/v6_y6/README.md`** + `docs/spec/v6_y6/CODE.md` + `docs/spec/v6_y6/TODO.md`；批 D 细化 `docs/spec/v6_y6/batch_d/`（README/TODO/CODE） | F6-13/F6-19 的落地形态（reed::Error / 请求内化 / job.rs 拆分 / C-5 / 联合 SFS / `--bed-file`）、与 roadmap 与 08 §6.2 的差异表、门与执行清单；批 D 的真值表/集成点/退出码契约 |
| 实施状态 | **`docs/rounds/v6/STATUS.md`**；缺陷修复轮各自的状态在 **`docs/rounds/v7/STATUS.md`**（Z0–Z4）、**`docs/rounds/v8/STATUS.md`**（P0–P2，含 D8-6）、**`docs/rounds/v9/STATUS.md`**（Q0–Q3，起始全部未修 + 复审基线与实测记录）、**`docs/rounds/v10/STATUS.md`**（R0–R3，起始全部未修 + 复审基线代码证据摘录）、**`docs/rounds/v11/STATUS.md`**（R0–R3，起始全部未修 + 复审基线实测/编译证据）、**`docs/rounds/v12/STATUS.md`**（R0–R3，起始全部未修 + 复审基线实测输出摘录与三条"已一致"基线）、**`docs/rounds/v13/STATUS.md`**（C0–C4，新功能轮：起始全部未开始 + 基线表已实测） | 历史：`docs/rounds/v5/STATUS.md`（X0–X7）、`v4/STATUS.md`（W0–W4）、`v3/STATUS.md`（V0–V6）、`v2/STATUS.md`（M0–M5b/R0–R8） |

## 历史文档（只作记录，不再修订）

| 文档 | 状态 |
|---|---|
| `docs/rounds/v1/design.md` | v1 设计（原 `DESIGN.md`）；除后续轮次明确推翻的章节外仍然有效 |
| `docs/spec/architecture.md` | v1 规范；除 docs/rounds/v2/03 明确推翻的章节外仍然有效 |
| `docs/spec/filter.md` | reed-filter 目标架构（v2，2026-10-05）；执行模型已按 docs/rounds/v2/03 与 v3 落地 |
| `docs/spec/naming.md` | 命名规范；稳定性状态已按 docs/rounds/v5/08 §6 更新 |
| `docs/spec/statistics-map.md` | 统计量到参考实现的映射；对拍结果见各轮 STATUS |
| `docs/rounds/v2/` | 第二轮复审与修订规范；被 v3（BCF/性能/清理）与 v4（参数方向）增量修订 |
| `docs/rounds/v3/` | 第三轮复审与 BCF/性能/清理；V0–V6 已完成，状态在 `v3/STATUS.md` |
| `docs/rounds/v4/` | 第四轮：复审（01）、参数语义规范回退版（02）、架构残留（03）、性能（04）、算法台账（05）、路线 W0–W6（06）、目录树（07）、状态；W0–W4 已落地 |
| `docs/rounds/v5/` | 第五轮：复审（01）、合并版 CLI 规范（02）、架构与前后端边界（03）、性能基线（04）、W5/W6 出口设计（05）、目录树执行版（06）、路线 X0–X7（07）、项目定位（08）、状态 |

## 每轮明确取代的章节（增量，不是整份作废）

| 被取代的条款 | 取代者 |
|---|---|
| v2 02 §3.1 `--mind` / §3.4 `--geno`（plink2 丢弃方向） | v4 02 §3：**名字保留**，语义精确定义（等位计、掩码后、精确边界） |
| v4 02 首版 §5.1–§5.4、§6、§7 的改名（`--max-missing*`） | v4 02（回退版）同名小节：撤回改名 |
| v2 02 §5.2 `--fst-pops A,B,C`（多群体合并 Fst） | v4 05 §2.1 / v5 05 §3.2：作废，改为 `--fst` + 布尔开关 `--fst-global`（默认成对） |
| v2 02（整体）与 v4 02（回退版） | **docs/spec/cli.md**（v5 02 合并版）是当前有效的 CLI 规范 |
| v4 05 §2–§2.1、06 的 W5/W6、07-layout | v5 05（出口设计）、07（X0–X7）、06（目录树执行版） |
| `docs/spec/architecture.md` D12（Python 绑定 M8 之后再做） | v5 08 §6.5：**不做**；项目定调为前后端形态，`reed-cli` 是唯一前端 |
| v2 05 §2.3"`--max-missing` 方向与 `--geno` 相反，不提供别名" | v4 02 §8：结论保留（不提供别名）；`--max-missing` 已不是 reed 的参数 |
| v3 04 §C3"`TagTable` 未做" | `docs/rounds/v3/STATUS.md` 更正为**已做**（v4 01 F4-2） |
| `docs/spec/architecture.md` §13.17 的漏斗 id 示例 | v4 02 §5.3 定名（`after_mask.missing`、`sample.mind`）；§13.17 已就地回改 |
| `docs/spec/architecture.md` §7.7「超出时报 `EngineError::MemBudgetExceeded`」 | **doc_v9 D9-7（待落地）**：改为背压表述——超出时 reader 阻塞，`MemBudgetExceeded` 只用于"单个 chunk 的足迹超过全部预算"。依据是 `docs/rounds/v2/03-architecture.md` §4.2，实现早已如此，spec 未回改 |
| `docs/spec/architecture.md` §7.5「所有阻塞点都用 `select!` 等待数据和 abort」 | **doc_v9 02 §2（待落地）**：追加规约 W2——任何线程的错误与 panic 必须先经 `FirstError` 记录（记录即触发 abort），不得绕过它直接返回。§7.5 的机制本身正确，缺的是"没有旁路"这一条 |
| `v8/STATUS.md:22-23`「`stitcher_reassembles_and_caps_fragments` 失败是预存在的，与本轮无关」 | **doc_v9 D9-8 更正**：预存在属实，但"无关"是误判——它是 D9-3（`NewlineStitcher` 缺 finish 语义）的症状，且暴露了 release-only 验收门对 `debug_assert` 的盲区 |
| README「Library use」的 `StatsRequest::to_args()` 写法、漏列 `run_reports`、及各驱动请求参数的按值/按引用不一 | **doc_v11 N3/N2/N6 更正**：`to_args(&group_names)`、驱动清单补 `run_reports`、`run_stats` 请求改按引用（README 示例纳入编译门） |
| `docs/rounds/v5/08-api-surface.md` §6.2 的后端公开面清单（未含 `run_reports`/报告根导出，未注 `_opened`/`OpenedInput` 例外） | **doc_v11 N6/N9 增补**：报告类型与 `run_reports` 置于根；`_opened` 与 `OpenedInput` 列为批准例外或 opaque 化 |
| `docs/spec/cli.md` §3.5 的 `--fill-tags` 行（承诺 AF/F_MISSING，实现只有 AC/AN/NS/MAF/HWE 且不校验标签） | **doc_v12 N10 更正（待落地）**：支持集收敛为已实现五项，AF/F_MISSING 标【待实现】；未知标签报 usage 错（bcftools `+fill-tags` 同款） |
| `docs/spec/cli.md` §3.5 的 `-x` 行与 §3.2 的 `-q/-Q` 行（未写默认陈旧清单；未写 cohort 门与警告） | **doc_v12 N11/N5 增补（待落地）**：`-x` 默认删 8 标签、`-x ''` 关闭；`-q/-Q` 在样本子集/掩码后按 GT 重算（bcftools `vcfview.c:341-350`、`filter.c:1471` 同款门），INFO 优先真正生效时打"陈旧 INFO"警告（兑现 cli.md:78 的承诺） |
| `docs/spec/cli.md` §5 的表达式支持清单（Flag 比较、向量 any、裸标签真值、INFO 数值解析差异未记） | **doc_v12 N6/N13 增补（待落地）**：Flag 按存在性 1/0（`doc/bcftools.txt:4067`）；Number≠1 按 any-element（`CMP_VECTORS`）；裸标签恒假或编译期拒绝【决】；新增与 htslib 的解析差异表（`+5`/`5.0`/`1e3`/`0x10`/Σac>an） |
| `docs/spec/cli.md` §3.2 的 `-m/-M` 行（`ALT='.'` 计 1 个 ALT、`-M 1` 因"0=关闭"而整条规则消失） | **doc_v12 N9 更正（待落地）**：`ALT='.'` → 0 个 ALT（bcftools `n_allele`）；`max_alt_alleles` 改 `Option<u32>`，`-M 1` 语义为"只留不变位点" |
| `docs/spec/architecture.md` §11.4 的 `SketchThenEmit` 触发条件（"有 `max_missing` 规则，且 sketch 估计大小 ≤ 内存预算"，未排除尾部表达式/行字段依赖） | **doc_v12 N1 增补（待落地）**：sketch 前置"无尾部表达式且 `-q/-Q` 不依赖 INFO"；`--mind + --indep-pairwise + 尾部表达式`编译期报错。§11.4 的三策略等价契约本身正确，缺的是这条前置 |

## 一致性守卫

`bash ci/check-doc-drift.sh` 校验"STATUS 声明 vs 代码事实"的断言：C3 无 noodles、无版本后缀类型、
`--geno`/`--mind` 是真实参数、`--max-missing` 只出现在报错文本里、`reed index` 拒绝 stdin、
`Verbosity` 单一实现、恰两处 `unsafe`、无 `sample.max_missing`/`--fst-pops`/`reed qc`、
前端只依赖门面、后端不打 stderr、门面无 `pub mod`、README 无 v0.2 示例、无 noodles re-export。
文档声明"未做"的项应在该脚本里有一条守卫。
doc_v11（门面 API）落地后追加其 01-review §3 的守卫：`Levels` 默认 `window=true`、
`run_stats` 请求按引用、`to_args` 可失败、根含 `run_reports`/`FilterReport`/`ReportsOutcome`、
无 `LdWindow`/`bgzf`/`header` 死再导出、README 库示例在编译门内。
doc_v12（INFO 兼容性）落地后追加其 01-review §3 的守卫 1–14：策略覆盖尊重尾部
表达式、`fill_info_set` 有 `has_gt` 门、`parse_int_vector` 引用 INT32 值域常量、
Flag/String 编码收口（无 `if v == b"."` 空串特判）、`MinAfInfo` 推入含 cohort 门
且 `-q` 有陈旧警告、`expr::check` 签名含类型信息、`any` 分支含 `cmp.apply`、
打包计数有上限 helper、head 的 ALT 计数走共享 helper、`--fill-tags` 未知标签
报 usage 错、`--help` 的 `-x` 条目含默认清单、AC/AN 读取单一实现、INFO 数值
解析单一实现、`line_edit.rs` 无 `LineRef::parse(line.line)`。
doc_v13（`convert` 子命令）已落地并追加其 01-review §6 的守卫 1–18：`Command`
恰四变体且 `--help` 列 convert、SFS 请求构建单一共享函数（stats/convert
双入口同源）、`push_sfs` 全仓唯一实现、`reed-core` 的 iupac 模块（AMBIG 62
项 / GEN_BIN 10 项）、convert 无 `--mind`、`--phylip` 默认关、resolve-iupac
确定性双线程测试存在、`run_convert` 请求按引用且根导出、
`gold/popgen/run_vcf2phylip.sh` 存在、convert 渲染层无临时文件 API、
行级 GT→等位索引解析全仓单一实现（phylo/012/hapmap/msa 四目标共用）、
`--012` 半调用置 -1 的单测存在、`--hapmap` 的 alleles 列来自 REF+ALT
声明集、AMBIG 未命中走跳过路径、hapmap 字母表收敛有测试、
`--vcf` 反向目标存在且与 VCF 输入目标互斥（门规则 6，根导出）、
`ambig_expand` 逆表与往返保真 e2e、hapmap→VCF 头无 `##contig` 且
`##reed_*` 可关（B13-11/12）。
