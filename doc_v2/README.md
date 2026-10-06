# reed doc_v2：复审 + 修订规范（2026-10-06）

> 本目录是 `docs/architecture-spec.md`（v1）之后的**第二版规范**。
> 它由两部分组成：对已完成的 M0–M5b 的复审，以及在复审基础上的修订规范（允许"爆破重构"）。
> 本目录不修改任何代码，所有改动建议都以规范形式给出，由实现者执行。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | M0–M5b 逐里程碑复审：实测性能、正确性对拍、缺陷清单（带 `file:line`、严重度、修法） |
| [02-cli.md](02-cli.md) | **CLI 规范（推翻 v1 §18）**：`reed filter / stats / index`；参数名取自 bcftools/plink2/vcftools/pixy 且语义一致；统计参数以科学方法命名（`--pi --dxy --fst wc84|hudson --theta --tajima-d --pbs --sfs --freq --hardy --het --missing --summary`），缺失数据默认采用 pixy 的方法，只有 `--variant-only` 一个声明开关 |
| [03-architecture.md](03-architecture.md) | 架构重构：删除 legacy、统一记录解码、crate 边界纠正、引擎修复、filter/stats/qc 驱动拆分、目标文件树 |
| [04-algorithms.md](04-algorithms.md) | **经典方法 + pixy 的缺失数据处理**：撤下 reed 自创语义；各软件在经典统计量上的异同（源码证据）；reed 的唯一定义与文献；窗口网格；对拍方案（软件只在测试中出现） |
| [05-naming.md](05-naming.md) | 命名规范增量：对 v1 §13 的修订（CLI 参数与别名规则、按科学方法命名的模块与枚举、输出文件名） |
| [06-roadmap.md](06-roadmap.md) | 修复与后续里程碑：P0 必修清单、R1–R8 重构步骤、每步验收门、流程纪律 |

## 与 v1 的关系

- v1 中**仍然有效**的部分：设计原则 §2（第 1 条除外：CLI 不再是"Spec 的薄翻译"，而是一等公民的参数语言）、
  数据模型 §5、BAPP §6.2、引擎骨架 §7、累加器与分箱 §8、统计算法 §9、QC 目录 §10、测试体系 §16、性能工程 §17。
- v1 中**被本目录取代**的部分：
  - §11.1 TOML Spec、§18 CLI → `02-cli.md`；
  - §4 文件树、§19 迁移路线 → `03-architecture.md`、`06-roadmap.md`；
  - §13.14 / §13.15 CLI 与配置键命名 → `05-naming.md`；
  - §20 决策 D6（TOML/JSON 规范文件）→ **废止**；D14 → 重新表述（见 02-cli.md §1）；D4（投影默认）→ **撤回**；D5（PBS 位点集合）→ 由 allel 语义决定（per-pair）；
  - v1 §5.8 清单 v2、§9 中所有 reed 自有估计器定义、DESIGN §4 的清单门 → **撤下**（04-algorithms.md §1）。

## 结论速览

| 维度 | 结论 |
|---|---|
| 性能 | **远超目标**。全量（6.6 GB、348 样本 × 3.3M 位点、56 核）：stats 2.65 s（目标 ≤ 30 s），filter 14.2 s（目标 ≤ 60 s；bcftools 约 40 min，plink2 约 2.5 min），GT 扫描 3.14 s，QC summary 2.38 s |
| 过滤正确性 | **通过**。1,768,759 个位点 × 348 样本的 GT、AC/AN/NS 与 bcftools gold 的 md5 完全一致；MAF/HWE 抽查区域逐字节一致；`bgzip -t` 通过；tbi 区域查询正常 |
| 统计正确性 | **未对齐主流**。π/dxy、SFS/θ/Tajima（complete-case）、窗口网格都是 reed 自定义语义，与 pixy/vcftools/allel 都无法逐数对上；另有输入契约缺陷：缺 contig 长度时编造、不校验位置顺序、群体文件首行被跳过、输出列用组下标 |
| 测试状态 | **红**。`reed-engine/tests/pipeline.rs::output_is_invariant_across_thread_counts` 失败，panic 位于 `reorder.rs:58`（5 次运行失败 4 次），正是复审发现的重排窗口缺陷 |
| 结构 | **双实现并存**。legacy（v0.2）与 v2 两套 reader、两套 filter、两个二进制同时存在；GT 文本解析器有 5 份；`reed-filter` 依赖 clap/toml/noodles，不再是纯策略层；`reed/src/filter.rs` 1699 行 |
| CLI | **推翻**。改为 bcftools/plink2/vcftools/pixy 风格的参数语言，删除 `--spec`/`--conf`/TOML（02-cli.md） |
| 方向 | **只做架构创新，不做算法创新**。统计量都是有文献的经典方法，缺失数据统一用 pixy（Korunes & Samuk 2021）的方法；reed 的价值是多核单遍、全字段保真、一条命令完成过滤加统计 |

## 决策状态

| 决策 | 状态 |
|---|---|
| CLI 不用 conf/TOML | 已定（用户，2026-10-06） |
| 只做架构创新，统计与主流对齐 | 已定（用户，2026-10-06） |
| CLI-1 命令划分 | 已定：`filter / stats / index`，统计合并进 `stats`（02-cli.md §1） |
| A-1 投影 n、A-2 同位置多记录 | 已消解（04-algorithms.md §3、§6） |
| 统计参数命名 | 已定（用户）：只写科学方法，不写软件名；`--compat` 撤回 |
| 缺失数据默认 | 已定（用户）：pixy 的方法；只有变异位点的输入用 `--variant-only` 补全分母（与 vcftools `--window-pi` 的做法相同） |

## 标记约定

【规】必须；【荐】推荐；【决】需要你拍板；【据】有实测或源码证据；**P0/P1/P2** 为修复优先级。
