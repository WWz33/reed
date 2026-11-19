# doc_v30：sgkit API 覆盖、科学符合性与 API 工程质量审计

> 立项 2026-11-11。回答三个问题：**（1）相对 sgkit 公开 API，reed 实现了多少；
> （2）已实现的部分在科学上是否站得住；（3）作为可被外部使用的 Rust API 库，工程上
> 是否"优雅"**。第（1）（2）问不要求与 Python/xarray 接口逐一对齐，CLI/库等价即可，
> 但公式与统计语义必须成立；第（3）问只审 API 形状，不改算法。
>
> 本轮是**审计轮**：只记录证据与缺口，不改代码。修复项落在
> [02-roadmap.md](02-roadmap.md)，实施与否由后续轮次决定。
>
> API 工程质量审计已并入本文件，不单独成轮：
> 第一部分 = V30-01…09，第二部分 = V30-Q01…Q36，第三部分 = V30-R01…R14（根因审查），
> 第四部分 = V30-T01…T20（工具化全库审计；2026-11-11 追加）。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | **§1–§5**：56 项对照表、科学符合性、真缺陷、实现但 sgkit 无 API 的能力；**§6–§13**：36 条 API 工程质量（V30-Q）；**§14–§16**：14 条根因审查（V30-R）；**§17–§20**：20 条工具化全库审计（V30-T）与既有条目修正表 |
| [02-roadmap.md](02-roadmap.md) | **第一部分**修复优先级；**第二部分** Q-A/Q-B/Q-C；**第三部分** R-A/R-B/R-C/R-D；**第四部分** T 的 S-A…S-E 批次与 8 条验收门；**第五部分** P0–P4 综合方案 |
| [03-synthesis.md](03-synthesis.md) | **综合诊断**：79 条 → 5 个根因 + 1 个产品裁决；P0–P4 根本方案、12 条工具验收门、关闭规则 |
| [STATUS.md](STATUS.md) | 审计状态与复现记录（第二部分见 §5，第三部分见 §6，第四部分见 §7，综合诊断见 §8） |
| [04-recheck.md](04-recheck.md) | **独立复核**（2026-11-11 晚）：三条修法错（V30-01 的 `e02`、V30-04 的"对齐 ADMIXTOOLS"、V30-R02 的"同源"门）；R01/R12/发布类条目的定性修正；**新缺陷 N-1 `--skip-empty` 截断窗口（数值错误）** 与 N-2…N-5；滑窗复核 |
| [04-implementation.md](04-implementation.md) | **实施记录**（2026-11-12）：P0 契约 ADR 与 CI 门、P1 R01/R02/R03/R04/R12、P2 R05/R06/R07/R08/R09/R10/R11；验证命令与结果；未完成清单 |
| [05-review-followup.md](05-review-followup.md) | **实施后复核与第二批规范**（2026-11-13）：确认 N-1、PI_HAT、`-I`、fill-tags 顺序、LD idle 和 phylo 口径等问题，并给出修复门 |
| [06-final-audit.md](06-final-audit.md) | **当前权威收口**（2026-11-19）：逐项确认已修复项、剩余工程风险、IO/API 优化建议、noodles 路线、main/dev 发布流程和需要拍板的事项 |

> **当前状态请以 [06-final-audit.md](06-final-audit.md) 为准。** `04-recheck.md` 和 `05-review-followup.md` 是修复前的证据与规范；`STATUS.md` 的早期表格是历史基线，后续追加章节及最终收口文件才反映当前工作树。

> **复核（[04-recheck.md](04-recheck.md)）**：V30-01 的现象成立但给出的 `e02=(p²+q²)²` 错（应为 `p⁴+q⁴+4p²q²`，且只有完整有限样本系数能与 plink 逐位相同）；V30-04 的"ADMIXTOOLS 用 Σ 尺度"与 `qpDstat.c`/`qp3Pop.c`/ADMIXTOOLS 2 源码不符；V30-R02 的 `.afreq`/`.stats` AF 差异是两边各自对齐 plink2 与 bcftools。另发现 **N-1：`--skip-empty` 让第一个非空窗口只算到 1 个位点**（真实数据上 π 差 67 倍，`stats/sink.rs:756-767`，v27 引入）。

## 结论（1）（2）：覆盖与科学符合性

**（1）覆盖**：sgkit 公开 API **56 项**（`docs/api.rst` 54 项 + `__init__.py` 额外导出的
`pca`、`regenie_loco_regression`）。reed **不能算实现了"大部分"**，但分布极不均匀：

| 判定 | 数量 | 占比 |
|---|---|---|
| 完整或等价 | 12 | 21% |
| 部分 / 换算法实现 | 13 | 23% |
| 未实现 | 27 | 48% |
| 不适用（Zarr/Dataset/BGEN 数据层） | 4 | 7% |

只看"硬调用 VCF 上的统计内核"（56 项剔除关联 4、系谱 8、概率/倍性 5、窗口/展示与模拟 4、
Zarr/Dataset N/A 4，剩 31 项），覆盖率 **25/31 ≈ 81%**：缺的 27 项里有 12 项属于
GWAS/regenie/genee 与系谱两块独立产品（reed 完全没有），5 项属于概率/混合倍性数据层，
其余是 display、window_by_variant 等边角。统计内核本身只剩 6 项未实现：
`pc_relate`、`Weir_Goudet_beta`、`Garud_H`、`individual_heterozygosity`、
`invert_relationship_matrix`、`count_call_alleles`。

**（2）科学符合性**：已实现项**绝大多数站得住**——公式对标 pixy / scikit-allel / plink2 /
vcftools / ADMIXTOOLS 等正典，差异多为有据可查的口径选择。但有 **3 条真实缺陷**：

- **V30-01（P0）** `reed distance --type pihat` 的期望 IBS 系数错位 —— 实测 60378 对里
  **60346 对**输出退化的 `Z0=0 Z1=1 Z2=0 PI_HAT=0.5`，与 plink 1.9 源码逐项不符。
  这不是 v17 记录的"缺 Nei 修正"，而是系数排列本身错误。
- **V30-02（P1）** `filter -I` 的文档承诺与行为相反：`filter.md` 写 "Leave INFO/AC and INFO/AN alone"，
  实测 `-I` 会把 AC/AN **删掉**（要保留必须同时 `-x ''`）；`--help` 的 "Do not update …" 措辞更弱，
  但同样不完整（无样本变化时也删）。
- **V30-03（P1）** 多等位在 `stats` 被静默折叠（`CollapseAlt`）、在 `popg` 被跳过、
  在 `filter --fill-tags` 报错 —— 同一份输入三个命令三种行为。

另有 5 条文档层问题（f-stat 的 Σ 量纲、`--f4-ratio` 角色序、f3\* 过度承诺、
MSA 缺失→REF、`identity_by_state` 命名），见 01 §2.2/§2.3。

## 结论（3）：API 工程质量（V30-Q01…Q36）

把"优雅"拆成 8 个可检查性质：可发现、可命名、可构造、可消费、可诊断、可演化、可验证、
可发布（定义与"本仓现状"表见 01 §6）。共记 **36 条**：**6 P0** / 17 P1 / 13 P2。
结论：**骨架是库的，脸面还不是。**

| 严重度 | 数量 | 含义 |
|---|---|---|
| P0 | 6 | 作为"外部可用库"不成立/会直接误导用户 |
| P1 | 17 | 设计不一致或契约缺失，造成误用与维护成本 |
| P2 | 13 | 文档/工程卫生 |

### 六条 P0

1. **V30-Q01 门面整树再导出**：`pub use reed_distance as distance` 等 11 个模块/整 crate
   别名，把 ≈410 个下层 pub 项变成门面 API；`reed::reduce` 全仓零引用。
2. **V30-Q02 公开签名含不可命名类型**：`StatOutput`、`BcfDict`、`QcReportAcc` 字段类型都不在
   门面清单；本仓 **13 个**集成测试已绕过门面直接 `use reed_stats::...`。
3. **V30-Q07 `PcaRequest`/`DapcRequest` 的 `threads`/`memory` 是死字段**：CLI 写入、驱动从不读取，
   真实预算走 `&ThreadBudget`；`PcaRequest::default()` 还会在当前目录静默写 `out.*`。
4. **V30-Q16 `run_distance` 缺失**：前端 472 行承担矩阵选择/bootstrap/写盘/错误分类，
   是 `v5/08 §6.1` 唯一没落地的命令，也是整 crate 别名的根因。
5. **V30-Q24 零用户级 doctest**：9/10 crate 的 lib.rs 零示例；全仓仅 3 个可运行 doctest、
   2 个 `# Examples`；examples/ 全是内部基准。
6. **V30-Q29 发布定位矛盾**：门面 `reed` `publish=false`；反而 distance/reduce/graph 未关发布，
   但 path 依赖缺 version，`cargo package` 必失败；全仓 0 个 repository/readme/keywords/docs.rs。

### 其他最值得先看的

- **V30-Q08** 请求全 `pub` 字段 + Default：`PopgenRequest::default()` 静默 0 行；
  `FilterSpec` 的 policy/alt 数可自相矛盾。
- **V30-Q11** 错误分类不可编程：`ErrorKind`/`FilterError` 非 `non_exhaustive`；`Foreign` 丢类型。
- **V30-Q18** `run_stats` 文件进文件出，与 popgen/filter 的"数据返回"不一致。
- **V30-Q31** `reed-reduce --no-default-features` 编译失败（E0063）——零依赖承诺是假的。
- **V30-Q32/Q33** CI 错位：main 与 PR→main 只 `cargo build`；公开树被 `publish.sh` 剥离
  tests/docs（main 259 文件 0 测试；dev 579 文件 65 测试）。

### 正面结论（修复时不要改坏）

0 个 `pub mod`、sealed `Error`、后端 0 处 `println!/eprintln!`、602 个 `#[must_use]`、
L0–L5 单向分层、`reed-stats`/`reed-filter` 纯策略无 I/O、98 条 doc-drift + 16 条 structure 守卫。

## 结论（4）：根因审查（V30-R01…R14）

按四项根因审查（语义归属 / 能力接口 / 生命周期与所有权 / 变更传播）把"同一条规则反复在多模块
修正、加一种输入要改好几层、边界判断不断复制"归到根因：**14 条，4 P0 / 10 P1**。

| 轴 | 条目 | 最重的根因 |
|---|---|---|
| 语义归属 | R01–R04 | **R01** GT/等位语义无主（5 份 ALT 计数、convert 自定"有数据"，半缺失判定与 `stats` 分叉）；**R02** AC/AN 第二份实现（debug panic / release 回绕 / 同命令两个 AC）；**R03** `pos0` 四处转换、三种哨兵语义 |
| 能力接口 | R05–R07 | **R05** "可重读"不是类型能力（stdin + two-pass 静默单趟）；**R06** 同一能力缺失 exit 1/1/2 三种类别；**R07** source 三种获取方式 + BCF 靠 `Option` 推断 |
| 生命周期与所有权 | R08–R11 | **R08** prune 线程失败路径 detach；**R09** `--memory` 不参与策略选择；**R10** `TempSpill` 只拥有路径；**R11** `OrderedOutput::adopting` 可覆盖 owner（写入静默丢失） |
| 变更传播 | R12–R14 | **R12（P0）** `reed-filter`（自称 L3）依赖 `reed-io`，违反 `architecture.md:1843`，是 R01/R02/R04 多份实现的共同原因；**R13** 一个 `--fill-tags` 标签 6 份枚举；**R14** 每条过滤规则要在 3 个 `EditPlan` 组装点重复 |

证据：R01/R02/R05/R06 本轮实测；R03/R04 为代码路径（端到端由子代理实测）；R08–R11 为代码推理
（故障注入测试尚未编写）；R12–R14 有 `git show --stat` 提交证据。详见 [STATUS.md](STATUS.md) §6。

## 结论（5）：工具化全库审计（V30-T01…T20）

六个子代理用已安装的六件工具覆盖 11 个 crate，主代理复核数字：**20 条，3 P0 / 12 P1 / 5 P2**。

**三条 P0：**

1. **T01 semver 门只覆盖 3/11**：`cargo semver-checks --workspace` 对 7 个 `publish=false` crate 打印
   `Skipping`，`reed-cli` 无 lib；实际检查的 `distance/graph/reduce` 恰是 0.1.0/2021/无 lints 的三个。
   显式 `-p reed` 能跑，说明不是工具限制，是发布策略让门失明。
2. **T02 公开签名仍含不可命名类型**：`StatOutput`、`BcfDict` 之外，**新增 `EngineError`** 出现在
   `write_fstats/write_ld/write_projected_sfs/write_joint_sfs` 的返回类型；`QcReportAcc` 已导出但字段类型未导出（半修复）。
3. **T03 门面整树再导出规模被低估**：源 `pub use` 30 条，展开 32 条；**20 条模块再导出**
   （13 crate/模块别名 + 7 writers），不是前文的 11 条；`reed::reduce/bgzf/header` 零内部引用；
   `cargo public-api` 不展开模块再导出，且没有 `ci/api` 快照。

**其他最值得先看的：** T04 三个"可发布" crate `cargo package` 全失败（两种根因）+ 元数据 11/11 缺；
T05 noodles 10 项 / serde 37 类型进入公共契约；T07 仅 `reed-reduce` 有 feature，`--no-default-features`
仍 E0063；T09 `reed-filter→reed-io` 唯一真实越层 + 2 处模块环；T12 rustdoc **5 crate/19 处**红；
T13 仅 6 个 doctest（3 执行）、强制 `missing_docs` 暴露 252 项；T16 CI 在 main/PR→main 只 build；
T17 lints 只覆盖 8/11，强制 reduce 有 28 条会失败。

**工具结论**：`cargo-check-external-types 0.5.0` 必须固定 `nightly-2026-03-20`（默认 nightly 的
rustdoc JSON v61 不兼容）；`cargo public-api` 只列再导出语句，可达面是下界；`cargo hack` 证明
`--all-features` 绿而空 feature 集红。修正表见 [01-review.md](01-review.md) §19（含删除 Q30 的
"reed-graph 无 rust-version"错误）。

## 综合诊断：79 条 → 5 个根因 + 1 个裁决

逐条 issue 不是工作分解。把 v30 的 79 条（9 + 36 + 14 + 20）按"删掉哪个所有者会让整簇症状
同时消失"归并，得到 5 个系统性根因 + 1 个未裁决的产品问题：

| 根因 | 缺的所有者 | 条数 |
|---|---|---:|
| **A 契约面** | 门面没有设计出来的 API | 17 |
| **B 语义** | 领域规则没有唯一类型/函数 | 15 |
| **C 状态与资源** | 能力/预算/临时文件/线程没有 owner | 18 |
| **D 演进** | 没有可执行的契约门 | 16 |
| **E 文档与结构** | 文档/规范不是派生品 | 9 |
| **F 产品裁决** | 未决定发布面与稳定面 | 4 |

**根本方案** = P0 契约裁决 → P1 语义内核 → P2 能力与状态类型 → P3 门面收口 + 六工具门 →
P4 派生文档与结构。详见 [03-synthesis.md](03-synthesis.md)（含 79→6 的完整归属表、因果链、
工具验收矩阵与关闭规则）。

## 边界

- **没有 sgkit 数值夹具**：所有"与 sgkit 等价"都是**公式层**结论；现有对拍是
  pixy/scikit-allel/plink2/vcftools/bcftools/ADMIXTOOLS。要对外声称兼容 sgkit，
  需另加 sgkit 夹具（见 02 V30-09）。
- 第二部分（API 工程质量）36 条中 **9 条本次实测复现**（V30-Q01/02/07/16/29/31/33/34/35：
  `rg` 零引用/零命中、`cargo check --no-default-features` 失败、`git ls-tree` 树计数、
  `git check-ignore`、`cargo doc -D warnings` exit 101），其余为确定性代码/清单证据；
  每条末尾标注证据级。见 [STATUS.md](STATUS.md) §5。
- 本轮不评估 `reed` 与 sgkit 的**功能取舍**对错（例如"要不要做 GWAS"），只记录事实。
- 不重开 v17 的 distance 结论（对拍矩阵），也不重开 v27 的多等位拍板（"项目不做多等位"）。
- 不裁决产品问题（是否公开发布、是否公开 tests/docs、是否保留三个 0.1.0 crate）。
- 第二部分只动 API 形状，不改任何统计公式/输出格式（v30 第一部分已审）。
- 第三部分的"建议归属"是选项，不是决定；其中 **R12（`reed-filter` 去掉 `reed-io` 依赖）需单独拍板**。
- 第三部分**未装工具**：`rust-analyzer`（仅有 rustup 代理）、`cargo-modules`、`cargo-public-api`、
  `cargo-semver-checks`、`cargo-check-external-types`、`cargo-hack`、`cargo-fuzz`、`cargo-deny`
  均未安装；本轮用 `rg` + `codegraph` + `cargo metadata` 替代，**公共 API 快照与 semver 基线仍未生成**。
  共变分析为 `git log --no-merges --name-only` 手工统计，未排除批量改名/格式化（可能高估）。
- 第四部分**工具链要求**：`cargo-check-external-types 0.5.0` 需固定 `nightly-2026-03-20`（默认 nightly
  的 rustdoc JSON v61 不兼容）；`cargo public-api` 不展开整模块再导出，故"可达面"是下界；
  `rust-analyzer diagnostics` 只有 WeakWarning，不能替代 clippy/rustc 门。
- 第四部分的子代理实测项（external-types 9 份表、未用依赖、CI 条件、RA 620 条诊断等）见
  [STATUS.md](STATUS.md) §7.4；本轮同样**未改代码或 CI**，未生成 `ci/api/` 快照。
