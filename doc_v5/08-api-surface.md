# 08 项目定位：现在像不像一个 API 库

> 用户问（2026-10-09 晚）："我们现在写的像是个 API 库么"。本文用代码事实回答，并给出两种自洽的终态供裁决。结论先行：**现在是一个 CLI 产品，库是它的内部结构；拿出去不像 API 库。** 这不是缺陷，是尚未做出的定位选择；但有几处"假装是库"的残留要清掉。
>
> **用户随后定调（同日）："类似于前后端分类的形态。"** 即 A 的加强版：后端有一条明确、可测的 API 边界，但不对外发布；前端只是这条边界的文本皮肤。细则与现状违规见 §6。

---

## 1. 事实

| 维度 | 现状（已核实） | 像库吗 |
|---|---|---|
| 消费者 | `reed` 门面库只有两个使用者：`reed-cli` 与 `crates/reed/tests/*`；工作区之外没有任何调用方 | 否 |
| 公开入口 | `run_filter`、`run_stats[_opened]`、`run_reports[_opened]`、`run_index`、`run_qc`、`run_qc_report`：吃路径 + 选项结构，返回报告。这是"把命令行搬进函数"，不是面向嵌入的 API | 半 |
| 公开面的收口 | `reed/src/lib.rs` 把 `bcf_input`、`filtered`、`input` 这些实现模块直接 `pub mod`，还 `pub use stats::writers::{pixy, qc, qc_report, reports}`（注释自称 "historical public paths"）；`reed-io` 把 `noodles::vcf` 整个 re-export，noodles 的版本因此成了公开 API 的一部分；reed-stats 30 个顶层 pub、reed-core 24 个，没有 curated 的 `pub use` 清单 | 否 |
| 文档与例子 | 7 个 `lib.rs` **零 doctest**；`examples/` 里的 `gt_scan`、`qc_summary`、`bcf_dump` 是 M2/M3 的验收驱动，不是给用户看的 | 否 |
| README | "Library (reed-core)" 示例用的是 `ReadOptions/WindowBuffer/estimators`，"Library (reed-filter)" 用 `exec/emitter/FilterSpec::builder`——这些 API 在 doc_v2/v3 已删除（现在只存在于 oracle）；CLI 示例写 `--groups --size --out out.tsv`、`reed qc`、"legacy `reed-filter` binary retained until 0.5.0"。整份 README 停在 v0.2 | 否，且有害 |
| 发布元数据 | 7 个 crate 可发布（只有 oracle `publish = false`），但没有 `repository/readme/keywords/categories`、没有 `[package.metadata.docs.rs]`、没有 CHANGELOG；`reed-core` 的 description 还是 v0.2 的"windowed genotype buffers"；`naming-conventions.md` 说 "reed-core 0.2.0 和 reed-io 0.2.0 是 published surfaces"，与现实不符 | 否 |
| 类型与错误 | `thiserror` + `#[non_exhaustive]`（4 处）+ `Locus` 定位；`FilterSpec`/`StatsRequest` 是类型化结构，`to_args()` 往返有测试；纯 crate（stats/filter 无 I/O 无线程）、L0–L5 单向分层、`#[doc(hidden)]` 纪律（6 处） | **是**：库的骨架在 |

一句话：**骨架是库的，皮肤是工具的，门脸是过期的。**

---

## 2. 为什么会这样（不是错误）

用户的目标自始至终是"多核把 stats 和过滤做快"，对标 bcftools / vcftools / plink2 ——这三者都是**工具**。把实现拆成 8 个 crate 是为了架构（纯度、可测、并行层与策略层分离、oracle 差分），不是为了给别人 `cargo add`。doc_v2 02 §7 写的"库 API 仍是类型化结构，CLI 是它的文本表示"，是内部纪律（往返性、可测性），不是对外承诺。所以今天的形态与目标一致；不一致的只有那些从 v0.2 时代留下的"库姿态"（README 示例、noodles re-export、historical public paths、可发布的 crate）。

---

## 3. 两种自洽的终态

| | **A. 工具优先（现状的诚实版）** | **B. 真正的 API 库** |
|---|---|---|
| 产品 | `reed` 二进制；crate 是内部结构 | 一个可 `cargo add` 的 `reed` crate（Python 绑定已裁决不做） |
| 公开面 | 门面只保留一组**嵌入 API**：`FilterSpec`/`StatsRequest` + `run_filter/run_stats/run_index` + 报告与错误类型；其余 `pub(crate)` 或 `#[doc(hidden)]` | 同左，但要 curated `pub use` 清单、每个公开类型有文档与 doctest、`examples/` 面向用户、不 re-export noodles |
| 发布 | 7 个内部 crate `publish = false`（oracle 已是）；`reed` 也不发布，或只打 git tag | `reed` 发布 crates.io；子 crate 要么 `publish = false` 要么作为 `reed-*` 内部包发布且声明无稳定性 |
| 稳定性 | 只对 **CLI 参数与输出格式**承诺稳定（这才是 bcftools/plink2 式的契约），版本号跟二进制 | semver：公开类型/函数变更 = minor/major；CHANGELOG；MSRV 政策；docs.rs 构建 |
| 文档 | README = 工具文档：安装、三条命令、与参照物的对拍结果、性能；删除所有库用法段 | README + rustdoc + examples + 迁移指南 |
| 成本 | 一个卫生 PR（§4） | 约一个里程碑，且应在 X4–X6（统计出口、多等位）稳定**之后**做，否则公开面会跟着改 |

**推荐 A**，理由三条：与目标一致；现在做 B 会把正在变的东西（`--fst-global`、SFS、f 统计、LD）锁成承诺；A 做完之后若想转 B，是打包工作不是重构（骨架已经是库的）。

---

## 4. A 的清单（并入 07 X1）

| # | 动作 | 验证 |
|---|---|---|
| L-1 | README 重写为 v0.3 现实：删 "Library (reed-core/reed-filter)"、`reed qc`、legacy binary、`--groups --size` 示例；写三条命令（引用 02-cli.md）、对拍结果、性能表、仓库结构（06） | 人工；`! grep -q 'WindowBuffer\|exec::execute\|--groups' README.md` 进守卫 |
| L-2 | 门面公开面收口：`bcf_input`、`filtered`、`input` → `pub(crate)`（测试经 `#[doc(hidden)] pub mod internal` 进入）；删 `pub use stats::writers::{pixy, qc, qc_report, reports}`（随 C-1/C-3）；`lib.rs` 顶部一段 crate 文档写明"嵌入 API 只有 `FilterSpec`/`StatsRequest`/`run_*`/报告/错误" | `cargo doc --no-deps -p reed` 的公开项与本行一致；`reed-cli` 与 e2e 测试仍编译 |
| L-3 | `reed-io` 不再 `pub use noodles::vcf as noodles_vcf`（下游 crate 在各自 Cargo.toml 里显式依赖 noodles） | 构建通过；`grep -rn noodles_vcf crates` 为空 |
| L-4 | 7 个内部 crate `publish = false`；`reed` 暂也 `publish = false`；crate description 更新（`reed-core` 不再是"windowed genotype buffers"，`reed` 不再提 `run_qc`） | `grep -c 'publish = false' crates/*/Cargo.toml` = 8 |
| L-5 | `naming-conventions.md` 头部"published surfaces"一句改为"内部 API，仅对 CLI 与输出格式承诺稳定"；稳定性承诺写进 02-cli.md §0 | 人工 |
| L-6 | `examples/` 三个验收驱动移到 `bench/suites/`（它们本来就是基准）或删除 | `cargo build --examples` 无目标 |
| L-7 | 对 `FilterSpec`/`StatsRequest`/`run_*` 各写一个 doctest（就是 02 §7.5 的例子用 API 写一遍），作为"嵌入 API 可用"的最低证明 | `cargo test --doc -p reed` ≥ 3 |

---

## 5. 若裁决 B，额外要做的

- 公开面冻结在 X4–X6 之后；`reed/src/lib.rs` 变成纯 `pub use` 清单，禁止 `pub mod` 实现模块；
- 每个公开类型 `#[non_exhaustive]`（已有 4 处，补齐 `FilterReport`/`StatsOutcome` 等）；错误类型保留 `source()` 链；
- `[package.metadata.docs.rs]`、`repository`、`readme`、`categories = ["science::bioinformatics"]`；CHANGELOG（Keep a Changelog）；`cargo semver-checks` 进 CI；
- API 级集成测试与 CLI 测试分开（`crates/reed/tests/api_*.rs` 只用公开项）；
- 外部消费者：Python 绑定已裁决不做；B 若启动，由 API 级集成测试与 doctest 充当第一个消费者。

---

## 6. 裁决：前后端形态【规】

### 6.1 分工

```
前端   reed-cli（唯一的前端；2026-10-09 裁决：不做 reed-py，也不做其他 Python 绑定）
       只做：参数解析与校验文案、单位、--help、verbosity 与 stderr 文案、退出码、--timing 落盘、provenance 的 argv、进度
契约   02-cli.md：参数语法（后端用类型表达 = FilterSpec / StatsRequest / IndexRequest，前端用文本表达）、
       输出文件格式、报告 JSON、错误分类（Usage / Runtime）
后端   reed（门面）→ reed-filter / reed-stats / reed-engine / reed-io / reed-core
       只做：计算、文件读写、结构化结果、typed 错误、tracing 事件；不认识 clap，不打 stderr，不定退出码
```

**依赖规则：前端只依赖门面 `reed`**（外加 clap、serde_json 等外部 crate）。门面 re-export 前端需要的一切。前端的 `Cargo.toml` 里出现任何 `reed-core / reed-io / reed-engine / reed-stats / reed-filter` 即违规。

### 6.2 后端公开面（`reed/src/lib.rs` 只有 `pub use` 清单，禁止 `pub mod` 实现模块）

| 类别 | 项 |
|---|---|
| 请求 | `FilterSpec`（含 `SiteRules`、`TailRules`、`MaskRules`、`RecomputeSpec`、`SampleSelection`、`Region`、`MissingRate`）、`StatsRequest`（含 `StatKind`、`Levels`、`WindowGrid`、`Populations`）、`IndexRequest`；表达式入口 `expr::parse` |
| 运行 | `run_filter(&FilterSpec, Input, Output, &ThreadBudget) -> Result<FilterReport>`、`run_stats(&StatsRequest, Input, &ThreadBudget) -> Result<StatsOutcome>`、`run_index(&IndexRequest, &Path, &ThreadBudget) -> Result<IndexReport>`；`plan::compile`、`assign_groups` 这类"把请求变成可执行计划"的步骤在 `run_*` **内部**，前端不碰 |
| 结果 | `FilterReport`（漏斗、`sample_rule`、`StageMetrics`）、`StatsOutcome`（窗口表、报告表、日志条目；**结果是数据，写文件是另一步**）、`IndexReport` |
| 写出 | `reed::writers::{pixy, plink2, vcftools, bcftools_stats, dadi, fsc2}`：把 `StatsOutcome` 写成 02 §7.4 的文件；CLI 调用它，Python 前端可以不调 |
| 错误 | 单一 `reed::Error`（包 Core / Io / Engine / Spec，保留 `source()` 链），`kind() -> ErrorKind::{Usage, Runtime}`；前端据此映射退出码 2 / 1，不再各自判断 |
| 预算 | `ThreadBudget`、`MemBudget`（从 reed-engine re-export） |
| 文本契约 | `FilterSpec::to_args()` / `StatsRequest::to_args()` 与 `parse` **留在后端**：provenance（`##reed_filterCommand`）必须与前端无关，`reed-py` 跑出的 header 要与 CLI 逐字节相同 |

### 6.3 现状违规（2026-10-09 核实）

| 违规 | 位置 | 处理 |
|---|---|---|
| 前端直接依赖下层 crate | `reed-cli`：`use reed_engine::budget::{MemBudget, ThreadBudget}`（1 处）；`reed_filter::{spec::*, expr::parse, plan::compile, funnel::FilterReport}`（19 处，含测试） | 门面 re-export；`reed-cli/Cargo.toml` 只留 `reed`（并入 03 C-8） |
| 后端直接打 stderr（5 处） | `reed/src/input.rs:268`；`reed/src/filter/mod.rs:672`（`--force-samples` 警告）；`reed-io/src/chunk/bgzf.rs:240,390`（缺 EOF 标记）；`reed-filter/src/plan/mod.rs:146`（陈旧 INFO 警告） | 改 `tracing::warn!`；前端的订阅器（W1 已有）按 verbosity 决定是否打印、打印成什么文案 |
| 前端里有后端逻辑 | `reed-cli/src/translate/filter.rs:310` 直接调 `reed_filter::plan::compile`；`translate/stats.rs` 自己 `assign_groups` 并拼 `StatKind` 列表 | 移入 `run_filter` / `run_stats`；前端只构造请求 |
| 门面泄漏实现 | `pub mod bcf_input / filtered / input`；writers 的 historical re-export | §4 L-2 |
| 前端渲染报告 | `reed-cli` 的 `eprint_report`（`translate/filter.rs:379`） | **正确**，渲染属前端，保持 |

### 6.4 守卫（进 `scripts/check-doc-drift.sh`，搬迁后 `ci/`）

1. `! grep -qE 'reed-(core|io|engine|stats|filter)' crates/reed-cli/Cargo.toml`
2. `! grep -rn 'eprintln!\|println!' crates/{reed,reed-core,reed-io,reed-engine,reed-stats,reed-filter}/src --include=*.rs` 为空
3. `grep -c '^pub mod' crates/reed/src/lib.rs` = 0
4. `cli_contract`：每个 clap 结构体都能从 `reed::spec` 的类型往返（现有往返测试扩到 stats 与 index）

### 6.5 第二个前端：不做（裁决 2026-10-09）

用户裁决**不做 `reed-py`**，也不做其他 Python 绑定。§6.1–6.4 的边界规则不依赖第二个前端而成立。"结果是数据、写文件是另一步"（§6.2 的 `StatsOutcome` 与 `reed::writers` 分离）仍然保留：理由是它让写出层可以单独测试、输出格式变动不触及计算，与是否有第二前端无关。v1 架构文档的 D12（"Python 绑定 M8 之后再做"）随之作废。
