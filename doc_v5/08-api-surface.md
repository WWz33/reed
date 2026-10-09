# 08 项目定位：现在像不像一个 API 库

> 用户问（2026-10-09 晚）："我们现在写的像是个 API 库么"。本文用代码事实回答，并给出两种自洽的终态供裁决。结论先行：**现在是一个 CLI 产品，库是它的内部结构；拿出去不像 API 库。** 这不是缺陷，是尚未做出的定位选择；但有几处"假装是库"的残留要清掉。

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
| 产品 | `reed` 二进制；crate 是内部结构 | 一个可 `cargo add` 的 `reed` crate（+ 可选 Python 绑定，v1 D12） |
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
- Python 绑定（PyO3）作为第一个外部消费者，倒逼公开面设计。
