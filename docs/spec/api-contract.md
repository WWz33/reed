# ADR: reed 公共 API 契约（P0）

> 状态：**接受并冻结（2026-11-12）**。依据 `docs/rounds/v30/03-synthesis.md` §4 P0。
> 发布包边界于 2026-11-19 按 [api-release-amendment.md](api-release-amendment.md) 修订；
> 门面稳定 API、MSRV、feature 与兼容承诺仍保持本 ADR 冻结内容。
> 本文件回答"发布什么、承诺什么、依赖什么、怎么验证"。变更必须另立 ADR，
> 不得在同一轮内悄悄改变；冻结期 6 个月（至 2027-05-12）。

## 1. 发布面

| 对象 | 决定 | 理由 |
|---|---|---|
| `reed`（门面，L4） | **发布**，crates.io | 唯一"设计出来的 API"；`run_*` + 请求/报告/错误 + `prelude` |
| `reed-core/io/engine/stats/filter/graph/distance/reduce` | **发布，0.5.0 lockstep** | `reed` 的正常依赖；必须先于门面按依赖顺序上架，作为实现支撑包，不扩大 `reed` facade 的稳定 API |
| `reed-cli` | `publish = false` | 二进制；CLI 输出另有 byte-identical 承诺 |
| `oracle` | `publish = false` | 测试/对拍工具 |

`reed` 的发布包由 Cargo 根据 crate 清单生成；测试和开发文档不是 facade 的发布契约。
发布前必须先按 `ci/check-package.sh` 输出的顺序上架全部正常 library 依赖，最后再上架
`reed` facade。`reed-cli` 与 `reed-oracle` 始终不发布。

## 2. 稳定面

- **承诺**：`reed` 根模块的公开项 + `reed::prelude`。这是 semver 的唯一对象。
- **`prelude` 的边界**：`reed::prelude` 只是根导出的便捷集合（`use reed::prelude::*;`），
  内容必须全部来自根导出；它不引入新的类型或额外承诺。
- **支持包**：8 个正常 library 依赖与 `reed` 使用 `0.5.x` lockstep 发布，满足 Cargo
  registry 依赖解析；它们的**公开 API 也是 lockstep 稳定面**（doc_v30 P1-7 第 2 选项）：
  支持包的删除/改签名必须按 semver bump，并由独立 API 快照与 semver 门验证。这样
  `reed` 公开签名中的支持包类型不会在未 bump 的情况下改变。
- **不发布**：`reed-cli` 与 `reed-oracle`；前者是二进制，后者是测试/对拍工具。
- **可命名性**：出现在 `reed` 公开签名里的每个类型，必须能从 `reed::*` 命名；
  不能命名的类型要么导出，要么改成 opaque（V30-Q02/T02）。
- **内部依赖**：`reed` 的公开签名不得直接使用 `reed-core` 等下层 crate 的类型名，
  除已在 `reed` 根显式 re-export 的项。

## 3. 外部依赖政策

默认 **opaque**：`noodles`、`serde`、`memmap2`、`rayon` 等第三方类型**不进入**
`reed` 的公开签名。确需暴露时只有两条路：

1. 在 `reed` 根显式 re-export，并登记到 `ci/external-types-allowlist.toml`；
2. 用 reed 自有类型包一层（推荐）。

现状（V30-T05：noodles 10 项、serde 37 个类型）属于待清偿债务，由 P3 的
`cargo-check-external-types` 门逐步收紧到 allowlist 之外为零。`serde` derive 产生的
`Serialize/Deserialize` 实现不视为签名暴露，但 `serde_json::Value` 之类的字段类型算。

## 4. MSRV / edition / feature

- **MSRV = 1.91**（`rust-version`，全 workspace 一致）；`edition = 2024`，`resolver = 3`。
- feature 只用于**可选后端**，不得改变默认语义：
  - `reed-reduce/faer-backend`（默认开）：纯 Rust 特征值后端；关闭后是零依赖回退，
    必须能编译（V30-Q31）。
  - 新增 feature 必须有 `cargo hack --feature-powerset` 计划与 `--no-default-features` 门。
- 同一输入下，feature 开关只允许影响性能，不允许影响输出字节。

## 5. 兼容承诺

| 面 | 承诺 | 门 |
|---|---|---|
| 库 API（`reed` 根 + `prelude`） | semver；删除/改签名按 major（0.x 按 minor） | `cargo semver-checks -p reed` |
| CLI stdout/stderr/输出文件 | 同一 0.x 系列内 **byte-identical**（除 `--version`/时间戳） | 黄金夹具 + 契约测试 |
| 公开文件格式（VCF/BCF/PHYLIP/012/…） | 只增不改；字段顺序与数值舍入视为契约 | 各命令的 parity/golden 测试 |
| 支持包 API（8 个 library crates） | 0.5.x lockstep；同样按 semver 管理 | 每个 crate 的 `cargo semver-checks` + `ci/api/<crate>.txt` 快照 |
| `docs/rounds/` | 历史记录，不是契约 | — |

## 6. 机器门

**本批已接入 CI（2026-11-12）**：

1. `cargo public-api -p <crate> -s > ci/api/<crate>.txt` 对 9 个发布 crate 入库，
   `ci/check-api.sh` 逐个做 diff；
2. `cargo check -p reed-reduce --no-default-features`（feature 空集必须能编译）；
3. dev/非 main 分支执行 fmt / release build / debug+release test / doctest /
   clippy `-D warnings` / `ci/check-doc-drift.sh` / `ci/check-structure.sh` 等完整门；
   main 和 PR→main 只验证公开树的 release build 与 MSRV source check，避免公开树
   缺少 dev-only docs/tests 时误报。具体分支边界以 `.github/workflows/ci.yml` 为准。

**P3 接入（2026-11-12）**：

4. **可命名性**：`ci/check-nameable.sh` 从快照提取下层 crate 的类型路径，
   要求每个类型都出现在 `crates/reed/tests/nameability.rs` 的 `reed::*` 导入里；
   编译失败即契约失败。`#[doc(hidden)] pub mod frontend` 是唯一的前端逃生门。
5. **外部类型**：`ci/check-external-types.sh`（固定 `nightly-2026-03-20`）+
   `ci/external-types-allowlist.toml`；只允许本 workspace 的模型 crate，
   `noodles`/`serde`/`rayon`/`memmap2` 进入签名即失败。
6. **feature**：`ci/check-features.sh`（`cargo hack --each-feature --workspace` +
   workspace 空 feature 集；包集合由 workspace 自动发现，新增 feature 不会漏门）。
7. **semver**：`ci/check-semver.sh` 对最新可达 tag 逐个跑 9 个发布 crate 的
   `cargo semver-checks -p <crate>`，**已设为 CI 必过**。workspace 已 bump 到
   `0.5.0`（0.x 的破坏版本），覆盖 P2 的 `RunContext` 签名变化、删除的
   `FilterOptions`/`PcaRequest`/`DapcRequest` 字段等；破坏性改动若不再 bump，
   该步会失败。
8. **文档派生**：`ci/check-docs.sh`（rustdoc `-D warnings`、每 crate doctest、
   `cargo modules orphans`、library.md 入口表对 `pub use … run_*` 的一致性）。
9. **package**：`reed` 及其全部正常 library 依赖均为 `0.5.0`、带 registry
   版本约束和可发布元数据；`reed-cli` 与 `reed-oracle` 保持 `publish = false`。
   `ci/check-package.sh` 验证 package 文件清单、正常 path 依赖的可发布性、版本约束
   和拓扑发布顺序。`cargo package -p reed --no-verify` 只有在这些前置 crate 已经发布
   到所选 registry 后才会通过；本地 gate 不把 registry 尚未包含它们误报为成功。

任何一条**已接入**的门在 CI 上失败，视为契约破坏，不得合并。

## 7. 功能边界（不是 semver 面）

以下能力**没有**完整实现，发布说明与文档不得写成已兼容：

- `.genome` 只保证 `Z0/Z1/Z2/PI_HAT/DST` 与 plink 1.9 对齐；`PPC/RATIO` 输出 `NA`，
  `RT/EZ/PHE` 是占位值（plink 的窗口计数与系谱推断不可从全位点统计复现）。
- `stats --summary` 不输出 `PSI` 与 `HWE` 两节；`--hardy`/`--het` 是替代路径。
- site-count (`WindowGrid::Sites`) 窗口仍返回 usage error，未实现。
- `--hwe P midp` 已实现 mid-p 调整（2026-11-23）；此前“可解析但拒绝”的边界解除。
