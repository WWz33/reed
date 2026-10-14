# doc_v11 实施状态

> 本页记录 `01-review.md` §1 缺陷清单与 `06-roadmap.md` R0–R3 的实施状态。
> 本轮为**落地轮**：复审（A11-1…A11-13）与修复规范（N1…N13）已在同一分支
> 实施完毕，每条附证据（测试名 / 守卫 / 变更点）。复审基线（修复前的
> 实测/编译证据摘录）见文末。

## 实施状态（R0–R3）

| 项 | 缺陷 | 严重度 | 状态 | 证据 |
|---|---|---|---|---|
| R1-a | A11-1 默认 `Levels` 产出 0 行 | **P0** | **已修** | 手写 `impl Default for Levels`（`request.rs`，window=true，删派生）；测试 `default_request_emits_window_rows`（api_surface.rs）；守卫 "the default Levels is the window level" |
| R1-b | A11-2 `run_stats` 按值取请求 | P1 | **已修** | `run_stats`/`run_stats_opened` 请求改 `&StatsRequest`（`driver.rs`）；41 处测试调用点、CLI 2 处同步；测试 `readme_library_example_compiles_and_runs`（借用后仍可用）；lib.rs 新增 doctest；守卫 "run_stats borrows its request" |
| R1-c | A11-3 `to_args` 对 `default_kinds` panic | P1 | **已修** | `to_args(&self, group_names) -> Result<Vec<String>, Error>`（`render.rs`），5 处 panic/assert 改 `Err(usage)`；测试 `to_args_errors_instead_of_panicking_on_unrenderable_kinds`；守卫 "StatsRequest::to_args errors instead of panicking"；README 句改 `to_args(&group_names)` |
| R1-d | A11-4 门面六种错误类型 | P1 | **已修** | `Levels::parse`/`assign_groups`/`assign_groups_in_memory` 返 `reed::Error`；`AfKind::parse`/`GtClassRule::parse` 由 `String` 改 `FilterError`；CLI 的 `to_string()` 桥接删除；测试 `facade_errors_converge_on_reed_error` |
| R1-e | A11-5 `Error` payload 不可命名 | P1 | **已修（②b，用户拍板）** | `pub struct Error(Box<ErrorInner>)` + 私有 `enum ErrorInner`（`error.rs`）；`index.rs:624` 改走 `From`；测试 `facade_errors_converge_on_reed_error` + `source_chain_reaches_the_payload`；守卫 "the facade Error is sealed"；下游探针确认 `match reed::Error::Usage(_)` 编译失败（预期） |
| R1-f | A11-6 `run_reports`/报告类型不在根 | P1 | **已修** | 根导出 `run_reports`/`run_reports_opened`/`ReportsOutcome`/`ReportRequest`/`FilterReport`（`lib.rs`）；测试 `run_reports_is_reachable_from_the_root` + 驱动签名检查；守卫 "run_reports and its types are at the facade root"；README 清单补 `run_reports` |
| R1-g | A11-7 `is_empty` 漏 `summary` | P2 | **已修** | `is_empty` 并入 `self.summary`（`reports.rs`）；测试 `summary_only_report_request_is_not_empty` |
| R2-a | A11-8 死再导出 | P2 | **已修（范围更正）** | 删 `LdWindow`、`reed_io::header`；**保留 `bgzf`**——首轮"0 引用"统计被 `.gitignore` 误导，`reed-cli/tests/{cli_contract,stdin_streaming}.rs` 有 6 处真实使用（见 01-review A11-8 的更正）；守卫 "no dead LdWindow/header re-exports" |
| R2-b | A11-9 `OpenedInput` 暴露传输内部 | P2 | **已修** | 三字段 `pub(crate)`，加 `facts()`/`bcf_dict()`/`source_mut()`（`source_mut` 为 `#[doc(hidden)]`，测试用）；CLI 与内部调用点改访问器；守卫 "OpenedInput fields are crate-internal" |
| R2-c | A11-10 测试专用 `ExecStrategy` 公开 | P3 | **已修** | `ExecStrategy`（`strategy.rs`）与 `FilterOptions.strategy`（`filter/mod.rs`）加 `#[doc(hidden)]` |
| R3-a | A11-11 缺内存态 `Populations` | P2 | **已修** | `Populations::Assignment(Vec<(String,String)>)` + `assign_groups_in_memory`（复用 `GroupTable::from_assignment`）；测试 `in_memory_assignment_matches_the_file_form`（与 File 形式逐字节一致）；README 新增示例 |
| R3-b | A11-12 `WindowRow.end0` 命名 | P3 | **已修** | `end0` → `end`（`request.rs`、`sink.rs`、`writers/windows.rs` 8 处、2 个测试）；文档写"1-based inclusive end (= 0-based exclusive end)" |
| R3-c | A11-13 `flag` helper 重复 | P3 | **已修** | 单一实现 `reed_core::args::push_flag`（新文件 `args.rs`），`stats/render.rs` 与 `spec/render.rs` 共用（后者沿用 `flag` 别名） |
| R3-d | 文档与守卫收口 | — | **已完成** | README「Library use」重写（借用、`run_reports`、`Assignment`、sealed error）；lib.rs 新增 A11-2 doctest；`ci/check-doc-drift.sh` 增 9 条 v11 守卫（55 条全绿） |

## R0–R3 进度

| 步骤 | 状态 | 说明 |
|---|---|---|
| R0 探针/红测试 | 已完成 | 新增 `crates/reed/tests/api_surface.rs`（7 测试，只 `use reed::…`），落地后全绿；复审期的临时复现见文末基线 |
| R1 P0+P1 行为契约（A11-1…7） | 已完成 | 见上表 |
| R2 公开面收口（A11-8/9/10） | 已完成 | 见上表（A11-8 范围按实测更正） |
| R3 缺失/命名/文档（A11-11/12/13 + 回改） | 已完成 | 见上表 |

## 验收门结果（2026-10-14）

| 命令 | 结果 |
|---|---|
| `cargo test --workspace`（debug） | 全绿（含 api_surface 7 项） |
| `cargo test --release --workspace` | 全绿 |
| `cargo clippy --workspace --all-targets -- -D warnings` | 0 警告 |
| `cargo test --doc -p reed` | 4 passed（新增 A11-2 doctest） |
| `cargo doc -p reed --no-deps` | 无警告 |
| `bash ci/check-doc-drift.sh` | 55 条全绿（含 9 条 v11 新增） |
| 下游探针（只依赖 `reed`） | README 示例编译并运行；`match reed::Error::Usage(_)` 按预期编译失败（封口生效） |

## 本轮未做（有意）

见 `06-roadmap.md` §4：不把下层 crate 的构造器错误全部统一；不删
`run_*_opened`；不删 `default_kinds`/`window_kinds`/`assign_groups`；
不删 `bgzf`/`chunk`/`metrics`/`budget`；不重排 v5/08 §6 的分工。

## 复审基线（2026-10-14，**修复前**）

复审 commit：`065175c`（doc_v10 全部落地），工作区干净。
以下为修复前证据摘录。证据通道：【实】实跑；【编】下游探针编译；【读】源码阅读。

| 缺陷 | 基线证据 | 通道 |
|---|---|---|
| A11-1 | 临时集成测试输出：`Levels { window: false, contig: false, genome: false }`；`windows.len() = 0, n_windows = 0` | 【实】 |
| A11-2 | 下游 crate 逐字编译 README:158-171：`error[E0382]: borrow of moved value: request` | 【编】 |
| A11-3 | 临时测试 `#[should_panic]`：`panicked at crates/reed/src/stats/render.rs:63: no --flag spelling for per-bp Pi` | 【实】 |
| A11-4 | 探针：`error[E0277]: the trait bound reed::Error: From<String> is not satisfied` | 【编】 |
| A11-5 | 探针：`error[E0425]: cannot find type EngineError / IoError in crate reed` | 【编】 |
| A11-6 | 探针：`E0432 unresolved import reed::run_reports / reed::FilterReport / reed::ReportsOutcome` | 【编】 |
| A11-7 | 探针程序运行：`summary-only is_empty() = true` | 【实】 |
| A11-8 | 首版：`grep -r` 报三者各 0 引用；**落地时复核更正**——`grep --no-ignore-files` 显示 `reed::bgzf` 有 6 处真实使用（CLI 集成测试），`header`/`LdWindow` 确为 0 | 【读】 |
| A11-9 | 探针：`fn touch(o:&reed::OpenedInput){ &o.facts; &o.bcf; &o.source }` 编译通过 → 字段公开 | 【编】 |
| A11-10 | 探针：`reed::plan::ExecStrategy::SinglePass` 解析成功 | 【编】 |
| A11-11 | 探针：`E0599: no variant ... named Assignment found for enum Populations` | 【编】 |
| A11-12 | `request.rs:231-232`：字段 `end0` vs 文档 "1-based inclusive end" | 【读】 |
| A11-13 | `stats/render.rs:15` 与 `spec/render.rs:19` 各一份 `fn flag` | 【读】 |

## 方法教训（记入本轮）

本仓 `.gitignore` 忽略 `crates/*/tests/`、`crates/reed/.gitignore` 忽略
`tests/`、根 `.gitignore` 忽略 `docs/`，而 ugrep 默认遵守 gitignore。因此
**任何"全仓引用统计"都必须带 `--no-ignore-files`**，否则会系统性漏掉
测试代码——A11-8 的首版结论（把 `bgzf` 判为死导出）就是这么错的，落地
编译时被 `cargo check --workspace` 的 6 个 `cannot find bgzf in reed` 抓出。
