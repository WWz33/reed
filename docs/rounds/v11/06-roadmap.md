# 06 修复路线（R0–R3）

> 原则：**先把公开面的坑钉成可跑的探针/红测试，再修 P0+P1（行为契约），
> 然后做公开面收口（删/移/藏），最后收尾文档与守卫**。修复对象是
> 01-review.md §1 的 A11-1…A11-13；规范依据见 02-architecture.md 的 N1…N13。
> 不重开 v5/08 §6 的分工裁决、不重开 v6–v10 的结构项。每步满足 §1 通用门。

---

## 1. 通用验收门（每一步都必须满足）【规】

1. **双 profile 全绿**：`cargo test --workspace`（debug）与
   `cargo test --release --workspace` 都通过；
   `cargo clippy --release --workspace --all-targets` 不新增警告。
2. **门面无 `pub mod` 不破**：`grep -c '^pub mod' crates/reed/src/lib.rs` 仍为 0
   （v5/08 §6.4 守卫 3）。
3. **红 → 绿**：A11-1／A11-2／A11-3／A11-7／A11-11 先提交复现（探针或测试），
   再提交修复转绿；测试名描述缺陷本身。
4. **下游探针门**：新增 `crates/reed/tests/api_surface.rs`，**只 `use reed::…`**，
   断言 A11-2/6 的可达性与 A11-4 的错误收敛（编译通过即证据）。
5. **README 示例编译门**：README「Library use」示例纳入编译（doctest 或
   `include` 进 `api_surface.rs`）——这是 A11-2 将来不再复发的根因门。
6. **同一 PR 内删除被替代的旧路径**（如 §8 的死再导出、§13 的重复 helper），
   不留并行代码（v2 06 §1 纪律 7 沿用）。
7. 相应条目标记完成并注明 commit；`docs/rounds/v11/STATUS.md` 同步；
   每落一步在 `ci/check-doc-drift.sh` 登记 01-review §3 的对应守卫。

---

## 2. 步骤

### R0 探针与红测试先行（把公开面的坑钉住）— 门的基础

> R0 不修代码。建立一个**只依赖 `reed` 的**测试目标 `crates/reed/tests/api_surface.rs`
> 与一个临时运行测试，把 A11-1/2/3/7/11 钉成红（或编译失败），后续每步
> 只让自己那条转绿。

| 项 | 内容 | 门 |
|---|---|---|
| R0-a | A11-1：`StatsRequest::default()` + kinds/pops/grid 运行，断言 `windows` 非空。**修复前：空** | 测试 `default_request_emits_window_rows`；debug 与 release 都红 |
| R0-b | A11-2：把 README 库示例 `include!`/抄进 `api_surface.rs`。**修复前：E0382 编译失败** | `cargo test -p reed --test api_surface` 编译即红 |
| R0-c | A11-3：`to_args(default_kinds(2), &names)`。**修复前：panic** | 现 `#[should_panic]` 在 R0 记录；R1 改为断言 `is_err()` |
| R0-d | A11-7：`ReportRequest{summary:true,..default}.is_empty()`。**修复前：true** | 断言 `== false`，修复前红 |
| R0-e | A11-11：`Populations::Assignment(...)`。**修复前：E0599** | R0 记为编译失败；R1 后转绿 |

### R1 P0+P1 行为契约修复（A11-1/2/3/4/5/6/7）— 主量

> 都是"一处默认值／一个签名／一个返回类型"级别（02 §1–§7）。
> A11-5 已定方案②b（`Error` struct 包装），迁移点只有 `error.rs` + `index.rs:624`。

| 项 | 内容 | 规范 | 门 |
|---|---|---|---|
| R1-a | `Levels` 手写 `Default`（window=true）；R0-a 转绿 | 02 §1 | 通用门 3 |
| R1-b | `run_stats[_opened]` 请求改 `&StatsRequest`；同步 CLI；R0-b 转绿 | 02 §2 | 探针门；README 门 |
| R1-c | `to_args` 返 `Result`；R0-c 转绿；README 回改 | 02 §3 | 守卫 #3 |
| R1-d | 门面构造器错误收敛到 `reed::Error`（含两个 `String` parse） | 02 §4 | 探针门（六源可 `.into()`） |
| R1-e | `Error` 变体封口（②b）：enum → `Error(Box<ErrorInner>)`，变体不可见；改 `error.rs` + `index.rs:624` | 02 §5 | 探针门（变体不可匹配）；`kind()`/`Display` 单测不变 |
| R1-f | 根补 `run_reports`/`ReportsOutcome`/`FilterReport`；README 补列 | 02 §6 | 探针门；守卫 #5 |
| R1-g | `is_empty` 并入 `summary`；R0-d 转绿 | 02 §7 | 断言 |

### R2 公开面收口（删 / 藏 / 移）（A11-8/9/10）

| 项 | 内容 | 规范 | 门 |
|---|---|---|---|
| R2-a | 删死再导出 `LdWindow`、`reed_io::header`（`bgzf` 有真实消费者，保留） | 02 §8 | 守卫 #6；workspace 编译 |
| R2-b | `OpenedInput` 字段私有化 + 受控构造（或 §6.2 批准例外注记） | 02 §9 | workspace 编译；`e2e_*` 经方法访问 |
| R2-c | `ExecStrategy`/`strategy` `#[doc(hidden)]`（或 `test-hooks` feature） | 02 §10 | `cargo doc` 公开项检查 |

### R3 缺失 API + 命名 + 收尾（A11-11/12/13）+ 文档

| 项 | 内容 | 规范 | 门 |
|---|---|---|---|
| R3-a | `Populations::Assignment` + `GroupTable::from_rows`；R0-e 转绿 | 02 §11 | File/Assignment 逐字节一致 |
| R3-b | `WindowRow.end0` → `end`（更新 writers） | 02 §12 | `grep '.end0'` 为空 |
| R3-c | 合并 `flag`/`f` helper（保留单份） | 02 §13 | `grep 'fn flag'` 命中 1 |
| R3-d | README / v5/08 §6.2 回改 + drift 守卫登记 + STATUS 收口 | 02 §14 | `bash ci/check-doc-drift.sh` 全绿 |

---

## 3. 顺序与依赖

```
R0（探针/红） → R1（P0+P1 契约） → R2（收口） → R3（缺失/命名/文档）
```

- R1 内各项基本互不依赖；R1-b（改签名）应与 R1-c（README 回改）同 PR。
- R2 依赖 R1（先稳签名，再删/藏），也可并行——R2-a/c 与 R1 无交集。
- **若排期紧张**：R1-a（A11-1，P0）必须最先；R2/R3 可单独成轮，
  但 A11-2 的 README 编译门（通用门 5）建议随 R1 一起落，否则同类漂移会复发。

## 4. 本轮不做

| 项 | 原因 |
|---|---|
| 把下层 crate 的构造器错误全部统一（`MissingRate`/`Region`/`WindowGrid`） | 01-review §2 核减：它们在各自 crate 是公开类型，错误类型可保持；只收敛门面自身构造器 |
| 删 `run_*_opened` | v6 Y6.5 的 stdin 流式入口，有真实前端依赖；只处理其 `OpenedInput` 内部量 |
| 删 `default_kinds`/`window_kinds`/`assign_groups` | 有真实消费者；A11-3 修 `to_args` 后即可组合 |
| 删 `chunk`/`metrics`/`budget` 再导出 | 有真实消费者（01-review §2） |
| 重排 `docs/rounds/v5/08 §6` 的分工（前后端边界） | v5/08 裁决仍有效；本轮只把公开面收回它，不改分工 |
