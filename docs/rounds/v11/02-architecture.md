# 02 修订规范（N1…N13）

> 每条 N 对应 01-review.md §1 的一条或多条缺陷（A11-x），给出**可执行的规范**
> 与**验收标准**。规范只约束"公开面该长成什么样"，不重开 v5/08 §6 的分工裁决。
> 标记：【规】必须；【荐】推荐；【决】需拍板。

---

## §1 请求默认值：`Levels` 不派生 `Default`（N1 ← A11-1）【规】

`Levels` 的"默认"是**语义**，不是 `bool` 的零值。删除
`#[derive(...)]` 里的 `Default`，手写：

```rust
impl Default for Levels {
    fn default() -> Self {
        Self { window: true, contig: false, genome: false }
    }
}
```

**验收**：新增测试 `default_request_emits_window_rows`——只补
`kinds`/`populations`/`grid` 的 `StatsRequest::default()` 运行后
`report.windows` 非空；把 `Default` 改回派生即红。

---

## §2 驱动签名对齐：`run_stats` 按引用取请求（N2 ← A11-2）【规】

`run_stats`／`run_stats_opened` 的 `request` 参数由 `StatsRequest` 改
`&StatsRequest`，与 `run_filter`（`&FilterSpec`）／`run_reports`（`&ReportRequest`）
一致。内部若需持有，`kinds` 等按需 `clone`（已用 `Arc` 的位置不变）。

**验收**：README「Library use」示例（含 `&request.kinds`）作为 doctest 或
下游编译门通过；`reed-cli` 调用点同步改为 `&request`。

---

## §3 渲染契约：`to_args` 可失败（N3 ← A11-3）【规】

`StatsRequest::to_args` 改为：

```rust
pub fn to_args(&self, group_names: &[String]) -> Result<Vec<String>, crate::error::Error>
```

把现在的 `assert!`/`panic!`（`render.rs:63/67/71/93`、混合 Fst/Pbs、
`Regions` grid、`LdRequest` 非默认）改成 `Err(Error::usage(...))`，
文案沿用现有 panic message。`#[must_use]` 保留。

同时回改 README:176 的 "`StatsRequest::to_args()`" 为
"`StatsRequest::to_args(&group_names)`"。

**验收**：`to_args(default_kinds(n), &names)` 返回 `Err`（而非 panic）；
正常窗口态仍返回原 flag 向量（现有往返测试扩到 `Result`）。

---

## §4 错误面收敛（N4 ← A11-4）【规】

1. 门面自身构造器返 `reed::Error`：`Levels::parse`、`assign_groups`
   （后者也可留在内部，若公开则返 `Error`）。
2. `String` 型 parse 改可并入 `reed::Error`：`AfKind::parse`、`GtClassRule::parse`
   返 `FilterError`（`Spec(String)`）或 `reed::Error`（二者择一，避免再添第七种）。
3. 下层 crate 的公开构造器（`MissingRate`/`Region`/`WindowGrid`）**可保持**原错误类型，
   只要 01-review §2 的核减成立；若要保持"单一错误类型"的对外表述，则在门面加薄包装。

**验收**：探针里六种来源都能 `.into::<reed::Error>()`（编译通过）；
`reed-cli` 去掉 `Levels::parse(...).map_err(|e| ...to_string())` 的字符串桥接。

---

## §5 `Error` 变体封口（N5 ← A11-5）【规】（已定：方案②b，用户 2026-10-14）

**裁决**：选②b——`Error` 由公开 enum 改为 struct 包装，变体对外不可见，
对外只留"分类 + 消息"。这正合 v5/08 §6.2 定的错误契约（"单一 `reed::Error` …
`kind() -> ErrorKind::{Usage, Runtime}`"），且与 cli.md §9 的两类出口码一致。

**前提（已核实，复审轮记录在案）**：`reed` 之外没有任何代码匹配 `Error` 的
变体（对 `match`/`if let`/`matches!` 全仓 grep 为空）；CLI 的错误出口只有
`eprintln!("reed: {err}")`（Display）与 `match err.kind()`（`main.rs:18`）。
因此封口对本仓**零破坏**。②a（只加 `#[non_exhaustive]`）两头不靠，**不采用**。

### 5.1 目标形态

`crates/reed/src/error.rs` 重写为：

```rust
use reed_core::error::CoreError;
use reed_engine::error::EngineError;
use reed_filter::error::FilterError;
use reed_io::error::IoError;

/// The single error type across every `run_*` entry point. Its variants are
/// sealed: callers classify with [`Error::kind`] and read with `Display`/`source`.
#[derive(Debug)]
pub struct Error(Box<ErrorInner>);

#[derive(Debug)]
enum ErrorInner {
    Core(CoreError),
    Io(IoError),
    Engine(EngineError),
    Spec(FilterError),
    Usage(String),
}

/// The exit-code class of an [`Error`] (cli.md §9: 2 = usage, 1 = runtime).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ErrorKind {
    /// Exit 2: the command line or input declaration is at fault.
    Usage,
    /// Exit 1: the run failed for reasons the command line cannot fix.
    Runtime,
}

impl Error {
    /// Wraps a message as a usage failure.
    #[must_use]
    pub fn usage(msg: impl Into<String>) -> Self {
        Self(Box::new(ErrorInner::Usage(msg.into())))
    }

    /// The exit-code class (cli.md §9). Request-construction failures are usage;
    /// transport, decode and engine failures are runtime.
    #[must_use]
    pub fn kind(&self) -> ErrorKind {
        match &*self.0 {
            ErrorInner::Usage(_)
            | ErrorInner::Core(_)
            | ErrorInner::Spec(FilterError::Spec(_) | FilterError::Format(_)) => ErrorKind::Usage,
            ErrorInner::Io(_)
            | ErrorInner::Spec(FilterError::Io(_) | FilterError::Decode { .. })
            | ErrorInner::Engine(_) => ErrorKind::Runtime,
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &*self.0 {
            ErrorInner::Core(e) => write!(f, "{e}"),
            ErrorInner::Io(e) => write!(f, "{e}"),
            ErrorInner::Engine(e) => write!(f, "{e}"),
            ErrorInner::Spec(e) => write!(f, "{e}"),
            ErrorInner::Usage(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &*self.0 {
            ErrorInner::Core(e) => Some(e),
            ErrorInner::Io(e) => Some(e),
            ErrorInner::Engine(e) => Some(e),
            ErrorInner::Spec(e) => Some(e),
            ErrorInner::Usage(_) => None,
        }
    }
}

impl From<CoreError> for Error {
    fn from(e: CoreError) -> Self { Self(Box::new(ErrorInner::Core(e))) }
}
impl From<IoError> for Error {
    fn from(e: IoError) -> Self { Self(Box::new(ErrorInner::Io(e))) }
}
impl From<EngineError> for Error {
    fn from(e: EngineError) -> Self { Self(Box::new(ErrorInner::Engine(e))) }
}
impl From<FilterError> for Error {
    fn from(e: FilterError) -> Self { Self(Box::new(ErrorInner::Spec(e))) }
}
/// Raw I/O (`?` in the CLI translation layer) enters the same `Io` class as
/// the typed transport errors.
impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self { Self(Box::new(ErrorInner::Io(IoError::from(e)))) }
}
```

对外**保留不变**：`usage()`、`kind()`、`ErrorKind`、五个 `From` impl、
`Display`、`source()`。对外**去掉**：五个公开变体（外部不能再 `match`
`Error::Io(..)` 这类写法，也不能再命名 payload 类型——这正是封口的目的）。

### 5.2 迁移点（只有两处）

封口的改动面极小，因为 reed 之外无人匹配变体：

1. **`error.rs` 本体**：上面 5.1 的全部。其 `#[cfg(test)]` 里的三处直接构造
   变体（`Error::Core(CoreError::UnknownSample{..})`、`Error::Engine(EngineError::Aborted)`、
   `Error::Spec(FilterError::Decode{..})`）改走 `From`：
   `Error::from(EngineError::Aborted)` 等。`kind()` 的断言一字不改。
2. **`crates/reed/src/index.rs:624`**：
   `crate::error::Error::Engine(EngineError::Job("BCF input ends inside a record (truncated file)".into()))`
   → `Error::from(EngineError::Job("…".into()))`。

### 5.3 验收

- `cargo doc --no-deps -p reed` 中 `Error` 显示为 struct（无变体页）；
  `use reed::IoError` 这类探针**继续失败**——这正是封口的预期，不再是缺陷。
- 下游 crate 写 `match err { reed::Error::Io(..) => .. }` 编译失败（变体不可见）。
- CLI 行为不变：双 profile 下 `cargo test --workspace` 全绿；`reed-cli` 的
  退出码（Usage=2/Runtime=1）与错误文案逐字节不变（`kind()`/`Display` 单测与
  e2e 不变）。
- `kind()` 的 6 个既有断言（error.rs 测试）原样通过。

### 5.4 后续若需更细分类

封口后若要"按错误类别结构化分支"（如 Io 才重试），**不要再把内部类型放回公开面**：
加命名谓词即可，如 `pub fn is_io(&self) -> bool` / `pub fn is_usage(&self) -> bool`。
这条留在公开面上的成本远低于暴露四个内部类型。


---

## §6 报告类型与驱动放置在根（N6 ← A11-6）【规】

`lib.rs` 根补充：

```rust
pub use crate::stats::writers::reports::{ReportsOutcome, run_reports, run_reports_opened};
pub use reed_filter::funnel::FilterReport;
```

（`run_reports` 现状只在 `reed::reports::*`；`FilterReport` 现状只在
`reed::funnel::*`。）README:172-175 的驱动清单补 `run_reports`。

**验收**：探针 `use reed::{run_reports, FilterReport, ReportsOutcome};` 通过；
`reed::reports::*` 旧路径保留（不破坏）。

---

## §7 报告谓词（N7 ← A11-7）【规/荐】

`ReportRequest::is_empty` 并入 `summary`：

```rust
!(self.freq || self.hardy || self.het || self.missing || self.depth || self.summary)
```

若拍板"该谓词无消费者，删"，则删除并移除 `#[must_use]`。

**验收**：`ReportRequest { summary: true, ..Default::default() }.is_empty() == false`。

---

## §8 删除死再导出（N8 ← A11-8）【规】

删 `lib.rs:78` 的 `pub use reed_stats::ld::LdWindow;` 与 `:98` 中
`reed_io::{bgzf, header}` 的 `header` 一项。**`bgzf` 保留**——首轮"0 引用"
的统计漏了 `crates/*/tests/`（`.gitignore` 所致，见 01-review §1 A11-8 的更正），
`reed-cli` 的两个集成测试用它写 BGZF 夹具。同步删 `lib.rs:98` 注释里
"the BGZF helpers the tests write fixtures with" 的措辞（改为 front end 与其
集成测试实际使用的说法）。

**验收**：`grep -rn --no-ignore-files 'reed::\(header\|LdWindow\)' crates` 为空；
`cargo check --workspace --all-targets` 通过。

---

## §9 `OpenedInput` 收口（N9 ← A11-9）【荐/决】

- **荐**：字段私有化，提供受控构造（内部用 `pub(crate)` + `#[doc(hidden)]` 入口），
  `run_*_opened` 的签名不变（仍收 `OpenedInput`）。`reed-cli` 只用 `open*` 的返回值，
  不受影响。
- **决**：若认为"传输词汇"应作为 §6.2 的批准例外，则维持现状，但在 §6.2 补一行
  "`OpenedInput` 及其传输类型是批准例外"，并保证这些类型可在门面命名。

**验收**：字段私有后 `cargo build --workspace --all-targets` 通过；
`e2e_*` 需要的字段经方法访问。

---

## §10 测试钩子不占公开面（N10 ← A11-10）【荐】

`FilterOptions.strategy` 与 `ExecStrategy` 加 `#[doc(hidden)]`；
进一步可 `#[cfg(any(test, feature = "test-hooks"))]`。CLI 的策略注入在
测试配置下保留。

**验收**：`cargo doc --no-deps -p reed` 的公开项不含 `ExecStrategy`；
`reed-cli` 测试仍通过。

---

## §11 内存态群体映射（N11 ← A11-11）【规/决】

`Populations` 增变体（名字待定，如 `Assignment`）：

```rust
pub enum Populations {
    All,
    File(PathBuf),
    Assignment(Vec<(String, String)>), // (sample, group)
}
```

`run_stats_opened` 对 `Assignment` 走与 `assign_groups` 相同的行规则
（`GroupTable` 复用一个 `from_rows`），组序仍取首次出现顺序。

**验收**：新增测试——同一映射分别以 `File` 与 `Assignment` 传入，
`StatsReport::windows` 逐字节一致。

---

## §12 `WindowRow.end0` 命名（N12 ← A11-12）【荐】

改名 `end0` → `end`（`request.rs:232`；更新 `windows.rs` 的 `row.end0` 引用），
文档改为 "1-based inclusive end (= 0-based exclusive end)"。

**验收**：`grep -rn '\.end0' crates` 为空；writers 测试全绿。

---

## §13 合并重复 helper（N13 ← A11-13）【荐】

`fn flag(args, name, value)`（`stats/render.rs:15`、`spec/render.rs:19`）合并到
一处共享实现（放 `reed-core` 的小工具模块或门面内部模块），两处调用改引用。
`fn f(value: f64)` 保留单份。

**验收**：`grep -rn 'fn flag' crates` 命中 1 处；两处渲染测试全绿。

---

## §14 spec / README 回改清单

| 文档 | 回改 |
|---|---|
| README「Library use」 | ① `run_stats` 示例改为 `&request`（或先取 `kinds`）；② `to_args` 句补 `&group_names`；③ 驱动清单补 `run_reports` |
| `docs/rounds/v5/08-api-surface.md` §6.2 | 补批准例外说明（`_opened`／`OpenedInput`），或在 v11 落定后于其上加"由 doc_v11 修订"注记 |
| `ci/check-doc-drift.sh` | 登记 01-review §3 的 7 条守卫 |
