# 01 门面 API 复审

> 复审对象：commit `065175c`（doc_v10 全部落地），工作区干净。
> 方法：门面公开项（`crates/reed/src/lib.rs` 的 41 项 `pub use`）通读；一个只依赖
> `reed` 的下游探针 crate 做编译期核对；一个临时集成测试做运行核对（跑完即删）。
> 基准契约：`docs/rounds/v5/08-api-surface.md` §6.2（后端公开面）、README「Library use」。
> 本轮不修改任何代码；修法以规范形式给出（02），实施顺序见 06。
> 证据通道：【实】实跑复现；【编】下游探针编译；【读】源码/文本（不可实跑）。

---

## 0. 总评

CLI（唯一前端）把门面的每一处粗糙都显式绕开了：`Levels::parse` 给默认值、
`run_stats_opened` 先开再传、`to_args` 只对窗口态调用、错误一律
`to_string()` 后再包 `Error::usage`。所以**在 CLI 这条路径上，本轮 13 条
缺陷没有一条会显形**。它们只在"换一个只依赖 `reed` 的调用方"时成片出现，
而 README 正是这样教人用库的——被教的那个例子，既不编译，也不产出结果。

分档见 README 的"结论速览"。以下逐条。

---

## 1. 缺陷清单

### A11-1 默认 `StatsRequest` 产出 0 行【实】【据】

- **严重度：P0**。后果形态：静默错误结果（空表，无报错）。
- **reed 现状**：
  - `Levels` 派生默认值：`#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]`
    （`request.rs:165`），三字段 `window/contig/genome` 全 `bool`，
    派生即全 `false`——而三处文档都反着写：`Levels` 的 `window` 字段注
    "One row per window (the default)"（`request.rs:167`），`StatsRequest.levels`
    注 "Output levels (default: window rows only)"（`request.rs:289`）。
  - sink 只在 `self.levels.window` 时出行：`emit_complete_windows`
    （`sink.rs:357`）、`close_contig` 的 Bp 分支（`sink.rs:405`）与 Regions
    分支（`sink.rs:427`）；`contig`/`genome` 行同理（`sink.rs:437`/`:190`）。
  - driver 把 `request.levels` 原样传给 sink，**无归一化**（`driver.rs:248`）。
- **实测**（临时集成测试，2 样本 / 100bp 窗口）：
  ```
  VERIFY levels default   = Levels { window: false, contig: false, genome: false }
  VERIFY windows.len()    = 0, n_windows = 0
  VERIFY render_tsv       = "contig\tstart\tend\tpi.pair.0\n"
  ```
- **正确先例（本仓库）**：CLI 从不依赖默认——`Levels::parse(&args.levels)`
  且 clap 默认 `--levels window`（`reed-cli/src/cli/mod.rs:338`）；
  测试也都显式给 `levels`（`e2e_stats.rs:82` 等）。问题只在**派生默认值**。
- **修法**：手写 `impl Default for Levels`，返回
  `Levels { window: true, contig: false, genome: false }`（删派生）。
  见 02 §1。
- **验收**：新增测试——`StatsRequest::default()` 只补 `kinds`/`populations`/`grid`
  后运行，断言 `!report.windows.is_empty()`（放回派生即红）。

### A11-2 `run_stats` 按值取请求，README 库示例不编译【编】

- **严重度：P1**。后果形态：对外承诺的用法编译失败。
- **reed 现状**：四个驱动的请求参数不统一——
  - `run_stats(path, request: StatsRequest, budget)` 按值（`driver.rs:28-35`）；
  - `run_filter(input, spec: &FilterSpec, opts: &FilterOptions)` 按引用（`filter/mod.rs:161`）；
  - `run_reports(input, request: &ReportRequest, budget)` 按引用（`reports.rs:733`）。
- **实测**（把 README「Library use」示例逐字拷进只依赖 `reed` 的下游 crate）：
  ```
  error[E0382]: borrow of moved value: `request`
    --> src/main.rs:11:58
    let report = reed::run_stats(Some("in.vcf.gz".as_ref()), request, ...)?;   // moved
    let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false); // borrowed after move
  ```
- **修法**：`run_stats`/`run_stats_opened` 改 `request: &StatsRequest`，
  与兄弟对齐；或至少把 README 示例改成先取 `kinds` 再 move。见 02 §2。
- **验收**：README 库示例进 doctest 门（`cargo test --doc`），编译通过。

### A11-3 `to_args` 对本仓 `default_kinds` panic【实】

- **严重度：P1**。后果形态：进程 panic（公开 helper 间不可组合）。
- **reed 现状**：`StatsRequest::to_args`（`render.rs:41`）对无 flag 拼写的 kind
  直接 `assert!`/`panic!`——`Pi/Dxy/ThetaW` 的 `per_bp=true`（`:63/:67/:71`）、
  `StatKind::Sfs`（`:93`）、混合 Fst/Pbs、`Regions` grid、`LdRequest` 的
  非默认参数。而公开 helper `default_kinds`（`request.rs:25`）**正好产出**
  `per_bp=true` 与 `Sfs` kind。
- **实测**：
  ```
  #[should_panic] to_args(default_kinds(2)):
  panicked at crates/reed/src/stats/render.rs:63: no --flag spelling for per-bp Pi
  ```
- **附带（文本契约漂移）**：README:176 写 "`StatsRequest::to_args()`"，实际
  签名需 `group_names: &[String]`（`render.rs:41`；名字型拼写 `--pbs A,B,C`
  需要它）。README 该句与实现不符。
- **修法**：让渲染可失败——`to_args(&self, group_names) -> Result<Vec<String>, Error>`
  （或加 `try_to_args`），把"无拼写"情形返 `Error::usage`；README 改为
  `to_args(&group_names)`。见 02 §3。
- **验收**：`to_args(default_kinds(n), &names)` 返回 `Err` 而非 panic；
  README 句子与签名一致。

### A11-4 门面横跨六种错误类型【编】【据】

- **严重度：P1**。后果形态：调用方被迫 `to_string()` 抹平错误分类。
- **reed 现状**：`run_*` 与 `Populations::group_names` 返 `reed::Error`
  （`request.rs:263`），但其余公开构造器各返不同错误：
  `Levels::parse` → `EngineError`（`request.rs:181`）、`assign_groups` →
  `EngineError`（`request.rs:784`）、`MissingRate::parse`／`Region::parse` →
  `FilterError`（`spec/missing_rate.rs:37`、`spec/region.rs:76`）、
  `WindowGrid::validate`／`ContigRegions::new` → `CoreError`（`grid.rs:90`/`:215`）、
  `AfKind::parse`／`GtClassRule::parse` → `String`（`spec/mod.rs:84`/`:290`）。
  CLI 的代价直接可见：`Levels::parse(...).map_err(|e| reed::Error::usage(e.to_string()))`
  （`reed-cli/src/translate/stats.rs:554`/`:774`）。
- **实测**（探针编译）：
  ```
  error[E0277]: the trait bound `reed::Error: From<String>` is not satisfied
    = help: `reed::Error` implements trait `From<T>`:
              From<reed_core::error::CoreError>  From<reed_engine::error::EngineError>
              From<reed_filter::error::FilterError>  From<reed_io::error::IoError>
              From<std::io::Error>
  ```
  即 `EngineError`/`FilterError` 尚可 `.into()` 进 `reed::Error`，
  `String`（两个 `parse`）不能。
- **修法**：门面自身的构造器（`Levels::parse`、`assign_groups`）返
  `reed::Error`；`String` 型 parse 改返 `FilterError::Spec` 或 `reed::Error`
  （二者择一，保持 `Error: From<_>` 全通）。见 02 §4。
- **验收**：六种来源都能 `.into()` 成 `reed::Error`（探针编译通过）；
  CLI 不再需要 `to_string()` 桥接。

### A11-5 `Error` 变体携带不可命名的内部类型【编】

- **严重度：P1**。后果形态：违反"下层 crate 是内部"的承诺。
- **reed 现状**：`pub enum Error { Core(CoreError), Io(IoError), Engine(EngineError),
  Spec(FilterError), Usage(String) }`，五变体全 `pub`（`error.rs:13-26`）；
  门面只再导出 `Error`/`ErrorKind`（`lib.rs:63`），**不导出**四个 payload 类型。
- **实测**（探针编译）：
  ```
  error[E0425]: cannot find type `EngineError` in crate `reed`
  error[E0425]: cannot find type `IoError` in crate `reed`
  ```
  只依赖 `reed` 的调用方能 `match`（类型靠推断）、能 `Display`，
  但**无法命名或存储**这些 payload 类型，也无法为它们写转换。README
  又说 "Only this facade carries a stability promise; the lower crates under
  `crates/` are internal"——可这些内部类型已随公开变体冻结进契约。
- **正是本仓自己的标准**：doc_v7 D7-22 把"pub 在私有模块、却没再导出、
  库调用方无法命名"判为缺陷；那批修了 `StatAcc`/`QcReportAcc`，
  `Error` 的 payload 漏了。
- **修法（已定，用户 2026-10-14：方案②b）**：`Error` 改 struct 包装
  （`pub struct Error(Box<ErrorInner>)`），变体对外不可见；只留
  `usage()`/`kind()`/`ErrorKind`/五个 `From`/`Display`/`source()`。改动面
  只有 `error.rs` 本体 + `index.rs:624` 一处（后者走 `From`）。**不采用** ②a
  （只加 `#[non_exhaustive]`：payload 仍在签名里、匹配变难，两头不靠）。
  完整代码与迁移点见 02 §5。
- **验收**：`cargo doc` 里 `Error` 无变体页；下游 `match err { reed::Error::Io(..) => .. }`
  编译失败；CLI 退出码与文案不变；`kind()` 单测原样过。

### A11-6 四驱动里 `run_reports` 与报告类型不在根，README 漏列【编】

- **严重度：P1**。后果形态：公开面按"四个驱动"承诺，根上只能发现三个。
- **reed 现状**：`run_filter`（`lib.rs:64`）、`run_stats`（`:67`）、`run_index`
  （`:65`）在根；`run_reports`／`ReportsOutcome` 只在 `reed::reports::*`
  （writers 模块再导出，`:82`）。报告类型同样分裂：`reed::StatsReport`／
  `reed::IndexOutcome`（根）对 `reed::funnel::FilterReport`／
  `reed::reports::ReportsOutcome`（子模块）。
- **实测**（探针编译）：
  ```
  error[E0432]: unresolved import `reed::run_reports`
  error[E0432]: unresolved import `reed::FilterReport`
  error[E0432]: unresolved import `reed::ReportsOutcome`
  （对照：`reed::reports::run_reports`、`reed::funnel::FilterReport`、
    `reed::reports::ReportsOutcome` 均解析成功）
  ```
- **附带（文档）**：README:172-175 的"re-exports everything needed"只列
  `run_filter`/`run_stats`/`run_index`，漏 `run_reports`（却把 `reed::reports`
  当写出模块列了）。
- **修法**：根补 `pub use ...::reports::{run_reports, run_reports_opened, ReportsOutcome};`，
  并把 `FilterReport` 也在根再导出（`run_filter` 的返回类型不该只住在叫
  `funnel` 的子模块）；README 驱动清单补 `run_reports`。见 02 §6。
- **验收**：探针能 `use reed::{run_reports, FilterReport, ReportsOutcome};`。

### A11-7 `ReportRequest::is_empty()` 漏 `summary`【实】

- **严重度：P2**。后果形态：谓词错误（`#[must_use]` 公开 API）。
- **reed 现状**：`is_empty` 只测五个字段（`reports.rs:65`），漏 `summary`。
- **实测**（探针程序运行）：
  ```
  summary-only is_empty() = true
  assert passed -> is_empty() returns true for a --summary-only request
  ```
- **范围**：全仓无调用者（CLI 直接按 flag 分支，grep 0 命中），
  所以是**潜伏缺陷**，不是现行 bug。
- **修法**：加 `|| self.summary`；若确无常驻消费者，删除该谓词。见 02 §7。
- **验收**：`ReportRequest { summary: true, ..default }.is_empty() == false`。

### A11-8 死再导出 `header`／`LdWindow`（`bgzf` 有真实消费者，已核减）【据】

- **严重度：P2（lean）**。后果形态：公开面背无关名字。
- **reed 现状**：`pub use reed_stats::ld::LdWindow;`（`lib.rs:78`）、
  `pub use reed_io::{bgzf, chunk, header};` 中的 `header`（`:98`）。
- **证据（2026-10-14 落地时更正）**：首轮复审用 `grep -r` 统计"全仓 0 引用"，
  但 ugrep 默认遵守 `.gitignore`，而 `crates/*/tests/` 被忽略——**统计漏了
  CLI 的集成测试**。改用 `grep --no-ignore-files` 复核后：
  - `reed::header`、`reed::LdWindow`：确实 0 引用 → 删除；
  - `reed::bgzf`：`reed-cli/tests/cli_contract.rs`（2 处）与
    `stdin_streaming.rs`（4 处）用它写 BGZF 夹具 → **保留**。
  编译器证实了这一点（删掉后 `cargo check --workspace` 报 6 个
  `cannot find bgzf in reed`）。
- **修法**：删 `LdWindow` 与 `header`，保留 `bgzf`/`chunk`。见 02 §8。
- **验收**：`grep -rn --no-ignore-files 'reed::LdWindow\|reed::header' crates`
  为空；`cargo check --workspace --all-targets` 通过。

> **方法教训（记入本轮）**：本仓 `.gitignore` 忽略 `crates/*/tests/` 与
> `docs/`，任何"全仓引用统计"都必须带 `--no-ignore-files`，否则结论会
> 系统性漏掉测试代码。A11-8 的首版结论就是这么错的。

### A11-9 `OpenedInput` 暴露传输层内部类型【编】

- **严重度：P2**。后果形态：内部类型进入稳定契约。
- **reed 现状**：`pub struct OpenedInput { pub source: Box<dyn ChunkSource>,
  pub facts: HeaderFacts, pub bcf: Option<Arc<BcfDict>> }`（`input.rs:24-32`），
  且 `run_stats_opened`（`driver.rs:43`）、`run_reports_opened`（`reports.rs:748`）
  公开接收它。
- **实测**（探针编译通过）：`fn touch(o: &reed::OpenedInput){ &o.facts; &o.bcf; &o.source }`
  ——三字段确实公开。`ChunkSource`/`HeaderFacts`/`BcfDict` 是 `reed_io`/`reed_engine`
  内部量，靠门面再导出才勉强可命名。
- **修法**：给 `OpenedInput` 一个 opaque 形态（字段私有 + 受控构造），
  保留 `_opened` 入口的流式能力；或把这一组"传输词汇"明确列为 §6.2 的批准例外。
  见 02 §9。
- **验收**：`OpenedInput` 字段非 `pub`，`_opened` 入口仍可用。

### A11-10 测试专用 `ExecStrategy` 公开【编】

- **严重度：P3**。后果形态：测试钩子进入稳定 API。
- **reed 现状**：`FilterOptions.strategy: Option<ExecStrategy>`（`filter/mod.rs:68`），
  字段注释写明 "Test-only strategy override"；`ExecStrategy` 经
  `reed::plan` 公开。CLI 用它做策略注入（`translate/filter.rs:349`）。
- **实测**：探针 `reed::plan::ExecStrategy::SinglePass` 解析成功。
- **修法**：`#[doc(hidden)]`（非破坏）或移入 `#[cfg(any(test, feature = "test-hooks"))]`。
  见 02 §10。
- **验收**：`cargo doc` 不列 `ExecStrategy`；CLI 测试仍可用。

### A11-11 缺内存态 `Populations`【编】

- **严重度：P2**。后果形态：程序化调用者被迫落临时文件。
- **reed 现状**：`enum Populations { All, File(PathBuf) }`（`request.rs:243-250`）——
  给群体命名只能给 TSV 路径；`assign_groups` 也只接路径（`request.rs:784`）。
- **实测**：探针 `reed::Populations::Assignment(...)` → `E0599: no variant ... named Assignment`。
- **背景**：README 自我定位是"other tools can drive the same statistics
  without shelling out"；一个手上已有 `sample→group` 映射的调用方却必须写盘。
- **修法**：加 `Populations::Assignment(Vec<(String, String)>)`（或等价内存表），
  `assign_groups` 加一个接内存表的兄弟（复用 `GroupTable` 的行规则）。见 02 §11。
- **验收**：探针能用内存映射跑通 `run_stats` 且与 `File` 版本逐字节一致。

### A11-12 `WindowRow.end0` 名字与文档不符【读】

- **严重度：P3**。后果形态：命名歧义。
- **reed 现状**：字段 `end0`（`request.rs:232`）文档写 "1-based inclusive end"；
  值来自 `interval.end.0`（0-based 半开闭端，`sink.rs:423`），数值恰等于
  1-based 闭端，故两种读法同值，但名字 `0` 后缀读作 0-based。
- **依据**：v2/05-naming.md §4 的命名纪律——名字应反映语义。
- **修法**：改名 `end`（更新 writers 的 `row.end0` 引用），或文档同时写明两种读法。
  见 02 §12。
- **验收**：字段名与文档一致；writers 全绿。

### A11-13 `render` 的 `flag`／`f` helper 重复【读】

- **严重度：P3（lean）**。后果形态：近似副本须同步维护。
- **reed 现状**：`fn flag(args, name, value)` 在 `stats/render.rs:15` 与
  `spec/render.rs:19` 各有一份；`fn f(value: f64) -> String` 只在后者（`:25`）。
- **修法**：合并到一处共享（如 `reed-core` 的小助手或门面 `render` 公共模块）。
  见 02 §13。
- **验收**：`flag` 单处定义；两处调用点编译通过。

---

## 2. 已核减 / 需说明的说法

| 说法 | 处置 | 依据 |
|---|---|---|
| "六种错误类型都要统一成 `reed::Error`" | **部分核减**：`MissingRate`/`Region`/`WindowGrid` 等是下层 crate 的公开类型，其构造器错误可保持下层类型；**本轮只要求门面自己的构造器**（`Levels::parse`、`assign_groups`）与 `String` 型 parse 收敛（A11-4） | §1 A11-4；门面 `From` 链仅缺 `String` |
| "`run_*_opened` 是多余 API，应删" | **不改**：它是 v6 Y6.5 为 stdin 流式加的（一次打开、多次消费），有真实前端依赖（`reed-cli` 5 处）。本轮只处理它携带的 `OpenedInput` 内部量（A11-9） | §1 A11-9；`reed-cli` 用法 |
| "`default_kinds`／`window_kinds` 不该公开" | **不改**：CLI 与测试都在用（`default_kinds` 2+2 处、`window_kinds` 2+8 处）；问题在它们与 `to_args` 不可组合（A11-3），修 `to_args` 即可 | §1 A11-3；用法计数 |
| "`chunk`／`metrics`／`budget` 也应删" | **核减**：三者有真实消费者（`chunk` 5、`metrics` 1、`budget` 7 处）；只有 `bgzf`/`header`/`LdWindow` 是死的（A11-8） | §1 A11-8；用法计数 |

---

## 3. 复核门（本轮新增，落地后由 `ci/check-doc-drift.sh` 登记）

1. **A11-1 守卫**：`request.rs` 含 `impl Default for Levels` 且 `window: true`
   （或 `Levels` 不再派生 `Default`）。
2. **A11-2 守卫**：`run_stats` 签名的请求参数为 `&StatsRequest`。
3. **A11-3 守卫**：`render.rs` 不含 `panic!("no --flag spelling`（to_args 已可失败）；
   README 的 `to_args` 句含 `group_names` 或改写。
4. **A11-4 守卫**：`! grep -q 'no --flag spelling for per-bp'`… 并入 #3；
   另加：门面可 `From<String>`（或两个 `String` parse 已改类型）。
5. **A11-6 守卫**：`lib.rs` 根含 `run_reports`、`ReportsOutcome`、`FilterReport`。
6. **A11-8 守卫**：`grep -rn 'LdWindow\|reed_io::bgzf\|reed_io::header' crates/reed/src/lib.rs` 为空。
7. **README 示例守卫**：一段把 README「Library use」示例编译的下游/doctest 门
   （现仓库 README 不在任何编译门内，A11-2 就是这样漏掉的）。
