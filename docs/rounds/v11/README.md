# reed doc_v11：门面 API 复审 + 修复规范（2026-10-14）

> 本目录是 `docs/rounds/v10/` 之后的**第十一轮规范**，主题是 **`reed` 门面
> （`crates/reed/src/lib.rs` 的公开面）的 API 复审**：拿 v5/08 §6.2 的
> "后端公开面"契约、README「Library use」的对外承诺、以及 `reed-cli`
> 与 `crates/reed/tests/*` 的实际用法，对照当前代码事实逐条核对。
>
> 复审对象：commit `065175c`（doc_v10 D10-1…D10-8 全部落地），工作区干净。
> 本轮**不修改任何代码**；所有修法以规范形式给出（02），实施顺序见 06。
> 方法：门面公开项通读 + 下游探针 crate 编译 + 临时集成测试实跑；每条
> 缺陷给出可复现的触发方式（§0 的证据通道标记）。

---

## 背景

v7–v10 修的是"行为"：错误路径、线程预算、终止语义、格式兼容。
`docs/rounds/v5/08-api-surface.md` §6 早已把**门面的公开面**定成契约
（§6.2 的请求／运行／结果／写出／错误／预算／文本契约清单），
但那一轮之后公开面又长了：`run_stats_opened`／`run_reports_opened`、
`default_kinds`／`window_kinds`、`assign_groups`、`LdWindow`、以及
为测试补的 `bgzf`／`header` 再导出，都没有被 §6.2 复核过。

本轮把这条边界重新量了一遍，发现的问题有一个共同模式：**门面在"CLI
这条唯一前端上"是对的——因为 CLI 每一处都显式绕开了默认值和缺失的类型；
一旦换个"只依赖 `reed` 的库调用方"，默认值、错误类型、报告类型放置、
README 示例就会成片失效。** README 自己举的那个库例子，今天既不编译、
就算编译也不产出窗口行。

分档：

| 档 | 缺陷 | 后果形态 |
|---|---|---|
| P0 | A11-1 `StatsRequest::default()` 产出 0 行 | **静默错误结果**（照 README 写库调用，得到空表） |
| P1 | A11-2 `run_stats` 按值取请求，README 示例不编译 | 对外承诺的用法编译失败 |
| P1 | A11-3 `to_args` 对本仓 `default_kinds` panic | 两个公开 helper 不能组合，进程 panic |
| P1 | A11-4 门面横跨六种错误类型 | 调用方被迫 `to_string()` 抹平错误 |
| P1 | A11-5 `Error` 变体携带不可命名的内部类型 | 违反"下层 crate 是内部"的承诺 |
| P1 | A11-6 `run_reports`／报告类型不在根，README 漏列 | 公开面按"四个驱动"承诺、只能从根发现三个 |
| P2 | A11-7 `ReportRequest::is_empty()` 漏 `summary` | 谓词错误（`#[must_use]` 公开 API） |
| P2 | A11-8 死再导出 `bgzf`／`header`／`LdWindow` | 公开面背无关名字（lean） |
| P2 | A11-9 `OpenedInput` 暴露传输层内部类型 | 内部类型进入稳定契约 |
| P2 | A11-11 缺内存态 `Populations` | 程序化调用者被迫落临时文件 |
| P3 | A11-10 测试专用 `ExecStrategy` 公开 | 测试钩子进入稳定 API |
| P3 | A11-12 `WindowRow.end0` 名不副实 | 命名歧义（v2 05 §4 命名纪律） |
| P3 | A11-13 `render` 的 `flag`／`f` 两处重复 | 可合并（lean） |

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 门面 API 复审：逐条缺陷 A11-1…A11-13（`file:line` + 严重度 + 触发方式 + 修法 + 验收）、已核减的说法、复核门 |
| [02-architecture.md](02-architecture.md) | 修订规范：请求默认值（N1）、驱动签名与放置（N2/N6）、渲染契约（N3）、错误面（N4/N5）、谓词（N7）、公开面收口（N8/N9/N10）、缺失 API（N11）、命名（N12/N13）；每条带验收标准 |
| [06-roadmap.md](06-roadmap.md) | 修复路线 R0–R3：探针/红测试先行、P0+P1 修复、公开面收口、文档与守卫收尾；通用验收门与"本轮不做"清单 |
| [STATUS.md](STATUS.md) | 实施状态（起始全部未修）+ 复审基线（每条缺陷的实测/编译证据摘录） |

## 结论速览

| 维度 | 结论 |
|---|---|
| 默认值（A11-1，P0） | `Levels` 用 `#[derive(Default)]` 得全 `false`，而字段文档与 `StatsRequest.levels` 文档都写"默认 window 行"；sink 只在 `levels.window` 时出行 → `StatsRequest::default()` 的 `windows` 为空。README 库示例正踩此坑 |
| 驱动一致性（A11-2/6，P1） | 四驱动里 `run_stats` 唯一按值取请求（README 示例因此 move 报错）；`run_reports` 唯一不在根、README 唯一漏列 |
| 渲染契约（A11-3，P1） | `to_args` 对 `per_bp`/`Sfs` 直接 `assert!`/`panic!`，而 `default_kinds` 恰好产出这些 → 两个公开 helper 不能组合；README 又写成 `StatsRequest::to_args()`（实际需 `group_names`） |
| 错误面（A11-4/5，P1） | 门面构造器横跨 `reed::Error`/`EngineError`/`FilterError`/`CoreError`/`String`/`io::Error` 六种；`Error` 的 `From` 链缺 `String`；`Error` 四个变体的 payload 类型（`CoreError` 等）在门面**不可命名**——正是 v7 D7-22 立的标准，那批漏了 `Error`。**A11-5 已定方案②b（变体封口）** |
| 公开面收口（A11-8/9/10，P2/P3） | `bgzf`/`header`/`LdWindow` 全仓 0 引用；`OpenedInput` 三个 `pub` 字段是传输内部量；`ExecStrategy` 是测试钩子。均属 v5/08 §6.2 收口范围 |
| 缺失 API（A11-11，P2） | `Populations` 只有 `All`／`File`，程序化调用者无内存态入口 |
| 命名（A11-12/13，P3） | `WindowRow.end0` 文档写"1-based 闭端"、名字读作 0-based；`flag`/`f` helper 在 stats 与 spec 两处重复 |

## 与 v2、v5 的关系

- v2/05-naming.md §4 立了"禁止版本/新旧后缀"的命名纪律；§3 曾规定库驱动路径为
  `reed::filter::run`／`reed::stats::run`／`reed::index::run`（"不写 `run_filter`"）。
  该路径条已被 v5/08 §6.2 取代（门面用 `run_*`），但 §4 的命名纪律仍然有效——
  A11-12 即援引它。
- v5/08 §6.2 是本轮的量尺：请求／运行／结果／写出／错误／预算／文本契约。
  A11-6/8/9/10 都是对它的**回归**（公开面长出 §6.2 未列、也不该列的东西）。
  本轮不重开 §6 的分工裁决，只把公开面收回到它。
- v5/08 §6.4 守卫 3（`lib.rs` 无 `pub mod`）当前成立（`grep -c '^pub mod' = 0`），
  本轮新增守卫登记于 01-review.md §3。

## 决策状态

| 决策 | 状态 |
|---|---|
| A11-5 `Error` 变体可见性 | **已定（用户，2026-10-14）：方案②b**——`Error` 由公开 enum 改 struct 包装，变体对外不可见，只留 `usage()`/`kind()`/`ErrorKind`/`From`/`Display`/`source()`。判据：cli.md §9 只承诺两类出口、全仓无外部变体匹配（封口零破坏）；完整代码与迁移点（`error.rs` + `index.rs:624`）见 02 §5 |

## 标记约定

【规】必须；【荐】推荐；【决】需要你拍板。
证据通道：**【实】**实跑复现（临时集成测试／探针程序运行）；**【编】**编译期复现
（下游探针 crate 的 `cargo build` 诊断）；**【读】**源码或文本证据（不可实跑，如命名）。
**P0/P1/P2/P3** 为修复优先级。
