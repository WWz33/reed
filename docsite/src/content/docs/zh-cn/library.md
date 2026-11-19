---
title: 库调用
---

`reed` crate（`crates/reed`）是公开 API，命令行是它的一层封装。

```toml
[dependencies]
reed = { path = "crates/reed" }
```

稳定性承诺只在 crate 根；`crates/` 下层的 crate 是内部实现，类型与签名随时可能变。错误用 `reed::Error` 表示，`Error::kind()` 分类，`Display`/`source()` 读消息。

## 稳定性

根模块条目与 `run_*` 入口遵循 semver（0.x：破坏性变更走 minor）；`#[doc(hidden)] reed::frontend` 与 `crates/` 下的其余 crate 均为内部实现。见 `docs/spec/api-contract.md`。

## 入口函数

crate 根每个子命令导出一个 `run_*`（`_opened` 变体接收已打开的 `OpenedInput`，
完整清单见 `ci/api/reed.txt`）：

| 入口 | 请求 | 报告 |
|---|---|---|
| `run_filter` | `FilterSpec` | `FilterReport` |
| `run_popgen` | `PopgenRequest` | `PopgenReport` |
| `run_stats` | `StatsReportRequest` | `StatsOutcome` |
| `run_rename` | `RenameSpec` | `RenameReport` |
| `run_index` | `WriteIndex`（参数） | `IndexOutcome` |
| `run_convert` | `ConvertRequest` | `ConvertReport` |
| `run_pca` | `PcaRequest` | `PcaReport` |
| `run_dapc` | `DapcRequest` | `DapcReport` |

签名模式统一：`fn run_*(input: Option<&Path>, request: &XxxRequest, ctx: &RunContext<'_>) -> Result<XxxReport>`。`RunContext`（`reed::RunContext`）借用本次运行的 `ThreadBudget` 与 `AbortSignal`：有取消源时用 `RunContext::new(&budget, &abort.signal)`（`Abort` 来自 `reed::abort_for(&budget)`），没有时用 `RunContext::idle(&budget)`。`ThreadBudget::default()` 即"全部核 + 物理内存一半"，与 CLI 缺省一致。这些名字也可用 `use reed::prelude::*;` 一次导入；prelude 只再导出根项（`run_*` 入口、请求、报告、错误、预算与 `RunContext`）。

## 群体统计

```rust
use reed::{PopgenRequest, Populations, StatKind};

let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::File("pops.txt".into()),
    ..PopgenRequest::default()
};
let budget = reed::ThreadBudget::default();
let ctx = reed::RunContext::idle(&budget);
let report = reed::run_popgen(Some("in.vcf.gz".as_ref()), &request, &ctx)?;
let tsv = reed::render_tsv(&report.windows, &request.kinds, false);
```

群体映射可以不落盘：

```rust
let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::Assignment(vec![
        ("sample1".into(), "popA".into()),
        ("sample2".into(), "popB".into()),
    ]),
    ..PopgenRequest::default()
};
```

`StatKind` 的变体覆盖 CLI 的全部统计量（`Pi`、`Dxy`、`Fst`、`WattersonTheta`、`TajimaD`、`Pbs`、f 统计量族、`Sfs`、LD 族）；CLI 拼不出来的请求（per-bp 统计、`Regions` 网格）在 `to_args()` 渲染时报错，不影响直接调用。

## 格式转换

```rust
use reed::{ConvertRequest, ConvertTarget, PhyloFormats, PhyloOptions};

let request = ConvertRequest {
    targets: vec![ConvertTarget::Phylo {
        formats: PhyloFormats { phylip: true, fasta: true, nexus: false },
        options: PhyloOptions::default(),
    }],
    prefix: "matrices".into(),
    ..ConvertRequest::default()
};
let budget = reed::ThreadBudget::default();
let ctx = reed::RunContext::idle(&budget);
let report = reed::run_convert(Some("in.vcf.gz".as_ref()), &request, &ctx)?;
```

## 请求与报告

| 请求类型 | 报告类型 |
|---|---|
| `PopgenRequest` | `PopgenReport` |
| `FilterSpec` | `FilterReport` |
| `StatsReportRequest` | `StatsOutcome` |
| `RenameSpec` | `RenameReport` |
| `ConvertRequest` | `ConvertReport` |
| （index 无请求结构） | `IndexOutcome` |

`render_tsv` 渲染窗口统计报告；其余写出模块只属于前端（位于 `#[doc(hidden)] reed::frontend`，不在稳定性承诺内）。

## 往返渲染

`PopgenRequest::to_args(&group_names)` 与 `FilterSpec::to_args()` 把请求渲染回 CLI 拼写，且 `parse(render(x)) == x`（有测试锁定）。用处：把库层构造的请求记进日志、或在库里拼好请求后交给 CLI 进程执行。
