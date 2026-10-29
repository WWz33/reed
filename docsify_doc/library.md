# 库调用

`reed` crate（`crates/reed`）是公开 API，命令行是它的一层封装。

```toml
[dependencies]
reed = { path = "crates/reed" }
```

稳定性承诺只在 crate 根；`crates/` 下层的 crate 是内部实现，类型与签名随时可能变。错误用 `reed::Error` 表示，`Error::kind()` 分类，`Display`/`source()` 读消息。

## 入口函数

crate 根导出六个入口，各自对应一个子命令：

```rust
run_filter, run_popgen, run_stats, run_rename, run_index, run_convert
```

签名模式统一：`fn run_*(input: Option<&Path>, request: &XxxRequest, budget: &ThreadBudget) -> Result<XxxReport>`。`ThreadBudget`（`reed::budget`）给线程数与内存预算，`Default` 即"全部核 + 物理内存一半"，与 CLI 缺省一致。

## 群体统计

```rust
use reed::{PopgenRequest, Populations, StatKind};

let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::File("pops.txt".into()),
    ..PopgenRequest::default()
};
let report = reed::run_popgen(Some("in.vcf.gz".as_ref()), &request,
                              &reed::budget::ThreadBudget::default())?;
let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false);
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
let report = reed::run_convert(Some("in.vcf.gz".as_ref()), &request,
                               &reed::budget::ThreadBudget::default())?;
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

写出模块也直接可用：`reed::windows`（窗口 TSV）、`reed::reports`、`reed::histogram`、`reed::sfs`、`reed::fstats`、`reed::ld`、`reed::summary`。

## 往返渲染

`PopgenRequest::to_args(&group_names)` 与 `FilterSpec::to_args()` 把请求渲染回 CLI 拼写，且 `parse(render(x)) == x`（有测试锁定）。用处：把库层构造的请求记进日志、或在库里拼好请求后交给 CLI 进程执行。
