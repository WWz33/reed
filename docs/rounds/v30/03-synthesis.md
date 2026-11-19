# 03 综合诊断：79 条 → 5 个根因 + 1 个裁决

> 2026-11-11。用户定调：v30 的 79 条（V30-01…09 9 条 + V30-Q01…Q36 36 条 +
> V30-R01…R14 14 条 + V30-T01…T20 20 条）**不应逐条设 issue**，要综合起来找根本方案。
>
> 方法：先按"症状簇"归并，再问"删掉哪个所有者，这一簇会同时消失"；最后用已安装的
> 六件工具 + `codegraph` 给每条根因配一个**机器可验证的不变量**。本文件是路线图，
> 不改代码；逐条 issue 降级为验收清单。

## 1. 结论先行

79 条不是 79 个问题。它们收敛为 **5 个"缺一个唯一所有者"**，外加 **1 个未裁决的产品问题**：

| 根因 | 缺的所有者 | 条数 | 一句话 |
|---|---|---:|---|
| **A 契约面** | 门面没有"设计出来的 API" | **17** | 下层实现的并集被当成公共 API |
| **B 语义** | 领域规则没有唯一类型/函数 | **15** | GT/等位/缺失/坐标/区域/INFO 各层重写 |
| **C 状态与资源** | 能力/预算/临时文件/线程没有 owner | **18** | 靠路径、`Option`、运行期分支维持 |
| **D 演进** | 没有可执行的契约门 | **16** | semver/public-api/external-types/feature/CI 缺位 |
| **E 文档与结构** | 文档/规范不是派生品 | **9** | 手写副本 + 内部票号 + 规格漂移 |
| **F 产品裁决** | 未决定"发布还是 git-only、稳定面是什么" | **4** | 决定 A/D 的目标形态 |

**根本方案 = 1 页契约裁决（F）+ 4 个所有者重构（A/B/C/D/E 归并到 4 个改动）+ 一套机器门。**
做完后逐条 issue 作为验收用例自动关闭；没有"第 80 条"。

## 2. 为什么是这 5 个根因（判定规则）

一条症状归入某根因，需同时满足：

1. **删除该根因会让症状按构造消失**（不是"顺便修好"）；
2. **同一根因下的多条症状共用同一个修复动作**；
3. **该根因有一个可用工具/编译门验证的不变量**。

按此规则得到下面 79 条的完整归属（每条只归一个主根因；跨根因的在正文标注）：

| 根因 | 成员 | 数量 |
|---|---|---:|
| **F 产品裁决** | Q29、Q33、T01、T04 | 4 |
| **A 契约面** | Q01、Q02、Q03、Q04、Q06、Q09、Q10、Q11、Q12、Q13、Q14、Q16、Q36、T02、T03、T05、T08 | 17 |
| **B 语义** | V30-01、V30-02、V30-03、Q05、Q23、R01、R02、R03、R04、R12、R13、R14、T09、T10、T11 | 15 |
| **C 状态与资源** | Q07、Q08、Q15、Q17、Q18、Q19、Q20、Q21、Q22、R05、R06、R07、R08、R09、R10、R11、T19、T20 | 18 |
| **D 演进** | Q24、Q26、Q28、Q30、Q31、Q32、Q34、Q35、V30-09、T06、T07、T12、T13、T16、T17、T18 | 16 |
| **E 文档与结构** | Q25、Q27、V30-04、V30-05、V30-06、V30-07、V30-08、T14、T15 | 9 |
| **合计** |  | **79** |

## 3. 因果链：症状不是并列的

```
F 产品裁决（publish=false / 稳定面未定义）
 ├─→ A 契约面失控（整 crate 再导出、不可命名类型、外部类型）
 │     └─→ D 演进门缺位（没有快照就无从 semver；workspace 跳过 7 个 crate）
 │           └─→ E 文档漂移（无 doc 门；README/library.md 手写）
 ├─→ B 语义无主（规则留在 L1/各消费方）
 │     ├─→ C 状态与资源隐式（每个 driver 自己处理 capability/预算/临时文件）
 │     └─→ A 的类型多路径（Region/Stand/Mat 随实现分散）
 └─→ D 的发布/feature 门缺位（可发布三 crate 打包全失败）
```

三条最有代表性的"一棵树"：

1. **GT/等位语义**（B）：`count_alts` 声明"唯一实现"但有 5 份（`site.rs:22`、`summary.rs:86`、
   `reports.rs:634`、`convert/job.rs:824`、`genotype.rs:263`），convert 又用
   `!raw.starts_with(b".")` 自定"有数据"（`convert/job.rs:851`）。同一根因长出的症状：
   R01（半缺失分叉）、Q05（Region 三份的同类问题）、T02/T05（类型不可命名/外部类型）、
   V30-03（多等位三路径）。修 5 个副本中的任意 1 个都不会消失。
2. **AC/AN**（B）：`filter::semantics::info_af_indexed` 是"唯一读取器"，但 `stats --summary`
   手写第二份（`summary.rs:294-328`），`.afreq` 又走 GT 计数（`reports.rs:843-857`）。
   一个位点、一条命令、两个 AC（0 vs 0.61），debug 下溢出 panic。症状：R02、Q02、Q08、Q14、
   R13（6 份标签枚举）、R14（3 个 `EditPlan` 组装点）、T09（filter→io）。
3. **能力/资源**（C）：`OpenedInput { source, path: Option<PathBuf> }` 把 stdin 与文件同构，
   于是 three-pass 选择无法询问"能重读吗"（R05），错误分类按路径分叉（R06），
   `source_mut`/`bcf.is_some()`/按路径重开成为三种入口（R07）；再往下是临时文件只拥有路径（R10）、
   预算在 8 处各自推导（R09）、线程 detach（R08/T19）、输出双 owner（R11/T20）。

## 4. 根本方案（P0–P4）

### P0 契约裁决（先做，1 页 ADR，不改代码）

**产出**：`docs/spec/api-contract.md`，必须回答：

1. **发布面**：`reed` 门面 + 哪些闭包发布？还是明确"仅 git 依赖"？`reed-cli`/`oracle` 永不发。
2. **稳定面**：只承诺 `reed` 根 + `prelude`；下层 crate 标注 `Internal`（不承诺 semver），
   还是逐 crate 承诺（则每个都要 semver 检查）。
3. **外部依赖政策**：noodles/serde/memmap2 是 public dependency（显式 re-export）还是 opaque（不进入签名）？
4. **MSRV / edition / feature 政策**：1.91、2024、feature 只用于可选后端（faer/serde/noodles）。
5. **兼容承诺**：CLI 输出 byte-identical；库 API semver；内部 crate 无承诺。

没有这一步，A 和 D 的目标形态无法确定（例如 semver 要不要覆盖 11 个 crate），任何重构都可能返工。

### P1 语义内核：让 B 只有一处定义（reed-core）

**动作**（把规则从 L1/消费方收回 L0，并给 L3 一个不依赖 L1 的视图）：

| 规则 | 唯一所有者 | 删除的副本 |
|---|---|---|
| GT 分类/半缺失 | `reed_core::genotype::{GtCode,GtClass,has_called_allele}` | convert 的 `present`、pca 的 dosage 特判 |
| ALT 计数 | `reed_core::site::AltCount` + `count_alts` | summary/reports/convert/genotype 的 4 份 |
| AC/AN/INFO | `reed_core::count::InfoCounts::from_info` | filter `semantics` 与 stats `summary` 的两份 |
| 坐标 | `reed_core::coord::{Pos0,Pos1,Interval1}` | 4 处 `pos0` 转换、三处裸 `+1/-1` |
| 区域 | `reed_core::region::Region::parse_in(text, contigs)` | filter/io/convert 三份 Region + 两份包含判定 |
| 多等位 | `reed_core::site::MultiallelicPolicy` | stats 折叠 / popg 跳过 / filter 报错三条路径 |
| 站点视图 | `reed_core::site::SiteFields` trait（文本/BCF 两实现） | 5 条 BCF native 入口、convert 自写 tokenizer |

**结构性结果**：`reed-filter` 只依赖 `reed-core`（删除 `reed-io` 依赖），
`cargo modules dependencies` 中 `reed-filter -> reed-io` 边消失；`reed` 的 driver 不再有
`bcf.is_some()` 分支。

**关闭**：B 的 15 条 + A 的 Q02/Q05/T05 中与语义类型相关的部分。
**工具门**：`cargo modules dependencies --lib -p reed-filter` 不含 `reed-io`；
`rg` 计数（ALT 5→1、pos0 4→1、Region 3→1、EditPlan 3→1）；半缺失/溢出/POS 哨兵的单测。

### P2 能力与状态类型：让 C 只有一处定义

**动作**：

1. **输入能力成类型**：`StreamedInput`（一次）/ `ReplayableInput::replay(regions)` /
   `SeekableInput`（`-r`，需 index）；`InputMeta { format, framing, bcf }` 由 `into_parts` 唯一带出；
   删除 `source_mut` 与按路径重开。
2. **运行上下文**：`RunContext { budget: &ThreadBudget, abort: &AbortSignal }` 统一传给所有 driver；
   删除 `PcaRequest/DapcRequest` 的 `threads/memory` 死字段。
3. **资源 owner**：`PrunePool`（Drop 先取消再 join，覆盖流式与 batch 两个 prune 引擎）；
   `TempFile` 持有 `File`（读写 `try_clone`，注册进 `AtomicFile` 的存活表）；
   `OutputTarget` 拆 `create` / `from_parts`（`adopting` 加断言或删除）。
4. **预算注入**：`ThreadBudget`（含 `MemBudget`）单向进入 `choose_for_mind`、`ConvertSink`、LD 矩阵；
   删除 `default_half_memory` 的 8 处重复推导。

**关闭**：C 的 18 条（含 Q07/Q08/Q15/Q17-Q22）。
**工具门**：`strace` 统计 `-r` 只 open 1 次；stdin + `--strategy two_pass_reparse` 报 usage；
故障注入后 prune 线程归零；symlink 同址不改写 victim；`--memory 64 --mind` 的策略反映 64 MiB；
`cargo hack --no-default-features` 绿。

### P3 契约面与机器门：让 A/D/F 有唯一入口

**动作**：

1. **门面**：逐项再导出（删除 13 条整 crate/模块别名 + 7 个 writers 整模块），只保留
   14 个 `run_*` + 请求/报告/错误/`prelude`；`#[doc(hidden)]` 只用于 `frontend` 模块。
2. **可命名性**：用 `cargo public-api` 的输出做编译门——公共签名出现的每个类型必须能从
   `reed::*` 命名；`StatOutput`/`EngineError`/`BcfDict` 要么导出要么 opaque。
3. **快照与 semver**：`ci/api/<crate>.txt` 入库，CI diff；`cargo semver-checks -p <crate>`
   逐 crate（或按 P0 的稳定面）；`cargo package --dry-run` 对发布目标。
4. **外部类型**：`cargo check-external-types`（固定 `nightly-2026-03-20`）+ allow-list；
   noodles/serde/memmap2 的政策由 P0 决定。
5. **feature/CI**：`cargo hack --each-feature` + `--no-default-features`；CI 的 main 与
   PR→main 同权跑 test/clippy/doc/fmt/deny/structure。

**关闭**：A 的 17 条 + F 的 4 条。
**工具门**：public-api diff、semver-checks、check-external-types、cargo-hack、cargo package、
CI 条件检查——六件工具全部变成门。

### P4 派生文档与结构：让 D/E 不再手写

**动作**：

1. **文档从代码派生**：每 crate ≥1 个可运行 doctest；门面 ≥1 个用户示例；
   `library.md` 的入口表由 `pub use`/`run_*` 生成；README 加最小库示例；
   `# Stability` 与三处口径统一；公开 rustdoc 的 `doc_vNN` 迁入 `# Implementation notes`。
2. **结构**：`RecomputeTag` 静态表（名字/ty/number/description/compute）驱动 filter/io/cli；
   `SurvivorEdit` 纯值替代 3 个 `EditPlan` 组装点；`convert/job.rs` 按格式拆分；
   architecture.md 补齐 distance/reduce/graph/oracle 的层级与 `reed-qc` 去向。
3. **门**：`RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --workspace`、
   `cargo test --doc --workspace`、`cargo modules orphans`、单文件行数上限、文档票号守卫。

**关闭**：D 的 16 条 + E 的 9 条。
**工具门**：rustdoc exit 0；每 crate doctest ≥1；`cargo modules orphans` 0；doc-drift 新增守卫。

## 5. 顺序与依赖（不要并行 79 条）

| 阶段 | 内容 | 依赖 | 为什么这个顺序 |
|---|---|---|---|
| **0** | P0 裁决 + 三个便宜门（public-api 快照、CI 修条件、doc/no-default 门） | 无 | 门先于重构，保护后续所有改动 |
| **1** | P1 语义内核 | P0 的稳定面决定公开类型 | 类型先于门面；否则门面导出的是错误类型 |
| **2** | P2 能力与状态 | P1 的 `SiteFields` 接口 | driver 简化后才能统一 capability |
| **3** | P3 门面与门（全量） | P1/P2 的类型稳定 | 逐项导出需要最终类型清单 |
| **4** | P4 文档与结构 | P3 的导出清单 | 入口表/示例依赖最终门面 |

**可并行**：P4 的票号清理、architecture.md 更新、`RecomputeTag` 表可与 P1 并行（不碰热路径）。
**不可并行**：`convert/job.rs` 拆分与 P1 的 GT 重构会冲突，先语义后拆分。

## 6. 验收：工具，而不是"感觉"

| 不变量 | 工具 | 命令 | 通过条件 |
|---|---|---|---|
| 公共面冻结 | cargo-public-api | `cargo public-api -p reed -s` + diff `ci/api/reed.txt` | 新增/删除必须显式更新 |
| 类型可命名 | cargo-public-api / rustdoc JSON | 编译门 + 探针 | 公共签名类型全部可从 `reed::*` 命名 |
| semver | cargo-semver-checks | `-p <crate> --baseline-rev <tag>` | 无 `Skipping`，无 breaking |
| 外部类型 | cargo-check-external-types | 固定 nightly + allow-list | 只剩批准项 |
| 分层 | cargo-modules | `dependencies --lib -p reed-filter` | 无 `reed-filter -> reed-io` 边 |
| 模块卫生 | cargo-modules | `orphans --lib -p <crate>` | 0 孤儿；模块 SCC = 0 |
| feature | cargo-hack | `--each-feature`、`--no-default-features` | 全绿 |
| 文档 | rustdoc / doctest | `-D warnings`、`test --doc --workspace` | exit 0；每 crate ≥1 |
| 结构 | 自写守卫 | 文件行数、票号、单实现计数 | 全部 0/1 |
| 资源 | 故障注入/strace | 取消、临时文件、open 次数 | 线程归零、victim 不变、open=1 |
| 发布 | cargo | `package --dry-run` | 发布目标全绿 |
| CI | workflow 检查 | main 与 PR→main | 都跑 test/clippy/doc |

## 7. 关闭规则：79 → 0

1. **按构造关闭**：P1 完成后，B 的 15 条中 12 条由"唯一实现 + 单测"直接关闭；
   剩余 3 条（V30-01/02/03）需要各自的数值对拍/行为测试。
2. **按门关闭**：P3 完成后，A 的 17 条中 12 条由 public-api/semver/external-types 门直接关闭；
   5 条（Q03/Q04/Q16/Q36/T08）需要显式清单/依赖清理。
3. **按类型关闭**：P2 完成后，C 的 18 条中 12 条由类型构造关闭；6 条（Q07/Q08/Q19-Q22）
   需要 API 形状调整。
4. **按文档门关闭**：P4 完成后，D/E 的 25 条中 20 条由 doc/doctest/structure 门关闭；
   其余是内容修正。
5. **未关闭的条目继续以测试用例存在**，不再作为"issue 列表"维护；每条必须有对应的
   失败测试或门命令。

## 8. 明确不做

- **不逐个修 79 条**：它们是验收清单，不是工作分解。
- **不做"每个参数一个 newtype"**：只对跨层语义（GT/ALT/AC-AN/坐标/区域/预算/能力）建类型。
- **不重写算法**：V30-01（PI_HAT 系数）等数值缺陷按原路线单独修，不混入结构重构。
- **不在 P0 前发布/打 tag**：semver 基线要等稳定面确定。
- **不为了分层而分层**：两处模块级环（engine abort↔budget、filter expr↔semantics）可接受，
  只在改动相邻代码时顺手解。

## 9. 风险与回退

| 风险 | 缓解 |
|---|---|
| P1 改语义内核触碰热路径 | 先加"半缺失/溢出/哨兵/多等位"四组夹具；gold 输出逐字节对拍 |
| P2 改 capability 影响所有命令 | 保留 `_opened` 入口作为过渡；CLI 输出 golden 不变 |
| P3 公开面收口造成下游破坏 | 快照 + semver 门先立；一次性 major 版本；`prelude` 提供迁移路径 |
| P4 文档门拖慢 CI | doc 门与 test 门并行；只对公开 crate `-D warnings` |
| P0 裁决反复 | 写成 ADR，冻结 6 个月；变更需新 ADR |

## 10. 本文件用到的工具证据

| 结论 | 证据 |
|---|---|
| 门面 20 条模块再导出 / 241 声明项 | `cargo public-api -p reed -s`（1120 行、32 `pub use`） |
| 8 个被再导出 crate ≈3,580 声明项 | `cargo public-api -p <crate> -s` 九份快照 |
| semver 只覆盖 3/11 | `cargo semver-checks --workspace -v` 的 7 条 `Skipping` |
| 外部类型 | `cargo +nightly-2026-03-20 check-external-types` 九份表 |
| filter→io 唯一越层 | `cargo modules dependencies --lib -p reed-filter` + `cargo metadata` |
| 2 处模块环 / 0 孤儿 | `cargo modules dependencies/orphans` 九份输出 |
| feature 唯一失败点 | `cargo hack check --workspace --each-feature`（14 plans / 1 fail） |
| 6 个 doctest / rustdoc 5×19 | `cargo test --doc --workspace`、`cargo doc -D warnings` |
| 概念散布 | `codegraph query Region`（2 个 struct）、`query count_alts`（1 个符号但 5 处内联副本） |
| 620 条 inactive_code | `rust-analyzer diagnostics .` |

**codegraph 概念散布核对**（`codegraph query <term>`，索引 274 文件 / 6,187 节点）：

| 查询 | 结果 | 含义 |
|---|---|---|
| `Region` | 2 个 `struct Region`（reed-io、reed-filter）+ `MsaRegion` | 三份区域类型 |
| `count_alts` | 1 个函数符号 | 声明唯一；5 处内联副本不在符号表里，需 `rg` 补 |
| `GtClass` | 1 个 enum + `GtClassRule` | 语义有类型，但 convert/pca 未都用 |
| `TempSpill` | 1 个 struct + `StdinSpill` | 两种 spill，均只拥有路径 |
| `ThreadBudget` | 1 个 struct + 7 个方法 | 预算有类型，但 8 处未使用它 |

结论：codegraph 能确认"同名类型/符号"的散布与集中；对"内联复制的规则"仍需 `rg` 补——两者互补。

> 原始输出见 [STATUS.md](STATUS.md) §7；逐条证据见 [01-review.md](01-review.md) 各部分。
