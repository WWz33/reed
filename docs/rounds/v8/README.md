# reed doc_v8：线程模型复审 + 并行预算统一规范（2026-10-13）

> 本目录是 `docs/rounds/v7/` 之后的**第八轮规范**。它由两部分组成：
> 对当前实现（v7 修复全部落地后）**线程模型与并行预算**的独立复审，
> 以及在复审基础上的**并行预算统一修复规范**。
> 本目录不修改任何代码，所有改动以规范形式给出，由实现者执行。
> 历史轮次 `docs/rounds/vN/` 只增不改；当前规范索引见 `docs/README.md`。

## 背景

`docs/spec/architecture.md` §7.3 规定：进程里**不允许出现第二个线程池**，rayon 只在
pass 之间的批量运算里使用，并且**复用同一个预算**。本轮以 codegraph 符号索引
（157 个 src 文件，2026-10-13 同步）+ 全仓库 grep 双向核对，枚举了全部线程创建点，
逐条对照该规范。结论：流水线本体合规，但**三个 rayon 辅助阶段与 stats LD sink
完全脱离 `ThreadBudget`**，规范承诺与实现之间存在系统性缺口，且此前任何一轮
（v2–v7）都未登记。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 线程模型复审：全部 7 个线程创建点台账、与规范逐条对照、缺陷清单 D8-1…D8-5（带 `file:line`、严重度、修法）、已降级/已消解的历史说法 |
| [06-roadmap.md](06-roadmap.md) | 修复路线 P0–P3：通用验收门、每步内容与门、确定性纪律 |
| [STATUS.md](STATUS.md) | 实施状态（起始全部未修） |

## 复审方法【据】

1. codegraph 符号索引全量枚举 `pipeline::run` 调用者与并行相关符号；
2. 全仓库 grep 七个模式：`thread::spawn` / `thread::scope` / `thread::Builder` /
   `par_iter` / `into_par_iter` / `par_bridge` / `rayon::join` / `rayon::scope`；
3. 反查不存在项：`ThreadPoolBuilder` / `rayon::ThreadPool` / noodles MT reader-writer /
   tokio / async-std / futures / `std::process::Command`（均不存在，见 01 §2）；
4. 逐函数检查签名是否接收 `ThreadBudget`、调用点是否传线程数。

## 结论速览

| 维度 | 结论 |
|---|---|
| 线程创建点 | 生产代码 **6 处** + 测试 1 处，全部登记在 01 §1；无隐藏线程池 |
| 流水线本体 | **合规**。1 reader + 1 reducer + N−2 workers，`threads<=2` 走顺序路径（F6-6 已落地，`is_sequential` 有调用者） |
| rayon 三阶段 | **失控**（D8-1）。`mid_stage` / `prune_matrices` / `decay_contig` 不接收 budget，用满 rayon 全局池（= 全部核），与 spec §7.3 直接冲突；`--threads 1` 时进程仍多线程 |
| stats LD sink | **无上限**（D8-2）。每 contig 一个 `thread::spawn`，`finish` 才 join，M 个 contig → M 个线程 + M 份常驻矩阵；filter pass0 的同一模式（D7-19）已修，此处未修 |
| 跨阶段超订 | **可能**（D8-3）。流水线 workers 与 rayon 全局池同时活跃，核数可被翻倍占用 |
| reader | 串行但**已非瓶颈**（mmap 后 2.2 → 0.35 s，v4 W3）；上轮"reader 与 gzip 并列"的权重作废 |
| 历史修正 | D7-19 **已修**（`max_prune` cap 在 `matrix.rs:941`），此前综述误引为未修；本目录澄清 |

## 决策状态

| 决策 | 状态 |
|---|---|
| rayon 必须复用 `ThreadBudget`，进程只有一个线程池 | 已定（spec §7.3，v1 起有效；本轮首次给出落地形态） |
| 顺序路径（threads≤2）必须覆盖辅助阶段 | 已定（spec §7.3 表格脚注的原有语义；本轮补实现规范） |
| stats LD sink 加 in-flight 上限 | 已定（对齐 D7-19 的修法） |

## 标记约定

【规】必须；【荐】推荐；【据】有实测或源码证据；**P0/P1/P2** 为修复优先级。
