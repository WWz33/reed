# 00 补丁下沉审计（本轮组织原则）

> 起点 `6fd061d`。本轮的每个条目不是"再加一个修复"，而是"把一类补丁从调用点
> 收回底层"。判据：把 236 个 commit 里的 63 个 fix 当作症状样本，逐条回答——
> ① 补丁最初打在哪一层？② 同类语义还有几份实现/几次修复？③ 只留下沉后的唯一
> 入口，这个补丁还能不能写出来？能写出来 = 没沉到底；写不出来 = 结构投入到位。

---

## 1. 问题陈述

"功能能跑，但不优雅"的具体形态不是代码风格，而是**补丁落点**：

- 同一个语义有 N 份实现，修一份漏 N−1 份（GT 位置、计数、缺失比较、AC/AF）；
- 修复是在调用点加 guard/assert/注释，而不是改类型或 API（S68 的 mid_stage
  `assert!`、`info_trusted: true` 硬编码）；
- 一个 bug 的修法在 N 个文件重复出现（`-O ''` panic 修 3 遍、
  `available_parallelism().expect` ×6）；
- 不变式靠"人肉同步"维持（counter 清单 core/io/filter 三处镜像；
  GRM 的 "single implementation" 注释已名存实亡）；
- 缺陷是"没有生产者/消费者"的死链路（`TRIG_HAPLOID_MISSING` 自环、
  `reed-reduce/src/linalg.rs` 死文件）。

补丁史已经给出结论：**反复出问题的不是缺补丁，是缺结构**（01 §3）。所以本轮
先做 P0 止血，再按本文件的判据把出血机制拆掉。

### 1.1 结构混杂是补丁的繁殖地

18 组大文件的职责审查（08）给出同一条因果链：

```text
职责混杂（文件/函数把多个阶段装在一起）
  → 同一语义在多个阶段各自实现/各自传参（语义多入口）
  → 修复只能打在调用点（guard/assert/注释），没有唯一所有者
  → 下一次改动只改一处，其余入口漂移
  → bug + 补丁密度上升
  → 再叠一层补丁，结构更混杂（正反馈）
```

所以文件混杂本身不是 bug；它是**语义多入口的土壤**。只拆文件不收语义 =
"把乱搬了个家"；只收语义不拆文件 = 下一个功能仍会把语义搅浑。顺序固定为：
**语义下沉 → 拆分（纯搬家）→ 常驻门**（08 §2–§3）。

---

## 2. 判据：什么算"散落补丁"

一条补丁满足任一条，即进入下沉候选：

| 标记 | 判据 | 本轮实例 |
|---|---|---|
| **P1** | 同一语义在 ≥2 处实现，修一处漏一处 | GT 位置 6 份、计数 5 份、缺失比较 3 份、AC/AF 2 份 |
| **P2** | 修复是在调用点加 guard/assert/注释，而不是改类型/API | S68 的 `assert!`；`from_line(..., true)`；`unwrap_or_else` 保底分组 |
| **P3** | 一个 bug 的修法在 N 个文件重复出现 | `-O ''` 三处；budget 构造七处；argv 捕获五处 |
| **P4** | 不变式靠人肉同步（镜像表/重复注释/手抄惯用法） | counter 三处；`Layout` 与 encode 各自 split FORMAT；kNN 比较子手抄 |
| **P5** | 没有生产者/消费者的死链路，或恒为常量的输出 | `TRIG_HAPLOID_MISSING`；`linalg.rs` 198 行；`max_residual` 被丢弃 |
| **P6** | 同一文件/函数承担 ≥3 种职责，且与补丁热点重合（结构混杂） | `filter/mod.rs run_filter` 686 行；`plan/mod.rs compile` 618 行；`convert/job.rs` 的 spill 寄生 |

**不下沉的重复**（刻意保留，见 01 §4.3）：oracle 的独立 `count_site`、
text/BCF 双 reader、reed-stats 与 reed-distance 的 per-pair 累积、word kernel。
保留条件是**有等价性测试**，不是"注释说它们一致"。

---

## 3. 动作：四选一

| 动作 | 条件 | 落点 / 产物 |
|---|---|---|
| **Sink（下沉）** | 语义有唯一所有者，且下层能表达该不变式 | reed-core（纯类型/语义）、reed-io（记录解码视图）、reed-engine（预算/指标）、reed-cli `params`（参数域）、facade（边界） |
| **Delete（删除）** | 无生产者、无消费者、输出恒为常量 | 幽灵计数器、`linalg.rs`、死分支 |
| **Keep（保留）** | 刻意重复（oracle/性能/双 reader） | 等价性测试 + 在 01 §4.3 登记 |
| **Verify（补验证）** | 代码已对，但验证装置不可信或缺失 | CI 门、KAT、对拍 fixture、常驻 grep |
| **Split（拆分）** | P6：混杂点与补丁热点重合，且拆分是下沉的载体 | 文件/函数边界；纯搬家 + no-output-change 门 + 结构守卫（08 §3） |

每一种动作都必须落到**一个层**和**一个验收门**；只写"重构一下"不算条目。

---

## 4. 本轮补丁族 → 下沉物

| 补丁族 | 反复证据 | 下沉物 | 目标层 | 条目 |
|---|---|---|---|---|
| GT 语义（位置/域/解码） | V22-R2-1 一 commit 修三处；位置 6 份（含 PCA 两份） | `Layout`/`gt_subfield`/`decode_gt` | reed-io | V25-01/02 |
| 行合法性 | 5 入口各自推导；B13-2 只对齐一对 | `RecordClass::classify_line` | reed-io | V25-03 |
| 缺失比较 / AC-AF | V23-59 只修 2/3；head/tail 两份 | `semantics::{compare_missing, info_af}` | reed-filter | V25-04 |
| merge 形状 | V23-28 修 3 处后仍开放 15 处 | `merge_zip_eq` | reed-stats | V25-05 |
| 参数校验 | 4 轮 15+ 补丁；`-O` 3 份、budget 7 份 | `params` newtype + `ThreadBudget::from_cli` | reed-cli | V25-06 |
| 计数/掩码 | 计数 5 份；mask 长度检查 5 处；counter 三处镜像 | `SiteCounts::from_codes`、`KeepMask`、counter 单一定义 | reed-core | V25-07 |
| 门面边界 | 直连 3 crate；退出码/时序/provenance 漂移 | re-export + `Error` 分类 + `StageMetrics::render_tsv` + `Provenance` | facade/core | V25-08 |
| 临时文件 | SEC-1 只修 stdin 一份 | `TempSpill` | core/io | V25-09 |
| Lanczos 观测 | 残差判据 + `STREAM_DEBUG` + `max_residual` 丢弃 | 混合尺度 + 告警 + `tracing` | reed-reduce | V25-10 |
| 结构混杂（P6） | 18 组判定（20 个文件）中 9 组判乱；8/9 乱组与补丁热点重合 | 文件/函数拆分（仅作语义下沉的载体） | 各层 | 08 §3 |
| 未归属小补丁 | R2/R11/R12/R14 | 并入对应条目，见 07 §1 | 见 03/06/04 | 本轮增补 |

**本轮新确认的漏项**：PCA 的 genotype collector 还有两份自写 GT 位置
（`crates/reed/src/genotype.rs:171-174` 的 `run_chunk`、`:247-248` 的
`run_seam`），且 `PcaJob::dosage`（`:106`）把越界 allele 静默变 NaN 而不是
`MalformedRecord`。这意味着 V25-02 原计划的"4→1"实际是 **6→1**，V25-01 的
域检查还多一条 PCA 路径；详见 02 §1–§2。

**本轮不下沉（v26 候选）**：同样满足 P1–P5，但规模或依赖关系不适合与正确性
主线同轮，登记在此防止"看不见就忘"：

| 补丁族 | 证据 | 为什么不在本轮 |
|---|---|---|
| stats `site_terms` | Hudson ×2、Tajima 常数 ×2、jackknife SE ×4、π/dxy 骨架【述】 | 规模最大、无正确性缺口；等 V25-01–04 稳定后单独对拍 |
| distance 骨架 | `decode_dose` ×6、GRM ×3、`MissingPolicy` ×2【述】 | 动 rayon 结构，需独立重跑线程不变性 |
| kNN 比较子 / `LastHit` | `total_cmp` 两轮 6+ 处、两个哨兵缓存各咬一次【述】 | 与 GT/记录语义正交，独立小轮 |
| `budget::admit` | 6 处检查、3 种公式、3 种超限行为【述】 | 与 V25-06 参数层相邻，但改的是执行点，留 v26 |
| R13 输出原子性 | 全仓 `fs::rename` 0 处【据】 | 数据安全项；与 V25-09 的 spill 安全分开做 |
| R8 / R9 | popg 裸路由；LD pruning 账目【述】 | 需先定契约/以 plink2 为准，不能顺手改 |

---

## 5. 下沉层的所有权

| 层 | 拥有的不变式 | 本轮下沉物 |
|---|---|---|
| reed-core | 纯类型/语义：`GtCode`、`SiteCounts`、`KeepMask`、`VariantClass`、`Provenance`、counter 定义 | V25-07、V25-08c |
| reed-io | 记录解码视图：`Layout`、`gt_subfield`、`decode_gt`、`classify_line`、样本列校验 | V25-01/02/03 |
| reed-engine | 预算、指标、中止语义：`thread_budget::from_cli`、`StageMetrics::render_tsv`、`budget::admit`（v26） | V25-06/08b |
| reed-filter | 过滤求值语义：`semantics` 模块 | V25-04 |
| reed-stats | 累加器合并律：`merge_zip_eq`；`site_terms`（v26） | V25-05 |
| reed-cli | 参数域：`params` newtype + value_parser | V25-06 |
| facade (`reed`) | 跨 crate 边界：re-export、`Error` 分类、`run_*_opened` | V25-08 |

规则：**不变式由能强制执行它的最低层拥有**。上层只做策略选择，不重写判别。

---

## 6. 验收：补丁复发测试

每个下沉项必须回答四问（写进对应条目的验收）：

1. **旧 bug 还能写出来吗？** 在唯一入口之外再造一份旧实现，编译器/可见性/CI
   是否阻止？（类型、私有字段、单点 API）
2. **旧补丁还有落点吗？** 原来的修复位置是否已不存在（grep 清零）？
3. **常驻门**：`ci/check-structure.sh` 是否有对应守卫？（GT 自写位置、裸 zip、
   mask `unwrap_or(true)`、core `_` 兜底、`available_parallelism` 绕过 budget、
   stats 手写 `total_cmp` 比较子）
4. **等价性**：保留下来的性能/双 reader 路径是否有与唯一入口的等价性测试？

一条不满足 = 没沉到底（07 §3.6 的自检问句）。

---

## 7. 与正确性的关系

F1–F10 里 8 条的根源都是"语义有多个实现点/契约无类型承载"。修 P0 是止血，
下沉是拆掉出血机制。顺序固定为：

1. 先修 P0（V25-01 → 02 → 03 → 04），failing test 先行；
2. 再收结构（05 → 06 → 07 → 08 → 09 → 10），每条填 "retired patch"；
3. 最后用常驻门锁住（`ci/check-structure.sh` 与各条目的 grep 验收）。

**边界**：本轮不改变统计口径与输出格式；唯一预期的行为变化是"静默错误 →
显式报错"（V25 README 决策表）。下沉若必须改用户可见数值，停下记录 STATUS，
拍板后再动。
