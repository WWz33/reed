# 08 大文件职责审查与结构因果（18 组 / 20 文件）

> 来源：2026-11-04 九个只读子代理 × 3 批，覆盖 18 组判定、20 个 ≥670 行文件
> （最后一组 `line_edit.rs`/`summary.rs`/`sink.rs` 合判），HEAD `6fd061d`。
> 本文回答两件事：① 子代理的判定是否属实；② "职责混杂 → 补丁热点 → bug"
> 的因果链是否成立。拆分建议一律**绑定 v25 语义项**，不做独立的美学重构。
> 标记：✓ 复核确认；△ 结论成立但路径/行数需修正；✗ 未成立。

---

## 1. 判定总表（已核对）

| 文件 | 行数 | 判定 | 核心问题 | 复核 |
|---|---|---|---|---|
| `crates/reed/src/convert/job.rs` | 2092 | **乱** | spill 基础设施（:380-709）+ FASTA 解析器（:1521-1566）两个独立单元寄生在 chunk-job 文件；`run_convert_opened` 312 行串装配/校验/渲染多阶段（:1679-1991）；gate 概念 6 处落地 | ✓ |
| `crates/reed/src/filter/mod.rs` | 1025 | **乱** | `run_filter` **686 行单体**（:162-848）：region 解析→spill→策略选择→pass0/mid/prune 编排→**BCF magic 手工字节越层**（:612-618）→输出推断→两条近乎逐字重复的 pass1 收尾段（:717-726 vs :770-790）→报告 | ✓ |
| `crates/reed-filter/src/plan/mod.rs` | 1484 | **乱** | `compile()` **618 行混 8 类工作**（:438-1056）；`tail_reads_record` 守卫镜像两份（:416 vs :977，注释自认 "mirrors"）；`expr_reads_gt_fields`（:195）跨域滞留 | ✓ |
| `crates/reed-io/src/bcf/encode_native.rs` | 1817 | **乱（中度）** | GT 解码器寄居编码文件（:1271-1355，无 `n_alt` 校验）、FORMAT 布局自切分（:439-441，GT 位置第 3 份）、typed 写原语库（:1071-1208）三块异质工作寄居——**V25-01/02 的震中** | ✓ |
| `crates/reed/src/stats/writers/reports.rs` | 1518 | **乱** | 6 工种一锅烩；手写第 5 份计数镜像（:566-583）+ 手抄 `stats_policies()` 字面量（:1372-1375）+ 直方图喂入职责外漏 | ✓（"可见性问题"✗，见 §6） |
| `crates/reed-cli/src/translate/popg.rs` | 860 | **乱** | 翻译层混决策层：fstats 谓词 ×3、invariant 梯子+warning ×2（:224-228 vs :403-410）、timing ×2（:210-211 vs :270-271）、`run_popg_cli` 275 行五合一（:20-295） | ✓ |
| `crates/reed/src/filter/pass0/matrix.rs` | 1199 | **乱（中等）** | `MatrixSink` 一个结构背 5 个 job（计数折叠/流簿记/矩阵生命周期/**预算核算**/**prune 线程池**）；`close_pending` 46 行（:982-1028）缠 3 种 job + **同函数内两个预算公式分支**（`ResidentBudget::admit` :1002 vs 裸比较 :1006） | ✓ |
| `crates/reed/src/stats/request.rs` | 937 | **乱（中度）** | 自称 "Request types" 实混 5 类：类型/默认策略/**CLI 文本解析**（`Levels::parse` :193）/**分组校验**（`assign_groups` :852）/**带 jackknife 的计算逻辑**（`FStatAcc` :625） | ✓（路径修正，见 §6） |
| `crates/reed-cli/src/translate/filter.rs` | 791 | 中度乱 | 混杂集中在 `run_filter_cli` 149 行（:358-507）：`-O` 三份之一 + `-W` 镜像 + 第 8 形状 budget（`available_parallelism().expect` :386-389）——V25-06 首要受益者 | ✓ |
| `crates/reed/src/rename/mod.rs` | 1430 | 乱（轻） | 388 行装配函数（`run_rename_opened` :386-774；header/samples/provenance/per-record plan/write 五段）；**互斥校验双镜像已漂移**：`validate` :774 缺 `set_id`+`set_missing_var_ids`（plan 侧 :603）与 `idspace_to`（plan 侧 :633）两条 | ✓ |
| `crates/reed-engine/src/pipeline.rs` | 898 | **不乱** | 但 :411-418 与 :730-736 确认**逐字手抄副本**；双路径保留是 v9 有记录的决策（有对拍门） | ✓ |
| `crates/reed-io/src/vcf/genotype_scan.rs` | 865 | **不乱** | `parse_int_sum` vs BCF `int_sum_column` 是**构造性共享**（`bcf/columns.rs:23` 复用 `INT_MISSING`+`clamp_int_sum`，有差分测试）——V25-01 落点就在此文件，勿另起第三份 | ✓ |
| `crates/reed/src/filter/record_pipeline.rs` | 800 | **不乱**（正面） | "策略差异收敛进 `PipelineOptions`，阶段零分支" 成立；仅 `run_stages` 253 行/19 参数（`:429`，`#[allow(clippy::too_many_arguments)]`）待拆 | ✓ |
| `crates/reed/src/stats/run.rs` | 800 | **不乱** | 单一 job 职责；仅行合法性双写（`:621`/`:681`，V25-03 落点）+ `stats_policies` 被 `reports.rs` 手抄（V25-07 落点） | ✓ |
| `crates/reed-filter/src/expr/mod.rs` | 1050 | **不乱** | 词法已外置（lexer/parser 子模块）；四个语义单点候选归 V25-04，**不要按职责拆文件** | ✓ |
| `crates/reed-cli/src/cli/mod.rs` | 1185 | **不乱** | "no translation logic" 声明成立；仅 `parse_optional_size` 校验后丢弃结果（V25-06 顺带删） | ✓ |
| `crates/reed/src/index.rs` | 819 | **不乱** | 两个嫌疑（tbi/csi 字节、区域查询混入）经实证**均不成立**——字节构造在 reed-io；`translate/index.rs` 是另一份 84 行薄封装 | ✓（路径修正） |
| `crates/reed-io/src/vcf/line_edit.rs` | 763 | **不乱** | 仅 `format_g6` 家族错置（顺手归位） | ✓ |
| `crates/reed/src/stats/writers/summary.rs` | 741 | **不乱** | — | ✓ |
| `crates/reed/src/stats/sink.rs` | 688 | **不乱** | LD 调度约 40% 可轻拆（非本轮） | ✓ |

> 行数为 `wc -l` 实测；函数跨度为"起始行到下一个顶层函数"的差值，与子代理
> 报的 683/612/312/268 差几行是计数口径（是否含签名/空行），结论不变。

---

## 2. 因果链：混杂是补丁的繁殖地

子代理的三条横贯结论里，第一条是核心：**补丁密度与混杂区重合**。把它展开成
可验证的因果链：

```text
职责混杂（文件/函数把多个阶段装在一起）
  → 同一语义在多个阶段各自实现/各自传参（语义多入口）
  → 修复只能打在调用点（guard/assert/注释），没有唯一所有者
  → 下一次改动只改一处，其余入口漂移
  → bug + 补丁密度上升（01 §3 的 9 组反复模式）
  → 再叠一层补丁，结构更混杂（正反馈）
```

这条链在 v25 的缺陷清单里有直接证据：

| 混杂点 | 语义多入口 | 补丁/缺陷证据 |
|---|---|---|
| `filter/mod.rs run_filter` 686 行 | GT 位置/掩码/策略三线在 pass1 内交织；两条 pass1 收尾段重复 | V25-02 要把 GT 位置 6→1；F4 的 head/tail 路由分叉在 plan+head 两处 |
| `encode_native.rs` GT 解码寄居 | GT 解码、FORMAT 切分、typed 编码三块同文件 | V25-01/02 震中；`push_gt_code` 无 `n_alt`；FORMAT 自切分是第 3 份 GT 位置 |
| `convert/job.rs` spill 寄生 | spill/FASTA/多阶段编排同文件 | SEC-1 的 O_EXCL 只修在 stdin 一份，convert spill 仍是可预测名（F9） |
| `plan/mod.rs compile` 618 行 | `tail_reads_record` 两份镜像；`expr_reads_gt_fields` 跨域 | F4：head 路由硬编码 `info_trusted=true`，同一量两种拼写两种队列 |
| `reports.rs` 6 工种 | 第 5 份计数镜像 + 手抄 policy 字面量 | V25-07：计数 5→1；policy 三处镜像 |
| `matrix.rs close_pending` | 矩阵生命周期/预算/prune 三 job 同函数 | V23-10：十 contig 各自"刚好不超"→总量超；两个预算公式分支 |
| `request.rs` 5 类混合 | CLI 文本解析/分组校验/计算逻辑同文件 | R11 退出码翻转；V25-06 参数层收口 |
| `popg.rs` 翻译混决策 | fstats 谓词/invariant/timing 各两份 | R8 裸路由绕过 `-o`/`--pop` 契约 |
| `rename/mod.rs` 装配+校验 | 校验双镜像，`validate` 少两条 | R11/R14 同类 CLI 早期校验缺口 |

**反事实检验**：把混杂点去掉，缺陷还能不能存在？

- GT 解码若只在 `bcf/gt_decode.rs` 一处并带 `n_alt`，就不会有 6 份位置实现、
  `DP:GT` 三种结果和 `0/9` 的静默 NaN——F1/F2 是"多入口"的必然，不是偶发；
- spill 若是共享的 `TempSpill`，SEC-1 的 O_EXCL 就不会只修在 stdin 一份，
  convert 的可预测名也不会存活——F9 是"两份实现"的必然；
- 计数若只有 `SiteCounts::from_codes`，`reports.rs` 的第 5 份镜像就不会漂移——
  V25-07 是"五份镜像"的必然；
- 路由若只有一棵语义树，`head.rs` 的 `info_trusted=true` 就没有落点——F4 是
  "两棵树"的必然。

这些不是"碰巧写错"，是结构允许它们被写错。**结论**：用户"很多 bug 热点也是
结构混乱造成"的判断**成立**，但因果方向要说准——文件混杂本身不直接产生 bug；
它产生的是**语义多入口**，多入口让补丁只能打在调用点，才产生 bug。所以：

- 只拆文件、不收语义 = "把乱搬了个家"，下一轮补丁仍会落在新的调用点；
- 只收语义、不拆文件 = 语义收口后仍可能被下一个功能重新搅浑；
- **正确顺序：语义下沉（V25-01…07）→ 文件拆分作为下沉的载体 → 常驻门锁住。**

---

## 3. 拆分优先级与 v25 绑定

| 优先 | 拆分 | 绑定 | 时机 | 验收 |
|---|---|---|---|---|
| **①** | `encode_native.rs` → `bcf/gt_decode.rs` + `bcf/typed.rs`（~320 行） | **V25-01/02** | 步 1–2 内：先纯搬家，再改域校验/GT 位置 | GT 代码只在新模块；6→1；`DP:GT`/`0/9` fixture 全绿 |
| **②** | `filter/mod.rs` → `strategy.rs`/`sample_subset.rs`/`output.rs`，pass1 合一（~300 行） | **V25-02/03** | 步 2 之后、步 3 之内 | pass1 只剩一条收尾段；BCF magic 下沉 `reed_io::bcf`；输出字节不变 |
| **③** | `convert/job.rs` → `convert/spill.rs`（**V25-09 TempSpill**）+ `reference.rs` + `resolve_targets`（~400 行） | **V25-09** | 步 9 | spill 用共享 `TempSpill`；攻击测试移植；输出字节不变 |
| ④ | `plan/mod.rs` 的 `compile()` 按 section 拆 4 个私有子模块（纯搬家） | V25-04 之后 | 步 4 后 | **已完成**（`7d235ab`）：`plan/compile/{head,mask,tail,finish}.rs`；语义不动 |
| ⑤ | `reports.rs`：计数收口（V25-07a）+ `stats_policies` 收口（V25-07）+ hist 喂入归位 `histogram.rs` | **V25-07** | 步 7 | 计数 5→1；policy 字面量消失；报告数值不变 |
| ⑥ | `request.rs` → `groups.rs`/`kinds.rs`；`FStatAcc` 下沉 reed-stats | v26 | — | 独立，随时可做 |
| ⑦ | `matrix.rs`：`ResidentBudget` 迁 reed-engine `budget`（**v26 `budget::admit` 的第一批**）+ `PrunePool` | v26 | — | 预算单一执行点；测试随迁 |
| ⑧ | `popg.rs` 谓词/invariant/timing 三归一（文件拆分留 v26） | **V25-06/08** | 步 6/8 顺带 | 同一条 timing 渲染；invariant 单点 |

**规则**：

1. **先纯搬家，再改语义**；同一 commit 不同时做"移动 + 行为变化"（否则 review
   无法区分回归来自哪一半）；
2. 每个拆分必须有 **no-output-change 门**（既有 golden/e2e 字节断言）和一条
   新结构守卫；
3. 拆分只在"它是语义下沉的载体"时进入本轮；否则列 v26，不单开"美化轮"；
4. 拆分后仍保留 00 §6 的四问：旧 bug 还能写出来吗、旧补丁还有落点吗、常驻门
   在哪、保留的重复有等价性测试吗。

---

## 4. 明确不拆

| 文件 | 理由 |
|---|---|
| `pipeline.rs` 双路径 | v9 有记录的决策（阻塞点/唤醒语义不同），有对拍门；只删 :411-418 与 :730-736 的逐字副本 |
| `genotype_scan.rs` | 三区段清晰；`parse_int_sum`/`int_sum_column` 已是构造性共享，V25-01 在此加域校验即可，勿另起第三份 |
| `expr/mod.rs` | 词法已外置；四个语义单点候选归 V25-04，按职责拆文件会把语义树拆散 |
| `cli/mod.rs` | "no translation logic" 成立；只删 `parse_optional_size` 的丢弃结果 |
| `index.rs` | 两个嫌疑实证不成立；字节构造在 reed-io |
| `record_pipeline.rs` | **正面对照**：策略差异收敛进 `PipelineOptions`，阶段零分支；只拆 `run_stages` 的参数（后续可做） |
| `line_edit.rs` / `summary.rs` / `sink.rs` | 基本不乱；只做错置函数归位/轻拆 |

**模板**：`record_pipeline.rs` 的 `PipelineOptions`（策略描述符 + 阶段零分支）、
`genotype_scan.rs`↔`bcf/columns.rs` 的构造性共享、`sink.rs`→`writers` 的格式化
委托——这是"不乱"文件共同的健康形态，拆分时照抄。

---

## 5. 与 v25 条目的融合

- **V25-01/02**：`encode_native.rs` 的 GT 解码搬进 `bcf/gt_decode.rs`，域校验、
  GT 位置、语法解码三事同居；这是防止 GT 代码第四次搬家的唯一时机。
- **V25-03**：`filter/mod.rs` 的 pass1 合一放在 `RecordClass` 收口之后；
  BCF magic 字节下沉 `reed_io::bcf`，`filter` 不再手写文件头。
- **V25-04**：`plan/mod.rs` 的 `compile()` 拆分放在语义收口之后，纯搬家；
  `expr_reads_gt_fields` 的归属由 V25-04 决定，拆分 PR 不先动语义。
- **V25-07**：`reports.rs` 的计数/policy 收口优先于文件拆分；计数第 5 份镜像
  必须先消失。
- **V25-09**：`convert/job.rs` 的 spill 先搬进 `convert/spill.rs` 并接
  `TempSpill`，再做加固——一次 PR 完成"搬家 + 共享 + 攻击测试"。
- **v26**：`request.rs`、`matrix.rs` 的拆分，以及 `popg.rs` 的**文件拆分**（其三处
  语义归一随 V25-06/08）不阻塞正确性主线；`matrix.rs` 的 `ResidentBudget` 是
  v26 `budget::admit` 的第一批迁移面。

---

## 6. 复核记录：哪些主张被修正

| 子代理主张 | 实测 | 处置 |
|---|---|---|
| `request.rs` 937 行 | `crates/reed/src/stats/request.rs` 937；`crates/reed-stats/src/request.rs` 是另一份 567 行 | 结论成立，路径写清 |
| `index.rs` 819 行 | `crates/reed/src/index.rs` 819；`translate/index.rs` 是 84 行薄封装 | 结论成立，路径写清 |
| `summary.rs` / `sink.rs` | `crates/reed/src/stats/writers/summary.rs` 741；`crates/reed/src/stats/sink.rs` 688 | 补全路径 |
| `stats_policies` "可见性问题" | **不成立**：`stats/mod.rs` 是 `pub mod run`，`stats_policies` 为 `pub(super)`，同 crate 的 `writers/reports.rs` 可达 | 改为"手抄字面量"（V25-07 收口） |
| `run_filter` 683 / `compile` 612 / `run_convert_opened` 312 / `run_popg_cli` 268 | 实测跨度 686 / 618 / 312 / 275（含签名/空行口径） | 结论不变，数字按实测 |
| 两条 pass1 分支"30 行逐字重复" | 是**近乎逐字**的收尾段（`open_pass1`/`run`/timing/`finish`/slot），不是逐字节全同 | 改为"近乎逐字" |
| `matrix.rs` 预算两公式 | `ResidentBudget::admit`（:1002）与裸比较（:1006）确为同函数两分支 | 确认 |
| `rename` 校验双镜像 | `validate` 缺 `set_id`+`set_missing_var_ids`、`idspace_to` 两条 plan 侧检查 | 确认，列入 V25-08 前置 |

**未独立复核**：`index.rs` 的 tbi/csi 与区域查询细节（子代理称在 reed-io）。
`compile` 的"8 类工作"、`run_popg_cli` 的"五合一"、`reports.rs` 的"6 工种"
是职责归纳，不是可数常量；验收以条目（计数 5→1、GT 6→1 等）为准，不拿这些
数字做门。两者不影响本轮结论。
