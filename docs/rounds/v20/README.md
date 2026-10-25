# reed doc_v20：全库评审 + must-fix 修复（2026-10-24）

> 本目录由两部分组成：对整个 workspace（11 个 crate，约 8.2 万行）的**只读评审**，
> 以及对评审 14 条 must-fix 的**先复现后修复**（测试名 `r1_*` … `r14_*`）。
> 评审由 11 个并行评审完成（每 crate 一个 + 跨 crate 接缝一个），完整报告存于
> 仓库外 `../../../reed-rs-review.md`；本目录只收录**逐行核实过**的结论。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 评审方法与总评；R1–R14 逐条：`file:line`、严重度、修复前失败测试的实测输出、修法；评审结论中被证实有出入的三处修正；核实为干净的部分 |
| [STATUS.md](STATUS.md) | 实施状态：14 项修复的 项/内容/验证 表（commit 列 = 未提交，与 doc_v18 WIP 同树）；被替换的两处旧锚点及其理由；偏差记录；全量验收 |
| [06-roadmap.md](06-roadmap.md) | 后续路线：评审的 ~66 条 should-fix 与 lean 清单按 P0/P1/P2 组织成 v21 候选，含每类的验收门 |

## 与 doc_v15/v17/v18 的关系

- doc_v15 也是"评审 + 修复"轮（40 findings，`01-review.md` + `02-architecture.md`），
  本轮沿用它的**先写失败测试再修复**纪律（v2 06 §1.2 起的规矩）。
- 不改架构、不改 CLI 形态、不动统计定义；doc_v18 的 RNG/xval/aboot WIP
  （工作树中未提交）不受影响，v20 修复与其共存于同一棵树。

## 结论速览

| 维度 | 结论 |
|---|---|
| 评审覆盖 | 11 个只读评审：reed-core/io/stats/reduce/engine/filter/distance/graph/oracle、reed 主 crate（输入层 + 命令层各一）、reed-cli + 工作区配置、跨 crate 接缝 |
| 公式正确性【据】 | 经典统计量手工验算全部正确：WC84 pair/multi（手算到精确分数）、Hudson、Tajima e1/e2、Watterson、pixy 补全、Patterson f2/f3/f4、Wigginton HWE、UPGMA/NJ、oracle 判定为真独立（零依赖、内嵌 scikit-allel 1.3.13 黄金值） |
| 确定性【据】 | rayon 归约顺序全部按索引固定；e2e 在 {1,2,3,7,56} 线程下字节一致；本轮 R13 的按行并行延续该契约（1 vs 4 线程 `to_bits()` 相等） |
| must-fix | 14 条全部属实，全部先复现（修复前失败测试）后修复；详见 [01-review.md](01-review.md) |
| 第二批（bug 类 should-fix） | 应用户指令**同轮追加**（不开新轮）：06-roadmap 中错误结果/崩溃/安全类 16 条（S28–S32/S36/S37/S39–S43/S46/S58/S68/S76/S79/SEC-1）同纪律修复完毕，详见 STATUS.md「第二批」 |
| 第三批（复审修复） | 两批完成后做只读复审：**R8 把正确的 t-SNE 梯度改错了**（对 KL 代价做中心差分定案），**S36 关于 ape 的断言与实测相反**（本机 ape 5.8.1 保留根节点），S42 漏了默认写入路径、S76 检查晚于分配、S58 让 `--set-gts` 的 `INFO/AC` 静默失效。16 条已在同轮串行修完，详见 STATUS.md「第三批」 |
| 修复后基线【据】 | 三批之后 `cargo test --workspace` exit 0（全部目标 0 失败，含 noodles 差分套件与线程确定性套件）；`cargo clippy --workspace --all-targets` 0 告警；`ci/check-doc-drift.sh` 通过 |
| should-fix 余量 | 验证面/CI 门/lean 等非 bug 类条目仍留在 [06-roadmap.md](06-roadmap.md)（已修条目画线标记） |

## 决策状态

| 决策 | 状态 |
|---|---|
| 14 条 must-fix 的修法 | 已定：最小改动；每条一个 `rN_*` 失败测试先行 |
| R9 单态列语义 | 已定：中心化 0（plink2/GCTA），与 flashpca 均值填充在 fixed-ALT 上**有意分歧**；CHANGELOG 已注明 |
| R10 stream + pairwise | 已定：拒绝而非静默降级（含 auto→stream 路径）；一处旧 e2e 锚定的正是该 bug，已改锚并注明 |
| R11 `--set-missing-var-ids` + `--update-name` | 已定：核实为合法组合（map 对具名 ID 仍可达），**不拒绝**，用测试钉住 |
| R13 范围 | 已定：本轮只做按行并行 + dosage 行预计算（数值逐位不变）；流式分块与 `--threads` 全链路接入留 backlog |
| 提交切分 | 未决（用户）：v20 改动与 doc_v18 WIP 同树未提交；按 STATUS 表所列文件 + `r*_` 测试名 + `doc_v20` 注释切分 |

## 标记约定

【规】必须；【荐】推荐；【决】需要拍板；【据】有实测或失败测试输出作证；
**P0/P1/P2** 为优先级。R1–R14 为本轮 must-fix 编号（与仓库外评审报告一致）。
