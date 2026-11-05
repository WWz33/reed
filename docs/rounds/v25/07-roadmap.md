# 07 顺序、验收门与待拍板

## 1. 顺序与验收门

| 步 | 条目 | 前置 | 验收门 |
|---|---|---|---|
| 0 | [00 补丁下沉审计](00-patch-sinking.md) | — | 每项填 "retired patch" 与 "sink target"；00 §6 四问有答案 |
| 1 | V25-01 GT 域 | failing fixtures 先行 | 四条路径（scan / encode / convert / pca）越界报错；oracle golden；全量绿（945 基线不降） |
| 2 | V25-02 GT 位置 | 步 1（域检查同居一处） | `DP:GT` fixture convert/filter/`-Ob`/pca（body+seam）一致；实现点 **6→1**；GT 解码搬入 `bcf/gt_decode.rs` |
| 3 | V25-03 RecordClass | — | 跨 seam 空行 e2e；CRLF 换行落在 chunk 边界；blanket-true 分支消失；策略表逐格单测；`filter/mod` pass1 合一 |
| 4 | V25-04 filter 语义 | R3 failing 复现 | INFO/AC 对拍 bcftools；三策略（single-pass/sketch/tail）存活位点集一致；实现点 3→1、2→1；R2 round-trip；S58 用例不再自相矛盾 |
| 5 | V25-05 merge_zip_eq | — | 15 处收敛；形状不匹配 panic 带定位；proptest 绿 |
| 6 | V25-06 params | — | usage 矩阵全绿；budget 构造 7→1；`.expect` 清零；常驻 grep 门 |
| 7 | V25-07 core 契约 | 步 3/4 后（消费者先稳） | 计数 5→1；restrict 报错；word kernel 等价性 proptest；counter 三处镜像→1；`reports.rs` 计数/policy 收口 |
| 8 | V25-08 facade | — | 退出码 e2e 表（含 R11）；timing 键名唯一；Cargo.toml 收口；R14 两列 |
| 9 | V25-09 TempSpill | — | 攻击测试移植通过；裸路径拼接清零；convert spill 搬入 `convert/spill.rs` |
| 10 | V25-10 Lanczos | V24-02 协议 | 轮数对比入 bench 日志；非收敛告警存在；STREAM_DEBUG 清零；R12 死文件删除 |
| 11 | `ci/check-structure.sh` 常驻门 | 步 1–10（守卫随条目落地，最后汇总） | 5 类 grep 门全绿并接入 CI |

- 每步独立提交（沿 v24 惯例：`fix(...)` / `refactor(...)` / `perf(...)` +
  `doc_v25` 编号），STATUS 随步更新；
- 步 1–4 是正确性主线，建议连续完成后再动 5–10 的结构项；
- 5 与 6 相互独立，可并行。

### 1.1 结构轨道（与语义主线绑定）

18 组判定（20 个文件）的判定与因果见 [08-large-files.md](08-large-files.md)。拆分只作为
语义下沉的载体，不单开"美化轮"：

| 拆分 | 绑定 | 时机 | 验收 |
|---|---|---|---|
| `encode_native.rs` → `bcf/gt_decode.rs` + `bcf/typed.rs` | V25-01/02 | 步 1–2 | GT 解码单点；6→1；`DP:GT`/`0/9` 全绿 |
| `filter/mod.rs` → pass1 合一 + BCF magic 下沉 `reed_io::bcf` | V25-02/03 | 步 2–3 | pass1 一条收尾段；`filter` 无 `BCF\x02\x02` 字面量；输出字节不变 |
| `convert/job.rs` spill → `convert/spill.rs` | V25-09 | 步 9 | `TempSpill` 共享；攻击测试移植 |
| `reports.rs` 计数/policy 收口 | V25-07 | 步 7 | 计数 5→1；policy 字面量消失 |
| `plan/mod.rs compile` 纯搬家 | V25-04 后 | 步 4 后 | 语义不动；`tail_reads_record` 单点 |
| `request.rs` / `matrix.rs` / `popg.rs` | v26 | — | 08 §3（`matrix.rs` 的 `ResidentBudget` 是 v26 `budget::admit` 第一批） |

三条纪律：

1. **先纯搬家，再改语义**；同一 commit 不同时做"移动 + 行为变化"；
2. 每个拆分必须有 **no-output-change 门**（既有 golden/e2e 字节断言）和一条
   新结构守卫（`ci/check-structure.sh`）；
3. 拆分后仍按 00 §6 四问验收；不满足就退回"只收语义、不拆文件"。

**明确不拆**：`pipeline.rs` 双路径（v9 决策）、`genotype_scan.rs`、
`expr/mod.rs`、`cli/mod.rs`、`index.rs`、`record_pipeline.rs`（正面对照）。

## 2. 待拍板（本轮明确不做，v26 候选）

R2/R11/R12/R14 已从"未归属"补入本轮（V25-04/08/10），不再列 v26；R10 本轮
拍板（03 §5）。下表其余项仍为本轮明确不做。

| 项 | 证据 | 规模 / 风险 |
|---|---|---|
| stats `site_terms` 层 | Hudson ×2、Tajima 常数 ×2、jackknife SE ×4、π/dxy 骨架孪生【述】 | 最大；纯公式唯一化，无正确性缺口；对拍驱动（沿 V24-05 模式） |
| distance 骨架 | `decode_dose` ×6、GRM ×3（`grm.rs` 的 "single implementation" 注释名存实亡）、`MissingPolicy` ×2【述】 | 中；动 rayon 结构需重跑线程不变性测试 |
| kNN 比较子 | `total_cmp` 修复史两轮 6+ 处、isomap ×2、LLE k 钳制【述】 | 小–中；抽 `knn_indices` + 统一比较子 |
| R8 popg 裸路由契约 | 零统计 flag 走 bare route，绕过 `-o`/`--pop` 契约；`-o ""` 退出码与 stats 不同【述】 | 小；需先定"裸路由是否保留" |
| R9 LD pruning 账目 | 同 POS 双记录：pruned 计数与 n_output 矛盾；乘法判定 vs 除法 r² 边界不同【述】 | 中；对拍 plink2 |
| R10 mask 归因顺序 | single-pass（表达式先）与 sketch（声明序）不一致，破坏 "byte-identical funnels" 契约【述】 | **本轮拍板**（03 §5）：A 并入 V25-04（荐，限单文件改动）/ B 保持现状并加测试，转 v26 |
| R13 输出原子性 | 全仓 `fs::rename` 0 处；长跑中断留半文件【述】 | 中；`write_report` tmp+rename，6+ 处收敛 |
| 内存预算集中核算 | v20 S76 → v21 M11 → v23 V23-10，6 处检查 3 种公式 3 种超限行为【述】 | 中；`budget::admit` 唯一执行点（`matrix.rs` 的 `ResidentBudget` 是第一批迁移面） |
| 大文件剩余拆分 | `request.rs` 5 类混合、`matrix.rs` `MatrixSink` 5 job、`popg.rs` 翻译混决策（08 §3 ⑥–⑧） | 中；不阻塞正确性主线；按 08 的顺序做 |
| 钉子户五项 | CI-1 data-gates 不构建被测二进制（历轮"已修复"结论的验证装置本身不可信）；S28 `N_MISSING` 口径分裂；`--hwe midp` 永远拒绝；reed-reduce 已知答案测试为空；MSRV 1.91 vs toolchain 1.99 无 CI 验证【述】 | 各自独立；CI-1 建议最先（它验证其他一切） |

## 3. 流程纪律（沿 v24）

1. **每条正确性项先写 failing 测试再修**（契约先行，沿用 V24-01）；动手前先填
   00 §4 的 "retired patch" 与 "sink target"，否则不知道该收回哪一层；
2. 不改统计口径与输出格式；实施中发现必须改用户可见数值 → 停下，STATUS 记录，
   拍板后再动；
3. 性能结论只认 V24-02 协议（warmup、交替顺序、≥5 次中位数、JSON 输出）；
4. 不把 v25 与 v24 的审计记录混写；v24 保持已提交状态；
5. 【述】行号实施时以 grep 重新锚定——子代理在 v24 落地窗口（01:09–02:55）
   内读取，可能有偏移；
6. 每个条目完成的自检问句：**"这个补丁类下次还会再发生吗？"**——还会，就没
   修到底（v22 R22-3 `LineRef::shape()` 是标杆：一次收口消 12 个散点条目）；
7. 【规】**验收 grep 固化为常驻 CI 门**：沿 `ci/check-doc-drift.sh` 模式建
   `ci/check-structure.sh`，至少五类守卫：merge 裸 zip、GT 位置自写查找、
   mask 的 `unwrap_or(true)`、core 枚举 `_` 兜底、`available_parallelism`
   绕过 budget / stats 手写 `total_cmp`。随步 1–10 落地，步 11 汇总接入 CI——
   一次性 grep 会衰减，常驻门不会；模式 1 的教训正是"按点名修、不按模式扫"。
