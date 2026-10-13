# 06 修复路线（P0–P3）

> 原则：**先恢复规范承诺（线程池唯一、顺序路径覆盖辅助阶段），再统一上限语义，最后补纪律防漂移**。
> 本轮的修复对象是 01-review.md §3 的缺陷清单；不重开 v7 的 Z 系列、v6 的 Y 系列结构项。
> 每一步都必须满足 §1 通用门。

---

## 1. 通用验收门（每一步都必须满足）【规】

1. `cargo test --release --workspace` 全绿；`cargo clippy --release --workspace --all-targets` 不新增警告。
2. **每个 P0/P1 缺陷先写出复现它的失败测试，再修复**（红 → 绿），测试名描述缺陷本身。
3. 确定性：threads ∈ {1, 2, 3, 7, 56} 输出 sha256 相同（filter 与 stats 都要测；
   rayon 三阶段的合并顺序保证此性质，本轮不得破坏）。
4. 性能：不回退超过 10%，并在 `bench/suites/e2e/results.tsv` 追加一行。
5. **运行时核占用**：`--threads N` 的进程在整个运行期间（含 rayon 阶段）的
   `/proc/self/status` `Threads` ≤ N+3（reader、reducer、主线程三个余量）；
   wall/user 比 ≤ N+ε。这是本轮新增的门，替代以往"只看输出正确"的验证。
6. 同一 PR 内删除被替代的旧实现，不留"以后再删"的并行代码。
7. 任何新增 `thread::spawn` / `par_iter` / `into_par_iter` 必须在 PR 描述里写明
   budget 接线（取数点 = `ThreadBudget::auxiliary()`）；由 CI grep 守卫。
8. 本目录相应条目标记为已完成并注明 commit；`docs/rounds/v8/STATUS.md` 记录状态。

---

## 2. 步骤

### P0 恢复线程池唯一性（D8-1、D8-2）

| 项 | 内容 | 门 |
|---|---|---|
| D8-1a | `ThreadBudget` 加 `pub fn auxiliary(&self) -> usize`（= `workers()`），rustdoc 写明"rayon 阶段与 per-contig 上限的唯一取数点" | 单元测试：`threads(8).auxiliary() == 6`、`threads(1).auxiliary() == 1` |
| D8-1b | `mid_stage::decide_chunk` 批次、`prune::prune_matrices`、`run::decay_contig` 的签名加 `budget: &ThreadBudget`；`is_sequential()` 时 `par_iter`/`into_par_iter` 换成串行 `iter`/`into_iter`，否则用 `rayon::ThreadPoolBuilder::new().num_threads(budget.auxiliary()).build()` 建**每调用一次**的临时池（`pool.install(...)`），不用全局池 | 复现测试先红：`threads1_rayon_stays_serial`（断言运行期间 `Threads` ≤ 3）；修复后转绿；threads ∈ {1,2,3,7,56} 输出 sha256 相同 |
| D8-1c | 三处调用点（`filter/mod.rs` 两处、stats `run.rs` 一处）传 `budget`；CLI 到 driver 的 budget 转发链不变 | 编译期保证（签名变更）；无行为回退 |
| D8-2 | `StatsSink` 加 `max_ld_jobs: usize`（= `budget.auxiliary()`），`finalize_ld` 在 `ld_jobs.len() >= max_ld_jobs` 时 join 最旧 job（对齐 `matrix.rs:941` 的模式）；`finish` 仍按 contig 顺序 join 剩余 job | 复现测试先红：`ld_sink_caps_in_flight_matrices`（100 contig × 小预算 → 峰值 ≤ `--memory`，`Threads` 有上界）；输出与无 cap 时逐字节相同 |

> P0 完成标志：`--threads 1 --mind`、`--threads 1 --ld-decay`、千级 contig 的 `--ld-decay`
> 三种场景下，进程线程数与峰值内存都有硬上界，且输出与并行模式逐字节相同。

### P1 统一上限语义（D8-3、D8-4）

| 项 | 内容 | 门 |
|---|---|---|
| D8-4 | `matrix.rs:1030` 的 `max_prune: opts.threads.max(1)` 改为从 `budget.auxiliary()` 取；`sink.rs` 的 `max_ld_jobs` 同取数点 | `grep -rn "opts.threads" crates/reed/src/filter/pass0/matrix.rs crates/reed/src/stats/sink.rs` 为空；cap 语义统一 |
| D8-3 | 全量 stats `--ld-decay --threads 8`、filter `--threads 8` 各跑 3 次中位数，记录 wall/user 比与 `Threads` 峰值 | wall/user ≤ 8+ε；`Threads` ≤ 11；数据进 `bench/suites/e2e/results.tsv` |

### P2 防漂移纪律（D8-5）

| 项 | 内容 | 门 |
|---|---|---|
| D8-5a | `docs/spec/architecture.md` §7.3 的"现状（2026-10-10，doc_v6 F6-6）"脚注整段删除（已过时），替换为指向本目录的一行 | spec 无过时脚注；`docs/README.md` 索引加 v8 行 |
| D8-5b | CI 加 grep 守卫：`crates/*/src` 里新增 `thread::spawn`/`par_iter`/`into_par_iter` 的 PR 必须同时改动 `docs/rounds/v8/01-review.md` §1 台账（在台账里登记新点） | 守卫脚本进 `ci/`；对现有 6 个生产点豁免（台账已登记） |

---

## 3. 顺序与依赖

```
P0（D8-1a → D8-1b → D8-1c → D8-2）→ P1（D8-4、D8-3）→ P2（D8-5a、D8-5b）
```

- D8-1a 是纯增量（加方法），先行；
- D8-1b 依赖 D8-1a 的取数点；D8-1c 依赖 D8-1b 的签名；
- D8-2 与 D8-1b 可同 PR 也可分开，但都用 `auxiliary()`；
- D8-4 依赖 D8-1a；D8-3 依赖 D8-1b/D8-2 落地后的新行为；
- D8-5b 的豁免清单以 P0 落地后的台账为准。

## 4. 本轮不做什么

- 不改 `pipeline.rs` 的 worker 结构、票据窗口、FrameChain（合规，不动）；
- 不给普通 gzip 输入做并行 inflate（spec §7.3 已声明只能串行，且会警告建议 bgzip）；
- 不动 D7-19 的 `max_prune` cap 逻辑本身（已修，只统一取数点）；
- 不重开 v6/v7 的任何结构项（Y 系列、Z 系列）。
