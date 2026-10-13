# 01 线程模型复审

> 复审对象：v7 修复全部落地后的当前实现。
> 方法：codegraph 符号索引（157 个 src 文件，2026-10-13 同步）+ 全仓库 grep 双向核对。
> 基准规范：`docs/spec/architecture.md` §7.3（线程预算）、§7.7（背压与内存预算）。

---

## 1. 线程创建点台账（全量）【据】

生产代码里创建线程的位置共 **6 处**，测试代码 1 处。此外无 `ThreadPoolBuilder`、
无 noodles MT reader/writer、无 tokio/async-std/futures、无进程 fork。

| # | 位置 | 原语 | 用途 | 是否受 `ThreadBudget` 约束 |
|---|---|---|---|---|
| T1 | `reed-engine/src/pipeline.rs:343` | `thread::scope` | 主流水线：1 reader + 1 reducer + N−2 workers | **是**（`budget.workers()`；`threads<=2` 走 `run_sequential`，`pipeline.rs:225`） |
| T2 | `reed/src/filter/mid_stage.rs:198` | rayon `par_iter` | `--mind` 等中间判定，逐 chunk 独立判定 | **否**（函数签名无 budget，走 rayon 全局池） |
| T3 | `reed/src/filter/pass0/prune.rs:187` | rayon `into_par_iter` | LD 剪枝，按 contig 一个 `LdPruner` | **否**（同上） |
| T4 | `reed/src/stats/run.rs:190` | rayon `par_iter` | LD decay，按 4096 位点块扫描 | **否**（`decay_contig` 签名无 budget） |
| T5 | `reed/src/stats/sink.rs:231` | `thread::spawn` | stats 每 contig 一个 decay 线程，`finish` 才 join | **否**（无上限，见 D8-2） |
| T6 | `reed/src/filter/pass0/matrix.rs:950` | `thread::spawn` | filter 流式剪枝每 contig 一个线程 | **部分**（有 `max_prune` cap，`matrix.rs:941`；但 cap 值 = `opts.threads`，与引擎 budget 不共享） |
| T7 | `reed-engine/src/abort.rs:244` | `thread::spawn` | 许可释放看门狗 | —（仅测试代码） |

### 1.1 rayon 三阶段的共同模式

T2/T3/T4 的函数签名与调用点均不接收、不转发任何线程数信息：

- `decide_chunk` 批次（mid_stage）：`chunks.par_iter().map(...)`，调用点 `filter/mod.rs` 不传 budget；
- `prune_matrices`（prune）：`matrices.into_par_iter()`；
- `decay_contig`（run）：`ranges.par_iter()`。

三者都注释自称"结果与线程数无关"（由合并顺序保证），这一性质成立；但**占多少核**
这一资源问题没有任何一层管理。

---

## 2. 与规范逐条对照

### 2.1 spec §7.3：「进程里不允许出现第二个线程池」

> 原文："进程里不允许出现第二个线程池（不用 noodles MT reader/writer；rayon 只在
> pass 之间的批量运算里使用，例如 sketch 的尾部判定，并且复用同一个预算）。"

**判定：未实现。** rayon 全局池（默认 = `available_parallelism`，56 核机器即 56 线程）
在 T2/T3/T4 活跃期间与流水线 workers 并存。规范允许 rayon 存在，但附加了
"复用同一个预算"的条件；该条件在代码里没有任何对应物（全仓库无
`ThreadPoolBuilder` / `rayon::ThreadPool` / `num_threads` 调用）。

### 2.2 spec §7.3：顺序模式承诺

> 原文（表格脚注）："N ≤ 2 时退化为 reader 和 worker 合一的顺序模式。
> 结果必须与并行模式逐字节相同（测试锁定）。"

**判定：部分实现。** 流水线本体合规（`pipeline.rs:225` 的 `is_sequential` 分支，
F6-6 已落地）；但 T2/T3/T4/T5 不看 `is_sequential`，`--threads 1` 时仍然多线程。
"one thread reads, inflates, computes and writes" 对辅助阶段不成立（D8-1）。

### 2.3 spec §7.7 / D7-19：in-flight 内存预算

**判定：filter pass0 已修，stats LD sink 未修。**
`matrix.rs:941` 的 `max_prune` cap 让 M 个 contig 不再驻留 M 份矩阵（D7-19 已修，
测试 `prune_multi_contig_is_thread_count_independent`）；但 `sink.rs:231` 的
`finalize_ld` 没有任何上限，M 个 contig → M 个线程 + M 份常驻 `ContigMatrix`（D8-2）。
两处是同一模式，v7 只修了一处。

---

## 3. 缺陷清单

严重度：**P0** = 资源失控 / 承诺失效；**P1** = 可观测的超订或预算突破；
**P2** = 文档/纪律缺口。

| # | 严重度 | 位置 | 缺陷 | 后果 | 修法 |
|---|---|---|---|---|---|
| D8-1 | **P0** | `reed/src/filter/mid_stage.rs:198`；`reed/src/filter/pass0/prune.rs:187`；`reed/src/stats/run.rs:190`（+ 各自调用点） | rayon 三阶段不接收 `ThreadBudget`，无 `ThreadPoolBuilder`，直接落在 rayon 全局池（= 全部核） | spec §7.3"复用同一个预算"未实现；`--threads 1` 时进程仍 56 线程；顺序模式承诺对辅助阶段失效 | 三阶段签名加 `&ThreadBudget`；`is_sequential()` 时走串行迭代器，否则在调用处用 `rayon::ThreadPoolBuilder::new().num_threads(budget.workers())` 建**每调用一次**的临时池（不用全局池），池线程数与流水线 workers 同预算 |
| D8-2 | **P0** | `reed/src/stats/sink.rs:223-243`（`finalize_ld` 每 contig 立即 `thread::spawn`，`ld_jobs` 到 `finish` 才 join） | 无并发上限；M 个 contig 可同时驻留 M 份 `ContigMatrix` 与 M 个 OS 线程 | `--memory` 预算被突破；contig 数大（千级 scaffold）时线程数与常驻内存同步失控；与 D7-19 同一模式但未被 v7 覆盖 | 对齐 `matrix.rs:941` 的修法：sink 加 `max_ld_jobs` 上限（= `budget.workers()`），达到上限时 join 最旧的 job；在 D8-1 落地后，cap 与 rayon 池共享同一 `workers()` 值 |
| D8-3 | **P1** | T1 与 T2/T3/T4 的组合 | 流水线 workers（N−2）与 rayon 全局池（N）在同一进程内可同时活跃；filter 的 pass0 流水线结束后立即进入 `prune_matrices`，stats 流水线结束后进入 `decay_contig`，时间上虽不重叠，但 mid_stage 的 rayon 判定在流水线**运行期间**发生 | CPU 占用瞬时可达 ~2N；在共享机器/调度器上表现为超订；`--threads N` 的实际核占用不可预期 | D8-1 的临时池把 rayon 阶段也压到 `workers()`；在 06 §1 通用门加"运行时核占用 ≤ `--threads`"的观测项（`/proc/self/status` 的 `Threads` 或 `time -v` 的 wall/user 比） |
| D8-4 | **P1** | `reed/src/filter/pass0/matrix.rs:1030`（`max_prune: opts.threads.max(1)`） | D7-19 的 cap 取的是 CLI 原始线程数，不是 `ThreadBudget::workers()`；与 T5/T6 的上限值三处各自为政 | `--threads 4` 时流水线 workers=2 但 prune cap=4，上限语义不统一；后续改动容易顾此失彼 | 三处上限统一为 `budget.workers()`；在 `reed-engine::budget` 加一个 `pub fn auxiliary(&self) -> usize`（= `workers()`）作为唯一取数点，三处都调它 |
| D8-5 | **P2** | `docs/rounds/v8/README.md` 之前的全部轮次 | "rayon 复用同一预算"自 v1 §7.3 起就是规范，但 v2–v7 的 review 均未核对辅助阶段的 budget 接线；stats LD sink 的无上限 spawn 与 D7-19 同模式却漏检 | 规范与实现的系统性漂移无人值守；同类缺口（同模式只修一处）反复出现 | 本轮登记；06 §1 通用门加"任何新增 `thread::spawn`/`par_iter` 必须在 PR 描述里写明 budget 接线"，由 `ci/check-doc-drift.sh` 同款思路的 grep 守卫 |
| D8-6 | **P0** | `reed-engine/src/budget.rs`（`run_auxiliary` 串行分支） + 三阶段闭包 | **D8-3 实测发现**（2026-10-13）：串行分支虽不建临时池，但闭包内的 `ranges.par_iter()` 仍会**惰性初始化 rayon 全局池**——`--threads 1` 实测 58 线程、cpu_ratio 5.5×，spec §7.3"one thread reads, inflates, computes and writes"名存实亡。strace `clone3` 76 次 + gdb 栈定位 | 顺序承诺完全失效；更广：任何 `run_auxiliary` 串行调用方的 par_iter 都逃逸到全局池 | `run_auxiliary` 闭包签名改 `FnOnce(bool)`，`parallel==false` 时闭包走 `iter()`/`into_iter()`；CI 守卫升级为"par_iter 必须在 `if parallel {` 分支内"（D8-5b 同步收紧） |

---

## 4. 已降级 / 已消解的历史说法（防止再被误引）

| 说法 | 出处 | 本轮核实结论 |
|---|---|---|
| "filter pass0 每 contig 线程无上限（D7-19）" | v7/01-review.md:75 | **已修**。`matrix.rs:941` 有 `max_prune` cap，达到上限 join 最旧 job；v7/STATUS.md:34 标"已修"。继续引用为"未修"是错的 |
| "reader 是串行瓶颈，与 gzip 解压并列" | 此前综述 | **权重作废**。v4/STATUS.md:69 W3 mmap 后全量 stats reader 2.16–2.20 → 0.35–0.36 s（−84%），wall 由 CPU 决定。reader 确实串行，但已不是主导开销 |
| "顺序模式未实现（F6-6），`is_sequential` 无调用者" | v6 时代 spec §7.3 脚注 | **已落地**。`pipeline.rs:225` `if budget.is_sequential()` 有调用者；spec 该脚注的"现状（2026-10-10，doc_v6 F6-6）"一段已过时，应由本目录取代 |
| "stats LD sink 无上限" | 此前综述标"未登记" | **属实且仍未修**，本轮登记为 D8-2 |

---

## 5. 复核门（每条缺陷的验证方式）

| # | 门 |
|---|---|
| D8-1 | `--threads 1 --mind` 与 `--threads 1 --ld-decay` 运行期间 `/proc/self/status` 的 `Threads` ≤ 3（reader+reducer+主线程）；`--threads 8` 时 rayon 阶段核占用 ≤ 6（= `workers()`）；threads ∈ {1,2,3,7,56} 输出 sha256 相同 |
| D8-2 | 100 个 contig × 小预算 → 峰值内存 ≤ `--memory`；运行期间 `Threads` 有上界（≤ `workers()`+3）；输出与无 cap 时逐字节相同（join 顺序 = contig 顺序不变） |
| D8-3 | 全量 stats `--ld-decay --threads 8` 的 wall/user 比 ≤ 8+ε；`mid_stage` 期间 `Threads` ≤ `threads`+3 |
| D8-4 | `grep -n "workers()\|opts.threads" crates/reed/src/filter/pass0/matrix.rs crates/reed/src/stats/sink.rs` 只剩 `budget.auxiliary()` 一个取数点；cap 值语义在 rustdoc 里写明 |
| D8-5 | 本目录存在；`docs/README.md` 索引更新；grep 守卫脚本进 CI |
| D8-6 | `--threads 1 --ld-decay` 运行期间 `/proc` Threads ≤ 3、cpu_ratio ≤ 1.3（实测 2 / 1.21）；`--threads 1 --indep-pairwise` Threads ≤ 2、ratio 1.01；输出与 threads=8 sha256 相同 |
