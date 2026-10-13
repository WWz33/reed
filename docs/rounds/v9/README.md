# reed doc_v9：引擎终止语义复审 + 错误路径修复规范（2026-10-13）

> 本目录是 `docs/rounds/v8/` 之后的**第九轮规范**。它由两部分组成：
> 对当前实现（v8 修复全部落地后）**引擎错误路径与终止语义**的独立复审，
> 以及在复审基础上的**修复规范**。
> 本目录不修改任何代码，所有改动以规范形式给出，由实现者执行。
> 历史轮次 `docs/rounds/vN/` 只增不改；当前规范索引见 `docs/README.md`。

## 背景

v6_y5 落地了流水线骨架（票据窗口、许可 drop guard、接缝上限、中止原子位），
v8 统一了并行预算。**两轮都只覆盖了"正常运行"与"worker 侧失败"两条路径**：
`worker_panic_returns_every_permit`、`source_errors_propagate`、
`memory_budget_stops_the_run`、`external_abort_stops_the_run` 四个既有验收测试
全部经由 `FirstError::record` 触发中止，因此都能正常返回。

本轮以"逐条阻塞点反查唤醒来源"的方法复审，发现**第三条路径——reducer 侧
（调用方线程）自己出错——完全没有唤醒机制**：`deliver!` 宏的 5 处 `?` 直接把
`Err` 抛出 scope 闭包，既不 `first.record`，也不触发 abort，而 `ticket_tx`
要等 `run_parallel` 返回才析构、`run_parallel` 又卡在 `thread::scope` 的 join 上。
结果是**任何 reducer 侧错误（输出磁盘满、索引遇到未排序输入、接缝记录超 1 GiB、
畸形接缝行）都会让整条流水线永久挂死**，而不是返回错误。已用 `eu-stack`
抓到 5 线程全部停在 `futex_wait` 的现场（01-review.md §1.2）。

同一族还有三个静默问题：通道发送失败被当成成功（放大 D9-1 为硬死锁）、
文本输入末条未换行记录被三条路径分别处理成"丢弃 / 报错 / 收录"、
以及带过滤的 `stats` 在接缝上静默跳过畸形记录。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 引擎错误路径复审：阻塞点—唤醒来源台账、逐条对照 spec §7.5/§7.6/§7.7、缺陷清单 D9-1…D9-8（带 `file:line`、严重度、修法）、已降级/已消解的说法、复核门 |
| [02-architecture.md](02-architecture.md) | 修订规范：引擎终止语义（唯一错误出口 + 唤醒不变量）、通道发送结果语义、输入截断契约、seam 与 chunk 体的错误一致性、指标语义、spec 回改清单、验收纪律 |
| [06-roadmap.md](06-roadmap.md) | 修复路线 Q0–Q3：通用验收门、每步内容与门、顺序与依赖、本轮不做什么 |
| [STATUS.md](STATUS.md) | 实施状态（起始全部未修）+ 复审基线 |

## 复审方法【据】

1. **阻塞点—唤醒来源台账**：枚举 `run_parallel` 里全部 6 个阻塞点，逐个列出
   "谁能唤醒它"，凡唤醒来源依赖某线程正常返回的即为可疑（01-review.md §1.1）；
2. **失败注入实测**：在 `crates/reed-engine/tests/` 下临时构造失败 sink /
   失败 `run_seam` / panic sink / 未换行末行四类夹具，分别在
   "chunk 数 ≫ 票据窗口"与"chunk 数 < 票据窗口"两种规模下运行，
   用 `timeout` 判定是否挂死；
3. **现场取证**：对挂死进程用 `/proc/<pid>/task/*/wchan` + `eu-stack -p`
   抓全线程栈，逐帧确认每个线程停在哪个阻塞点；
4. **对照实验**：同一个失败 sink，只改变 chunk 数（3 vs 2788），
   证明挂死与否完全由"reader 是否已先到 EOF"决定；
5. **双 profile 交叉**：`cargo test`（debug）与 `cargo test --release`
   分别运行，暴露 `debug_assert` 类不变量对验收门不可见的问题。

复审基线：commit `5ef4454`（doc_v8 D8-1…D8-6 全部落地），工作区干净。

## 结论速览

| 维度 | 结论 |
|---|---|
| reducer 侧错误 | **挂死**（D9-1，P0）。实测 2788 chunk × `budget(4)`：sink 报错后进程 45 s / 90 s / 120 s 均不返回；`eu-stack` 显示 reader 停在 `select{signal, ticket_rx}`、2 个 worker 停在 `select{signal, chunk_rx}`、测试线程（= reducer）停在 `thread::scope` 的 join，5 线程全 `futex_wait` |
| 同族变体 | **同样挂死**（D9-1b/1c）。`run_seam` 返回 `Err` → 挂死；reducer 侧 panic → panic 消息已打印但进程永不返回（reducer 没有 `catch_unwind`，worker 有） |
| 为什么测试没抓到 | **夹具规模**（01-review.md §1.3）。对照实验：同一失败 sink，3 个 chunk（在 `budget(4)` 的 4 张票内）→ 971 µs 正常返回 `Err(Job("sink exploded"))`；2788 个 chunk → 永久挂死。既有 e2e 夹具全部落在票据窗口内 |
| 通道发送结果 | **失败当成功**（D9-2，P1）。`send(...) -> _ => true` 把 `SendError` 也算成已发送，worker 不 break 继续消费 chunk；这是把 D9-1 从"错误"放大成"硬死锁"的机制 |
| 输入截断 | **三条路径三种行为**（D9-3，P0）。实测 50 条记录去掉末行 `\n` → filter/stats 交付 **49** 条（顺序与并行路径一致，静默丢弃）；`reed index` **收录**该记录（`IndexJob` 解析 tail + `IndexSink::finish` 以 `eof_voff` 收尾）；BCF 路径**报错**（"BCF input ends inside a record"）。同一文件产生"索引里有、过滤输出里没有"的记录 |
| 接缝畸形记录 | **静默跳过**（D9-4，P1）。`FilteredJob::run_seam` 见 `LineRef::parse` 失败即 `Ok(default())`；同一条记录落在 chunk 体内则报 `"record line has fewer than 9 columns"`。是否致命取决于 chunk 对齐 |
| 指标语义 | **两处矛盾**（D9-5/6，P2）。`in_flight_peak_bytes` 文档写"压缩字节"，实测 = 压缩 × `MEM_FACTOR(8)`（实测比值恰为 8.000），且顺序路径恒为 0；`job_time` 文档写"所有 worker 的 run_chunk/run_seam"，实测同输入同 273 ms 墙钟下顺序路径 `job_time=273.12 ms`、并行路径 `job_time=49.2 µs`（reducer 上的 seam 时间漏计约 5600 倍） |
| spec 漂移 | **§7.7 是过时的一方**（D9-7，P2）。`docs/spec/architecture.md:841` 写"超出时报 `MemBudgetExceeded`"，但 `docs/rounds/v2/03-architecture.md` §4.2 早已改为背压（"超出预算时 reader 阻塞，而不是报错"），实现遵循的是 v2 03 §4.2。**应回改 spec，不改代码** |
| 验收纪律 | **release-only 门有盲区**（D9-8，P1）。`stitcher_reassembles_and_caps_fragments` 自 Y5（`111dbde`）加入起从未在 debug 下绿过，但 `cargo test --release` 下 `debug_assert` 被编译掉、测试通过。v8/STATUS.md 将其登记为"预存在、与本轮无关"——本轮更正：它正是 D9-3 缺失的 finish 语义，测试最后一行 `push_head(None)` 就是在够一个不存在的"流结束"API |

## 决策状态

| 决策 | 状态 |
|---|---|
| 引擎只有一个错误出口：任何线程的错误都必须经 `FirstError::record` 触发 abort | 已定（spec §7.5 的原有语义；本轮首次给出 reducer 侧的落地形态，02 §2） |
| 文本输入以未换行的残行结束时，行为与 BCF 路径一致（报错） | **【决】需拍板**：报错（推荐，与 BCF 一致）／警告后保留该记录／维持现状静默丢弃。见 02 §4.3 |
| `in_flight_peak_bytes` 的语义取"估算足迹"（改文档）还是"压缩字节"（改代码） | **【决】需拍板**：推荐改文档 + 顺序路径填真实值（02 §6.1） |
| spec §7.7 的 `MemBudgetExceeded` 表述回改为背压 | 已定（实现遵循 v2 03 §4.2，spec 是过时的一方；02 §7） |
| 验收门必须包含一次 debug 运行 | 已定（对齐 v2 06 §4 纪律 1"测试红着不提交"；02 §8） |

## 与前轮的关系

- **不重开** v6 的 Y 系列、v7 的 Z 系列、v8 的 P 系列结构项；
- v8 判定"流水线本体合规"**继续有效**——本轮发现的是**错误路径**缺口，
  不是线程预算缺口；v8/01 §1 的线程创建点台账无需修订；
- v6_y5 §4 已登记的 `FrameChain::wait_inner` 10 ms 轮询偏差**不重复计数**
  （见 01-review.md §4）；
- v8/STATUS.md:22-23 关于红测试"与本轮无关"的判断由本轮**更正**（D9-8）。

## 标记约定

【规】必须；【荐】推荐；【决】需要拍板；【据】有实测或源码证据；**P0/P1/P2** 为修复优先级。

严重度定义（本轮口径）：
**P0** = 挂死、静默丢数据或产生自相矛盾的输出；
**P1** = 可观测的错误行为、或让 P0 逃过验收的流程缺口；
**P2** = 文档/指标语义矛盾（不影响输出正确性）。
