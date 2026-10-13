# 01 引擎错误路径与终止语义复审

> 复审对象：commit `5ef4454`（doc_v8 D8-1…D8-6 全部落地后的当前实现），工作区干净。
> 方法：阻塞点—唤醒来源台账 + 失败注入实测 + `eu-stack` 现场取证 + 双 profile 交叉。
> 基准规范：`docs/spec/architecture.md` §7.5（中止与错误）、§7.6（确定性契约）、
> §7.7（内存预算）、§7.10（阶段指标）；`docs/rounds/v2/03-architecture.md` §4（引擎）。
> 本轮不修改任何代码；所有实测夹具是可复现的临时测试，修法见 02/06。

---

## 0. 总评

v6_y5 与 v8 把**正常路径**与**worker 侧失败路径**做得很扎实：票据窗口、许可
drop guard（`PermitGuard` 的析构是唯一释放路径）、接缝 1 GiB 上限、
`FirstError::record → trigger` 的即时唤醒、rayon 三阶段接 `ThreadBudget`，
这些都有对应验收测试且实测通过。

问题集中在**一条从未被覆盖的失败路径**：

1. **reducer 侧（调用方线程）自己出错时没有任何唤醒机制**（D9-1，P0）。
   `deliver!` 宏的 5 处 `?` 把 `Err` 直接抛出 `thread::scope` 闭包，
   既不 `first.record`，也不 drop 任何 sender；`pipeline.rs:521-524` 的注释
   "Errors return immediately — the dropped senders/receivers unblock every
   other stage" **与代码不符**。实测挂死，`eu-stack` 现场见 §1.2。
2. **通道发送失败被当成成功**（D9-2，P1）。worker 与 reader 的
   `send(...) -> _ => true` 把 `SendError` 也算已发送，这是把 D9-1
   从"一个错误"放大成"硬死锁"的机制。
3. **输入截断的契约在三条路径上互不相同**（D9-3，P0）。文本 filter/stats
   静默丢末条记录、`reed index` 收录它、BCF 路径报错。
4. **seam 与 chunk 体的错误判定不一致**（D9-4，P1）。同一条畸形记录，
   落在体内是致命错误，跨接缝则被静默跳过。
5. **指标语义与文档/spec 矛盾**（D9-5/6/7，P2）。
6. **验收门只跑 release，`debug_assert` 类不变量失效无人值守**（D9-8，P1）。
   一个自 Y5 起就红的测试被 v8 登记为"预存在、与本轮无关"，实际是 D9-3 的症状。

D9-1 的现实触发条件**全部是运维级常见事件**，不需要构造畸形输入：

| 触发源 | 位置 | 场景 |
|---|---|---|
| `IndexSink::push` 未排序 | `reed/src/index.rs:402-405` | 对未按位置排序的大 VCF 建索引 → `"positions are not ascending"` |
| `OutputSink` 写失败 | filter/stats 的 sink `accept` | 输出盘 ENOSPC、管道 EPIPE（`reed filter … \| head`） |
| 接缝记录超上限 | `pipeline.rs:139-142` | 无换行的病态输入累积到 1 GiB |
| `run_seam` 报错 | `pipeline.rs:535` | 畸形行恰好跨 chunk 接缝（含 D9-4 修好之后的正常报错路径） |
| reducer panic | `pipeline.rs:531-548` | sink/job 的任何 `expect`/索引越界；reducer 无 `catch_unwind` |

也就是说：**"100 GB 文件 + 输出盘写满"今天是永久挂死，而不是报错退出。**

---

## 1. 实测结果【据】

以下全部在 commit `5ef4454`、debug profile、`--test-threads=1` 下测得。
夹具：`fixture_text(20_000)` 条 VCF 记录，`make_chunks(&text, 300)` 切成
**2788 个 chunk**；`budget(4)` → `workers()=2`、`tickets()=4`。

### 1.1 阻塞点—唤醒来源台账

`run_parallel` 共 6 个阻塞点。"唤醒来源"列出**所有**能让它返回的事件：

| # | 线程 | 位置 | 阻塞点 | 唤醒来源 | 判定 |
|---|---|---|---|---|---|
| B1 | reader | `pipeline.rs:363-366` | `select!{ recv(signal), recv(ticket_rx) }` | ① abort 触发 ② 取到票 ③ `ticket_tx` 全部析构 | ⚠ ③ 依赖 `run_parallel` 返回，而它卡在 join |
| B2 | reader | `:398` | `permits.acquire(footprint, signal)` | ① abort ② 许可释放 | ✅ condvar，abort 已注册（`abort.rs:31`、`budget.rs:169-187`） |
| B3 | reader | `:407-410` | `select!{ recv(signal), send(chunk_tx) }` | ① abort ② worker 取走 chunk ③ `chunk_rx` 全部析构 | ⚠ 同 B1；且发送失败被当成成功（D9-2） |
| B4 | worker | `:440-446` | `select!{ recv(signal), recv(chunk_rx) }` | ① abort ② reader 送来 chunk ③ `chunk_tx` 全部析构 | ⚠ ③ 依赖 reader 退出 |
| B5 | worker | `:498-507` | `select!{ recv(signal), send(result_tx) }` | ① abort ② reducer 取走结果 ③ `result_rx` 析构 | ❌ reducer 出错时 ②③ 都不发生（D9-1/D9-2） |
| B6 | reducer | `:565-580` | `select!{ recv(signal), recv(result_rx) }` | ① abort ② worker 送来结果 ③ `result_tx` 全部析构 | ✅ ③ 成立（`:515` drop 了本体，只剩 worker clone） |

**关键观察**：B1/B3/B4/B5 的"安全"唤醒来源都依赖别的线程正常退出，
唯一的**无条件**唤醒来源是 abort（①）。而 abort 只由
`FirstError::record`（`abort.rs:160-166`）或外部 `handle.trigger()` 触发。

**reducer 的错误路径两者都不做**，于是 ① 不成立；
`ticket_tx` 是 `run_parallel` 的局部变量（`:321`），只在 `:544` 被闭包借用发送，
闭包 `return Err` 后借用结束但**本体要等 `run_parallel` 返回才析构**，
而 `run_parallel` 正卡在 `thread::scope` 的 join（`:343`/`:603`）——③ 也不成立。
循环等待成立。

### 1.2 挂死现场（`eu-stack`，2788 chunk × `budget(4)`）

失败 sink 在第 5 个 chunk 返回 `Err`：

```
running 1 test
test f1_reducer_sink_error_returns_instead_of_hanging ... F1: n_chunks = 2788
（45 s / 90 s / 120 s 超时均无返回）
```

`/proc/<pid>/task/*` 共 5 个线程，全部 `State: S`、`wchan: futex_wait_queue`、
`syscall_nr: 202`（futex）。逐线程栈（符号已还原）：

| TID | 角色 | 栈顶关键帧 |
|---|---|---|
| 233520 | test harness | `std::sync::mpsc … CompletedTest::recv` → `Thread::park`（等测试完成） |
| 233521 | **测试线程 = reducer** | `pipeline::run_parallel` 的 `std::thread::scope` → `Thread::park`（**join 永不完成**） |
| 233522 | reader | `crossbeam_channel::select::run_select` → `park`（B1） |
| 233523 | worker | `crossbeam_channel::select::run_select` → `park`（B4/B5） |
| 233524 | worker | `crossbeam_channel::select::run_select` → `park`（B4/B5） |

线程数与 `budget(4)` 的布局（1 reader + 1 reducer + 2 worker + harness）完全吻合。
注意 reducer 栈顶是 **scope 的 join**，不是 `deliver!`——它已经"返回"了错误，
只是这个返回永远走不出 scope。

### 1.3 对照实验：为什么既有测试抓不到

同一个失败 sink，只改 chunk 数：

| 夹具 | chunk 数 | 票据窗口（`budget(4)`） | 结果 |
|---|---|---|---|
| `fixture_text(2_000)` / `make_chunks(_, 40_000)` | **3** | 4 | ✅ **971 µs 正常返回** `Err(Job("sink exploded"))` |
| `fixture_text(20_000)` / `make_chunks(_, 300)` | **2788** | 4 | ❌ **永久挂死** |

差别只在：3 个 chunk 时 reader 在 reducer 出错前就到达 EOF 并 break（`:380`），
`chunk_tx` 析构 → worker 的 B4 走 ③ 退出 → `result_tx` 析构 → join 完成。
2788 个 chunk 时 reader 早在第 4 张票上就阻塞于 B1，永远到不了 EOF。

**既有 e2e / 引擎验收测试的夹具全部落在票据窗口内**，所以这条路径从未被执行到。
`reed-engine/tests/pipeline.rs` 的 4 个错误测试
（`worker_panic_becomes_worker_panicked`、`source_errors_propagate`、
`memory_budget_stops_the_run`、`external_abort_stops_the_run`）
全部经由 `FirstError::record` 或外部 `trigger`，即台账里的 ①，因此都能返回。

### 1.4 同族变体实测

| 变体 | 注入点 | 结果 |
|---|---|---|
| D9-1b | `run_seam` 返回 `Err`（`pipeline.rs:535`） | ❌ 挂死（`timeout 40` → exit 124，测试名已打印、无返回） |
| D9-1c | reducer 侧 sink `panic!` | ❌ 挂死：panic 消息 `"reducer sink panicked"` **已打印**，但进程永不返回。reducer 循环（`:549-591`）没有 `catch_unwind`，worker 有（`:450`）——两侧不对称 |

### 1.5 输入截断实测（D9-3）

`fixture_text(50)` 后 `text.pop()` 去掉末行 `\n`，`make_chunks(_, 200)` → 10 chunk：

```
F3: threads=1 delivered 49 of 50 records (n_chunks=10)
assertion failed: left=49 right=50
```

顺序（threads=1）与并行（threads=4）路径**都丢**最后一条，且两者互相一致——
所以 `sequential_and_parallel_outputs_agree`（S4 门）与
`output_is_invariant_across_thread_counts`（确定性门）**都抓不到它**：
它们比较的是两条路径之间，而不是与输入真值之间。

三条路径的行为对照（源码核实）：

| 路径 | 末条未换行记录 | 依据 |
|---|---|---|
| 文本 filter / stats | **静默丢弃** | `NewlineStitcher` 无 finish/flush（`pipeline.rs:132-182`）；并行 drain 结束 `:592-600`、顺序结束 `:705-710` 都只调 `sink.finish()`；`owned_body_range` 把 body 止于最后一个 `\n`（`chunk/mod.rs:246-251`） |
| `reed index`（文本 VCF） | **收录并建索引** | `IndexJob::run_chunk` 解析 `view.tail`（`index.rs:188-194`）；`IndexSink::accept` push 它（`:439-441`）；`finish` 用 `eof_voff` 收尾（`:446-449`） |
| BCF（filter/stats/index） | **报错** | `BcfRunner::finish` 的 `expect_complete`（`bcf_input.rs:119-133`）、`index.rs:546-556`：`"BCF input ends inside a record (truncated file)"` |

后果：同一个文件，`reed index` 产出的索引**指向一条 `reed filter` 输出里不存在的记录**。

### 1.6 指标实测（D9-5 / D9-6）

**(a) `in_flight_peak_bytes` 不是压缩字节**

单 chunk 夹具，`max_compressed = 149`：

```
F5: peak=1192  max_compressed=149  peak/max=8.000  MEM_FACTOR=8
```

比值恰为 8.000 = `MEM_FACTOR`（`budget.rs:118`）。
`metrics.rs:13` 的文档写的是 "Peak in-flight compressed bytes"。
另外顺序路径把它**硬编码为 0**（`pipeline.rs:715`），
所以 `--threads 1` 下这个指标永远没有信息量。

**(b) reducer 上的 seam 时间在并行路径完全漏计**

`run_seam` 内 `sleep(2ms)`，134 个 seam-rich chunk，两条路径同输入：

```
threads=1  wall=273.33 ms  job_time=273.12 ms      reduce_time=7.5 µs
threads=4  wall=273.83 ms  job_time=49.2 µs        reduce_time=3.4 µs
```

墙钟相同（273 ms，全是 seam 的 sleep），但并行路径的 `job_time` 只有
49.2 µs——**约 5600 倍的 seam 工作量在指标里不存在**，
且没有落到任何其它计数器（`reduce_time` 3.4 µs 也不含它）。
原因：`run_seam` 在 reducer 上执行（`pipeline.rs:535`）而 `job_ns` 只在
worker 里累加（`:480`）；seam 的 `sink.accept`（`:536`）也在 `reduce_ns`
的计时区间（`:538-540`）之外。
顺序路径反而把 `run_seam` 计进了 `job_time`（`:678-680`）——
所以 `metrics.rs:34-36` 的文档对顺序路径成立、对并行路径不成立，
两条路径的同一字段语义不同。

---

## 2. 与规范逐条对照

### 2.1 spec §7.5：「所有阻塞点都用 `select!` 同时等待数据和 abort，中止即时生效」

> 原文（`architecture.md:820-823`）："`AbortSignal`：一个只用来关闭的 crossbeam 通道。
> 所有阻塞点都用 `select!` 同时等待数据和 abort，中止即时生效，不再用 20 ms 轮询。
> `FirstError`：首个错误胜出…"

**判定：对 reducer 的错误路径未实现。** §7.5 的机制本身是对的（§1.1 台账里
① 一栏全部成立），但它的前提是**错误必须先进 `FirstError`**。
`deliver!` 的 5 处 `?` 绕过了 `FirstError`，于是 abort 从未触发，
6 个阻塞点里有 4 个（B1/B3/B4/B5）失去了唯一的无条件唤醒来源。
这不是"轮询 vs 即时"的问题，而是"根本没有信号"。

`pipeline.rs:521-524` 的注释是这条规范在代码里的**错误落地表述**：

```
// Reducer (this thread): deliver in chunk order, chain-stitching
// seam records between chunks (§6.2 step 4, §7.6 rule 3). Errors
// return immediately — the dropped senders/receivers unblock every
// other stage, and `run` reports the recorded cause.
```

三处与代码不符：① 出错路径没有 drop 任何 sender；② "the recorded cause"
不存在（从未 record）；③ 因此 `run` 不报告任何东西——它不返回。

### 2.2 spec §7.6 规则 5：「同一版本、同一输入、同一 Spec ⇒ 输出对线程数逐字节一致」

**判定：对 D9-1 成立但无意义，对 D9-3 是盲区。**
D9-3 在顺序与并行路径上**同样丢**末条记录（§1.5），所以
`sequential_and_parallel_outputs_agree` 与
`output_is_invariant_across_thread_counts` 两个门都通过。
§7.6 只约束"路径之间一致"，**不约束"与输入真值一致"**——
本轮要求补一条真值门（02 §4.4、06 §1 门 4）。

### 2.3 spec §7.7：「超出时报 `EngineError::MemBudgetExceeded`」

> 原文（`architecture.md:838-841`）："引擎负责运行时监控（在途字节计数），
> 超出时报 `EngineError::MemBudgetExceeded`，不会默默 OOM。"

**判定：spec 是过时的一方，实现正确。**
`docs/rounds/v2/03-architecture.md` §4.2 早已把它改成背压：

> "令牌改为按字节计…超出预算时 reader **阻塞，而不是报错**。
> `MemBudgetExceeded` 只保留给'单个 chunk 本身就超过预算'这种无法继续的情况。"

实现与此完全一致（`budget.rs:107-109`、`:145`、`pipeline.rs:390-396`），
且 v2 06 R3 的门"内存上界测试（在途字节 ≤ 预算）"由
`memory_budget_bounds_in_flight_bytes` 落地并通过。
**应回改 `architecture.md:841`**，不改代码（02 §7）。

同一节的 §7.4（`:811-818`）描述的票据背压与实现一致，无需改动。

### 2.4 spec §7.10：「每个阶段记录忙碌时间、空闲时间、字节数、记录数」

**判定：部分失效。** §7.10 的承诺是"每个阶段"可观测，但 §1.6(b) 显示
**reducer 阶段的 job 时间（seam）完全没有计数器**：
它既不在 `job_time`（只在 worker 累加），也不在 `reduce_time`（计时区间之外）。
在 seam 密集的输入上（小 chunk、宽记录），reducer 的真实负载对
`--profile run.json` 不可见，"reducer-bound"因此不可判定——
这正是 doc_v4 03 A4-5 引入 reader/reducer 计时想消除的那类推断。

### 2.5 v6_y5 §4：`FrameChain::wait_inner` 的 10 ms 轮询

**判定：已登记的已知偏差，本轮不重复计数。**
`docs/spec/v6_y5/README.md:38-41`、`:167-168` 已明确登记理由
（要同时等 abort/poison/finish 三个来源）与中止延迟上界（≤10 ms）。
见 §4。

---

## 3. 缺陷清单

严重度口径见 README「标记约定」。

| # | 严重度 | 位置 | 缺陷 | 触发条件 / 后果 | 修法 |
|---|---|---|---|---|---|
| D9-1 | **P0** | `reed-engine/src/pipeline.rs:531-548`（`deliver!` 的 `?` 在 `:534`/`:535`/`:536`/`:541`/`:545`）、热循环 `:549-591`、注释 `:521-524` | reducer 侧错误直接 `return Err` 出 scope 闭包，不经 `FirstError::record`、不触发 abort、不 drop 任何 sender；`ticket_tx`（`:321`）要等 `run_parallel` 返回才析构，而它卡在 `thread::scope` join | 输出盘 ENOSPC / 管道 EPIPE / 索引遇未排序输入（`index.rs:402-405`）/ 接缝记录超 1 GiB（`:139-142`）/ `run_seam` 报错 → **整条流水线永久挂死**。实测 2788 chunk × `budget(4)` 120 s 不返回；`eu-stack` 5 线程全 `futex_wait`（§1.2） | reducer 的错误出口统一为一个 helper：`first.record(err)` → `drop(result_rx)` → `return Err(Aborted)`；`run` 已优先报告记录的因（`:231-233`）。同时把 `:521-524` 的注释改成实际语义。详见 02 §2 |
| D9-1b | **P0** | 同上，`:535` | `run_seam` 的错误走同一条无唤醒路径 | 畸形行跨接缝、或 D9-4 修好后正常的报错路径 → 挂死（实测 exit 124） | 同 D9-1（同一个出口即覆盖） |
| D9-1c | **P0** | `pipeline.rs:549-591`（reducer 循环无 `catch_unwind`；worker 侧 `:450` 有） | reducer panic 两侧不对称：worker panic 被转成 `WorkerPanicked` 并 record，reducer panic 直接 unwind 出闭包 | sink/job 的任何 `expect`、索引越界 → panic 消息打印后**进程永不返回**（实测）。比 D9-1 更糟：用户看到了 panic 却等不到退出码 | reducer 循环包 `catch_unwind(AssertUnwindSafe(...))`，unwind 时 `first.record(ReducerPanicked{…})` 后走 D9-1 的出口。错误变体见 02 §2.4 |
| D9-2 | **P1** | `pipeline.rs:498-510`（worker→result）、`:407-416`（reader→chunk） | `send(tx, …) -> _ => true` 把 `Result` 绑到 `_`，`SendError` 也算 `sent = true` | 对端已 drop 时：worker 不 break、继续取 chunk 做无用功；reader 把"chunk 没进管道"当成成功。在 abort 路径上被 reader 停摆掩盖，但它是把 D9-1 变成**硬死锁**（worker 永不自行退出）的机制 | 两处都改 `-> res => res.is_ok()`；`if !sent { break }` 的注释补上"对端断开也是 break 条件"。详见 02 §3 |
| D9-3 | **P0** | `NewlineStitcher` 无 finish/flush（`pipeline.rs:132-182`）；并行 drain 结束 `:592-600`；顺序结束 `:705-710` | 流结束时残留的 fragment 被静默丢弃 | 末行无 `\n` 的 VCF：**实测 50 条只交付 49 条**（顺序与并行一致，所以 §7.6 的门全绿）。同时 `reed index` **收录**该记录（`index.rs:188-194` + `:446-449`），BCF 路径**报错**（`bcf_input.rs:119-133`、`index.rs:546-556`）→ 索引指向过滤输出里不存在的记录；文本截断被隐藏 | 加 `NewlineStitcher::finish(self) -> Option<Vec<u8>>`；两条路径在 drain 后、`sink.finish()` 前检查，残留即 `EngineError::job("input ends inside a record")`，与 BCF 一致。**行为变更，需拍板**（02 §4.3） |
| D9-4 | **P1** | `reed/src/filtered.rs:327-329` | `run_seam` 见 `LineRef::parse` 失败即 `return Ok(J::Partial::default())` | 带 site/threshold 过滤的 `reed stats` 遇到 <9 列的行：跨接缝则静默跳过，落在 chunk 体内则报 `"record line has fewer than 9 columns"`（`record_pipeline.rs:165-166`、`:254-255`）。**是否致命取决于 chunk 对齐**，即取决于线程数与文件大小——违反确定性契约的意图 | 删掉这个提前返回，让 `filter_line`/`scan_line` 走 `pipeline.process` 报同一个错。详见 02 §5 |
| D9-5 | **P2** | `metrics.rs:13-14`；`pipeline.rs:389`、`:715`；`budget.rs:118`、`:137-139` | `in_flight_peak_bytes` 文档写"压缩字节"，实际存估算足迹 = 压缩 × 8（**实测比值 8.000**）；顺序路径恒为 0 | `--profile run.json` 与 §17.3 基准日志里这个字段被误读为压缩在途量，差 8 倍；`--threads 1` 下恒为 0，无信息量 | 改文档为"估算足迹（压缩 × `MEM_FACTOR`）"，字段更名 `in_flight_peak_footprint_bytes`；顺序路径填真实的单 chunk 足迹峰值。**语义取舍需拍板**（02 §6.1） |
| D9-6 | **P2** | `metrics.rs:22-23`、`:34-36`；`pipeline.rs:535-536`、`:538-540`、`:678-680` | `job_time` 文档写"所有 worker 的 run_chunk/run_seam"，但并行路径 `run_seam` 在 reducer 上执行且**不计时**；seam 的 `sink.accept`（`:536`）也在 `reduce_ns` 区间（`:538-540`）之外 | 实测同输入同 273 ms 墙钟：顺序 `job_time=273.12 ms`、并行 `job_time=49.2 µs`（漏计约 5600 倍）。同一字段在两条路径语义不同；seam 密集输入上"reducer-bound"不可判定，违背 §7.10 | 并行路径给 `run_seam` 与 seam 的 `accept` 各加计时：seam 的 job 时间并入 `job_time`（或新增 `seam_time`），seam 的 accept 并入 `reduce_time`；`metrics.rs` 文档写明"跨两条路径同语义"。详见 02 §6.2 |
| D9-7 | **P2** | `docs/spec/architecture.md:841`（§7.7） | spec 写"超出时报 `MemBudgetExceeded`"，实现是背压阻塞、只有单个超大 chunk 才失败 | spec 与实现矛盾；实现遵循的是 `v2/03-architecture.md` §4.2 的后续修订，spec 未回改 → 后来者按 spec 改代码会**破坏**已落地的背压与 v2 06 R3 的门 | **回改 spec**（不改代码）：§7.7 改为"背压为主，`MemBudgetExceeded` 只用于单个 chunk 超过全部预算"，并注明依据 v2 03 §4.2。详见 02 §7 |
| D9-8 | **P1** | `pipeline.rs:765-777`（测试）、`:157-160`（`debug_assert`）；`v8/STATUS.md:22-23`；`v2/06-roadmap.md:10` 与 `v8/06-roadmap.md:11` 的门 | `stitcher_reassembles_and_caps_fragments` 最后一行在 fragment 非空时调 `push_head(None)`，撞 `debug_assert`。该测试 Y5（`111dbde`）加入、`debug_assert` M3（`04c1f1a`）就有 → **从未在 debug 下绿过**；但验收门写的是 `cargo test --release`，release 下 `debug_assert` 被编译掉、测试通过 | ① 一个红测试在仓库里存活了两轮，违反 v2 06 §4 纪律 1"测试红着不提交"；② **所有 `debug_assert` 类不变量对验收门不可见**；③ v8 把它登记为"预存在、与本轮无关"是误判——它正是 D9-3 的症状：测试最后一行就是在够一个不存在的"流结束"API | ① 按 D9-3 加 `finish()`，把测试最后一行改成 `stitcher.finish()` 的正确断言（02 §4.2）；② 验收门增加一次 debug 运行（02 §8）；③ 更正 v8/STATUS.md 的归因（本目录已登记） |

### 3.1 D9-1 的循环等待（形式化）

```
reducer 出错 → return Err（未 record、未 trigger、未 drop sender）
   ↓
abort 从未触发            ⇒ B1/B3/B4/B5 失去唯一无条件唤醒源
ticket_tx 未析构          ⇒ B1 的 ③ 不成立（reader 阻塞在票上）
reader 阻塞               ⇒ chunk_tx 未析构 ⇒ B4 的 ③ 不成立
worker 阻塞在 B4/B5       ⇒ result_tx 未析构（但 reducer 已不再 recv）
   ↓
thread::scope join 等待 reader + workers ⇒ 永不完成
   ↓
run_parallel 永不返回 ⇒ ticket_tx 永不析构  ← 回到起点
```

打破环的最小改动只有一个：**让 reducer 的错误进入 `FirstError`**（触发 abort，
即台账里的 ①）。`drop(result_rx)` 是加速项（让 B5 的 ③ 也成立），
不是必要条件——但它让 D9-2 修好后 worker 能立即 break，避免多做无用功。

### 3.2 修法在代码里已有先例【据】

正常 drain 路径（`pipeline.rs:592-594`）的注释与代码是：

```rust
// Draining finished; dropping the receiver unblocks any worker
// still sending, so the scope join below cannot deadlock.
drop(result_rx);
```

即**作者完全理解"断开 `result_rx` 才能避免 scope join 死锁"这个机制**，
并在成功路径与三个 abort 出口（`:551`、`:561`、`:587`）上正确应用了它。
D9-1 的成因不是概念缺失，而是**错误路径漏掉了同一个动作**——
`deliver!` 的 5 处 `?` 是 5 个绕过这层处理的旁路出口。

这条先例同时决定了修法的形状：不需要新机制，只要把已有的
"`record` + `drop(result_rx)` + `return Err(Aborted)`"三件套收成一个出口，
让 5 处 `?` 与 `sink.finish()?` 都走它（02 §2.2）。
`:551`/`:561`/`:587` 三处也应一并收进同一个 helper——它们本身是对的
（`first.failed()` 为真意味着 abort 已触发），只是重复了三遍。

---

## 4. 已降级 / 已消解的说法（防止再被误引）

| 说法 | 出处 | 本轮核实结论 |
|---|---|---|
| "`FrameChain::wait_inner` 的 10 ms 轮询违反 §7.5" | `architecture.md:820-823` | **已登记的已知偏差，不是新缺陷**。`v6_y5/README.md:38-41` 与 `:167-168` 已写明理由（同时等 abort/poison/finish 三源）与上界（中止延迟 ≤10 ms，每 chunk 一次而非每记录一次）。本轮不重复计数，也不要求本轮修 |
| "`stitcher_reassembles_and_caps_fragments` 失败是预存在的、与本轮无关" | `v8/STATUS.md:22-23` | **归因错误**。它确实是预存在的（干净 HEAD `e50f6fd` 与本轮 `5ef4454` 同样失败），但**不是无关**：它是 D9-3（`NewlineStitcher` 缺 finish 语义）的直接症状，且暴露了 D9-8（release-only 门的盲区）。v8 的判断在"不阻塞 v8 交付"上是对的，在"与引擎缺陷无关"上是错的 |
| "reducer 出错会立即返回，dropped senders 解锁其它阶段" | `pipeline.rs:521-524` 注释 | **与代码不符**（§2.1）。出错路径没有 drop 任何 sender，`run` 也不返回。该注释必须随 D9-1 一起改 |
| "spec §7.7 的 `MemBudgetExceeded` 是实现的既定行为" | `architecture.md:841` | **spec 是过时的一方**。`v2/03-architecture.md` §4.2 已改为背压，实现与 v2 一致并有 `memory_budget_bounds_in_flight_bytes` 门。修 spec，不修代码 |
| "`debug_assert` 保护了 stitcher 的首块不变量" | `pipeline.rs:157-160` | **生产不可达，保护是虚的**。`push_head(None)` 只在 chunk 0 发生：`split_text` 的 `head: None` 分支要求 `chunk.id == 0` 或 `first_in_block_offset > 0`（`chunk/mod.rs:213-243`），而两个源都只给 chunk 0 设非零偏移（`region.rs:146` 的 `mem::take`、`bgzf.rs:461-465` 的 `next_chunk_id == 1` 即"刚发出 id 0"）。所以它只会被测试触发——正是 D9-8 |
| "确定性门（§7.6）覆盖了输出正确性" | `architecture.md:830-836` | **只覆盖路径间一致，不覆盖真值**。D9-3 在两条路径上同样丢记录，两个确定性门全绿（§1.5、§2.2）。需补真值门 |

---

## 5. 复核门（每条缺陷的验证方式）

测试名即 06-roadmap.md 里要求先红后绿的用例名。所有挂死类测试**必须**带
`timeout` 语义（用 `Duration` 上界断言，不用 `#[should_panic]`），
否则红的时候会把 CI 挂住——这本身是 D9-1 的教训。

| # | 门 |
|---|---|
| D9-1 | `reducer_sink_error_returns_instead_of_hanging`：2788 chunk × `budget(4)`，sink 在第 5 个 chunk 报错 → 断言 **5 s 内返回**且 `Err` 携带 sink 的原始错误（不是 `Aborted`）。修复前先红（挂死）。同规模再跑 `run_seam` 报错版 `reducer_seam_error_returns_instead_of_hanging`（D9-1b） |
| D9-1c | `reducer_panic_returns_instead_of_hanging`：sink 在第 5 个 chunk `panic!` → 断言 5 s 内返回，错误为 `ReducerPanicked`（或复用的 panic 变体），且 `abort.signal.is_aborted()` |
| D9-1（规模不变量） | 上述三个测试各跑两个规模：chunk 数 ∈ {3（窗口内）, 2788（窗口外）}，**两者都必须在 5 s 内返回同一个错误类型**。这是 §1.3 对照实验的固化，防止将来只在窗口内测试 |
| D9-2 | `worker_send_on_disconnected_result_channel_breaks`：reducer 提前 drop `result_rx` 后，worker 必须在有限步内退出（断言 worker 处理的 chunk 数有上界，而不是"跑完全部 2788 个"）。reader 侧同型：`reader_send_on_disconnected_chunk_channel_breaks` |
| D9-3 | `unterminated_final_record_is_rejected`：50 条记录去掉末行 `\n` → 顺序与并行都返回 `Err`，消息含 "ends inside a record"。`terminated_final_record_is_accepted`：带 `\n` 的正常输入不受影响（仍交付 50 条）。`index_and_filter_agree_on_truncated_input`：同一截断文件跑 `reed index` 与 `reed filter`，两者**都失败**（当前是一个成功一个静默丢） |
| D9-3（真值门，新增） | `delivered_records_equal_input_records`：对以 `\n` 结束的输入，`sink_report` 必须等于用串行扫描数出的记录数——把 §2.2 的盲区补上（现有门只比较路径之间） |
| D9-4 | `malformed_seam_record_is_fatal`：构造一条 <9 列的记录**强制跨接缝**（小 chunk），断言带过滤的 stats 报错且消息与体内路径逐字相同（`"record line has fewer than 9 columns"`）。`malformed_record_fatality_is_chunk_alignment_independent`：同一条行放在体内 vs 跨接缝，两种情况的错误类型与消息必须相同 |
| D9-5 | `in_flight_peak_matches_its_documented_semantics`：断言 peak 等于该 run 中"单 chunk 估算足迹"的最大在途和，并与 rustdoc 声明的口径一致（改文档后按足迹口径断言）。顺序路径断言 peak **非 0**（当前恒 0） |
| D9-6 | `seam_time_is_counted_on_both_paths`：`run_seam` 内 sleep 2 ms × N 个 seam，断言并行路径 `job_time ≥ N × 2 ms × 0.8`（当前 49 µs），且顺序/并行两条路径的 `job_time` 落在同一量级（当前差 5600 倍）。`reduce_time_covers_seam_accept`：seam 的 accept 时间计入 `reduce_time` |
| D9-7 | `grep -n "MemBudgetExceeded" docs/spec/architecture.md` 的上下文出现"背压"/"单个 chunk"；`ci/check-doc-drift.sh` 加一条守卫：spec §7.7 不得再声称"超出即报错" |
| D9-8 | `cargo test -p reed-engine --lib stitcher`（**debug**）全绿；`cargo test --release -p reed-engine --lib stitcher` 全绿；06 §1 通用门新增"至少一次 debug 运行"，`ci/` 脚本体现；`v8/STATUS.md:22-23` 的归因由本目录更正（v8 只增不改，故在 v9 登记） |
| 全局 | `cargo test --workspace`（debug）与 `cargo test --release --workspace` 均全绿；threads ∈ {1,2,3,7,56} 输出 sha256 相同（filter 与 stats）；`bash ci/check-doc-drift.sh` 全绿；`bench/suites/e2e/results.tsv` 追加一行证明性能未回退 |
