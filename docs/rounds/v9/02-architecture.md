# 02 修订规范：引擎终止语义、截断契约与指标语义

> 目标：**任何错误都必须让 `run()` 返回**；**同一份输入在 filter / stats / index
> 三条路径上必须有同一个"是否完整"的判定**；**指标字段名与其实际口径一致**。
> 本轮是缺陷修复轮，不重开 v6/v7/v8 的结构项（Y/Z/P 系列），
> 不改票据窗口、许可 drop guard、`FrameChain`、rayon 预算接线——
> 这些经 v8 复审判定合规，本轮实测也未发现问题。
>
> 标记：【规】必须；【荐】推荐；【决】需拍板；【据】有实测或源码证据。

---

## 1. 现状结构问题

| # | 问题 | 证据 |
|---|---|---|
| S1 | **错误出口不唯一**。引擎有四种失败形态，只有两种进了 `FirstError`：worker 的 job 错误与 panic（`pipeline.rs:494-497`）、reader 的传输错误与预算超限（`:381-384`、`:391-395`）会 record；reducer 的 job 错误、sink 错误、seam 超限、panic（`:534-545`、`:549-591`）**全部绕过** | 01-review §1.1 台账、§3.1 循环等待 |
| S2 | **panic 处理两侧不对称**。worker 有 `catch_unwind` 并转成 `WorkerPanicked`（`:450`、`:484-492`），reducer 没有 | 01-review §1.4 D9-1c 实测 |
| S3 | **通道发送的 `Result` 被丢弃**。两处 `send(...) -> _ => true` 把断开当成功（`:409`、`:506`） | D9-2 |
| S4 | **`NewlineStitcher` 没有生命周期终点**。只有 `push_head`/`push_tail`，没有 finish/flush，流结束时残留 fragment 无处可去，只能被 drop（`:132-182`、`:592-600`、`:705-710`） | D9-3；也是 D9-8 红测试的根因 |
| S5 | **"截断"的判定分散在三处且互不相同**。引擎（文本）不判定、BCF runner 判定（`bcf_input.rs:119-133`）、index 的 BCF 分支判定（`index.rs:546-556`）；文本 index 走的是"收录 tail"（`index.rs:188-194`） | D9-3 三路对照表 |
| S6 | **seam 路径的错误判定比 body 路径宽松**。`FilteredJob::run_seam` 提前 `Ok(default())`（`filtered.rs:327-329`），而 `RecordPipeline::process`/`process_scanned` 对同样的解析失败报错（`record_pipeline.rs:165-166`、`:254-255`） | D9-4 |
| S7 | **指标口径与文档/spec 三处不一致**：peak 字节（×8）、seam 时间（漏计）、§7.7 的预算策略（spec 过时） | D9-5/6/7；01-review §1.6 实测 |
| S8 | **验收门只跑 release**，`debug_assert` 保护的不变量在门下不可见；一个红测试存活两轮并被误归因为"无关" | D9-8 |

---

## 2. 引擎终止语义：唯一错误出口【规】（修 D9-1 / D9-1b / D9-1c）

### 2.1 唤醒不变量

> **规约 W1**：`run_parallel` 里每一个可能无限期阻塞的线程，
> 都必须存在一个**不依赖其它线程正常返回**的唤醒来源。

现状违反 W1：台账（01-review §1.1）里 B1/B3/B4/B5 的无条件唤醒来源只有
abort 一个，而 abort 只由 `FirstError::record` 或外部 `trigger` 触发。

> **规约 W2**：任何线程（reader / worker / reducer）产生的任何错误与任何 panic，
> 在离开该线程之前**必须先经 `FirstError::record`**。
> `record` 内部会 `trigger`（`abort.rs:160-166`），因此 W2 ⇒ W1。

W2 是本轮的核心规范。它把 §7.5"中止即时生效"从"机制存在"收紧为
"**没有旁路**"：`?` 直接出闭包就是旁路。

### 2.2 reducer 的错误出口

reducer 的错误出口统一为一个 helper（命名建议 `bail`），
`deliver!` 的 5 处 `?` 与热循环的 `sink.finish()?` 全部经它：

```rust
// 【规】reducer 的唯一错误出口：先 record（触发 abort，唤醒所有阻塞点），
// 再断开 result 通道（让 worker 的 send 立即失败），最后返回 Aborted。
// run() 会优先报告 record 进去的因（pipeline.rs:231-233），
// 所以调用方看到的是 sink/job 的原始错误，不是 Aborted。
fn bail(first: &FirstError, result_rx: &mut Option<Receiver<..>>, err: EngineError)
    -> EngineError
{
    first.record(err);          // W2：触发 abort → B1/B3/B4/B5 的 ① 成立
    result_rx.take();           // 加速项：B5 的 ③ 也成立，worker 立即 break
    EngineError::Aborted
}
```

要点【规】：

1. **`record` 必须在 `return` 之前**。顺序反了就退化成现状。
2. **返回 `Aborted` 而不是原始错误**。`run()` 的
   `if let Some(err) = FirstError::take_shared(&first) { return Err(err) }`
   （`:231-233`）保证调用方拿到的是 record 进去的**因**；
   reducer 返回的 `Aborted` 是**果**，会被因覆盖。这条既有语义不需要改。
3. **`drop(result_rx)` 是加速项而非必要条件**。W2 的 abort 已足够打破环
   （01-review §3.1），drop 只是让 worker 少做无用功——
   但它与 D9-2 配套时才生效，见 §3。
4. **`ticket_tx` 不需要显式归还**。abort 触发后 reader 在 B1 的 ① 分支 break，
   不再需要票。不要为了"归还票"去改票据语义。
5. **收编已有的三处 abort 出口**【据】：`:551`、`:561`、`:587` 已经是
   `drop(result_rx); return Err(Aborted)` 的正确形状（前提是 `first.failed()`
   或 signal 已触发，即 W2 已被别的线程满足），只是重复了三遍。
   它们应与 `bail` 共用同一个出口，差别只在于不需要 `record`
   （错误已被别人记录）。正常 drain 路径 `:592-594` 的注释
   "dropping the receiver unblocks any worker still sending, so the scope
   join below cannot deadlock"说明这个机制作者本来就清楚——
   D9-1 是**错误路径漏掉了它**，不是设计缺失，所以修法不需要新机制。

### 2.3 需要改的具体位置【据】

| 位置 | 现状 | 改为 | 会挂死？ |
|---|---|---|---|
| `pipeline.rs:534` | `stitcher.push_head(...)?` | 经 `bail`（seam 超 1 GiB 的错误在此产生，`:139-142`） | **是** |
| `:535` | `job.run_seam(...)?` | 经 `bail`（D9-1b，实测挂死） | **是** |
| `:536` | `sink.accept(seam)?` | 经 `bail` | **是** |
| `:541` | `accepted?`（body 的 accept） | 经 `bail`（ENOSPC / EPIPE / 未排序索引在此产生，实测挂死） | **是** |
| `:545` | `stitcher.push_tail(...)?` | 经 `bail` | **是** |
| `:599` | `sink.finish()?` | 经 `bail`（**仅为满足 W2**，见下方说明） | **否** |
| `:521-524` | 注释"the dropped senders/receivers unblock every other stage" | 【规】改写为实际语义：错误经 `FirstError` 触发 abort 唤醒各阶段；出错路径没有 drop 任何 sender | — |

**关于 `:599` 为什么不挂**【据】：到达 `:594` 的唯一途径是 `break 'outer`，
而它只在 `recv(result_rx)` 返回 `Err`（所有 `result_tx` clone 已析构）时发生。
`chunk_tx` 的唯一存活副本在 reader 手里（`:349` clone 进去、`:420` drop 本体），
所以 worker 在没有 abort 的情况下**只能**在 `chunk_rx` 断开后退出，
即 reader 已退出。因此 `sink.finish()` 时 reader 与 worker 都已结束，
scope join 必然完成，错误也能正常经 `outcome?`（`:234`）返回调用方。

它仍然要经 `bail` 的理由是 W2 的一致性，后果是具体但次要的：
`Abort` 是调用方与引擎共享的（`:222-224` 注释），
`FrameChain` 的等待者靠 `abort.is_aborted()` 逃生
（`bcf_input.rs:84`、`index.rs:540`、`frame_chain.rs` 的 wait 循环）。
`sink.finish()` 的错误不 record ⇒ 这个共享 handle 从不触发 ⇒
调用方侧任何以 `is_aborted()` 为条件的清理/退出路径看到的是"一切正常"。
不要在文档或 commit message 里把它写成挂死。

`deliver!` 是宏，改造时【荐】把宏体收进一个返回 `Result` 的闭包，
在闭包外统一 `if let Err(e) = step { return Err(bail(...)) }`，
避免 5 处各写一遍 record 逻辑（重复即漏改风险）。

### 2.4 reducer panic

> 【规】reducer 的热循环与 drain 循环必须包在
> `catch_unwind(AssertUnwindSafe(...))` 里，与 worker 侧（`:450`）对称。
> unwind 时 `first.record(...)` 后走 §2.2 的出口。

错误变体【决】：`EngineError` 当前只有
`WorkerPanicked { chunk, message }`（`error.rs:19-25`）。两个选项：

- **A（推荐）**：新增 `ReducerPanicked { chunk: ChunkId, message: String }`。
  语义更准（reducer panic 与 worker panic 的可恢复性、定位含义都不同：
  前者来自 sink/输出层，后者来自 inflate/job），且 `#[non_exhaustive]`
  已声明（`error.rs:6`），新增变体不是破坏性变更。
- **B**：复用 `WorkerPanicked`，`chunk` 填正在交付的 id。
  改动小，但错误信息会误导（"worker panicked" 而 panic 的是 reducer/sink）。

两者都必须保证：panic 消息**逐字保留**（沿用 `panic_message`，`:748-756`），
且 `run` 返回它而不是裸 `Aborted`。

### 2.5 顺序路径

`run_sequential`（`:627-727`）跑在调用方线程上，没有 scope、没有通道，
错误直接 `?` 返回即可——**不受 W1/W2 约束**。但两件事必须与并行路径一致：

1. panic 行为【荐】：同样包 `catch_unwind`，返回同一个 `ReducerPanicked`，
   使"输出与错误对线程数一致"（§7.6 规则 5）覆盖错误类型，不只是输出字节；
2. §4 的截断检查必须同样落地（`:705-710`）——这是 D9-3 的一半。

---

## 3. 通道发送结果的语义【规】（修 D9-2）

> **规约 C1**：`select!` 里 `send(tx, msg)` 分支的返回值**必须**反映发送结果。
> `SendError`（对端已全部 drop）必须与 abort 同等对待：break 出循环。

两处改动：

| 位置 | 现状 | 改为 |
|---|---|---|
| `pipeline.rs:498-507`（worker → result） | `send(result_tx, WorkerMsg{..}) -> _ => true` | `-> res => res.is_ok()` |
| `:407-410`（reader → chunk） | `send(chunk_tx, InFlightChunk{..}) -> _ => true` | `-> res => res.is_ok()` |

要点：

1. **`if !sent { break }` 的控制流不变**，只是 `sent` 的取值变正确。
   worker 的 `:508-510` 注释"Unsent message drops with its guard"仍然成立——
   `SendError` 会把 `WorkerMsg` 连同 `PermitGuard` 还回来，许可不泄漏
   （这正是 S2 drop guard 设计的收益，本轮不需要改）。
2. **reader 侧的 `InFlightChunk` 同理**：`SendError` 归还 chunk 与 guard，
   `:413-416` 的注释"the guard is still a local and drops here"要补上
   "或对端断开时由 `SendError` 归还"。
3. C1 与 W2 的关系【据】：C1 单独不足以修 D9-1（abort 没触发时
   `result_rx` 也没 drop，`send` 不会失败），W2 单独足以修 D9-1
   （abort 唤醒 B5 的 ①）。但**两者都要做**：W2 保证唤醒，
   C1 保证唤醒后 worker 不做无用功、并且不会把"没发出去"当成"发出去了"
   而继续消费 chunk。D9-2 的复核门测的就是这个上界。

---

## 4. 输入截断契约【规】（修 D9-3 / D9-8）

### 4.1 现状：三条路径三种判定【据】

| 路径 | 末条未换行记录 | 依据 |
|---|---|---|
| 文本 filter / stats | 静默丢弃 | `NewlineStitcher` 无 finish；`owned_body_range` 止于最后一个 `\n`（`chunk/mod.rs:246-251`）；drain 后只调 `sink.finish()` |
| `reed index`（文本 VCF） | **收录并建索引** | `index.rs:188-194`（解析 tail）、`:439-441`（push）、`:446-449`（`eof_voff` 收尾） |
| BCF（filter/stats/index） | 报错 | `bcf_input.rs:119-133`、`index.rs:546-556` |

实测：50 条记录去掉末行 `\n` → 交付 **49** 条（顺序与并行一致）。

### 4.2 `NewlineStitcher::finish`【规】

```rust
impl NewlineStitcher {
    /// 流结束后取出残留的 fragment（`None` = 输入以记录边界结束）。
    /// 调用它即消耗 stitcher：fragment 不再参与任何后续 push。
    pub fn finish(self) -> Option<Vec<u8>> { /* std::mem::take 语义 */ }
}
```

同时【规】收紧 `push_head(None)` 的契约，把 D9-8 的 `debug_assert` 从
"生产不可达的虚保护"变成"有明确语义的断言"：

- `push_head(None)` 的含义是"本 chunk 没有 head fragment"
  = "上一个 chunk 以行边界结束" = **不应有残留 fragment**；
- 该断言在生产不可达（01-review §4：`split_text` 的 `head: None` 只在
  chunk 0，或 `first_in_block_offset > 0`，而两个源都只给 chunk 0 设非零偏移）；
- **流结束不是 `push_head(None)`**，是 `finish()`。
  D9-8 那个红测试的最后一行正是把两者混淆了，必须改成：

```rust
// 现状（红，且语义错）：
//   assert!(stitcher.push_head(None).unwrap().is_none());
// 改为：
assert_eq!(stitcher.finish().unwrap(), b"xy");   // 残留 fragment 必须交出来
```

### 4.3 引擎的截断判定【决】

> 【规】两条路径（并行 `:592-600`、顺序 `:705-710`）都必须在 drain 之后、
> `sink.finish()` 之前调用 `stitcher.finish()`；有残留即按下述决策处理。

**这是用户可见的行为变更**（今天静默成功的文件将来会失败），需要拍板：

| 选项 | 行为 | 评价 |
|---|---|---|
| **A（推荐）** | `Err(EngineError::job("input ends inside a record (truncated file)"))` | 与 BCF 路径逐字一致（`bcf_input.rs:129`）；三条路径统一；符合 htslib 对截断 BGZF 的处理方向；缺点：以前"能出结果"的截断文件现在报错 |
| B | 把残留 fragment 当最后一条记录交付 + stderr 警告 | 不丢数据，但与 index 的 `eof_voff` 语义、与 BCF 路径都不一致；且"警告 + 继续"违反 v2 06 §4 纪律 3"由输入数据决定的条件一律返回错误" |
| C | 维持现状（丢弃） | 不可接受：filter 输出与 index 内容自相矛盾（§4.1） |

选 A 时的连带影响【据】：

1. `reed index` 对**文本** VCF 的截断输入也会开始报错
   （引擎级检查，`IndexJob` 无法绕过）。这是**期望的**——
   它使文本 index 与 BCF index（`index.rs:546-556`）行为一致；
2. `IndexJob` 解析 tail 的逻辑（`:188-194`）**不需要改**：
   它处理的是 chunk 内部的 tail（会被下一个 chunk 的 head 补完），
   与"流末尾的残留 fragment"是两回事。引擎的 `finish()` 检查只在
   **最后一个 chunk 之后**触发；
3. 需要一条 e2e 夹具：末行无 `\n` 的 BGZF VCF，三个命令都必须失败且消息一致。

### 4.4 真值门（补 §7.6 的盲区）【规】

现有确定性门只比较"路径之间"，D9-3 因此在两条路径上同样丢记录也全绿。
新增一条**与输入真值比较**的门：

> 对以 `\n` 结束的输入，`RunReport::sink_report`（或等价的记录计数）
> 必须等于对该输入做单线程串行扫描得到的记录数。

这条门与 §4.3 的选项无关（A/B/C 都要满足"不静默丢"），
应作为 `reed-engine/tests/pipeline.rs` 的常驻验收测试。

---

## 5. seam 与 chunk 体的错误一致性【规】（修 D9-4）

> **规约 E1**：同一条记录，无论落在 chunk 体内还是跨 chunk 接缝，
> 必须产生**同一个错误类型与同一条消息**。
> 记录的致命性与 chunk 对齐无关，因而与线程数、文件大小无关。

E1 是 §7.6 确定性契约在**错误行为**上的延伸：§7.6 规则 5 只承诺输出字节一致，
但"某条记录是致命错误还是被跳过"同样不能依赖 chunk 边界。

现状违反【据】：

```rust
// filtered.rs:327-329  —— seam 路径：静默跳过
if LineRef::parse(line).is_none() {
    return Ok(J::Partial::default());
}

// record_pipeline.rs:165-166 / :254-255  —— body 路径：报错
let Some(lr) = LineRef::parse(line) else {
    return Err(EngineError::job("record line has fewer than 9 columns"));
};
```

`LineRef::parse` 在少于 9 个 tab 时返回 `None`（`line.rs:19-31`），
两处判的是同一个条件，结论相反。

修法【规】：**删除 `filtered.rs:327-329` 的提前返回**，
让 seam 记录走 `scan_line`/`filter_line` → `pipeline.process`/`process_scanned`，
由它们报同一个错。删除后 `run_seam` 的两个"filter 丢弃了 seam 记录"分支
（`:332-335`、`:341-344`）**保留不动**——那两个是过滤器**正常**丢弃
（记录合法但不通过），与畸形记录是两回事，不要一起删。

【据】注意 `self.plan.is_none()` 的早退分支（`:323-325`）不受影响：
无过滤计划时直接转发给 inner job，inner 自己判定畸形。

---

## 6. 指标语义【规】（修 D9-5 / D9-6）

### 6.1 `in_flight_peak_bytes`

实测：`peak / max_compressed = 8.000`，恰为 `MEM_FACTOR`（`budget.rs:118`）。
`metrics.rs:13` 写的是 "Peak in-flight compressed bytes"。顺序路径恒为 0
（`pipeline.rs:715`）。

| 选项 | 内容 | 评价 |
|---|---|---|
| **A（推荐）** | 承认它是**估算足迹**：字段更名 `in_flight_peak_footprint_bytes`，rustdoc 写明"= Σ 在途 chunk 的 `MemBudget::estimate(compressed_len)`，即压缩字节 × `MEM_FACTOR`"。顺序路径填真实的单 chunk 足迹峰值（它同时只有一个 chunk 在途，峰值 = 该 run 的最大 `estimate`），不再恒 0 | 不改运行时行为、不改预算语义；`memory_budget_bounds_in_flight_bytes` 的断言（peak ≤ budget）本来就是按足迹口径写的（`pipeline.rs:389` 的 footprint 与 `budget.mem.bytes` 同口径），更名后名实相符 |
| B | 改代码让它真的存压缩字节 | 会让 `memory_budget_bounds_in_flight_bytes` 的门变弱（压缩字节 ≪ 预算，断言失去意义），且与 `BytePermits` 的许可口径脱钩。**不推荐** |

【规】无论选哪个，`metrics.rs:13-14` 的 rustdoc 必须与
`budget.rs:107-109`（MemBudget 的背压语义）、`:118`（MEM_FACTOR）三处口径一致，
并在 rustdoc 里交叉引用 `MEM_FACTOR`。

### 6.2 `job_time` 与 `reduce_time`

实测【据】：同输入、同 273 ms 墙钟，顺序 `job_time=273.12 ms`、
并行 `job_time=49.2 µs`（`reduce_time` 3.4 µs 也不含它）。
seam 工作量在并行路径的指标里**不存在**。

> 【规】`job_time` 的定义改为：**所有 `run_chunk` 与 `run_seam` 调用的累计时间，
> 不论执行线程是 worker 还是 reducer**。`reduce_time` 改为：
> **所有 `OrderedSink::accept` 调用的累计时间，包含 seam 记录的那一次**。

落地要求：

1. 并行路径给 `deliver!` 里的 `job.run_seam`（`:535`）加计时，累加进
   `job_ns`（与 worker 侧 `:480` 同一个计数器）；
2. seam 的 `sink.accept`（`:536`）纳入 `reduce_ns` 的计时区间
   （当前区间是 `:538-540`，只覆盖 body 的 accept）；
3. `metrics.rs:34-36` 的 rustdoc 删掉 "across all workers"，
   改为 "across all workers **and the reducer**（seam 在 reducer 上执行）"；
   `:22-23` 的 `reduce_time` 同理写明含 seam；
4. **两条路径必须同语义**：顺序路径已经把 `run_seam` 计入 `job_time`
   （`:678-680`），是对的，不动；改的是并行路径；
5. 可选【荐】：若不想让 `job_time` 混入 reducer 时间，则新增
   `seam_time` 字段（`job_time` 保持"worker only"，`seam_time` 单列），
   但必须在 rustdoc 与 §7.10 的说明里写清两者相加才是"全部 job 时间"。
   **不允许**维持现状（既不计入 `job_time` 也不计入任何其它字段）。

### 6.3 `worker_busy` 的一致性

`metrics.rs:11-12` 定义 `worker_busy` = "inflate + job across all workers"。
并行路径由 `worker_busy_ns`（`:474-475`）测得，顺序路径是
`inflate_time + job_time`（`:714`）。§6.2 落地后顺序路径的 `job_time`
含 seam 而 `worker_busy` 也随之含 seam——这与"worker busy"的字面义冲突。

> 【规】顺序路径的 `worker_busy` 明确定义为"同一份工作的串行耗时"，
> 在 rustdoc 里写明它在顺序模式下含 seam（因为顺序模式没有独立 reducer 线程）。
> 现有断言 `worker_busy >= inflate_time + job_time`
> （`tests/pipeline.rs:775-777`）在两条路径上都必须继续成立。

---

## 7. spec 回改清单【规】（修 D9-7）

`docs/spec/architecture.md` 是 v1 规范，按 `docs/README.md` 的规则
"后一轮只修订它明确点到的章节"。本轮点到两处：

| 位置 | 现文 | 改为 |
|---|---|---|
| `:838-841`（§7.7） | "引擎负责运行时监控（在途字节计数），超出时报 `EngineError::MemBudgetExceeded`，不会默默 OOM。" | "引擎负责运行时监控（在途**估算足迹**计数，= 压缩字节 × `MEM_FACTOR`）。预算以**背压**方式执行：超出时 reader 阻塞等待许可归还，**不报错**；`EngineError::MemBudgetExceeded` 只保留给'单个 chunk 的足迹就超过全部预算'这种无法继续的情况（依据 `docs/rounds/v2/03-architecture.md` §4.2，本轮 doc_v9 D9-7 回改）。不会默默 OOM。" |
| `:820-823`（§7.5） | "所有阻塞点都用 `select!` 同时等待数据和 abort，中止即时生效" | 追加一句：**"任何线程的错误与 panic 都必须先经 `FirstError` 记录（记录即触发 abort），不允许绕过它直接返回——否则其余阻塞点失去唯一唤醒来源（doc_v9 §2 规约 W2）。"** |

另外两处**只登记不改**：

- §7.6（`:826-836`）：规则 5 的"输出对线程数逐字节一致"不覆盖
  "与输入真值一致"，本轮以 02 §4.4 的真值门补，不改 spec 条文；
- §7.10（`:843-846`）：`StageMetrics` 的字段语义在 `metrics.rs` 的 rustdoc 里，
  spec 只说"记录忙碌时间、空闲时间、字节数、记录数"，本轮 §6 改 rustdoc 即可。

`ci/check-doc-drift.sh` 加一条守卫【规】：
spec §7.7 不得再出现"超出时报 `MemBudgetExceeded`"这类"超出即失败"表述
（对齐 README「一致性守卫」的做法：文档声明的事实应有一条守卫）。

---

## 8. 验收纪律【规】（修 D9-8）

v2 06-roadmap.md §4 纪律 1 是"测试红着不提交"。本轮发现它有一个盲区：
**门只跑 release**（`v2/06:10`、`v8/06:11`），而 `debug_assert!` 在 release 下
被编译掉。结果 `stitcher_reassembles_and_caps_fragments` 自 Y5 加入起
在 debug 下一直红、在 release 下一直绿，存活两轮无人发现。

> 【规】通用验收门必须包含**至少一次 debug profile 的
> `cargo test --workspace`**（即不带 `--release`）。
> 顺序建议：debug 先跑（快，且能抓 `debug_assert`），release 后跑（性能门需要）。

> 【规】`debug_assert!` 只能用于**可证明生产不可达**的不变量，
> 且必须在旁边写明"为什么不可达"（本轮 `pipeline.rs:157-160` 缺这句话，
> 以致测试作者以为 `push_head(None)` 可以表达"流结束"）。
> 由输入数据决定的条件一律返回错误——这条 v2 06 §4 纪律 3 已有，本轮重申并补
> "debug_assert 需附不可达理由"。

> 【规】挂死类缺陷的验收测试**必须自带时间上界断言**
> （`assert!(elapsed < Duration::from_secs(5))`），不能用
> `#[should_panic]` 或依赖 CI 超时——否则测试红的时候会把 CI 挂住，
> 而这正是 D9-1 的表现形式。01-review §5 的每个 D9-1 门都已按此写。

> 【荐】`ci/check-doc-drift.sh` 增一条守卫：
> `crates/*/src` 里新增的 `debug_assert!` 必须在同一行或上一行有
> `invariant:` 前缀的理由（现有代码已是这个风格，如
> `pipeline.rs:159`、`abort.rs:161`），把惯例固化成守卫。

---

## 9. 本轮不做什么

- 不改票据窗口、许可 drop guard、`PermitGuard` 的析构语义（v6_y5 S2 已验证，
  本轮实测 `in_flight_final_bytes == 0` 仍成立）；
- 不改 `FrameChain` 的 10 ms 轮询（v6_y5 §4 已登记的偏差，01-review §4）；
- 不改 rayon 三阶段的预算接线（v8 D8-1…D8-6 已落地）；
- 不把 `run_sequential` 合并进 `run_parallel`
  （`v2/03-architecture.md` §4.3 曾建议"删除 `run_sequential`"，
  但 v6_y5 S4 选择了"保留独立顺序路径 + 逐字节对拍门"并已落地；
  重开它超出缺陷修复轮的范围）；
- 不动 `IndexJob` 解析 chunk 内 tail 的逻辑（§4.3 要点 2）；
- 不重开 v6 的 Y 系列、v7 的 Z 系列、v8 的 P 系列。
