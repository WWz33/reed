# doc_v6 · Y5：引擎健壮性（spec，2026-10-11）

> 范围：`docs/rounds/v6/06-roadmap.md` 的 **Y5** 行，即 `docs/rounds/v6/03-architecture.md`
> §4 的 **A6-1…A6-4** 四项加一项计量收口。本文是落地规格；实现记录写进
> `docs/rounds/v6/STATUS.md` 的 Y5 段。
>
> 与 roadmap 附录的差异在 §3 逐条列出（吸收 / 改写 / 不做），实现按本文执行。

## 1. 现状（本树 `v0.3.0-21-g2d5299c`）

| 点 | 代码 | 缺陷 |
|---|---|---|
| 中止信号 | `reed-engine/src/abort.rs:35-38`：`is_aborted = !rx.is_empty()`；容量 `available_parallelism()+4`（`abort.rs:63-67`）；`trigger` 只灌通道 | token 被一个 `select!` 消费后信号回假；`--threads` > 核数时等待者多于 token，唤醒不完整 |
| 许可归还 | `pipeline.rs:355-368` worker 失败路径直接 `break`，手里的 chunk footprint（`pipeline.rs:326-327`）不 release | `in_flight` 与许可池双双泄漏；reader 可能卡在 `BytePermits::acquire` 的 10 ms 轮询上（`budget.rs:118-138`） |
| 轮询 | `BytePermits::acquire` 用 `wait_timeout(10ms)` 轮询中止；`FrameChain::wait_inner` 同样（`frame_chain.rs:185-188`） | 违反 `docs/spec/architecture.md:820`「中止即时生效，不再用 20 ms 轮询」 |
| 接缝上限 | 文本路径 `stitch_seam`/`finish_fragment`（`pipeline.rs:522-566`）无上限；BCF 有 1 GiB（`bcf_input.rs:26-28,181`） | 无换行的病态输入把 reducer 的 fragment 撑到文件大小 |
| 顺序模式 | `run()` 无条件 `run_parallel`（`pipeline.rs:155`）；`is_sequential` 无调用者（`budget.rs:46-50`）；规范 `architecture.md:804` 承诺 `threads ≤ 2` 退化顺序且逐字节相同 | 规范偏差（F6-6） |
| 计量 | `StageMetrics`（`metrics.rs:8-24`）无 reader/reducer 等待、无 inflate/job/deflate 拆分；deflate 落点在 job 的 `BlockWriter`（`block_writer.rs:177`） | 「reader 受限」与 deflate 占比靠推断（A6-4） |

## 2. 五项规格

### S1 中止：原子位是真值 + 按预算定容 + 许可唤醒（A6-1a）

**真值来源**：`AbortSignal.flag: Arc<AtomicBool>`。`is_aborted()` = `flag.load(Acquire)`。
token 只负责唤醒 `select!`，被消费不影响真值。

**容量**：`abort_channel_for(threads)` = `bounded(threads + 4)`，覆盖 reader + workers +
reducer + 余量。新增 `abort_for(&ThreadBudget) -> Abort`；保留无参 `abort()`（测试与
无预算调用者，容量按 `available_parallelism()+4` 不变）。

**即时唤醒（修 §7.5 的轮询）**：
- `BytePermits` 可被注册进 `AbortHandle`：`AbortHandle` 持
  `Mutex<Vec<Weak<BytePermits>>>`，`register(&self, &Arc<BytePermits>)` 在创建时调用一次。
- `trigger()`：`flag.store(true, Release)` → `while try_send` → 对每个 `Weak::upgrade`
  成功者 `notify_all()`。
- `BytePermits::acquire` 改 `cv.wait(guard)` 无超时等待；循环内先查 `flag`。
  许可释放的 `notify_all` 照旧。中止延迟 = 唤醒延迟（微秒级），不再有 10 ms 轮询。
- `FrameChain` 已有自己的 poison/finish/`wait_timeout(10ms)`。本轮**不改**：它必须同时
  等 abort、poison、finish 三个来源，改为纯 condvar 需要把三者的通知路径一起重排；
  10 ms 轮询只发生在跨 chunk 的等待者上（每个 chunk 一次，不是每记录一次），
  登记为已知偏差（§4）。中止时 chain 的等待者靠 `abort.is_aborted()` 在 ≤10 ms 内退出，
  满足 2 s 门。

**错误定位**：`FirstError` 语义不变（首错胜出）。`run()` 收口不变。

### S2 在途许可：drop guard 是唯一出口（A6-1b）

`pipeline.rs`：

```rust
struct PermitGuard {
    permits: Arc<BytePermits>,
    shared: Arc<Shared>,
    bytes: u64,
}
impl PermitGuard {
    fn new(permits: Arc<BytePermits>, shared: Arc<Shared>, bytes: u64) -> Self {
        shared.record_in_flight(bytes);
        Self { permits, shared, bytes }
    }
}
impl Drop for PermitGuard {
    fn drop(&mut self) {
        self.permits.release(self.bytes);
        self.shared.release_in_flight(self.bytes);
    }
}
```

- `InFlightChunk { chunk, guard }`（替代 `{ chunk, footprint }`）；`WorkerMsg { id, guard, .. }`。
- reader：`acquire` 成功后 `PermitGuard::new(...)`，删掉三处手写 `release`（早退、abort
  分支、send 失败分支）——`select!` 的 send 操作数惰性求值，abort 分支胜出时 guard
  仍在栈上，`break` 即 drop。
- worker：`let InFlightChunk { chunk, guard } = in_flight;` 一次解构；失败分支 `break`
  即归还；成功把 guard 移进 `WorkerMsg`。
- reducer：`accept` 后 `drop(msg)` 归还，删掉显式 `release_in_flight`/`permits.release`。
- abort 早退：`result_rx` drop 连带 ring 里所有 guard 释放，reader 不再卡 `acquire`。
- `Shared.in_flight` 提供 `in_flight_final_bytes()`，作为门断言（S5、§5 门表）。

**采用 roadmap 附录的 permits 注册，采用本文的 `Arc` 形态**（附录是 `PermitGuard<'a>`）。
理由：guard 要进 `WorkerMsg<P>`、`ReorderBuffer`、channel，`'a` 会沿这些类型扩散；
每 chunk 两次原子增减（十万级 chunk）可忽略。偏差登记见 §3。

### S3 接缝上限：`NewlineStitcher` + 共享上限（A6-2）

- `pub const MAX_RECORD_BYTES: usize = 1 << 30;` 上移到 `reed-engine::pipeline`；
  `crates/reed/src/bcf_input.rs` 的私有常量改为 `use` 引用（`bcf_input.rs:181` 同步）。
- `NewlineStitcher`（`pipeline.rs`，pub）：

```rust
#[derive(Default)]
pub struct NewlineStitcher { fragment: Vec<u8> }
impl NewlineStitcher {
    fn grow(&mut self, bytes: &[u8]) -> Result<(), EngineError> { /* 超限 → EngineError::job */ }
    pub fn push_head(&mut self, head: Option<&[u8]>) -> Result<Option<Vec<u8>>, EngineError>;
    pub fn push_tail(&mut self, tail: Option<&[u8]>) -> Result<(), EngineError>;
}
```

- 并行 reducer 与顺序路径共用它；删 `stitch_seam`/`finish_fragment` 两个自由函数。
- 超限错误文案：`record exceeds {MAX_RECORD_BYTES} bytes across chunk seams`。
- **不抽 `SeamStitcher` trait**（§3）。

### S4 顺序模式（A6-3）

`run()` 按 `budget.is_sequential()` 分派 `run_sequential` / `run_parallel`。

```rust
fn run_sequential<J, S>(source, job, sink, budget, abort: &AbortSignal)
    -> Result<(u64, u64, StageMetrics), EngineError>
```

- 单线程循环：`next_chunk` → 单 chunk 超 `budget.mem` 报 `MemBudgetExceeded`（与并行同
  语义）→ `inflate_chunk` → 按 `framing()`：
  - `Newline`：`stitcher.push_head` 闭合则 `run_seam` + `accept`；`run_chunk` + `accept`；
    `push_tail`。
  - `Raw`（BCF）：整块交 `run_chunk`（链在 job 内）。
- 每步查 `abort.is_aborted()` 早退 `Aborted`；错误原样上抛（`run()` 仍走
  `FirstError` 收口，与并行一致）。
- 填 `reader_time`/`inflate_time`/`job_time`/`reduce_time`；`reader_wait`/`reducer_wait`
  为 0（无通道）；`in_flight_peak_bytes`/`in_flight_final_bytes` 为 0。
- 文档同步：`budget.rs:12-15` 注释与 `docs/spec/architecture.md:804-806` 的 F6-6 现状段
  改为「已实现」。

### S5 计量（A6-4）

`StageMetrics` 增字段（保持 `Copy`，`merge` 同步）：

| 字段 | 语义 | 测点 |
|---|---|---|
| `reader_wait` | reader 阻塞在票据 `select!`、`acquire`、发送 `select!` 的累计 | `pipeline.rs` reader 三处 |
| `reducer_wait` | reducer 阻塞在结果 `recv` 的累计 | reducer 底部 `select!` |
| `inflate_time` | 全部 worker 的 `inflate_chunk` 累计 | worker 闭包 |
| `job_time` | 全部 worker 的 `run_chunk`/`run_seam` 累计（deflate 是其子集） | worker 闭包 |
| `deflate_time` | 输出压缩累计，sink 上报 | 见下 |
| `in_flight_final_bytes` | 收尾时的在途字节（门断言用） | `Shared.in_flight` |

`worker_busy` 保持 = `inflate_time + job_time`（原语义：worker 忙碌总和）。

**deflate 的落点**：deflate 在 job 的 `BlockWriter`（`block_writer.rs:177`），不在引擎。
采用「随输出流到 sink」的形态，不动 job：

- `ChunkOutput` 增 `deflate_nanos: u64`；`BlockWriter` 累计，`finish_chunk` 时
  `std::mem::take` 出去（`block_writer.rs`）。
- `OrderedOutput` 累计 `deflate_nanos`（`accept_chunk` 处），提供 `deflate_time()`。
- `OrderedSink` trait 增 `fn deflate_time(&self) -> Duration { Duration::ZERO }`；
  `OutputSink`（`reed/src/filter/mod.rs:849`）转发 `self.out.deflate_time()`。
  引擎在 `sink.finish()` **之前**读一次（finish 消费 sink）。

**CLI 出口**：`reed-filter/src/funnel.rs` 的 `TimingReport` 增
`reader_wait_secs`/`reducer_wait_secs`/`inflate_secs`/`job_secs`/`deflate_secs`；
`timing_text` 增 `reader_wait_ms`/`reduce_wait_ms`/`inflate_ms`/`job_ms`/`deflate_ms`；
`reed/src/filter/mod.rs:666` 的构造补齐。stats/reports 只填结构体字段，渲染留 Y6。

## 3. 与 06-roadmap 附录的差异（逐条）

| 点 | 附录写法 | 本文 | 理由 |
|---|---|---|---|
| S1 中止 | flag + 容量 `budget.threads+4` + `trigger` `notify_all` permits | **吸收** | §7.5 要求不轮询；附录形态正确 |
| S2 guard | `PermitGuard<'a>{ &'a BytePermits, &'a AtomicU64 }` | 改 `Arc` 持有 | `'a` 会扩散进 `WorkerMsg`/`ReorderBuffer`/channel；原子开销可忽略 |
| S3 接缝 | 抽 `SeamStitcher` trait，BCF 的 `FrameChain` 也实现，引擎按 framing 选用 | 只做 `NewlineStitcher` + 共享上限 | BCF 走 `RecordFraming::Raw`，接缝在 job 内的链上，reducer 不碰；trait 只有一个实现类型。合一需要把链的所有权搬进引擎，超出 A6-2 的范围 |
| S4 顺序 | `run_sequential(&Abort) -> RunReport`，只写换行路径 | `&AbortSignal` + 返回 tuple；含 `Raw` 分支 | ① 顺序模式 BCF 输入仍走 Raw；② 直接返回 `RunReport` 会绕过 `run` 的 FirstError 收口 |
| S5 计量 | `ChunkJob::stage_nanos() -> [u64;3]` 回报 inflate/job/deflate | 引擎测 inflate/job；deflate 挂 `ChunkOutput` → sink | `inflate_chunk` 在引擎 worker（`pipeline.rs:330`），job 看不到；deflate 在 per-worker 的 `Scratch` 里，job 是 `&J` 共享，`stage_nanos` 报 deflate 需给每个 job 加共享原子。sink 上报零改 job，且覆盖所有走 `BlockWriter` 的出口 |

## 4. 明确不做（登记进 STATUS 偏差表）

1. `FrameChain::wait_inner` 的 10 ms 轮询：它要同时等 abort/poison/finish 三个来源，
   本轮不重排通知路径；中止延迟 ≤10 ms，满足门。
2. `SeamStitcher` trait 与 BCF 接缝合一：见 §3 S3。
3. `ChunkJob::stage_nanos()` 形态：见 §3 S5。
4. Y6（前后端收尾、C-5、联合 SFS、`--bed-file`/`--sites-file`）、Y7（CI fmt/doc/deny/MSRV、
   fuzz）不在本轮。

## 5. 门（与 06-roadmap Y5 行一致）

| 门 | 命令 / 断言 |
|---|---|
| 外部中止及时 | `threads = 8×核数` + 极小 `--mem`，外部 trigger 后 2 s 内返回 `Aborted`（新测试） |
| panic 后归零 | `worker_panic_becomes_worker_panicked` 之后 `in_flight_final_bytes == 0` |
| 跨 ≥10 chunk 的记录 | 10 B chunk × 100 B 行，fold == 串行扫描 |
| 无换行不 OOM | `NewlineStitcher` 单测：超限返回 `EngineError::Job`，不无限增长 |
| 顺序 = 并行 | `budget(1)`/`budget(2)` 与 `budget(8)` 的 fold 相同 |
| `--timing` 新字段 | `reader_wait_ms`/`reduce_wait_ms`/`inflate_ms`/`job_ms` 存在且非零（全量或子集） |
| 不变性 | gold 正文 md5 `7aedc59ac798729a2674e6e7ef19adfc`；`check_parity` PASS；子集 263,420 / 381,672；报告 8 件 md5 不变；三策略与跨格式 md5 不变 |
| 四件套 | `cargo test --release --workspace` 全绿；fmt、clippy（CI 命令）、rustdoc 退出 0；守卫 29/29 |
| 性能 | filter 全量 user 不回退；记录进 `bench/suites/e2e/results.tsv`（含批次起始负载与中位数） |

## 6. TODO（实现顺序）

见 [TODO.md](TODO.md)。批 1 引擎代码一个提交；批 2 测试与并行验证；批 3 文档收口与提交。
