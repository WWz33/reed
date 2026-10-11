# Y5 难点处参考实现（敲定版）

> 与 [README.md](README.md) §2 一一对应。这里给的是对照着改的骨架，不是最终成品；
> 类型名与现有代码一致。
>
> **注意 S1**：本文件采用 roadmap 附录的形态（`trigger` 唤醒注册的 `BytePermits`，
> `acquire` 无超时等待）。早前讨论里「保留 10 ms 轮询」的版本**已否决**——
> `docs/spec/architecture.md:820` 要求「中止即时生效，不再用轮询」。

## 1. 中止信号：原子位是真值，通道只负责唤醒，许可池被唤醒

`abort.rs`。`is_aborted` 改读 `AtomicBool` 后，token 被 `select!` 消费不再让信号回假；
容量按 `--threads` 定，超核数也够每个等待者拿一个 token；注册的许可池在 `trigger`
时被 `notify_all`，`acquire` 因此可以无超时等待。

```rust
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, Weak};

use crate::budget::BytePermits;
use crate::budget::ThreadBudget;

#[derive(Clone, Debug)]
pub struct AbortHandle {
    tx: Sender<()>,
    flag: Arc<AtomicBool>,
    /// Byte-permit pools blocked in `acquire`; notified on trigger so the
    /// abort is observed without a polling timeout (§7.5).
    permits: Arc<Mutex<Vec<Weak<BytePermits>>>>,
}

impl AbortHandle {
    /// Registers a permit pool to wake on abort. Call once per run, when the
    /// pool is created.
    pub fn register(&self, permits: &Arc<BytePermits>) {
        self.permits
            .lock()
            .expect("invariant: permit registry not poisoned")
            .push(Arc::downgrade(permits));
    }

    /// The flag is the source of truth; the channel only wakes `select!`
    /// waiters, so a full channel is fine. The registry wakes permit waiters.
    pub fn trigger(&self) {
        self.flag.store(true, Ordering::Release);
        while self.tx.try_send(()).is_ok() {}
        let mut registry = self
            .permits
            .lock()
            .expect("invariant: permit registry not poisoned");
        registry.retain(|weak| {
            let Some(permits) = weak.upgrade() else {
                return false;
            };
            permits.notify_abort();
            true
        });
    }
}

#[derive(Clone, Debug)]
pub struct AbortSignal {
    rx: Receiver<()>,
    flag: Arc<AtomicBool>,
}

impl AbortSignal {
    #[must_use]
    pub fn is_aborted(&self) -> bool {
        self.flag.load(Ordering::Acquire)
    }
    // raced_by_abort / receiver unchanged: `select!` still needs a token to wake.
}

/// Capacity covers every thread that may block on the signal in one run
/// (reader + workers + reducer + margin). `--threads` may exceed the cores.
#[must_use]
pub fn abort_channel_for(threads: usize) -> (AbortHandle, AbortSignal) {
    let (tx, rx) = crossbeam_channel::bounded::<()>(threads + 4);
    let flag = Arc::new(AtomicBool::new(false));
    (
        AbortHandle { tx, flag: Arc::clone(&flag), permits: Arc::new(Mutex::new(Vec::new())) },
        AbortSignal { rx, flag },
    )
}

#[must_use]
pub fn abort_for(budget: &ThreadBudget) -> Abort {
    let (handle, signal) = abort_channel_for(budget.threads);
    Abort { handle, signal }
}

/// Default for tests and callers without a budget.
#[must_use]
pub fn abort_channel() -> (AbortHandle, AbortSignal) {
    let n = std::thread::available_parallelism().map_or(8, std::num::NonZeroUsize::get);
    abort_channel_for(n)
}
```

`budget.rs` 的两处：

```rust
impl AbortHandle { /* ... */ }  // register 见上

impl BytePermits {
    /// Blocks until `bytes` permits are available (or abort fires); returns
    /// `false` on abort. The wait is a plain condvar wait, not a poll: the
    /// abort handle notifies this pool on `trigger`.
    pub fn acquire(&self, bytes: u64, abort: &crate::abort::AbortSignal) -> bool {
        let mut guard = self.available.lock().expect("invariant: permits not poisoned");
        loop {
            if abort.is_aborted() {
                return false;
            }
            if *guard >= bytes {
                *guard -= bytes;
                return true;
            }
            guard = self.cv.wait(guard).expect("invariant: permits not poisoned");
        }
    }

    /// Wakes every waiter so it re-checks the abort flag.
    pub fn notify_abort(&self) {
        self.cv.notify_all();
    }
}
```

注册点：`run_parallel` 创建 `permits` 后 `abort.handle.register(&permits)`（因此
`run_parallel` 收 `&Abort` 而不是 `&AbortSignal`）。

## 2. 许可归还：drop guard 是唯一出口

`pipeline.rs` 最难的一处。要点有两个：`select!` 的 send 操作数是惰性求值的，abort
分支赢了 guard 会留在栈上自然 drop；`InFlightChunk` 必须一次性解构，不能先拿走
`chunk` 再用 `guard`。

```rust
/// Owns one in-flight chunk's memory permit and its in-flight accounting.
/// Dropped, never explicitly released: on in-order delivery, on worker
/// failure, or when an aborted pipeline drops the pending message.
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

struct InFlightChunk {
    chunk: ChunkBytes,
    guard: PermitGuard,
}
```

Reader 循环里，`acquire` 之后到发送之前的所有手写 `permits.release` 全删：

```rust
let wait_started = Instant::now();
let acquired = permits.acquire(footprint, abort);
reader_wait_ns.fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::AcqRel);
if !acquired || first.failed() {
    break;
}
let guard = PermitGuard::new(Arc::clone(&permits), Arc::clone(&shared), footprint);
let sent = crossbeam_channel::select! {
    recv(abort.receiver()) -> _ => false,
    send(chunk_tx, InFlightChunk { chunk, guard }) -> _ => true,
};
if !sent || first.failed() {
    break; // abort won: `guard` is still a local, dropped here.
}
```

Worker：一次解构，失败即 drop，成功把 guard 移进消息。

```rust
let in_flight = match crossbeam_channel::select! {
    recv(abort.receiver()) -> _ => break,
    recv(chunk_rx) -> msg => match msg { Ok(x) => x, Err(_) => break },
};
let InFlightChunk { chunk, guard } = in_flight;
// ... inflate + job（见 §5）...
if let Some(err) = failure {
    first.record(err);
    break; // `guard` drops: permit + in_flight returned on the error path.
}
let sent = crossbeam_channel::select! {
    recv(abort.receiver()) -> _ => false,
    send(result_tx, WorkerMsg { id: chunk.id, guard, partial, head, tail }) -> _ => true,
};
if !sent {
    break; // unsent message drops with its guard.
}
```

Reducer：`msg` 持 guard，交付后才 drop，删掉 `release_in_flight`/`permits.release` 两行。

```rust
sink.accept(id, &mut msg.partial)?;   // in-order delivery
let tail = msg.tail.take();
drop(msg);                            // permits + in_flight return here
let _ = ticket_tx.send(());
stitcher.push_tail(tail.as_deref())?;
delivered += 1;
```

abort 早退时 `result_rx` 被 drop，ring 里所有 `WorkerMsg` 连带 guard 释放，reader 不再卡
`acquire`。`run_parallel` 收尾加 `in_flight_final_bytes: shared.in_flight.load(...)`，
panic 门就能断言 0。

## 3. 接缝上限：一个带界的状态机

`pipeline.rs`。parallel reducer 和 sequential 共用同一份。把 `MAX_RECORD_BYTES` 上移到这里，
`bcf_input.rs:28` 改成引用。

```rust
/// The largest record the engine reassembles across chunk seams. Shared with
/// the BCF frame chain, which enforces the same ceiling on carried prefixes.
pub const MAX_RECORD_BYTES: usize = 1 << 30;

#[derive(Default)]
pub struct NewlineStitcher {
    fragment: Vec<u8>,
}

impl NewlineStitcher {
    fn grow(&mut self, bytes: &[u8]) -> Result<(), EngineError> {
        if self.fragment.len() + bytes.len() > MAX_RECORD_BYTES {
            return Err(EngineError::job(format!(
                "record exceeds {MAX_RECORD_BYTES} bytes across chunk seams"
            )));
        }
        self.fragment.extend_from_slice(bytes);
        Ok(())
    }

    /// Absorbs a chunk's head fragment; `Some(line)` when the boundary line
    /// closes here. `None` head means the first chunk.
    pub fn push_head(&mut self, head: Option<&[u8]>) -> Result<Option<Vec<u8>>, EngineError> {
        let Some(head) = head else {
            debug_assert!(self.fragment.is_empty(), "invariant: fragment before the first chunk");
            return Ok(None);
        };
        self.grow(head)?;
        Ok(self.fragment.ends_with(b"\n").then(|| std::mem::take(&mut self.fragment)))
    }

    pub fn push_tail(&mut self, tail: Option<&[u8]>) -> Result<(), EngineError> {
        if let Some(tail) = tail {
            self.grow(tail)?;
        }
        Ok(())
    }
}
```

BCF 的 `FrameChain` 已有同等上限（`bcf_input.rs:181`），保持不动：它走
`RecordFraming::Raw`，接缝不经过 reducer，硬套一个 trait 只会多一个单实现抽象。

## 4. 顺序模式：复用 stitcher，分派在 `run`

`pipeline.rs`。早退与错误语义和并行保持一致，单 chunk 超预算照样报错。

```rust
fn run_sequential<J, S>(
    source: &mut dyn ChunkSource,
    job: &J,
    sink: S,
    budget: &ThreadBudget,
    abort: &AbortSignal,
) -> Result<(u64, u64, StageMetrics), EngineError>
where
    J: ChunkJob + 'static,
    S: OrderedSink<J::Partial> + 'static,
{
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut inflate_scratch = Vec::new();
    let mut scratch = J::Scratch::default();
    let mut seam_scratch = J::Scratch::default();
    let mut stitcher = NewlineStitcher::default();
    let mut sink = sink;
    let (mut reader_time, mut inflate_time, mut job_time, mut reduce_time) =
        (Duration::ZERO, Duration::ZERO, Duration::ZERO, Duration::ZERO);
    let (mut chunks_read, mut delivered) = (0_u64, 0_u64);

    loop {
        if abort.is_aborted() {
            return Err(EngineError::Aborted);
        }
        let t = Instant::now();
        let chunk = source.next_chunk()?;
        reader_time += t.elapsed();
        let Some(chunk) = chunk else { break };
        chunks_read += 1;

        let footprint = budget.mem.estimate(chunk.compressed.len() as u64);
        if footprint > budget.mem.bytes {
            return Err(EngineError::MemBudgetExceeded {
                peak_bytes: footprint,
                budget_bytes: budget.mem.bytes,
            });
        }

        let t = Instant::now();
        let view = inflate_chunk(&chunk, &mut inflater, &mut inflate_scratch)?;
        inflate_time += t.elapsed();

        match job.framing() {
            RecordFraming::Newline => {
                if let Some(seam) = stitcher.push_head(view.head_bytes(view.text))? {
                    let t = Instant::now();
                    let mut partial = job.run_seam(chunk.id, &seam, &mut seam_scratch)?;
                    job_time += t.elapsed();
                    let t = Instant::now();
                    sink.accept(chunk.id, &mut partial)?;
                    reduce_time += t.elapsed();
                }
                let t = Instant::now();
                let mut partial = job.run_chunk(&chunk, view.body_bytes(view.text), &mut scratch)?;
                job_time += t.elapsed();
                let t = Instant::now();
                sink.accept(chunk.id, &mut partial)?;
                reduce_time += t.elapsed();
                stitcher.push_tail(view.tail_bytes(view.text))?;
            }
            // The job walks the length prefixes; no seam assembly here.
            RecordFraming::Raw => {
                let t = Instant::now();
                let mut partial = job.run_chunk(&chunk, view.text, &mut scratch)?;
                job_time += t.elapsed();
                let t = Instant::now();
                sink.accept(chunk.id, &mut partial)?;
                reduce_time += t.elapsed();
            }
        }
        delivered += 1;
    }

    let deflate_time = sink.deflate_time();
    let sink_report = sink.finish()?;
    Ok((delivered, sink_report, StageMetrics {
        n_chunks_read: chunks_read,
        worker_busy: inflate_time + job_time,
        in_flight_peak_bytes: 0,
        reader_time,
        reduce_time,
        chain_wait: job.wait_time(),
        reader_wait: Duration::ZERO,
        reducer_wait: Duration::ZERO,
        inflate_time,
        job_time,
        deflate_time,
        in_flight_final_bytes: 0,
    }))
}
```

分派只改 `run` 一处；两条路径共用同一个 `FirstError` 收口。

```rust
let outcome = if budget.is_sequential() {
    run_sequential(source, job, sink, budget, &abort.signal)
} else {
    run_parallel(source, job, sink, budget, &first, abort)
};
```

## 5. 计量：三段拆分 + 两处等待

`metrics.rs` 加字段（仍 `Copy`）：

```rust
pub struct StageMetrics {
    pub n_chunks_read: u64,
    pub worker_busy: Duration,        // = inflate_time + job_time（worker 侧）
    pub in_flight_peak_bytes: u64,
    pub reader_time: Duration,        // next_chunk
    pub reduce_time: Duration,        // accept
    pub chain_wait: Duration,
    pub reader_wait: Duration,        // blocked on ticket/permits/send
    pub reducer_wait: Duration,       // blocked on recv
    pub inflate_time: Duration,       // engine-measured
    pub job_time: Duration,           // engine-measured (deflate is a subset)
    pub deflate_time: Duration,       // sink-reported, rides ChunkOutput
    pub in_flight_final_bytes: u64,
}
```

worker 闭包里拆两段（`inflate` 由引擎测，`job` 是 `run_chunk` 墙钟）：

```rust
let inflate_started = Instant::now();
let view = inflate_chunk(&chunk, &mut inflater, &mut inflate_scratch)?;
let inflate_elapsed = inflate_started.elapsed();
// framing copy ...
let job_started = Instant::now();
let partial = job.run_chunk(&chunk, body, &mut scratch)?;
let job_elapsed = job_started.elapsed();
// catch_unwind 外层再把两段累加进原子
```

reader 的三处阻塞（票据 select、`acquire`、send select）累加进 `reader_wait_ns`；
reducer 底部的 `select! recv` 包进 `reducer_wait_ns`。

deflate 的测点在 `BlockWriter::cut`（`block_writer.rs:177`），让它随 `ChunkOutput` 流到
sink，不需要 job 侧原子：

```rust
// block_writer.rs
pub struct ChunkOutput { pub blocks: Vec<OutBlock>, pub raw: Vec<u8>, pub deflate_nanos: u64 }
pub struct BlockWriter { /* ... */ deflate_nanos: u64 }

// cut()
let t = Instant::now();
deflater.deflate_block(&self.pending[..cut], &mut compressed)
    .map_err(|e| IoError::Io(format!("bgzf deflate: {e}")))?;
self.deflate_nanos += t.elapsed().as_nanos() as u64;

// finish_chunk()
Ok(ChunkOutput {
    blocks: std::mem::take(&mut self.blocks),
    raw: std::mem::take(&mut self.raw),
    deflate_nanos: std::mem::take(&mut self.deflate_nanos),
})
```

```rust
// ordered_output.rs
pub struct OrderedOutput { /* ... */ deflate_nanos: u64 }

pub fn accept_chunk(&mut self, chunk: &ChunkOutput) -> Result<(), IoError> {
    self.deflate_nanos += chunk.deflate_nanos;
    // ... existing ...
}
pub fn deflate_time(&self) -> Duration { Duration::from_nanos(self.deflate_nanos) }
```

```rust
// pipeline.rs trait，默认 0，与 wait_time 同一形态
pub trait OrderedSink<P>: Send {
    // accept / finish ...
    /// Cumulative output-compression time this sink observed (A6-4).
    fn deflate_time(&self) -> Duration { Duration::ZERO }
}

// reed/src/filter/mod.rs
impl OrderedSink<FilterPartial> for OutputSink {
    fn deflate_time(&self) -> Duration { self.out.deflate_time() }
    // ...
}
```

引擎在 `sink.finish()` 之前读一次（finish 会消费 sink），顺序不能颠倒：

```rust
let deflate_time = sink.deflate_time();
let sink_report = sink.finish()?;
```

`funnel.rs` 的 `TimingReport` 加 `reader_wait_secs`/`reducer_wait_secs`/`inflate_secs`/
`job_secs`/`deflate_secs`，`timing_text` 加同名毫秒行，`filter/mod.rs:666` 的构造补字段。
stats/reports 只补结构体字段，渲染留 Y6。

## 落地顺序与坑

1. `abort.rs` + `abort_for` 调用点（`filter/mod.rs:589,637`、`index.rs:511`、
   `stats/job.rs:1887`、`reports.rs:798`、`matrix.rs:889`）。
2. `PermitGuard`，先删 reader/reducer 的手写 release 再编译，编译器会指出漏点。
3. `ChunkOutput` 加字段，唯一构造点是 `finish_chunk` 和 `block_writer.rs` 的测试字面量，
   一并改。
4. `NewlineStitcher` + `MAX_RECORD_BYTES` 上移，`bcf_input.rs` 改引用。
5. `run_sequential` + 分派。
6. metrics 与 sink 管道。

两个易错点：`select!` 的 send 操作数惰性求值，abort 分支胜出时 guard 仍在栈上，靠
`break` drop，不要手写 release；`InFlightChunk` 必须
`let InFlightChunk { chunk, guard } = in_flight;` 一次解构，否则 partial move 编译不过。
