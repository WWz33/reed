# 06 路线（Y0–Y7）

> 通用门沿用各轮：测试全绿、gold 对拍（正文 md5 `7aedc59ac798729a2674e6e7ef19adfc` + `check_parity` PASS）、子集 263,420 / 381,672、三策略与跨格式 md5、跨线程 sha256、性能不回退并入 `bench/suites/e2e/results.tsv`（含批次开始负载与中位数，04 §9）、同 PR 删旧实现、更新 `docs/rounds/v6/STATUS.md`（每行引用测试名或可复现命令）。难点都带代码示意（02 §2.2、03 §7–9）。

```
Y0 lint 门变绿 + 漂移/死代码归零 ──► Y1 全局分析矩阵化（剪枝/衰减/--mind 共用 pass0）──► Y6 前后端收尾 + 数据即结果
                                   ├─► Y2 I/O 热路径            ├─► Y3 过滤热路径 ──► Y4 -Ob/BCF 记录模型
                                   └─► Y5 引擎健壮性            Y7 CI 与测试（随 Y0 起，贯穿）
```

| 步 | 内容 | 门 |
|---|---|---|
| **Y0** | 05 E6-1（修 `tail.rs:75-76`；测试目标 pedantic 修或限定 `-D warnings` 到 lib/bin）、E6-7（05 §4 八处漂移，含 `ld_prune.rs:7-9`）、E6-9（删 `field_plan.rs`、`CountTable`/`count_site_group`、`math::hypergeom`、`ChunkBytesExport`、`from_half_memory_for_tests`；`is_sequential` 随 Y5 决定）、03 A6-23（`ContigId` 统一 newtype）、v5 STATUS 的 clippy 陈述更正；**05 §7–§8（交付性与命名）**：`LICENSE-MIT`/`LICENSE-APACHE`、`v0.3.0` tag + CHANGELOG、`--version` 带 commit/日期（`build.rs`）、README 加 Installation/Citation/输入要求/退出码、`pixy` 词汇出生产代码（`writers/windows.rs`、`WindowLayoutOptions`、`window_kinds`）、测试文件改名（`m7_fst`/`qc_merge_laws`/`individual_joint`）、`cli.rs`→`cli/mod.rs`、`ponytail`→`TODO(#)`、4 个公开参数名去缩写 | `cargo clippy --release --workspace --all-targets -- -D warnings` 退出 0；`RUSTDOCFLAGS=-D warnings cargo doc --workspace --no-deps` 通过；守卫 + 新增"`allow(unsafe_code)` 数量"随删除同步；`ls LICENSE-MIT LICENSE-APACHE`；`reed --version` 含 commit；`! grep -rniE "\bpixy" crates/*/src`（测试除外）；`! grep -rn ponytail crates` |
| **Y1** | 02 §2.2 / 03 §7：`ContigMatrix` + `MatrixJob`/`MatrixSink`（reducer 只拼接）；`prune_all` contig 并行（严格 plink2 顺序语义，含 kb 窗口按 `LdPruneNextWindow` 对齐，F6-10）；`ld_decay` 位置区块并行；`--mind` 从矩阵平面按列 popcount；删 `MissingJob`/`SketchJob`/`LdPruneJob`；`--ld-decay` 不再要求 `--pop`（F6-9）；内存估计进 `--memory`，超预算按 contig 分段 | 剪枝集合与 plink2 逐 ID 相同：夹具、区域（277）、全量（235,403）、kb 窗口（新增 `50kb` 对拍，对称差 0）；衰减 r² 与现实现逐对相同（区域 427,152 对）；`--mind` 输出与现 sketch/two-pass 逐字节相同；1/8/56 线程相同；**剪枝 ≤ 3 s、衰减 ≤ 10 s**；在途峰值 ≤ 预算 |
| **Y2** | 03 A6-5/6/8/10/11/12 + A6-7：deflate/inflate 常驻缓冲、`memrchr`、`GtTable` 进程级 `Arc`、`format_g6` 栈缓冲、`Locus` 不分配、索引 contig 名缓存；**流式源偏移漂移修复 + 流式 vs mmap 元数据等价测试**；reports 双峰诊断（04 §7） | 全部输出 md5 不变；filter 全量 user **≤ 420 s**（−40）；stats ≤ 68 s；`cat x.vcf.gz \| reed stats` 与文件路径输出相同；reports 3 次极差 ≤ 3% |
| **Y3** | 03 A6-15…20 + §9：`stale` 预算进计划、HWE 走 `HweCache`、列式掩码（整数比较出位图）、一记录一次 INFO 片段索引供表达式与 head 共用、谓词预编译、`emit_survivor` 合并两套闭包、`masked` 整字清零 | 子集 263,420 / 381,672 不变；gold md5 不变；`stats_site_filters_match_the_filter_then_stats_composition` 通过；filter 全量 user 再 −30 s；`--set-gts .` 表达式路径无每样本分配（基准） |
| **Y4** | 03 §8：`Record` 模型 + `RecordSource`/`RecordSink`（VcfText / Bcf）；`-Ob` 与 BCF 输入走模型，未修改 FORMAT 列原字节拷贝；A6-14 两处渲染差异对齐 | `-Ob` 读回正文 ≡ `-Oz` 正文；noodles 编码器 oracle 字节相等；全量 `-Ob` user **≤ 470 s**（现 521）；BCF→BCF 全量 ≤ 9.5 s；VCF4.4 相位与 `-0.0` 夹具与 `bcftools view -H` 一致 |
| **Y5** | 03 A6-1…4：中止通道按 `budget.threads` 定容 + `AtomicBool`；在途 chunk drop guard 释放许可；`MAX_RECORD_BYTES` 共享 + `SeamStitcher`；`StageMetrics` 加 `reader_wait`/`reducer_wait` 与 inflate/job/deflate 三段；顺序模式实现或改规范；计数内核若保留则 release 断言 + 带检查构造器（否则已在 Y0 删） | 新测试：threads = 8× 核数 + 小预算下外部中止 2 s 内返回 `Aborted`；worker panic 后 `in_flight == 0`；跨 ≥ 10 chunk 的记录重建；无换行输入报错；`--timing` 新字段非零；顺序路径 sha256 与并行相同 |
| **Y6** | 前后端收尾与数据即结果：`plan::compile`/`assign_groups`/`pixy_kinds` 移入 `run_filter`/`run_stats`（`translate/stats.rs:330,438,507,660`）；单一 `reed::Error` + `kind()` 定退出码；stats/index 往返测试；`StatsSink` 只产 `WindowRow`/partial，文本由 `reed::writers`；`stats/job.rs` 拆 `request/run/sink`（1,936 行）；C-5 拆 `individual.rs`；联合 SFS 接 CLI；`--bed-file`/`--sites-file`（`WindowGrid::Regions`） | 输出 md5 不变；`cargo doc -p reed` 公开项 = 08 §6.2 清单；`cli_contract` 三条命令往返；联合谱与 easySFS 二维对拍；BED 窗口与 pixy `--bed_file` 逐列一致 |
| **Y7** | 05 E6-2…6、E6-8：CI 加 fmt / doc / deny / MSRV；基准进 CI（`run_baseline.sh`，回退 > 5% 阻断）；引擎边界与两种源的测试；索引内容测试；`cargo-fuzz` 四个目标；公开树守卫或公开 `docs/spec/` | CI 全绿且含基准行；fuzz 目标各跑 10 分钟无崩溃 |

顺序：Y0 当天；Y1 单独一个 PR（它改变三条命令的 pass0，是本轮的主项）；Y2/Y3/Y5 可并行、各自独立 PR；Y4 在 Y3 之后；Y6 的六项各自 PR；Y7 从 Y0 起随每个 PR 补。

---

## 附：Y5 代码示意（非最终代码；签名取自 `reed-engine/src/{abort,budget,pipeline}.rs` 现状）

现状：`AbortSignal::is_aborted` = "通道非空"（`abort.rs:35-38`），容量按 `available_parallelism()+4`（`abort.rs:63-67`）；`BytePermits::acquire` 每 10 ms 轮询一次中止（`budget.rs:124-143`）；worker 失败路径 `break` 时不归还许可（`pipeline.rs:355-368`）；文本接缝 `stitch_seam`/`finish_fragment`（`pipeline.rs:533-574`）无上限，BCF 的上限在 `bcf_input.rs:28`。

```rust
// reed-engine/src/abort.rs —— 原子位 + 按预算定容的唤醒通道
pub struct AbortSignal { flag: Arc<AtomicBool>, rx: Receiver<()> }
impl AbortSignal {
    pub fn is_aborted(&self) -> bool { self.flag.load(Ordering::Acquire) }   // 不再依赖"通道里还有 token"
    pub const fn receiver(&self) -> &Receiver<()> { &self.rx }               // select! 仍靠通道唤醒
}
impl AbortHandle {
    pub fn trigger(&self) {
        self.flag.store(true, Ordering::Release);
        while self.tx.try_send(()).is_ok() {}
        for p in self.permits.iter().filter_map(Weak::upgrade) { p.notify_all(); }   // 唤醒卡在许可上的 reader
    }
}
/// 容量 = 本次运行可能阻塞在信号上的线程数（reader + workers + reducer + 余量），来自 ThreadBudget，不是核数。
pub fn abort_channel_for(budget: &ThreadBudget) -> (AbortHandle, AbortSignal) {
    let (tx, rx) = crossbeam_channel::bounded::<()>(budget.threads + 4);
    let flag = Arc::new(AtomicBool::new(false));
    (AbortHandle { tx, flag: Arc::clone(&flag), permits: Vec::new() }, AbortSignal { flag, rx })
}

// reed-engine/src/budget.rs —— 许可改为 condvar 等待（无轮询），并用 drop guard 保证归还
impl BytePermits {
    pub fn acquire(&self, bytes: u64, abort: &AbortSignal) -> bool {
        let mut g = self.available.lock().expect("invariant: permits not poisoned");
        loop {
            if abort.is_aborted() { return false; }
            if *g >= bytes { *g -= bytes; return true; }
            g = self.cv.wait(g).expect("invariant: permits not poisoned");   // trigger() 会 notify_all
        }
    }
}
pub struct PermitGuard<'a> { permits: &'a BytePermits, in_flight: &'a AtomicU64, bytes: u64 }
impl Drop for PermitGuard<'_> {
    fn drop(&mut self) { self.permits.release(self.bytes); self.in_flight.fetch_sub(self.bytes, Ordering::AcqRel); }
}
// pipeline.rs：InFlightChunk { chunk, guard: PermitGuard }；WorkerMsg 把 guard 一路带到 reducer，
// 交付后 drop（语义不变：按序交付才还许可）；worker 出错/panic 时 guard 随栈展开归还，in_flight 归零。

// reed-engine/src/pipeline.rs —— 接缝拼接的唯一抽象 + 共享上限
pub const MAX_RECORD_BYTES: usize = 1 << 30;        // bcf_input.rs:28 改为引用这里
pub trait SeamStitcher: Send {
    /// 吞下下一 chunk 的头片段；拼成完整记录时返回它。
    fn push_head(&mut self, head: &[u8]) -> Result<Option<Vec<u8>>, EngineError>;
    /// 吞下本 chunk 的尾片段。
    fn push_tail(&mut self, tail: &[u8]) -> Result<(), EngineError>;
}
pub struct NewlineStitcher { fragment: Vec<u8> }
impl NewlineStitcher {
    fn grow(&mut self, bytes: &[u8]) -> Result<(), EngineError> {
        if self.fragment.len() + bytes.len() > MAX_RECORD_BYTES {
            return Err(EngineError::job(format!("record exceeds {MAX_RECORD_BYTES} bytes across chunk seams")));
        }
        self.fragment.extend_from_slice(bytes);
        Ok(())
    }
}
impl SeamStitcher for NewlineStitcher {
    fn push_head(&mut self, head: &[u8]) -> Result<Option<Vec<u8>>, EngineError> {
        self.grow(head)?;
        Ok(self.fragment.ends_with(b"\n").then(|| std::mem::take(&mut self.fragment)))
    }
    fn push_tail(&mut self, tail: &[u8]) -> Result<(), EngineError> { self.grow(tail) }
}
// BCF 的长度前缀链实现同一 trait（FrameChain 内部不变）；引擎按 ChunkJob::framing() 选用其一，
// 删除 pipeline.rs 的 stitch_seam/finish_fragment 两套手写逻辑。

// reed-engine/src/metrics.rs —— 补齐计量，让"reader 受限 / deflate 占比"成为测量而不是推断
pub struct StageMetrics {
    /* 现有字段 */
    pub reader_wait: Duration,     // reader 阻塞在票据 + 许可上的累计
    pub reducer_wait: Duration,    // reducer 阻塞在 recv 上的累计
    pub job_stages: [Duration; 3], // ChunkJob::stage_nanos() 的和：inflate / job / deflate
}
pub trait ChunkJob: Sync { /* 现有 */ fn stage_nanos(&self) -> [u64; 3] { [0; 3] } }

// 顺序模式（规范 docs/spec/architecture.md:804）：
pub fn run_sequential<J, S>(source: &mut dyn ChunkSource, job: &J, mut sink: S, abort: &Abort) -> Result<RunReport, EngineError> {
    let mut stitcher = NewlineStitcher::default(); let mut scratch = J::Scratch::default();
    while let Some(chunk) = source.next_chunk()? {
        if abort.signal.is_aborted() { return Err(EngineError::Aborted); }
        let view = inflate_chunk(&chunk, &mut inflater, &mut buf)?;
        if let Some(seam) = view.head_bytes(view.text).map(|h| stitcher.push_head(h)).transpose()?.flatten() {
            let mut p = job.run_seam(chunk.id, &seam, &mut scratch)?; sink.accept(chunk.id, &mut p)?;
        }
        let mut p = job.run_chunk(&chunk, view.body_bytes(view.text), &mut scratch)?; sink.accept(chunk.id, &mut p)?;
        if let Some(t) = view.tail_bytes(view.text) { stitcher.push_tail(t)?; }
    }
    /* sink.finish() … */
}
// run()：budget.is_sequential() ? run_sequential : run_parallel；测试断言顺序路径与 1/N worker 的 sha256 相同。
```

门（与 Y5 行一致）：`threads = 8 × 核数` + 极小 `--memory` 下外部中止 2 s 内返回 `Aborted`；`worker_panic_becomes_worker_panicked` 之后 `in_flight == 0`；10 字节 chunk × 100 字节行（记录跨 ≥ 10 chunk）精确重建；无换行输入报错而非耗尽内存；`--timing` 的新字段非零；顺序路径输出 sha256 与并行相同。
