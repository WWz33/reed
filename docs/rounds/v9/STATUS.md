# doc_v9 实施状态

> 本页记录 `01-review.md` §3 缺陷清单与 `06-roadmap.md` Q0–Q3 的实施状态。
> 每条附证据（测试名 / 守卫 / 变更点）。

## 实施状态（Q0–Q3）

| 项 | 严重度 | 状态 | 证据 |
|---|---|---|---|
| D9-1 reducer 侧错误挂死 | **P0** | **已修** | `pipeline.rs`：reducer 的唯一错误出口 `bail(first, err)`（`record` → 返回 `Aborted`，`run` 优先报告记录的因）；`deliver!` 的 5 处 `?` 与 `sink.finish()` 都经它。测试 `reducer_sink_error_returns_instead_of_hanging`、`reducer_first_chunk_sink_error_returns`（修复前 2788 chunk × `budget(4)` 120 s 不返回，修复后 <5 s 返回 sink 的原始错误） |
| D9-1b `run_seam` 错误挂死 | **P0** | **已修** | 同一出口覆盖；`reducer_seam_error_returns_instead_of_hanging` |
| D9-1c reducer panic 挂死 | **P0** | **已修** | 整个 reducer 体（含 `accept`/`run_seam`/`finish`）包 `catch_unwind(AssertUnwindSafe(..))`，unwind 转 `EngineError::ReducerPanicked { chunk, message }`（`error.rs` 新增变体，`#[non_exhaustive]` 非破坏性）；测试 `reducer_panic_returns_instead_of_hanging`、`reducer_finish_error_returns_instead_of_hanging` |
| D9-2 发送失败当成功 | P1 | **已修** | `worker→result`（原 `:498`）与 `reader→chunk`（原 `:407`）的 `-> _ => true` 改为 `-> res.is_ok()`；守卫 "pipeline sends check their result" |
| D9-3 末条未换行记录：丢弃/收录/报错三路不一 | **P0** | **已修** | `NewlineStitcher::finish(self) -> Option<Vec<u8>>`；并行与顺序两条路径在 drain 后检查，残留即 `EngineError::job(TRUNCATED_INPUT)`（消息与 BCF 路径逐字一致）。测试 `unterminated_final_record_is_rejected`（并行/顺序 × 两种 chunk 划分）、`terminated_input_delivers_every_record`（防过度收紧，兼真值门）、`stitcher_reassembles_and_caps_fragments`（红测试转绿） |
| D9-4 `FilteredJob::run_seam` 静默跳过畸形记录 | P1 | **已修** | 删除 `filtered.rs` 的 `LineRef::parse().is_none() → Ok(default)` 提前返回；malformed seam 现在与 body 走同一个 `process`，报同一条错误。e2e `malformed_record_under_site_filter_is_fatal` |
| D9-5 `in_flight_peak_bytes` 口径 | P2 | **已修（选项 A）** | rustdoc 改为"估算足迹 = 压缩 × `MEM_FACTOR`"；顺序路径填真实峰值（原硬编码 0）。测试 `in_flight_peak_is_the_estimated_footprint`（两条路径 peak>0 且为 `MEM_FACTOR` 的整数倍）。**未更名**（字段名保留，只改语义文档——更名会波及 bench/profile 消费方，超出本轮范围） |
| D9-6 reducer 上的 seam 时间漏计 | P2 | **已修** | 并行 `deliver!` 给 `run_seam` 计时并入 `job_time`、seam 的 `accept` 并入 `reduce_time`；`metrics.rs` 的 `job_time`/`reduce_time` rustdoc 改为跨两条路径同语义。测试 `seam_time_is_counted_on_the_parallel_path`（修复前并行 ~49 µs，修复后 >50 ms） |
| D9-7 spec §7.7 过时 | P2 | **已修** | `docs/spec/architecture.md` §7.7 回改为背压表述（依据 v2/03 §4.2）；§7.5 追加规约 W2。守卫 "spec §7.7 describes backpressure, not an error" |
| D9-8 release-only 门盲区 | P1 | **已修** | CI 增加 debug `cargo test --workspace`（`ci.yml`，在 release 之前）；红测试已随 D9-3 转绿；`v8/STATUS.md:22-23` 的归因更正登记在本轮（v8 只增不改） |

## 本轮未做（有意）

| 项 | 原因 |
|---|---|
| Q3-d `debug_assert!` 必须附 `invariant:` 理由的 CI 守卫 | 【荐】项。现有 20 处 `debug_assert` 中约 10 处无理由，加守卫需顺带改 10 个无关文件，属范围外；D9-8 的 CI debug 运行已覆盖"断言类失败不可见"这一根因。**shortcut：待有人批量补齐理由时再上守卫** |
| D9-5 的字段更名 `in_flight_peak_footprint_bytes` | 更名会波及 `--profile`/bench 消费方；本轮只改语义文档，字段名保留 |
| `run_sequential` 合并进 `run_parallel` | v2/03 §4.3 的旧建议，v6_y5 S4 已选择保留独立顺序路径 + 对拍门；不在缺陷修复轮范围 |

## 复审基线（2026-10-13，**修复前**）

复审 commit：`5ef4454`（doc_v8 全部落地）。以下为修复前实测，记录在此以说明每条门修的是什么。

| 实验 | 夹具 | 修复前 |
|---|---|---|
| reducer sink 错误（窗口外） | 2788 chunk × `budget(4)` | **挂死**：`timeout 45/90/120` 无返回；`eu-stack` 5 线程全 `futex_wait`（reducer 停在 `thread::scope` join，reader 停在 `crossbeam select`，2 worker 同） |
| reducer sink 错误（窗口内，对照） | 3 chunk × `budget(4)` | 971 µs 正常返回 → 证明挂死只由"reader 是否先到 EOF"决定 |
| `run_seam` 错误 | 2788 chunk × `budget(4)` | **挂死**（exit 124） |
| reducer panic | 同上 | **挂死**（panic 消息已打印，进程永不返回） |
| 末条未换行记录 | 50 条去掉末行 `\n` | 交付 **49/50**（顺序与并行一致） |
| `in_flight_peak_bytes` | 单 chunk，max_compressed=149 | `peak = 1192` = 149 × **8.000**；顺序路径 `peak = 0` |
| seam 时间 | `run_seam` sleep 2 ms × 134 seam | threads=1 `job_time` 273.12 ms；threads=4 `job_time` **49.2 µs** |
| stitcher 红测试 | `cargo test -p reed-engine --lib` | debug 红（`pipeline.rs:157`）/ release 绿 |

## 测试与守卫基线

| 命令 | 结果 |
|---|---|
| `cargo test -p reed-engine`（debug） | **39 通过 0 失败**（lib 16 + 集成 23，含 9 个新增 v9 测试） |
| `cargo test --workspace`（debug，CI 新增门） | **52 个 test target 全绿** |
| `cargo test --release --workspace` | 全绿 |
| `cargo test -p reed --test e2e_stats` | 12 通过（含新增 `malformed_record_under_site_filter_is_fatal`） |
| `cargo clippy --release --workspace --all-targets -- -D warnings` | 0 错误 |
| `bash ci/check-doc-drift.sh` | 38 条全绿（含 4 条新增 v9 守卫） |
