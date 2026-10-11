# Y5 TODO（执行清单）

> 规格见 [README.md](README.md)，难点参考实现见 [CODE.md](CODE.md)。
> 每项验证结果已写进 `docs/rounds/v6/STATUS.md` 的 Y5 段。全部完成（2026-10-11）。

## 批 1 · 引擎代码（一个提交 `Y5 — abort/permit/seam/sequential/metrics`）

- [x] **S1a** `abort.rs`：`AbortSignal`/`AbortHandle` 共享 `Arc<AtomicBool>`；
      `is_aborted` 读 flag；`trigger` 先 store 再灌通道再 `notify_all` 注册的 permits。
- [x] **S1b** `abort.rs`：`abort_channel_for(threads)`、`abort_for(&ThreadBudget)`；
      `abort_channel()` 保留（测试用）。
- [x] **S1c** `budget.rs`：`AbortHandle::register(&Arc<BytePermits>)`；
      `BytePermits::acquire` 改无超时 `cv.wait`（循环内查 flag）。
- [x] **S1d** 调用点换 `abort_for`：`reed/src/filter/mod.rs:589,637`、
      `reed/src/index.rs:511`、`reed/src/stats/job.rs:1887`、
      `reed/src/stats/writers/reports.rs:798`、`reed/src/filter/pass0/matrix.rs:889`。
- [x] **S2** `pipeline.rs`：`PermitGuard`；`InFlightChunk`/`WorkerMsg` 持 guard；
      reader 删三处手写 release；worker 一次解构；reducer 交付后 drop；
      `Shared` 加 `in_flight_final_bytes()`。
- [x] **S3** `pipeline.rs`：`MAX_RECORD_BYTES` pub + `NewlineStitcher`；
      reducer 改用它，删 `stitch_seam`/`finish_fragment`；
      `bcf_input.rs` 改引用共享常量。
- [x] **S4** `pipeline.rs`：`run_sequential`；`run()` 按 `is_sequential()` 分派；
      `budget.rs` 与 `architecture.md:804` 的 F6-6 现状段改写。
- [x] **S5a** `metrics.rs`：六个新字段 + `merge`。
- [x] **S5b** `pipeline.rs`：reader/reducer 等待、inflate/job 两段、最终在途字节的原子。
- [x] **S5c** `block_writer.rs`：`ChunkOutput.deflate_nanos` + `BlockWriter` 累计。
- [x] **S5d** `ordered_output.rs`：累计 + `deflate_time()`；
      `OrderedSink::deflate_time()` 默认 0；`OutputSink` 转发；引擎 finish 前读取。
- [x] **S5e** `funnel.rs`：`TimingReport` 五字段 + `timing_text` 五行；
      `filter/mod.rs:666` 构造补齐。
- [x] **构建** `cargo build --release --workspace --all-targets` 0 警告；`cargo fmt`。

## 批 2 · 测试与并行验证

- [x] **T0** 新测试（`crates/reed-engine/tests/pipeline.rs` 与 `pipeline.rs` 内联）：
      - [x] `sequential_and_parallel_outputs_agree`
      - [x] `stitcher_reassembles_records_across_ten_chunks`
      - [x] `stitcher_rejects_a_record_past_the_ceiling`、`stitcher_reassembles_and_caps_fragments`
      - [x] `abort_at_eight_times_the_cores_returns_without_leaking_permits`
      - [x] `worker_panic_returns_every_permit`、`permit_guard_returns_its_bytes_on_drop`
      - [x] `stage_metrics_split_inflate_and_job_time`
      - [x] `consuming_the_token_does_not_clear_the_flag`、`trigger_wakes_registered_permit_waiters`
- [x] **T1–T5 并行**：
      - [x] 子集数字门 + 正文 md5：263,420 / 381,672，`fef4cb07…` / `371d3a0d…`
      - [x] 子集跨格式（`-Ob` 读回、BCF 入）== `-Oz` 正文；`--mind` t1/t56 == `fcc6494c…`
      - [x] 全量 reports 8 件 md5 == Y0–Y4 记录
      - [x] `gold/popgen/run_plink2_ld.sh` PASS（r² + prune）
      - [x] 全量 gold `run_reed.sh` 正文 md5 `7aedc59a…` + `check_parity.py` PASS
- [x] **T6** `cargo test --release --workspace`：398 通过 / 0 失败。
- [x] **T7** 串行：clippy（CI 命令）、rustdoc、`ci/check-doc-drift.sh` 29/29 ok。
- [x] **T8** 子集 `-Oz` 计时 ×3（1.118/1.771/1.121 s）+ `--timing` 新字段冒烟；
      行写入 `bench/suites/e2e/results.tsv`。

## 批 3 · 收口

- [x] **D1** `docs/rounds/v6/STATUS.md` Y5 段（含偏差表）。
- [x] **D2** `CHANGELOG.md` Unreleased 补条目。
- [x] **D3** `docs/README.md` 索引补 `docs/spec/v6_y5/`。
- [ ] **D4** `git commit`。
