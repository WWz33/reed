# 06 修复路线（Q0–Q3）

> 原则：**先让错误能返回（消除挂死），再让错误判定一致（消除静默丢数据），
> 最后让文档/指标与实现一致（消除漂移），并补上让前面这些逃过验收的门**。
> 本轮的修复对象是 01-review.md §3 的缺陷清单；
> 不重开 v6 的 Y 系列、v7 的 Z 系列、v8 的 P 系列结构项。
> 规范依据见 02-architecture.md 的对应小节。每一步都必须满足 §1 通用门。

---

## 1. 通用验收门（每一步都必须满足）【规】

1. **双 profile 全绿**：`cargo test --workspace`（debug）**与**
   `cargo test --release --workspace` 都通过；
   `cargo clippy --release --workspace --all-targets` 不新增警告。
   debug 运行是本轮新增的门（02 §8，修 D9-8）——`debug_assert` 类不变量
   在 release 下被编译掉，只跑 release 有盲区。
2. **每个 P0/P1 缺陷先写出复现它的失败测试，再修复**（红 → 绿），
   测试名描述缺陷本身。
3. **挂死类测试必须自带时间上界断言**（`assert!(elapsed < Duration::from_secs(5))`），
   不得依赖 CI 超时或 `#[should_panic]`——否则测试红的时候会把 CI 挂住
   （这正是 D9-1 的表现形式，02 §8）。
4. **规模不变量**：每个错误路径测试都要在两种夹具规模下跑——
   chunk 数 **< 票据窗口**（reader 先到 EOF）与 **≫ 票据窗口**
   （reader 阻塞在票上）。两者的结果类型与时间上界必须相同。
   这是 01-review §1.3 对照实验的固化；既有测试全部只覆盖前者。
5. **确定性**：threads ∈ {1, 2, 3, 7, 56} 输出 sha256 相同（filter 与 stats 都测）。
6. **真值门（本轮新增）**：对以 `\n` 结束的输入，交付的记录数必须等于
   串行扫描数出的记录数（02 §4.4）。这条补 §7.6"只比较路径之间"的盲区。
7. **性能**：不回退超过 10%，并在 `bench/suites/e2e/results.tsv` 追加一行
   （含 threads=8 与 threads=1 各一次，因为 §6.2 给 reducer 加了计时）。
8. **同一 PR 内删除被替代的旧实现**，不留"以后再删"的并行代码
   （v2 06 §1 纪律 7 沿用）。
9. 本目录相应条目标记为已完成并注明 commit；`docs/rounds/v9/STATUS.md` 记录状态。

---

## 2. 步骤

### Q0 让错误能返回（D9-1 / D9-1b / D9-1c / D9-2）— **P0，最高优先**

> Q0 是本存在的理由：在它落地前，"输出盘写满"或"索引遇到未排序输入"
> 会让进程永久挂死而不是退出。其余各步都排在它后面。

| 项 | 内容 | 门 |
|---|---|---|
| Q0-a | 实现 reducer 的唯一错误出口 `bail`（02 §2.2）：`first.record(err)` → `drop(result_rx)` → `return Err(Aborted)`。把 `deliver!` 的 5 处 `?`（`pipeline.rs:534`/`:535`/`:536`/`:541`/`:545`）改为经它；`sink.finish()?`（`:599`）同样经它，但**它不是挂死点**（02 §2.3 已证明到达该处时 reader 与 worker 必然已退出），改它只为满足规约 W2 的一致性。`:551`/`:561`/`:587` 三处已有的正确形状一并收编 | 复现测试先红：`reducer_sink_error_returns_instead_of_hanging`（2788 chunk × `budget(4)`，sink 在第 5 个 chunk 报错）。修复前**挂死**（实测 120 s 不返回）；修复后 5 s 内返回，且 `Err` 携带 sink 的原始错误而非 `Aborted`（验证 `run` 的 `take_shared` 优先报告因，`pipeline.rs:231-233`）。通用门 4 的两种规模都要跑 |
| Q0-b | `run_seam` 的错误路径（D9-1b） | `reducer_seam_error_returns_instead_of_hanging`：同上规模，`run_seam` 恒返回 `Err`。修复前挂死（实测 exit 124） |
| Q0-c | reducer 循环包 `catch_unwind(AssertUnwindSafe(...))`（02 §2.4），unwind 时 record 后走 `bail`；错误变体按 02 §2.4 的【决】定（推荐新增 `ReducerPanicked`） | `reducer_panic_returns_instead_of_hanging`：sink 在第 5 个 chunk `panic!`。修复前 panic 消息已打印但进程永不返回（实测）；修复后 5 s 内返回 `ReducerPanicked`，消息逐字保留，且 `abort.signal.is_aborted()` |
| Q0-d | 两处 `send(...) -> _ => true` 改为 `-> res => res.is_ok()`（02 §3）：`pipeline.rs:498-507`、`:407-410`；同步修正 `:413-416`、`:508-510` 的注释 | `worker_send_on_disconnected_result_channel_breaks`：断言 worker 处理的 chunk 数有上界（不是跑完全部 2788 个）。`reader_send_on_disconnected_chunk_channel_breaks`：同型。许可不泄漏由既有 `permit_guard_returns_its_bytes_on_drop` 与 `in_flight_final_bytes == 0` 继续覆盖 |
| Q0-e | 改写 `pipeline.rs:521-524` 的错误注释（01-review §2.1 指出三处与代码不符） | 注释描述与 `bail` 的实际行为一致；grep 守卫：不得再出现 "the dropped senders/receivers unblock every other stage" |

> **Q0 完成标志**：reducer 侧的 job 错误、sink 错误、seam 超限、panic
> 四种形态，在 chunk 数 ≫ 票据窗口时都能在 5 s 内返回**原始错误**，
> 且 `in_flight_final_bytes == 0`（许可无泄漏）、`abort.signal.is_aborted()`。

### Q1 让错误判定一致（D9-3 / D9-4）— **P0 / P1**

| 项 | 内容 | 门 |
|---|---|---|
| Q1-a | `NewlineStitcher::finish(self) -> Option<Vec<u8>>`（02 §4.2） | 单元测试：以行边界结束 → `None`；有残留 → `Some(fragment)` 且内容正确 |
| Q1-b | 按 02 §4.3 的【决】落地截断判定（推荐选项 A：报错，消息与 BCF 路径一致 `"input ends inside a record (truncated file)"`）。**两条路径都要改**：并行 `:592-600`、顺序 `:705-710`，都在 drain 后、`sink.finish()` 前检查 | 复现测试先红：`unterminated_final_record_is_rejected`（50 条记录去掉末行 `\n`）。修复前**交付 49 条且返回 Ok**（实测）；修复后顺序与并行都返回 `Err`。`terminated_final_record_is_accepted`：带 `\n` 的正常输入仍交付 50 条（防过度收紧）。e2e：末行无 `\n` 的 BGZF VCF 夹具，`reed filter` / `reed stats` / `reed index` **三者都失败且消息一致** |
| Q1-c | 真值门（02 §4.4、通用门 6） | `delivered_records_equal_input_records`：对以 `\n` 结束的输入，`sink_report` == 串行扫描的记录数。作为常驻验收测试 |
| Q1-d | 删除 `filtered.rs:327-329` 的提前 `Ok(default())`（02 §5 规约 E1）。**保留** `:332-335`、`:341-344` 那两个"过滤器正常丢弃"分支——它们与畸形记录是两回事 | 复现测试先红：`malformed_seam_record_is_fatal`（构造 <9 列的记录**强制跨接缝**）。`malformed_record_fatality_is_chunk_alignment_independent`：同一条行放在体内 vs 跨接缝，错误类型与消息**逐字相同**（都必须是 `"record line has fewer than 9 columns"`）。threads ∈ {1,2,3,7,56} 下结论一致 |

> **Q1 完成标志**：同一份输入，filter / stats / index 三条路径对
> "记录是否完整、是否合法"的判定完全一致；
> 一条畸形记录的致命性不再取决于它是否恰好跨 chunk 接缝。

### Q2 让文档与指标一致（D9-5 / D9-6 / D9-7）— **P2**

| 项 | 内容 | 门 |
|---|---|---|
| Q2-a | 按 02 §6.1 的【决】定 `in_flight_peak_bytes` 口径（推荐选项 A：更名 `in_flight_peak_footprint_bytes` + rustdoc 写明 = 压缩 × `MEM_FACTOR` + 顺序路径填真实峰值而非恒 0） | `in_flight_peak_matches_its_documented_semantics`：peak 等于该 run 的最大在途估算足迹和；顺序路径 peak **非 0**（当前恒 0，`pipeline.rs:715`）。既有 `memory_budget_bounds_in_flight_bytes` 必须继续通过（peak ≤ budget） |
| Q2-b | 按 02 §6.2 给并行路径的 `run_seam`（`:535`）与 seam 的 `sink.accept`（`:536`）加计时，分别并入 `job_ns` 与 `reduce_ns`；`metrics.rs:22-23`、`:34-36` 的 rustdoc 改为跨两条路径同语义 | `seam_time_is_counted_on_both_paths`：`run_seam` sleep 2 ms × N 个 seam，断言并行路径 `job_time ≥ N × 2 ms × 0.8`（修复前实测 49.2 µs，而顺序路径 273.12 ms）。`reduce_time_covers_seam_accept`。既有断言 `worker_busy >= inflate_time + job_time` 两条路径都继续成立（02 §6.3） |
| Q2-c | 回改 `docs/spec/architecture.md:838-841`（§7.7）为背压表述，`:820-823`（§7.5）追加规约 W2（02 §7） | `grep -n "MemBudgetExceeded" docs/spec/architecture.md` 的上下文出现"背压"/"单个 chunk"；`ci/check-doc-drift.sh` 新增守卫：spec §7.7 不得再声称"超出即报错" |

### Q3 防漂移纪律（D9-8）— **P1**

| 项 | 内容 | 门 |
|---|---|---|
| Q3-a | 修好 D9-8 的红测试：把 `pipeline.rs:776` 的 `push_head(None)` 改为 `stitcher.finish()` 的正确断言（02 §4.2）；给 `:157-160` 的 `debug_assert` 补"为什么生产不可达"的理由（02 §8） | `cargo test -p reed-engine --lib stitcher`（**debug**）全绿；release 同样全绿。这是通用门 1 首次真正生效 |
| Q3-b | 通用门 1 的 debug 运行写进 `ci/`（02 §8） | CI 脚本含一次不带 `--release` 的 `cargo test --workspace`；对现有全绿 |
| Q3-c | 更正 `v8/STATUS.md:22-23` 的归因（v8 只增不改，故在本目录登记）：该红测试不是"与本轮无关"，而是 D9-3 的症状 + D9-8 的实例 | 01-review §4 已登记；`docs/README.md` 的 v9 行指向它 |
| Q3-d | 【荐】`ci/check-doc-drift.sh` 增守卫：`crates/*/src` 里的 `debug_assert!` 必须附 `invariant:` 理由（02 §8） | 守卫脚本进 `ci/`；对现有代码全绿 |

---

## 3. 顺序与依赖

```
Q0（a → b → c → d → e）─► Q1（a → b → c，d 可并行）─► Q2（a、b、c 互相独立）─► Q3（a 依赖 Q1-a）
```

- **Q0-a 是全部的前提**：在它落地前，Q1-b 的截断报错会在大文件上**挂死**
  而不是返回错误（因为报错点就在 reducer 的 `deliver!`/`finish` 里）。
  这一点必须写进 PR 描述——Q1-b 不能先于 Q0-a 合并；
- Q0-b/Q0-c 依赖 Q0-a 的 `bail` 出口；Q0-d 独立但【荐】同 PR
  （它让 Q0-a 的唤醒更快、并消除无用功）；Q0-e 依赖 Q0-a 定稿；
- Q1-a 是 Q1-b 的前提；Q1-c 独立（但【荐】与 Q1-b 同 PR，
  因为它就是防止 Q1-b 过度收紧的安全网）；Q1-d 与 Q1-a/b 完全独立，可并行；
- Q2 三项互相独立，可分三个 PR；
- Q3-a 依赖 Q1-a（`finish()` 存在才能改测试断言）；Q3-b/c/d 独立。

**最小可交付子集**：若只做一件事，做 **Q0-a**（消除挂死）。
若做两件，加 **Q1-b**（消除静默丢记录）——但注意二者的顺序不可颠倒。

---

## 4. 本轮不做什么

- 不改票据窗口、许可 drop guard、`PermitGuard` 的析构语义（v6_y5 S2 已验证，
  本轮实测 `in_flight_final_bytes == 0` 仍成立）；
- 不改 `FrameChain` 的 10 ms 轮询（`v6_y5/README.md` §4 已登记的已知偏差，
  中止延迟 ≤10 ms 满足 §7.5 的门）；
- 不改 rayon 三阶段与 per-contig 上限的预算接线（v8 D8-1…D8-6 已落地，
  v8 的线程创建点台账无需修订）；
- 不把 `run_sequential` 合并进 `run_parallel`
  （`v2/03-architecture.md` §4.3 曾建议删除，v6_y5 S4 选择了保留 + 逐字节对拍门；
  重开它超出缺陷修复轮范围）；
- 不动 `IndexJob` 解析 chunk 内 tail 的逻辑（`index.rs:188-194`）——
  它处理的是会被下一个 chunk 的 head 补完的 tail，与流末尾的残留 fragment 无关
  （02 §4.3 要点 2）；
- 不新增统计语义、不改 CLI 参数（v2 06 §4 纪律 4 沿用）。
