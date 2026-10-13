# doc_v8 实施状态

> 本页记录 `01-review.md` §3 缺陷清单与 `06-roadmap.md` P0–P2 的实施状态。
> 每条附证据（测试名 / 守卫 / 变更点）。

| 项 | 严重度 | 状态 | 证据 |
|---|---|---|---|
| D8-1a `ThreadBudget::auxiliary()` + `run_auxiliary()` | P0 | **已修** | `reed-engine/src/budget.rs`：`auxiliary()`（顺序=1，否则 `workers()`）与 `run_auxiliary`（每调用一次的 `ThreadPoolBuilder` 临时池；顺序直接调闭包）；测试 `auxiliary_stages_share_the_worker_budget`、`auxiliary_pool_stays_inside_the_budget`（断言池内 `current_num_threads()==auxiliary()`） |
| D8-1b rayon 三阶段接 budget | P0 | **已修** | `mid_stage.rs`（par_iter 包进 `run_auxiliary`）、`prune.rs::prune_matrices`（into_par_iter 同）、`run.rs::decay_contig`（par_iter 同）——三处均不再落 rayon 全局池；顺序预算在调用线程执行 |
| D8-1c 三处调用点传 budget | P0 | **已修** | `filter/mod.rs`：budget 构造提前到 pass0 之前，`mid_stage(..., &budget)`、`prune_matrices(..., &budget)`；`stats/driver.rs`：`StatsSink { max_ld_jobs: budget.auxiliary(), ld_budget: *budget }`，`finalize_ld` 闭包内传给 `decay_contig` |
| D8-2 stats LD sink in-flight 上限 | P0 | **已修** | `sink.rs::finalize_ld`：`ld_jobs.len() >= max_ld_jobs` 时 join 最旧（新字段 `ld_merged` 累积 cap 期 join 的结果，`finish` 以它为合并基底——contig 顺序不变，输出不变）；测试 `ld_sink_caps_in_flight_decay_jobs`（10 contig × threads {1,8}：30 对全在、两线程数结果相等；threads=1 时 auxiliary=1，每个 contig 都走 cap-join 路径）。该测试先红：抓出第一版 cap 丢弃已 join 结果的真 bug（3/30） |
| D8-3 运行时核占用观测门 | P1 | **已修** | `bench/suites/e2e/results.tsv` 2026-10-13 四行（全量 clean.vcf.gz 6.6 GB，各 3 次中位数）：`--ld-decay` threads=8 wall 25.8 s / ratio 16 / Threads max 9，threads=1 wall 342.8 s / ratio 1.21 / Threads max 2；`--indep-pairwise` threads=8 wall 34.1 s / ratio 6.0 / Threads max 9，threads=1 wall 203.2 s / ratio 1.01 / Threads max 2；`--mind + --indep-pairwise` threads=8 Threads max 8、threads=1 Threads max 1。输出 sha256 跨线程数相同。全部 ≤ N+1 门 |
| D8-6 串行分支的全局池逃逸 | P0 | **已修**（D8-3 实测中发现） | `run_auxiliary` 闭包改 `FnOnce(bool)`：`parallel==false` 时三阶段走 `iter()`/`into_iter()`，不再惰性拉起 rayon 全局池。修复前 `--threads 1 --ld-decay` 实测 58 线程（strace `clone3` 76 次定位：sink 线程里的 `par_iter` 初始化全局池）、wall 76 s；修复后 Threads max 2、ratio 1.21、wall 343 s（真实串行耗时），user CPU 不变。守卫升级为"par_iter 必须在 `if parallel {` 分支内"。测试 `auxiliary_pool_stays_inside_the_budget`、`e2e_ld` 11/11 全绿 |
| D8-4 per-contig 上限统一取数 | P1 | **已修** | `matrix.rs` `max_prune: budget.auxiliary()`（原 `opts.threads.max(1)`）；sink `max_ld_jobs` 同取数点；守卫 `ci/check-doc-drift.sh` "per-contig caps use budget.auxiliary" |
| D8-5a spec §7.3 过时脚注清理 | P2 | **已修** | `docs/spec/architecture.md` §7.3：删除 F6-6"现状"段，写入 doc_v8 落地形态（临时池 + `auxiliary()` + 顺序承诺覆盖辅助阶段） |
| D8-5b CI grep 守卫 | P2 | **已修** | `ci/check-doc-drift.sh` 新增两条："rayon parallelism goes through run_auxiliary"（生产 src 里含 par_iter/into_par_iter 的文件必须同时含 run_auxiliary）、"per-contig caps use budget.auxiliary" |

## 基线说明

- `cargo test -p reed`：100 通过 0 失败；`e2e_ld` 11/11（含新验收测试；D8-6 修复后重跑仍全绿）。
- `cargo test -p reed-stats -p reed-filter -p reed-core -p reed-io`：236 通过 0 失败。
- `cargo test -p reed-engine --lib`：15 通过，1 失败（`pipeline::tests::stitcher_reassembles_and_caps_fragments`，
  "invariant: fragment before the first chunk"）——**预存在**：干净 HEAD e50f6fd 上同样失败，与本轮无关。
- `cargo clippy --workspace --all-targets`：`reed-filter` 2 个 `field_reassign_with_default` 错误——**预存在**（stash 后同样出现）；本轮改动无新增警告（`mid_stage` 加 `#[allow(clippy::too_many_arguments)]`，沿用文件内既有做法）。
- `bash ci/check-doc-drift.sh`：全绿（含 v8 新增两条守卫）。
