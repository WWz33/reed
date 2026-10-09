# doc_v4 实施状态

> 对照 `06-roadmap.md` 的 W0–W6。每项注明 commit 与验证方法（测试名或可复现命令）；偏差单独登记。格式沿用 `doc_v3/STATUS.md`。

## 已完成

### W0 参数语义统一（02）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| `CallRate` 类型（03 §3.2） | `reed_filter::spec::CallRate`：十进制量化为有理数 `q/10^9`，`keeps(called, total)` 全 u128 整数判定，`complement()`/`render()`（最短十进制）/serde 往返 | 本提交 | 单测 `call_rate::tests::{exact_boundary_decisions, render_round_trips_and_is_shortest, complement_is_exact_decimal, serde_round_trips, out_of_range_is_rejected}`；`eval::tail::tests::call_rate_boundaries_are_exact`（`X=0.8,N=10` 等号成立、`X=0.7` 的 `1−X` 浮点陷阱） |
| 位点级参数 | `--max-missing FLOAT`（取代 `--geno`）、`--max-missing-count INT`（新）；分子改**等位计**（`an ≥ X·2N`、`2N−an ≤ C`）；`TailRule::{MinCallRate, MaxMissingCount}`；漏斗 id `after_mask.max_missing`/`after_mask.max_missing_count` | 本提交 | 子集（`gold/work/subset.Chr01-03.vcf.gz`，`--min-dp 5 --min-gq 20`）：`--max-missing 0.9` → **263,420**、`0.8` → **381,672**（== 旧二进制 `--geno 0.1/0.2` 的条数）；`--max-missing-count 69` → 263,420（与 0.9 边界一致）；`e2e_missing_boundary`（N=10/7、`X∈{0,0.05,0.1,0.5,0.7,0.714,0.72,0.75,0.8,0.857,0.9,1}`）与 vcftools 逐位点一致，无半调用用例 == plink2 `--geno 1−X` |
| 样本级参数 | `--max-missing-sample FLOAT`（取代 `--mind`）只属 `reed filter`；pass0（`MissingJob`/sketch 的 cohort 通道）由"完全缺失基因型计数"改为"已调用等位计数"（半调用 1、掩码 0）；`FilterSpec.min_sample_call_rate`（W0 时名 `FilterSpecV2`，W1 改名） | 本提交 | `e2e_filter::sample_call_rate_counts_post_mask_over_head_survivors`（三策略；头部失败记录不入分母；`0.5` 恰在界上保留、`0.500000001` 剔除；报告 `sample_rule.n_dropped`）；`e2e_missing_boundary::sample_call_rate_matches_plink2_mind` 与 plink2 `--vcf-min-gq 20 --mind 0.1` 保留集合一致（{s0,s2}，s2 恰在 0.9 界上） |
| 删除旧参数（02 §5.2） | `--geno X`/`--mind X` 在解析期报错（退出码 2）并给出精确十进制换算（`--geno 0.1`→`--max-missing 0.9`、`--geno 0.8`→`--max-missing 0.2`）；**不提供**静默换算的隐藏别名 | 本提交 | CLI 单测 `removed_flags_error_with_the_conversion`；手工：`filter --geno 0.1`、`stats --mind 0.2` 均 rc=2 且提示含目标参数 |
| `reed stats` 的样本级参数 | 隐藏参数 `--max-missing-sample`/`--mind` 一律拒绝并提示先用 `reed filter`（02 §3.2） | 本提交 | 手工：`reed stats x.vcf --max-missing-sample 0.9` rc=2，提示"run `reed filter --max-missing-sample X` first" |
| spec 字段名 / `to_args` / provenance | `TailRules.max_missing`→`min_call_rate`、新增 `max_missing_alleles`、`FilterSpecV2.mind`→`min_sample_call_rate`；`to_args()` 渲染 `--max-missing/--max-missing-count/--max-missing-sample`；`##reed_filterCommand` 随之 | 本提交 | `translate::filter::tests::{spec_to_args_round_trips_through_the_flags, sample_call_rate_round_trips, static_subset_and_sample_rate_round_trip}`、`spec::render::tests::{renders_a_combined_sample_selection, renders_the_documented_flags}`；gold 输出的 `##reed_filterCommand` 含 `--max-missing 0.9` |
| gold 脚本换算（02 §6） | `gold/scripts/run_reed.sh` 用与 `run_vcftools.sh` 相同的 awk 精确换算 `1−SITE_MISS_MAX`（`gold/` 未纳入版本控制，改动在工作区） | 本提交 | 全量 `run_reed.sh` 正文 md5 **7aedc59ac798729a2674e6e7ef19adfc**（与 V0–V6 记录相同）；`check_parity.py out/reed.vcf.gz --gold out/full/gold.vcf.gz` **PASS**（1,768,759 位点，jaccard 1.0，基因型不一致 0） |

W0 门（02 §7）六条全部通过；`cargo test --workspace` 全绿，`cargo build --release --workspace --all-targets` 0 警告，`cargo clippy --release --workspace --all-targets` 0 error。

### W1 CLI / 文档一致性（03 §3.3、§3.6；01 F4-2/3/4/9/11/12/13、A4-7）

| 项 | 内容 | 验证 |
|---|---|---|
| `--verbosity` 统一（F4-3） | 新增 `cli::verbosity::Verbosity`：0 silent / 1 summary（默认）/ 2 timing / 3 debug，>3 解析期报错（退出码 2）；filter 的 stderr 输出从驱动移到 CLI（`reed/src/filter/report.rs` 的 `eprint_report` 删除、`FilterOptions.quiet` 删除，`FilterReport.timing_text()` 供 level 2 与 `--timing FILE` 共用）；filter level 1 = 一行摘要 + 漏斗表 + `sample.max_missing` 行；stats/index 同一映射（stats 的 `let _ = verbosity` 删除，改按级别打印，level 2 打印 open/run 阶段计时）；level 3 安装 tracing 订阅器（`RUST_LOG` 优先，默认 `debug`），filter/stats/index 记 debug 事件 | `crates/reed-cli/tests/cli_contract.rs::verbosity_levels_are_uniform`；`cli::verbosity::tests::levels_map_and_above_three_is_rejected` |
| `reed index` 拒绝 stdin（F4-4） | `input == "-"` 报 "reed index needs a seekable file: stdin is not supported（索引存文件尾偏移且需重读）"；`--help` 的 `<INPUT>` 改为 "file only" | `cli_contract::index_rejects_stdin_with_a_clear_error`（断言无 "No such file"） |
| `<INPUT>` 帮助统一（F4-13） | `cli.rs` 的 `INPUT_HELP` 常量，filter 与 stats 共用（`#[arg(help = INPUT_HELP)]`） | `cli_contract::input_help_is_shared` |
| 去 `-g` 短名（F4-9） | `stats --pop` 去掉 `-g`（保留可见别名 `--populations`、隐藏 `--groups` 一个版本） | `cli_contract::help_lists_the_renamed_flags_only` |
| `FilterSpecV2` → `FilterSpec`（F4-12） | 版本后缀一并清掉：`FilterSpecV2`→`FilterSpec`、`SiteRulesV2`→`SiteRules`；相关注释里的 "v2 spec"/TOML 措辞更新 | `scripts/check-doc-drift.sh` 的 grep 守卫；`cargo build --release --workspace --all-targets` |
| `doc_v3/STATUS.md` C3 更正（F4-2） | "未完成"表 C3 行改为**已做**并附证据（`reed-core/src/tag_table.rs`、`reed_io::header::tag_table()` 接入 `plan::compile`、`reed-filter/Cargo.toml` 无 noodles） | 同上 grep 守卫（`! grep -q noodles crates/reed-filter/Cargo.toml`） |
| `docs/INDEX.md`（A4-7） | 新增：当前规范（按主题）与历史文档一览、每轮**明确取代的章节**表、守卫脚本说明；`doc_v2/README.md`、`doc_v3/README.md` 顶部加"后继修订"说明；`docs/architecture-spec.md` §13.17 的漏斗 id 示例更新 | 人工；文件在 `docs/INDEX.md` |
| `results.tsv` 补登（F4-11、04 §4.5） | 表头新增 `load_1m`、`median_of`、`timing_path`；补登 V0–V6 数字 14 行（负载当时未记录 → `unrecorded`，来源 `doc_v3/STATUS.md`）、2026-10-09 复审基线 7 行（负载 0.3）、W0 子集门 1 行（3 次中位数 1.244 s，负载 4.5–8.5 忙碌）；timing 文件放 `bench/e2e/timing/` | 每行 18 列（`awk -F'\t'` 校验）；子集行含 263,420 与 381,672 |
| grep 守卫（03 §3.6） | 新增 `scripts/check-doc-drift.sh`：C3 无 noodles、无 `reed-qc` crate、无版本后缀类型、`--max-missing` 家族存在、`--geno` 报错文本存在、index 拒绝 stdin、`Verbosity` 单一实现 | `bash scripts/check-doc-drift.sh` 全部 ok |
| 额外修正：输出索引改为 opt-in（非文档列出） | `WriteIndex::parse(None)` 此前返回 `Csi`，即不给 `-W` 也写索引；`-Oz -o /dev/null` 因此报 `/dev/null.csi: Permission denied`（旧二进制同样复现，属既有缺陷）。改为 bcftools 语义：`-W`/`-W=csi` 才写 | `write::index::tests::flavor_parsing_and_suffixes`；手工 `filter subset --max-missing 0.9 -Oz6 -o /dev/null` 通过 |
| 额外修正：错误渲染（非文档列出） | `main` 改为 `eprintln!("reed: {err}")`（Display）；此前 `String` 错误经 `Result` 的 Debug 打印带引号（`Error: "…"`） | `cli_contract` 的 index 用例；手工 `reed index -` 输出 `reed: reed index needs a seekable file: …` |

### W2 三处「多付的工作」与计时（O-1/O-2/O-5；01 F4-6/F4-7/F4-14）

| 项 | 内容 | 验证 |
|---|---|---|
| O-1 报告一次读取（F4-6） | `ReportRequest.summary`：QC summary 的累加器（`QcAccumulator`）并入 report job 的 partial——每 chunk 取走上交、sink 按序合并，`reed stats --freq…--summary` 单遍完成；stdin 双族不再落盘（`StdinSpill` 分支删除） | 全量 8 个输出文件与旧两遍实现的 **md5 全同**（afreq/hardy/vmiss/ldepth.mean/idepth/het/smiss/stats）；同一时间窗 A/B：旧 5.25–5.42 s / 231.6 s user → 新 3.62–3.67 s / 178.0–178.7 s（wall −31%，CPU −23%）；子集 0.88 → 0.64 s |
| O-2 传 codes（F4-7） | `FilteredJob::with_records` + `ScannedRecord{line, codes}` + `StatsJob::run_records`：掩码后的 `GtCode` 直接交给内层，跳过 `edit_vcf_line` 的整行重写与内层的 LineRef/SiteView/GT 重扫；文本路径（`FilteredJob::new`）保留给需要字节的消费者 | `stats_site_filters_match_the_filter_then_stats_composition`（CLI，`stats --X` ≡ `filter --X` + `stats` 逐字节）；子集上两条路径的 `_pi.txt` md5 同为 `438b1a1e…`；全量 `--min-dp 5 --min-gq 20 --max-missing 0.9 --maf 0.05 --pi` 3.00 / 3.04 / 3.00 s（负载 1.04 / 5.37 / 5.37）→ 门 ≤3.0 s 达标 |
| O-5 计时扩展（A4-5） | `StageMetrics` 新增 `n_chunks_read`（此前恒 0）、`reader_time`（`next_chunk` 累计）、`reduce_time`（`accept` 累计）、`chain_wait`（`ChunkJob::wait_time`；BCF 链在 `FrameChain::wait` 内累计）与 `merge`（多遍求和、在途峰值取最大）；filter 的 `--timing`/`--verbosity 2` 与 stats 的 `--verbosity 2` 都打印（`chunks_read`/`reader_ms`/`reduce_ms`/`chain_wait_ms`/`worker_busy_ms`） | 子集 VCF：reader 535 ms / reduce 64 ms / chain_wait 0；子集 BCF：reader 613 ms / chain_wait 72 ms；全量五统计：VCF reader 2.16–2.20 s、BCF 2.28–2.43 s（reader 占 wall 的 87%，即 W3 的靶子） |
| F4-14 看门狗 | `STALL_LIMIT` 由「单 cell 等待 120 s」改为「自上次任何 cell 发布以来的全局静默 120 s」，错误文本带等待的 chunk id、等待时长与最后发布的 chunk id；`FrameChain::wait_time()` 供 O-5 读取 | `frame_chain::tests::stall_watchdog_reports_the_global_silence`（缩短限时断言错误文本与发布后恢复）；全部 BCF e2e 通过 |

## 未完成 / 偏差

| 项 | 状态 | 说明 |
|---|---|---|
| 输出索引的默认值 | 已修正（非文档列出） | 见 W1 表末两行的"额外修正" |
| `sample.max_missing` 的形态 | 偏差（已登记） | 现实现没有样本级漏斗条目（报告只有 head/tail/mask 的**记录**计数），故 `sample.max_missing` 不是 `FilterReport.funnel` 的元素，而是报告新增的 `sample_rule: {rule, threshold, n_dropped}`（`reed-filter/src/funnel.rs`）；`after_mask.max_missing`、`after_mask.max_missing_count` 按 §5.3 作为 tail 漏斗 id |
| W3–W6 | 待做 | 见 `06-roadmap.md` |

## 基准（2026-10-09 复审基线，56 线程，负载 0.3，热缓存）

| 命令 | wall | user |
|---|---|---|
| filter gold `-Oz6`（VCF 输入） | 8.8 s | 455 s |
| filter gold `-Ob` | 9.9 s | 519 s |
| filter gold，BCF 输入 | 12.1 s | 630 s |
| stats 五统计 VCF / BCF | 2.7 / 3.1 s | 71 / 56 s |
| stats 边过滤 `--pi` | 3.4 s | 170 s |
| reports 两族 | 5.3 s | 232 s |
| index VCF / BCF | 2.4 / 2.2 s | 77 / 46 s |
