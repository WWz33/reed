# doc_v4 实施状态

> 对照 `06-roadmap.md` 的 W0–W6。每项注明 commit 与验证方法（测试名或可复现命令）；偏差单独登记。格式沿用 `docs/rounds/v3/STATUS.md`。

## 已完成

### W0 语义与精确判定（02 首版；**CLI 拼写部分已被 W0b 取代**）

> 本表记录 W0 的实现层成果：精确有理数判定、等位计分子、掩码后求值、头部存活记录作分母、
> 漏斗机制、`to_args()` 往返、`sample_rule` 报告形态。**凡涉及 `--max-missing` 家族 CLI 拼写的行，
> 已由下一条 W0b 取代（实现层全部保留）**；原行保留并标注，以示来路。

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| ~~`CallRate` 类型~~ **[W0b：改名 `MissingRate`，判定式改缺失方向]** | `reed_filter::spec::CallRate`：十进制量化为有理数 `q/10^9`，`keeps(called, total)` 全 u128 整数判定，`complement()`/`render()`（最短十进制）/serde 往返 | 本提交 | 单测 `call_rate::tests::{exact_boundary_decisions, render_round_trips_and_is_shortest, complement_is_exact_decimal, serde_round_trips, out_of_range_is_rejected}`；`eval::tail::tests::call_rate_boundaries_are_exact`（`X=0.8,N=10` 等号成立、`X=0.7` 的 `1−X` 浮点陷阱） |
| ~~位点级参数~~ **[W0b 取代：`--geno X`]** | `--max-missing FLOAT`（取代 `--geno`）、`--max-missing-count INT`（新）；分子改**等位计**（`an ≥ X·2N`、`2N−an ≤ C`）；`TailRule::{MinCallRate, MaxMissingCount}`；漏斗 id `after_mask.max_missing`/`after_mask.max_missing_count` | 本提交 | 子集（`gold/work/subset.Chr01-03.vcf.gz`，`--min-dp 5 --min-gq 20`）：`--max-missing 0.9` → **263,420**、`0.8` → **381,672**（== 旧二进制 `--geno 0.1/0.2` 的条数）；`--max-missing-count 69` → 263,420（与 0.9 边界一致）；`e2e_missing_boundary`（N=10/7、`X∈{0,0.05,0.1,0.5,0.7,0.714,0.72,0.75,0.8,0.857,0.9,1}`）与 vcftools 逐位点一致，无半调用用例 == plink2 `--geno 1−X` |
| ~~样本级参数~~ **[W0b 取代：`--mind X`]** | `--max-missing-sample FLOAT`（取代 `--mind`）只属 `reed filter`；pass0（`MissingJob`/sketch 的 cohort 通道）由"完全缺失基因型计数"改为"已调用等位计数"（半调用 1、掩码 0）；`FilterSpec.min_sample_call_rate`（W0 时名 `FilterSpecV2`，W1 改名） | 本提交 | `e2e_filter::sample_call_rate_counts_post_mask_over_head_survivors`（三策略；头部失败记录不入分母；`0.5` 恰在界上保留、`0.500000001` 剔除；报告 `sample_rule.n_dropped`）；`e2e_missing_boundary::sample_call_rate_matches_plink2_mind` 与 plink2 `--vcf-min-gq 20 --mind 0.1` 保留集合一致（{s0,s2}，s2 恰在 0.9 界上） |
| ~~删除旧参数~~ **[W0b 反转：`--geno/--mind` 恢复为正式参数，`--max-missing*` 改为报错]** | `--geno X`/`--mind X` 在解析期报错（退出码 2）并给出精确十进制换算（`--geno 0.1`→`--max-missing 0.9`、`--geno 0.8`→`--max-missing 0.2`）；**不提供**静默换算的隐藏别名 | 本提交 | CLI 单测 `removed_flags_error_with_the_conversion`；手工：`filter --geno 0.1`、`stats --mind 0.2` 均 rc=2 且提示含目标参数 |
| `reed stats` 的样本级参数 **[W0b：改为拒绝 `--mind`]** | 隐藏参数 `--max-missing-sample`/`--mind` 一律拒绝并提示先用 `reed filter`（02 §3.2） | 本提交 | 手工：`reed stats x.vcf --max-missing-sample 0.9` rc=2，提示"run `reed filter --max-missing-sample X` first" |
| ~~spec 字段名 / `to_args` / provenance~~ **[W0b 回退：`max_missing`/`mind`，删 `max_missing_alleles`]** | `TailRules.max_missing`→`min_call_rate`、新增 `max_missing_alleles`、`FilterSpecV2.mind`→`min_sample_call_rate`；`to_args()` 渲染 `--max-missing/--max-missing-count/--max-missing-sample`；`##reed_filterCommand` 随之 | 本提交 | `translate::filter::tests::{spec_to_args_round_trips_through_the_flags, sample_call_rate_round_trips, static_subset_and_sample_rate_round_trip}`、`spec::render::tests::{renders_a_combined_sample_selection, renders_the_documented_flags}`；gold 输出的 `##reed_filterCommand` 含 `--max-missing 0.9` |
| ~~gold 脚本换算~~ **[W0b：删 awk，直接 `--geno "$SITE_MISS_MAX"`]** | `gold/scripts/run_reed.sh` 用与 `run_vcftools.sh` 相同的 awk 精确换算 `1−SITE_MISS_MAX`（`gold/` 未纳入版本控制，改动在工作区） | 本提交 | 全量 `run_reed.sh` 正文 md5 **7aedc59ac798729a2674e6e7ef19adfc**（与 V0–V6 记录相同）；`check_parity.py out/reed.vcf.gz --gold out/full/gold.vcf.gz` **PASS**（1,768,759 位点，jaccard 1.0，基因型不一致 0） |

W0 门（02 §7）六条全部通过；`cargo test --workspace` 全绿，`cargo build --release --workspace --all-targets` 0 警告，`cargo clippy --release --workspace --all-targets` 0 error。

### W0b CLI 拼写回退（02 回退版 §5/§9；用户第二次裁决）

| 项 | 内容 | 验证 |
|---|---|---|
| 恢复的正式参数 | 可见的 `--geno X`（位点缺失率 ≤ X）与 `--mind X`（样本缺失率 ≤ X，仅 filter），值即用户写的缺失比例，CLI 路径无 `1−x` 算术 | `translate::filter::tests::{spec_to_args_round_trips_through_the_flags, mind_round_trips, static_subset_and_mind_round_trip}`、`spec::render::tests::{renders_a_combined_sample_selection, renders_the_documented_flags}`；`cli_contract::help_lists_the_renamed_flags_only` |
| 撤回的参数 | `--max-missing X`/`--max-missing-sample X` 解析期报错并给**精确十进制**换算（`--max-missing 0.9` → `use --geno 0.1`、`--max-missing-sample 0.9` → `--mind 0.1`）；`--max-missing-count C` 报错并说明本族只有两个参数、计数上限写作 `C/2N`；均退出码 2，无静默别名 | `translate::filter::tests::withdrawn_flags_error_with_the_conversion`（四例）；手工：`filter --max-missing 0.9`、`--max-missing-count 3` 的报错文本 |
| 阈值类型与判定式 | `spec::MissingRate`（`tolerates(missing, total)`：`missing·10^9 ≤ q·total`，u128 精确），`complement()`/`render()` 保留供报错文本；位点 `2N−an ≤ X·2N`、样本 `2L−called ≤ X·2L`（分子仍是**等位计**） | `missing_rate::tests::{exact_boundary_decisions, complement_is_exact_decimal, render_round_trips_and_is_shortest, serde_round_trips}`；`eval::tail::tests::missing_rate_boundaries_are_exact`（`X=0.2,N=10,missing=2` 等号保留、`missing=3` 剔除；`X=0.3,missing=3` 保留） |
| spec 字段 / 漏斗 id / 报告 | `TailRules.max_missing`、`FilterSpec.mind`（删 `max_missing_alleles`、删 `TailRule::MaxMissingCount`）；漏斗 id `after_mask.missing`、`sample.mind`（常量集中在 `funnel.rs`）；`sample_rule.threshold` = 用户写的 X，`rule` = `sample.mind` | `ci/check-doc-drift.sh`（守卫按 §7.5 翻转：`after_mask.missing` 在 `funnel.rs`、`--max-missing` 不出现在代码文本里）；`e2e_filter::mind_counts_post_mask_over_head_survivors`（`rule == "sample.mind"`、`threshold == "0.5"`、`n_dropped == 1`） |
| 数字与 gold | 子集 `--min-dp 5 --min-gq 20 --geno 0.1` → **263,420**、`--geno 0.2` → **381,672**（== §2.3 三家基线）；gold 全量正文 md5 `7aedc59ac798729a2674e6e7ef19adfc` 不变，`check_parity` PASS，`##reed_filterCommand` 含 `--geno 0.1` | `e2e_missing_boundary`（N=10/7、含/不含半调用；`X ∈ {0,0.1,0.2,0.3,0.5,0.9,1}` 与 0.1428572/0.2857143/0.4285715 等重复小数边界）：reed `--geno X` == plink2 `--geno X` == vcftools `--max-missing (1−X)`（换算用 `complement().render()` 精确十进制）；样本级 `--min-gq 20 --mind 0.1` == plink2 `--vcf-min-gq 20 --mind 0.1`（保留 {s0,s2}，s2 恰在界上） |
| gold 脚本 | `gold/scripts/run_reed.sh` 删除 W0 的 awk 换算，恢复 `--geno "$SITE_MISS_MAX"`（与 `run_plink2.sh` 同拼写） | 全量 `run_reed.sh` 正文 md5 与 `check_parity` 见上 |

### W1 CLI / 文档一致性（03 §3.3、§3.6；01 F4-2/3/4/9/11/12/13、A4-7）

| 项 | 内容 | 验证 |
|---|---|---|
| `--verbosity` 统一（F4-3） | 新增 `cli::verbosity::Verbosity`：0 silent / 1 summary（默认）/ 2 timing / 3 debug，>3 解析期报错（退出码 2）；filter 的 stderr 输出从驱动移到 CLI（`reed/src/filter/report.rs` 的 `eprint_report` 删除、`FilterOptions.quiet` 删除，`FilterReport.timing_text()` 供 level 2 与 `--timing FILE` 共用）；filter level 1 = 一行摘要 + 漏斗表 + `sample.mind` 行；stats/index 同一映射（stats 的 `let _ = verbosity` 删除，改按级别打印，level 2 打印 open/run 阶段计时）；level 3 安装 tracing 订阅器（`RUST_LOG` 优先，默认 `debug`），filter/stats/index 记 debug 事件 | `crates/reed-cli/tests/cli_contract.rs::verbosity_levels_are_uniform`；`cli::verbosity::tests::levels_map_and_above_three_is_rejected` |
| `reed index` 拒绝 stdin（F4-4） | `input == "-"` 报 "reed index needs a seekable file: stdin is not supported（索引存文件尾偏移且需重读）"；`--help` 的 `<INPUT>` 改为 "file only" | `cli_contract::index_rejects_stdin_with_a_clear_error`（断言无 "No such file"） |
| `<INPUT>` 帮助统一（F4-13） | `cli.rs` 的 `INPUT_HELP` 常量，filter 与 stats 共用（`#[arg(help = INPUT_HELP)]`） | `cli_contract::input_help_is_shared` |
| 去 `-g` 短名（F4-9） | `stats --pop` 去掉 `-g`（保留可见别名 `--populations`、隐藏 `--groups` 一个版本） | `cli_contract::help_lists_the_renamed_flags_only` |
| `FilterSpecV2` → `FilterSpec`（F4-12） | 版本后缀一并清掉：`FilterSpecV2`→`FilterSpec`、`SiteRulesV2`→`SiteRules`；相关注释里的 "v2 spec"/TOML 措辞更新 | `ci/check-doc-drift.sh` 的 grep 守卫；`cargo build --release --workspace --all-targets` |
| `docs/rounds/v3/STATUS.md` C3 更正（F4-2） | "未完成"表 C3 行改为**已做**并附证据（`reed-core/src/tag_table.rs`、`reed_io::header::tag_table()` 接入 `plan::compile`、`reed-filter/Cargo.toml` 无 noodles） | 同上 grep 守卫（`! grep -q noodles crates/reed-filter/Cargo.toml`） |
| `docs/README.md`（A4-7） | 新增：当前规范（按主题）与历史文档一览、每轮**明确取代的章节**表、守卫脚本说明；`docs/rounds/v2/README.md`、`docs/rounds/v3/README.md` 顶部加"后继修订"说明；`docs/spec/architecture.md` §13.17 的漏斗 id 示例更新 | 人工；文件在 `docs/README.md` |
| `results.tsv` 补登（F4-11、04 §4.5） | 表头新增 `load_1m`、`median_of`、`timing_path`；补登 V0–V6 数字 14 行（负载当时未记录 → `unrecorded`，来源 `docs/rounds/v3/STATUS.md`）、2026-10-09 复审基线 7 行（负载 0.3）、W0 子集门 1 行（3 次中位数 1.244 s，负载 4.5–8.5 忙碌）；timing 文件放 `bench/suites/e2e/timing/` | 每行 18 列（`awk -F'\t'` 校验）；子集行含 263,420 与 381,672 |
| `##reed_command` 与环境信息（F4-8） | 新增 `--no-command`：只省略 `##reed_command`（记录的 argv 含二进制路径/`--threads`），其余 provenance 行保留；`--no-version` 仍省略全部 `##reed_*`。文档化"正文确定、header 含环境信息" | `cli_contract::no_command_drops_only_the_argv_line` |
| 报告模式的 `-t/-T`（F4-10） | `ReportRequest` 增加 `targets`/`targets_exclude`，report job 对每条记录先过 `region_allows`；报告模式不再报 "not supported yet"，`--summary` 也走同一次读取（QC 累加器在 report partial 里），`-r` 的粗区间同样被精确谓词过滤 | `cli_contract::report_mode_accepts_targets`（`-t Chr01:1-5` → `.afreq` 5 行、`.stats` 记录数 5）；`e2e_reports`/`e2e_bcf` 全绿 |
| grep 守卫（03 §3.6） | 新增 `ci/check-doc-drift.sh`：C3 无 noodles、无 `reed-qc` crate、无版本后缀类型、`--geno/--mind` 存在、`--max-missing` 只出现在报错文本、index 拒绝 stdin、`Verbosity` 单一实现 | `bash ci/check-doc-drift.sh` 全部 ok |
| 额外修正：输出索引改为 opt-in（非文档列出） | `WriteIndex::parse(None)` 此前返回 `Csi`，即不给 `-W` 也写索引；`-Oz -o /dev/null` 因此报 `/dev/null.csi: Permission denied`（旧二进制同样复现，属既有缺陷）。改为 bcftools 语义：`-W`/`-W=csi` 才写 | `write::index::tests::flavor_parsing_and_suffixes`；手工 `filter subset --geno 0.1 -Oz6 -o /dev/null` 通过 |
| 额外修正：错误渲染（非文档列出） | `main` 改为 `eprintln!("reed: {err}")`（Display）；此前 `String` 错误经 `Result` 的 Debug 打印带引号（`Error: "…"`） | `cli_contract` 的 index 用例；手工 `reed index -` 输出 `reed: reed index needs a seekable file: …` |

### W2 三处「多付的工作」与计时（O-1/O-2/O-5；01 F4-6/F4-7/F4-14）

| 项 | 内容 | 验证 |
|---|---|---|
| O-1 报告一次读取（F4-6） | `ReportRequest.summary`：QC summary 的累加器（`QcAccumulator`）并入 report job 的 partial——每 chunk 取走上交、sink 按序合并，`reed stats --freq…--summary` 单遍完成；stdin 双族不再落盘（`StdinSpill` 分支删除） | 全量 8 个输出文件与旧两遍实现的 **md5 全同**（afreq/hardy/vmiss/ldepth.mean/idepth/het/smiss/stats）；同一时间窗 A/B：旧 5.25–5.42 s / 231.6 s user → 新 3.62–3.67 s / 178.0–178.7 s（wall −31%，CPU −23%）；子集 0.88 → 0.64 s |
| O-2 传 codes（F4-7） | `FilteredJob::with_records` + `ScannedRecord{line, codes}` + `StatsJob::run_records`：掩码后的 `GtCode` 直接交给内层，跳过 `edit_vcf_line` 的整行重写与内层的 LineRef/SiteView/GT 重扫；文本路径（`FilteredJob::new`）保留给需要字节的消费者 | `stats_site_filters_match_the_filter_then_stats_composition`（CLI，`stats --X` ≡ `filter --X` + `stats` 逐字节）；子集上两条路径的 `_pi.txt` md5 同为 `438b1a1e…`；全量 `--min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 --pi` 3.00 / 3.04 / 3.00 s（负载 1.04 / 5.37 / 5.37）→ 门 ≤3.0 s 达标 |
| O-5 计时扩展（A4-5） | `StageMetrics` 新增 `n_chunks_read`（此前恒 0）、`reader_time`（`next_chunk` 累计）、`reduce_time`（`accept` 累计）、`chain_wait`（`ChunkJob::wait_time`；BCF 链在 `FrameChain::wait` 内累计）与 `merge`（多遍求和、在途峰值取最大）；filter 的 `--timing`/`--verbosity 2` 与 stats 的 `--verbosity 2` 都打印（`chunks_read`/`reader_ms`/`reduce_ms`/`chain_wait_ms`/`worker_busy_ms`） | 子集 VCF：reader 535 ms / reduce 64 ms / chain_wait 0；子集 BCF：reader 613 ms / chain_wait 72 ms；全量五统计：VCF reader 2.16–2.20 s、BCF 2.28–2.43 s（reader 占 wall 的 87%，即 W3 的靶子） |
| F4-14 看门狗 | `STALL_LIMIT` 由「单 cell 等待 120 s」改为「自上次任何 cell 发布以来的全局静默 120 s」，错误文本带等待的 chunk id、等待时长与最后发布的 chunk id；`FrameChain::wait_time()` 供 O-5 读取 | `frame_chain::tests::stall_watchdog_reports_the_global_silence`（缩短限时断言错误文本与发布后恢复）；全部 BCF e2e 通过 |

### W3 reader 吞吐（O-3；03 §3.4）

| 步 | 内容 | 验证 |
|---|---|---|
| 1 fadvise + 8 MiB 缓冲 | **未做（被第 3 步取代）**：mmap 路径没有 read 系统调用，缓冲层不存在；fadvise 需要新依赖且无独立可测收益 | 见第 3 步的门 |
| 2 整 chunk 直读 | `BgzfChunkReader` 的 `append_block` 改为 `fill_to_target`/`read_into`：按 slab 直接读进 chunk 缓冲（不再经 `BufReader` 的二次拷贝），空成员就地 `drain`，`carry` 承接跨 chunk 的块尾；`input.rs` 的文件路径不再包 `BufReader` | 全量五统计：reader VCF 2.16–2.20 → 2.04–2.15 s、BCF 2.28–2.43 → 2.21–2.27 s（约 −3%：真正的开销是每 chunk 新分配 1 MiB 缓冲的缺页，不是缓冲拷贝） |
| 3 mmap 零拷贝 | 新增 `MmapChunkSource`：文件输入映射为 `Arc<Mmap>`，chunk 是映射里的一个区间（新类型 `ChunkPayload::{Owned, Mapped}`，`Deref<Target=[u8]>` 保持所有调用点不变）；切分规则与流式 reader 相同（整块直到 target），空成员留在区间内（`(len, 0)`，索引偏移保持连续）；非普通文件/映射失败回退流式 | 全量五统计：reader VCF 2.16–2.20 → **0.35–0.36 s**、BCF 2.28–2.43 → **0.35–0.37 s**（−84%）；wall VCF 2.44 → **1.71–1.75 s**、BCF 2.57 → **1.47 s**（**BCF ≤ VCF** ✓，W3 的门） |
| 一致性 | 输出与 reader 无关 | 10 个 stats 输出文件与 W2 时逐字节相同；gold 全量正文 md5 `7aedc59ac798729a2674e6e7ef19adfc` 不变 + `check_parity` PASS；子集 CSI sha256 与流式 reader 完全相同（`7048abbd…`）；`cargo test --workspace` 362 通过 |
| 安全性（偏差） | `memmap2::Mmap::map` 是 `unsafe`，工作区 `unsafe_code = deny`；在 `MmapChunkSource::open` 上加了 `#[allow(unsafe_code)]` 并写明契约：只读映射、生命周期覆盖全部 chunk、输入不被写入、非普通文件回退流式。**更正（doc_v5 03 §2）**：工作区共有**两处**允许——本处与 `reed-core/src/count/kernel.rs` 的 AVX-512 计数内核（模块级，自 doc_v1 §14.4 起） | 评审点：两处有意的策略例外，见"未完成 / 偏差"表 |

### W4 BCF 原生过滤输入（O-4；01 F4-5）

| 项 | 内容 | 验证 |
|---|---|---|
| 原生判定输入 | `reed-io`：`render_site_line`（固定列 + 完整 INFO + 占位 FORMAT + 每样本 `.`，保住 `N_SAMPLES`）与 `decode_format_columns`（计划内 Number=1 FORMAT 列的类型化解码；缺失用与文本扫描器相同的哨兵 `INT_MISSING`/`NaN`，负值同 `parse_int` 视为缺失）；GT 复用既有 64K 表 `decode_record_gt` | 新增 `crates/reed-io/tests/bcf_native_input.rs`：fixture 上站点列/INFO/`N_SAMPLES`/DP 列/GT codes 与「整行渲染 + 文本扫描」逐个相同 |
| 管线两段式接口 | `RecordPipeline::process_scanned` + `NativeScan`：跳过文本扫描器，用原生 codes/列走同一套 mask/count/tail/body；`run_stages` 由 `process` 与 `process_scanned` 共用（不是两份实现） | 策略等价与 gold 门；`cargo test --release --workspace` 364 通过 |
| 输出侧仍渲染 | 单遍 `FilterJob::run_bcf_records`：只有存活记录才 `render_record` 全量渲染 → `edit_vcf_line` → 写块（与文本路径共用 edit/索引）；pass0 `MissingJob::run_bcf_records` 同构（`--mind` 的 BCF 输入全程不渲染） | 见下表门 |
| 适用条件与回退 | 仅当 BCF 输入 ∧ 掩码不含 `--set-gts .` FMT 表达式 ∧ 计划内 FORMAT 列均为 Number=1 时走原生，否则回退渲染路径 | `e2e_bcf` 7 条；子集 `--set-gts . -e 'FMT/DP>20'` 正文 md5 `ae4c8d79…` 与 pre-W4 二进制相同（2.4 s） |
| 门（全量 gold，56 线程，负载 2.15） | BCF 输入：wall **8.95 / 9.06 / 9.36 s（中位 9.06）≤ 9.5 s**；user **631.6 → 468.3–469.5 s（−26%）**；worker_busy **639.0 → 472.6–480.9 s（−25%）** | 全量正文 md5 `7aedc59ac798729a2674e6e7ef19adfc` 与 VCF 输入相同；子集 VCF/BCF 正文 md5 相同（`ca00a0de…`；`--mind` 版 `7a2ba1c5…`）；`-Ob` 读回正文 md5 相同（`fef4cb07…`） |

### W5 统计接入（**未做**；前置调查与 fixture 生成器已完成）

| 项 | 状态 | 调查结论（可复现） |
|---|---|---|
| 投影 SFS（05 §2 第 1 项） | 未做，**参照链已打通** | `easySFS` 不在 PyPI，已从 GitHub 取得（`git clone https://github.com/isaacovercast/easySFS`）；其依赖 `dadi` 已装（`pip install --user dadi` → 2.4.4，`numpy`/`scipy` 已在本机）。在 `gold/popgen/make_fixtures.py` 的 `sfs.vcf`（4 样本 1 群体，ALT 计数 1,1,2,2,3,3,4,4,5,6,7）上 `easySFS.py -i sfs.vcf -p sfs.pops --proj 6 -o out --unfolded --dtype int -a -f` 的输出已实测固定：dadi 文件三行 = `7 unfolded` / `1 3 3 2 1 1 0` / mask `1 0 0 0 0 0 1`；fsc2 文件 = `1 observations. No. of demes…` + `1\t6` + 数据行。**待做的改动**：(a) `ProjectedSfsAcc::push_weighted` 只累加 k∈1..n，丢弃 k=0/k=n 的非分离质量，而 easySFS 的数组是 k=0..n，需先保留该质量（`bins` 加两个槽位）；(b) SFS 不是窗口标量，套不进现有 `StatKind`/窗口行模型，需新增 genome-wide 输出路径（`StatsRequest.sfs` + sink 收尾写出）；(c) dadi/fsc2 两个写出格式 + easySFS 的 mask 策略（variant-only 输入时 mask[0]=mask[n]=1；给 `--total-length` 时把单态质量加到 k=0）；(d) `gold/popgen/run_easysfs.sh` + compare 脚本 |
| 联合 SFS / 多群体 WC84 / LD / f 统计 / block jackknife | 未做 | 累加器都在（`#[doc(hidden)]`）；参照物均已就绪：scikit-allel 1.3.13、plink2、vcftools、admixture。每项需要 CLI 出口 + **输出形状设计**（文档未固定：`--fst-pops` 的矩阵形态、LD/pruning 的列、f 统计的 block jackknife 表），属功能设计工作 |
| fixture 生成器（05 §4） | **已交付** | `gold/popgen/make_fixtures.py`：固定、可手算、无随机数；输出 `missing_rate.vcf`（缺失阈值边界 + 半调用）、`sfs.vcf`+`sfs.pops`（投影 SFS）、`multiallelic.vcf`（W6 的逐等位重算）、`fstats.vcf`+`fstats.pops`、`ld.vcf`+`ld.pops`；全部经 `bcftools view -h` 校验，且 `sfs.vcf` 已跑通 easySFS（见上）。`gold/` 未纳入版本控制，改动在工作区 |


## 未完成 / 偏差

| 项 | 状态 | 说明 |
|---|---|---|
| 输出索引的默认值 | 已修正（非文档列出） | 见 W1 表末两行的"额外修正" |
| reports 两族 ≤ 3.2 s（W2 门） | **未达标（已登记）** | 实测 3.58–3.65 s（负载 18.5；CPU 181 s，wall 下限 ≈ 181/54 = 3.35 s）——瓶颈已从 reader 变成 CPU：QC summary 用**自己的扫描器**再解一遍 GT（策略与 reports 相同但带 `tracking_alleles`）。下一步：让 reports job 在 `--summary` 时用带 tracking 的扫描器，并把 QC 的累加部分拆成消费 `ScanOutput` 的入口（预计省 ~25–30 CPU-s，可望 ≤3.2 s）。旧两遍实现为 5.35 s / 232 s CPU，本轮已 −22% CPU |
| `unsafe_code` 策略 | **两处例外（已登记，doc_v5 03 §2 更正）** | `MmapChunkSource::open` 的 `Mmap::map` 与 `reed-core/src/count/kernel.rs` 的 AVX-512 计数内核（模块级）；契约写在调用点，回退路径保证非普通文件不受影响。其余 crate 仍为 `deny` |
| `sample.mind` 的形态 | 偏差（已登记，02 回退版 §5.3 接受） | 现实现没有样本级漏斗条目（报告只有 head/tail/mask 的**记录**计数），故 `sample.mind` 不是 `FilterReport.funnel` 的元素，而是报告新增的 `sample_rule: {rule, threshold, n_dropped}`（常量在 `reed-filter/src/funnel.rs`）；`after_mask.missing` 作为 tail 漏斗 id |
| W5（统计接入：SFS → `--fst-pops` → f 统计 → LD） | **未做**（05 §3 的"报告模式 `-t/-T`"已随 F4-10 完成） | 累加器都已存在且 `#[doc(hidden)]`（`ProjectedSfsAcc`/`JointSfsAcc`/`WcMultiAcc`/`LdDecayAcc`/`F2Acc`…），缺 CLI 出口、请求种类、写出格式与 `gold/popgen` 对拍脚本（easySFS、plink2 LD、scikit-allel f 统计）；每项都要「逐位/逐行对拍通过后才去掉 `#[doc(hidden)]`」，属功能开发 |
| W6（多等位逐等位重算） | **未做** | 现状是**显式拒绝**（默认双等位门；显式允许多等位且需 AC/MAF/HWE 重算时报错），不会静默写错；实现需要 `fill_tags` 按 Number=A 输出逐 ALT 向量（MAF 取第二常见等位）并与 `bcftools +fill-tags` 在多等位 fixture 上逐字节对拍 |

## 最终验证（2026-10-09，W0–W3 全部落地之后）

| 检查 | 结果 |
|---|---|
| `cargo test --release --workspace` | **361 通过 / 0 失败** |
| `cargo build --release --workspace --all-targets` | **0 警告** |
| `cargo clippy --release --workspace --all-targets` | 0 error（pedantic 提示沿用 doc_v3 C7 的豁免政策） |
| 01 发现清单 | **F4-1…F4-14 全部处理**（F4-5 = W4 已完成） |
| `bash ci/check-doc-drift.sh` | 全部 ok（C3 无 noodles、无 `reed-qc`、无版本后缀类型、`--geno/--mind` 存在、`--max-missing` 只在报错文本、index 拒绝 stdin、`Verbosity` 单一实现） |
| gold 全量（W0b + mmap reader 之后重跑） | 正文 md5 `7aedc59ac798729a2674e6e7ef19adfc` 不变；`##reed_filterCommand` 含 `--geno 0.1`；`check_parity.py` 对 `out/full/gold.vcf.gz` **PASS**（1,768,759 位点，jaccard 1.0，基因型不一致 0） |
| 子集数字门 | `--min-dp 5 --min-gq 20 --geno 0.1` → 263,420、`--geno 0.2` → 381,672（== 三家基线） |
| 参考工具对拍 | `e2e_missing_boundary`：plink2 `--geno X` 与 vcftools `--max-missing (1−X)`（精确十进制换算）在 N=10/7、含/不含半调用、含重复小数边界上与 reed 一致；样本级 == plink2 `--vcf-min-gq 20 --mind 0.1` |
| 跨格式/跨 reader 一致 | stats 五统计 VCF 与 BCF 的 5 个输出文件逐字节相同；mmap 与流式 reader 的 10 个 stats 输出、CSI sha256（`7048abbd…`）、gold 正文全部相同 |
| 性能（见 `bench/suites/e2e/results.tsv`，每行带负载与中位数） | stats 五统计：VCF 2.44 → **1.71 s**、BCF 2.57 → **1.47 s**（reader 2.2 → 0.35 s，mmap）；reports 两族 5.35 → **3.58 s**（CPU 232 → 181 s，门 ≤3.2 s 未达，见上表）；边过滤 `--pi` 3.4 → **2.99 s**（门 ≤3.0 s 达标）；filter gold `-Oz6` 8.8 → 8.75 s（不回退） |

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
