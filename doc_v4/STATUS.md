# doc_v4 实施状态

> 对照 `06-roadmap.md` 的 W0–W6。每项注明 commit 与验证方法（测试名或可复现命令）；偏差单独登记。格式沿用 `doc_v3/STATUS.md`。

## 已完成

### W0 参数语义统一（02）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| `CallRate` 类型（03 §3.2） | `reed_filter::spec::CallRate`：十进制量化为有理数 `q/10^9`，`keeps(called, total)` 全 u128 整数判定，`complement()`/`render()`（最短十进制）/serde 往返 | 本提交 | 单测 `call_rate::tests::{exact_boundary_decisions, render_round_trips_and_is_shortest, complement_is_exact_decimal, serde_round_trips, out_of_range_is_rejected}`；`eval::tail::tests::call_rate_boundaries_are_exact`（`X=0.8,N=10` 等号成立、`X=0.7` 的 `1−X` 浮点陷阱） |
| 位点级参数 | `--max-missing FLOAT`（取代 `--geno`）、`--max-missing-count INT`（新）；分子改**等位计**（`an ≥ X·2N`、`2N−an ≤ C`）；`TailRule::{MinCallRate, MaxMissingCount}`；漏斗 id `after_mask.max_missing`/`after_mask.max_missing_count` | 本提交 | 子集（`gold/work/subset.Chr01-03.vcf.gz`，`--min-dp 5 --min-gq 20`）：`--max-missing 0.9` → **263,420**、`0.8` → **381,672**（== 旧二进制 `--geno 0.1/0.2` 的条数）；`--max-missing-count 69` → 263,420（与 0.9 边界一致）；`e2e_missing_boundary`（N=10/7、`X∈{0,0.05,0.1,0.5,0.7,0.714,0.72,0.75,0.8,0.857,0.9,1}`）与 vcftools 逐位点一致，无半调用用例 == plink2 `--geno 1−X` |
| 样本级参数 | `--max-missing-sample FLOAT`（取代 `--mind`）只属 `reed filter`；pass0（`MissingJob`/sketch 的 cohort 通道）由"完全缺失基因型计数"改为"已调用等位计数"（半调用 1、掩码 0）；`FilterSpecV2.min_sample_call_rate` | 本提交 | `e2e_filter::sample_call_rate_counts_post_mask_over_head_survivors`（三策略；头部失败记录不入分母；`0.5` 恰在界上保留、`0.500000001` 剔除；报告 `sample_rule.n_dropped`）；`e2e_missing_boundary::sample_call_rate_matches_plink2_mind` 与 plink2 `--vcf-min-gq 20 --mind 0.1` 保留集合一致（{s0,s2}，s2 恰在 0.9 界上） |
| 删除旧参数（02 §5.2） | `--geno X`/`--mind X` 在解析期报错（退出码 2）并给出精确十进制换算（`--geno 0.1`→`--max-missing 0.9`、`--geno 0.8`→`--max-missing 0.2`）；**不提供**静默换算的隐藏别名 | 本提交 | CLI 单测 `removed_flags_error_with_the_conversion`；手工：`filter --geno 0.1`、`stats --mind 0.2` 均 rc=2 且提示含目标参数 |
| `reed stats` 的样本级参数 | 隐藏参数 `--max-missing-sample`/`--mind` 一律拒绝并提示先用 `reed filter`（02 §3.2） | 本提交 | 手工：`reed stats x.vcf --max-missing-sample 0.9` rc=2，提示"run `reed filter --max-missing-sample X` first" |
| spec 字段名 / `to_args` / provenance | `TailRules.max_missing`→`min_call_rate`、新增 `max_missing_alleles`、`FilterSpecV2.mind`→`min_sample_call_rate`；`to_args()` 渲染 `--max-missing/--max-missing-count/--max-missing-sample`；`##reed_filterCommand` 随之 | 本提交 | `translate::filter::tests::{spec_to_args_round_trips_through_the_flags, sample_call_rate_round_trips, static_subset_and_sample_rate_round_trip}`、`spec::render::tests::{renders_a_combined_sample_selection, renders_the_documented_flags}`；gold 输出的 `##reed_filterCommand` 含 `--max-missing 0.9` |
| gold 脚本换算（02 §6） | `gold/scripts/run_reed.sh` 用与 `run_vcftools.sh` 相同的 awk 精确换算 `1−SITE_MISS_MAX`（`gold/` 未纳入版本控制，改动在工作区） | 本提交 | 全量 `run_reed.sh` 正文 md5 **7aedc59ac798729a2674e6e7ef19adfc**（与 V0–V6 记录相同）；`check_parity.py out/reed.vcf.gz --gold out/full/gold.vcf.gz` **PASS**（1,768,759 位点，jaccard 1.0，基因型不一致 0） |

W0 门（02 §7）六条全部通过；`cargo test --workspace` 全绿，`cargo build --release --workspace --all-targets` 0 警告，`cargo clippy --release --workspace --all-targets` 0 error。

## 未完成 / 偏差

| 项 | 状态 | 说明 |
|---|---|---|
| `sample.max_missing` 的形态 | 偏差（已登记） | 现实现没有样本级漏斗条目（报告只有 head/tail/mask 的**记录**计数），故 `sample.max_missing` 不是 `FilterReport.funnel` 的元素，而是报告新增的 `sample_rule: {rule, threshold, n_dropped}`（`reed-filter/src/funnel.rs`）；`after_mask.max_missing`、`after_mask.max_missing_count` 按 §5.3 作为 tail 漏斗 id |
| W1–W6 | 待做 | 见 `06-roadmap.md` |

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
