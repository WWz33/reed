# 实施状态（doc_v22，2026-11-03）

> 纪律：每条先写复现缺陷的失败测试，再最小修复；修复后跑全量测试与真实数据
> 回归。测试名以 `v22_` / 语义名标注。

## 修复表

| # | 层 | 文件 | 内容 | 测试 |
|---|---|---|---|---|
| V22-1 | io | `vcf/line.rs` | `LineRef` 接受 8 列 sites-only；`format_text`/`samples_text` 返回空；`n_sample_cols`=0；`samples_start`=行尾 | `line::tests::sites_only_lines_parse_without_format_or_samples` |
| V22-2 | io | `vcf/genotype_scan.rs` | `n_samples == 0` 时跳过样本循环 | `genotype_scan::tests::sites_only_and_trailing_empty_columns_decode` |
| V22-3 | io | `vcf/line_edit.rs` | 无 FORMAT/样本区时输出 8 列 | CLI `v22_accepts_sites_only_haploid_and_trailing_tabs`（断言每行 8 列） |
| V22-4 | io | `vcf/genotype_scan.rs` | 超出 header 样本数的**空**列忽略；非空仍 malformed | 同上（`\t\t` 与 `0/0:7` 两分支） |
| V22-5 | filter | `reed-filter/plan/mod.rs` | header 仅声明 Integer `FORMAT/AD` 时 DP 掩码规则合法 | `e2e_filter::min_dp_falls_back_to_ad_when_dp_is_absent` |
| V22-6 | filter | `reed/src/filter/record_pipeline.rs` | 文本扫描器按编译计划的显式 `IntSpec` 建整数列 | 同上 |
| V22-7 | io | `bcf/columns.rs`（R22-2 后） | BCF 原生 `decode_format_columns` 同语义（主键标量、回退键求和） | `e2e_bcf::bcf_dp_mask_falls_back_to_ad` |
| V22-8 | io | `vcf/format_layout.rs` `vcf/genotype_scan.rs` | `IntSpec { key, fallback }`、`parse_int_sum`、`FormatLayoutCache::new_with_specs` | `format_layout::tests::fallback_key_extends_last_needed_and_resolves_per_record`；`genotype_scan::tests::depth_column_falls_back_to_ad_sum`、`parse_int_sum_matches_htslib_missing_rules` |
| V22-9 | stats | `stats/writers/{reports,summary}.rs`、`histogram.rs` | `--summary` DP/PSC、`--depth`、`--hist dp` 共用深度列；`dp > 0` 才入直方图 | `e2e_reports::summary_depth_falls_back_to_ad`；`summary_idd_qual_psc_match_bcftools` |
| V22-10 | filter/cli | `reed-cli/src/cli/mod.rs`、`translate/filter.rs` | `--ploidy require-diploid\|allow-haploid` | CLI `v22_accepts_sites_only_haploid_and_trailing_tabs` |
| V22-11 | 全层 | `index/genotype/rename/stats/convert/filter` | "fewer than 9 columns" 的错误信息与判断改为 8 列下限 | 现有 R6/D9-4 测试更新；全量绿 |
| V22-12 | io/bcf | `bcf/render.rs`（`render_sample_part`） | `n_fmt == 0` 渲染 FORMAT `.` + 每样本 `.`（htslib `vcf_format` 口径；修前为空 FORMAT 列 + 多一个空样本列） | `bcf::render::tests::fmt_no_format_keys_render_as_missing`；真实 `bcftools view -Ob` 的 `FORMAT=.` BCF 上 `reed filter -Ov` 数据行与 `bcftools view` 逐字节相同 |

## 文档同步

- `docs/spec/cli.md`
  - §3.3 `--min-dp/--max-dp`：深度 = `FMT/DP`，记录无 DP 键时 `FMT/AD` 之和；
    AD 回退标注为 reed 扩展。
  - §7.2 `--depth` / `--summary`：同一深度口径（`calc_sample_depth`，只计
    `dp > 0`）。
- `crates/reed-stats/src/summary/report.rs` 模块注释：AD 回退改为"AD 之和"
  （原 "maximum non-ref AD" 从未实现且口径错误）。
- `CHANGELOG.md`：本轮条目（doc_v22）。

## 验收记录

| 检查 | 结果 |
|---|---|
| `ci/test-parallel.sh`（全量工作区） | exit 0；867 passed / 0 failed（56 个测试二进制，15 s 墙钟；V22-R2 新增 7 项） |
| `cargo clippy --workspace --all-targets` | 0 告警 |
| `ci/check-doc-drift.sh` | exit 0；97 项全部 ok（编码器改名后的路径已同步） |
| `cargo doc -p reed-io --no-deps` | 0 warning（`super::encode_noodles` 断链已修） |
| 真实 BSA（chr1:1-200000）`--min-dp 10` | masked=1285 = 1111 原调用 + 174 原缺失；kept=3075 |
| 真实 BSA `--summary` DP/PSC | 与 `bcftools stats -s -` 逐行相同 |
| `subvcf.vcf.gz` `--summary` DP/PSC | 与 bcftools 逐行相同（有 DP 的文件无行为变化） |
| sites-only `filter -Oz/-Ob`、`index`、`stats --freq --summary` | 全部通过；输出 8 列 |
| BCF 输入的同语义（`e2e_bcf`） | 与 VCF 路径逐字节相同的掩码结果 |
| BCF `FORMAT=.`（`n_fmt == 0`） | `reed filter -Ov fmt_dot.bcf` 数据行 == `bcftools view -Ov fmt_dot.bcf`（10 列，逐字节） |

## R22 架构收口（完成，2026-11-03）

| ID | 结果 | 内容 | 测试 / 验收 |
|---|---|---|---|
| R22-1 | 完成 | `GenotypeMaskRule.depth`；`FilterPlan.int_fields: Vec<IntField>`；`IntSpec` 拥有型；BCF `decode_format_columns(&[IntSpec], …)`；删除 `format_layout::spec_for_key` | `depth_policy_tests` 4 项；全量回归 |
| R22-2 | 完成 | `bcf/mod.rs` 模块地图（读/写/适配）；`render.rs` 的列解码拆到 `columns.rs`；`encode.rs` → `encode_noodles.rs`、`native_encode.rs` → `encode_native.rs`；`ci/check-doc-drift.sh` 的旧路径同步 | 全量测试与对拍绿 |
| R22-3 | 完成 | `LineRef::shape()`（`SitesOnly` / `FormatOnly` / `WithSamples`）与 `validate_sample_cols(expected)`；`GenotypeScanner`、`PcaJob`、`ConvertJob`、`line_edit` 共用同一契约，不再各自数 tab；FORMAT-only 记录在 header 声明样本时按 bcftools 语义报 malformed；尾随空列仍忽略 | `vcf::line::tests::{shape_and_slices_match_a_naive_split,validate_sample_cols_is_total_and_matches_the_contract}`（各 2,000 cases）、`vcf::line_edit::tests::column_shape_is_preserved`、`vcf::genotype_scan::tests::sites_only_and_trailing_empty_columns_decode`、CLI `v22_accepts_sites_only_haploid_and_trailing_tabs` |
| R22-4 | 完成 | `crates/reed-io/tests/format_columns_vs_bcf.rs`：随机 FORMAT 布局 / AD 宽度（含空向量、全缺失、短向量）/ DP 有无 → native BCF 编码 → 原生列解码 == 原文扫描 == 渲染回文本后扫描；列覆盖深度回退、AD 求和、严格 DP、GQ | 500 cases，0 差异 |
| R22-5 | 完成（文档） | `reed-io/src/header.rs` 已是 `TagTable` 单点构造，无需新抽象 | — |
| R22-7 | 完成 | `bcf_vs_noodles` 最小子集、`reed-graph` 深树直构、`ci/test-parallel.sh`、gold 并发与默认 `subvcf` 输入 | 见上表验收记录 |

## V22-R 并行评审修复（2026-11-03）

评审明细与证据见 [02-review.md](02-review.md)。18 项中 17 项属实并修复；
#2（`.gitignore` 吞新文件）按用户要求跳过。

| ID | 层 | 修复 | 测试 / 验证 |
|---|---|---|---|
| V22-R1 | io/bcf | `render_sample_part` 在 `n_fmt == 0 && n_sample == 0` 时不写 FORMAT/样本区；`render_site_line` 同样 | `bcf::render::tests::sites_only_bcf_stays_eight_columns`；真实 sites-only BCF 与 `bcftools view` 均 8 列 |
| V22-R3 | ci | 构建失败、`xargs` 失败、日志数不符、doctest 失败均 exit 1 | `JOBS=abc` 与假 cargo 构建设障均实测 exit 1 |
| V22-R4 | stats | `push_hist` 的 FMT/GQ 样本切片 `take(n_samples)` | `e2e_hist::hist_fmt_ignores_trailing_empty_columns` |
| V22-R5 | gold | `run_pixy.sh` 仅在输入含 `ALT="."` 时加入 `pixy_default` | variant-only → 3 commands / 0 failed；all-sites → 4 commands / 0 failed |
| V22-R6 | spec | `to_args` 渲染 `--ploidy allow-haploid`；DP 掩码仅在 `depth=true` 时渲染为 `--min-dp/--max-dp` | `spec::render::tests::{renders_allow_haploid,depth_less_dp_mask_is_not_rendered_as_min_dp}`、CLI `ploidy_allow_haploid_round_trips`；真实输出头含 `##reed_filterCommand=--ploidy allow-haploid` |
| V22-R7 | filter | `depth=true` 时无论 DP 是否声明都校验 AD 必须 Integer | `depth_policy_tests::{float_ad_fallback_is_rejected_even_when_dp_is_declared,float_dp_depth_mask_is_rejected}`；真实 Float-AD 输入编译期报错 |
| V22-R8 | rename | 删除 sites-only + record-plane 的拒绝守卫；CHROM/ID 只读前 5 列 | CLI `v22_accepts_sites_only_haploid_and_trailing_tabs`（sites-only `--set-id`）；真实 `--rename-chrs` 输出 chrX |
| V22-R9 | cli | 新增 `--min-dp` 对 `GT:AD` 的真实 CLI 测试 | `cli_contract::v22_cli_min_dp_falls_back_to_ad` |
| V22-R10 | test | `REED_SUBVCF_BCF` 覆盖 + 相对工作区上级目录解析 | `bcf_vs_noodles` 全绿；文档改“工作区上级目录” |
| V22-R11 | io/bcf | FORMAT 描述符遍历收敛到 `record::format_fields`，columns/render 共用 | `cargo test -p reed-io` 全绿；proptest 仍 500 例 |
| V22-R12 | docs | v22 README/01-review 的 `spec_for_key`、`native_encode.rs`、`bcf/render.rs` 路径改为最终名 | 文档与 STATUS 一致 |
| V22-R13 | test | `--hist dp` 的 AD 回退 + 零深度排除直接断言 | `e2e_hist::hist_dp_falls_back_to_ad_and_skips_zero` |
| V22-R14 | filter | `impl From<&IntField> for IntSpec`，删除 5 份映射拷贝 | 全量编译/测试 |
| V22-R15 | io | 删除 `IntMode`/`IntSpec::sum`、`FormatLayoutCache::new`、`SampleColumnsError::found/expected/Display`；主键统一标量、回退统一求和 | 全量测试；`IntSpec` 只剩 `key`/`fallback` |
| V22-R16 | io | 扫描器去掉与 `validate_sample_cols` 重复的尾随空列/计数分支，按 `0..n_samples` 解码 | `sites_only_and_trailing_empty_columns_decode`、2,000-case proptest |
| V22-R17 | cli test | v22 CLI 测试断言 haploid/trailing 的实际行内容与列数 | `v22_accepts_sites_only_haploid_and_trailing_tabs` |
| V22-R18 | docs | `e2e_stats` “9 required” → 8；测试计数 860；`encode_native` 断链修复 | `cargo doc -p reed-io --no-deps` 0 warning |
| #2 | — | `.gitignore`/提交项按要求跳过 | 提交时对新文件 `git add -f` |

## V22-R2 第二轮评审修复（2026-11-03，完成）

明细与复现证据见 [03-review-pass2.md](03-review-pass2.md)。按 A–D 批完成；
E 批（有符号整数列）按设计留到后续轮。

| ID | 层 | 修复 | 测试 / 验证 |
|---|---|---|---|
| V22-R2-1 | io/filter | `IntLayout.index/fallback_index`、`Layout.gt_index/float_indices`、`gt_subfield_index` 全部 `u8` → `usize`；`filtered.rs` 的重复实现改为复用 `record_pipeline::gt_subfield_index`（`pub(crate)`）；`line_edit` 的 GT 位置同改 | `vcf::format_layout::tests::keys_past_index_255_keep_their_position`；`vcf::genotype_scan::tests::keys_past_index_255_decode`（GQ@256、AD@257、GT@258）；`vcf::line_edit::tests::sample_pruning_and_gt_masking`（GT@256 输出 `./.`） |
| V22-R2-2 | io/bcf/stats | 整数文本解析收敛到 `info::parse_info_int`（可选 `+`/`-`）再限 `0..=i32::MAX`：超出/负数 → `INT_MISSING`（htslib 口径）；`parse_int_sum`/BCF `int_sum_column` 用 `u64` 累加并 `clamp_int_sum` 到 `INT_MISSING-1` | `genotype_scan::tests::{parse_int_uses_the_htslib_i32_domain,int_sum_never_collides_with_the_missing_sentinel}`；`format_columns_vs_bcf::boundary_ints_match_text_and_native`（`i32::MAX`、`i32::MAX+1`、三个最大值求和、`+50`、负数/超范围分量） |
| V22-R2-3 | io | `parse_float` 的 `inf`/`nan` 语义为 htslib `strtod` 兼容，文档固化；R2-7 的 NaN 哨兵设计选项记录在 03 §8 | `genotype_scan::tests::parse_float_accepts_htslib_inf_and_nan`（含 `parse_int(b"inf") == INT_MISSING`，GQ 是 Integer） |
| V22-R2-4 | cli/stats | `--hist dp` 的 header gate 接受 `FORMAT/DP` 或 `FORMAT/AD`（二者皆无才报 usage）；类型由 `validate_report_tags` 统一校验 | CLI `v22_r2_hist_dp_ad_only_and_stats_type_checks`（AD-only 出 `10/30` 两箱、无 0 箱） |
| V22-R2-5 | filter/stats | `--hist` 的 `parse_int_token` 改为调用扫描器同一个 `parse_int`；INFO/DP（`--summary` 的 `dp_sites`）同改。`+50` 两条路径同值，超 `i32` 两边都 missing | 上述 parser 测试 + `boundary_ints_match_text_and_native` |
| V22-R2-6 | stats | `reports::validate_report_tags` 在 `run_stats_opened` 入口校验：深度消费者需 Integer `DP`（无 DP 时 Integer `AD`），`--hist gq` 需 Integer `GQ`，`--hist dp-site` 需 Integer `INFO/DP`；非 Integer → usage（exit 2）。`--summary/--depth` 在两者皆无时仍允许（深度段为空） | CLI `v22_r2_hist_dp_ad_only_and_stats_type_checks` 的 Float DP/AD/GQ 三 case（含 `--summary`） |
| V22-R2-7 | — | 设计选项记录（03 §8）：`Option<f32>` 列 vs BCF 位型 vs 维持 NaN 哨兵；本轮不改，文本/BCF 两条路径已一致地无法区分 | 文档；`parse_float` 的 doc 注释 |
| V22-R2-8 | — | 负数/有符号整数列留到后续轮；当前负数按 missing（与 `u32` 列一致），已在 parser 文档与测试声明 | 文档；`parse_int` 负数用例 |
| V22-R2-9 | test | 三形态（无 header / 部分记录 / 全部记录无 GQ）经复核不是 bug；既有 `e2e` 行为由本轮 CLI 类型校验补强 | 见 03 §1；GQ Float 的 usage 用例 |

### 验收记录（V22-R2）

| 检查 | 结果 |
|---|---|
| `ci/test-parallel.sh`（全量工作区） | exit 0；867 passed / 0 failed（56 个测试二进制，15 s 墙钟） |
| `cargo clippy --workspace --all-targets` | 0 告警 |
| `ci/check-doc-drift.sh` | exit 0；全部 ok |
| `cargo doc -p reed-io --no-deps` | 0 warning |
| 真实 `subvcf.vcf.gz` `--summary --hist dp --depth` | DP/PSC/GDB 与 `bcftools stats -s -` 逐行相同；`.dp.hist` 最小箱 3（无 0 深度箱） |


## 未做 / 残余

- R22-6（`MalformedRecord` 的结构化 FORMAT/样本上下文）按路线图暂缓。
- R22-3 的行为对齐：9 列 FORMAT-only 记录在 header 声明样本时此前会被扫描器
  当作"一个缺失样本"接受，现与 bcftools 一致地报
  `sample column count disagrees with the header`；`PCA`/`convert` 对尾随空列
  的处理也随之改为与扫描器一致（仅忽略超出 header 样本数的空列）。
- `--min-dp/--max-dp` 的 AD 回退是 reed 扩展：plink2/vcftools 在无 DP 时不做
  回退。仅在"记录 FORMAT 无 DP 键"时生效，有 DP 的 callset 行为不变。
- `reed filter` 的默认 `require-diploid` 保持不动；需要单倍体数据时显式传
  `--ploidy allow-haploid`（`stats` 一直允许半调用）。
- 多等位记录仍受默认 biallelic gate 限制（plink2 风格，非本轮范围）。