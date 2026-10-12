# STATUS — doc_v7 缺陷修复

> 本页记录 `01-review.md` §2 缺陷清单与 `06-roadmap.md` Z0–Z4 的实施状态。
> 每条附证据（测试名 / 对拍命令 / 提交）。起始状态全部未修，现全部有结论。

## 结论总览

29 条发现：**24 条改为真修**（含 5 条 P0、9 条 P1、其余 P2 代码/测试改动），
**5 条经核实为"非 bug / 文档或设计说明"**（D7-20、D7-23、D7-25、D7-26、D7-27），
逐条给出了独立提交。

## 缺陷状态

| # | 级 | 状态 | 证据 |
|---|---|---|---|
| D7-1 | P0 | 已修 | `sketch_path_rejects_bcf_output` |
| D7-2 | P0 | 已修 | `indep_pairwise_stdin_spills_without_mind` |
| D7-3 | P0 | 已修 | `plain_and_gzip_inputs_are_refused` |
| D7-4 | P0 | 已修 | 单测 + **plink2 `--keep` 对拍**（Chr01-03：86/86、363/363） |
| D7-5 | P0 | 已修 | `ld_cohort_excludes_removed_samples` |
| D7-6 | P0 | 已修 | `report_mode_rejects_unapplied_filters` + `--sites-file` 生效 |
| D7-7 | P0 | 已修 | **bcftools 1.24 逐行一致**（`summary_idd_qual_psc_match_bcftools`） |
| D7-8 | P1 | 已修 | `f4_ratio_accepts_five_populations` |
| D7-9 | P1 | 已修 | `fstats_match_scikit_allel`（estimate/se==z；admixtools sum 尺度） |
| D7-10 | P1 | 已修 | `f4_ratio_is_the_ratio_of_sums` |
| D7-11 | P1 | 已修 | `skipped_multiallelic_is_present_for_variant_only` |
| D7-12 | P1 | 已修 | `missing_sample_file_errors` |
| D7-13 | P1 | 已修 | `skip_empty_and_levels_apply_to_windowed_stats` |
| D7-14 | P1 | 已修 | `half_calls_are_not_counted_as_genotypes` |
| D7-15 | P1 | 已修 | `bcf_input_refuses_tbi` |
| D7-16 | P1 | 已修 | `truncated_bcf_errors_without_panicking` |
| D7-17 | P1 | 已修 | `value_total`→`value_avg`（CLI 已同时输出 avg/raw，pixy 布局） |
| D7-18 | P1 | 已修 | `indep_pairwise_accepts_integer_r2` |
| D7-19 | P1 | 已修 | `prune_multi_contig_is_thread_count_independent` |
| D7-20 | — | **非 bug** | 多群体 Fst 与 scikit-allel 一致（代码注释 + 文档） |
| D7-21 | P1 | 已修 | `variant_only_denominator_does_not_overflow` |
| D7-22 | P1 | 已修 | 重导出 7 类型；`FilterOptions` 默认压缩 6；`to_args` 补 af_kind/ld_prune |
| D7-23 | P2 | **设计意图** | variant-only 单一估计器（文档说明） |
| D7-24 | P2 | 已修 | `missing_input_names_the_file` |
| D7-25 | P2 | **已知限制（文档）** | 同位点共享剪枝决定（vcftools/plink2 位置键约定） |
| D7-26 | P2 | 已修 | het.rs 文档随 D7-14 重写 |
| D7-27 | P2 | **非 bug + 已修** | `h_bar` 与 allel `count_het` 一致；`pi.rs:57` 溢出随 D7-21 修 |
| D7-28 | P2 | 已修 | `merge_and_reset_carry_the_ld_prune_count` |
| D7-29 | P2 | 部分 | 移除未用的 `assert_cmd`；`--summary` 覆盖缺口由 D7-7 的夹具补 |

## 未覆盖边界（待补测，见 06 §4）

- BGZF 虚拟偏移非零路径（`VoMap::voff`）——未测。
- tabix/CSI 随机访问边界（位置 ≥ 2^29、contig ≥ 2^16、CSI 大区间）——未测。
- gold 管道覆盖扩展（indel / 缺失 QUAL / PSC 全列）——D7-7 已在集成测试层覆盖，管道脚本层面未改。

## 已知既有失败（非本轮引入）

- `reed-engine` 的 `pipeline::tests::stitcher_reassembles_and_caps_fragments`：Y5 seam 上限
  （commit 111dbde），不在 v7 清单内。
