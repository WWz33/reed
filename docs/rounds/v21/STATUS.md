# 实施状态（2026-10-28）

> 纪律（沿 v20）：每条**先写复现缺陷的失败测试**（红），再最小修复（绿）。
> 测试名 `mN_*`（本轮 must-fix）/ `sN_*`（评审 should-fix 编号）。
> commit 列＝未提交（与本树既有 WIP 同树）。

## Must fix（M1–M11）

| 项 | 文件 | 内容 | 验证（修前红→修后绿） |
|---|---|---|---|
| M1 | reed `genotype.rs`、`pca/mod.rs` | `read_genotype_matrix` 的 engine 调用改走 `InputJob`（其余四个入口都包装；裸 `PcaJob` 用 `Newline` 帧把你 BCF 二进制按 `\n` 切）。全命令（pca/dapc/distance）受益 | `e2e_bcf::m1_genotype_matrix_is_format_independent`（修前 BCF 报 malformed record）；真实 `sub.bcf` 上 pca/dapc/distance 通过 |
| M2 | reed-io `chunk/region.rs` | `IndexQuery` 加 `has_reference_names`/`tid_of`；`query_range` 优先用**索引自身名字表**取 tid，仅在索引无名字表（contig-id CSI）时退回头顺序；索引命名而头也声明却无记录的 contig → 空（非错） | `chunk::region::tests::m2_query_range_resolves_tids_through_the_index_names`；真实 `.tbi`/重建 `.csi`/干净 `.csi` 三形态 |
| M3 | reed-io `chunk/mod.rs`、`chunk/region.rs` | `ChunkBytes` 加 `truncate_to`：区域源的**末尾块**在索引结束虚拟偏移（记录边界）处截断膨胀文本；范围终点恰在块边界时不再读该块。消除 `input ends inside a record` | `chunk::region::tests::m3_chunk_truncates_to_the_range_end_inside_the_last_block`（修前 `truncate_to==None` 断言失败）、`m3_block_aligned_range_end_needs_no_truncation`；真实 `subvcf` 40 区间 36→0 失败 |
| M4 | reed-cli `translate/common.rs`、`translate/filter.rs` | 新增 `canonical_region_line`/`region_tokens`：`-r/-R/-t/-T` 的文件按行解析，≥3 列=BED（0-based 半开→1-based 闭），2 列=`CHROM POS`，1 列=原样；`^` 排除保留；空/注释行跳过 | `translate::common::tests::m4_*` ×2；真实 `-T/-R real.bed` == bcftools 124 |
| M5 | reed-filter `spec/region.rs` | `Region::parse` 接受裸位置 `CHROM:POS`（`start==end`） | `spec::region::tests::m5_bare_position_is_a_single_site`、`m5_non_numeric_position_is_rejected`；真实 `-t Chr05:16004641` 命中 |
| M6 | reed-io `vcf/site.rs` | `ContigTable::id_of` 快路径加 `last.1 != ContigId(u32::MAX)` 哨兵守卫（空 CHROM 不再命中种子） | `vcf::site::tests::m6_empty_chrom_never_hits_the_last_hit_sentinel`；真实首行空 CHROM + `-W` 由 rc 134 变正常报错 |
| M7 | reed-io `vcf/line_edit.rs`；reed `filter/mod.rs` | 选中零样本时记录去掉 FORMAT 与样本列（8 列，对齐 8 列头）；`build_keep` 后零样本 warn 一次 | `vcf::line_edit::tests::m7_zero_selected_samples_drops_format_and_samples`；真实 `--remove all` 头/记录均 8 列（VCF 与 BCF） |
| M8 | reed `rename/mod.rs` | `--clear-ids` 的互斥块补上 `--set-missing-var-ids`（并镜像进 `validate()`，I/O 前即拒） | 真实 `--clear-ids --set-missing-var-ids` 由"ID 全变 ." 变 usage 错 |
| M9 | reed-graph `tree.rs`、`bootstrap.rs` | UPGMA 合并进**较低槽位**（`(older, newer)`），对齐 R `hclust`；改正错误注释；更新随之变化的 s36 支持向量锚点 | `tree::tests::upgma_tied_distances_follow_hclust_slot_order`（修前根高 0.25→断言 1/1/4/3 失败）；整数矩阵穷举 older 25/25 vs newer 19/25 |
| M10 | reed-core `group.rs`、reed `stats/request.rs`、reed-cli `translate/{popg,convert}.rs` | 删掉第二个解析器 `names_from_tsv`；`Populations::group_names` 改收 `sample_names`+`force_samples`，直接走运行期同一 `assign_groups`；`WindowedContext` 带上 header 样本名 | `e2e_stats::m10_group_names_match_the_runs_forced_assignment`；真实破损 pop 文件由 panic/错标变清晰 usage 错，五组场景正确文件与破损文件输出一致 |
| M11 | reed `pca/mod.rs`、`dapc/mod.rs`、reed-cli `translate/distance.rs` | 门控按路径计真实峰值：`Plain=2nm`、`Stream=3nm`、`DenseMeanImpute=2nm+n²`、`DensePairwise=2nm+2n²`；抽出共享 `chosen_solver`/`matrix_use` | `pca::tests::m11_dense_paths_charge_the_grm_against_the_budget`；真实 2000×50 在 16 MiB 下 dense 由"通过后栈溢出"变 usage 拒绝 |

## Should fix（bug 类，沿用评审 S 编号）

| 项 | 文件 | 内容 |
|---|---|---|
| S11 | reed-cli `translate/filter.rs` | `--hwe` 第二个 token 非 `midp` 报用法错（原来静默忽略） |
| S12 | reed-filter `plan/mod.rs` | `--hwe midp` 的拒绝去掉内部里程碑码 "M7"，改面向用户的说明 |
| S13 | reed-filter `plan/mod.rs` | `--maf > --max-maf`、`--mac > --max-mac` 报用法错（原来静默空文件） |
| S14 | reed-cli `translate/filter.rs` | `--indep-pairwise` 的 r² 必须有限且 ∈[0,1]，window/step ≥ 1 |
| S15 | reed-cli `translate/{filter,rename}.rs` | `-O` 级别越界改用法错（原来 exit 1 双 `io error` 前缀） |
| S16 | reed-cli `translate/{filter,rename}.rs` | `-W` 在 stdout 或非 BGZF 输出时报用法错（CLI 层，避免伤库调用方） |
| S17 | reed-cli `translate/popg.rs` | `--proj/--sfs-int/--block-size/--pbs-normed` 缺伴随项报用法错 |
| S18 | reed-reduce `xval.rs`、reed-cli `translate/pca.rs` | `n_rep == 0` 库层 `XvalError::NoRepetitions` + CLI usage 错 |
| S19 | reed-cli `translate/{pca,stats,convert}.rs` | 空/空白 `-o` 前缀报用法错（原来在 CWD 写隐藏文件） |
| S20 | reed `rename/mod.rs` | 空 `--set-id/--set-missing-var-ids` 模板报错（原来写空 ID 列） |
| S21 | reed `stats/writers/windows.rs` | 空窗口 θ 打 `NA`（与 π/dxy 一致；原来打 0） |
| S22 | reed `stats/request.rs`、`stats/driver.rs` | 零分母的 D/f4-ratio 仍出一行（NA），`FStatAcc::row` 改不可失败返回 |
| S29 | reed-core `site.rs`(已有分类)、reed-filter `spec/mod.rs`/`plan/mod.rs`/`eval/head.rs`、reed-cli | 新增 `SiteRules.just_acgt` + `HeadRule::JustAcgt`；裸 `--snps-only` = 任意 1-bp 替换，`=just-acgt` 限 A/C/G/T |
| S30 | reed-cli `translate/filter.rs` | `--max-mac 0` 报用法错（0 是内部 off 哨兵，plink2 语义不同） |
| S32 | reed-cli `translate/{common,filter}.rs` | `-S` + `--keep` 同给报错（原来忽略后者） |
| S33 | reed `rename/mod.rs`、reed-cli `translate/filter.rs` | 空 `--add-meta`/`-s`/`--fill-tags` 报错 |
| S34 | reed-cli `cli/mod.rs`、`translate/{pca,distance}.rs` | `pca`/`distance` 加 `--verbosity`（0–3），摘要走 `shows_summary()` |

## 全量验收

- `cargo test --workspace`：**exit 0，839 个测试通过、0 失败**（含 373 s 的
  noodles 差分套件与线程确定性套件）。
- `cargo clippy --workspace --all-targets`：**0 告警**。
- 数值锚点变动仅一处且可归因：`reed-graph` 的 s36 支持向量（旧值即旧 tie 规则
  的产物，已在测试内注明 doc_v21 M9）。
- 真实数据终检：M1–M11 + S 系列的关键命令在 `subvcf.vcf.gz` / `clean.vcf.gz` /
  `clean.bcf` 上结果与 bcftools/参照一致。

## 偏差记录

1. **M9 的真实矩阵残余**：本机 348×348 矩阵上 reed 与 R `hclust` 的合并高度序列
   逐位一致，但 clade 组成仍差 100/237。归因：R 的最近邻扫描在并列时的取序与
   全对扫描不同——槽位规则（本轮修的）已对齐，扫描策略未对齐。记录为残余，
   roadmap 跟踪。
2. **M11 的门控单位**：`--memory` 是 MiB（`m << 20`）；报告 #10 的文字复现需按
   MiB 读，本轮据此重测。
3. **S16 放在 CLI 层**：库的 `FilterOptions::default().write_index` 是 `Csi`
   （库默认要求索引），库级"非 BGZF 即拒"会伤到未显式指定压缩的库调用方，
   故只在 CLI 翻译层拒绝。
4. **S29 语义改动**：原 `default_missing_value="just-acgt"` 使裸 `--snps-only`
   等价于 `=just-acgt`；改为裸=任意 1-bp、`=just-acgt` 才限 ACGT，与 bcftools/
   plink2 对齐（`docs/spec/cli.md` 同步）。

## 未完成（转 06-roadmap.md）

- 评审与文档中其余 should-fix（验证面、CI 门、测试金字塔）、nice-to-have 全部不动。
