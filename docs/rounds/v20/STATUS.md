# 实施状态（2026-10-24）

> 对照 01-review.md 的 R1–R14（第一批）与 06-roadmap.md 的 S 系列（第二批）。
> 全部条目按 v2 06 §1.2 纪律执行：先写复现缺陷的失败测试（红），再最小修复
> （绿），测试名 `rN_*` / `sN_*` / `sec1_*` 描述缺陷本身。
> commit 列 = **未提交**：v20 改动与 doc_v18 的 WIP 同树，由用户决定切分
> （按各表"文件"列 + 测试名 + 源码内 `doc_v20` 注释定位）。

## 已完成

| 项 | 文件 | 内容 | 验证 |
|---|---|---|---|
| R1 | reed-filter `eval/head.rs`、`plan/mod.rs`(doc) | `-f` 按 `;` 拆分逐值匹配（零分配）；`.` 仅显式列出时匹配；枚举文档改正 | `r1_filter_list_matches_any_semicolon_separated_value`（修前 FAILED）；`-p reed-filter` 30/30；`-p reed --test e2e_filter` 26/26 |
| R2 | reed-io `vcf/format_layout.rs` | `intern` 快路径加 `last.1 != usize::MAX` 哨兵守卫；空 FORMAT → 无 GT 布局 | `r2_empty_format_column_does_not_panic`（修前 panic `len is 0 … index 18446744073709551615`）；`-p reed-io` 74+全部集成套件 |
| R3 | reed `src/index.rs` | `IndexSink::push` 顶部 tid 边界检查，报含 tid 的 job 错误（覆盖 `add` 与排序错误路径） | `r3_sink_rejects_tid_beyond_the_header_contigs`（修前 panic `len is 1 … index is 7`） |
| R4 | reed-cli `translate/{common,filter,rename,convert}.rs` | 共享 `split_output_type`；空 `-O` 报 usage 并列出合法拼写；convert 的 1–9 级别校验保留 | 三个 `r4_empty_output_type_is_a_usage_error`（修前 panic `end byte index 1 …`） |
| R5 | reed-cli `translate/{common,filter,rename,convert,index}.rs`、`tests/cli_contract.rs` | 共享 `refuse_output_eq_input`（canonicalize/absolute，`-o -` 跳过，symlink 缺口文档注明）；index 在 exists/`-f` 门之前 | 四个 `r5_*`（修前实测：filter/rename/convert/index 全部 exit 0 且毁掉输入副本；修后 exit 2 + 消息 + 输入逐字节不变）；`-p reed-cli` 86/86 |
| R6 | reed `src/genotype.rs`、`pca/mod.rs`、`dapc/mod.rs`(1 行) | `malformed_line` 镜像 `IoError::MalformedRecord`（列数不符 / 不可解析行）；空行、`#` 仍跳过；`sites` 死 Vec → `n_sites` 计数器，POS 解析随之删除 | `r6_*` ×3（修前走 Ok 分支）；`-p reed` 20 目标全绿 |
| R7 | reed-stats `src/sfs.rs` | 众数锚定双向递推（k* u128 求值）；i64 分母、归一化、1e-15 尾部保留不变 | `r7_projection_weights_do_not_overflow_at_large_n`（修前 `== Vec::new()`）+ `r7_weights_match_exact_hypergeometric`（8 组小参数 vs u128 精确 pmf，1e-12）；既有锚点逐位不变；`-p reed-stats` 59/59 |
| R8 | reed-reduce `src/extra.rs` | 梯度提取为 `tsne_grad`（行为保持）后改 `* qij`；无补偿常量需调 | `r8_tsne_gradient_uses_normalized_q`（修前 got/expected = sum_q = 1.741254914 精确）；`-p reed-reduce` 65/65（`tsne_separates_blobs` 不变通过） |
| R9 | reed-distance `src/standard.rs`、`matrices.rs`(测试)；reed-reduce `src/pca.rs`(测试夹具) | 退化列写中心化 0；文档注明与 flashpca 的有意分歧 | `r9_fixed_alt_column_standardizes_to_zero_not_the_mean`（修前 `left: 2.0`）+ GRM 等价断言；`grm_diagonal` 夹具换真正变化列（旧夹具 100% 退化列，旧锚点即 bug，测试内注明） |
| R10 | reed `src/pca/mod.rs`、`dapc/mod.rs`；reed-cli `tests/e2e_ordination.rs` | `solve()` 变 `Result`；`chosen` 处拒绝 stream+pairwise（含 auto）；假文档句删除 | `r10_*` ×2（修前返回 Ok solver=stream）；`e2e_ordination::dapc_var_agrees_between_solvers` 旧锚即 bug，改锚 `--missing meanimpute` 并注明；e2e_pca 6/6、e2e_ordination 23/23 |
| R11 | reed `src/rename/mod.rs`、`tests/e2e_rename.rs` | 互斥块新增 `--set-id` × `--update-name` 一臂 | `r11_set_id_conflicts_with_update_name`（修前 Ok 报告）+ `r11_set_missing_var_ids_keeps_update_name_reachable`（合法组合钉住）；`e2e_rename` 18/18 |
| R12 | reed `src/convert/writers.rs`(+`mod tests`) | `pos0 < pending` 报 `input is not contig-sorted`（沿用 stats sink 措辞）；接缝同位点 last-wins 保留 | `r12_write_msa_refuses_regressing_frames`（修前 20 bp 写出 28 碱基，输出留档在测试注释）+ `r12_write_msa_keeps_seam_duplicates_and_ascending_blobs`；`e2e_convert` 17/17 |
| R13 | reed-distance `src/matrices.rs`、`Cargo.toml`；`Cargo.lock`(1 行) | 全估计量按行 rayon 并行 + 串行 scatter + dosage 行预计算；浮点不跨线程；死 `miss` 向量删除 | `r13_estimator_results_thread_invariant`（1 vs 4 线程池 `to_bits()` 相等，含缺失 + 退化列）；25 个 reed-distance 手算锚点全部不变；e2e_ordination 23/23（含 `--tree/--bootstrap` 重建矩阵） |
| R14 | reed-core `src/grid.rs` | span 成员表（构造期输入序）+ `partition_point` 二分 + 成员精确过滤 | `r14_stab_matches_linear_reference`（旧线性实现作 oracle 穷举；并抓出成员表设计的嵌套区间错误）+ `r14_stab_large_bed_fast`（200k×2000 探测毫秒级）；`-p reed-core` 106/106；`e2e_stats::bed_file_windows_match_pixy_layout` 通过 |

## 第二批（已完成，2026-10-24 续）

> 应用户指令与 v20 同轮（不开新轮）：06-roadmap.md 中的 **bug 类**条目——
> 错误结果 / 崩溃 / 安全，不含 CI 门、验证类小项与 lean。编号沿用评审 S 编号。
> 测试名 `sN_*` / `sec1_*`。

| 项 | 文件 | 内容 | 验证 |
|---|---|---|---|
| S58 | reed-filter `expr/mod.rs`；reed `filter/record_pipeline.rs`×2 行 | `info_or_gt` 按 `info_trusted` 门控（对齐 `af_of`/bcftools vcfview.c:341-350）：cohort 变化后表达式 AC/AN/AF/MAC/MAF 从基因型重算；`SiteValues::from_line` 改收 `info_trusted` 参数，删两处事后覆盖；直接 `INFO/x` 读取与 head 规则不受影响（head 规则仅在 cohort 未变时编译） | `s58_untrusted_info_recomputes_from_genotypes`（修前返回陈旧 MAF 0.5）；`-p reed-filter` 31/31；e2e_filter 26 + e2e_info 17 |
| S28 | reed-reduce `umap.rs` | 学习率去掉 `.max(0.25)` 下限，抽出 `sgd_alpha`，照 umap-learn 线性衰减到 0；注释改引 `layouts.py` | `s28_alpha_schedule_decays_linearly_to_zero`（末轮恰 0.0）；CLI 无 umap 出口，无锚点变动 |
| S29 | reed-reduce `extra.rs`、`grm_reduce.rs`、`lib.rs` | 困惑度搜索非收敛行也写入最终权重（不再静默清零）；`tsne` 返回 `Result` + 类型化 `TsneError::{PerplexityBelowOne, PerplexityTooLarge}`（`<1 或 ≥ n−1` 拒绝，不静默钳制）；`grm_tsne` 透传 | `s29_impossible_perplexity_is_rejected`（修前 Ok 噪声）、`s29_nonconvergent_rows_keep_their_final_beta_weights`（几何比值断言：主导亲和对 ≥4× 近于零亲和对——纯斥力给不出该几何）、`s29_converging_case_embedding_is_unchanged`（精确锚 73.6595911798…） |
| S31 | reed-reduce `umap.rs:192`、`graph.rs:191`、`extra.rs:338` | 三处 kNN 排序 `partial_cmp().unwrap()` → `f64::total_cmp`（NaN 排尾、确定性）；graph.rs:326（isomap）非评审点名点，留 backlog | `s31_*` ×3（NaN 不再 panic、有限序正确、NaN 不入 k 窗） |
| S32 | reed-reduce `xval.rs` | 组 id 计数前校验 `groups[i] < ng`，返回类型化 `LdaError::GroupIndexOutOfRange` | `s32_group_id_out_of_range_is_a_typed_error_not_a_panic`（修前 panic） |
| S36 | reed-graph `bootstrap.rs` | `branch_keys` 跳过全 1 键（根/中心节点不再恒 100% 支持）；两构造器根都是最大节点 id，`to_newick` 标签索引天然对齐；support 数量对齐 ape（n−2/n−3）；ape-parity 测试去掉根的 workaround | `s36_*` ×3（修前 UPGMA 根 `)100;`，其余 support 逐一不变）+ `clade_counts_match_ape` 断言根键不存在；`-p reed-graph` 35/35 |
| S37 | reed-cli `translate/distance.rs`、`Cargo.toml`(+rayon) | discordance 的 bootstrap 闭包改捕已构建的 `DiscordanceOptions`（修前恒 default：Count 树配 Probability(40) 的 support）；`--threads` 经 scoped `ThreadPoolBuilder.install` 接通估计量与 `graph::support`（R13 收尾） | e2e 三测：`s37_bootstrap_discordance_respects_mismatch_count`（修前 92/33/25 → 修后 58/25/42）、`s37_bootstrap_error_probability_reaches_the_replicates`（修前 E=1 与 E=40 support 相同）、`s37_thread_count_does_not_change_distance_output`（1 vs 4 线程 `.nwk`/`.gtcheck` 逐字节同）；e2e_ordination 26/26 |
| S76 | reed `pca/mod.rs` | 剂量矩阵分配前对照 `--memory`：超限报 usage（含所需/预算 GiB、指向 `reed filter --indep-pairwise`）；`run_pca_opened` 改道 `read_genotype_matrix`（顺带消掉 P2 的 ~40 行重复块，dapc/distance 同路继承检查） | `s76_dosage_matrix_beyond_the_budget_is_a_usage_error`（修前 Ok）+ 舒适预算对照；e2e_pca 6/6、e2e_ordination 26/26 |
| S68 | reed `filter/record_pipeline.rs`、`filter/mid_stage.rs` | record_pipeline 五处 `keep.get(i).unwrap_or(true)` → 直接索引 + `debug_assert_eq!(keep.len(), n)`；mid_stage 加位宽锚定断言；契约钉死："短 keep 是 debug 期捕获的编程错误，必须全长" | `s68_*` ×2（修前"test did not panic"，sketch/单遍等价 e2e 守卫保持绿） |
| S46 | reed-stats `fst.rs` | **评审有误**：multi 在 r=2 同样跳过（`n_C=2n₁n₂/(n₁+n₂)` 任一 nᵢ 取 0 即 0，既有 `n_c≤0` 门先触发）；修法改两条路径一起对齐 an∈{0,1} → keep-c 分支（pair 触发条件放宽为 `n_a<2‖n_b<2`；multi keep-c 门扩为 `an∋0 ‖ (n∋0 ∧ n_c≤0)`，r≥3 的 `an=[1,6,6]` 类行为不变）；oracle 冻结旧跳过行为，parity 测试按"oracle 分量 + 角落位点 keep-c 项回加"调整并注明 | `s46_wc84_pair_keeps_c_for_lone_half_call_site`（修前 tally (0,1)）、`s46_wc84_pair_equals_multi_at_r2`（5 位点混合缺失，<1e-12）；既有锚点 0.125/0.16981… 不变；e2e_fst_global/e2e_fstats 全绿 |
| S79 | oracle `estimators/pbs.rs` | 模块文档改为与代码一致的 per-pair 规则（三群交集的说法是假的，`:235-236` 本有正确描述；`n_sites`/`n_skipped` 字段文档一并改正）；oracle 代码零改动 | `s79_pbs_accumulates_per_pair_not_intersection`（C 全缺位点 ab 对仍积累，`fst_ab`=1/27 而非交集规则的 1/5，且等于独立 `fst(A,B)`——旧文档否认的等式） |
| S39 | reed-io `bcf/native_encode.rs` | FORMAT 整数越界由钳制（DP=3e9 → 2147483647）改为 missing 哨兵 + `warn_once_out_of_range`（对齐同文件 INFO 路径与 htslib vcf.c:4004-4013）；越界值不再把列撑宽 | `s39_out_of_range_format_int_is_missing_not_clamped`（修前 2147483647，INFO 对照本就 `.`） |
| S40 | reed-io `chunk/mod.rs` | `split_text` 变 `Result`：`scan_start > text.len()` 报 `MalformedBgzf`（损坏索引虚拟偏移），不再切片 panic；界内偏移行为不变 | `s40_*` ×2（修前 `range start index 64 out of range for slice of length 5`） |
| S43 | reed-io `bcf/native_encode.rs`、`bcf/render.rs` | POS > i32 范围（文本→BCF）与负 pos0（BCF→文本/BCF）由静默钳制改为 `MalformedRecord`（POS=2³¹ 边界仍合法，已测）；`reed/src/index.rs:337/409` 的同型 `unwrap_or(0)` 非本条范围，留 backlog | `s43_*` ×3（修前 Ok 且写 2147483648 / POS=1） |
| S42 | reed-io `vcf/site.rs`（新共享 `indexed_span_len`）+ `bcf/native_encode.rs`；reed `filter/emit.rs`、`index.rs` | 符号等位 + INFO/END 的索引跨度 = `END−POS+1`，三处 span 站点收敛为一个 `SiteView::indexed_span_len()`（非符号记录保持 REF 长，SNP/indel 不变）；BCF 索引路径经编码器值自动继承。**测试设计要点**：原建议夹具（DEL 100-500 + 查 300）不判别——CSI 叶 bin 宽 16384,300 与 100 同桶;实际用跨 bin 夹具（DEL 100 END=30000 + SNP 20000,查 20000,对照 noodles-csi reg2bin 源码） | `s42_indexed_span_len_uses_end_for_symbolic_alts`（11 例）、`s42_symbolic_end_extends_bcf_rlen`（rlen 字节 401）、`s42_region_query_finds_spanning_symbolic_del`（e2e_index：修前查不到 DEL，修后命中；e2e_bcf 12、e2e_filter 26 全绿，正常数据字节不变） |
| SEC-1 | reed `input.rs` | stdin 落盘改 `create_new`（O_EXCL：文件与符号链接都拒、不跟随不截断）+ 计数重试上限 16；命名方案不变（排他性本身即防预置,无需随机名）;Drop 清理经核安全（unlink 不跟随符号链接） | `sec1_*` ×3（修前受害者文件与符号链接目标**均被截断为空**、耗尽时覆盖 -1 号候选；修后全不触碰）；stdin_streaming 13/13 |

## 全量验收

**第一批（R1–R14，2026-10-24）**：`cargo test --workspace` exit 0，全部目标
0 失败（57 个测试二进制，含 371.8 s 的 noodles 差分套件与线程确定性套件）；
既有锚点仅两处变动，均为"旧锚即 bug"（见 R9/R10 行）。

**第二批（S 系列 16 条，2026-10-24 续）**：逐 crate 全绿（reed-filter 31、
reed-reduce 73、reed-graph 35、reed-stats 10 套件、reed-oracle 全部、
reed-io 82+集成、reed 20 目标、reed-cli 5 目标含 e2e_ordination 26）；
`cargo clippy --all-targets` 各改动 crate 无新告警。数值锚点变动仅三处且
均可归因：S37 的 discordance bootstrap support（旧值 = 错误度量）、
S46 的 oracle parity 角落（oracle 冻结旧跳过行为，测试内注明回加项）、
S79 无数值变动。工作区级全量复跑见下方"最终门"。

**最终门（两批之后，2026-10-24）**：`cargo test --workspace` **exit 0，
全部测试目标 0 失败**（含 noodles 差分套件与线程确定性套件）；
`ci/check-doc-drift.sh` 通过。两批合计 30 条修复（R1–R14 + S 系列 16 条），
全部测试先行，与 doc_v18 WIP 同树未提交。

- 并行修复期间的 cargo 以 `flock /tmp/reed-v20.lock` 串行,最终全量跑无锁。
- 其余全部数值锚点（手算 fixture、gold 对拍、oracle parity）零变动。

## 偏差记录

1. ~~**R13 `--threads` 未全链路接入**~~ → **已由第二批 S37 补齐**（scoped
   pool 接通估计量与 `graph::support`，1 vs 4 线程输出逐字节一致锁定）。
   R13 的流式 sites 分块、n×m 冗余临时全消仍留 backlog。
2. **R5 的符号链接缺口**：输出不存在时用 `std::path::absolute`（词法
   归一,不解析符号链接）——输出目录经 symlink 指向输入目录的形态不设防,
   已在 `refuse_output_eq_input` 文档注明。
3. **R6 顺带删除 `GenotypeBlock.sites`**：评审 backlog 第 67 条（死内存
   ~1 GB @ 10⁷ 位点）与本条修复面重叠,一并完成;`GenotypeBlock` 为 crate
   私有类型,无 API 兼容问题。
4. **rustfmt**：工作树存在**先于本轮**的 fmt 漂移（用户 WIP 文件 +
   rustfmt 版本差异）,两批新增代码未引入新 diff,未代跑 `cargo fmt`。
5. **S46 与冻结 oracle 的有意分歧**：oracle 实现旧的 an=1 跳过；parity
   测试按"oracle 分量 + 角落位点 keep-c 项回加"调整（等价于 oracle 自己
   的 an=0 分支算术），每窗口仍受门控——冻结库不改源码，改的是期望值的
   推导并注明 doc_v20 S46。
6. **S28–S32 的中断运行残留**：第二批开工时树里已有一次被中断运行的
   半成品（doc_v20 S* 注释 + s28–s32 测试齐全但 s29 断言错误）；该代理
   对照 HEAD 逐条审计、修复断言（几何比值替代绝对尺度）、删除遗留调试
   探针后全绿——最终状态经过完整复验。

## 未完成（转 06-roadmap.md）

- should-fix 剩余条目（第二批已清掉 bug 类 16 条）：验证面/CLI 批、
  stats 批余项、core 批、filter 测试批、input 批余项、测试金字塔、lean
  清单——仍按优先级留在 06-roadmap.md。
- R2 同型隐患：`ContigTable::id_of` 空 CHROM 哨兵（01-review R2 段）。
- 第二批新发现的 backlog：
  - `reed-reduce/src/graph.rs:326`（isomap 的 kNN 选择）仍是
    `partial_cmp().unwrap()`——第四个 NaN 排序点，非评审点名。
  - `reed/src/index.rs:337/409`：BCF 索引路径的负 pos0 仍 `unwrap_or(0)`
    （外来畸形 BCF 直接索引时记 POS=1）——S43 只修了 reed-io 两侧。
  - t-SNE 在 `learning_rate: 200`、微小 n 时坐标可涨到 ~1e154（HEAD 既有、
    保形；步长控制在 roadmap"明确不做"里）。
  - `warn_stale_info` 的文案对尾部表达式已过度保守（S58 后尾部与
    `--maf` 同口径重算），措辞微调属验证面批。
- CI 层（fmt/deny 缺失、data-gates 不构建被测二进制）属 roadmap P0,
  但改动面在 `.github/` 与脚本,留待后续批次。

## 第三批：对前两批的复审修复（2026-10-25）

> 两批修复完成后做了一次只读复审（3 个复审代理 + 主代理亲自验证关键数值）。
> 结论：**R8 把正确的梯度改错了**，另有几条修复不彻底或参考实现说法与实测
> 相反。以下全部由主代理串行修复，每条都有可失败的回归测试。
> 判定手段：对 t-SNE 的 KL 代价函数做中心差分（独立判据），用本机 R 5.8.1 +
> ape 实跑 `prop.clades`/`write.tree`，对照工作区 htslib 源码核对 `get_rlen`
> 与端粒坐标。

| 项 | 文件 | 内容 | 验证 |
|---|---|---|---|
| R8 回退 | reed-reduce `extra.rs` | **改回原始核** `4(p−q)·w·(y_i−y_j)`（前一批改成了 `·q`，梯度小了 `sum_q≈n²` 倍，n≈180 时嵌入停在初始化量级）。新增 `tsne_kl_cost` 作为测试专用代价函数 | `r8_tsne_gradient_matches_the_kl_cost`（有限差分：原始核误差 5.7e-11/相对 3.9e-10，归一核 8.5e-2/0.58）。**原测试是 `tsne_grad` 的逐行复写，恒真，已删除** |
| S29 加深 | reed-reduce `extra.rs` | beta 从 `1/mean(d²)` 起搜（剂量距离 d² 是 O(位点数)，固定 beta=1 时 `exp(-d²)` 下溢，整行清零）；`sum==0` 不再 `break` 而是压低上界继续搜 | `s29_genotype_scale_distances_do_not_zero_every_row`（12 样本、d²~1e2–1e3，修前 sum\\|y\\| 塌到初始化量级）+ 收敛锚点按新行为重钉（2.264348394311807e3） |
| S42 补漏 | reed `filter/single_pass.rs`、`rename/mod.rs`；reed-io `vcf/site.rs` | 默认单遍写入器与 rename 也走 `indexed_span_len`（前一批只改了 sketch 的 `emit.rs`，`reed filter -W` 不带 `--mind` 时仍写 rlen=1）；span 规则改为 **htslib `get_rlen`：END 与 ALT 形状无关**（`ALT=.` 的 gVCF 块、`G,<DEL>` 多等位都延伸） | 新增 e2e `s42_filter_write_index_extends_symbolic_deletions_to_end`（修前只查到 SNP，已用临时回退验证）；11 例单元测试改为 htslib 规则 |
| S36 回退 | reed-graph `bootstrap.rs` | **恢复根节点支持**：本机 ape 5.8.1 实测 `prop.clades` 返回 Nnode 个值（含整棵树），`write.tree` 也写根标签；前一批丢弃根并声称 ape 同样丢弃，是**事实错误**。support 长度回到 `n_internal` | `s36_root_clade_is_a_branch_key`、`s36_bootstrap_newick_labels_the_root`、`support_length_matches_ape_clade_count`、`s36_supports_are_the_ape_vectors`；`clade_counts_match_ape` 断言根键存在 |
| S76 加深 | reed `genotype.rs`、`pca/mod.rs` | 预算上限移入 `CollectSink::accept`（逐块检查，超限即中止扫描）；事后检查改按峰值 `2×矩阵` 计（块 + 矩阵 + dense 的标准化副本） | `s76_*` ×2（修前扫描到 96 KB 才拒绝；现扫描中即拒） |
| S58 补漏 | reed-filter `plan/mod.rs` | `Expr::Info(name)` 且是内置字段（`INFO/AC` 等）纳入 `expr_reads_gt_fields`——`--set-gts` 掩码阶段没有基因型值，前一批让它编译通过并静默求值为 Missing | e2e `format_expressions_need_set_gts` 新增一段（修前返回 Ok 报告、掩码 20 个样本槽，已用临时回退验证） |
| S43 修正 | reed-io `bcf/native_encode.rs`、`bcf/render.rs` | **pos0 = −1 是 htslib 的端粒哨兵**（`v->pos == UINT32_MAX` → POS 0，vcf.c:1649），必须原样传递；只有 < −1 才是畸形。前一批把 −1 也拒绝，会拒掉 bcftools 自己写的 BCF | `s43_telomere_pos0_is_copied_through`（新增）+ 两个旧测试改断言 |
| S39 文案 | reed-io `bcf/native_encode.rs` | 越界警告按列族命名（FORMAT 值不再说 "Integer INFO value"） | 既有 `s39_*` 保持绿 |
| R1 语义 | reed-filter `eval/head.rs` | 列出的 `.` 只匹配整列为 `.`（bcftools 把 `.` 映射为 filter id −1，`LowQual;.` 不应通过 `-f .`，实测 bcftools 1.24）；`spec/mod.rs` 与 `docs/spec/architecture.md` 的"`.` 不算通过"旧文案同步改正 | `r1_*` 增 3 个断言 |
| R13 内存 | reed-distance `matrices.rs` | 提升的 dosage 由 f64 改 i8（值是 0/1/2），循环内加宽——2000×1M 省 16 GB；数值逐位不变 | 25 个手算锚点全绿，`r13_estimator_results_thread_invariant` 保持 |
| R10 提前 | reed `pca/mod.rs` | `stream + pairwise` 的拒绝移到读输入之前（`n_kept` 由 header + keep mask 即可算出，auto 的选择也在这里判定） | `r10_*` ×2 保持绿 |
| R11 提前 | reed `rename/mod.rs` | `--set-id` × `--update-name` 的拒绝同时放进 `validate()`（I/O 之前），文件错误不再抢先 | `r11_*` ×2 保持绿 |
| R6 补回 | reed `genotype.rs` | 删掉 `sites` 时丢掉的 POS 校验补回（非数字 POS 报 `MalformedRecord`），`run_chunk`/`run_seam` 各一处 | `r6_run_chunk_errors_on_unparseable_pos`（新增） |
| R12 补漏 | reed `convert/writers.rs` | 乱序检查移到窗口过滤之前（窗口外的回退也拒绝）；`flush_msa_site` 的 `min` 钳制换成 `debug_assert`（死代码） | `r12_regression_outside_the_window_is_still_refused`（修前通过，已用临时禁用验证） |
| S68 一致 | reed `filter/mid_stage.rs` | `debug_assert` 改 `assert`：release 下两条策略不再一个静默丢尾样本、一个 panic | `s68_*` ×2 保持绿 |
| R9 夹具 | reed-reduce `pca.rs` | 对角断言收紧为 `\|diag − 8/6\| < 1e-12`（0.6 的带宽同时接受正确值与 R9 通胀值） | `grm_diagonal_is_one_when_complete` 通过 |

**第三批验收**：`cargo test --workspace` exit 0（全部目标 0 失败）；
`cargo clippy --workspace --all-targets` 0 告警；`ci/check-doc-drift.sh` 通过。
两处旧锚点被证实是 bug 的产物（R8 的 73.66、S29 的 2.2773e3），已在测试注释
里写明来历；S36 的旧测试断言的是**错误的 ape 行为**，已按实测重写。

**仍未处理**（记录在 06-roadmap）：R1 的 `-f` 列表元素含 `;` 静默不匹配、
S58 在 AN==0 时 AC/MAC/MAF 与 `--max-maf` 的口径差异、S29 的 roadmap 文案
（`(n−1)/3` vs 代码的 `n−1`）、R10 的默认 `pairwise`+`auto` 在大 N 下报错
（是否改默认值需拍板）、S42 用输入 INFO 而非写出记录算跨度（`-x INFO/END`
会过度包含）、`e2e_ordination.rs`/`e2e_pca.rs` 被 `.gitignore` 忽略（`crates/*/tests/`
规则，仓库级问题）。
