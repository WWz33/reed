# doc_v23 实施记录：P0 崩溃/数据损坏与拍板落地（2026-11-03）

> 起点：HEAD `8b54588`（doc_v22 R2）。本文件记录**代码改动**，与
> [01-verification.md](01-verification.md)（复核证据）和
> [STATUS.md](STATUS.md)（条目状态）配套。
>
> 纪律：每条先写会失败的测试，再最小修复，最后全量回归。测试名以语义命名
> （`doc_v23 V23-xx` 写在注释里），便于 review 时逐条对照。

## 落地范围

| 组 | 条目 | 依据 |
|---|---|---|
| P0 崩溃 | V23-01、V23-02、V23-03、V23-05 | 复核实测 SIGABRT/panic |
| P0 静默错数据 | V23-04、V23-06 | 负浮点丢符号、DAPC off-by-one |
| P1 | V23-08、V23-09（拍板 A） | NaN 排序 panic、跨 backend 变号 |
| 拍板 | V23-14（C）、V23-21（A）、V23-50（A）、V23-47（A） | README 拍板表 |
| 第二批 | V23-59（决策 A）、V23-60/61（决策 B）、V23-10/11/12（决策 C ③） | 讨论后拍板 |
| 第三批 | V23-13/15/16/17/18/22/27/28/29/30/32/33（决策 C ④） | 低成本确定性项 |
| 第四批 | V23-38/40/42/44/46/48（决策 C ⑤⑥） | 性能与测试护栏 |
| 第五批 | V23-24/25 | 两个低成本 P2 |
| 丢弃 | V23-07 | 拍板：不接受多等位 |

未落地：V23-26（需要 header-aware 区域解析并同步 CLI seeking）、V23-34…37/39/41
（等 profile 的性能/精简）、V23-47 的剩余 oracle 项（WC84 全缺失实测、hobs/hexp）。
V23-19/20/23/43/45 已随相关条目落地或被覆盖。

## 实现说明（按文件）

### reed-io：BCF 描述符（V23-01）

`bcf/record.rs`

- `typed()` 的 count==15 分支不再调 `typed_int()`，改调新的私有
  `extended_count()`：读**一个**描述符字节，要求其 count 半字节为 1
  （`bcf_enc_size` 的编码恒为 `1<<4 | type`）、类型为 INT8/16/32，然后读值。
  二次 `15` 直接 `MalformedRecord("second-level extension")`。
- `typed_int()` 不再经 `typed()`：自己读描述符并要求 count 半字节为 1。
  这正是 htslib `bcf_dec_typed_int1` 的形态（读 type nibble + 值，不递归）。
- 两个函数都不再有互递归路径，栈深度变为常数。

测试（`bcf::record::tests`）：

| 测试 | 覆盖 |
|---|---|
| `second_level_extension_is_malformed_not_a_stack_overflow` | 2/8/64/4096 个 `0xF1`：两条入口都返回错误 |
| `extended_count_rejects_non_integer_and_non_unit` | 扩展描述符为 CHAR/INT16(非 1 count)/count 6：拒绝 |
| `typed_int_reads_a_single_descriptor_without_extending` | 正常 INT8、count 2、`0xF1` 三种边界 |

### reed-io：CHAR FORMAT（V23-02）

- `bcf/render.rs::fmt_gt` 返回 `Result<(), &'static str>`，首行拒绝非整数
  类型（与 htslib `bcf_format_gt_v2` 的 `default: return -2` 对齐）；
  `render_sample_part` 把它映射成带 locus 的 malformed。
- `bcf/gt.rs` 的两个内部调用点：表构建用 `expect`（INT8 恒成立），
  `decode_slow` 用 `map_err` 保留 locus。
- `bcf/columns.rs::float_element_column` 对非 FLOAT/非整数返回 NaN。
- `bcf/encode_native.rs::encode_gt_column_bcf` 的守卫从 `width == 0` 收紧为
  "必须是整数类型"。

测试：`char_typed_gt_is_malformed_not_an_out_of_bounds_panic`（函数级）、
`render_record_rejects_char_typed_gt`（record 级，断言错误含 `Chr01:20`）、
`char_typed_float_column_element_is_missing_not_a_panic`。

### reed-io：尾部空列（V23-03）

`vcf/line_edit.rs::edit_samples` 增加 `sample_idx < n_samples` 上界，样本数
从 plan 的 bitset 读取（`SampleMask::n_samples` / `Bitset::len`，两者都由
header 构造），因此 `edit_vcf_line` 的公开签名不变、三个调用点不动。

测试 `trailing_empty_columns_do_not_overrun_the_header_sample_count`：
两个尾随空列 × 样本子集、× GT 掩码，外加单尾随 tab 的对照。

### reed-io：负浮点（V23-04）

`format_g6_into` 指数分支在 mantissa 前补 `-`（`wide`/`scientific` 由
`v.abs()` 生成，定点分支本来就格式化 `v` 本身）。

测试 `g6_negative_scientific_values_keep_their_sign`（9 个值，含
`-0.000099999`、`-999999.5`、`-0.0` 与正数对照）与
`negative_info_float_survives_an_info_edit`（未触碰标签逐字节保真；
`info_set` 重算走渲染器）。

实跑：`bcftools view -Ob` 生成的负 INFO 浮点 BCF，`reed filter -Ov` 输出
`XF=-1.5e-07`，与 `bcftools view -Ov` 逐字节相同。

### reed-filter：表达式深度（V23-05）

`expr/parser.rs`

- `MAX_NESTING = 128`：`(` 与 `!` 每层递增，超限返回
  `"expression nesting exceeds 128 levels"`。
- `MAX_CHAIN = 2_000`：单条 `||`/`&&` 链的项数上限（左深树的 Drop 与求值
  同样递归）。
- `Parser::nested()` 负责增减计数；`unary` 的 `!` 与 `primary` 的 `LParen`
  两个递归点都经它。

标定过程（`/tmp/v23probe2`，独立 probe crate，不改仓库）：debug 构建下
`parse + Drop` 在 2 MiB 线程栈上 ~430 层溢栈、1 MiB 上 ~200 层溢栈；2,000
项链在 512 KiB 栈上仍可 parse+drop。拍板文档写的是"如 500"，实测 500 在
2 MiB 栈上不安全，故按实测取 128（约 2× 余量）。

测试：`nesting_is_bounded_at_the_limit`（127/128/129 + 128/129 个 `!` +
10,000 层）、`logic_chains_are_bounded`（2,000/2,001）。

### reed-core：R 平局（V23-06）

`math/r_compat.rs` 的平局判定改为精确 `frac == 0.5`；文件级
`#![allow(clippy::float_cmp)]` 并注明"精确比较是语义要求"。

测试 `values_near_a_half_are_not_ties`：复核报告的三个实测值
（`1.4999999999999998`、`2.5000000000000004`、`0.35*90`）与
0.05…0.95 × 组大小 1…500 的全扫描（对照精确判定）。

### reed-reduce：NaN 排序（V23-08）

`graph.rs::isomap`、`grm_reduce.rs::knn_from_distances`、
`grm_reduce.rs::isomap_from_distances` 的 `partial_cmp(..).unwrap()` 改
`total_cmp`。同文件两处已有 finite 过滤（`auto_gamma` 的 `d > 0.0`、
`umap` 的 `is_finite()`）也顺手改成 `total_cmp`，与仓库风格统一。

复核的同类位置核查结论：`grm_distances` 的 `.max(0.0)` 会把 NaN 映射为 0，
所以 GRM 入口实际不会喂给这两个排序 NaN——但函数级仍应按 NaN 安全对待，
测试据此写成"函数级 fixture"，并另加一条钉住 clamp 行为的
`grm_distances_clamp_nan_to_zero`。

测试：`isomap_sort_tolerates_nan_rows`、
`knn_and_isomap_sorts_tolerate_nan_distances`。

### reed-reduce：特征向量符号（V23-09，拍板 A）

`backend.rs::SymEigResult::top_k` 与 `eig.rs::SymEig::top_k` 对每列做
"最大 |分量| 取正"（严格大于才替换，故并列取行序最先；全零列 sign 为 +1，
保持 0 而非 -0）。

测试：`top_k_sign_is_normalized_across_backends`（规范化断言 + faer 与
fallback 逐值对比，容差 1e-9）、`top_k_leaves_an_all_zero_column_alone`、
`top_k_tie_rule_takes_the_first_row`。

影响：`.eigenvec`/embedding 坐标的符号会变（拍板已接受）。全量测试与
`e2e_pca`、`ordination` 等既有锚点未发现依赖绝对符号的断言。

### reed：LD pairs 流式化（V23-14，拍板 C）

- `stats/request.rs`：`LdRequest` 新增 `pairs_path: Option<PathBuf>`（库调用
  方默认 `None` 时行为不变）；`LdPlan` 同步携带。
- `stats/sink.rs`：`StatsSink::pairs_file` 惰性 `File`；`join_ld`（reducer
  线程、contig 顺序、`finish` 与中途 cap join 共用）对每个作业的 rows 调
  `write_pairs`，随后 `set.pairs.clear()`，report 的 `pairs` 保持为空。
- `stats/writers/ld.rs`：导出 `PAIRS_HEADER`、`render_pair_row`，让流式与
  内存两条路径共用同一份字节定义。
- `reed-cli/translate/popg.rs`：`build_ld_request` 增加 `out_prefix` 参数，
  `--ld-pairs` 且有可用前缀时设置 `<prefix>.ld`。

测试：`streamed_pairs_file_is_byte_identical_to_the_accumulated_report`
（流式 vs 内存渲染逐字节 + 无重复写 + decay 不变）、
`streamed_pairs_file_is_not_created_without_pairs`（单一位点不建文件）。

实跑：`reed popg ld.vcf --ld-decay 1000 --ld-pairs -o out` 产出
`out.ld`（4 行，含表头）与 `out.ld.decay`。

### reed-cli：前缀拒绝（V23-21，拍板 A；含 V23-31）

`translate/common.rs` 新增 `require_output_prefix(prefix, command)`：

- `prefix == "-"` → usage（"writes PREFIX files: -o - is not stdout"）
- `prefix.trim().is_empty()` → usage（隐藏文件）

调用点：`stats`（替换原空检查）、`convert`、`pca`、`distance`（新增）、
`popg`（从 run 之后移到参数校验阶段，删除原来的 `prefix == "-"` 后置检查）。
`filter` 的 stdout 语义不变。

测试：`prefix_commands_refuse_dash_and_empty_output`——5 条命令 ×
`-`/`""`/`"   "` 共 15 组，断言 exit 2、stderr 点名 `-o`，且工作目录不出现
以 `-` 或 `.` 开头的输出；末尾用真实前缀跑一次 `stats --freq` 作对照。

### 文档（V23-50、V23-47）

- `docs/spec/cli.md` §5：新增 `AN = 0` 的 AF/MAF NaN 与 IEEE 比较条目。
- `crates/reed-filter/src/expr/mod.rs::field_value` 文档注释记录同一决定。
- `crates/oracle/src/estimators/spectrum.rs` 模块文档：complete-case 口径、
  与 scikit-allel 的差异及其理由。
- `CHANGELOG.md`：本轮条目。

测试：`af_and_maf_are_nan_at_an_zero_and_compare_by_ieee`（四条比较符 ×
全缺失/已调用两种位点）、oracle
`sfs_theta_and_tajima_are_complete_case_per_site`（3 位点含 2 个部分缺失，
断言 `n_sites=1`、`n_skipped_missing=2`、bin 与 θ_W）、
`pi_is_pairwise_complete_where_the_spectrum_is_complete_case`（同一 fixture
上 π 的 pair 口径与谱系口径不同）。

## 第二批：决策 A/B 与内存三项（2026-11-03）

> 拍板依据见 [STATUS.md](STATUS.md) §"决策点（第二批）"。这一批的代码在
> 拍板后先行完成，本条记录补齐逐项实现、测试与文档。

### reed-filter：缺失值的算符表（V23-59；含 V23-23）

bcftools `filter.c:4449-4455` 的 `missing_logic` 是 `{0,0,1}`/`{0,1,1}`：

| 操作数状态 | `==` | `!=` | `< <= > >=` |
|---|---|---|---|
| 一侧缺失 | 假 | **真** | 假 |
| 两侧都缺失 | **真** | **真** | 假 |

三处入口统一到这张表：

- `expr/mod.rs::Operand::any_pair`：缺失对不再 `continue`，按
  `(a_missing, b_missing)` 查表；真则短路，否则继续。`Eq`/`Ne` 两侧缺失时
  为真，一侧缺失只有 `Ne` 为真。`> <` 等仍为假。
- `eval/head.rs::missing_value_passes`（新）：INFO field 规则的
  `on_missing = pass` 表示"按算符表"，`fail` 是显式一律丢弃；`keep`/`mask`
  在 `plan/mod.rs` 被拒绝，防御性按 `pass` 处理。tag 缺席与 `TAG=.` 同一条
  路径（V23-23）。
- `eval/head.rs` 的 `Number=A` 分支：`.` 元素同样按算符表计入 `any`/`all`。
  这一处是补齐时自查发现的：第一版只写了 `all = false`，与注释和 bcftools
  都不符；现在 `!=` 在 `.` 元素上通过，其余算符失败，`all` 量化器要求每个
  元素各自通过。

测试（`doc_v23 V23-59` 注释）：

| 测试 | 覆盖 |
|---|---|
| `missing_operands_follow_the_bcftools_operator_table` | `DP=.`/`DP` 缺席 × 六算符 + 已调用对照 |
| `both_sides_missing_pass_for_eq_and_ne_only` | 字面量 `.` 与缺失 tag：`==`/`!=` 真，`>`/`<` 假 |
| `missing_info_values_follow_the_bcftools_operator_table` | head 规则：`. `与缺席 × `!=`/`==`/`>`/`<`，`fail` 覆盖 |
| `number_a_missing_elements_follow_the_operator_table` | `AD=3,.`/`AD=5,.` × `!=`/`==`/`>` × `any`/`all` + `fail` |

实测（bcftools 1.24-2，`DP=.`@10、`DP` 缺席@20、`DP=5`@30）：`-i 'DP!=3'`
两边都是 `[10 20 30]`，`-e 'DP!=3'` 两边都是空。

### reed-filter / reed-cli：区域两族（V23-60/61；含 V23-19/20）

- `spec/region.rs`：`RegionRule` 新增 `indexed: Vec<String>`（`-r/-R`），与
  `include`/`exclude`（`-t/-T`，`^` 为排除）分开；`is_empty` 计入三组。
- `plan/mod.rs::compile`：`indexed` 先编译成一条 `HeadRule::Regions`
  （exclude 空），`-t/-T` 再一条；两条 AND，即 bcftools 的跨族交集。
- `cli/filter_flags.rs`：四个 `Option<String>` 改为 `Vec<String>`，允许重复。
- `translate/common.rs`/`translate/filter.rs`：同族累加（并集），跨族分别
  进入 `indexed` 与 `include/exclude`。
- `spec/render.rs::to_args`：`indexed` 渲染回 `-r`，`include/exclude` 仍渲染
  为一条 `-t`（排除项带 `^`），round-trip 保真。

同族重复取并集是**有意分歧**：bcftools 把 `-r`/`-R` 写进同一个槽位，后写
覆盖且 `is_file` 不复位（`-R f -r chr1` 把 `chr1` 当文件名）。顺序相关的
静默丢数据不固化成契约，已登记在 `docs/spec/cli.md` 的分歧表。

测试：`region_flag_combinations_match_bcftools`（9 组：单族、跨族交集、
不相交、`^` 排除、`-T` 文件、重复 `-t` 两种顺序、双 include），外加
`-r`+`-R` 并集。fixture 与实测表（STATUS 决策 B）一致。

### reed：pass0 驻留预算累计（V23-10）

`filter/pass0/matrix.rs` 新增 `ResidentBudget::admit`：`MatrixSink` 记录
`resident_bytes`，每接收一个 contig 矩阵用 `saturating_add` 累计并与预算
比较；拒绝的矩阵不计数，所以 `matrices` 与计数器保持同步。旧检查是
`matrix.size_bytes() > budget`，十个各占 90% 预算的 contig 会逐个通过、
全部驻留。

测试：`admits_by_running_total_not_per_matrix`（6+6 > 10、拒绝后不计数、
恰好填满）、`accepts_an_exact_fill`（`>` 而非 `>=`）、
`overflow_saturates_instead_of_wrapping`；e2e
`sketch_budget_covers_the_sum_of_every_resident_contig`（6 contig × 40 位点，
16 KiB 预算；无预算对照成功）。

### reed：`--reference` 流式读取（V23-11）

`convert/job.rs::read_reference` 从 `fs::read` + 整文件按 `\n` 切分改为
`BufReader::with_capacity(1 MiB)` + `split(b'\n')`。峰值从"文件全量 + 每个
contig 的 `Arc`"降为"当前 contig 的 `Vec` 与 `Arc`（短暂两份）+ 其余
`Arc`"。解析规则不变：`\r` 去掉、空行跳过、`>` 后的描述截到空白。

测试：`msa_reference_parses_crlf_blanks_and_multi_record_fastas`（CRLF、空行、
带描述的 header、两条记录；断言 FASTA 逐字节）。

### reed-distance：GRM 原地内积（V23-12）

`matrices::grm` 删掉 `n × m` 的 `x` 缓冲（2000×1e6 时 16 GB），直接读
`st.w`；内积循环同时累计 `called_both` 作为 `Pairwise` 分母。`w` 在未调用
格为 0，所以"缺一个就加 0"与旧零填充缓冲完全等价；方差标准化用 `vi*vj`，
协方差用 `vi*vj*sd²`。

测试：`reference_grm`（旧算法的逐字重写）在四组
`(MeanImpute|Pairwise) × (variance|full)` 下与 `grm` 逐值对拍。

## 第三批：决策 C ④ 低成本确定性项（进行中）

### reed-stats：bp 窗口的 `pos_log` 截断（V23-13）

`ld.rs::LdPruner` 的 `pos_log` 记录每个 arrival 的位置，供 bp 窗口推进时
按 `bp >= trigger − W` 找新起点。它原来只随 contig 增长、从不缩短；复核
已修正为"单 contig 内"（`finish_contig` 会 `clear()`），但大染色体 × 每
线程仍是 O(位点数) 的内存。

改动：新增 `pos_base`（`pos_log[0]` 的 arrival 下标）。`evaluate_window_bp`
先按 `arrival − pos_base` 索引搜索，`evaluate_window` 把 `window_start`
推进到 `new_start` 后，把 `pos_log` 的前缀 `[0, window_start − pos_base)`
`drain` 掉并同步 `pos_base`。搜索区间本来就从 `window_start + 1` 开始，
被丢弃的前缀不可能再被读取；`finish_contig` 的 `usize::MAX` 路径不截断
（随后的 `clear()` 负责）。

测试（`ld::v23_pos_log_tests`）：

| 测试 | 覆盖 |
|---|---|
| `bp_window_drops_the_dead_prefix` | 10,000 个 1 kb 间隔位点、100 bp 窗口：每次 push 后 `pos_log.len() <= 2`，全部保留 |
| `bp_window_keeps_the_in_span_entries` | 3 kb 窗口：日志只保留 span 内到达项（`<= 5`），结束 contig 后 kept == seen |

### reed-reduce：`tql2` 不收敛不再静默（V23-15）

`eig.rs::tql2` 的 `iter > 100 { break; }` 让未变换的尾块照常返回，PCA 会
静默得到错误特征值。现在 `tql2` 返回 `bool`，`sym_eig` 用
`assert!(converged, "tql2: ... did not converge ...")` 终止。EISPACK 的
100 次上限远高于正常的 2–5 次，命中意味着输入退化（NaN/Inf）或矩阵病态；
继续返回错误结果比报错更糟。

测试：`non_convergence_is_reported` 通过 `tql2_with_limit(.., 0)` 直接触
发上限分支（普通对称矩阵在 100 次内必收敛，无法自然构造失败样本），断言
返回 `false`。

### reed-engine：先取 cell，再判看门狗（V23-16）

`frame_chain.rs::wait_inner` 原来先检查全局静默是否超过 `stall_limit`，
再取已就绪的 cell。慢存储停顿 >120 s 后恢复、且 cell 已经发布时，健康
运行仍被判"chain stalled"。现在循环内先 `take` cell，拿到就返回；看门狗
只对"cell 未就绪"的等待生效。

测试：`a_ready_cell_is_served_despite_a_stale_publish_clock`（1 ms 上限、
发布 cell 1 后 sleep 5 ms，`wait(1)` 仍返回 cell）。

### reed-distance：`standardize` 只把有限值当观测（V23-17）

`standard.rs` 的两处 `is_nan()` 改为 `is_finite()`：Inf/-Inf 既不计入
`sum`/`cnt`/`called`，也不参与标准化（写 0），与 NaN 同一条缺失路径。
此前 Inf 被当作已调用，`mean`/`sd` 变成 Inf/NaN，整条链静默全 NaN。

测试：`infinite_dosages_are_missing_not_observations`（`[0,1,+Inf,-Inf]`：
mean=0.5、sd=√0.375、两个 Inf 不算 called、`w` 有限）。

### reed-reduce：核 PCA 负特征值截零（V23-18）

`graph.rs:120`、`grm_reduce.rs:71/288` 的 `vals[c].abs().sqrt()` 改为
`kernel_scale(vals[c]) = vals[c].max(0.0).sqrt()`。PSD 核上的负特征值来自
数值误差，`abs()` 会给出"虚数模长"的坐标；截零才与"该方向不存在"一致。

测试：`negative_kernel_eigenvalues_scale_to_zero`（-1e-12/-3 → 0 且非
NaN，4 → 2）。

### reed：LD decay 的 `None` 约定（V23-22）

`LdPlan::decay_edges` 原来对 `decay_max_bp = None` 合成 `[0, bin)`，而
sink 的 pair 跨度缺省是 100_000，两处缺省不一致；pairs-only 的库调用会
得到一条"看似完整"的单行 decay。现在 `None` 返回空 edges，`LdDecayAcc`
的 bin 数就是 `edges.len()`（空 = 0），push 全部 no-op，`write_ld` 因此
不写 `.ld.decay`——与 CLI（要求 `--ld-decay`）同一约定。

测试：`no_edges_means_no_decay_bins`。

### reed-stats：decay 的 `bin_of`（V23-27）

`LdDecayAcc::bin_of` 在距离小于首条 edge 时返回 `Err(0).saturating_sub(1)`
= 0，把该 pair 静默并进 bin 0。现在 `bin_of`/`bin_index` 返回
`Option<usize>`（`Err(0) → None`），`push_bin` 收到 `None` 直接跳过。
仓库调用方的 edges 总从 0 起步，行为不变；库 API 不再误报。

测试：`a_distance_before_the_first_edge_has_no_bin`（edges `[100,200]`：
50 → None、150 → 0、250 → 1）。

### reed-stats：merge 不变式（V23-28）

三处 `debug_assert!` 在 release 是空操作，随后 `zip` 会截断较长的一侧：

- `StatAcc::merge`：不同 variant → 无条件 `assert_eq!`，`_` 分支改
  `unreachable!`（原来 `_ => {}` 静默丢数据）。
- `ContigBins::merge_bins`：bin 宽度 + 每 bin 累加器数量都无条件断言。
- `BpBlockJackknife::merge`：block 大小无条件断言。

测试：`merge_bins_rejects_a_different_width`、
`merge_bins_rejects_a_different_accumulator_count`、
`merge_rejects_different_accumulator_kinds`、
`bp_merge_rejects_a_different_block_size`（均 `should_panic`）。

### reed-io：诊断与空 ALT（V23-29、V23-30）

- `IoError::MalformedRecord.reason` 与 `MalformedBgzf.reason` 从
  `&'static str` 改为 `String`（`encode_noodles` 的 `Box::leak` 随之删除）。
  `genotype_scan` 的样本列错误改为 `format!("... (found {found}, expected
  {expected})")`，两种方向（过少/多余）各自带计数。
- `count_alts(b"")` 返回 0（原来 commas+1 = 1）；`.` 与空列同为零。

测试：`column_count_mismatch_is_an_error` 增加 `found 1`/`expected 3`
断言；`empty_alt_has_zero_alleles`（`""`/`.` → 0，`A` → 1，`A,G` → 2，
尾随逗号仍是 2）。

### reed-cli / reed-graph：表类 `--format`、路径与 tree 文档（V23-32）

- `DistanceCli.format` 改为 `Option<String>`（无 clap 默认值）。表类
  （`king`/`pihat`/`phi`/`ajk`/`discordance`）显式给 `--format` 是 usage
  错误；矩阵类缺省仍为 `triangle`。
- `common.rs` 的 `--contig-lengths` 读取错误带文件名（`{path}: {e}`）。
- `tree.rs` 模块文档：内部节点 `n_tips + 1..=2 * n_tips - 1`（最大
  id = 2n−1）。

测试：`format_is_rejected_for_table_distance_kinds`（表类 exit 2 且点名
`--format`，矩阵 `ibs` 对照成功）、`contig_lengths_error_names_the_file`
（`popg`，stderr 含缺失路径）。

### reed-reduce：diffusion 相对阈值与 backend 文档（V23-33）

- `graph.rs` 的平稳模式判定改为 `is_stationary_mode(v, λ_max) =
  |v − 1| <= 1e-9 * λ_max`，不再用绝对 1e-9。
- `pca::backend()` 文档改为"编译期 backend 常量"，并指向
  `sym_eig_top_k` 返回的 per-call backend。

测试：`stationary_tolerance_is_relative_to_the_spectrum`（λ_max=1e-6 时
绝对 1e-9 会误跳过的模式不再被跳过）。

## 第四批：性能与测试护栏（V23-38/40/42/44/46/48）

### reed-reduce：PCA 的 faer GEMM（V23-38）

`pca.rs` 的 mean-impute GRM 从 `Grm::mean_impute`（标量三重循环）改为
`mean_impute_grm`：`faer-backend` 下走 `backend::gram_divide`（此前全仓
零调用），无 faer 构建仍走标量路径。分层原因：`reed-distance` 是
`reed-reduce` 的下游，不能反向调用 `gram_divide`，所以委托放在 PCA 层而
不是 `Grm` 本身；`Grm::mean_impute` 的公开行为不变。

测试：`faer_gram_matches_the_scalar_mean_impute`（6×4 fixture，faer 与
标量逐值 1e-12、对称性）。

### reed-cli：样本来源校验单点化（V23-40）

`common::split_sample_sources` 收拢 `-s` 两种极性、`-S`/`--keep`、
`--remove` 的解析与冲突规则（重复同极性、两个 keep、两个 remove、空
列表）。`filter::spec_from_flags` 用它生成 `SampleSelect`，
`common::parse_sample_keep` 用它解析文件与掩码，删除两份约 40 行的重复
实现。错误文案统一为信息更全的一版（含 "combine the names"）。

### 测试护栏（V23-42/44/46/48）

- **V23-42**：`pruned_header_preserves_order_and_selection`（4 样本
  `[T,F,T,T]` → `s1,s3,s4`）、`pruned_header_treats_a_short_keep_as_kept`
  （`keep=[false]` 的短列表按文档保留尾部）。
- **V23-44**：`every_nucleotide_subset_has_its_iupac_letter`——15 个非空
  {A,C,G,T} 子集逐一断言字母且必须在 62 项表中；每个表键必须是排序去重
  形式。`pretty_sequence` 的 `n≠10` 行为固定为逐整数（R 可给出 0.5 步长，
  `usize` 返回类型无法表达；`xval` 只用 `n=10` 的实测表）。
- **V23-46**：`wc84_three_populations_accept_a_half_called_group`——r=3、
  一组为孤立 `./1`，用 allel 1.3.13 实测的 `n = an // 2` 口径对拍
  （`ref_wc84_floored_n`），池化 Fst 为 `-0.70271714922049`；
  `patterson_accumulators_merge_split_partials`——F2/F3/F4 按 2+2 拆分
  后 merge，值与单遍完全一致（含 F4 的分子/分母和）。
- **V23-48**：`faer_is_stable_across_one_and_four_threads`——10×10 随机
  对称阵，faer 在 1/4 线程池与 fallback 的特征值/特征向量逐值 1e-9
  一致（补充既有 `top_k_sign_is_normalized_across_backends`）。

## 第五批：两个低成本 P2（V23-24、V23-25）

### reed-filter：`to_args` 的两字符算符（V23-24）

`spec/render.rs` 的掩码渲染先用 `starts_with('<')` 选 `--min-dp`，但取值用
`trim_start_matches(['<', '>'])`，`<= 10` 只剩 `"= 10"`，provenance 命令行
不可重放。现在按 `<=`/`>=`/`<`/`>` 的顺序整段剥离算符再 trim。

测试：`two_character_mask_bounds_render_replayable_args`（`<= 10` →
`--min-dp 10`，`>= 20` → `--max-dp 20`）。

### reed-filter：空 `SampleSelect`（V23-25）

`plan::compile` 用 `spec.sample_select.is_some()` 判定 cohort 变化，JSON 的
`"sample_select": {}` 会让 `-q/-Q` 从 INFO AC/AN 短路静默切到基因型重算。
现在用已有 `SampleSelect::is_empty()`：只有真正选择样本才算变化。

测试：`empty_sample_select_keeps_the_info_af_short_circuit`（空选择保留
`HeadRule::MinAfInfo`，真实 `keep` 移除它）。

### 仍待拍板

- **V23-26（HLA contig 的 `:`）**：需要把 header contig 表带进
  `Region::parse`（或新增 header-aware 解析），并同步 CLI 的索引 seeking
  路径；不在本轮拍板范围。
- V23-34…37、39、41（profile 先行的性能/精简项）与 V23-47 的剩余 oracle
  项（WC84 全缺失实测、hobs/hexp 口径）继续挂起。

## 验收

| 检查 | 结果 |
|---|---|
| `ci/test-parallel.sh`（全量工作区） | exit 0；**931 passed / 0 failed**（56 二进制，15 s 墙钟；v22 基线 867，本轮 +64） |
| `cargo clippy --workspace --all-targets -- -D warnings` | 0 告警 |
| `ci/check-doc-drift.sh` | exit 0；97 项全 ok |
| 真实 `subvcf.vcf.gz` 的 `--summary` | DP/PSC 与 `bcftools stats -s -` 逐行相同 |
| 真实 BCF（`bcftools view -Ob`，负 INFO 浮点） | `reed filter -Ov` 数据行 == `bcftools view -Ov`（`XF=-1.5e-07`） |
| 尾部空列 + `-s` | 修复前 panic（复核记录），修复后 exit 0、输出单样本列 |
| 损坏 BCF（`0xF1` 填充 indiv） | 修复前 SIGABRT；修复后 exit 1 + `MalformedRecord` 诊断 |
| `popg --ld-pairs` 流式 | `out.ld` 4 行 + `out.ld.decay`；无 pair 输入不建 `.ld` |
| 五命令 `-o -`/空前缀 | 全部 exit 2，前置报错，无 `-`/`.` 开头文件 |

## 自查修正（实施后）

对最高风险项做了一次只读自查，抓到并修掉一个真实缺陷：

- **V23-03 的第一版修复仍会写回尾部空列**：`edit_samples` 用 `sample_idx <
  n_samples` 提前结束循环，但循环末尾的"整段逐字节拷贝"
  （`samples[start..]`）把超出 header 的空白列一并写出。第一版测试用
  `trim_end()` 断言，恰好把多出来的尾部 tab 掩盖了。现在先把样本区截到
  `sample_region_end(samples, n_samples)`，并改用**字节精确**断言（含
  "无样本编辑时仍逐字节保真"的对照）。

其余自查结论（无需改动）：`nested()` 用闭包包住子解析，内层任何早退都会
递减深度计数；`write_pairs` 对空 rows 直接返回，所以"先出现无 pair 的
contig、后出现有 pair 的 contig"仍按需创建（新增
`streamed_pairs_skip_a_pairless_leading_contig` 钉住）；`top_k` 的 sign
不影响特征值；`popg -o -` 现在 4 ms 返回（此前要跑完整轮）。

第二批自查又抓到一处：

- **V23-59 的 `Number=A` 分支只改了注释**：补文档时发现 `.` 元素仍走
  `all = false`，`!=` 在该元素上不通过，与同一函数的标量分支和 bcftools
  都相反。改为 `missing_value_passes` 计入 `any`/`all`，并补
  `number_a_missing_elements_follow_the_operator_table`。这是"边写文档边
  核对代码"直接暴露的缺陷。

## 未做的事（刻意）

- **V23-07（AC 多等位）**：拍板丢弃，未改代码。
- **B 组剩余项**（V23-13 `pos_log` 截断、V23-15 tql2、V23-16 frame_chain、
  V23-17 `standardize`、V23-18 kernel PCA 负特征值、V23-22 `decay_max_bp`
  约定）：仍待拍板。
- **性能/精简（V23-34…41）与其余测试缺口（V23-42…46、48、49）**：未动。
- 本轮未重写历史、未改 `Cargo.toml` 版本、未动 gold 基线。