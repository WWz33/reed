# doc_v15 实施状态

> 本页记录 F1–F5 的落地状态与证据。阶段二核实基线（40 条判决、复现
> 证据、产物路径）见 [01-review.md](01-review.md) §0/§7。每条完成以
> "修复前的复现命令重跑不再复现 + 回归测试绿"为证据。

## 实施状态（F1–F5）

| 项 | 内容 | 状态 | 证据 |
|---|---|---|---|
| F1 | 机械修复 10 条（D15-1/5/21/22/27/30/33/13/14/39） | **完成** | `e2e_filter`（Float 掩码回归）、`cli_contract::v15_usage_gates`；复现夹具全部转不复现：Float GQ 掩码、三态跳过、`--levels ""` exit 2、`--step`+`--bed-file` clap 冲突、debug 无 panic（POS=0 计数跳过）、`--idspace-to` 互斥、`--pbs A,B,A` 点名、无 `-o` 0.006s 快失败、EPIPE 由 panic(101)→exit 1 "Broken pipe"、`-o -` 走 stdout 且不再生成字面文件 `-`（filter/rename/convert/popg 裸路由） |
| F2 | 统计口径 7 条（D15-2/11/12/31/18/9/10） | **完成** | D15-2：`precac.vcf --min-dp 100 -q 0.2` 保留→丢弃（与无 INFO 重算路径一致；`info_trusted: !cohort_changed` 单一来源，驱动只留 subset 运行期清空）+ `e2e_filter::min_af_recomputes_after_a_mask_changes_the_cohort`；D15-11：Hudson/PBS 半调用窗口 NA→0.0，`fst_multi::hudson_zero_numerator_is_zero_not_na`；D15-12/31：oracle 三组门移除、n_sites 全计、f2/f3/f4 移植偏差校正式与 `Σ` 口径（`families.rs` 手算值经 scikit-allel 1.3.13 复核：f2=0.0、f3=0.5、f4=−0.25、dstat=−1.0），`contracts.rs` PBS 审计改为 (1,1)；D15-18：`--pi` 与 `--pi --f2` 输出逐字节一致、全变异文件 + `--pi --f2` 报检查错（固化测试 `cli_contract::f2_does_not_switch_the_pi_convention`）；D15-9：检查已下沉共享的 `build_sfs_plan(sfs, group_names, group_sizes)`，**popg 与 convert 同一实现**；`popg|convert --sfs --proj 40` 均报 `exceeds population A's 4 chromosomes`，`--proj 4` 对照成功（见「审计回补」）；D15-10：重复分配报 `sample "s1" is assigned to both "A" and "B"`；`cargo test -p reed-stats`/`-p reed-oracle` 全绿 |
| F3 | 格式与窗口 11 条（D15-3/4/6/8/15/19/20/24/25/26/29） | **完成** | D15-3：`boundary.bcf`（H=65280）4/4 条与 bcftools 一致，H=65279/65281/130560 对照不回归；D15-4：`--phylip --012` 的 `.012.pos` 与单独 `--012` 逐字节相同；D15-6：实现取"translates 侧补 formats-empty 的 Phylo 目标"（不改库门面），bin-only 行序 s3 首、未知名报错、无 outgroup 路径不变（job.rs 空格式校验放宽为 `options.outgroup.is_none()` 时才拒绝）；D15-8：**每三元组一份文件**（单三元组 `PREFIX_pbs.txt` 逐字节不变——sha256 与改前相同；多三元组 `PREFIX_pbs.<A>-<B>-<C>.txt` 各持自身表头，`--pbs A,B,C --pbs C,D,A` 的两份文件表头/行一致；依据 = cli.md 自身“可重复”承诺 + §5 每三元组列布局，命名规则已明文写入 cli.md §7；回归测试 `e2e_stats::multiple_pbs_triples_get_their_own_files`）；D15-15：裸 popg 全变异文件报检查错、`--variant-only` 可跑、含不变位点文件不变；D15-19：**延迟发射一位**（窗口需出现严格更后位点才 emit），跨块同 POS 重复落入仍开窗口；seam 夹具输出与同块对照**逐行完全一致**（窗口 [1,100] 两份拷贝、窗口和=contig）；依据 = VCFv4.3 Fixed fields 原文“It is permitted to have multiple records with the same POS”；守卫保留为"不得低于已发射水位线"的不变量检查；D15-20：尾窗输出、diffs/comps 窗口和 = contig 行；D15-24：idx.bcf CHROM 与 bcftools 一致（`entry_idx` 单一实现，tag/contig 共用）；D15-25：撞名报错、无输出；D15-26：MSA 300bp REF 输出宽 1299（原 1043）；D15-29：`cat *.hmp.txt.gz | reed convert -` 成功（2 records） |
| F4 | 规模 4 条（D15-17/28/16/7） | **基本完成**（D15-7 降级登记） | D15-17：`SiteKeepSet::contains` 改 partition_point + 二分（`sort_index` 在两条装配路径调用）；1 万 contig、30 万位点单线程 **8.708s → 0.800s（≈11×）**，保留位点数不变（68,706）；D15-28：**两阶段落地**（`rename/dedup.rs`：串行定序 + 按字节封口的批次由 sink 管线池并行 deflate，见 02 §D15-28）。worker 去掉冗余 block push，并把索引键随 partial 带走（sink 不再重复解析）；无 `-W` 时 sink 不解析键；`dedup_line` 改 `Cow`（`.` ID / 未重复 ID 常用路径零分配）；3M 记录 `-Oz`：**1.394s → 0.23s**（对照 `--idspace-to _` 0.091s，15.4× → **2.5×**），`-Oz -W` **0.820 → 0.515s**（对照 -W 0.299s）；decompressed 输出与改前逐字节一致；块布局只由累计字节决定，压缩流与 sidecar 对任意 `--threads` **逐字节相同**；D15-16：`track_chunks` 分支补同款预算检查、`sketch_size_bytes` 补 mask 行（`(4+R)`→`(5+R)`）；三种模式在 `--memory 1` 下报错文案与行为不回归；D15-7：**有界形态落地（三目标）**——hapmap rows、`--012` cells+pos、phylo cols/bin_cols/used 均由 sink 在合并前接管，两级落盘：单 blob 超 `SPILL_THRESHOLD`=4 MiB 即落临时文件，sink 缓冲总量超 `SPILL_TOTAL_BUDGET`=8 MiB 时按固定序强制落盘（覆盖"多点位×多样本"下小 blob 不越单阈值的形态）；渲染时流回，`RowsBlob` Drop 清理（成功/失败/中止三路径）。实测 800k×50 夹具（1 线程，容器含引擎开销）：hapmap **123.4 → 15.6 MB**、`--012` **70.8 → 26.0 MB**、`--phylip` ≤28.1 MB；输出与改前**逐字节一致**（012 三件套 + hapmap `cmp` 通过）、无 spill 残留。测试：4 个单元（append/visit/落盘清理/writer 流回）+ 2 个跨阈值 e2e（hapmap；`--012`+phylip 合并夹具）。**MSA 已并入（同一轮追加）**：MSA 站点按 contig 以帧流（`pos0 u64 LE ‖ width u32 LE ‖ vals`）落盘，渲染改为**逐样本一趟帧遍历**（区域尺寸级的缓冲不再存在）；800k×50 夹具 1 线程峰值 **151.4 → 15.8 MB**（运行期可见 spill 文件、结束清理），输出与旧版**逐字节一致**；测试 `msa_frames_spill_past_the_threshold`（60k 位点×80 样本跨 4 MiB 阈值，逐字节公式对照）。渲染的每样本趟**按样本 seek 只读该样本切片**（不做整帧重读，延续 v13 D13-2 的禁止重复 IO；守卫 `ci/check-doc-drift.sh`）。**"输出早现"决策：保持现状**——输出仍在成功后写；cli.md 无承诺，运行期直写会打破"失败不留半成品"；如需早现按 `.partial`+rename 另开（记于 02 §D15-7 与 README"本轮不做"）|
| F5 | 收口 | **部分** | 回归测试已落：`e2e_filter`（D15-1/D15-2）、`fst_multi`（D15-11）、`cli_contract::v15_usage_gates`（D15-21/22/25/27/30/33）、oracle `families.rs`/`contracts.rs`（D15-12/31，修复前会红）；`cargo test --workspace` 全绿（56 个测试二进制 / 571 项，2026-10-19 冻结态实测）；`cargo clippy --workspace --all-targets` 零告警；cli.md §7 `--pbs` 行已改；**已补**：D15-36 死代码删除（`JointSfsAcc` 及其测试文件、`PattersonKind::Naive` 及其测试——`f2_site`/`f3_site`/`F2Acc`/`F3Acc` 随之简化为校正式单一形态、`Cursor::skip`、`SiteView::info_int/info_float`、`HeaderScan::n_samples`+`count_samples`、`WindowGrid::bin_of_pos0`（含集成探针）、`keep_mode` 死参、`header_level` 弃值、`sink.rs` 死 `or(take)`、`convert` 空臂）、D15-37（活的 `summary::report::SummaryAcc` 增 `merge`+`PartialEq/Eq`，`QcReportAcc::merge` 改为委托，死孪生结构与其单测删除，`summary::SummaryAcc` 再导出活类型，merge-law 测试改指活类型并改为字段生成）、D15-38（采用“不假装跟踪”变体：`PairAcc::add_partial` 用于 fst/PBS 审计，`missing` 只在有 cohort 总数的 Pi 路径维护，字段文档写明）、D15-32（`notify_abort` 持锁后 notify，丢失唤醒窗口关闭）、D15-35 第一处（`run.rs` 每站点 contig clone 去除）；另补固化测试：D15-3（`boundary.bcf`/`control.bcf` 夹具入 `crates/reed-io/tests/fixtures/bcf/` + `e2e_bcf::bcf_header_at_a_block_boundary_decodes`）、D15-9 的 convert 面（`cli_contract::v15_usage_gates` 追加）。**已补（2026-10-19 复核批次）**：D15-35 第二处（`push_sfs` 权重缓冲入 `StatsScratch.sfs_weights`）；回归测试固化——D15-5（`count_matrix...` 夹具加"无 allele-2 使用的三态行"）、D15-10/D15-13/D15-15（`cli_contract::v15_usage_gates` 与 `popg_validates_before_running`）、D15-14（`popg_broken_pipe_is_not_a_panic`）、D15-24（`idx.bcf`/`idx_control.bcf` 夹具 + `e2e_bcf::bcf_contig_idx_wins_over_line_order`）、D15-26（`msa_long_deletion_is_not_truncated`）、D15-29（`gzipped_hapmap_on_stdin_decodes`）。**已补（2026-10-19 收尾）**：D15-16 单测（`sketch_model_includes_the_mask_plane` + `choose_for_mind_switches_at_the_estimate`）；D15-17 bench 场景 `bench/suites/e2e/ld_prune_contigs.sh`（10 万 contig 收尾实测 10k=0.816s）。**收尾后状态**：40 条无未落地的必须修/应当修项（D15-16 为预防性、已按预防性形态落地）。唯一仍在跟进的是**性能类**：D15-17 属 bench 范畴（场景脚本 `bench/suites/e2e/ld_prune_contigs.sh` 已加入，10 万 contig 收尾实测 10k=0.816s）；D15-28 的 `-W` 串行残留与大 ID 集分配见「观察」。D15-19 的固化测试：`stats::sink::tests::same_pos_seam_duplicate_lands_in_its_window` 直接喂同 POS 双段并断言窗口含两份拷贝、窗口和=contig 行；seam 大夹具亦从"报错"翻为"与同块逐行一致"；仓内先例注"真接缝回归需 >1 MiB 夹具"见 `e2e_stats::malformed_record_under_site_filter_is_fatal` 的 shortcut 注。**已补（2026-10-19 收尾二批）**：D15-40（reader 线程纳入 `catch_unwind`，新增 `EngineError::ReaderPanicked` + 测试 `reed-engine/tests/pipeline.rs::reader_panic_becomes_reader_panicked`）；D15-39（`cli_contract::dash_output_streams_to_stdout`）；D15-23（popg 去掉 `--tmp-dir`，stats 帮助改为诚实表述，cli.md §2/D-3 同步）；D15-35 第四处（`format_g6_into` 栈缓冲零分配：3M 次调用 **242.9 → 182.9 ns/call**，累加值不变）；D15-28 两阶段（见 F4 行）。**gold 脚本（本机全可跑）**：`run_convert_parity` PASS、`run_split_parity` PASS（22 cases 逐字节）、`run_hist` PASS、`run_allel_fstats` PASS（最大相对差 0；该脚本入口命令仍是 v14 拆分前的 `reed stats`，本批坐正为 `reed popg`）。**已固化**：D15-20 窗口算术（`e2e_stats::trailing_skipped_site_advances_the_windows`：尾窗存在性 + 窗口和=contig + 尾窗 594 / contig 2994 的分母断言） |

## 与 02 修复规范的偏差（实现选择）

1. **D15-8（最终）**：按 02 的默认方案 (a) 落地——单三元组保持 `PREFIX_pbs.txt`（逐字节不变），多三元组每份 `PREFIX_pbs.<A>-<B>-<C>.txt` 各持自身表头；依据 = cli.md 自身的“可重复”承诺与 §5 的每三元组列布局，命名规则已明文写入 cli.md §7（先前的“拒绝”实现已按 02 默认回正）。
2. **D15-6**：规范写"outgroup 提升到 ConvertRequest 级"，实现改为 **translate 侧在有 `--outgroup` 时也推一个 formats-empty 的 `ConvertTarget::Phylo`**（job.rs 空格式校验放宽为携带 outgroup 时允许）。理由：零库门面改动，行为等价（bin-only 的行序与名字校验都生效）。
3. **D15-14**：`print!` 之外的 `eprintln!`（timing/metrics 路径）未接线；stdout 主路径已修复。
4. **D15-28（补充）**：实现引入 sink 自持的压缩池，触及 `docs/spec/architecture.md` §7.3
   的"进程里不允许出现第二个线程池"【规】——已按 doc_v8 D8-1 的先例在该节点补一条**预算内**
   例外（池大小 = `auxiliary()`、随 sink 生命周期、在途窗口受内存预算约束；worker 在该路径
   先背压空闲，池占其空位）。此为该规格的**明文修订**，实现与实测见 02 §D15-28。
5. **D15-19（最终）**：改为延迟发射（窗口需严格更后位点才 emit），跨块同 POS 重复落入仍开窗口；seam 夹具与同块对照逐行一致。依据 = VCFv4.3 Fixed fields 原文（“Positions are sorted numerically, in increasing order … It is permitted to have multiple records with the same POS”）。

## 审计回补（doc 偏移核查，2026-10-18 追加）

对 02 修复规范逐条做了"与项目自身权威语义对照"的核查，发现并修复 3 处偏移/缺口：

1. **D15-9 检查放错层（真实缺口）**：规范写明 convert 共享 `build_sfs_plan`、同一缺口，但初版把检查放在 `run_popgen_opened`——`reed convert --sfs --proj 40` 仍静默全零（实测）。已把校验移入共享 `build_sfs_plan`（签名改为 `group_names + group_sizes`，同时做数量/下限/上限三查），driver 与 convert 两处调用点同步；实测 convert 报错、`--proj 4` 对照成功。
2. **`to_args` 的 all-false `Levels` 不可拼写**：D15-21 后 CLI 拒绝空 `--levels`，而 `to_args` 对全 false 状态会渲染 `--levels ""`，往返断裂。已让 `to_args` 对该状态返回用法错（与 Msa-无-reference 同款"不可拼写"处理）；`api_surface` 往返测试全绿。与 v11 N1（Default = window）不冲突。
3. **popg 窗口路由的 `-o -`**：裸路由把 `-` 映射 stdout 后，窗口路由仍把 `-` 当 PREFIX（会生成 `-_pi.txt` 之类）。已显式拒绝（`-o - is not stdout`），行为一致化。

核查后**无偏移**的项（对照结论）：

- **D15-5**：v13 §3.8 明文"位点门：双等位（n_alleles ≤ 2），对齐 vft-conv 的 `>2` 判定"——`raw_alt_count > 1` 正是"总等位 >2 即跳"，与自家规范一致；`ALT='.'` 单等位记录保留的语义也保持。
- **D15-19**：收紧后的守卫在"窗口已输出、位点无法表示"时才报错（等位低于水位线）；同一 chunk 内重复、或重复位置尚未越过水位线（窗口仍开着）的合法形态不受影响——与 sink 注释"不得回卷进已输出窗口"的意图一致。
- **D15-24**：`BcfDict::contig()` 的全部读侧消费者（render、native_encode 的 `contig()`）都经该入口；`native_encode` 的写入侧 map 是独立的行序索引，不受影响。
- **D15-17**：`SiteKeepSet` 仅两处构造（prune.rs 单条、matrix.rs 装配点），均已 `sort_index`；`contains` 对重复 contig 条目用 partition_point+take_while，仍正确。
- **D15-2**：`info_trusted = !cohort_changed` 与 53c8aeb 的"仅真正读 INFO 才告警"语义一致（cohort 变化时不再读 INFO，也不该告警）。

曾列为 **owner 决策**的两处，2026-10-19 按"须有明文规范依据、不作凭空假设"的原则**已定案**（实现与依据见上表 F3 行与「偏差」1/4）：

- **D15-8**：**每三元组一份文件**（单三元组文件名逐字节不变）；依据 = cli.md 自身的"可重复"承诺 + §5 每三元组列布局，命名规则已明文写入 cli.md §7。
- **D15-19**：**窗口延迟发射**（窗口需出现严格更后位点才 emit）；依据 = VCFv4.3 Fixed fields 原文"It is permitted to have multiple records with the same POS"。

## 观察（超出本轮 40 条，登记跟进）

- **D15-20 残留（复核排除）**：曾记录"末条为被跳过的三态位点时分母比窗口和少 1 个位置"，2026-10-19 逐项复核为**误读**。跳过位点按 D7-11 既"在位"（计入 `n_present_file`）又"不计比较"（不入 `comparable`），该处理在窗口与 contig 两级一致：尾窗含真实变异的夹具实测分母 3×600 + 600 + 594 = **2994 = contig 分母**（594 = 98 缺位×6 + 该窗 1 个真实位点×6；早前误把每个 100bp 窗都按 600 计）。非缺陷，无需修复；保留记录以防再误读。
- **D15-28 残留（两阶段落地后）**：`-W` 运行中键传递与 CSI 记账仍在 sink 串行（后者 `--idspace-to` 对照同样支付）；大 ID 集下 `HashSet<Vec<u8>>` 的每唯一 ID 堆分配未优化（1M 全唯一 ID 夹具 0.371s vs 对照 0.065s）。

## 修复后的验证门（每批，实测）

| 门 | 结果 |
|---|---|
| 复现夹具重跑 | F1–F4 全部转不复现（证据见上表；夹具仍留在 `/tmp/reed-verify`、`/tmp/verify-bcf`、`/tmp/verify-oracle`、`/tmp/verify-misc`） |
| `cargo test --workspace` | 全绿（56 个测试二进制 / 571 项；含 oracle_parity、accumulator_laws、e2e_* 与 372s 的 `bcf_vs_noodles` 差分套件） |
| `cargo clippy --workspace --all-targets` | 零告警 |
| 输出字节对照 | `--012` 组合修复后与单独一致；`--rm-dup-ids` decompressed 输出一致（仅 `##reed_command` 路径不同），**压缩流与 sidecar 跨 `--threads` 逐字节相同**（`e2e_rename::dedup_batches_keep_order_index_and_bytes`）；LD 剪枝保留数不变；`-q/-Q` 修复改变了"掩码+陈旧 INFO"的预期输出（即修复目标本身） |

## 变更文件索引（本轮）

| crate | 文件 |
|---|---|
| reed | `filter/record_pipeline.rs`、`filter/mod.rs`、`filter/pass0/prune.rs`、`filter/pass0/matrix.rs`、`stats/request.rs`、`stats/run.rs`、`stats/sink.rs`、`stats/driver.rs`、`convert/job.rs`、`convert/hapmap.rs`、`convert/writers.rs`、`rename/mod.rs`、`rename/dedup.rs` |
| reed-filter | `plan/mod.rs`、`plan/strategy.rs` |
| reed-io | `chunk/bgzf.rs`、`bcf/dict.rs`、`vcf/line_edit.rs`（D15-35#4 零分配 G6）；`write/{block_writer,ordered_output}.rs` 未改（复用） |
| reed-stats | `fst.rs`、`pbs.rs` |
| oracle | `estimators/pbs.rs`、`estimators/fstats.rs` |
| reed-cli | `units.rs`、`cli/mod.rs`、`translate/{popg,convert,rename,filter}.rs`、`tests/cli_contract.rs` |
| reed-engine | `error.rs`（`ReaderPanicked`）、`pipeline.rs`（reader `catch_unwind`）、`tests/pipeline.rs` |
| 测试 | `crates/reed/tests/e2e_filter.rs`、`crates/reed-stats/tests/fst_multi.rs`、`crates/oracle/tests/{families,contracts}.rs` |
| 文档 | `docs/spec/cli.md`（--pbs、--tmp-dir）、`docs/README.md`、`docs/rounds/v15/*`、`CHANGELOG.md` |
| gold | `gold/popgen/run_allel_fstats.sh`（入口命令 stats→popg，v14 拆分漏改） |
| F5 批次 | `reed-stats`：`sfs.rs`、`lib.rs`、`patterson.rs`、`pi.rs`、`summary/{mod,report}.rs`、`tests/{patterson_fstats,summary_merge_laws}.rs`（`tests/sfs_joint.rs` 删除）；`reed-io`：`bcf/record.rs`、`vcf/site.rs`、`chunk/bgzf.rs`、`tests/fixtures/bcf/{boundary,control,idx,idx_control}.bcf`；`reed-core`：`grid.rs`、`tests/grid.rs`；`reed-engine`：`budget.rs`；`reed`：`tests/e2e_bcf.rs`（另见上表 F5 行）|