# v27 状态（2026-11-08 立项；2026-11-09 实施第三至十批）

> 本页跟踪 v27 文档项。状态只在代码改动与验收证据齐备后更新；审查发现本身不计「完成」。

## 总览

| 优先级 | 条数 | 已修 | 待修 |
|---|---:|---:|---:|
| P0 | 3 | 3 | 0 |
| P1 | 17 | 17 | 0 |
| P2 | 19 | 19 | 0 |
| **合计** | **39** | **39** | **0** |

## P0

| ID | 初审状态 | 实施状态 | 证据 / 验收 |
|---|---|---|---|
| V27-01 裸路由 `-o` 覆盖输入 | confirmed | **已修**（`translate/popg.rs` + `cli/mod.rs`） | 修复前：exit 0，输入 md5 `2c96ad7b…` → `60815d1e…`，`--verbosity 0` 无任何输出。修法：裸路由（`writes_files == false`）在 `prelude` 前调用 `refuse_output_eq_input`；`-o` 帮助文本改为「窗口统计时是 PREFIX，否则是单个 TSV 的路径」。验收：`cli_contract::popg_bare_route_refuses_output_naming_the_input` —— `-o <输入>` → exit 2、错误串 `would overwrite the input`、输入 md5 不变；窗口路由 `-o <输入>` 仍 exit 0 且只写 `in.vcf_pi.txt`。 |
| V27-02 失败运行发布截断 `.ld` | confirmed | **已修**（`sink.rs` + `driver.rs`） | 修复前：截断 BCF（合法 BGZF）→ exit 1，`o.ld` 从 4498501 行/`bbc8e12a` 变成 4489507 行/`089ee2a6`（覆盖上一轮好文件）。修法：`SinkOutput` 新增 `pending_commit: Option<AtomicFile>`；`sink.finish()` 只 flush 不 commit，`driver` 在 `job.finish(...)` 完整性校验**之后**才 `commit()`；提前返回时 `AtomicFile` 的 `Drop` 删除临时文件。验收：`e2e_ld::failed_run_leaves_the_previous_ld_untouched`（800 位点跨多个 BGZF 块、裁掉末 8 字节重压；**修复前实测失败**）。手工：截断运行 exit 1 后旧 `.ld` md5 与行数完全不变、无临时残留。 |
| V27-03 Hudson/PBS 分母混入 n=1 位点 | confirmed | **已修**（`fst.rs` + `pbs.rs`） | `mpd_ab` 改为两组均 `n>=2` 才进 num/den，跟 pixy。`fst_multi::hudson_two_sites_match_pixy_one_third` 实测 1/3；单点低覆盖改为 `None`。CLI 低覆盖夹具与 pixy 同值 `0.33333333333333337`；`oracle_parity` 用独立 pixy 参照，WC84 与整数审计仍对拍旧 oracle。 |

## P1

| ID | 初审状态 | 实施状态 | 证据 / 验收 |
|---|---|---|---|
| V27-04 LD 单元素尾块漏配对 | confirmed | **已修**（`reed-stats/src/matrix.rs`） | `decay_range` 只在 `start == end` 时早退；单元素范围继续与前序位点配对。验收：`reed-stats::one_site_decay_range_keeps_pairs_to_the_preceding_sites`（0 → 4096 对）、`stats::run::ld_block_boundaries_keep_all_pairs_including_a_single_site_tail`（4095/4096/4097/4098/8193）、CLI `ld_pairs_keep_the_single_site_tail_block`。真实 CLI 4097 位点输出 8,390,656 对 = 4097×4096/2，`.ld` 行数 == decay `n_pairs` 总和。反向验证：恢复 `< 2` → 单元素测试 0 vs 4096。 |
| V27-05 空结果保留旧输出 | confirmed | **已修**（`driver.rs` + `sink.rs` + `writers/ld.rs`） | 请求决定输出集合，不再由非空内容决定：SFS/fstats 在 sink 构造时初始化零值累加器；LD 按请求生成全部 0/NA bins；`--ld-pairs` 在运行开始创建临时表并写表头，完整性检查通过后才 commit。验收：`empty_popg_outputs_replace_the_previous_run_and_are_listed`；库级 `an_empty_sfs_run_returns_and_writes_a_zero_spectrum`、`an_empty_ld_run_publishes_header_only_pairs_and_zero_decay`。反向验证：恢复 `sfs: None` 与延迟建 `.ld` 均失败。 |
| V27-06 `step > size` 误报乱序 | confirmed | **已修**（`sink.rs` seam 校验） | 输入顺序只与 `cur.last_pos0` 比较，`next_start` 仅作窗口 frontier。修复前真实 1.3 MiB BGZF 夹具（40000 位点、--threads 1）报 `33224 followed by 33225`；修复后输出首窗 `no_sites=100`。验收：`sparse_windows_accept_sorted_sites_across_a_real_chunk_seam` + 3 个 sink 单测（空隙合法、同 POS 跨 seam 不漏、真乱序仍拒绝）。反向验证：恢复 `.max(cur.next_start)` 后失败。 |
| V27-07 全基因组统计仍要 contig 长度 | confirmed | **已修**（`driver.rs` + `sink.rs`） | 无窗口级 kind 时使用 `GridMode::None`，`close_contig` 不再读取长度，span 取观测范围；`--bed-file` 仍要求窗口统计。验收：`genome_wide_popg_statistics_do_not_require_contig_lengths`（SFS/f2/LD 三类）及库级无长度测试。反向验证：强制旧 Bp 网格后报 `contig chr1 has no ##contig length`。 |
| V27-08 SFS 组名路径穿越 | confirmed | **已修**（`writers/sfs.rs`） | 修复前：组名 `../../ESCAPED` → `sub/ESCAPED-2.sfs`，四级上跳直达 `/tmp/ESCAPED-6.sfs`；popg 与 `convert --sfs` 同缺陷。修法：新增 `stem()`（`[A-Za-z0-9_.-]` 白名单，首尾点与 `..` 串映射为 `_`）与 `unique_stems()`（同 stem 追加 `-2`/`-3`），四个 writer 统一走它。验收：`e2e_sfs::sfs_group_names_cannot_escape_the_output_directory` —— 文件全在 `-o` 树内、路径不含 `..` 组件；**用未消毒的 `stem` 反向验证**会写到 `/tmp/ESCAPED-6.sfs`。 |
| V27-09 SFS 一维/联合谱撞名 | confirmed | **已修**（`writers/sfs.rs`） | 修复前：组 `A`+组 `4`、`--proj 4` → 一维 `A-4.sfs` 被联合谱覆盖、stderr 把同一路径列两次。修法：`claim(path, seen)` 对本次运行的每个输出路径去重，重复即报 `EngineError::Job`（不静默挑一个赢）。验收：`e2e_sfs::sfs_one_dimensional_and_joint_names_cannot_collide` 断言 `two SFS outputs would write …/A-4.sfs`。正常组名输出与修复前**逐字节相同**。 |
| V27-10 variant-only 空窗口写 NA | confirmed | **已修** | 窗口 101-200、201-300 的 `avg_pi` 与 `avg_watterson_theta` 全 NA；`pi.rs:165` 实际返回 `Some(0.0)`。 |
| V27-11 genome `bp` 随 `--step` 变 | confirmed | **已修**（`WindowRow.denom_bp` + sink） | `bp` 保持公开跨度；variant-only 另用 `denom_bp`。contig/genome 取声明长度，无长度回退观测 span；BED 取裁剪后的区间并集。验收：`contig_and_genome_variant_only_denominators_do_not_follow_step`、`bed_contig_and_genome_denominators_use_the_clipped_region_union`，step 100 vs 1M 数值相同。 |
| V27-12 variant-only `\|S\|` 定义 | confirmed | **已修**（`pi.rs` occupied 计数 + 删除 skipped-present） | 全缺失记录不再占 `|S|`；多等位过滤后也不再占。尾窗 3/594 改为 3/600，contig 13/2994 改为 13/3000。全缺失审计列仍保留；与 pixy 实测 `no_sites`、π、count_missing 全一致。 |
| V27-13 π/θ 的 `no_sites` 矛盾 | confirmed | **已修**（`pi.rs`/`dxy.rs`；θ 已同口径） | π 按该群体有调用计数，dxy 按两组交集计数。CLI `popg_matches_pixy_on_pairwise_occupancy_and_missing_audits` 与 pixy 实测 `A no_sites=1, B=2, dxy=1`，`count_missing=6/16` 全对齐。 |
| V27-14 `--sites-file` 与 `-t` 成并集 | confirmed | **已修**（`common.rs` + `filter/spec/region.rs`） | 修复前：仅 sites-file = 1 条；加 `-t chr1:300-500` = **4 条**。修法：`--sites-file` 落 `SiteSelector.sites` 独立族，与 region 族 AND。验收：`cli_contract::popg_and_stats_intersect_the_site_families`（popg 窗口/裸路由 + `stats --freq`）；真实数据 `--sites-file`(3 位点) + `-t Chr05:1-100` = **0 条**，`-t` 覆盖时 = 3 条。 |
| V27-15 `-r`+`-t` 成并集 | confirmed | **已修**（`SiteSelector` + `common.rs` + `run.rs` + `reports.rs`） | 修复前：popg 2 条、stats 2 行 vs `filter`/`bcftools` 0 条；违反 `cli.md` §5 的 V23-61。修法：新增 `reed_filter::spec::region::SiteSelector`（四族 `allows` 逐族 AND），`PopgenRequest`/`StatsReportRequest` 三字段收敛为 `select`；`-r` 只作 seek + 自家族。验收：单元 `site_selector_intersects_the_families` + CLI `popg_and_stats_intersect_the_site_families`。真实数据三方一致：1416 / 492 / 0（reed popg == reed filter == bcftools view）。 |
| V27-16 LD 线程平方增长 | confirmed | **已修**（单一 `LdExecutor` + 总预算拆分） | 每 run 只建一个 Rayon 池，contig 结果按序归并；`--threads` 拆给输入 pipeline 与 LD，不再各拿满额。实测 8 contig × 20000 位点：1/2/3/4/8/16/56 的线程峰值为 **1/1/3/4/8/16/56**，均 ≤ 预算（旧 56 → 441）。`ld_worker_share` 单测覆盖 1..64；反向恢复满额分配后失败。 |
| V27-17 LD 任务出错后不取消不 join | confirmed | **已修**（`LdExecutor::drop` + abort） | 失败 drop 触发共享 `AbortSignal`，块间检查退出，pool drop 等待 worker；成功 finish 不触发 abort。`e2e_ld::ld_error_cancels_workers_and_leaves_no_auxiliary_threads`：库级 `run_popgen` 返回 Err 后线程回到基线。 |
| V27-18 窗口累加器内存 | confirmed | **已修**（流式窗口会话 + 标量输出 + raw rows） | CLI 每窗完成即写临时文件，成功后统一发布；`WindowRow.outputs` 保存标量/审计值，库 API 默认保留完整 `accs`。Bp worker 传 `SiteRow`（`SmallVec`）而非每站点完整 `StatAcc`，sink 按位置分组逐行推进水位线。实测 2 样本 100k 个 1bp 窗口峰值 **约 25.0 MiB**（旧同类 4000 位点已达 669 MiB）；256 样本 50k 窗口 **79.6 MiB**；输出逐字节不变。反向恢复整段末尾水位线 → 首窗 `1 vs 100` 失败。 |
| V27-19 `--bed-file` 内存随 chunk 放大 | confirmed | **已修**（按 contig 共享 `RegionAccs`） | `StatsJob.region_accs` 按 contig 持有共享累加器，sink 采用同一 `Arc` 并只合并一次，跨 worker/seam 不重复；50k BED 区间、64 样本、--threads 8/--memory 64 峰值 **49.0 MiB**（旧 `--pi` 253 MiB、θ/Tajima 2.65 GiB）。跨 chunk/线程输出逐字节一致，gold `run_pixy_bed.sh` PASS。 |
| V27-20 `.ld` 不在 stderr 文件列表 | confirmed | **已修**（`PopgenReport.pairs_path` + CLI 汇总） | sink 把实际发布的流式 pair path 通过 `PopgenReport` 返回，CLI 不再猜路径；零配对表同样列入。验收：`empty_popg_outputs_replace_the_previous_run_and_are_listed` 的 stderr 断言。 |

## P2

| ID | 初审状态 | 实施状态 | 证据 / 验收 |
|---|---|---|---|
| V27-21 `--levels bogus` 静默通过 | confirmed | **已修** | sfs/ld 路由 exit 0；`--bed-file --levels window` 反而误报；`--f2 --bed-file` 完全不生效。 |
| V27-22 单群体 `--fst` 静默成功 | confirmed | **已修** | exit 0、无文件、`-> ` 空列表。 |
| V27-23 BED 未覆盖 contig 输出非法行 | confirmed | **已修** | `chr2 1 0`（`pos2 < pos1`）；debug 构建 panic exit 101。 |
| V27-24 `--contig-lengths` 被静默丢弃 | confirmed | **已修** | header 1000 + fai 400 → 仍输出 `401-500`；两处合并顺序相反。 |
| V27-25 不变位点检查的退出码 | confirmed | **已修** | exit 1，而 `cli.md:342` 写 exit 2；`cli_contract.rs:719` 断言 Some(1)。 |
| V27-26 `.log` methods 写死 π/dxy | confirmed | **已修** | `--theta` 日志仍印 pi/dxy；SFS 丢弃计数无消费者；裸 `--sfs` 不写 log。 |
| V27-27 `no_snps` 是记录数 | confirmed | **已修**（`FstAcc`/`PbsAcc::push_with_variant`） | 只在 `all_variant` 时递增 `no_snps`；默认 `push` 保留旧矩阵 API。验收：`fst_and_pbs_no_snps_count_only_variant_sites`、`e2e_pixy`（1 变异 + 2 不变 → 1）。 |
| V27-28 多等位绕过升序检查 | confirmed | **已修** | 双等位尾置 → exit 1；多等位尾置 → exit 0；`--memory`/BGZF 布局改变退出码。 |
| V27-29 SIGINT 留临时 `.ld` | confirmed | **已修** | 261 MB 残留；下次运行不清理。 |
| V27-30 `--skip-empty` 逐窗空转 | confirmed | **已修** | 100 Mb 单点、`--window 1` → 737 ms；`--levels genome` 同样 737 ms。 |
| V27-31 `RegionSet::stab` 退化 | confirmed | **已修** | 10000 区间相邻 85.7 ms vs 有隙 0.41 ms（209×）；每翻倍 ×4。 |
| V27-32 `closed` 线性扫描 | confirmed | **已修** | 40000 contig → 1.695 s（每翻倍 ~4.5×）。 |
| V27-33 `finish_bins` 每窗克隆 | confirmed | **已修** | 1 bin/窗 0.1288 s vs 100 bins/窗 0.3497 s（2.7×）。 |
| V27-34 jackknife 每 push 分配 | confirmed | **已优化（观察项）** | `Arc<[u8]>` 缓存 last block key，同块不再 `to_vec()`；不改变 SE/网格语义，端到端性能影响仍不作缺陷声称。 |
| V27-35 jackknife 块网格与加权 | confirmed | **已修** | 按 ADMIXTOOLS 当前块首位点跨度分块，非 `pos/block_bp` 或 `pos-prev_pos`。共享位置网格 + 每统计量有效位点加权；有序归并后统一分块，跳过不加权、合法零贡献计权。少于两个有效块 SE/Z=NA。跨所有切分点/线程数一致；真实 v7.0.2 D 自然估计与 SE 全精度一致到 2e-12。 |
| V27-36 ThetaAcc merge 形状 | code / 常规 CLI 不可达 | **已修（库 API 防御性）** | `ThetaAcc::merge` 校验 group/n_chrom、直方图内部等长及未 compact，再补零尾逐类相加。合法 `(3 vs 7)` 零尾在修复前 panic，现正反合并与整体累加的直方图和所有 θ 输出一致；异群体、异染色体数、compact、非法内部形状均拒绝。全局 `merge_zip_eq` 不放宽；θ 公式与正常 CLI 输出不变。 |
| V27-37 `--threads 1` 多一线程 | confirmed | **已修** | 无辅助池时 contig 在 sink 串行计算，实测 `--threads 1` 峰值 1；`--threads 2` 也保持 1。 |
| V27-38 裸命令忽略 `--variant-only` | confirmed | **已修** | `render_tsv` 使用 layout 策略；pair 做完整 comparable-pair 补全，bp 为 `Σsite_mean/L`；θ total=raw、bp=raw/L。缺失夹具裸 pair == 窗口 avg，pair/bp 定义分离；库级宽表与 CLI 回归通过。 |
| V27-39 默认 `-M 2` 门缺失 / 跳过不可见 | confirmed | **已修**（`StatsJob::site_prefix` + `popg.rs` + 删除 skipped-present） | 多等位在创建 segment/累加器前剔除并计数；`-M > 2` exit 2，默认与 `-M 2` 输出逐字节相同；窗口 stderr 和 `.log` 均报 `multiallelic skipped`。D7-11 作废；`stats` QC 的 collapsed-ALT 语义不受影响。 |

## 关联的既有决策（需在实施时同步更新）

| 既有 ID | 影响 | 动作 |
|---|---|---|
| D15-11 | V27-03 推翻其「分母累加」部分 | 拆测试为「单点 0.0」+「双位点 1/3」 |
| D7-11 | **V27-39 使其作废** | 多等位不再进入统计后「skipped present」机制整体删除；`n_present_file` ≡ `n_occupied`（pixy 口径） |
| D7-23 | V27-38 已重新定案 | `per_bp:false` 在 variant-only 下做完整配对补全；`per_bp:true` 是独立 `Σsite_mean/L`。裸 θ total=raw、bp=raw/L，不再重复或名不副实。 |
| V23-14 | V27-05/V27-20 的成因 | **已恢复其流式所有权并补齐空输出**：`pairs_path` 由 sink → driver → CLI 返回；不是回退到内存累积 |
| V23-61 | **V27-15 曾违反，第二批已恢复** | `SiteSelector` 把 `indexed` 与 `include`/`exclude` 分成独立族并逐族 AND，与 `reed filter` 的 `RegionRule` 结构一致；`cli.md` §7.3 新增了组合语义行 |
| V26-08 | V27-02 违反 | 把完整性校验移到 commit 之前 |
| D15-20 | V27-11 是其未覆盖的面 | 补 `--step > --window`、`--skip-empty`、BED 下的「窗口和 = contig 行」断言 |

## 验证边界（2026-11-09 第十批实施后更新）

- 第一批只跑了受影响 crate 的测试（36 个测试二进制约 499 项）；**第二批补齐了全 workspace**
  （当时 68 个测试二进制 1074 项全绿，含 `reed-io`/`reed-engine`/`reed-filter`）。
  **第八批后为 1111 项全绿，第九批后为 1116 项全绿，第十批后为 1126 项全绿**；各批均通过 `clippy --workspace --all-targets` 0 告警、
  `check-structure.sh` 16/16、gold `run_easysfs.sh` PASS，第三至五批另跑 `run_plink2_ld.sh` PASS；第五批另跑 `run_vcftools_fst.sh` PASS 与 pixy 小夹具逐列对拍；第六批另做线程峰值矩阵与库级取消验证；第七批另做流式内存与失败清理验证；第八批另做稀疏窗口/RegionSet 性能与 CLI 契约验证。
- 第十批另有 `run_allel_fstats.sh` point estimate PASS、真实 ADMIXTOOLS v7.0.2 动态网格及 D 自然估计/SE PASS，
  以及 easySFS、plink2 LD、vcftools Fst、pixy BED 全部 PASS。参考 `qpDstat` 默认打印 bias-corrected estimate，
  gold 使用其未改公式的全精度自然估计/SE 诊断；不把 reed 的自然估计与校正估计混为一谈。
- **未跑**：完整的外部工具矩阵和真实数据（`clean.bcf` 7 GB）性能矩阵；当前 gold 是仓库可用脚本与新增 ADMIXTOOLS 对拍。
- `check-doc-drift.sh` 原有的 **1 条既有 DRIFT**（根 `README.md` 缺 `--threads 上限 1024`，doc_v26 V26-03 遗留）
  已在第九批文档收尾补齐；当前 **98 ok / 0 DRIFT**，无行为变化。
- 未实测：`--min-dp`/`--min-gq` 掩码路径下的 variant-only 分母（与 V27-12 同源，预计同结论）。
- `--bed-file` 与 `--sites-file` 的组合语义：第二批在**谓词层**已确定（`--sites-file` 属 `SiteSelector.sites`，与 `--bed-file` 的窗口网格是两个正交概念——后者决定窗口，前者决定哪些位点进入统计；回归测试覆盖了 `--sites-file` + `-t`，但**未**构造 `--bed-file` + `--sites-file` 的端到端夹具）。
- V27-35 的 f2/f3/f4/f4-ratio 保持 reed 的 Σnum 或 Σnum/Σden 量纲；不以不同角色映射或 bias-corrected estimate 宣称全量 ADMIXTOOLS 输出相同。

## 实施批次

| 批次 | 条目 | 状态 | 门 |
|---|---|---|---|
| 第 1 批（2026-11-08） | V27-01、V27-02、V27-08、V27-09 | **已落地** | 见下 |
| 第 2 批（2026-11-08） | V27-14 + V27-15（区域谓词 AND） | **已落地** | 见下 |
| 第 3 批（2026-11-09） | V27-05、V27-06、V27-07 + V27-20 | **已落地** | 见下 |
| 第 4 批（2026-11-09） | V27-04（LD 单元素尾块） | **已落地** | 见下 |
| 第 5 批（2026-11-09） | V27-03/11/12/13/27/39（pixy 口径 + 多等位过滤期剔除） | **已落地** | 见下 |
| 第 6 批（2026-11-09） | V27-16、V27-17、V27-37（LD 线程与取消所有权） | **已落地** | V27-18/19 内存项另列 |
| 其余 P2 | V27-21…V27-38（除 V27-27/37/39 已修） | **已全部落地** | 第十批 V27-35/38 完成 |

## 2026-11-08 第一批实施记录

| ID | 改动 | 修复前复现 | 回归测试 | 反向验证 |
|---|---|---|---|---|
| V27-01 | `translate/popg.rs`（裸路由接 `refuse_output_eq_input`）、`cli/mod.rs`（`-o` 帮助） | `-o <输入>` exit 0，输入 md5 `2c96ad7b…` → `60815d1e…` | `cli_contract::popg_bare_route_refuses_output_naming_the_input` | n/a（新检查，无旧行为可对照） |
| V27-02 | `stats/sink.rs`（`pending_commit`）、`stats/driver.rs`（commit 移到 `job.finish` 之后） | 截断 BCF → exit 1 且 `.ld` 被覆盖成 4489507 行 | `e2e_ld::failed_run_leaves_the_previous_ld_untouched` | **已做**：`git stash` 掉修复后测试失败（`.ld` 内容不匹配） |
| V27-08 | `writers/sfs.rs`（`stem` + `unique_stems`） | `../../../../ESCAPED` → `/tmp/ESCAPED-6.sfs` | `e2e_sfs::sfs_group_names_cannot_escape_the_output_directory` | **已做**：把 `stem` 换回 `name.to_string()` 后测试失败并打印 `/tmp/ESCAPED-6.sfs escaped the tree` |
| V27-09 | 同上 + `claim(path, seen)` | 组 `A`+`4` → 一维谱被联合谱覆盖 | `e2e_sfs::sfs_one_dimensional_and_joint_names_cannot_collide` | **已做**：`git stash` 掉修复后测试 FAILED |

## 2026-11-08 第二批实施记录

| ID | 改动 | 修复前复现 | 回归测试 | 反向验证 |
|---|---|---|---|---|
| V27-15 | 新增 `reed-filter/src/spec/region.rs::SiteSelector`（`indexed`/`include`/`exclude`/`sites` + `allows` 逐族 AND，`is_empty`）；`PopgenRequest`/`StatsReportRequest`/`StatsJob`/`ReportJob` 的三个扁平字段 → 一个 `select`；`common.rs` 的 `-r` 只填 `select.indexed`（不再 push 进 target 列表） | popg `-r chr1:201-400 -t chr1:1-60` → 2 条；`filter`/`bcftools` 0 条；`stats --freq` 2 行 | `reed-filter::spec::region::site_selector_intersects_the_families`、`cli_contract::popg_and_stats_intersect_the_site_families` | **已做**：把 `allows` 的三处 `&&` 改回 `\|\|` → 单元测试 FAILED（`cross-family` 断言）、CLI 测试 FAILED（5 vs 0） |
| V27-14 | `common.rs` 的 `--sites-file` 落 `select.sites`（每行 `(contig, Some(pos), Some(pos))`），不再 append 到 `-t` 的 include 列表 | 仅 sites-file = 1 条；+ `-t chr1:300-500` = 4 条 | 同上（`sites_and_t` / `sites_and_t_hit` / `s_*` 四例） | 同上（同一个 `\|\|` 反向验证同时打破两条） |

**真实数据对拍**（`subvcf.vcf.gz`，记录 16004641–16998874）：

| 参数 | reed popg | reed filter | bcftools view |
|---|---:|---:|---:|
| `-r Chr05:16000000-16500000` | 1416 | 1416 | 1416 |
| `-r Chr05:16000000-16500000 -t Chr05:16400000-16600000` | 492 | 492 | 492 |
| `-r Chr05:16000000-16100000 -t Chr05:16900000-17000000` | 0 | 0 | 0 |

**第二批门**：`cargo clippy --workspace --all-targets` 0 告警；`cargo test --workspace` **1074 项 / 68 个测试二进制全绿**；
`check-structure.sh` 16/16；`check-doc-drift.sh` 97 ok + 1 条既有 DRIFT；`gold/popgen/run_easysfs.sh` PASS
（unfolded 8.635e-16 / folded 4.363e-16）。pixy gold 门无法跑通是**既有环境/夹具问题**
（pixy 对 `subvcf.vcf.gz` 报 `No genotype records found ... cannot infer ploidy`：header 声明 20 个 contig 而记录只在 Chr05），
与第二批改动无关。

**第一批时未做（已于第三批完成）**：V27-05（空结果输出会话）——它与 V27-02 共用「运行开始确定输出集合 + 统一发布」的设计，
单独改 `.ld` 会在 `_sfs/`/`.fstats` 留下同类缺陷，故与第 2/3 批一起做。

**过程教训（已固化进测试）**：V27-08 的测试初版只上跳三级（停在测试自己的目录内）且断言共享路径
`/tmp/ESCAPED-6.sfs`，两次误报——先被自己的验证运行污染，后又因上跳层数不足而无法区分修复前后。
最终版用 pid 唯一组名 + 自清理 + 断言「路径含 `..` 组件」与「树外文件不存在」两个条件。

## 2026-11-09 第三批实施记录

| ID | 改动 | 修复前复现 | 回归测试 | 反向验证 |
|---|---|---|---|---|
| V27-06 | `sink.rs` 跨块 seam 只比较上一位点；`next_start` 不再冒充输入顺序约束 | 40000 位点、长 ID 的 1.3 MiB BGZF 输入，`--window 100 --step 1M --threads 1` → `33224 followed by 33225` | `cli_contract::sparse_windows_accept_sorted_sites_across_a_real_chunk_seam`；sink 三个 seam 单测 | 恢复 `.max(cur.next_start)` → `a_site_between_sparse_windows_is_legal_across_a_chunk_seam` FAILED |
| V27-07 | `driver.rs` 在无窗口级 kind 时选择 `GridMode::None`；`sink.rs` 接受该模式并不再要求长度 | 无长度 header 的 SFS/f2/LD 全部 exit 1 | CLI 三类请求 + 库级无长度 SFS/f2/LD | 强制旧 Bp 分支 → `genome_wide_sfs_does_not_require_contig_lengths` FAILED |
| V27-05 | 请求级零值 SFS/fstats 初始化；LD 按请求产生全部 bins；pair 表运行开始创建，成功后统一发布 | 同前缀零配对运行保留旧 `.ld`、`.ld.decay`、`.fstats`、`_sfs/` | CLI 空结果覆盖旧文件；库级零谱、空 f2、空 LD | 恢复 `sfs: None` → 空 SFS 测试 FAILED；恢复延迟创建 `.ld` → 空 LD 测试 FAILED |
| V27-20 | `SinkOutput.pairs_path` → `PopgenReport.pairs_path` → CLI 文件清单 | 汇总只列 `.ld.decay` | CLI 成功汇总含 `out.ld` | 同一空 LD 反向验证覆盖 |

**第三批验收**：`cargo test --workspace` **1085 项 / 68 个测试二进制全绿**；
`clippy --workspace --all-targets` 0 告警；`check-structure.sh` 16/16；
`check-doc-drift.sh` 97 ok + 1 条既有 DRIFT；gold `run_easysfs.sh` PASS；
`run_plink2_ld.sh` PASS（15 pairs，max rel 1.020e-06；prune 3 kept / 3 pruned）。
`failed_run_leaves_the_previous_ld_untouched` 仍通过，证明提前创建 pair 临时文件没有回退 V27-02 的失败安全。

**已改变的旧测试**：`streamed_pairs_file_is_not_created_without_pairs` 更名为
`streamed_pairs_file_is_header_only_without_pairs`，断言从「无文件」改为「表头文件」；
`write_ld` 新增显式 `pairs: bool`，区分 decay-only 与 pair 请求，避免靠空 Vec 猜输出集合。

## 2026-11-09 第四批实施记录

| ID | 改动 | 修复前复现 | 回归测试 | 反向验证 |
|---|---|---|---|---|
| V27-04 | `reed-stats/src/matrix.rs::decay_range` 将 `len < 2` 早退改为仅 `start == end`；单元素块继续扫描前序 lookback | `decay_range(4096..4097)` 返回 0 对，应返回 4096；4096k+1 整族受影响 | 库级单元素测试；`stats::run` 的 4095/4096/4097/4098/8193 分块测试；CLI 4097 位点 `.ld` 与 decay 总数一致 | 恢复 `end.saturating_sub(start) < 2` → 单元素测试 **0 vs 4096** FAILED |

**第四批验收**：全 workspace **1088 项 / 68 个测试二进制全绿**；
`clippy --workspace --all-targets` 0 告警；gold `run_plink2_ld.sh` PASS
（15 pairs，max rel 1.020e-06；prune 3 kept / 3 pruned）。
真实 CLI 4097 位点 `.ld` = **8,390,656** 对，等于 `4097×4096/2`，与 decay `n_pairs` 总和一致。

## 2026-11-09 第五批实施记录

| ID | 改动 | 修复前复现 | 回归测试 | 反向/参照验证 |
|---|---|---|---|---|
| V27-03 | `fst.rs`/`pbs.rs` 的 Hudson between 项要求两组 `n>=2`，num/den 同进出 | 同一低覆盖夹具 allel/reed 0.2222、pixy 1/3 | `hudson_two_sites_match_pixy_one_third`；CLI pixy 低覆盖夹具 | 直接运行 pixy 2.2.4：Fst 1/3、`no_snps` 与 reed 一致；旧 oracle 不再充当 Hudson 数值权威 |
| V27-27 | `push_with_variant` 只对变异记录累加 `no_snps` | 1 变异 + 2 不变：reed 3、pixy 1 | `fst_and_pbs_no_snps_count_only_variant_sites`；`e2e_pixy` | pixy 输出同值；WC84/审计路径不变 |
| V27-12 | `PiAcc` 只把有调用记录计为 `|S|`；全缺失仍保留 missing 审计 | variant-only 分母把全缺失/多等位算作 present | occupied 单测、trailing-window 新分母、CLI 缺失夹具 | pixy π/no_sites/count_missing 逐列一致 |
| V27-13 | `DxyAcc` 只在两组都有调用时计 `no_sites` | π/θ/dxy 的 `no_sites` 不同口径 | `no_sites_is_the_pairwise_intersection`；CLI  pairwise occupancy | pixy `A=1, B=2, dxy=1`、missing 6/16 全一致 |
| V27-39 | gate 移到 segment 创建之前；`-M > 2` usage；删除 skipped-present；stderr/log 报计数 | 多等位进统计后才跳过，窗口不可见 | CLI `popg_rejects_multiallelic_mode_and_filters_multiallelic_records` | 默认 == `-M 2` 逐字节；`-M 3` exit 2；`multiallelic skipped: 1` |
| V27-11 | `WindowRow` 新增 `denom_bp`；contig/genome 声明长度与 BED 并集 | genome π 随 `--step` 改变 | step 100/1M 一致性测试；BED 150bp 并集测试 | `bp` 公开跨度不变，denominator 独立 |

**第五批验收**：全 workspace **1096 项 / 68 个测试二进制全绿**；
`clippy --workspace --all-targets` 0 告警；gold `run_easysfs.sh`、`run_plink2_ld.sh`、
`run_vcftools_fst.sh` 全 PASS。另以 pixy 2.2.4 运行低覆盖、全缺失两个小夹具，
π/dxy/Fst、`no_sites`、`no_snps`、`count_missing` 全部一致（仅浮点打印精度差异）。

**测试更新的理由**：D15-11 的旧 Hudson 断言改为 pixy 的 `None` 与双位点 `1/3`；
`oracle_parity` 只对 Hudson/PBS 改用独立 pixy 参照，WC84 和整数审计仍对拍 frozen oracle。
`trailing_skipped_site_advances_the_windows` 的分母由 D7-11 的 594/2994 改为无多等位 present 的 600/3000。

## 2026-11-09 第六批实施记录

| ID | 改动 | 修复前复现 | 回归测试 | 反向验证 |
|---|---|---|---|---|
| V27-16 | 新增每 run 一个 `LdExecutor`；contig 结果按序归并；`ld_worker_share` 把 worker 预算在 pipeline 与 LD 之间拆分 | 8 contig、`--threads 56` → 441 线程 | `thread_allocation_tests`；真实线程峰值 1/2/3/4/8/16/56 = 1/1/3/4/8/16/56 | 恢复满额分配 → 分配单测 FAILED |
| V27-17 | `LdExecutor::drop` 失败时 trigger abort；`decay_contig_with_abort` 块间检查；pool drop join；成功不触发 | 库 Err 后继续烧 CPU | `e2e_ld::ld_error_cancels_workers_and_leaves_no_auxiliary_threads` | 同测覆盖失败路径；成功 LD gold 不受影响 |
| V27-37 | 无辅助池时在 sink 串行执行 contig | `--threads 1` 峰值 2 | 实测峰值 1/1（threads 1/2） | 包含在线程矩阵中 |

**第六批验收**：全 workspace **1098 项 / 68 个测试二进制全绿**；
`clippy --workspace --all-targets` 0 告警；gold `run_plink2_ld.sh` PASS；
真实线程峰值全部 ≤ `--threads`。系统级 `ulimit -u 64` 未复测（线程池创建失败已映射为 `EngineError::Job`，非 panic）。
V27-18/V27-19 的窗口累加器与 BED 内存预算问题**尚未修**，另列后续批次。

## 2026-11-09 第七批实施记录

| ID | 改动 | 修复前复现 | 回归测试 | 反向验证 |
|---|---|---|---|---|
| V27-18 | 新增 `WindowOutputSession`（每窗临时文件 + 完整性后发布）；`StatOutput` 标量快照；Bp raw `SiteRow<SmallVec>`；逐位置/同 POS 分组发窗 | `--window 1` 大窗口完整累加器堆积，2 样本 4000 位点 669 MiB | `streaming_window_output_counts_rows_and_cleans_failed_sessions`、sparse-window seam、θ/Tajima 输出逐字节对比 | 恢复整段末尾水位线 → 首窗 **1 vs 100** FAILED；失败会话临时文件清理通过 |
| V27-19 | `StatsJob.region_accs` 按 contig 共享；sink 采用同一 Arc，避免重复 merge | 50k BED × 64 样本 × 8 线程，旧内存 253 MiB（π）/2653 MiB（θ） | `bed_accumulators_merge_once_across_chunks_and_workers`（>1 MiB 跨 chunk）、gold BED | 1/2/8 线程输出 md5 一致；峰值 **49.0 MiB**，且与修复前输出一致 |

**第七批验收**：全 workspace **1101 项 / 68 个测试二进制全绿**；
`clippy --workspace --all-targets` 0 告警；gold `run_easysfs.sh`、`run_plink2_ld.sh`、`run_pixy_bed.sh` 全 PASS。
实测 100k 个 1bp 窗口（θ+Tajima）峰值 **约 25.0 MiB**；256 样本 50k 窗口 **79.6 MiB**；
50k BED 区间 64 样本 **49.0 MiB**。仍保留每行/每站点必要数据与分配开销，未宣称整个进程 RSS 为常数；
`--memory` 继续只约束输入在途块，不能解释为全进程硬上限。

## 2026-11-09 第八批实施记录

本批完成 **13 条修复 + 1 条观察项优化**：V27-21、22、23、24、25、26、28、29、30、31、32、33、34，以及 V27-10。
每条均在 01-review 保留原始证据，并在上表记录实施与验收。

| 类别 | 条目 | 验收 |
|---|---|---|
| CLI/网格 | V27-21/22/23/24 | 非法 levels、单群体、BED 空 contig、长度优先级全部 CLI 回归通过 |
| 退出码/日志 | V27-25/26 | 不变位点 exit 2；动态 methods 与裸 SFS/LD `.log` 通过 |
| 输入契约/临时文件 | V27-28/29 | 多等位乱序仍拒绝；旧临时文件显式告警不误删 |
| 性能 | V27-30/31/32/33 | 稀疏窗口约 737 ms → 6–7 ms；RegionSet 索引；HashSet；scratch 重用 |
| 观察项 | V27-34 | `Arc<[u8]>` 缓存 last block key，不改变统计结果 |

**第八批验收**：全 workspace **1111 项 / 68 个测试二进制全绿**；
`clippy --workspace --all-targets` 0 告警；结构门 16/16；doc-drift 97 ok + 1 条既有 DRIFT。
旧 `e2e_pixy::writer_emits_pixy_columns_with_na_for_empty_windows` 的 variant-only 空窗断言由 NA 改为 0（V27-10 拍板语义），Require/Bypass 的 NA 不变。

## 2026-11-09 第九批实施记录

| ID | 改动 | 失败先行 | 回归测试 | 边界 |
|---|---|---|---|---|
| V27-36 | `ThetaAcc::merge` 把 Vec 长度与请求身份分开：相同 group/n_chrom、直方图内部等长且未 compact 时，省略的高 n 零类别补齐后逐类相加；不修改全局 `merge_zip_eq` | 合法零尾短表示与长表示 `(3 vs 7)`，修复前 panic `ThetaAcc sites: merge operands have different shapes (3 vs 7)` | `merge_accepts_equivalent_zero_trimmed_histogram_tails`：正反合并与整体累加的每个类别、raw/weighted θ、位点/变异/等位计数完全一致；4 个拒绝测试覆盖异 group、异 n_chrom、compact、内部非法形状 | 常规 CLI 未复现该 panic，不声称端到端正确性事故；θ 公式与正常 CLI 数值不变。Tajima 内嵌 θ 自动获得同一合并语义 |

**第九批验收**：全 workspace **1116 项 / 68 个测试二进制全绿**；`clippy --workspace --all-targets` 0 告警；格式检查通过。
结构门 **16/16**；`check-doc-drift.sh` **98 ok / 0 DRIFT**。原先根 README 的 `--threads` 上限遗漏已补齐，公开中英文手册同步为 1–1024；只修文档，不改变预算代码。


## 2026-11-10 第十批实施记录

| ID | 改动 | 失败先行与判别 | 验收 |
|---|---|---|---|
| V27-35 | 固定物理桶改为 ADMIXTOOLS 观察位点动态块；共享位置网格，按有效位点加权；`usable` 区分跳过与合法零贡献；有序 sink 归并后计算 | `0,5,10,15` 区分当前块首位点与相邻游走；8:2 权重独立公式探针；稀疏小于/恰好块长边界；共享/owned 所有切分点；缺 grid 误用拒绝 | `run_admixtools_jackknife.sh` 真实 v7.0.2 D 自然估计 `-0.099648300117233329`、SE `0.042291544590272186` 与 reed 一致到 2e-12；fstats 多线程一致；point estimate 不变 |
| V27-38 | 裸路由使用 layout invariant；`per_bp:false` 做完整配对补全，`per_bp:true` 保留 `Σsite_mean/L`；裸 θ total 输出 raw、bp 为 raw/L | A=`0/1 ./.`、B=`0/0 0/0`、L=100：π pair `1/595`、π bp `0.01`、dxy pair `1/398`、dxy bp `0.005`、θ total `1`、θ bp `0.01` | CLI 裸 pair 与窗口 avg 相等，pair/bp 不再冒充同一定义；库级宽表测试、普通 θ 期望同步为 raw |

**第十批验收**：全 workspace **1126 项 / 68 个测试二进制全绿**；`clippy --workspace --all-targets` 0 告警；
结构门 **16/16**；`check-doc-drift.sh` **98 ok / 0 DRIFT**。
`run_allel_fstats.sh`、`run_admixtools_jackknife.sh`、`run_easysfs.sh`、`run_plink2_ld.sh`、
`run_vcftools_fst.sh`、`run_pixy_bed.sh` 均 PASS。未跑完整外部矩阵及 7 GB 性能矩阵，见「验证边界」。

## 完成状态

V27-01…V27-39 **全部 39 项已实施**。V27-36 是库 API 防御性修复；V27-35 只改动态块/SE，保留 reed 自然估计量纲；
V27-38 明确裸表物理量与 θ raw/bp 列契约，不擅自把 pair 与 bp 强行做成不同数值。
