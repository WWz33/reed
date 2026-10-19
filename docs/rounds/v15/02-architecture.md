# doc_v15 / 02 修复规范

> 每条给出：落地改动（文件/位置/形态）、语义影响、回归测试（必须能失败）、
> 需要回改的规范章节。编号沿用 01 的 D15-N。修复纪律：**每条必须附带
> 一条在修复前会失败的回归测试**（D15-34 证明现有套件对 15 条必须修全
> 部失明）；凡改变 cli.md 承诺行为者同轮回改 spec。

## F1 机械修复（一行级，先落地）

### D15-1 Float 掩码
- 改动：`crates/reed/src/filter/record_pipeline.rs:509`：
  `|| !(x.is_nan() && mask_on_missing)` → `|| (x.is_nan() && !mask_on_missing)`
  （与 :537 的 Integer 形式对齐）。
- 测试：`e2e_filter.rs` 新增 Float `GQ` 头部夹具 + `--min-gq 20`，断言
  低 GQ 基因型变 `./.`（hand fixture 见 /tmp/reed-verify/gq_float.vcf）。

### D15-5 `--012` 等位门差一
- 改动：`crates/reed/src/convert/job.rs:755` `raw_alt_count(lr) > 2` →
  `> 1`（注释语义：>2 个等位 = REF+1 ALT 即跳过，vcftools 行为）。
- 测试：`e2e_convert.rs` 新增 `A G,T`（无 allele-2 使用）夹具，断言该
  位点不在 `.012.pos`；保留一条 ALT='.' 的既有行为。

### D15-21 `--levels ""`
- 改动：`crates/reed/src/stats/request.rs` `Levels::parse`：无 token
  命中时报用法错（`--levels lists no level (want window, contig, genome)`），
  与 `--hist ""` 的既有先例一致。
- 测试：`e2e_stats.rs` 或 cli 契约测试断言 exit 2 且信息含 `--levels`。

### D15-22 `--step` 与 `--bed-file` 互斥
- 改动：`crates/reed-cli/src/cli/mod.rs:558-560` 给 `step` 增加
  `conflicts_with = "bed_file"`（`--window` 已有）。
- 测试：cli_contract 断言 `popg --bed-file x.bed --step 500` 解析期报错。

### D15-27 POS=0
- 改动：`crates/reed/src/convert/job.rs:439-443`（parse_pos1 消费处）：
  `pos1 == 0` 按 malformed 记录处理（跳过 + 计数），去掉 `pos1 - 1`
  的下溢路径。
- 测试：debug 构建下 `POS=0` fixture 不得 panic，输出不含该位点、
  计数入 malformed。

### D15-30 `--idspace-to` 互斥漏项
- 改动：`crates/reed/src/rename/mod.rs:603-607` 的互斥检查加
  `set_missing_var_ids.is_some()`。
- 测试：`e2e_rename.rs` 断言两标志同时给出为用法错（exit 2）。

### D15-33 `--pbs` 重复组
- 改动：`crates/reed-cli/src/translate/popg.rs` `parse_pbs_triples`：
  三个名字互异校验，重复时报错点名（`--pbs A,B,A: population A repeats`）。
- 测试：cli_contract / e2e 断言报错含 `--pbs` 与点名，不再落进
  "invalid window grid" 文案。

### D15-13 `-o` 校验顺序
- 改动：`crates/reed-cli/src/translate/popg.rs:109-116`：把 `args.out`
  检查移到 `run_popgen_opened` 之前（与 `stats` 的预检同位）。
- 测试：既有 stdin_streaming/契约测试补一条 300k 级冒烟不由本轮引入；
  单元级断言顺序（错误在打开输入/计算之前返回）。

### D15-14 EPIPE
- 改动：`translate/popg.rs:225` `print!("{text}")` →
  `io::stdout().write_all(text.as_bytes())?`（走既有 `Error` 通道，
  EPIPE 变为普通运行错误 exit 1，不再 panic）。
- 测试：`stdin_streaming.rs` 或 e2e 冒烟：600 contig 夹具管 `head -1`，
  断言非 101 并含 broken pipe 文案。夹具生成见 /tmp/reed-verify/manycontigs.vcf。

### D15-39 `-o -`
- 改动：filter/rename/convert 的输出打开处把 `-` 识别为 stdout。若
  下游 writer 只支持路径流，最小正确形态是**显式拒绝**（用法错：
  `-o -: stdout output is not supported`）并在 help 写明；若 BlockWriter
  已有 stdout 分支则直接接线。落地时以现有 writer 能力定，二者都消灭
  "字面文件 `-`"。
- 测试：断言 `-o -` 不再产生名为 `-` 的文件。

## F2 统计口径

### D15-2 mask 后 `-q/-Q`
- 改动：`crates/reed-filter/src/plan/mod.rs`（compile 收口处，约 :897）
  返回 `info_trusted: !cohort_changed`；删除
  `crates/reed/src/filter/mod.rs:399-401` 的驱动清空块（单一来源）。
- 影响：掩码与 `--set-gts .` 之后的 cohort 变化使 `-q/-Q` 一律按掩码后
  基因型重算（cli.md §3.2 已明文）。
- 测试：e2e 以 /tmp/reed-verify/precac.vcf 为形：`--min-dp 100 -q 0.2`
  必须丢弃（掩码后 AF=1/6）；另加 `--set-gts .` 变体。

### D15-11 Hudson/PBS 的 `num_any` 门槛
- 改动：`crates/reed-stats/src/fst.rs:216-222` 删 `|| !self.hudson_num_any`；
  `pbs.rs:51-57` 删对应 `|| !self.num_any`（保留 `den != 0`）。
- 语义：对齐 allel `nansum`（全 NaN 分子 → 0.0）与 oracle。
- 测试：reed-stats 单测 + e2e `f11d.vcf` 形（B 组 1 条被调用等位）断言
  输出 0.0 而非 NA。

### D15-12 oracle PBS 口径
- 改动：`crates/oracle/src/estimators/pbs.rs:227-231` 去掉三组交集门
  （每对独立累积），`n_sites` 计所有站点；与 shipped/allel 同口径。
- 测试：oracle_parity 增加零调用组用例（把 §1 的 C 全 `./.` 夹具固化），
  修复前会红。

### D15-31 oracle fstats
- 改动：`crates/oracle/src/estimators/fstats.rs` 的 `f2/f3/f4/d_stat`
  移植偏差校正式（与 `reed-stats/src/patterson.rs:44-65` 同式）：
  `(pA−pB)² − p̂A(1−p̂A)/(anA−1) − p̂B(1−p̂B)/(anB−1)`；跳过 `an<2`。
- 测试：oracle_parity 增加 F2/F3/F4 条目（kinds() 现无）。

### D15-18 `--f2/--ld-decay` 口径切换
- 改动：`crates/reed-cli/src/translate/popg.rs:321-336`：仅当请求里
  **没有** π/dxy/θ 类 kind 时才强制 VariantOnly；有窗口统计时保留
  Require/Bypass，组合不合法则显式报错。
- 测试：`pop8.vcf` 形 e2e：`--pi --f2` 的 π 值与单独 `--pi` 一致；
  全变异文件 + `--pi --f2` 不再静默通过（转为检查报错）。

### D15-9 `--proj` 上限
- 改动：`crates/reed/src/stats/driver.rs:66-91`（group_sizes 已知后）：
  `projections[g] > 2 * group_sizes[g]` 报用法错，点名群体与两数。
- 测试：/tmp/reed-verify 的 4 染色体夹具断言 exit 2 且信息含 `--proj`；
  `--proj 4` 保持成功。

### D15-10 重复分配
- 改动：`crates/reed/src/stats/request.rs:875-895` `group_masks`：样本已
  有槽位时报案名错（`sample X is assigned to both A and B`）。
- 测试：e2e 用 `pops_dup.tsv` 形断言报错。

## F3 格式与窗口边界

### D15-3 BCF 块边界
- 改动：`crates/reed-io/src/chunk/bgzf.rs:269-281` `scan_header`：BCF 且
  `end == text.len()` 时不返回，读下一块（`read_block` 已跳过空成员）后
  重查 NUL；carry 变为该块、`first_in_block_offset = 1`（无 NUL 则 0）；
  流终止则保持现状。
- 测试：固化 /tmp/verify-bcf/boundary.bcf 的两块手工夹具
  （H=65280 与 H=130560 各一），断言记录数与 `bcftools view` 一致。

### D15-4 组合门吞位点
- 改动：`crates/reed/src/convert/job.rs:457-526`：把 phylo 段抽成独立
  方法（三个早退 `return` 方法内化），各目标独立求值。
- 测试：`e2e_convert.rs` 固化 `--phylip --012` 夹具（含 indel），断言
  `.012.pos` 与单独 `--012` 相同。

### D15-6 outgroup
- 改动：`outgroup` 提升到 `ConvertRequest` 级（对每个比对目标生效）；
  bin-only 计划读取它，未知名字走既有校验（job.rs:1407 同款）。
- 测试：`--nexus-binary --outgroup s3` 行序 s3 首位；`nosuch` 报错。

### D15-8 多 `--pbs` 三元组
- 改动：落地时按 cli.md 的可重复性承诺二选一：(a) 每三元组单独文件
  （文件名带组名后缀）；(b) 保留单文件但表头改为不误导。默认取 (a)，
  与 `--pbs` 可重复的文档一致；实现前读 cli.md §7/§help 与
  `gold/popgen/compare_pbs.py` 的解析器再钉死格式。
- 测试：双三元组夹具断言两个文件、各自表头与行一致。

### D15-15 裸 popg 检查
- 改动：`crates/reed-cli/src/translate/popg.rs:179-196`：route 2 与
  windowed 路由同样计算 invariant 模式（variant_only/bypass/Require）
  并传 `layout: Some(...)`。
- 测试：全变异文件裸跑 `--variant-only` 与裸跑输出不一致（前者按
  variant-only），不带 flag 时按 Require 报错（与 `--pi` 一致）。

### D15-19/D15-20 窗口接缝与尾窗
- 改动 a（D15-20）：`crates/reed/src/stats/run.rs:435-444` 的
  multiallelic 分支同时更新 `segment.first_pos0/last_pos0`（min/max）。
- 改动 b（D15-19，最终实现）：`emit_complete_windows` 的发射条件改为"出现**严格更后**位点才
  emit"（`watermark = cur.last_pos0`），窗口在其末位被越过前保持开放，跨块的同 POS 重复因此落入
  仍开窗口；接缝守卫保留为"不得低于已发射水位线"的不变量检查（`first_pos0 <
  cur.last_pos0.max(cur.next_start)`）。依据 VCFv4.3 Fixed fields 原文"It is permitted to have
  multiple records with the same POS"。回归测试：`stats::sink::tests::same_pos_seam_duplicate_lands_in_its_window`。
- 测试：固化 /tmp/verify-misc/A 与 /tmp/verify-bcf/seam.vcf.gz 两个
  夹具：窗口和必须等于 contig 行（修复前 seam 少 −1/−4/−6，A 少
  一条尾窗）。

### D15-24 contig IDX
- 改动：`crates/reed-io/src/bcf/dict.rs:46-51,105-110`：`##contig` 行按
  与 tag 行相同的 `IDX=` 规则解析（越界跳过、稀疏映射）。
- 测试：固化 /tmp/verify-bcf/idx.bcf，断言 CHROM 与 bcftools 一致。

### D15-25 contig 改名撞名
- 改动：`crates/reed/src/rename/mod.rs:429-439`：contig_map 构建后对
  最终名做重名检查，报错（镜像 `add_sample_pair` 的样本检查）。
- 测试：`e2e_rename.rs` 断言 `chr1→X, chr2→X` exit 2，无输出。

### D15-26 MSA 长度
- 改动：`crates/reed/src/convert/job.rs:180-183,246,947-968`：
  `msa_lens` → `Vec<u32>`、`MsaSite::width` → `u32`，去掉 `as u8`。
- 测试：REF 300bp 删除夹具断言输出宽度 1299、无跨样本错位。

### D15-29 gzip hapmap stdin
- 改动：`crates/reed/src/convert/hapmap.rs:493-515`：stdin 分支嗅探前两
  字节 `1f 8b` 时套 `MultiGzDecoder` 并把嗅探字节链回。
- 测试：`cat *.hmp.txt.gz | reed convert - --vcf` 成功，计数与路径读相同。

## F4 规模项

### D15-17 SiteKeepSet
- 改动：`crates/reed/src/filter/pass0/prune.rs:38-46`：pass1 前构建
  `HashMap<&[u8], &[u64]>`（或按名排序 + 二分）。
- 测试：性能门（1 万 contig 夹具单线程目标 ≥5× 改善）＋既有 e2e 全绿。

### D15-28 `--rm-dup-ids` 双重压缩（按两阶段完整形态落地）
- 形态（`crates/reed/src/rename/dedup.rs` + `mod.rs` 的 sink）：
  **阶段一（串行，有序 sink）**只做全局定序必需的事——按记录序给重复 ID
  加后缀、把改写后的字节填进**按字节封口**的批次（目标 1 MiB，整记录）；
  **阶段二**由 sink 私有管线池（`budget.auxiliary()` 个线程）并行
  deflate 已封口的批次，按批次序回交 `OrderedOutput`；批内切块固定
  65,280 字节，锚点按"记录起始所在块"重新定位。在途批次窗口 = 池大小
  + 2，且不超过运行内存预算的 1/4（下限 2 批）——池的批次缓冲也被
  `--memory` 兜住。
- 关键决定：
  - **封口只看累计字节**（不看时间/线程数），因此块布局与压缩字节流对
    任意 `--threads` 逐字节一致（测试断言含 sidecar 索引）；
  - worker 侧在 `rm_dup_ids` 下不再做无用 block push，但把**索引键随
    partial 带走**（一记录一键、按序；worker 本已解析该记录），有序 sink
    不再重复解析——键流错位按错误退出，不静默错配；
  - 无索引（无 `-W`）时 sink 不解析键（与 worker 路径同款的 `need_key`
    短路）；
  - 输出不足一批（小文件）时在 `finish` **内联 deflate**，不为一个小
    作业起线程池；`:308` 的 `dedup_line` 改 `Cow`，`.` ID 与未重复 ID
    的常用路径零分配；
  - 池随 sink 生命周期：`Drop` 先断开作业队列再 join，中止路径不悬挂。
- 实测（3M 记录 `-Oz`，tmpfs、min-of-5）：`--rm-dup-ids` **1.394s →
  0.23s**（对照 `--idspace-to _` 0.091s：15.4× → **2.5×**）；`-Oz -W`
  0.820 → **0.515s**（对照 -W 0.299s）；1M 全唯一 ID（HashSet 主导）
  0.371s。decompressed 输出与改前**逐字节一致**（仅 `##reed_command`
  行按命令行不同）。
- 测试：`dedup.rs` 单元（块切分/锚点重定位/往返、在途窗口的预算算式、中止路径 Drop 不悬挂）；e2e
  `dedup_batches_keep_order_index_and_bytes`（80k 记录跨多批 + `-W`：
  后缀跨批全局有序、索引查询与扫描一致、**压缩流与 sidecar 跨线程数
  逐字节相同**）；e2e `dedup_single_batch_matches_the_plain_run`（单批
  内联路径，含顺序引擎 `--threads 2`）。
- 残留：`-W` 时的键解析与 CSI 记账仍在 sink 串行（后者与 `--idspace-to`
  对照同样支付），大 ID 集下的 `HashSet<Vec<u8>>` 分配未优化。

### D15-7 convert 驻留内存（架构级；已按完整形态落地）
- 形态：每个目标的 merged partial 由 sink 在合并前接管并落临时文件
  （追加式）：hapmap rows、`--012` cells+pos、phylo cols/bin_cols/used、
  MSA 帧流（`pos0 u64 LE ‖ width u32 LE ‖ vals`）；两级阈值（单缓冲
  4 MiB / sink 总量 8 MiB），render 时流回，`RowsBlob` 的 Drop 清理
  成功/失败/中止三路径。
- MSA 渲染：逐样本一趟帧遍历，**按样本 seek 只读该样本切片**——不做
  O(样本数×数据) 的整帧重读（v13 D13-2 的"禁止重复 IO"决定在 v15 的
  延续，守卫见 ci/check-doc-drift.sh）。
- 实测（800k×50，1 线程）：hapmap 123.4→15.6 MB、`--012` 70.8→26.0 MB、
  MSA 151.4→15.8 MB，三者输出与改前逐字节一致。
- 运行期"输出早现"：**保持现状**（输出仍在成功后写）；cli.md 对此无
  承诺，运行期直写会打破"失败不留半成品"的既有行为；如日后要做，按
  `.partial` + 成功 rename 另开。
- 测试：跨阈值 e2e（hapmap；`--012`+phylip；MSA 60k×80）＋ 4 个单元
  （append/visit/writer 流回/清理）。

### D15-16（预防性）
- 改动：`pass0/matrix.rs:929-935` 的 `track_chunks` 分支在 push 前做
  与流式分支同款的体积检查（超预算报 `MemBudgetExceeded`）；
  `strategy.rs:48-53` 的 `sketch_size_bytes` 补 mask-plane 行。
- 测试：既有 --memory 门不回退；新增单测覆盖 track_chunks 分支的预算
  超限路径。

## F5 收口

- D15-36 死代码删除（`JointSfsAcc`、`bin_of_pos0`、`Cursor::skip`、
  `SiteView::info_*`、`HeaderScan::n_samples`、`keep_mode` 死参、
  `header_level` 弃值、`sink.rs:186`、`convert/job.rs:390-397` 空臂）；
- D15-37 merge-law 测试改指 `summary/report.rs` 的 `SummaryAcc`；
- D15-38 两处 `add` 调用去死字段或传真实 total；
- D15-32 触发路径加锁（预防性）；
- D15-35 只做已标点的两处（`run.rs:492` 去 clone、`push_sfs` 权重缓冲）；
- D15-34 每个必须修的回归测试清单随 F1–F3 落地；
- cli.md 回改点：§3.2（mask 后重算的驱动条件已明，补"含 `--set-gts .`"
  的落地注）、§5（裸 popg 的布局行为）、§7（`--pbs` 多三元组输出形态、
  `--proj` 上限报错）；`--help` 文案（D15-14 broken pipe 语义、
  `-o -`）；CHANGELOG 0.4.1（缺陷修复）；docs/README.md 索引行。

## 通用门（每批）

1. 对应 crate 的 `cargo test` 全绿；
2. 该批每条发现的最小复现从"复现"转为"不复现"（用 §7 夹具重跑）；
3. `cargo clippy --workspace --all-targets` 零告警；
4. 涉及输出字节的（convert/rename/filter 输出、pbs 布局）与修复前逐字节
   对照：除本条修复的差异外不得变化。