# doc_v15 / 01 核实台账：全库评审 40 条的实证复核

> 阶段一为 6 子代理静态评审（功能区划分见 README），本文件记录阶段二的
> 逐条实证核实。核实基线 HEAD `53c8aeb`，release 二进制
> `/tmp/reed-target/release/reed`（debug 用于溢出类对照），参照真值
> bcftools 1.24 / scikit-allel 1.3.13 / plink2（仅读其语义）。
> 判决用语：**实锤**（真实命令/构造成文件复现）｜**代码层核实**（逐行
> 确认、运行路径未触发或无法触发）｜**未复现**（构造后行为与指控不符）。
> 复现产物索引见 §7。

## 0. 判决总表

| # | 标题 | 判决 |
|---|---|---|
| D15-1 | Float FORMAT 掩码整体失效（`record_pipeline.rs:506-514`） | 实锤 |
| D15-2 | 掩码/`--set-gts .` 后 `-q/-Q` 仍读陈旧 INFO/AC,AN | 实锤（论据修正，§5-4） |
| D15-3 | BCF 头恰落 BGZF 块边界 → 记录流整体错位一字节 | 实锤 |
| D15-4 | `--phylip` 组合运行吞掉 `--012/--hapmap/--msa` 的位点 | 实锤 |
| D15-5 | `--012` 双等位门差一：三态位点被保留 | 实锤 |
| D15-6 | 仅 `--nexus-binary` 时 `--outgroup` 被忽略（未知名字也静默） | 实锤 |
| D15-7 | convert 全矩阵驻留内存（输出在 99.3% 时刻才出现） | 实锤 |
| D15-8 | 多 `--pbs` 三元组共用一行表头 | 实锤 |
| D15-9 | `--proj` 超过群体染色体数 → 全零 SFS、exit 0 | 实锤 |
| D15-10 | pop 文件同一样本两组 → 静默取后者 | 实锤 |
| D15-11 | Hudson Fst/PBS 打印 NA 而 oracle/allel 为 0.0 | 实锤（条件修正，§5-1） |
| D15-12 | oracle PBS 口径与 shipped/allel 不同，parity 测试靠种子运气全绿 | 实锤 |
| D15-13 | popg 跑完全量后才报缺少 `-o` | 实锤 |
| D15-14 | 裸 `popg | head` EPIPE panic（exit 101） | 实锤 |
| D15-15 | 裸 popg 忽略 `--variant-only/--bypass` 且跳过不变位点检查 | 实锤 |
| D15-16 | `--mind --indep-pairwise` 绕过 `--memory` | **未复现**（描述修正，§5-2） |
| D15-17 | `SiteKeepSet::contains` 每条记录 O(contig 数) | 实锤（24×） |
| D15-18 | `--f2/--ld-decay` 静默把 π/dxy/θ 切到 variant-only 口径 | 实锤 |
| D15-19 | 接缝处重复 POS 丢失于窗口行、计入 contig/genome 行 | 实锤（位置修正，§5-3） |
| D15-20 | 最后一条双等位位点之后的窗口永不输出但计入总数 | 实锤 |
| D15-21 | `--levels ""` 只写表头、exit 0 | 实锤 |
| D15-22 | `--bed-file` 下 `--step` 被静默丢弃 | 实锤 |
| D15-23 | `--tmp-dir` 在 stats/popg 上被接受但无人使用 | 实锤 |
| D15-24 | BcfDict 忽略 `##contig` 的 `IDX=` | 实锤 |
| D15-25 | `--rename-chrs` 允许两 contig 改成同名（非法 VCF） | 实锤 |
| D15-26 | MSA 长度 u8 截断（>255 破坏比对、跨样本错位） | 实锤 |
| D15-27 | POS=0 下溢：debug panic / release 输出 POS 0 | 实锤 |
| D15-28 | `--rm-dup-ids` 双重压缩（被动路径慢 15×） | 实锤 |
| D15-29 | gzip hapmap 走 stdin 被拒（路径读正常） | 实锤 |
| D15-30 | `--idspace-to` 与 `--set-missing-var-ids` 组合半生效 | 实锤 |
| D15-31 | oracle 的 f2/f3/f4 是朴素矩估计，与 shipped/allel 不同 | 实锤 |
| D15-32 | abort 通知丢失唤醒窗口 | 代码层核实（无外部触发源） |
| D15-33 | 重复组 PBS 报错文案错配（"invalid window grid"） | 实锤 |
| D15-34 | 测试缺口：`-q/-Q` 零 e2e、Float 掩码零夹具、窗口重叠零用例、搬迁重定向 1/27 | 实锤（grep） |
| D15-35 | 热循环分配（WcMultiAcc × 4、render/native_encode、g6、run.rs:492 等） | 代码层核实 |
| D15-36 | 死代码（JointSfsAcc、bin_of_pos0、Cursor::skip 等） | 代码层核实 |
| D15-37 | merge-law 测试锁死结构 `summary::SummaryAcc` | 代码层核实 |
| D15-38 | Fst/PBS 的 `missing` 恒为 0（第三参传了第二参） | 代码层核实 |
| D15-39 | `-o -` 生成字面文件 `-` | 实锤 |
| D15-40 | reader 线程不在 `catch_unwind` 下 | 代码层核实 |

## 1. D15-1…D15-15（必须修）逐条证据

**D15-1 Float 掩码失效。** `gq_int.vcf`（GQ=Integer）`--min-gq 20` →
`s1 0/1:5` 掩成 `./.:5`；`gq_float.vcf`（同数据、GQ=Float）同命令 →
原样输出 `0/1:5`。skip 子句 `!(x.is_nan() && mask_on_missing)` 对任何
非 NaN 值恒真；Integer 分支（:537）是正确形式。一行修复。

**D15-2 陈旧 INFO。** 判别夹具：掩码前 AF=2/8=0.25、掩码后 AF=1/6≈0.167
（s2 因 DP<100 被掩），阈值 `-q 0.2`。头部带 `AC=2;AN=8`：保留（用了
陈旧 0.25，且触发"may be stale"警告）；同基因型、头部无 AC/AN：丢弃
（重算正确）。规范 §3.2 明言掩码改变 cohort 后按掩码后基因型重算；
`plan.info_trusted` 的驱动清空条件漏了 `mask_if/set_gts`。

**D15-3 BCF 块边界。** 手构 BCF（`bcftools view -Ob` 原文件 + 一行可调
`##PAD` + 手工 BGZF 分块）：头文本（不含 NUL）恰 65280 字节时，reed 报
`BCF input ends inside a record`，输出只剩 `#CHROM` 头；`H=65279/65281`
与 `H=130560`（2 块）中仅 65280、130560（≡0 mod 65280）失败；bcftools 对
全部夹具 4/4 条正常。首记录 `l_shared=0x23` 位于裸偏移 65281，reed 从
65280（NUL 处）起读——与 `bgzf.rs:276` 的 `text.get(end)==Some(&0)` 在
跨块时为假吻合。

**D15-4 phylo 组合门吞位点。** `conv.vcf`（SNP/indel/SNP，3 样本）：
`--012` 单独 → `.pos` 3 行（100/200/300）；`--phylip --012` → 2 行
（200 缺失）。`push_line` 的 phylo 早退 `return` 跳过其余目标。

**D15-5 012 差一。** `A G,T`（无样本用 allele 2）→ `.pos` 含 pos 100
（被保留，按 REF/首 ALT 编码）；改 `0/2` → 整条记录被丢弃并打印
"malformed: a genotype index beyond the ALT alleles"。代码门 `>2`，注释
自述语义是"超过 2 个等位（REF+1ALT）即跳过"。

**D15-6 outgroup。** 4 样本：`--nexus-binary --outgroup s3` 行序
s1,s2,s3,s4（log 却记 `outgroup: s3`）；`--nexus --nexus-binary
--outgroup s3` 行序 s3,s1,s2,s4。`--outgroup nosuch` 在仅二进制模式
exit 0、输出与合法名字逐字节相同。

**D15-7 驻留内存。** 单线程峰值 RSS：`--012` 50k→32 MB、200k→24 MB、
800k→71 MB（约 77.7 B/位点）；`--hapmap` 800k×50 → 257 MB（约
141 B/位点，接近磁盘产物 118 MB）；运行 800k `--hapmap` 时输出文件在
5.152s/5.188s（99.3%）才首次出现。无 `--threads 1` 时另有引擎 chunk
arena 固定开销（300k→100 MB 量级）。

**D15-8 pbs 表头。** `--pbs A,B,C --pbs C,D,A` 的 `PREFIX_pbs.txt` 表头
为 `A B C chromosome ...`，第二行数据以 `C D A` 开头。行自带组名可自救，
按表头列位解析必错。

**D15-9 proj 超限。** 每群体 2 个二倍体样本（≤4 染色体）：
`--sfs --proj 40` → `dadi/A-40.sfs` 谱和 = 0.0（41 个全零），
fastsimcoal2 同样全零，exit 0、无警告；`--proj 4` 对照谱和 = 9.0。
`build_sfs_plan` 只查 n≥2 与每群体一个值。

**D15-10 重复分配。** `pops_dup.tsv`（s1→A、s1→B、s2→A、s3/s4→B）：
A 行 count_comparisons 54→9（只剩 1 个二倍体样本），B 行 13/54→40/135
（s1 归 B）。无任何报错。

**D15-11 Hudson/PBS NA。** B 组每站点恰 1 条被调用等位（1×`./0` +
其余 `./.`）：reed `--fst hudson` 窗口 = NA；allel `hudson_fst` 同数据
= 0.0（nansum 语义）；oracle = Some(0.0)。`--pbs` 同：reed 的 fst_12=NA
而 allel pbs=0.0。**触发条件修正见 §5-1。**

**D15-12 oracle PBS 口径。** C 组某位点全 `./.`、A/B 已调用：reed PBS
= 0.13118213223374553（与 allel 及项目 gold `compare_pbs.py` 完全一致，
max|diff|=0）；oracle = 0.11383935323980057、n_sites=5 vs 6（临时
Rust 测试直测对比，测后已删）。重放 `oracle_parity` 种子：996 位点、
0 次零调用组事件 → 现有断言从未被走到；`cargo test -p reed-stats
--test oracle_parity` 绿。

**D15-13 缺 `-o` 的检查顺序。** 300 万条记录、`--threads 1`：无 `-o`
0.8075s（随后 `window-layout mode writes files: pass -o PREFIX`），带
`-o` 0.8078s；clap 级用法错误基准 0.0059s。完整计算确实先发生。

**D15-14 EPIPE。** 600 contig、stdout 164054 字节：
`reed popg manycontigs.vcf --pop pops4.tsv | head -1` → exit 101，
`failed printing to stdout: Broken pipe (os error 32)`（
`translate/popg.rs:225` 的 `print!`）。≤64 KB 输出时因管道缓冲不触发。

**D15-15 裸 popg 丢检查。** 全变异文件（无 `ALT="."` 记录）：裸
`popg` exit 0 直接算；加 `--variant-only` 输出逐字节相同（标志无效）；
窗口路由 `--pi` 报 `input has no invariant sites ... declare
--variant-only`。裸路由 `layout: None` 使 sink 的检查与模式都不生效。

## 2. D15-16…D15-35（应当修）证据

- **D15-16 未复现**：三种组合在 `--memory 1` 下全部同样报
  `memory budget exceeded: in-flight 8388608 bytes > budget 1048576`
  （exit 1），`--memory 8` 下全部成功；找不到 (a)/(b) 行为分叉的配置。
  代码层：`matrix.rs:929-935` 的 `track_chunks` 分支确实先于体积检查
  return，但引擎按块 footprint 检查（pipeline.rs:408）对所有模式先触发；
  `mod.rs:245-252` 的 `choose_for_mind` 用 `default_half_memory()`
  （物理内存一半，本机 ≈257 GB），与 `--memory` 无关。**原描述前提错误，
  见 §5-2。**
- **D15-17 24×**：30 万位点、`--indep-pairwise 50 5 0.999`：1 contig
  0.146s vs 10000 contig 0.995s（默认线程）；`--threads 1` 下
  0.359s vs 8.708s。对照条件注意：两者保留位点数不同（35,726 vs
  68,706），1.9× 的记录差不足以解释 24×。
- **D15-18**：`pop8.vcf`（含 3 条不变记录）同一窗口 `--pi` = 0.37037；
  `--pi --f2 A,B` = 0.00167（diffs/comps 列不变）。全变异文件上
  `--pi --f2` 不再触发不变位点检查（`--pi` 单独触发）。
- **D15-19**：同文本 VCF 的分块 A/B——跨块（块边界恰在首个 POS 100 记录
  之后）：窗口 [1,100] no_sites=1/diffs=4/comps=6，窗口和比 contig 少
  −1/−4/−6（恰为重复记录贡献）；同块：窗口=contig=genome 一致；
  `nodup`/`seam99`（POS 99 重复）一致。**触发位置修正见 §5-3。**
- **D15-20**：双等位止于 250、三态在 500：窗口行止于 [201,300]，contig
  行 9/1794 与窗口和 0.005 不一致；仅含三态记录的 contig 只出 1 个窗口
  行；`--threads 1` 逐字节相同。
- **D15-21**：`--levels ""` → `out_pi.txt` 1 行（表头）、exit 0、log
  `windows: 0`。
- **D15-22**：`--bed-file w.bed --step 500` 与不加 `--step` 逐字节相同。
- **D15-23**：`--tmp-dir /nonexistent-dir` 正常完成；目录未被创建/触碰。
- **D15-24**：交换两 contig 的 IDX：bcftools `chr2 chr2 chr1 chr1`，
  reed `chr1 chr1 chr2 chr2`（exit 0）。
- **D15-25**：`chr1→X, chr2→X`：输出两行重复
  `##contig=<ID=X,length=1000>`，两 contig 记录合并；bcftools 宽容可读，
  文件非法。
- **D15-26**：REF 300bp 删除 → 输出宽度 1043（应为 1299，300 显示 44）；
  ALT 256bp → 宽度 1000（应 1255），且 0/0 样本取到 ALT 第 2 个碱基
  （跨样本错位）。exit 0。
- **D15-27**：debug `attempt to subtract with overflow`（job.rs:443，
  worker panic → exit 1）；release 接受记录并 `.012.pos` 写出 `chr1 0`。
- **D15-28**：3M 记录 `-Oz`：`--idspace-to '_'` 0.090s vs
  `--rm-dup-ids` 1.375s（15×），输出体积相同（5.71 MB）。
- **D15-29**：`cat hm.hmp.txt.gz | reed convert - --vcf` 报
  `not a hapmap file`；同文件按路径正常（2 records）。
- **D15-30**：ID 输出 `x:y → x_y` 而缺失 ID → `a:b`（半生效，无互斥
  报错）。
- **D15-31**：`--f2 A,B` = −0.125 = allel `patterson_f2` 之和（gold
  `compare_fstats.py` PASS）；oracle 朴素均值 = 0.078125；parity 测试
  的 kinds() 无 F2/F3/F4 条目。
- **D15-32 代码层**：`budget.rs:169-187` 的 waiter 在池锁下查 abort 后
  睡眠，`abort.rs` 的 `trigger`/`notify_abort` 不持该锁 `notify_all`；
  当前无外部触发源（无信号接线），内部错误路径有 permit drop 兜底。
- **D15-33**：`--pbs A,B,A` → `invalid stats request: invalid window
  grid: size and step must both be positive`。
- **D15-34 测试缺口（grep 实测）**：`e2e_filter.rs` 无 `min_af` 用例；
  Float FORMAT 夹具（XH，:619）只用于宽记录往返、未挂掩码规则；
  `e2e_stats.rs` 网格 size=100/step=100，全库无 step<size；
  moved-flag 单测（`translate/stats.rs:284`）只测 `--pi` 1/27。
- **D15-35 代码层**：`WcMultiAcc::push` 4 个 Vec/位点（fst.rs:287-331）、
  jackknife 每站点 clone contig 名（jackknife.rs:196-207）、
  render/native_encode 每条记录新建 Vec（render.rs:201/223/307、
  native_encode.rs:370/433）、`format_g6_into` 每浮点 2 次 String
  （line_edit.rs:279-305）、`push_fstats` 每站点 clone contig
  （run.rs:492）、single_pass/mid_stage 每记录 Vec。

## 3. D15-36…D15-40（建议）证据

- **D15-36 代码层**：`JointSfsAcc` 全库无生产调用方；`PattersonKind::
  Naive` 无 CLI 构造；`WindowGrid::bin_of_pos0` 仅自身测试调用且与
  `BpBinning` 的绝对 bin 约定冲突；`Cursor::skip`、`SiteView::info_*`、
  `HeaderScan::n_samples`、`keep_mode` 死参、`header_level` 弃值、
  `sink.rs:186` 死 `or(take)`、`convert/job.rs:390-397` 空臂。
- **D15-37 代码层**：`summary_merge_laws.rs` 导入 `summary::SummaryAcc`
  （mod.rs:20），产品路径为 `summary/report.rs:58` 的 `SummaryAcc`
  （被 `QcReportAcc` 与 reed writers 使用）。
- **D15-38 代码层**：`self.within_a.add(diffs, comps_a, comps_a)` —— 第三
  参传第二参，`missing = total − comparable` 恒 0。
- **D15-39 实锤**：`reed filter ren.vcf -o -` 在 cwd 生成 762 字节的
  字面文件 `-`。
- **D15-40 代码层**：worker（pipeline.rs:468）与 reducer（:589）有
  `catch_unwind`，reader 闭包内无；拟构造的 `ChunkSource` panic 未实测。

## 4. 判决汇总

实锤 32（必须修 15 全在内）｜代码层核实 6（D15-32/35/36/37/38/40）｜
未复现 1（D15-16，机制真实、描述有误）。无"证据不足"项。

## 5. 对原报告的修正

1. **D15-11 触发条件**：B 组 2–3 个 `./0` 样本**不会**触发（每组被调用
   等位 n_B≥2 时分子会累积）；须构造每站点恰好 1 条被调用等位（一个
   `./0`，其余 `./.`）。复核已按修正后的夹具实锤。
2. **D15-16 前提**：`choose_for_mind` 比较的是**物理内存的一半**，不是
   `--memory`；`--memory` 只进引擎的按块 footprint 检查。所有可达配置
   下三种模式的成败与报错完全一致——原"运行时不对称"不成立。降级为
   预防性修复（在 `track_chunks` 分支加上预算检查与体积模型补 mask 行）。
3. **D15-19 触发位置**：`--window 100 --step 100` 下丢失发生在
   **POS 100（pos0 99）**——bin 必须严格低于 `next_start/width` 才会被
   `drop_bins_below` 丢弃；POS 99（pos0 98）的窗口经 `close_contig`
   延后输出、仍能看到两份拷贝，不触发。
4. **D15-2 论据**：bcftools 无单命令等价式（`filter` 不能同时 `-e` 掩码
   与 `-i AF`）；两步等价管线（`filter -S .` 后再 `view -q`）同样保留
   该位点（它读文件内陈旧标签）。因此论据改为：**违反 cli.md §3.2 的
   明文承诺（"一旦…掩码…改变了 cohort，改按（掩码后的）基因型重算"）
   且与 reed 自家无 INFO 重算路径的结果不一致**。
5. **D15-17 对照口径**：两夹具保留位点数 35,726 vs 68,706（1.9×），
   24× 差距不可全部归因于记录数；若需更严格的可比口径，应对齐保留数
   后再测。

## 6. 测试缺口实测（D15-34 展开）

| 缺口 | grep 证据 | 复现的必须修 |
|---|---|---|
| `-q/-Q` 零 e2e | `grep min_af crates/reed/tests/e2e_filter.rs` = 空 | D15-2 |
| Float FORMAT 未被掩码规则覆盖 | XH 夹具仅宽记录往返（:619-） | D15-1 |
| 窗口重叠（step<size）零用例 | `e2e_stats.rs:86` size=100/step=100 唯一网格 | D15-19/20 |
| 搬迁重定向 1/27 | `translate/stats.rs:284` 仅 `--pi` 一条 | —（守卫缺口） |
| BCF 块边界零用例 | 无 `H≡0 (mod 65280)` 夹具 | D15-3 |
| 多 `--pbs` 三元组零用例 | 现有 pbs e2e 均单三元组 | D15-8 |

## 7. 复现产物索引

| 路径 | 内容 |
|---|---|
| `/tmp/reed-verify/`（gq_*.vcf、ac*.vcf、prec*.vcf、conv.vcf、tri_*.vcf、og.vcf、pos0.vcf、ren.vcf、pop8*.vcf、big300k*.vcf、big3m.vcf、d28_1m.vcf、d28_1m_ids.vcf（D15-28 分解：全 `.` ID vs 全唯一 ID）、manycontigs.vcf、ldA/ldB.vcf、c29/ 等） | 阶段二 CLI 级夹具与输出 |
| `/tmp/verify-bcf/`（build_fixtures.py、seam_build.py、run_all.sh、run_seam.sh、boundary*.bcf、idx*.bcf、seam*.vcf.gz） | BCF 块边界、contig IDX、接缝重复 POS 的手工 BGZF 夹具与脚本 |
| `/tmp/verify-oracle/`（f11*.vcf、f12.vcf、f31.vcf、check*.py、scan_seeds.py、/tmp/allelvenv） | oracle/allel 对照夹具与种子扫描 |
| `/tmp/verify-misc/`（A/ B/ C/ D/ 子目录） | 窗口/multi-only、MSA u8、内存缩放、--memory 对照 |
| 复核二进制 | `/tmp/reed-target/{release,debug}/reed`（0.4.0，v0.3.0-105-g53c8aeb） |