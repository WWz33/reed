# 01 评审：全库 must-fix（R1–R14）

> 评审对象：分支 `dev` 工作树（HEAD `5414326` + doc_v18 未提交 WIP）。
> 方法：11 个只读评审（不编译不改码）覆盖全部 crate 与跨 crate 接缝；
> 候选 must-fix 逐条由复审者**读源码行级核实**后，再以"先写失败测试"的方式
> 在 v20 分支上复现。本文件只收录复现成功的条目；每条附修复前的实测输出【据】。
> 完整报告（80 条发现 + 干净项清单）在仓库外 `../../../reed-rs-review.md`。

---

## 0. 总评

骨架是好的：经典统计量全部手工验算正确（WC84 手算到精确分数 0.125 与
9/53；Tajima n=4 的 0.0742269619…；pixy variant-only 的 3/60、6/160），
oracle 是真独立参照（零依赖、scikit-allel 1.3.13 黄金值内嵌），确定性契约
（线程数无关、字节级 e2e 锁定）在 reed-engine 里靠票据窗口构造性成立，
错误链在 58 处 `map_err` 里只有 3 处丢源（已列 backlog）。

must-fix 集中在四类：

1. **静默错果**（R1、R6、R7、R8、R9、R10、R12）：过滤丢弃该保留的记录、
   截断文件产出错误矩阵、SFS 质量静默归零、t-SNE 不是 t-SNE、GRM 整体
   膨胀、PCA 换了估计量、MSA 错位。
2. **崩溃 / 数据丢失**（R2、R3、R4、R5）：两条畸形输入 panic 路径、一个
   flag panic、四个子命令 `-o` == 输入时 exit 0 并毁掉输入。
3. **规模**（R13、R14）：距离矩阵族串行 O(n²·m) 且 `--threads` 被忽略；
   `--bed-file` 每位点线性扫全部区间。
4. **一条契约违背**（R11）：rename 模块自己声明的"拼写错误不得静默跳过"
   被 `--set-id` 早返回绕过。

---

## 1. 逐条：缺陷、复现、修法

格式：**R#（crate，严重度）标题** `file:line` — 问题；修前实测【据】；修法。

### 崩溃 / 数据丢失（P0）

**R2（reed-io，P0）空 FORMAT 列命中缓存哨兵，scanner 越界 panic**
`crates/reed-io/src/vcf/format_layout.rs:77-79` + `vcf/genotype_scan.rs:146`
— `FormatLayoutCache::new` 播种 `last: (Vec::new(), usize::MAX)`；`intern(b"")`
满足 `self.last.0 == format`（`Vec::new() == b""`）返回 `usize::MAX`，
`genotype_scan.rs:147` 用它直接索引布局表。一条 9 列 + 空 FORMAT 的坏行
（`…\t.\t\t0/0`）通过 `LineRef::parse` 后中止整跑。
修前【据】：`r2_empty_format_column_does_not_panic` panic 于
`format_layout.rs:126`：`index out of bounds: the len is 0 but the index is 18446744073709551615`。
修法：快路径加哨兵守卫 `self.last.1 != usize::MAX`；空 FORMAT 走慢路径
成为无 GT 布局（全部样本解码为缺失）。
同型隐患（未触发，列 backlog）：`ContigTable::id_of`（`vcf/site.rs:61`）
对空 CHROM 列会返回伪 `ContigId(u32::MAX)`。

**R3（reed，P0）BCF 记录 tid 越界直接索引 `contig_names`，panic 逃出
顺序路径** `crates/reed/src/index.rs:449`（tid 来自文件字节，`:321-330`
只查负值）— 文本路径按名字解析 contig 并报错；BCF 路径信任数值 tid。
`--threads 1/2` 走 `run_sequential`（无 `catch_unwind`），panic 逃出承诺
返回 `Result` 的 `run_index`。
修前【据】：`r3_sink_rejects_tid_beyond_the_header_contigs` panic 于
`index.rs:449`：`the len is 1 but the index is 7`。
修法：`IndexSink::push` 顶部边界检查（覆盖 `:433` 错误路径与 `add`），
报 `BCF record references contig id {tid} not declared in the header`。

**R4（reed-cli，P0）`-O ''` 三处 `split_at(1)` panic（exit 101）**
`translate/filter.rs:365`、`rename.rs:56`、`convert.rs:57` — clap 接受
空串，三处解析器直接 `split_at(1)`。
修前【据】：三个 `r4_empty_output_type_is_a_usage_error` 均 panic：
`end byte index 1 is out of bounds for string of length 0`。
修法：共享 `split_output_type(spec, spellings)`（common.rs），空串报
usage，错误消息列出合法拼写。

**R5（reed-cli，P0）`-o` == 输入：先读后截断，四个子命令 exit 0 毁掉输入**
全仓库无路径相等检查（grep 证实）；输出 `File::create` 晚于输入打开。
修前【据】（临时目录副本上实测）：`filter -o 同名 -Oz` exit 0、输入被
替换（242→608 B）；`rename` exit 0（564→697 B）；`convert` exit 0
（hapmap 变成 `##fileformat=VCFv4.3` 开头）；`index -f -o 同名` exit 0、
612 B BGZF VCF 被 120 B CSI 覆盖。filter 大文件形态下是被信号杀（缓冲
读未完成即截断）。
修法：共享 `refuse_output_eq_input`（canonicalize 输入；输出存在则
canonicalize、否则 `std::path::absolute`；`-o -` 跳过；符号链接缺口
在文档注明），接入 filter/rename/convert/index 四处；convert 在 `--vcf`
分发前检查（prefix == 输入同类错误）；index 在 exists/`-f` 门**之前**
（`-f` 不得成为绕过手段）。pca/distance/popg 写 `PREFIX.*` 后缀文件，
不覆写输入本身，不在本条范围。

### 静默错果（P0）

**R1（reed-filter，P0）`-f/--apply-filters` 匹配整个 FILTER 列**
`crates/reed-filter/src/eval/head.rs:86` — `set.iter().any(|s| s ==
line.filter_text())` 拿列表项与**原始整列**比较；README 承诺"contains a
value from LIST"，bcftools 按 `;` 拆分逐值匹配。`FILTER=LowQual;LowMQ`
（`bcftools filter -m +` 的常规产物）被 `-f LowQual` 静默丢弃。
修前【据】：`r1_filter_list_matches_any_semicolon_separated_value` 断言
`require_filter_passes("LowQual;LowMQ", &["LowQual"])` 失败。
修法：列按 `;` 拆分、任一元素等于任一列表项即通过（零分配子切片迭代）；
保留原正确语义：`.` 仅在显式列出时匹配 `.` 列、单值列精确匹配、
PASS 与具名值互不匹配。`HeadRule::RequireFilter` 的旧文档（"`.` never
passes"）本就与 README 漂移，一并改正。

**R6（reed，P0）pca/distance 静默跳过畸形记录**
`crates/reed/src/genotype.rs:101-112`（run_chunk）、`:179-190`（run_seam）
— `LineRef::parse` 失败、样本列数不符、POS 不可解析三处 `continue`，
无错误无计数；同一条坏行在 filter/stats 里是硬错误、在 convert 里计数、
在 pca/distance 里静默消失。截断的 WGS 文件产出数值错误的矩阵。
修前【据】：`r6_run_chunk_errors_on_sample_column_mismatch` 等 3 个测试
期望 Err，实测走 Ok 分支（"must fail the chunk, not be skipped"）。
修法：`malformed_line(line, reason) -> EngineError`（镜像 scanner 的
`IoError::MalformedRecord`，locus 取 CHROM:POS，POS 不可解析时按
`parse_pos1` 惯例置 0）；空行与 `#` 行仍跳过。顺带把死内存
`sites: Vec<(Vec<u8>,u64)>`（消费者只读 `.len()`，grep 证实）收敛为
`n_sites` 计数器，POS 解析随之消亡——即评审 backlog 第 67 条。

**R7（reed-stats，P0）投影 SFS 权重溢出 inf，位点计数但质量归零**
`crates/reed-stats/src/sfs.rs:111-149` — 超几何递推从左端 `k_lo` 以
`w=1` 起乘、最后归一化；未归一峰值 ~e^{n·H(p)}，n≳1024（p≈0.5）时
`sum=inf` → 返回空 `Vec` → `push_weighted` 仍 `n_sites += 1`、质量为零。
`--proj 1500` 在千样本上是现实用法，无上游上限、无告警。
修前【据】：`r7_projection_weights_do_not_overflow_at_large_n` 用
(10000, 5000, 1024)（峰值 ln w ≈ 765 > ln f64::MAX ≈ 709.8）与
(3000, 1500, 1500)：`projection_weights(10000,5000,1024) == Vec::new()`。
修法：众数锚定 `k* = ⌊(alt+1)(n+1)/(m+2)⌋`（u128 求值，m 大时 u64 乘积
溢出），`w=1` 向上用原比值、向下用因式化逆比值双向递推（离众数比值
≤1，任意 n 不溢出）；保留 i64 分母技巧、归一化与 1e-15 尾部保留。
小 n 数值不变：`r7_weights_match_exact_hypergeometric` 对 8 组小参数与
u128 精确二项 pmf 在 1e-12 内一致；既有手算锚点（1/6,4/6,1/6）逐位不变。

**R8（reed-reduce，P0）t-SNE 梯度乘未归一核，整体大 `sum_q≈O(n²)` 倍**
`crates/reed-reduce/src/extra.rs:141` — `mult = 4·(p_ij − q_ij) ·
q.get(i,j)`，末因子应为归一 `qij`（van der Maaten & Hinton 2008 两处都
是归一值）。平衡点不变但有限迭代轨迹不是 t-SNE，有效学习率随 n² 增长。
修前【据】：`r8_tsne_gradient_uses_normalized_q`（梯度提取为可测的
`tsne_grad` 后）实测 got/expected 比值 = 1.741254914 = 该用例的
`sum_q`，精确到全部输出位数。
修法：一个 token——`* qij`。全路径核查：`q` 只在梯度处被消费两处；
无 early-exaggeration，momentum/gain/lr 均为参考实现原值，无需补偿调整。
既有 `tsne_separates_blobs`（定性）修前修后都过——正是它区分不了
sum_q 倍的错误，新增单测才是真守卫。

**R9（reed-distance，P0）单态（fixed-ALT）列把均值写进标准化矩阵**
`crates/reed-distance/src/standard.rs:145-152` — `sd ≤ 1e-9` 分支写
`mu` 而非中心化 0。全 REF 列 mu=0 无感；全 ALT 列（p=1, mu=2，`--keep`
子集后常规出现）给每个 GRM 元素加 `4k/m`（k=1000/m=100k → 每对 +0.04）
并给 PCA 加秩一常数方向。plink2/GCTA 写 0；锁定测试用全 REF 列（均值
与 0 重合）看不到它。
修前【据】：`r9_fixed_alt_column_standardizes_to_zero_not_the_mean`
实测 `left: 2.0 right: 0.0`。
修法：退化分支写 `0.0`，文档注明与 flashpca 均值填充的有意分歧；
测试并断言含全 ALT 列的 GRM 与同位全 REF 列逐位相同。
波及：`reed-reduce/src/pca.rs::grm_diagonal_is_one_when_complete` 的旧
夹具**每列都退化**（`3i mod 3 = 0`），旧锚点 8/6 全是 R9 通胀——夹具
换为真正变化的 `(i+j)%3`（对角恰 8/6，断言带不变），测试内注明理由。

**R10（reed，P0）stream 求解器静默忽略 `--missing pairwise`**
`crates/reed/src/pca/mod.rs:283-330` — `--missing` 先验证、只有 dense
分支消费；`if chosen == "stream"` 一律 mean-impute，**auto 在大 N 内存
超限时必然走到这里**——恰是用户最关心的场景。模块文档 `:6-8` "both are
checked against plink2 in the e2e tests" 为不实陈述（无此测试）。
修前【据】：`r10_stream_solver_refuses_missing_pairwise` 与
`r10_auto_solver_refuses_missing_pairwise_when_memory_forces_stream`
（100 样本 + 60 KB 预算迫使 auto→stream）实测返回 `Ok`（solver=stream）。
修法：`solve()` 返回 `Result<Components>`；在 `chosen` 判定处拒绝
`stream + pairwise`（覆盖 auto），usage 错误给出三个出路（`--solver
dense` / `--missing meanimpute` / 提高内存）；删除不实文档句。
波及：`e2e_ordination::dapc_var_agrees_between_solvers` 旧锚定的正是
该静默降级（CLI 默认 `--missing pairwise` + `--solver stream`）——两处
调用改锚 `--missing meanimpute` 并注明（该测试比较的是 solver 间 `var`
归一一致性，与缺失策略无关）。

**R12（reed，P0）`convert --msa` 对重访 contig 游标倒退、位点重复发射**
`crates/reed/src/convert/writers.rs:292-381` — 帧低于 `start0` 被跳过，
但 `pos0 < cursor` 不被拒绝：`flush_msa_site` 钳 `from = cursor.min(pos0)`
不 panic、静默写出重复位点 + 倒退游标。乱序输入（chr1…,chr2…,chr1…）
其余命令全部拒绝，convert 静默产出错位 FASTA。
修前【据】：`r12_write_msa_refuses_regressing_frames`（pos0 5,10 后重访
3；20 bp 参考序列、1 样本）实测 Ok，写出 28 个碱基：
`>c1:1-20_s0\nAAAAAGAAAAGTAAAAAAAAAAAAAAAA\n`——pos 3 落在 pos 10 之后，
尾部 `seq[4..20]` 发射两次。
修法：进入 `pending` 判重之前，`pos0 < pending` 即报
`contig {} site {} follows {}: input is not contig-sorted`（沿用 stats
sink 的措辞）。语义边界用第二个测试钉死：接缝处同位点（blob1 尾 10、
blob2 头 10）保持 last-wins；已冲洗位点的重访（…5,10,12 后再来 10）
`10 < pending=12` 报错；blob 内升序由构造保证。

### 规模（P0/P1）

**R13（reed-distance，P0）估计量族串行 O(n²·m)、矩阵驻留 2–4 份、
`--threads/--memory` 解析后即弃** `matrices.rs:146-172, 255-287, 326-336,
402-411, 540-549, 699-708` + `translate/distance.rs:65-83` — n=2,000 ×
m=1M ≈ 34 GB 驻留 + 单核 ~2×10¹² 标量运算（小时级），`--bootstrap N`
再乘 N；`allele_ct_distance` 在对循环内重算 O(n·m) 的 per-sample dosage。
修前【据】：`grep -rn "par_iter|rayon" crates/reed-distance/` = 0 条
（结构性证明，性能 bug 无法以失败测试复现）。
修法（本轮范围）：全部估计量改 `(0..n).into_par_iter().map(row i)`、
行缓冲按索引序回收、**串行 scatter** 执行原有 `set()`/镜像写——每
(i,j) 条目的运算序列与串行完全一致，浮点不跨线程（无归约），任意池大小
逐位一致；`allele_ct_distance` 的 dosage 行提出对循环外（同操作数同
次序）；grm 死 `miss` 向量删除；`rayon.workspace = true` 入
reed-distance。守卫：`r13_estimator_results_thread_invariant`（含缺失
数据 + 退化列的矩阵，1 vs 4 线程池，grm/dist/king 的值与计数
`to_bits()` 相等）。未做（backlog）：driver 侧 `--threads` 接入
（scoped pool 一行）、sites 分块流式、n×m 冗余临时全消。

**R14（reed-core，P0）`RegionSet::stab` 每位点线性扫全部区间**
`crates/reed-core/src/grid.rs:252-263`（自带 `TODO(#1)`）— 外显子 BED
（~10⁵ 区间）× WGS（~10⁷ 位点）≈ 10¹² 次比较，`popg --bed-file` 实际
跑不完；构造时已建好的有序合并副本 `merged` 只有 `accessible_bp` 在用。
修法：合并时每个 span 记录成员输入下标（构造期排成输入序）；`stab`
对 `merged` `partition_point` 二分找唯一候选 span，再对**该 span 的成员**
套用原精确谓词过滤——注意 span 成员 ≠ 覆盖（嵌套区间只覆盖自己的子区间），
等价测试正是靠这一点抓出了"直接返回成员表"的错误设计。
守卫【据】：`r14_stab_matches_linear_reference`（旧线性实现原样作
oracle，对空/单/无序重叠/嵌套同界/相邻接缝/长跨多组输入 × 每区间
{start−1,start,end−1,end,end+1} 及 0/u64::MAX 穷举）；
`r14_stab_large_bed_fast`（200k 离散 span × 2000 次探测 < 10 s 上限，
旧实现 ~4×10⁸ 次比较，实测新代码毫秒级）。复杂度 O(log n + k)，
最坏（全部区间互相重叠的单 mega-span）退化为旧 O(n)——那是输出本身
就是 O(n) 的病态输入。

### 契约违背（P1）

**R11（reed，P1）`rename --set-id` 早返回静默禁用 `--update-name`**
`crates/reed/src/rename/mod.rs:1199-1207` — `new_id()` 先返回 `--set-id`
模板结果，从不查 `plan.update_name`；互斥块 `:609-630` 不拒绝该组合；
收尾的"拼写错误不得静默跳过"检查在 `set_id.is_some()` 时显式跳过
（`:746-753`）。
修前【据】：`r11_set_id_conflicts_with_update_name` 实测
`Ok(RenameReport { n_records: 3, n_ids: 3, … })`——exit 0，映射文件
从未被读过。
修法：互斥块新增一臂（沿用其消息风格）：`--set-id cannot be combined
with --update-name (the template rewrites every ID, so the mapping would
never be applied)`。
**不拒绝** `--set-missing-var-ids` + `--update-name`：核实为合法——早返回
只对 `.` ID 生效，map 对具名 ID 仍可达，防拼写检查也仍武装（跳过条件是
`set_id.is_none() && !clear_ids`）；假设映射键 `"."` 会撞上收尾检查而非
静默丢弃。`r11_set_missing_var_ids_keeps_update_name_reachable` 钉住
（rs1→new1 应用、`.`→v3 填充、拼错的键仍报错）。

---

## 2. 对评审结论的三处修正（以实测为准）

1. **R7 触发参数**：评审首选例 (m=30000, alt=9000, n=1500)（p=0.3）实测
   **不溢出**——左端起乘的 ln 增长只有 ≈548 < 709.8，该例只走到文档化的
   1e-15 尾部保留。真实触发区是评审括号里的推论：p≈0.5 且 n ≳ 1024
   （ln 增长 ≈ n·ln2）。复现测试按真实边界取参。
2. **R5 严重度上修**：四个子命令在实测里都是 **exit 0 报告成功**（小文件
   缓冲读先完成，截断不影响本跑），比"读中截断报错"更糟——操作者得不到
   任何信号；仅 filter 在大文件形态下被信号杀。
3. **R11 范围收窄**：`--set-missing-var-ids` 与 `--update-name` 的组合
   合法且有用（见上），评审建议一并拒绝是过度的。

另：R14 建议稿里"成员表即答案"的设计被等价测试证伪（嵌套区间），最终
实现是二分定 span + 成员精确过滤。

## 3. 核实为干净的部分（不再列为发现）

- **公式**：π/dxy/Hudson/WC84（pair 与 multi 手算复验，含 an=1 角落在
  backlog 的口径差异）、Watterson 分等位类、Tajima e1/e2、pixy
  variant-only 分母补全、plink `--het` F、Patterson f2/f3/f4 偏差校正、
  折叠/非折叠 SFS、Wigginton HWE 递推、UPGMA Lance-Williams/NJ Q 准则、
  Newick 转义（逗号/冒号/引号/空格/括号 `''`）。
- **确定性**：三处 rayon 归约全按索引序（LD_BLOCK=4096、per-contig、
  per-chunk）；chunk 边界只依赖 BGZF 块；Lanczos/bootstrap 种子为常量；
  e2e 在 {1,2,3,7,56} 线程字节一致；本轮 R13 延续该契约。
- **reed-engine**：票据窗口使重排溢出构造性不可能；三个 panic 面全被
  捕获并类型化；20 测试锁死接缝缝合与 reducer 死锁路径。
- **oracle**：`Cargo.toml` 零依赖，与 reed-core/reed-stats 无共享代码；
  期望值锚在外部且固定（scikit-allel 1.3.13 数值逐字内嵌 + Levene-Haldane
  枚举 + R 夹具），parity 测试是真证据。
- **错误链**：`reed::Error` 保源；`IoError` 带 `contig:pos`；58 处
  `map_err(|_|…)` 全部是用描述性文本替换无信息 payload（distance.rs
  三处 GraphError 例外，backlog）。
