# Changelog

All notable changes to reed are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the
lockstep workspace versions, and the CLI flags plus the output file formats
are the stability promise (`docs/spec/cli.md` §0).

## [Unreleased] — doc_v20 全库评审 must-fix 修复（2026-10-24）

11 个只读评审覆盖全部 crate（报告在仓库外 `../reed-rs-review.md`），14 个
must-fix 全部先写失败测试复现、再修复（测试名 `r1_*` … `r14_*`）。明细与
验证记录：`docs/rounds/v20/README.md`。

### Fixed

- **filter**：`-f/--apply-filters` 之前匹配整个 FILTER 列，`LowQual;LowMQ`
  被 `-f LowQual` 静默丢弃；现在按 `;` 拆分逐值匹配（bcftools 语义，R1）。
- **io**：空 FORMAT 列命中格式缓存哨兵 `(空, usize::MAX)` 导致 scanner
  越界 panic、整跑中止；现在按无 GT 布局处理（R2）。
- **index**：BCF 记录的 contig id 越界直接索引 `contig_names` panic（顺序
  路径无 catch_unwind）；现在报 job 错误（R3）。
- **cli**：`-O ''` 在 filter/rename/convert 三处 `split_at(1)` panic
  （exit 101）；现在报 usage 错误（R4）。
- **cli**：`-o` 与输入同路径时先读后截断，四个子命令（filter/rename/
  convert/index）实测 exit 0 并毁掉输入；现在拒绝并保持输入不动（R5）。
- **pca/distance**：畸形记录（列数不符、不可解析）被静默 `continue`，
  截断文件产出错误矩阵零提示；现在按 scanner 同款 `MalformedRecord` 报错，
  顺带把每位点死内存 `sites: Vec<(Vec<u8>,u64)>` 收成计数器（R6）。
- **stats**：投影 SFS 权重递推从左端锚定，`--proj ≳ 1024`（p≈0.5）时
  溢出为 inf 并静默丢弃位点质量；改为众数锚定双向递推，任意 n 不溢出，
  小 n 数值逐位不变（R7）。
- **reduce**：t-SNE 梯度用了未归一核（整体大 `sum_q≈O(n²)` 倍，输出不是
  t-SNE）；改为归一 `q_ij`（R8）。
- **distance**：单态（fixed-ALT）列把均值而非 0 写进标准化矩阵，给每个
  GRM 元素加 `4k/m` 并给 PCA 加秩一常数方向；改为中心化 0（plink2/GCTA
  约定，与 flashpca 的均值填充在 fixed-ALT 上有意分歧）（R9）。
- **pca**：stream 求解器静默忽略 `--missing pairwise`（auto 在大 N 时必然
  走到）；现在拒绝该组合并给出三个出路，模块文档中"已对 plink2 做 e2e
  校验"的不实陈述已删（R10）。
- **rename**：`--set-id` 早返回静默禁用 `--update-name`；现在互斥拒绝
  （`--set-missing-var-ids` 与 `--update-name` 本就兼容，保持可用）（R11）。
- **convert**：`--msa` 对重访 contig（乱序输入）游标倒退、位点重复发射，
  静默产出错位 FASTA；现在报 "input is not contig-sorted"，接缝处同位点
  last-wins 语义保留（R12）。
- **distance**：估计量族全程串行 O(n²·m) 且矩阵驻留 2–4 份、
  `--threads/--memory` 被忽略；现在按行 rayon 并行 + dosage 行预计算，
  任意线程数结果逐位一致（1 vs 4 线程 `to_bits()` 相等锁定）。流式分块
  与 `--threads` 全链路接入留作 backlog（R13）。
- **core**：`RegionSet::stab` 每位点线性扫描全部 BED 区间（自带 TODO），
  外显子 BED × WGS 实际跑不完；改为合并区间二分 + 成员精确过滤，
  等价性由旧线性实现做 oracle 的穷举测试锁定，200k 区间 × 2000 查询
  从分钟级降到毫秒级（R14）。

### Fixed（第二批：should-fix 中的 bug 类 16 条，同轮追加）

- **filter**：cohort 变化后（`--keep`/`--mind`/掩码），表达式
  AC/AN/AF/MAC/MAF 仍读全群体陈旧 INFO，与同 run 的 `--maf` 口径分裂；
  现按 `info_trusted` 门控从基因型重算（S58）。
- **reduce**：UMAP 学习率去掉 1/4 下限，照 umap-learn 线性衰减到 0（S28）；
  t-SNE 困惑度搜索非收敛行不再静默清零，`perplexity <1 或 ≥ n−1` 报
  类型化错误（`tsne` 变 `Result`）（S29）；三处 kNN 排序改
  `f64::total_cmp`，NaN 基因型不再 panic（S31）；`xval` 组 id 越界
  返回类型化 LDA 错误而非 panic（S32）。
- **graph/cli**：bootstrap 根节点不再恒标 100% 支持（跳过全 1 键，
  support 数量对齐 ape 的 n−2/n−3）（S36）；`--tree --bootstrap` 的
  discordance 重复矩阵改用用户的 `--mismatch-count/--error-probability`
  （修前恒默认参——Count 树配 Probability(40) 的 support）（S37）；
  同批把 `--threads` 经 scoped pool 接通 distance 估计量与 bootstrap
  （R13 收尾，1 vs 4 线程输出逐字节一致锁定）。
- **pca**：剂量矩阵分配前对照 `--memory`，超限报 usage 并指向
  `reed filter --indep-pairwise`（修前 5000×3M 直接 OOM）；`run_pca_opened`
  收敛进 `read_genotype_matrix`（S76）。
- **filter**：keep-mask 长度契约统一——record_pipeline 五处短掩码容忍
  改直接索引 + debug 断言，与 sketch 位平面路径口径一致（S68）。
- **stats**：WC84 单群 an∈{0,1}（如单个半调用）时 pair 与 multi 两条
  路径统一走 keep-c 分支（修前 pair 整位点丢弃，`--fst wc84` 与
  `--fst-global` 在该角落分歧；评审"multi 保留"的说法在 r=2 也不成立，
  一并修正）（S46）。
- **io**：FORMAT 整数越界由钳制改 missing 哨兵（对齐 INFO 路径与
  htslib）（S39）；损坏索引虚拟偏移越界报 `MalformedBgzf` 而非切片
  panic（S40）；POS 超出 BCF i32 范围与负 pos0 报 `MalformedRecord`
  而非静默钳制（S43）；符号等位 + INFO/END 的索引跨度改用 END（三处
  span 站点收敛为共享 `indexed_span_len`，区域查询不再漏跨 bin 的
  符号缺失）（S42）。
- **oracle**：PBS 模块文档改为与代码一致的 per-pair 规则（"三群交集"
  的说法不实；代码零改动，契约测试钉死）（S79）。
### Fixed（第三批：对前两批的复审修复）

- **reduce**：**撤回 R8 的梯度"修正"**——末因子应为原始核 `w` 而非归一
  `q`（对 KL 代价做中心差分验证：原始核误差 1e-10，归一核 0.58），前一版
  让梯度小了 `sum_q≈n²` 倍、几百样本时嵌入停在初始化量级；测试改为对代价
  函数做有限差分（旧测试是梯度自身的复写，恒真）。同批：困惑度搜索的
  beta 从 `1/mean(d²)` 起步（真实剂量距离下 `exp(-d²)` 下溢导致整行清零）。
- **io**：符号跨度改为 htslib 的 `get_rlen` 规则——INFO/END 与 ALT 形状
  无关（`ALT=.` 的 gVCF 块、多等位同样延伸），并补上默认单遍写入器与
  `rename -W` 两条漏掉的路径（此前只有 sketch 写入器正确）；BCF 的
  `pos0 = −1` 恢复为 htslib 的端粒哨兵（渲染 POS 0、原样传递），只有
  < −1 报错；FORMAT 越界警告不再自称 "INFO"。
- **graph**：**撤回 S36**——本机 ape 5.8.1 实测 `prop.clades`/`write.tree`
  保留根节点（Nnode 个值），前一版丢弃根并声称 ape 也丢弃是事实错误；
  支持向量恢复为每个内部节点一项，根恒 100%。
- **pca/filter**：`--memory` 上限移入逐块收集（超限即中止扫描，不再先攒
  满）；`stream + pairwise` 与 `--set-id`×`--update-name` 的拒绝移到读输入
  之前；`INFO/AC` 等带前缀的内置字段在 `--set-gts` 掩码表达式中编译期
  拒绝（此前静默求值为 Missing）。
- **filter/distance/convert**：`-f .` 只匹配整列为 `.`（对齐 bcftools 的
  filter id −1 语义）；dosage 提升缓存由 f64 改 i8（省 16 GB @ 2000×1M，
  数值逐位不变）；MSA 乱序检查覆盖窗口外帧；keep-mask 断言在 release 也生效。
- **安全**：stdin 落盘临时文件改 `create_new`（O_EXCL）+ 有界重试——
  修前共享主机上可预测路径可被预置文件/符号链接劫持并截断受害者文件
  （实测复现）（SEC-1）。

## [Unreleased] — doc_v18 随机数基础设施（2026-10-21）

`xvalDapc` 与 `aboot` 的共同前置，加上两个消费者本身。没有可分割的 RNG，
任何并行重抽样的结果都会依赖线程调度，违反 `architecture.md` §7.6 的
确定性契约。

### Added

- **`reed-core::rng::Rng`**：PCG32（O'Neill 2014），可分割。
  - `from_seed` / `from_seed_stream` / **`derive(index)`**：按**重复序号**
    派生独立子流，因此第 `i` 次重复在 1 线程与 64 线程下取到同一批随机数。
  - `below` 用**拒绝采样**而非 `% n`：模偏在置换检验里是系统性误差。
  - `shuffle`（Fisher–Yates）、`sample_n`（Floyd 算法，代价
    `O(min(k, n-k))`，`poppr::genotype_curve` 用的就是它）、`resample_n`
    （有放回）、`multinomial`（`poppr` 参数化 bootstrap 用的；负权重直接
    拒绝，否则会静默偏置）。
  - 与独立的 PCG32 参考实现逐位对拍，参考序列锁进测试。

### Changed

- `reed-reduce` 里两处私有 xorshift（`umap.rs`、`extra.rs`）迁到共享的
  `Rng`。UMAP/t-SNE/NMF/PPCA 的随机序列因此改变；它们的测试只断言结构性
  质（分离度、对称性、重构误差），没有锁死数值，所以全部通过。

### Added — `xvalDapc`（RNG 的第一个消费者）

- **`reed-core::math::r_compat`**：R 的 `round()`（银行家舍入）与
  `pretty()` 步长表。`xvalDapc` 的候选 PC 列表与训练集大小都依赖这两者，
  而 Rust 的 `f64::round` 是 half-away-from-zero，不等价。
- **`reed-reduce::xval`**：按群体分层抽样、每档 `n_rep` 次重复、以
  `sqrt(mean((success-1)^2))` 选最优（照抄 adegenet 的非标准 RMSE），
  并列取最大 `n_pca`。PCA 在全数据上算一次、不随重复重算——重算会把验证
  样本泄漏进主成分。
  - 每组训练集大小是 `round(training.set * n_g)`，**不夹取到 `n_g-1`**，
    与参照一致（该情形在 `training.set` 被压到 `(popmin-1)/popmin` 后不可达，
    但去掉夹取让算术与参照逐项相同）。
  - 随机基线是 **`1/ng`**，即 `adegenet` 标签置换基线的期望值；不是
    `Σ p_g²`（后者是按先验猜的准确率，在不等组数据上高估 17%）。
- **`reed pca --dapc --xval`**：`--xval-rep`、`--xval-training`、
  `--xval-result`、`--xval-seed`、`--xval-pca`；写 `PREFIX.dapc.xval`。

### Added — `aboot`（bootstrap 树支持率）

- **`reed-graph::bootstrap`**：按**位点**（列）有放回重采样、每次重算距离
  与树、统计每个分支在重复中出现的比例。有根树按 clade 计，无根树按
  bipartition 计（规范化到含 tip 0 的一侧，避免同一划分记两次）。
- **`reed distance --tree upgma|nj --bootstrap N`**：`--bootstrap-seed`、
  `--bootstrap-cutoff`；支持率写进 `PREFIX.nwk` 的节点标签。
- 与 `ape::prop.clades` 的确定性对拍：固定三个重采样索引，计数
  `3 3 3 2 1` 逐项一致。

### Verified

| 检查 | 结果 |
|---|---|
| PCG32 参考序列（seed=42，前 5 个 u32） | 与独立实现逐位一致 |
| bootstrap clade 计数 vs `ape::prop.clades` | 逐项一致（3/3/3/2/1） |
| 候选 n.pca 列表 vs `pretty()` | 7 个范围逐项一致 |
| `round(0.9*n)` 银行家舍入 vs R | 一致 |
| 弱信号夹具成功率 vs adegenet | 同为"PC 越多越差"，量级一致 |
| 串行 vs 4 线程分片（256 次重复） | 逐元素相同 |
| 重复 `i` 是否受其它重复影响 | 不受（烧掉 1000 次抽样后仍相同） |
| `below` 均匀性（70k 次，7 桶） | 每桶偏差 < 5% |
| `sample_n` 位置均匀性（30k 次，10 选 3） | 每位偏差 < 5% |
| `sample_n` 抽样预算（20 选 5） | 恰好 5 次，不是 20 次 |
| `derive` 流间相关（2 万条流） | 相邻流 r = -0.0005 |
| `multinomial` 权重跟随（50k 次） | 每类偏差 < 5% |

## [Unreleased] — doc_v17b 判别分析与距离树（2026-10-21）

面向自然群体的一轮：收**有监督判别**（DAPC）与**由距离矩阵建树**
（UPGMA/NJ/newick）两族。

### Added

- **`reed-reduce::lda`**：Fisher 线性判别分析，照 `MASS:::lda.default` 的
  `method="moment"` 分支（`adegenet::dapc` 调的就是它）。
  - `fit(x, groups, levels, tol)` / `Lda::project(x, dimen)` →
    `LdaProjection { coords, posterior, assign }`。
  - 单边 Jacobi SVD（工作区原本没有 SVD）；不形成 `A^T A`，精度更好。
  - 与 `MASS::lda` 对拍：`means` 2.2e-15、`posterior` 4.9e-16、`assign`
    逐样本一致；`scaling` 按 `|dot|` 归符号后 1.0e-15。
- **`crates/reed-graph`**：新 crate，只依赖 `reed-core`（L3 纪律）。
  - `Tree`（ape `phylo` 形状）、`to_newick`（`ape::write.tree` 格式，
    显式栈，深度上限 10 000）。
  - `upgma`：`hclust(method="average")` 的 Lance–Williams 递推，含并列
    扫描序与存活槽（取较新簇）规则——两者都是逐位一致的必要条件。
  - `nj`：`ape::nj` 的 `C_nj` 语义（`Q=(n−2)d−S_i−S_j`，边数 `2n−3`，
    负枝长原样保留）。
  - 与 `hclust` 对拍：合并高度 4.9e-13；与 `ape::nj` 对拍：RF=0、
    13 条枝长逐位相同。
- **`reed pca --dapc --pop FILE [--n-da N]`**：PCA → LDA，写
  `PREFIX.dapc.ind`（`#IID POP ASSIGN LD1..LDk`）、`PREFIX.dapc.eig`、
  `PREFIX.dapc.var`。`--n-da` 默认 `#{svd² > 1e-10}`，照 adegenet 的
  `lda.dim`。assign 与 `adegenet::dapc` 逐样本一致。
- **`reed distance --tree upgma|nj [--digits N]`**：由已算好的矩阵建树，
  写 `PREFIX.nwk`（不写矩阵文件）。
  - `--tree` 只接受**距离**矩阵：对角线非零的相似度/亲缘矩阵（`grm` 默认值、
    `cov`、`king`、`pihat`、`phi`、`ibs`）被拒绝并提示改用 `1-ibs`。
  - `--digits` 是 **R `format()` 语义的有效数字**（默认 10），不是小数位；
    `0` 被拒绝（会把所有枝长舍成 0）。
  - UPGMA 的枝长按 `ape::as.phylo.hclust` 的口径折半：`hclust$height` 是
    合并距离（节点高度的 2 倍），不折半会让每条枝长变成正确值的 2 倍。
    超度量输入下树的 patristic 距离与输入矩阵逐对相等，有性质测试锁死。
- 门面：`reed::run_dapc` / `run_dapc_opened` / `DapcRequest` /
  `DapcReport`；`reed::Mat` 转出（`reed-graph` 与 CLI 共用）。

### Notes

- `docs/rounds/v17b/`（README / 01-research / 02-spec / 06-roadmap /
  STATUS）。核心决定：DAPC 的 PCA 沿用 reed 口径而非 `ade4::dudi.pca`
  的 `li = U·d`；两者每轴差一个常数因子，而 LDA 对特征空间整体缩放不变，
  故 `posterior`/`assign` 一致、`ind.coord` 仅尺度不同（01 §1.5 有推导）。
- 明确不做：`find.clusters`、`snapclust`、`dist.genpop` 的五种距离、
  全部绘图层、Bruvo/SSR/GenAlEx。（`xvalDapc` 与 `aboot` 后来在 doc_v18
  完成。）

## [Unreleased] — doc_v17 个体间距离矩阵族（2026-10-20）

新增 `reed-distance` crate 与 `reed distance` 命令；PCA 的 GRM 构建迁入
该 crate，成为所有矩阵的统一底座。

### Added

- **`crates/reed-distance`**：10 种个体×个体矩阵，每种复刻一个参照工具。
  - `allele_ct_distance` / `ibs` / `one_minus_ibs`（plink `--distance`）
  - `grm` / `covariance`（plink2 `--make-rel`、GCTA）
  - `king`（plink2 `--make-king`，KING-robust）
  - `pi_hat`（plink `--genome`）
  - `phi`（vcftools `--relatedness2`，Manichaikul）
  - `ajk`（vcftools `--relatedness`，Yang）
  - `discordance`（bcftools `gtcheck`）
  - 4 种输出格式：`triangle`（含对角线）、`square`、`square0`、`bin`
- **`reed distance`**：`--type`、`--format`、`--stand`、`--mean-impute`、
  `--unbounded`/`--nudge`、`--mismatch-count`/`--error-probability`。
- **`reed::read_genotype_matrix`**：VCF→剂量矩阵的共享桥接，PCA 与
  distance 共用。
- **`docs/rounds/v17/`**：三子代理的源码调研台账（file:line + 公式）、
  规范、路线、状态。

### Changed

- `Mat` 从 `reed-reduce::linalg` 下沉到 `reed-core::linalg`，打破
  `reed-reduce` ↔ `reed-distance` 的循环依赖。
- `Grm`/`Stand`/`Standardized`/`standardize` 迁入 `reed-distance`；
  `reed-reduce::pca` 改为 re-export，v16 API 不变（PCA 数值不变：
  HapMap3 λ1 = 26.200944934）。

### Verified

| 矩阵 | 参照 | 结果 |
|---|---|---|
| `grm` | plink `--make-rel` | max abs 5e-6 |
| `king` | plink2 `--make-king` | max abs 5e-7 |
| `allele-ct` | plink `--distance` | max rel 5e-6 |
| `phi` | vcftools `--relatedness2` | 915849 对，计数零差异 |
| `discordance`（count） | bcftools `gtcheck -E 0` | 逐对相同 |

### Known issues

- **`pihat` 未对齐 plink**：plink 用 PLINK 1.07 的 Nei 偏差修正，reed 用
  理论期望值。输出保留但标记未对齐。
- `discordance` 的 HWE 列、`--genome` 的 DST/PPC/RATIO 未实现。
- `ajk` 与 `cov` 未对拍。

## [Unreleased] — doc_v16 CLI 接入：`reed pca`（2026-10-20）

把 `reed-reduce` 接进命令面，新增 `reed pca`。

### Added

- **`reed pca in.vcf.gz --pcs 20 -o out`**：样本关系矩阵（GRM）的主成分。
  - `--stand binom2|binom`（默认 binom2）、`--missing pairwise|meanimpute`
    （默认 pairwise，对齐 plink2/GCTA）、`--solver auto|dense|stream`
    （auto 按内存预算选）、`-o PREFIX` 写 plink 兼容的
    `PREFIX.eigenvec`/`.eigenval`/`.pve`。
  - `--samples/-S/--keep/--remove` 样本子集；`--threads`/`--memory` 复用
    既有预算机制。
  - 门面新增 `reed::run_pca` / `reed::run_pca_opened` / `PcaRequest` /
    `PcaReport`（`run_pca_opened` 让 CLI 只开一次输入，stdin 可流式）。
- **`crates/reed-cli/tests/e2e_pca.rs`**：6 个用例——文件与表头形状、
  binom/binom2 精确 2 倍、无缺失时两种策略一致、样本子集、5 个用法错
  （exit 2）、dense 与 stream 一致（rel < 1e-6）。

### Verified

HapMap3（957×14266，k=10）与 plink2 逐位一致：

| stand + missing | reed λ1 | plink2 |
|---|---|---|
| binom2 + meanimpute | 26.200944934 | 26.2009 |
| binom2 + pairwise | 26.280628904 | 26.2806 |
| binom + meanimpute | 52.401889869 | = binom2 × 2 |

`cargo test --workspace`：59 个测试二进制全绿。

### Not included

`--missing em`、`--use-gl`、`sd/center/none` 标准化、`--window`（Local
PCA）、`--tw`（Tracy–Widom）、`--indep-pairwise`（上游预处理）。理由见
`docs/rounds/v16/06-roadmap.md` §P2 与 `docs/spec/cli.md` §11b。

## [Unreleased] — doc_v16 降维族（2026-10-20）

新增 crate `reed-reduce`（零外部依赖，3119 行，29 个单元测试），并完成
PCA 求解器与降维族的源码级调研。**未接入 CLI**；命令面接入见
`docs/rounds/v16/06-roadmap.md` P2。

### Added

- **`crates/reed-reduce`**：零依赖的降维原型库。
  - `pca`：binom2 标准化、GRM（mean-impute / pairwise）、稠密 PCA、
    scores/loadings。HapMap3（957×14266）上 λ1/2/3 =
    26.200944934 / 23.394076211 / 6.932494547，与 plink2 exact 逐位一致。
  - `lanczos`：矩阵-free 块子空间迭代（top-k），与稠密解特征值一致。
  - `graph`：PCoA、Kernel PCA（linear/RBF/IBS）、Laplacian eigenmaps、
    Diffusion map、Isomap。
  - `grm_reduce`：上述五种方法的 GRM 驱动入口（缺失策略在 GRM 构建时
    确定一次，全部方法继承）。
  - `umap`：fuzzy simplicial set + SGD，按 umap-learn 0.5.12 逐行实现；
    `find_ab_params(1.0, 0.1)` = (1.576943618, 0.895060710) 与参考实现
    差 < 1e-7。
  - `extra`：t-SNE、NMF、LLE、PPCA-EM（缺失值原生）。
  - `eig`：稠密对称特征分解（EISPACK `tred2` + `tql2`），以及**部分谱
    求解器**（Householder 三对角化 + Sturm 二分 + 逆迭代），只提取前 k
    个特征对。N=4000 时比稠密路径快 2×（170.6 s vs 334 s），λ 逐位一致，
    残差 < 1e-13。
- **`docs/rounds/v16/`**：README、01-research（12 条判决 + 产物索引）、
  02-spec（API 与精度契约）、06-roadmap（P1–P4）、STATUS。
- **`docs/pca-algorithms.md` §8–§15**：flashpca `stand` 五选项、跨工具
  对拍、GRM 角色、求解器族谱、速度/精度实测、纯 Rust 实现路线。

### Verified

- 与 numpy `eigh` 相对误差 4.5e-16（机器精度）。
- 与 plink2 / GCTA / plink 1.9 / flashpca 在同集合同缺失策略下数值一致
  （2e-6 = 文本位数极限）。
- 与 sklearn 对拍：PCA 2.91 vs 2.91、KernelPCA 1.80 vs 1.80、
  Isomap 1.93 vs 1.95、NMF 重构误差 0.1907 vs 0.1905。

### Known limitations

- 块子空间迭代在 HapMap3 上 `max_rel_residual ≈ 1e-2`（特征值已收敛，
  Ritz 向量残差未达生产级）；稠密路径是精确参照。
- 稠密路径 O(N³) 单线程（4000×20000 需 334 s），因零依赖约束无 BLAS。
- UMAP 用精确 kNN、t-SNE 无 Barnes-Hut、Isomap 用 Floyd-Warshall。

## [Unreleased] — doc_v15 缺陷修复（2026-10-18）

对一份全库评审（6 子代理静态评审 + 逐条实证复核，40 条发现）的修复；
规范与证据见 `docs/rounds/v15/`。

### Fixed

- **filter**：Float 类型 FORMAT 字段的基因型掩码不再整体失效（D15-1）；
  掩码与 `--set-gts .` 改变 cohort 后 `-q/-Q` 一律按掩码后基因型重算，
  不再读陈旧 INFO/AC,AN（D15-2，cli.md §3.2）；`--proj` 超过群体染色体
  数、pop 文件把同一样本分配进两个群体、`--levels` 空列表均报用法错并
  点名（D15-9/10/21）；`--bed-file` 与 `--step` 解析期互斥、`-o -` 写
  stdout 而非字面文件 `-`（D15-22/39）。
- **popg**：缺少 `-o` 在打开输入前即报错（D15-13）；stdout 管道关闭由
  panic 改为普通运行错误（D15-14）；裸 popg 不再忽略
  `--variant-only/--bypass-invariant-check`，不变位点检查按窗口路由同款
  执行（D15-15）；`--f2/--ld-decay` 不再静默把 π/dxy/θ 切到 variant-only
  口径（D15-18）；窗口发射改为"严格更后位点才 emit"，跨 BGZF 块接缝的**合法同 POS 重复**（VCFv4.3："It is permitted to have multiple records with the same POS"）落入仍开窗口，窗口和与 contig 行逐项一致（D15-19）；
  尾随三态记录之后的窗口正常输出（D15-20）；`--pbs` **可重复**：单三元组写 `PREFIX_pbs.txt`（逐字节不变），多三元组每份
  `PREFIX_pbs.<A>-<B>-<C>.txt` 各持自身表头（D15-8，命名规则见 cli.md §7）。
- **convert**：`--phylip` 组合运行不再吞掉 `--012/--hapmap/--msa` 的位点
  （D15-4）；`--012` 双等位门修正，三态位点按 vcftools 语义跳过（D15-5）；
  仅 `--nexus-binary` 时 `--outgroup` 生效且未知名报错（D15-6）；MSA 长
  等位不再以 u8 截断（D15-26）；`POS=0` 按 malformed 跳过（D15-27）；
  gzip hapmap 可走 stdin（D15-29）。
- **rename**：`--rename-chrs` 撞名报错而非写出重复 `##contig` 的非法
  VCF（D15-25）；`--idspace-to` 与 `--set-missing-var-ids` 互斥（D15-30）；
  `--rm-dup-ids` 改为两阶段（串行定序 + 批次并行 deflate：worker 不再做
  被丢弃的压缩、索引键随记录并行算好，有序 sink 按字节封口批次并行压缩）：
  3M 记录 `-Oz` 1.39s → 0.23s（对照 15.4× → 2.5×），压缩流与 sidecar 对
  任意 `--threads` 逐字节相同，decompressed 输出与改前逐字节一致（D15-28）。
- **engine**：输入源线程的 panic 与 worker/reducer 一样收敛到单一错误
  出口，报 `reader panicked: …`（`EngineError::ReaderPanicked`；此前该闭包
  不在 `catch_unwind` 下，panic 会绕过引擎的错误通道，D15-40）。
- **io**：BCF 头恰落 BGZF 块边界（头文本长度 ≡ 0 mod 65280）不再整体
  错位一字节（D15-3）；`##contig` 行的显式 `IDX=` 按 htslib 规则生效
  （D15-24）。
- **stats**：Hudson Fst/PBS 在窗口内无可用分子时输出 0.0 而非 NA
  （scikit-allel `nansum` 口径，D15-11）；oracle 的 PBS 与 f2/f3/f4 改为
  与 shipped/allel 同口径（D15-12/31）。
- **性能**：`--indep-pairwise` 的保留位点查询由逐 contig 线性扫描改为
  二分（1 万 contig、30 万位点单线程 ≈11×，D15-17）；`--rm-dup-ids` 的
  去重输出由"sink 串行重压"改为两阶段并行（3M 记录 `-Oz` 1.39s →
  0.23s，D15-28）；GT 浮点渲染 `format_g6_into` 改为栈缓冲零分配
  （242.9 → 182.9 ns/次，D15-35）；`sketch_size_bytes`
  补上 mask 位平面行，`--mind`+`--indep-pairwise` 的 sketch 分支加预算
  检查（D15-16，预防性）。

### Changed

- `reed popg` 不再接受 `--tmp-dir`（该命令不落临时文件）；`reed stats`
  保留该参数并在帮助中写明"单趟报告从不落盘，此参数仅为兼容接受"
  （doc_v15 D15-23，cli.md §2 同步）。
- 内部收敛：Patterson 只保留校正式（`PattersonKind::Naive` 与相关死结构删除，
  `f2_site`/`f3_site` 去掉 kind 参数）；QC summary 只剩一个实现
  （`summary::SummaryAcc` 再导出 `report` 的活结构，`QcReportAcc::merge` 委托）；
  fst/PBS 审计改用 `PairAcc::add_partial`，不再假装跟踪 missing（D15-38）；
  `notify_abort` 持锁通知关闭丢失唤醒窗口（D15-32）；`--indep-pairwise` 的
  保留位点查询与 convert/popg 的 `--proj` 校验收口到共享层。

### Changed

- `convert` 的输出缓冲改为两级落盘（doc_v15 D15-7）：单缓冲超 4 MiB 或
  sink 缓冲总量超 8 MiB 即落临时文件、渲染时流回；覆盖 hapmap 行、
  `--012` 矩阵与 `.pos`、phylo 列与 `used_sites`、MSA 帧流（渲染改为逐样本遍历，
  并按样本 seek 只读该样本切片——不做整帧重读，延续 v13 D13-2 的禁止重复 IO）。
  800k 位点/50 样本夹具
  单线程峰值内存：hapmap 123.4 → 15.6 MB、`--012` 70.8 → 26.0 MB、MSA 151.4 → 15.8 MB，
  输出逐字节不变；临时文件在成功、失败与中止路径均清理。MSA 的
  per-contig 缓冲仍驻留内存（留待单独设计）。

### Removed

- 死代码：`JointSfsAcc`（连同 `tests/sfs_joint.rs`）、`PattersonKind::Naive`、
  `Cursor::skip`、`SiteView::info_int/info_float`、`HeaderScan::n_samples`、
  `WindowGrid::bin_of_pos0`、`keep_mode` 死参、`header_level` 弃值、stats sink 的
  死 `or(take)`、convert 的空臂（doc_v15 D15-36）。

## [0.4.0] - 2026-10-17

The doc_v14 command-surface round: `reed stats` is split into the QC command
(`stats`) and the windowed statistics command (`popg`); `--hist` and
`reed rename` are new; `reed --help` groups its commands seqkit-style. The
split is byte-identical for every existing output
(`gold/popgen/run_split_parity.sh`, 22 cases).

### Changed

- **`reed stats` is now the QC-report command** (`--freq/--hardy/--het/
  --missing/--depth/--summary`); the windowed population-genetics statistics
  moved to **`reed popg`** (alias `popgen`; doc_v14 §2). The moved flags
  (`--pi --dxy --theta --tajima-d --fst --fst-global --pbs --pbs-normed
  --f2 --f3 --f4 --dstat --f4-ratio --block-size --ld-decay --ld-pairs
  --sfs --proj --sfs-int --variant-only --bypass-invariant-check
  --contig-lengths --window --step --levels --skip-empty --bed-file`) on
  `reed stats` are parse-time errors naming the exact rewrite (doc_v14
  D14-2) — never silent aliases.
- A bare `reed stats` (no report flag) is a usage error naming the flags.
- The library facade renamed its drivers to follow the commands:
  `run_stats` → `run_popgen` (`StatsRequest` → `PopgenRequest`,
  `StatsReport` → `PopgenReport`), `run_reports` → `run_stats`
  (`ReportRequest` → `StatsReportRequest`, `ReportsOutcome` →
  `StatsOutcome`). The library has no stability promise (cli.md §0.9).
- `reed --help` lists the commands in groups (Filtering, Statistics,
  Transform, File Operations, Miscellaneous); command paths are unchanged
  (doc_v14 §1).

### Added

- Expressions (`-i/-e`) accept negative number literals (`MQRankSum<-12.5`,
  the GATK hard-filter shape; bcftools' unary minus, `filter.c:3902`). A
  `-` before a digit starts the number; a leading `+` stays an error, as in
  bcftools, which lexes it as the addition operator (`filter.c:332`); `-`
  inside identifiers stays part of the name, and there is no arithmetic, so
  `A - B` remains an error.
- `reed popg`: the windowed statistics (π, dxy, Fst, θ, Tajima's D, PBS,
  SFS, Patterson f-statistics, LD decay), with the 0.3.x flags and
  byte-identical outputs.
- `reed stats --hist LIST`: histograms from the same pass (doc_v14 §3) —
  built-ins `qual`, `dp`, `dp-site`, `gq`, `idl`, `af`, plus generic
  `INFO/<TAG>` and `FMT/<TAG>` requests; one `PREFIX.<m>.hist` TSV per
  metric with exact mean/median/p95 for integer metrics.
- `reed stats --singletons`: the per-allele singleton/private-doubleton
  table in vcftools' layout (`CHROM POS SINGLETON/DOUBLETON ALLELE INDV`,
  `PREFIX.singletons`), with a POP column under `--pop`.
- `reed stats --pop FILE` (with `--force-samples`): stratified QC
  (doc_v14 §2.2) — `--freq`/`--missing` write `PREFIX.afreq.pop` and
  `PREFIX.vmiss.pop` (one row per population, first-appearance order),
  `--singletons` gains a POP column (`.` for a carrier outside every
  group; counting still runs over the whole cohort). The stratified values
  equal a per-population plink2 `--freq/--missing` run
  (`gold/popgen/run_popqc.sh`); the cohort files are unchanged. `--pop`
  without one of `--freq/--missing/--singletons` is a usage error.
  Multiallelic ALT frequencies stay collapsed (the `.afreq` convention).
- `reed rename` (alias `reheader`): identifier edits without decoding
  genotypes (doc_v14 §4) — samples (`-s/-N/--replace-sample-name`), contigs
  (`--rename-chrs`, `--chr-style plain|ucsc|mt-is-m`), variant IDs
  (`--set-id`, `--set-missing-var-ids`, `--update-name`, `--idspace-to`,
  `--rm-dup-ids`, `--clear-ids`), header lines (`-h -f --add-meta
  --remove-meta`); VCF byte editing with `-Ov/-Oz`.

### Fixed

- The `-q/-Q` stale-INFO warning fired on every run whose header declares
  `INFO/AC` and `INFO/AN` with an unchanged cohort (every GATK header does),
  asserting "-q/-Q reads INFO/AC,INFO/AN" of commands that never passed
  `-q/-Q`. It now fires only when a `nref`-form threshold actually takes
  the INFO path (spec cli.md §3.2: 取到时).
- `reed rename` wrote a blank line at every chunk boundary (the seam record
  arrives with its terminator and got a second one), so every output over
  one chunk carried stray empty lines; with `--rm-dup-ids` the same doubled
  newline made the dedup pass fail outright (`job failed: dedup: unparsable
  record`). The suffix probe no longer restarts at `_2` for every duplicate,
  so a run of identical IDs is linear instead of quadratic.
- `reed rename` ran its record planes on inputs whose CHROM is not declared
  in `##contig`: a contig-less VCF (legal VCF) and a sites-only 8-column VCF
  both failed mid-run. Plain output now passes them through, `-W` says what
  the header is missing, and the record planes are refused up front on
  sites-only input.
- `reed rename -h FILE -W` built the index from the *input* header's contig
  order while writing the replacement header, so region queries could
  resolve to the wrong contig; the index now follows the written header.
- `reed rename -f FILE` appended the missing `##contig` lines after the
  `#CHROM` line (not a valid VCF header); they are inserted before it.
- `--rename-chrs` and `--chr-style` overwrote each other instead of
  stacking (the style was applied to the pre-rename name and discarded the
  explicit rename).
- A CRLF header's last sample kept its `\r`, so sample renames counted the
  edit but silently did not apply it.
- `--replace-sample-name` used to accept regex syntax and silently do
  nothing (or strip the wrong text); the flag is documented as literal, so
  metacharacters are now a usage error.
- `--update-name` did not check that a mapped name exists (the one
  validation-table row that was missing): a typo renamed nothing and exited
  0; the run now fails naming the key.
- `--no-command` was parsed and ignored by `reed rename` (the
  `##reed_command` line was written anyway).
- `--hist` float tags were labelled one 0.1 step high (`iqual` divided by
  10 instead of `(iqual−1)/10`), so `SC=1.0` printed as the bin `1.1 1.2`;
  a missing QUAL became a real `0.0` bin instead of `#missing`; the `af`
  bins used the wrong denominator and dropped indels; integer metrics
  passed the overflow bucket off as an exact value (`DP 600` → row `501`,
  mean 501). Bins are now labelled in the quantity's own units, capped
  observations print `>max` and stay out of the exact quantiles, and only
  occupied bins are written (a single huge value used to materialize a
  half-million-row file).
- `--hist gq`/`dp`/`dp-site` did not check that the header declares the
  field and wrote an all-zero histogram with exit 0; an all-empty
  `--hist ""` list did the same. Both are usage errors now.
- `--hist FMT/<TAG>` counted one missing observation per record when the
  FORMAT key was absent (it is one per sample) and never split
  Number=A/G/R subfields, so `FMT/AD` looked empty.
- `StatsOutcome::hist` (the library's histogram output) was always empty.
- `reed popg --ld-pairs` without `--ld-decay` exited 0 and produced no
  output; it is a usage error now.
- INFO values survive the BCF round trip: an Integer outside the 32-bit
  range becomes missing instead of wrapping (`3000000000` used to be
  written as `-1294967296`; htslib and VCFv4.5 §1.3 both forbid the value),
  a Flag carrying a value keeps it (`DB=1` lost it, and `DB=0` — bcftools'
  "clear the flag" — came back as a set flag), and a String value of `.` is
  no longer turned into an empty string (doc_v12 D12-3/4/5).
- `-q/-Q` follow the genotypes once the sample set changed: after a static
  subset, a mask or `--mind`, the input's `INFO/AC,AN` describe a different
  cohort, so the decision is recomputed from the kept genotypes like
  bcftools does (vcfview.c:341-350); the promised stale-INFO warning is
  emitted on the INFO path (doc_v12 D12-6).
- `--mind` no longer evaluates a tail `-i/-e` expression blind: the sketch
  path has no record bytes, so `--mind` with such an expression takes the
  two-pass path, and `--mind --indep-pairwise` with one is a usage error
  instead of silently dropping every record (doc_v12 D12-1).
- Expressions see the header's tag types: `DB=1`/`DB=0` test a flag's
  presence (they used to match nothing), a multi-valued INFO tag compares
  any element (`DP4=0,0,30,40` at `DP4>10` used to miss), and a lone tag
  with no comparison is a compile-time error instead of a silent
  constant (bcftools evaluates it to false; doc_v12 D12-7).
- One INFO number grammar for the filter and the writer: `XI=1e3` no longer
  filters as 1000 while writing `.` to BCF, and an unparsable numeric tag
  is missing rather than compared byte-wise (`XI=abc` used to pass
  `XI>4`; doc_v12 D12-14).
- A GT-less record keeps its input `AC/AN` under a sample subset (it used
  to be recomputed from nothing as `AC=0;AN=0`; bcftools' `bcf_get_fmt`
  null check, vcfview.c:343; doc_v12 D12-2).
- `ALT='.'` counts zero ALT alleles, so `-m 2`/`-M 1`/`N_ALT=0` match
  bcftools on invariant sites; `-M 1` means "invariant only" instead of
  being dropped as "off" (doc_v12 D12-10).
- Records whose INFO/FORMAT counts do not fit the BCF header fields are
  refused instead of packed with the overflow leaking into the next field
  (65536 INFO keys or 300 FORMAT keys produced an unreadable file;
  doc_v12 D12-9).
- `--fill-tags` validates its list: an unsupported tag (AF, F_MISSING, a
  typo) is a usage error naming the supported set, like `bcftools
  +fill-tags`; it used to be a silent no-op that the default `-x` then
  turned into a deletion (doc_v12 D12-11).
- `[[site.info]]`'s `quantifier = "any"` requires an element to satisfy the
  comparison (it was true whenever an element parsed; doc_v12 D12-8).
- An INFO `AC` vector whose length does not match the ALT count falls back
  to the genotype counts (htslib's `bcf_calc_ac` check, vcfutils.c:56-64;
  doc_v12 D12-13).

### Added

- `reed convert` (docs/rounds/v13): one parallel pass over the input, any
  combination of targets. Phylogenetic matrices with vcf2phylip 2.9
  semantics (Ortiz 2019, DOI 10.5281/zenodo.2540861) — `--phylip`,
  `--fasta`, `--nexus`, `--nexus-binary` (SNAPP), plus `--outgroup`,
  `--min-samples-locus`, `--resolve-iupac` (deterministic, seed 0),
  `--write-used-sites`; the vcftools `--012` triple; the TASSEL hapmap
  (`.hmp.txt`); a reference-driven whole-site MSA with vcf2msa semantics
  (`--msa --reference`, `--region`, `--flank`, `--collapse-suffix`); and
  `--sfs`, byte-identical to `reed popg --sfs` (one shared request builder,
  one SFS fold). Outputs are independent of `--threads`; parity is gated by
  `gold/popgen/run_convert_parity.sh` (vcf2phylip five files, vcftools
  `--012`, TASSEL read-back, vcf2msa sequences). Deliberate divergences from
  the reference tools are registered in `docs/spec/cli.md` §8b: malformed
  records are skipped and counted, `--012` writes `-1` for half calls and
  haploids (vcftools writes 1), the hapmap `alleles` column lists the
  declared REF/ALT alleles, the MSA pads indels deterministically instead of
  calling an external aligner, and an allele set outside the IUPAC table
  skips the record where vcf2phylip crashes.
- `reed-core`'s `iupac` module: the IUPAC ambiguity table (62 entries), the
  binary-NEXUS genotype table, and the inverse lookup (`ambig_expand`: `R` →
  `AG`), the single implementation shared by every convert target.
- `reed convert --vcf`: the hapmap → VCF reverse direction (v13 C4), so the
  two formats convert both ways. The input is a hapmap (`.hmp.txt`, one- or
  two-letter cells auto-detected); `chrom`/`pos`/`rs#`/`alleles` are copied
  verbatim (`ID` = `rs#`, REF/ALT in declared order), each IUPAC cell expands
  to a `GT` (`R` → `0/1`), only `N`/`?` read as missing, and `alleles=NA`
  becomes an all-missing record. `-O v` (default) or `-Oz[level]` writes
  plain text or BGZF, `-o` omitted streams to stdout; the header carries no
  `##contig` (TASSEL writes none) and its `##reed_*` provenance lines are
  controlled by `--no-version`/`--no-command`. Exclusive with every
  VCF-input target and with sample subsetting; malformed rows are skipped and
  counted. On the lossless domain (biallelic diploid SNPs, no missing) a VCF
  → hapmap → VCF round trip reproduces every record, and TASSEL 5 reads
  reed's VCF back (`gold/popgen/run_convert_parity.sh` gate 6).
- `-g/--genotype [^](hom|het|miss)` (bcftools `-g`) for `reed filter` and
  `reed popg`: a post-mask site rule over the genotype classes
  (`TailRules::gt_class` → `TailRule::GtClass`). A half-call (`0/.`) reads
  as missing (bcftools `bcf_gt_type` UNKN), so `SiteCounts` gains
  `n_half_called`, counted by `count_site`, the plane sketch and the
  single-pass pipeline. Repeated flags AND (bcftools keeps only the last
  one); exact duplicates collapse and opposite polarities (`hom` + `^hom`)
  are a compile-time error. Funnel ids are `after_mask.gt_class.(hom|het|
  miss|^hom|^het|^miss)`. Diploid parity against bcftools 1.24 PASS over
  every call shape; haploid sites stay policy-gated (`require_diploid`
  rejects, `allow_haploid` counts them as half calls, so `-g hom`/`^miss`
  diverge from bcftools there by the same convention as `--maf`/`--mac`).
  A record whose FORMAT lacks the GT field skips every `-g` rule and
  survives (`SiteCounts::has_gt`; bcftools' `bcf_get_fmt` null check does
  the same) — distinct from a GT whose values are all `./.`. Class names
  are matched exactly, case-insensitively (bcftools `strcasecmp`; no
  surrounding whitespace is trimmed). Like bcftools, the rules evaluate
  after subsetting, over the kept samples only.

### Changed

- WC84 Fst (`--fst wc84` pairwise and `--fst-global`): a site where a
  population has no called alleles now contributes only the finite
  Weir-Cockerham `c` component (h̄/2 over the typed populations) instead of
  being skipped — pixy / scikit-allel `nansum` semantics, for both the
  pairwise and the r-population forms (doc_v7 D7-20). vcftools and plink2
  drop such sites (previously reed's behavior); hierfstat recomputes them
  with the typed populations only. Parity: pairwise == pixy to 13
  significant digits on a masked Chr01 with 1–3% all-missing sites
  (0.0086203 / 0.0020595 / −0.0005171), global == scikit-allel r=3
  (0.003237); the vcftools gold gate still passes (the benchmark data has
  no such sites, and values there are bit-identical).

- The engine's abort signal is an atomic flag (doc_v6 A6-1, Y5): `is_aborted`
  no longer flips back when a `select!` waiter consumes a channel token, the
  wake-up channel is sized by `--threads` (which may exceed the core count),
  and registered byte-permit pools are notified on trigger, so
  `BytePermits::acquire` waits on a condvar instead of polling every 10 ms
  (`docs/spec/architecture.md` §7.5). Every engine entry point now builds its
  abort channel from the run's thread budget (`abort_for`).
- In-flight chunk permits live in a drop guard (A6-1, Y5): a worker failure
  or panic, an aborted send, or a dropped pending message all return the
  permit and the in-flight accounting, so the reader can no longer block
  forever on a leaked permit. `StageMetrics::in_flight_final_bytes` makes the
  residue assertable.
- `threads <= 2` now takes a real sequential path (A6-3, Y5): one thread
  reads, inflates, computes and writes, with the same seam stitcher and sink
  as the parallel path; output is byte-identical (locked by a test). The
  `is_sequential` selector documented in §7.3 finally has a caller.
- The seam stitcher is bounded (A6-2, Y5): a record reassembled across chunk
  seams may not exceed `reed_engine::pipeline::MAX_RECORD_BYTES` (1 GiB,
  shared with the BCF frame chain), so a newline-free input fails the run
  instead of growing the reducer's fragment to the file size.
- `StageMetrics` gains `reader_wait`, `reducer_wait`, `inflate_time`,
  `job_time`, `deflate_time` and `in_flight_final_bytes` (A6-4, Y5). Output
  deflate time rides `ChunkOutput` to the sink (`OrderedSink::deflate_time`),
  so the three worker stages are measured instead of inferred; `--timing`
  and `--verbosity 2` print `reader_wait_ms`, `reduce_wait_ms`, `inflate_ms`,
  `job_ms` and `deflate_ms`.

- One pass0 serves every filter-side global analysis (doc_v6 03 §7, Y1e):
  the `MatrixJob`/`MatrixSink` pair fills per-contig site-major plane
  matrices that `--mind`, `--indep-pairwise` and the sketch path's mid
  stage all read back. The three previous normalizations
  (`MissingJob`, `SketchJob`, `LdPruneJob`) and `reed_filter::sketch` are
  deleted. `--mind` counts come from the miss-plane column popcounts
  (`called(s) = 2·sites − pop(slot0_miss ∨ mask) − pop(slot1_miss ∨ mask)`);
  `--mind` + `--indep-pairwise` now combine in one run (they previously
  could not share the first open).
- The pruning runs after every site rule: without `--mind` the pass0
  evaluates the tail rules inline and the pruner sees only post-tail
  sites (plink2's load-then-prune order, 02-cli §6 step 6b); with `--mind`
  the mid stage's survivor bitmap filters the pruner. The report's
  `ld_prune` counts are post-tail (`n_seen` shrinks when tail rules drop
  sites).
- `PipelineOptions::scan_head_failures` and `RecordView::head_passed` are
  gone (their only user was the sketch pass).
- The window-layout output writers are `reed::windows` (`WindowLayoutOptions`,
  `window_kinds`) instead of the previous vendor-named module (doc_v6 05 §7.2,
  E6-13); the output bytes are unchanged.
- The unsound AVX-512 counting kernel and `CountTable` were deleted
  (doc_v6 03 §5, A6-22); the mmap is now the workspace's only audited
  `unsafe` site.
- `ContigId` is a single newtype shared by every layer (doc_v6 A6-23).
- The whole tree is `cargo fmt`-clean under the pinned 1.99 toolchain, and the
  CI clippy command (`--all-targets -- -D warnings`) passes (doc_v6 05 §1–§2).

- I/O hot paths reuse their buffers (doc_v6 03 §3, A6-5/6/8, Y2): the
  deflater keeps a resident output scratch instead of zero-filling the
  worst-case bound per block, `inflate_chunk` walks the chunk's block table
  into exact-size windows (no second header parse, no per-block 64 KiB
  memset), and the chunk tail scan uses `memrchr`. The filter gold run
  drops from 460 to 400 CPU-s.
- The streaming BGZF source keeps empty members in the chunk (A6-7, Y2):
  dropping them used to shift every later chunk's `compressed_start`, so
  virtual offsets computed from stdin/pipe input could be wrong. The
  streaming and mmap sources now produce identical chunk metadata, locked
  by a test with an empty-member fixture.
- Per-record allocations left the hot paths (A6-11/12, Y2): the error
  `Locus` is built from borrowed bytes only when an error carries it
  (`LocusSource`), INFO edits write into the output buffer directly, and
  the index writer no longer builds a contig-name `String` per record on
  the (default) CSI path. The BCF GT decode table is one `Arc` per run
  instead of one per worker (A6-10).
- The filter's survivor emit is one method (`emit_survivor`) for the text
  and native BCF paths, the stale-tag list builds once per chunk, and the
  mask rules run columnar — one contiguous pass per rule over the scanned
  FORMAT columns, integer thresholds compared in the integer domain
  (A6-15/17/20, Y3). The `--set-gts .` expression path reuses the
  pipeline's scanned FORMAT table in place (no per-sample allocation) and
  expressions compile to operator enums at parse time (A6-18). INFO reads
  go through a per-record fragment index shared by the typed head rules
  and the expressions (`reed_io::vcf::info::InfoIndex`, A6-18). The
  exact-test HWE (`--fill-tags HWE`, `--hwe`, `--hardy`) is memoized per
  worker through `HweCache` (A6-16); the emitted p-values are unchanged.
  Full-gold CPU moved 398.2 → 395.6 s; the measurements behind the
  roadmap's larger estimate are recorded in `docs/rounds/v6/STATUS.md`
  (Y3 deviation).

### Added

- Streaming per-contig pruning (Y1f): each contig's matrix feeds its own
  prune thread as soon as the stream closes it, so one matrix is resident
  at a time; a single contig over the `--memory` budget fails with
  `MemBudgetExceeded`.

- `LICENSE-MIT` and `LICENSE-APACHE` (the manifests declared
  `MIT OR Apache-2.0` with no license text until now).
- `reed --version` now prints the source state the binary was built from
  (`reed 0.3.0 (v0.3.0-12-g8709f1f 2026-10-10, rustc 1.99.0)`), injected by
  `build.rs` (E6-11).
- README sections: Installation, Citation, Input requirements, Exit codes
  (E6-12).

## [0.3.0] - 2026-10-10

The current development line (the W0–W6 / X0–X7 rounds in `docs/rounds/`).

### Added

- `reed stats`: π, dxy, Fst (`wc84`/`hudson`, pairwise and global), Watterson
  θ, Tajima's D, PBS, projected SFS, f-statistics (`--f2/--f3/--f4/--dstat/--f4-ratio`),
  LD (`--ld-decay`, `--ld-pairs`), per-site/per-sample reports
  (`--freq/--hardy/--het/--missing/--depth`) and the bcftools-stats
  compatible `--summary`, all in one pass.
- `reed filter`: bcftools-view-style hard filtering, genotype masking,
  INFO recomputation (`--fill-tags`), `-q/-Q` with the `:type` suffixes,
  LD pruning (`--indep-pairwise`) as a site rule, index writing (`-W`/`--write-index`).
- `reed index`: CSI (default) and TBI, parallel.
- Window layout output (`PREFIX_pi.txt`, `PREFIX_dxy.txt`, …) verified
  column-for-column against pixy 2.2.4; `docs/rounds/v5/05-algorithms.md`.
- BCF read and write (native, no text round-trip for the filter decision),
  BGZF chunk streaming with mmap and streaming sources.
- The deterministic chunk engine: ticket window, ordered reducer, byte-budget
  backpressure; output byte-identical across thread counts and strategies.

### Verified against

bcftools, plink2, vcftools, pixy, scikit-allel, easySFS, admixtools-style
f-statistics; see `gold/` and the `Verification` section of the README.

## [0.2.0]

The pre-migration architecture (`docs/rounds/v1/design.md`), frozen as the
differential-testing oracle (`crates/oracle`). Not a published release; the
history is recorded in `docs/rounds/v2/STATUS.md`.
