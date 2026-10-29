# distance

## 用法

```
reed distance [OPTIONS] --type KIND -o PREFIX <INPUT>
```

个体对个体的矩阵：等位差距离、IBS、GRM、协方差、KING、PI_HAT、相关度、不一致率；`--tree upgma|nj` 从距离矩阵建 Newick 树，`--bootstrap` 重采样位点给枝支持度。默认 `--type grm`、`--format triangle`。

输入门与 `pca` 相同：只统计双等位、二倍体 GT 的位点，无已调用 ALT 的位点跳过。

## 参数

### 矩阵

```
<INPUT>                     VCF/VCF.gz/BCF；- 读 stdin
    --type KIND             矩阵类型，默认 grm。可选：allele-ct（别名 dist）、ibs、
                            1-ibs、grm（别名 rel）、cov、king、pihat（别名 genome）、
                            phi（别名 relatedness2）、ajk（别名 relatedness）、
                            discordance（别名 gtcheck）。公式见下节
    --format FMT            矩阵型输出的布局：triangle（默认，下三角含对角）、
                            square、square0（方阵，右上三角置 0）、bin（小端 f64）
    --stand binom2|binom    标准化口径，同 pca，默认 binom2
    --mean-impute           仅 grm/cov：分母用全部位点数，不用两样本的实际联合调用数
    --unbounded             仅 pihat：关掉 PI_HAT 的 [0,1] 裁剪
    --nudge                 仅 pihat：PI_HAT² < Z2 时把系数吸附到最近的可行配置
    --mismatch-count        仅 discordance：报告失配数，而不是概率得分
    --error-probability INT 仅 discordance：Phred 尺度的基因型错误率，默认 40
```

### 建树与自展

```
    --tree upgma|nj         建树而不是写矩阵，写 PREFIX.nwk。upgma 算术平均连锁
                            （有根），nj 邻接法（无根）
    --digits INT            Newick 枝长的小数位数，默认 10；0 报错
                            （全部枝长会变成零）
    --bootstrap INT         自展：有放回重采样位点这么多次，支持度写成内节点标签
    --bootstrap-seed INT    自展随机种子，默认 42
    --bootstrap-cutoff FLOAT
                            支持度低于这个百分数的枝标签写空，默认 0（全保留）
```

### 样本与输出

```
-s, --samples [^]LIST       样本子集，逗号分隔，^ 前缀取补
-S, --samples-file FILE     样本名单文件，一行一个
    --keep FILE             同 -S，另一种拼写
    --remove FILE           要丢掉的样本名单
-o, --out PREFIX            输出前缀
    --threads INT           工作线程数，默认全部可用核
    --memory MB             内存预算，MiB，默认物理内存一半
    --verbosity INT         0 静默、1 摘要（默认）、2 阶段耗时、3 调试
```

## 矩阵类型与公式

标准先用 `--stand binom2`（`√(2p(1−p))`）或 `--stand binom`（`√(p(1−p))`）做，`p` 按已调用剂量估计。部分类型需要原始剂量，从标准化值反推：`dose = round(W·sd + mean)`，取值 0/1/2，缺失记 -1。

| `--type` | 公式 | 输出文件 |
|---|---|---|
| `allele-ct`（别名 `dist`） | 联合调用位点上的等位差之和 `Σ|x_i − x_j|`，缺失用 HWE 权重修正：权重 `w = p(1−p)(p²−p+1)`，分母是权重和减去各自缺失及双方同缺的权重，最后乘回总权重 | `PREFIX.dist` |
| `ibs` | `1 − d/2m`，d 是 allele-ct 距离，m 位点数 | `PREFIX.mibs` |
| `1-ibs` | `d/2m` | `PREFIX.mdist` |
| `grm`（别名 `rel`） | `A_ij = Σ W_is·W_js / 分母`；分母默认是两样本联合调用位点数，`--mean-impute` 改为总位点数 m。对角线是 `ΣW²/m`，不强行置 1 | `PREFIX.rel` |
| `cov` | 同 grm，但取消标准化里的除以 sd（保留中心化），得到协方差尺度 | `PREFIX.rel` |
| `king` | 联合调用位点上数四类基因型组合（ibs0、双杂合、het1-hom2、het2-hom1），`φ = 0.5 − (4·ibs0 + h1h2 + h2h1) / (4·(hethet + min(h1h2, h2h1)))`；分母为零记 NaN | `PREFIX.king` |
| `pihat`（别名 `genome`） | 逐对解 IBS 系数：`z0 = ibs0/(e00·n)`，`z1 = ibs1/(e11·n)`，`z2 = (oo − n(z0·e02 + z1·e12))/n`，期望系数按位点等位频率累计；裁剪到 [0,1]（`--unbounded` 关掉裁剪；`--nudge` 在 `π² < z2` 时把系数吸附到 `z0=(1−π)², z1=2π(1−π), z2=π²`）；`PI_HAT = z1/2 + z2` | `PREFIX.genome` |
| `phi`（别名 `relatedness2`） | `φ = (N_AaAa − 2·N_AAaa)/(N_Aa[i] + N_Aa[j])`，分母用两样本各自的杂合位点数之和（边缘计数）；对角线 0.5 | `PREFIX.relatedness2` |
| `ajk`（别名 `relatedness`） | 只用频率在开区间内的位点：`A_ij = mean (a−2p)(b−2p)/(2p(1−p))`；对角线用同一分子加 `x² − (1+2p)x + 2p²` 项 | `PREFIX.relatedness` |
| `discordance`（别名 `gtcheck`） | 基因型错误模型 `e = 10^(−E/10)`（`--error-probability` 默认 40），逐位点对每个可能的真基因型取 `cost(a,G)+cost(b,G)` 的最小值后求和；`--mismatch-count` 改为直接数失配数 | `PREFIX.gtcheck` |

矩阵型输出（dist、mibs、mdist、rel、king）附 `PREFIX.<suffix>.id` 边车（两列样本名）。表格型输出（genome、relatedness2、relatedness、gtcheck）按各自表格布局写，`.genome` 的 DST/PPC/RATIO 列写占位 0。

## 输出布局

`--format` 只作用于矩阵型输出：

| 格式 | 布局 |
|---|---|
| `triangle`（默认） | 下三角含对角线：第 i 行是 `d(i,0)..d(i,i)`，9 位小数，tab 分隔 |
| `square` | 完整 N×N 方阵 |
| `square0` | 方阵，右上三角置 0 |
| `bin` | 无表头的 n² 个小端 f64，二进制直接可读 |

## 建树

`--tree` 只接受对角线为零的距离型矩阵。相似度类型（grm、cov、king、ibs、pihat、phi）会被拒绝，错误信息列出可用的距离类型：`1-ibs`、`allele-ct`、`ajk`、`discordance`。

| 方法 | 算法 |
|---|---|
| `nj` | 邻接法：Q 矩阵 `Q(i,j) = (n−2)d(i,j) − S_i − S_j` 取最小者合并，新枝长 `d(i,u) = d(i,j)/2 + (S_i−S_j)/(2(n−2))`，归约 `d(u,k) = (d_ik + d_jk − d_ij)/2`；负枝长保留不截断；最后三个分类群并到一个中心节点；结果无根，共 `2n−3` 条边 |
| `upgma` | 算术平均连锁（Lance-Williams 更新），合并距离对半分成枝长；超度量输入精确还原；结果有根 |

建树前校验距离阵：方阵、有限值、对角线为零、对称（容差 1e-9）。Newick 输出的尾零剪掉，标签里的保留字符加单引号。

## 自展支持度

- 支持度按分支计：每个内节点的后代尖端集合（无根树取包含尖端 0 的二分）在重复中出现的比例，写成百分数；根节点的标签是 100。
- 每次重复重新标准化、重跑同一个 `--type`、重建同一种树。
- 失败的重复（如重采样后无可用位点）跳过并计数，支持度分母相应减小。

## 例子

```bash
# 遗传关系矩阵
reed distance clean.vcf.gz --type grm -o rel

# 协方差尺度（不除以 sd）
reed distance clean.vcf.gz --type cov -o cov

# 亲缘度矩阵（输入应是常染色体位点）
reed distance clean.vcf.gz --type king -o kin

# 相关度表
reed distance clean.vcf.gz --type phi -o rela
reed distance clean.vcf.gz --type ajk -o yang

# 距离矩阵与 IBS
reed distance clean.vcf.gz --type allele-ct -o d
reed distance clean.vcf.gz --type 1-ibs -o ibsd

# 不一致率，Phred 错误率 30
reed distance clean.vcf.gz --type discordance --error-probability 30 -o disc

# 完整方阵与二进制输出
reed distance clean.vcf.gz --type grm --format square -o rel
reed distance clean.vcf.gz --type 1-ibs --format bin -o ibsd

# NJ 树，100 次自展，只留支持度 ≥ 50 的标签
reed distance clean.vcf.gz --type 1-ibs --tree nj \
    --bootstrap 100 --bootstrap-cutoff 50 -o tree
```

`rel.rel`（`--format triangle` 默认布局，第 i 行到第 i 列）与 `rel.rel.id`：

```
1.310800633
-0.206295975	1.437578015
-0.296796369	-0.411552149	1.284675298

s1	s1
s2	s2
s3	s3
```

`rela.relatedness2`（`--type phi`，完整矩阵两个方向都写）：

```
INDV1	INDV2	N_AaAa	N_AAaa	N1_Aa	N2_Aa	RELATEDNESS_PHI
s1	s1	37	0	37	37	0.500000000
s1	s2	12	26	37	31	-0.588235294
s1	s3	6	29	37	30	-0.776119403
```

`kin.genome`（`--type pihat`，严格下三角；DST/PPC/RATIO 是占位列）：

```
FID1	IID1	FID2	IID2	RT	EZ	Z0	Z1	Z2	PI_HAT	PHE	DST	PPC	RATIO
s2	s2	s1	s1	OT	NA	1.0000	0.0000	0.0000	0.0000	-9	0	0	0
s3	s3	s1	s1	OT	NA	1.0000	0.0000	0.0000	0.0000	-9	0	0	0
```

`tree.nwk`，内节点的数字是自展支持度：

```
(((s1:0.2010492714,s4:0.2172579309)35:0.01874847015,
  ((s2:0.2185298336,s6:0.2104283332)25:0.01339241676,
    (s5:0.182356414,s9:0.1875419981)60:0.03474948761)20:0.01475450946)20:0.008743806348,
 ((s3:0.1950385303,s10:0.1764491831)40:0.01806552434,s8:0.2309362055)20:0.009008248687,
 s7:0.2161362545)100;
```
