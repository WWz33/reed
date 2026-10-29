# pca

## 用法

```
reed pca [OPTIONS] -o PREFIX <INPUT>
```

样本关系矩阵（GRM）的主成分。默认输出前 10 个主成分，用 `binom2` 标准化与 `pairwise` 缺失处理。只统计双等位、二倍体 GT 的位点：多等位位点、单倍体与半调用、没有已调用 ALT 的位点都跳过；一个可用位点都没有时报 `no biallelic sites with a diploid GT were found`。

## 参数

### 主成分

```
<INPUT>                     VCF/VCF.gz/BCF；- 读 stdin
    --pcs INT               输出前几个主成分，默认 10，最小 1
    --stand binom2|binom    标准化口径：binom2 用 √(2p(1−p))（默认），
                            binom 用 √(p(1−p))。两种口径成分相同、特征值差 2 倍
    --missing pairwise|meanimpute
                            缺失基因型处理：pairwise（默认）按每对样本的实际联合
                            调用数算分母；meanimpute 把缺失填成位点均值，分母用
                            全部位点
    --solver auto|dense|stream
                            求解器。auto（默认）在 N×N 矩阵装得下时走 dense
                            （全谱分解），否则走 stream（不物化矩阵的块 Lanczos）
```

### DAPC 与交叉验证

```
    --dapc                  在主成分上做判别分析，需要 --pop；写 PREFIX.dapc.ind、
                            PREFIX.dapc.eig、PREFIX.dapc.var
    --pop FILE              DAPC 的群体文件，样本<TAB>群体，无表头。
                            别名 --populations
    --n-da INT              保留几个判别轴，默认保留全部显著轴（svd² > 1e-10）
    --force-samples         群体文件里的未知样本只警告并跳过；默认报错
    --xval                  交叉验证选保留的 PC 数；写 PREFIX.dapc.xval。需要 --dapc
    --xval-rep INT          每个候选重复几次，默认 30
    --xval-training FLOAT   训练集比例，默认 0.9；最小的群体留不满时自动下调
    --xval-result overall|group-mean
                            正确率口径：overall 总体正确率（默认）、
                            group-mean 各群体等权平均
    --xval-seed INT         随机种子，默认 42
    --xval-pca LIST         候选 PC 数，逗号分隔；默认按固定刻度生成约 10 个候选
                            并去重
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

`--missing pairwise` 与 `--solver stream` 不兼容，同给报错并给出三个出路：换 `dense`、换 `meanimpute`、或提高 `--memory`。内存不够时错误信息会建议先跑 `reed filter --indep-pairwise` 减位点。

## 算法

数据流分三步。

1. 标准化。 每个位点估计等位频率 `p`（分母只计已调用剂量），把剂量 `x ∈ {0,1,2}` 变换为

```
--stand binom2（默认）:  W = (x − 2p) / √(2p(1−p))
--stand binom:           W = (x − 2p) / √(p(1−p))
```

缺失剂量在标准化空间记 0（即该位点的均值）；单态或固定 ALT 的退化列整列置 0，不进入分子。

2. 关系矩阵。

```
--missing pairwise（默认）:  A_ij = Σ_{s: i,j 都已调用} W_is·W_js / (联合调用位点数)
--missing meanimpute:        A_ij = Σ_s W_is·W_js / m           （m = 位点总数）
```

`pairwise` 的分子分母都只数两位点都已调用的位点，每个样本对用自己的分母；`meanimpute` 因为缺失已填成均值，分母就是全体位点。两套分母算出来的 `A_ij` 不会相同：`meanimpute` 等于拿均值补缺失，把关系往均值拉；`pairwise` 对覆盖低的样本对分母小、波动大。跟别的软件对坐标时，用同一种缺失处理口径才有可比性。

3. 特征求解。`dense` 物化 N×N 矩阵做稠密对称特征分解，取全谱；`stream` 按块扫描基因型迭代求解，固定种子，结果可复现。

## 输出

| 文件 | 内容 |
|---|---|
| `PREFIX.eigenvec` | 表头 `#FID IID PC1..PCk`，tab 分隔；FID 与 IID 都用样本名；特征向量单位范数，9 位小数。特征向量整体符号不定，跨运行比较请按方向比 |
| `PREFIX.eigenval` | 特征值，降序，每行一个，9 位小数（= GRM 的特征值） |
| `PREFIX.pve` | 每个成分占全谱的比例：λᵢ 除以全谱和（`stream` 模式下用 `Σw²/m` 作全谱和，等价于迹） |

```
PREFIX.eigenvec
#FID	IID	PC1	PC2	PC3
s1	s1	0.253837528	-0.199416559	0.047029802
s2	s2	0.498821822	-0.283700482	-0.197574393
s3	s3	-0.238191065	0.372053836	-0.082970657
…

PREFIX.eigenval
2.912527165
2.560773428
2.291419510

PREFIX.pve
0.218853435
0.192421917
0.172182096
```

## DAPC 判别分析

`--dapc` 在主成分之上做判别分析：先取前 `n_pca` 个成分（默认 `min(样本数−1, 50)`，特征值小于 1e-14 的成分不计入秩），样本得分取 `U·√λ`，再对得分做线性判别（组内白化 + 组均值 SVD）。

输出：

| 文件 | 内容 |
|---|---|
| `PREFIX.dapc.ind` | 表头 `#IID POP ASSIGN LD1..LDn_da`：样本、真实群体、判别指派的群体、判别坐标 |
| `PREFIX.dapc.eig` | 判别轴特征值（每行一个） |
| `PREFIX.dapc.var` | `#n_pca n_da var` 一行：保留的 PC 数、判别轴数、这些 PC 承载的全谱方差比例 |

```
PREFIX.dapc.ind
#IID	POP	ASSIGN	LD1
s1	popA	popA	3.802819344
s2	popA	popA	1.383581045
s3	popA	popA	3.283469062
…

PREFIX.dapc.var            PREFIX.dapc.eig
#n_pca	n_da	var          75.856959395
5	1	0.813187618
```

## 交叉验证

`--xval` 自动选保留的 PC 数。PCA 在全量数据上只算一次；每次重复按固定种子派生随机流，从每个群体按训练比例抽训练样本，训练集上重拟合判别，在留出集上算正确率。重复结束后按候选 `n_pca` 汇总。

评价量是 `√(mean((成功率 − 1)²))`，取误差最小的候选；平局取最大的候选。随机基线 `1/群体数` 一并报告。结果写 `PREFIX.dapc.xval`：`#n_pca mean_success` 每候选一行，后接 `#best_n_pca`、`#best_mean_success_n_pca`、`#random_chance` 三行摘要。

```
PREFIX.dapc.xval
#n_pca	mean_success
1	0.710000000
2	0.600000000
3	0.400000000
4	0.615000000
5	0.840000000
6	0.735000000
7	0.495000000
8	0.355000000

#best_n_pca	5
#best_mean_success_n_pca	5
#random_chance	0.500000000
```

选中的 PC 数覆盖 `--pcs` 并重新拟合最终 DAPC，同时打印一行 note：

```
reed pca --dapc --xval: best n_pca=5 (mean success 0.8400), highest mean success at n_pca=5, random chance 0.5000
  note: --pcs 10 was overridden by the cross-validated choice (5)
```

## 注意事项

- LD 剪枝是上游步骤：先 `reed filter --indep-pairwise`，把剪枝后的文件喂给 `pca`。本命令不收 `--indep-pairwise`。
- 关联大效应位点会污染 PCs；若做关联分析，常见做法是用剔除目标位点后的 PCs 作协变量。
- 特征向量符号任意：比较两批 PCs 请用方向（内积绝对值），逐元素差没有意义。

## 例子

```bash
# 输出前 20 个主成分
reed pca pruned.vcf.gz --pcs 20 -o pcs

# 换标准化口径与缺失处理
reed pca pruned.vcf.gz --stand binom --missing meanimpute -o pcs

# 强制流式求解，限制内存
reed pca pruned.vcf.gz --solver stream --memory 2000 -o pcs

# 只用一部分样本
reed pca pruned.vcf.gz -S samples.txt -o pcs

# DAPC：保留 5 个 PC 做判别分析
reed pca pruned.vcf.gz --dapc --pop pops.txt --pcs 5 -o dapc

# DAPC + 交叉验证选 PC 数（100 次重复）
reed pca pruned.vcf.gz --dapc --pop pops.txt --xval --xval-rep 100 -o dapc
```
