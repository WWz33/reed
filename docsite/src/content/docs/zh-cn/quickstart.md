---
title: 快速开始
---

用一条完整的分析线走完 reed 的主要命令：索引、QC、过滤、LD 剪枝、群体统计、PCA、建树、格式转换。

## 准备数据

输入是一个 VCF 或 BCF。本页假设：

- `raw.vcf.gz`：原始 calls，包含不变位点；
- `pops.txt`：群体文件，`样本<TAB>群体` 两列，无表头：

```
sampleA	popA
sampleB	popA
sampleC	popB
sampleD	popB
```

群体文件里出现 header 没有的样本名会报错；想让不认识的样本只给警告并跳过，加 `--force-samples`。

## 第一步：索引

```bash
reed index raw.vcf.gz
```

写 `raw.vcf.gz.csi`。后面 `filter -r/-R` 的区域查询会用到。要 TBI 就加 `-t`。

## 第二步：QC，先看数据长什么样

```bash
mkdir -p qc
reed stats raw.vcf.gz --freq --hardy --het --missing --depth --summary -o qc/raw
reed stats raw.vcf.gz --hist qual,dp-site,gq -o qc/raw
```

第一行给出位点等位频率、HWE 精确检验、每样本近交系数、缺失率、平均深度和一份汇总。第二行给 QUAL、位点深度、GQ 的直方图，用来选过滤阈值。输出都以 `qc/raw` 为前缀（`qc/raw.afreq`、`qc/raw.hardy`、`qc/raw.dp-site.hist` 等）。

直方图文件长这样，整数指标末尾附均值、中位数和 95 分位：

```
bin_lower	bin_upper	count	percent
2	2	126	2.520000
3	3	47	0.940000
4	4	58	1.160000
5	5	70	1.400000
…
#mean	31.012800
#median	30
#p95	58
```

## 第三步：硬过滤

位点阈值和基因型掩码可以放在一条命令里。求值顺序是：样本子集 → 位点类型 → 表达式 → 基因型掩码 → 样本缺失率 → 位点缺失率与 HWE 与频率阈值。基因型先被掩码，再参与位点过滤：

```bash
reed filter raw.vcf.gz -o clean.vcf.gz -Oz6 \
    --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 --hwe 1e-6 -W
```

跑完 stderr 给一行摘要和过滤漏斗，每条规则过后还剩多少位点：

```
reed filter: 500 -> 131 records (10 of 10 samples), strategy single_pass, 0.0s
funnel:
  site.biallelic_gate          500
  after_mask.missing           177
  after_mask.maf               131
  after_mask.hwe               131
```

- `--min-dp 5 --min-gq 20`：DP 低于 5 或 GQ 低于 20 的基因型置 `./.`。DP 或 GQ 本身缺失的基因型不受影响。
- `--geno 0.1`：掩码后位点缺失等位比例超过 10% 的位点丢弃。
- `--maf 0.05`：次要等位频率低于 5% 的位点丢弃。不变位点（MAF=0）也会被它滤掉，滤完的文件喂 `popg` 时要记得 `--variant-only`。
- `--hwe 1e-6`：HWE 精确检验 p 小于阈值的位点丢弃。
- `-Oz6`：BGZF 压缩，级别 6。
- `-W`：边写边建索引（需要 BGZF 输出和 `-o`）。

`-W` 不给就不写索引。`-O` 还接受 `v`（纯文本 VCF）、`b`（BCF）、`u`（未压缩 BCF）。

表达式过滤是另一条路，适合"按某个 INFO 值挑位点"这类规则：

```bash
# 丢掉 QUAL 低于 30 或 MAF 低于 0.05 的位点
reed filter raw.vcf.gz -i 'QUAL>30 && MAF>=0.05' -Oz -o clean2.vcf.gz

# 深度不足的基因型置为缺失，位点本身保留
reed filter raw.vcf.gz -e 'FMT/DP<5' --set-gts . -Oz -o masked.vcf.gz
```

表达式语法见[过滤表达式](../expressions/)。

## 第四步：LD 剪枝

PCA 之前先剪 LD。`--indep-pairwise` 直接输出剪枝后的 VCF，一条命令完成选位点和提取：

```bash
reed filter clean.vcf.gz --indep-pairwise 50 5 0.2 -Oz6 -o pruned.vcf.gz
```

三个数是窗口、步长和 r² 上限：窗口 50 个变异位点，每次前进 5 个，窗口内 r² 严格大于 0.2 的位点对，保留次要等位频率更高的那个。窗口写 `50kb` 就按碱基距离算（此时步长必须是 1）。

## 第五步：窗口化群体统计

```bash
reed popg clean.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out
```

一次开多个统计量，每个写一个文件（`out_pi.txt`、`out_dxy.txt`、`out_fst.txt`、`out_watterson_theta.txt`、`out_tajima_d.txt`）。

π 和 dxy 的分母要靠不变位点，`popg` 默认假设文件里有不变位点（全基因组 calls），启动时会检查。`--maf` 这类频率过滤会把不变位点也滤掉，所以过滤后的文件直接跑会收到这个错误：

```
reed: job failed: input has no invariant sites (122 records, none with ALT = "."):
declare --variant-only for an all-variant file, or --bypass-invariant-check to
compute over present sites only (both overestimate π/dxy denominators)
```

两种典型情况：

- 全基因组 calls：直接跑。空窗口默认也输出、写 `NA`，不想看到就加 `--skip-empty`。
- 过滤后只剩变异位点的文件（如 `clean.vcf.gz`）：加 `--variant-only`，窗口内文件不含的位置按"全部样本已调用的单态位点"计入分母。

两条路都不走、坚持只算文件里有的位点，用 `--bypass-invariant-check`，reed 会警告 π 和 dxy 被高估。

`out_pi.txt` 的行：

```
pop	chromosome	window_pos_1	window_pos_2	avg_pi	no_sites	count_diffs	count_comparisons	count_missing
popA	chr1	1	10000	0.0003778491492837536	9	170	320	85
popB	chr1	1	10000	0.00040893523488217553	9	184	354	51
```

PBS、f 统计量和 LD 衰减：

```bash
# PBS：A 是焦点群体，--pbs 可重复给多个三元组
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --pbs popA,popB,popC --pbs popA,popB,popD --window 50k -o out

# f 统计量：块刀切标准误，默认 5 Mb 块，--block-size 可改
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --f4 popA,popB;popC,popD --dstat popA,popB;popC,popD -o out

# LD 衰减：r² 按 1 kb 分箱，默认看到 100 kb；--ld-pairs 另给逐对表
reed popg clean.vcf.gz --pop pops.txt --variant-only --ld-decay 200k -o out

# 投影 SFS：dadi 和 fastsimcoal2 两种布局；两群体给两个投影数
reed popg clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out
```

群体名必须和 `pops.txt` 里的一致；f 统计量的群体组合写法是 `A,B`、`C;A,B`、`A,B;C,D`、`A,B;C,D;X`。

## 第六步：PCA 与判别分析

```bash
reed pca pruned.vcf.gz --pcs 20 -o pcs
```

写 `pcs.eigenvec`（表头 `#FID IID PC1..PCk`）、`pcs.eigenval`、`pcs.pve`：

```
#FID	IID	PC1	PC2	PC3
s1	s1	0.253837528	-0.199416559	0.047029802
s2	s2	0.498821822	-0.283700482	-0.197574393
```

带群体标签的判别分析，以及用交叉验证选保留的 PC 数：

```bash
reed pca pruned.vcf.gz --dapc --pop pops.txt --pcs 5 -o dapc

# 交叉验证：默认 30 次重复、训练集 90%、固定随机种子
reed pca pruned.vcf.gz --dapc --pop pops.txt --xval -o dapc
```

交叉验证会打印它选的 `n_pca` 并用它覆盖 `--pcs`：

```
reed pca --dapc --xval: best n_pca=5 (mean success 0.8167), highest mean success at n_pca=5, random chance 0.5000
  note: --pcs 10 was overridden by the cross-validated choice (5)
```

候选曲线和摘要写进 `dapc.dapc.xval`。

## 第七步：个体间矩阵与建树

```bash
reed distance clean.vcf.gz --type grm -o rel        # rel.rel + rel.rel.id
reed distance clean.vcf.gz --type king -o kin       # kin.king + kin.king.id
reed distance clean.vcf.gz --type phi -o rela       # rela.relatedness2
```

从距离型矩阵建树，加自展支持度：

```bash
reed distance clean.vcf.gz --type 1-ibs --tree nj --bootstrap 100 -o tree
```

`tree.nwk` 里内节点的数字是自展支持度：

```
(((s1:0.2010,s4:0.2173)35:0.0187,((s2:0.2185,s6:0.2104)25:0.0134,
  (s5:0.1824,s9:0.1875)60:0.0347)20:0.0148)20:0.0087,
 ((s3:0.1950,s10:0.1764)40:0.0181,s8:0.2309)20:0.0090,s7:0.2161)100;
```

`--tree upgma|nj` 只接受对角线为零的距离型矩阵；`grm`、`king`、`ibs` 这类相似度会被拒绝并提示改用 `1-ibs`、`allele-ct`、`ajk` 或 `discordance`。随机种子默认 42，结果可复现。

## 第八步：格式转换

```bash
# 系统发育与基因型矩阵，一次读入多目标输出
reed convert clean.vcf.gz --phylip --fasta --nexus --nexus-binary --012 --hapmap -o matrices

# 参考驱动的 MSA
reed convert clean.vcf.gz --msa --reference ref.fa --region chr1:1-100000 -o aln

# 投影 SFS（与 popg --sfs 输出一致）
reed convert clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out

# hapmap 反向转回 VCF
reed convert panels.hmp.txt --vcf -Oz -o panels.vcf.gz
```
