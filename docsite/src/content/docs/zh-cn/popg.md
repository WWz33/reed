---
title: popg
---

## 用法

```
reed popg [OPTIONS] --pop FILE --pi --fst ... -o PREFIX <INPUT>
```

窗口化群体遗传统计。一次运行可以同时开多个统计量，每个写自己的文件；另有一份 `PREFIX.log` 记录命令、版本、每个统计量的方法、文献与过滤计数。别名 `popgen`。默认假设输入含不变位点，没有就拒绝运行。

## 参数

### 群体与统计量

```
    --pop FILE              样本<TAB>群体 两列，无表头。只给 --sfs 或 LD 类统计时
                            可以不给，全部样本合成一个名为 all 的组。
                            别名 --populations
    --pi                    每群体 π（Nei & Li 1979）
    --dxy                   每群体对 dxy（Nei 1987）
    --theta                 每群体 Watterson θ（Watterson 1975）
    --tajima-d              每群体 Tajima's D（Tajima 1989）
    --fst [wc84|hudson]     全部群体对的 Fst，默认 wc84（Weir & Cockerham 1984）；
                            hudson 走 Hudson 1992 / Bhatia 2013 的比值和
    --fst-global            与 --fst 同给：报告合并的 r 群体 WC84 全局 Fst，
                            写 PREFIX_fst_global.txt
    --pbs A,B,C             三元组 PBS，A 是焦点群体（Yi et al. 2010）；可重复
    --pbs-normed            报告 PBSn1（Malaspinas et al. 2016），而不是原始 PBS
    --f2 A,B                Patterson f2(A, B)；可重复
    --f3 C;A,B              Patterson f3(C; A, B)；可重复
    --f4 A,B;C,D            Patterson f4(A, B; C, D)；可重复
    --dstat A,B;C,D         Patterson's D(A, B; C, D)；可重复
    --f4-ratio A,B;C,D;X    f4-ratio，即 alpha 估计；可重复
    --block-size SIZE       f 统计量块刀切的块长，默认 5M
    --ld-decay [MAX_BP]     LD 衰减：MAX_BP 内的成对 r²，按 1 kb 分箱，默认 100k
    --ld-pairs              另写逐对 LD 表 PREFIX.ld；需与 --ld-decay 同给
    --sfs [folded]          投影 SFS，全基因组；写 PREFIX_sfs/dadi 与
                            PREFIX_sfs/fastsimcoal2。默认展开（REF 为祖先），
                            --sfs folded 折叠。需要 --proj
    --proj N[,N...]         投影到的染色体数。一个值广播到所有群体；
                            多个值按群体首现顺序一一对应
    --sfs-int               投影值取整
```

### 窗口与位点

```
    --window SIZE           窗口大小，默认 100,000 bp
    --step SIZE             步长，默认等于窗口大小（给 0 也按窗口大小处理）
    --bed-file FILE         BED3（0-based 半开）自定义窗口，每个区间一窗，
                            与 --window/--step 互斥
    --sites-file FILE       CHROM<TAB>POS（1-based），只统计列出的位点
    --levels LIST           报告哪些行：window,contig,genome，默认 window；
                            contig 行折叠整条 contig，genome 行折叠全部
    --skip-empty            丢弃没有位点落入的窗口（默认输出这些窗口，统计量写 NA）
    --contig-lengths FILE   header 缺 ##contig 行时，从 .fai 风格的
                            名字<TAB>长度 文件补长度；仍缺失则报错，不编造
```

### 不变位点

```
    --variant-only          声明输入只含变异位点：窗口内文件不含的位置按
                            "全部样本已调用的单态位点"补进 π/dxy/θ 的分母
    --bypass-invariant-check
                            只算文件里有的位点，打印警告 π and dxy are overestimated
```

两者都不给时，输入不含 `ALT = "."` 的记录就运行失败，错误信息给出两条出路。只请求 SFS 或 f 统计量时自动按 `--variant-only` 处理。

### 输入与区域

```
<INPUT>                     VCF/VCF.gz/BCF；- 读 stdin
-r, --regions LIST          只要这些区域，逗号分隔；需要输入旁边的索引
-R, --regions-file FILE     区域清单文件；需要索引
-t, --targets LIST          同样的区域写法，流式，不需要索引
-T, --targets-file FILE     区域清单文件；行首 ^ 表示排除。不需要索引
```

### 样本与位点

```
-s, --samples [^]LIST       样本子集，逗号分隔，^ 前缀取补
-S, --samples-file FILE     样本名单文件，一行一个
    --keep FILE             同 -S，另一种拼写
    --remove FILE           要丢掉的样本名单
    --force-samples         群体文件里的未知样本只警告并跳过；默认报错
-v, --types LIST            保留这些变异类型：snps、mnps、indels、other、ref
-V, --exclude-types LIST    丢掉这些变异类型
    --snps-only[=just-acgt] 只留 SNP；=just-acgt 要求每个 REF/ALT 都在 A/C/G/T 内
-m, --min-alleles INT       等位数下限，REF+ALT 计
-M, --max-alleles INT       等位数上限，REF+ALT 计
-f, --apply-filters LIST    FILTER 列在 LIST 里的记录才留
-i, --include EXPR          保留满足表达式的记录，可重复
-e, --exclude EXPR          丢掉满足表达式的记录，可重复
-q, --min-af FLOAT[:type]   等位频率下限，:type 同 filter
-Q, --max-af FLOAT[:type]   等位频率上限
    --geno FLOAT            缺失等位比例超过 FLOAT 的位点丢掉（掩码之后）
    --maf FLOAT             次要等位频率下限
    --max-maf FLOAT         次要等位频率上限
    --mac INT               次要等位计数下限
    --max-mac INT           次要等位计数上限；0 是内部的"关闭"值，会被拒绝
    --hwe P [midp]          HWE 精确检验 p 值下限
-g, --genotype [^](hom|het|miss)
                            按基因型组成选位点，可重复，之间是 AND
    --min-dp INT            基因型掩码，FMT/DP 下限。别名 --vcf-min-dp
    --max-dp INT            基因型掩码，FMT/DP 上限。别名 --vcf-max-dp
    --min-gq INT            基因型掩码，FMT/GQ 下限。别名 --vcf-min-gq
    --set-gts .             配合 -i/-e 按样本掩码，只接受 .
    --indep-pairwise W[kb] [S] R2
                            本命令不收：报错并要求先跑 reed filter --indep-pairwise，
                            再把剪枝后的文件喂进来
```

### 输出与全局参数

```
-o, --output PREFIX         输出前缀。窗口模式必须给，-o - 报错；
                            裸命令（无统计量）把汇总 TSV 写 stdout。别名 --out
    --threads INT           工作线程数，默认全部可用核
    --memory MB             在途内存预算，MiB，默认物理内存一半
    --verbosity INT         0 静默、1 摘要（默认）、2 阶段耗时、3 调试
```

群体 id 按群体文件首现顺序编号；重复群体名、空群体、一个样本分给两个群体都是用法错。

## 统计量与算法

### π（`--pi`，Nei & Li 1979）

对群体内每条染色体对，位点上的差异数是 `n_ref·n_alt`，比较数是 `C(n,2)`（`n` 是该群体在该位点已调用的染色体数）。窗口值

```
π = Σ位点的差异数 / Σ位点的比较数
```

分母只累计两条染色体在该位点都已调用的对（成对完整），缺失对的个数记进输出列 `count_missing`。窗口内无任何可比较位点时写 `NA`。

### dxy（`--dxy`，Nei 1987）

群体 A、B 间每条跨群体染色体对的差异数是 `a_ref·b_alt + a_alt·b_ref`，比较数是 `(a_ref+a_alt)·(b_ref+b_alt)`。窗口值同样是跨位点的比值和。`--variant-only` 时，比较数追加 `(窗口长度 − 文件位点数)·2N_A·2N_B`。

### Fst（`--fst [wc84|hudson]`，默认 wc84）

两个估计器共用一次位点遍历，全部输出到 `PREFIX_fst.txt`。

WC84（Weir & Cockerham 1984，r = 2 成对）：按等位基因拷贝数算方差分量。记 `p̄` 为合并等位频率，`h̄` 为平均杂合度，`n̄` 为平均基因拷贝数（各群体已调用拷贝数除以 2 后的均值），

```
s² = Σ nᵢ(pᵢ−p̄)² / (n̄(r−1))
a = (n̄/n_C)·(s² − (p̄(1−p̄) − (r−1)s²/r − h̄/4)/(n̄−1))
b = (n̄/(n̄−1))·(p̄(1−p̄) − (r−1)s²/r − (2n̄−1)h̄/(4n̄))
c = h̄/2
窗口值 = Σa / (Σa + Σb + Σc)
```

`n_C = (Σnᵢ − Σnᵢ²/Σnᵢ)/(r−1)`。某群体不足 2 条已调用染色体时该位点只贡献 `c` 分量。

Hudson（Hudson 1992；Bhatia et al. 2013 的比值和形式）：分子是"组间平均差减去组内平均差"，分母是组间平均差：

```
组内  mpd_a = n_ref_a·n_alt_a / C(n_a,2)      （要求 n_a ≥ 2）
组间  mpd_ab = (a_ref·b_alt + a_alt·b_ref)/(n_a·n_b)
Fst  = Σ(mpd_ab − (mpd_a+mpd_b)/2) / Σ mpd_ab
```

某位点缺分子项时按 0 计入分子和，分母和独立累计。

全局 Fst（`--fst-global`）三个以上群体时算合并的 r 群体 WC84（上面的公式取一般 r）。Hudson 没有 r 群体形式，`--fst-global` 与 `--fst hudson` 同给报错；两群体时全局值等于成对值，只提示不报错。

### Watterson θ（`--theta`，Watterson 1975）

位点按观测等位数 `n`（该群体已调用染色体数）分箱，`0 < alt < n` 的分离位点计一个"突变"。窗口值

```
θ_W = (Σ 位点突变数/a₁(n)) / n_sites,   a₁(n) = Σ_{k=1}^{n−1} 1/k
```

调和数表预计算到 2¹⁶，更大时按补偿求和直接累加。缺失数据通过每位点自己的 `n` 进入公式，位点不会被整体丢弃。`--variant-only` 时分母是窗口长度。固定 ALT 的位点贡献 0。

### Tajima's D（`--tajima-d`，Tajima 1989）

`raw_pi` 是 π 的原始和，`raw_theta` 是上节的 θ_W 原始值，`S` 是窗口总突变数。方差常数在窗口平均等位数处取值（四舍六入五成双取整）：

```
a₁ = Σ1/k,  a₂ = Σ1/k²   （k = 1..n−1）
b₁ = (n+1)/(3(n−1)),      b₂ = 2(n²+n+3)/(9n(n−1))
c₁ = b₁ − 1/a₁,           c₂ = b₂ − (n+2)/(a₁n) + a₂/a₁²
e₁ = c₁/a₁,               e₂ = c₂/(a₁²+a₂)
variance = e₁·S + e₂·S(S−1),   D = (raw_pi − raw_theta)/√variance
```

无位点、无突变或平均等位数 < 2 的窗口写 `NA`。输出附 `tajima_d_stdev`。

### PBS（`--pbs A,B,C`，Yi et al. 2010）

输入是三个成对的窗口 Hudson Fst（每对独立决定自己的覆盖样本，不做三方交集）。Fst 先截断到 `[0, 0.99999]`，再取 `t(f) = −ln(1−f)`：

```
PBS = (t_AB + t_AC − t_BC)/2      （A 是焦点群体）
```

`--pbs` 可重复，每三元组一份文件 `PREFIX_pbs.<A>-<B>-<C>.txt`；单三元组写 `PREFIX_pbs.txt`。输出附三个未截断的成对 Fst 列。

`--pbs-normed` 把 PBS 换成 PBSn1（Malaspinas 2016）：`PBS/(1 + (t_AB+t_AC+t_BC)/2)`。

### f 统计量（Patterson et al. 2012 的无偏校正形式）

等位频率 `p = n_alt/an`（按 ALT 算），`ĥ = n_ref·n_alt/(an(an−1))` 是有限样本杂合度校正项：

| 统计量 | 每位点贡献 |
|---|---|
| `--f2 A,B` | `(pA−pB)² − ĥA/anA − ĥB/anB` |
| `--f3 C;A,B` | `T = (pC−pA)(pC−pB) − ĥC/anC`；归一化变体另给 `ΣT/ΣB`，`B = 2ĥC` |
| `--f4 A,B;C,D` | `(pA−pB)(pC−pD)` |
| `--dstat A,B;C,D` | `Σ(pA−pB)(pC−pD) / Σ(pA+pB−2pA·pB)(pC+pD−2pC·pD)` |
| `--f4-ratio A,B;C,D;X` | 含 X 的 f4 分子和 ÷ 不含 X 的 f4 分子和 |

所有 f 统计量的标准误用块刀切：基因组被切成固定长度的块（默认 5,000,000 bp，`--block-size` 可改），逐块删除后算伪值

```
v_k = (Σnum − num_k)/(Σden − den_k),   SE = √((n−1)/n · Σ(v_k − v̄)²),   z = estimate/SE
```

块按 `(contig, pos/块长)` 键控，多线程分块合并安全。输出 `PREFIX.fstats`，列 `stat pops estimate se z n_blocks n_snps`；`--f2/--f3/...` 都可以重复给，同一份文件汇总。零分母写 `NA`。

### 投影 SFS（`--sfs [folded] --proj N[,N…]`）

全基因组（不开窗）。每群体一维谱，两个群体起另写每对的二维联合谱。投影用超几何分布

```
P(k | n, m, alt) = C(alt,k)·C(m−alt,n−k)/C(m,n)
```

（m 是该位点实际调用数，n 是投影目标），按众数锚定的比值递归计算，相对权重小于峰值 1e-15 截断。输出两套布局：

- `PREFIX_sfs/dadi/<pop>-<n>.sfs`（一维）、`PREFIX_sfs/dadi/<A>-<B>.sfs`（二维），三行：维度行、值行、mask 行；
- `PREFIX_sfs/fastsimcoal2/<pop>_MAFpop0.obs`（折叠）/ `<pop>_DAFpop0.obs`（展开）、`<base>_jointMAFpop<idxB>_<idxA>.obs`（二维，转置布局）。

mask 行按谱的结构性零位（k=0 与展开谱顶、折叠谱上三角）标记。

### LD 衰减（`--ld-decay [MAX_BP]`、`--ld-pairs`）

r² 在未定相的 0/1/2 剂量上算（未定相 Rogers-Huff r² 的整数实现），两位点取都已调用的样本（成对剔除）。`MAX_BP` 默认 100,000；距离按 1 kb 固定分箱，每箱记平均 r² 与对数。输出 `PREFIX.ld.decay`，列 `pop bin_start bin_end n_pairs mean_r2`：先全体（pop = `all`），多群体时每群体一组。`--ld-pairs` 另写逐对表 `PREFIX.ld`（列 `#CHROM_A POS_A ID_A CHROM_B POS_B ID_B UNPHASED_R2`）。

## 窗口机制

窗口网格锚定在步长的整数倍上：第一个窗口从包含首个位点之前的最近锚点开始，窗口末端裁到 contig 长度。输出坐标是 1-based 闭区间的 `window_pos_1`、`window_pos_2`。

## 输出文件与列

| 文件 | 列 |
|---|---|
| `PREFIX_pi.txt` | `pop chromosome window_pos_1 window_pos_2 avg_pi no_sites count_diffs count_comparisons count_missing` |
| `PREFIX_dxy.txt` | `pop1 pop2 chromosome window_pos_1 window_pos_2 avg_dxy no_sites count_diffs count_comparisons count_missing` |
| `PREFIX_fst.txt` | `pop1 pop2 chromosome window_pos_1 window_pos_2 avg_wc_fst \| avg_hudson_fst no_snps` |
| `PREFIX_fst_global.txt` | `pops chromosome window_pos_1 window_pos_2 avg_wc_fst no_snps` |
| `PREFIX_watterson_theta.txt` | `pop chromosome window_pos_1 window_pos_2 avg_watterson_theta no_sites raw_watterson_theta no_var_sites weighted_no_sites` |
| `PREFIX_tajima_d.txt` | `pop chromosome window_pos_1 window_pos_2 tajima_d no_sites raw_pi raw_watterson_theta tajima_d_stdev` |
| `PREFIX_pbs.txt`（多三元组 `PREFIX_pbs.<A-B-C>.txt`） | `popA popB popC chromosome window_pos_1 window_pos_2 pbs no_snps fst_12 fst_13 fst_23` |
| `PREFIX.fstats` | `stat pops estimate se z n_blocks n_snps` |
| `PREFIX.ld.decay` | `pop bin_start bin_end n_pairs mean_r2` |
| `PREFIX.ld` | `#CHROM_A POS_A ID_A CHROM_B POS_B ID_B UNPHASED_R2` |

`no_sites` 是窗口内文件包含的位点数（含被跳过的多等位记录）；`no_snps` 是推给该统计量的全部位点数。缺失值写 `NA`。Fst 出现负值是正常的：样本量小的窗口里方差分量估计可以低于零，解读时通常当 0 看。

## 例子

```bash
# 多统计量一次跑完（全基因组 calls 输入）
reed popg allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k --threads 56 -o out

# 变异位点输入 + 半重叠窗口
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --fst hudson --window 100k --step 50k -o out

# 全局 WC84 Fst
reed popg clean.vcf.gz --pop pops.txt --variant-only --fst --fst-global --window 1M -o out

# PBS：焦点群体 popA 对两个对照
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --pbs popA,popB,popC --pbs-normed --window 50k -o out

# f4-ratio 检验
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --f4-ratio popA,popB;popC,popD;popX --block-size 2M -o out

# 自定义窗口（BED）与位点清单
reed popg clean.vcf.gz --pop pops.txt --variant-only --pi --bed-file windows.bed -o out
reed popg clean.vcf.gz --pop pops.txt --variant-only --pi --sites-file sites.txt -o out

# 边过滤边统计
reed popg raw.vcf.gz --pop pops.txt --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 \
    --pi --window 100k -o out
```

`out_pi.txt`：

```
pop	chromosome	window_pos_1	window_pos_2	avg_pi	no_sites	count_diffs	count_comparisons	count_missing
popA	chr1	1	10000	0.0003778491492837536	9	170	320	85
popB	chr1	1	10000	0.00040893523488217553	9	184	354	51
popA	chr1	10001	20000	0.0005067049510407453	11	228	461	34
popB	chr1	10001	20000	0.00046229207969204235	11	208	427	68
```

`out_tajima_d.txt`：

```
pop	chromosome	window_pos_1	window_pos_2	tajima_d	no_sites	raw_pi	raw_watterson_theta	tajima_d_stdev
popA	chr1	1	10000	2.0866267869882233	9	4.843650793650793	3.3423177051634436	0.7195024514442903
popB	chr1	1	10000	1.9896866772136814	9	4.70952380952381	3.2779393676625217	0.7195024514442903
```

`out_fst.txt`：

```
pop1	pop2	chromosome	window_pos_1	window_pos_2	avg_wc_fst	no_snps
popA	popB	chr1	1	50000	-0.03926597324664331	63
popA	popB	chr1	50001	100000	-0.026026648105796968	59
```

`out.fstats`（f 统计量与 D 一行一个检验）：

```
stat	pops	estimate	se	z	n_blocks	n_snps
f4	popA,popB,popC,popD	0.4027777777777775	NA	NA	0	122
dstat	popA,popB,popC,popD	0.013145965548504073	NA	NA	0	122
```

上面这份是 122 个位点的小演示数据：整个 contig 落进一个块，块数 0，标准误写 `NA`。真实数据要跨出几十个块，`se`/`z` 才有内容。

`out.ld.decay`（先全体 `all`，多群体时每群体一组）：

```
pop	bin_start	bin_end	n_pairs	mean_r2
all	0	1000	118	0.15795037970268636
all	1000	2000	144	0.15256733855146404
all	2000	3000	155	0.10808903029910441
all	3000	4000	144	0.13556198515852427
```

## 校验与常见错误

- `--proj` 需要 `--sfs`；`--sfs` 需要 `--proj N`；投影值个数必须等于群体文件组数。
- `--block-size` 需要 f 统计量；`--pbs-normed` 需要 `--pbs`；`--ld-pairs` 需要 `--ld-decay`。
- `--pbs` 的三元组必须三个互异的已知群体名。
- f 统计量组合写法：`A,B` / `C;A,B`（也收 `C;A;B`）/ `A,B;C,D` / `A,B;C,D;X`；群体名必须都在 `--pop` 文件里。
- 窗口模式缺 `-o`、或 `-o -`，报错（窗口模式只写文件不写 stdout）；裸命令（无统计量）时输出汇总 TSV 到 stdout。
