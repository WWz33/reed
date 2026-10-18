# reed 统计族参考实现地图(Statistics Reference Map)

> 产生方式:2026-10-04,11 个并行子代理逐一审读 14 个参考项目源码(codegraph 辅助),每条 file:line 均经抽查核对(6 处代表引用全部属实)。本文是**真实记录**——只登记源码里实际存在的东西,不登记印象。
> 许可证红线:仅 MIT/Apache/BSD 源可翻译(vcftools=LGPL、plink-ng/GCTA/flashpca/SNPGenie=GPL 只读);HWE/SFS 等算法规则可从论文重写。

## 1. 按族的跨项目对照总表

标记:`✅=已实现并测试`、`🔶=已实现但见分歧注记`、`空=参考项目有而 reed 未实现`、`—=全行业参考项目都罕见/缺`。

| 族 | 统计量 | reed | vcftools(LGPL) | scikit-allel(MIT) | pixy(MIT) | sgkit(Apache) | plink2(GPL) | vcflib(MIT) | snputils(BSD) | GCTA(GPL) | flashpca(GPL) | vcfstats(MIT) | SNPGenie(GPL) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 频率 | allele counts / AF | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | — | ✅(宏) | — |
| 频率 | SFS(折叠/展开) | ✅ | — | ✅ `sf.py` | — | (S 隐含) | — | — | — | — | — | — | — |
| 多样性 | π(per-site/window) | ✅ | ✅ | ✅ | ✅ | ✅ | — | ✅ `sequenceDiversity` | — | — | — | — | ✅(逐核苷酸) |
| 多样性 | Watterson θ / Tajima's D | ✅ | ✅θ+D | ✅ | — | ✅D(无θ导出) | — | — | — | — | — | — | ✅D(R 版) |
| 多样性 | πN/πS、dN/dS | — | — | — | — | — | — | — | — | — | — | — | ✅ |
| 分化 | dxy | ✅ | — | ✅ | ✅ | ✅ | — | — | — | — | — | — | ✅(组间) |
| 分化 | Fst·Hudson/Bhatia | ✅ | — | ✅ | — | ✅ | ✅(引 allel) | — | ✅ | — | — | — | — |
| 分化 | Fst·WC84 | ✅ | ✅ | ✅ | ✅ | — | ✅(引 allel) | ✅ | ✅ | ✅ | — | — | — |
| 分化 | Fst·Nei86 | ✅ | — | — | — | ✅ | — | — | — | — | — | — | — |
| 分化 | Fst·Patterson / Tsallis | — | — | ✅Patterson | — | — | — | — | ✅Tsallis | — | — | — | — |
| 分化 | PBS | ✅ | — | ✅ | — | ✅ | — | — | — | — | — | — | — |
| HWE | HWE 精确检验 | ✅ | ✅(+卡方) | — | — | ✅ | ✅(+mid-p) | — | ✅(SNPObject) | — | — | — | — |
| HWE | Hobs/Hexp/Fis | ✅ | ✅(--het) | ✅ `hw.py` | — | ✅ | ✅(--het) | ✅ `popStats` | — | ✅(--ibc) | — | — | — |
| LD | r²(Rogers–Huff/PLINK) | ✅ | ✅ | ✅ | — | ✅ | ✅(+相位三次方程) | ✅ `vcfld` | — | ✅ | — | — | — |
| f-统计 | f2/f3/f4/D | ✅(无偏+块刀切) | — | ✅(无偏校正) | — | — | — | ✅ abba-baba | ✅(无偏校正+刀切) | — | — | — | — |
| 选择 | EHH/iHS/XPEHH/nSL | — | — | ✅ | — | — | — | ✅(iHS/xpEHH/EHH) | — | — | — | — | — |
| 单倍型 | Garud H1/H12 | — | — | ✅ | — | ✅ | — | — | — | — | — | — | — |
| 结构 | PCA | — | ✅(dgeev,有缺失bug) | ✅ | — | ✅ | ✅ | — | ✅(processing) | ✅(GRM谱) | ✅(Lanczos) | — | — |
| 结构 | GRM/kinship/relatedness | — | ✅(Yang/KING) | — | — | ✅ | ✅(KING) | — | — | ✅(GRM 正典) | — | — | — |
| 结构 | ROH | — | ✅(HMM,实验) | ✅(HMM×2) | — | — | ✅(滑窗计数) | — | — | — | — | — | — |
| 密码子 | (同 πN/πS 行) | — | — | — | — | — | — | — | — | — | — | — | ✅ |

## 2. 逐族金标准认定与定义分歧

### 2.1 Fst —— 四种估计器 × 六家实现,分歧最大的一族

| 实现 | 位置 | 语义要点 |
|---|---|---|
| **WC84 · scikit-allel** | `allel/stats/fst.py:136-241`(s²:`:202-206`;n_C:`:179`;a/b/c:`:225-241`) | `s² = Σnᵢ(pᵢ−p̄)²/(n̄(r−1))`;`n = an//2` **整数地板**;a/b/c 逐位点×逐等位基因,总和比 `Σa/Σ(a+b+c)` |
| **WC84 · vcftools** | `variant_file_output.cpp:3564-3620`(a/b/c `:3599-3620`) | 同公式;`n[i] = N_hom + 0.5·N_het`(基因拷贝数);窗口版分 bin 后 `Σa/Σall` |
| **WC84 · vcflib** | `src/wcFst.cpp:274-311` | 同公式 r=2 特例;`n̄` 用**二倍体个体数**(WC84 用基因拷贝 → 有限样本校正项 1/(n̄−1) 略有出入);逐位点输出,"负值合法" |
| **WC84 · plink2** | `plink2_misc.cc:13205-13510, 14288-14364`(注释直引 `allel/stats/fst.py`) | `b` 的 (2n̄−1)h̄/(4n̄) 写作 `(0.5−0.5/n_total)·h̄` —— 与 allel **代数恒等**;带 **block jackknife** |
| **WC84 · GCTA** | `main/popu_genet.cpp:244-257` | 与我们修复后的实现逐项一致(`d_buf = p̄q̄ − (r−1)s²/r` 等);单 SNP 输出 |
| **WC84 · snputils** | `stats/fstats.py:1538-1563`(经 `_prepare_weir_cockerham_inputs` :353) | 用观测杂合度;r=2;比值之和 jackknife |
| **Hudson · scikit-allel** | `allel/stats/fst.py:247,317-323` | `within=(mpd_A+mpd_B)/2`,`between=mpd_AB`,`num=between−within`,`den=between`,**逐位点归一化求和**;位点上任一项无定义(fill=NaN)即由 `nansum` 跳过 |
| **Hudson · plink2** | `plink2_misc.cc:13533, 14288-14364` | 同 allel(注释自证);`fst_numer = between − within`,**block jackknife** SE |
| **Hudson · snputils** | `stats/fstats.py:1384, 523-572, 1487-1493` | 加**无偏校正**:`π_x = 2p(1−p)·n/(n−1)`(有限样本),块长加权 jackknife |
| **Hudson · sgkit** | `popgen.py:341-346` + `popgen_numba_fns.py:60-79` | **分歧矩阵表述**:`1 − ((d_ii+d_jj)/2)/d_ij`(与我们/allel 的 mpd 形式代数等价);负值不 clip |
| **Nei86 · sgkit** | `popgen_numba_fns.py:90-110` | `den = d_ii + 2d_ij + d_jj`,`Fst = 1 − 2(d_ii+d_jj)/den`(自证"tskit 同款") |
| **Patterson · allel** | `fst.py:327`(num=`patterson_f2`,den=num+ĥ_A+ĥ_B,Patterson 2012 App.A) | 唯一实现 |
| **Tsallis · snputils** | `fstats.py:1564-1601`(+ `_tsallis_entropy_bernoulli` :625) | q=2 退化为 (H_T−H_S)/H_T;唯一实现 |

**裁决**:Hudson 的**正典**是 scikit-allel(plink2 源码注释自证引用它);WC84 的正典是 allel + plink2 + GCTA 三家交叉一致;reed 现有 Hudson/WC84/Nei86 已分别对齐正典(2026-10-04 修复后,`tests/unequal_groups.rs` 逐位验证)。**落库计划**:采纳 snputils 的有限样本校正 Hudson 变体(`π_x·n/(n−1)`)为第四估计器;Tsallis 列远期。

### 2.2 π / dxy —— 三种分母语义并存(本质分歧,不是谁对谁错)

| 语义 | 实现 | 位置 | 分母 |
|---|---|---|---|
| 变异位点内成对差异(vcftools 式,有偏) | vcftools windowed-pi | `variant_file_output.cpp:4065-4283`(分母 `:4267-4272`) | **单态位点假设 2N_kept 条染色体**(`N_MONOMORPHIC = window_size − N_variant_sites`,`:4252`),mask 可改;`N_site_mismatches==0` 的窗口**不输出**(`:4245`) |
| 窗口 bp 分母(allel 式) | scikit-allel / reed | `diversity.py:283-288`:`n_bases = stop−start+1`,或 `is_accessible` 计数 | 位点均值之和 ÷ bp;位级掩码是正道(v0.3) |
| 逐样本对 + 不变位点贡献(pixy 式,最严格) | pixy | `stats/summary.py:161-256`(π),`:257-358`(dxy);不变位点贡献 `invariant_contribution.pi`(`inv.comps`)入分母 | 逐群体、逐样本对;要求全位点 VCF |

**裁决**:reed 已封两种(PerPairComparable / AccessibleBp);vcftools 式不再实现(有偏且被 pixy 论文批评);pixy 语义依赖全位点输入(管线重跑后解锁)。

### 2.3 Tajima's D / θ_W —— 公式统一,守卫分歧

- 常数块 a1/a2/b1/b2/c1/c2/e1/e2 三家逐字一致:allel `diversity.py:951-962`、vcftools `variant_file_output.cpp:3953-3968`、reed `spectrum.rs`(已对拍)。
- **n 的取法分歧**:vcftools 用全局 `2·N_kept_individuals`(`:3955`);allel 用 `n = an.max()`;reed 用 complete-case 的组内 2×样本数。
- **min_sites**:allel 默认 3(`:944-946`);sgkit **无此参数**(IEEE 除零自决,`popgen.py:496-499` 注释);reed 为显式参数。
- **Tajima 变体**:vcftools 的 `pi = 2·Σp(1−p)·n/(n−1)`(等价变形);SGPGenie 另有 R 版。

### 2.4 HWE 精确检验 —— 同一算法(Wigginton 2005),精度分级

| 实现 | 位置 | 增值 |
|---|---|---|
| vcftools | `entry.cpp:18-101`(SNPHWE,作者原代码逐字移植)+ 卡方 `variant_file_output.cpp:354-365` | 基础版 |
| sgkit | `hwe.py:18-105`(numba fastmath) | 基础版 |
| reed | `estimators/hwe.rs`(自 sgkit 可译移植) | 基础版;独立枚举交叉验证 |
| **plink2** | `include/plink2_stats.cc:1841-1876`(HweLnP)+ chrX 专用 `:3973` | **mid-p**(Graffelman & Moreno 2013)+ log 域 + double-double 近平局处理 |
| snputils | `snp/genobj/snpobj.py:3463`(hwe_pvalue) | 基础版 |

**裁决**:reed 已有基础版(正典);**mid-p 选项列为 v0.3**(plink2 证明其价值,规则可从 Graffelman & Moreno 2013 重写)。

### 2.5 LD r² —— 计算一致,缺失语义分歧

| 实现 | 位置 | 缺失语义 |
|---|---|---|
| scikit-allel | `ld.py:13`(cython `opt/stats.pyx:90`) | **缺失不剔除**:-1 直接参与相关计算(文档未承诺剔除) |
| reed | `estimators/ld.rs` | **成对剔除**(pairwise-complete) |
| vcftools | `calc_geno_r2` `variant_file_output.cpp:1023-1075` | 成对剔除 |
| plink2 | `plink2_ld.cc:7569-7601`(整数精确);相位未知 het-het 对解三次方程 `:5538` | 成对剔除;127 位精确乘积防 FP 漂移 |
| GCTA | `src/LD.cpp:58-69`(dsyrk) | **均值填充** |

**裁决**:reed 的成对剔除最保守;若要与 allel 数值对拍需提供"含缺失模式"选项(v0.3)。plink2 的整数精确 r² 与相位三次方程是 v0.3 值得封的变体。

### 2.6 f2/f3/f4/D —— reed 的朴素实现需要升级

| 实现 | 位置 | 校正 |
|---|---|---|
| scikit-allel | `stats/admixture.py`:`h_hat`:11、`patterson_f2`:43、`f3`:90、`patterson_d`:147;jackknife `average_patterson_d`:373 | **Patterson 2012 无偏 ĥ 校正**(f2 减 ĥ 项;f4 共享群体协方差权重) |
| snputils | `fstats.py`:f2:668、f3:812、f4:945、D:1087、f4_ratio:1222 | 同类校正 + **五种 jackknife 助手**(:85/:125/:175/:188)+ `BlockJackknifeResult`(est/se/z/p) |
| vcflib | `abba-baba.cpp:27-46` | 仅位点分类计数 |
| reed | `stats/writers/fstats.rs` + `reed-stats::patterson` | **无偏校正 + bp 块刀切**(`se/z/n_blocks` 列) |

**裁决**:正典 = snputils(BSD-3,五助手+数据类可直接移植)+ allel admixture.py(MIT,Patterson 校正)。reed v0.3:加 `apply_correction` 旗标(f2:−p(1−p)/(n−1)/群;f3:−目标项;f4:重叠权重)+ block jackknife(块长加权,ADMIXTOOLS 语义)。

### 2.7 PBS —— clip 上界的实现分歧

- scikit-allel:`selection.py:1322-1323` clip 到 **[0, 0.99999]**;两定义 Yi2010(`:1331`)/PBSn1(Malaspinas 2016,`:1335` 归一化)。
- sgkit:`popgen.py:583-586` clip 到 **[0, 1−ε]**;`_pbs_cohorts` numba 同款归一化(`popgen_numba_fns.py:150-152`)。两家测试对拍通过是因为典型 Fst 远离上界。
- reed:已对齐 **allel 的 0.99999**(2026-10-04 审计修复)。

### 2.8 PCA / GRM / 相关性

| 实现 | 位置 | 方法 |
|---|---|---|
| flashpca | `randompca.cpp:121-218` + `svdwide.cpp:21-68` | **隐式 Lanczos**(Spectra 对 `y=XXᵀx` 算子),流式 .bed 分块;binom2 scaler(`util.cpp:111-134`,Price 2006);缺失=标准化空间均值填充 |
| plink2 | `plink2_matrix_calc.cc:7179`(dsyevr)/`:7289`(Halko 2011 随机,移植自 EIGENSOFT kjg_fpca) | 精确 GRM+dsyevr / 随机 SVD 双模式;2-bit 流式块 |
| GCTA | `src/GRM.cpp:933-1038,1356-1520`(dsyrk+容斥 N) | GRM 三实现;PCA=GRM 谱分解(`main/grm.cpp:676-682`) |
| scikit-allel | `decomposition.py:8,123` + `preprocessing.py:102`(PattersonScaler) | sklearn SVD / randomized_svd |
| sgkit | `pca.py:138`(tsqr/randomized)+ `pc_relate.py`(Conomos 2016)+ `grm.py:293-320`(VanRaden/Endelman-Jannink 开关) | 缺失不支持(hard error) |
| vcftools | `variant_file_output.cpp:4871-5042` | dgeev + **缺失数据行错位 bug**(仅全无缺失可用) |

### 2.9 ROH —— 两大家族

- HMM 族:vcftools LROH(`:4375-4644`,forward-backward,实验性)、scikit-allel `roh.py:10`(multinomial HMM)/`:149`(Poisson HMM)。
- 计数族:plink2 `plink2_misc.cc:11520-11680`(滑窗 het/missing 阈值)。

### 2.10 其他

- **KING/relatedness**:plink2 `plink2_matrix_calc.cc:1558-1580`(`φ = 0.5 − (4·ibs0 + het1hom2 + het2hom1)/(4·smaller_het_ct)`);vcftools Yang(`:4842-4858`)/Manichaikul(`:4752`,已知与原文有偏差);sgkit pc_relate(Conomos 2016)。
- **Garud H**:sgkit `popgen.py:610-632` + allel `selection.py:797`(DJBX33A 单倍型哈希)。
- **SFS**:allel `sf.py:34-99`(folded:`mac = amin(ac,axis=1)`,bincount;fold_sfs 折半合并)。联合(二维)谱:dadi `Spectrum.from_data_dict` 的逐位点外积累加 + 全局折叠(按 `i+j` 对 `(n1+n2)/2` 镜像,对角 0.5 分裂),easySFS `dadi_twoD_sfs_combinations` 是文件布局参照。
- **πN/πS/dN/dS**:SNPGenie 唯一参考(Nei-Gojobori 分数位点表 `within_group.pl:1625-1721` + 路径平均差异 `:1728+`;GPL——规则从 Nelson et al. 2015 论文重写)。
- **传输层**:htslib `kfunc.c`(lgamma/erfc/betai/Fisher)、`vcfutils.h:91`(bcf_calc_ac)、`vcf.h:245-258`(bcf1_t 惰性)、CSI bin 数学 `hts.h:899`;noodles 惰性 Series(`variant/record/samples/series.rs:13`);rust-htslib **Record 构造即 `bcf_unpack(BCF_UN_ALL)`**(`src/bcf/record.rs:186`)——我们基准结果的机制级解释。
- **vcfstats**:`formula.py`(Y~X Lark 语法)、`macros.py`(三类装饰器注册)、`instance.py`(单遍流 + 末端聚合)——统计即插拔定义的架构范本;无窗口机制。

## 3. 金标准总认定(reed 的对拍优先级)

1. **Fst·Hudson**:scikit-allel `hudson_fst`(plink2 源码注释自证引用)——reed 已对拍。
2. **Fst·WC84**:allel + plink2 + GCTA 三家一致(相互独立)——reed 修复后已对拍。
3. **π/dxy**:scikit-allel(bp 语义)+ pixy(逐样本对语义)双正典——reed 双估计器已对拍。
4. **Tajima/θ_W**:allel 常数块(教科书逐字)——已对拍。
5. **HWE**:Wigginton 算法三家同源;**mid-p 升级版看 plink2**。
6. **PBS**:allel(normed 两版)——已对拍;注意与 sgkit 的 clip 上界分歧。
7. **f-统计**:snputils(校正+刀切,BSD-3 直接移植)> allel admixture.py —— reed 已落地无偏校正 + bp 块刀切(与 snputils 同类)。
8. **LD**:allel cython 核;plink2 的整数精确/相位三次方程为速度与相位变体参考。
9. **GRM**:GCTA(Yang 2010 正典);PCA:Halko 2011(plink2/flashpca 两家)。
10. **πN/πS**:SNPGenie 论文规则(未来密码子模块)。

## 4. 落库计划(不改代码,按优先级)

1. **v0.2**:spectrum 族接 availability 门(θ_W 拒绝 VO(Unknown));六种窗口级结果补 contig/start/end;`DetectedFacts` 扫描。
2. **v0.3-a(f-统计升级)**:Patterson 无偏校正(`apply_correction` 旗标,allel/snputils 语义)+ block jackknife(snputils 五助手移植,BSD-3)。
3. **v0.3-b(HWE mid-p)**:按 Graffelman & Moreno 2013 重写(参考 plink2 行为)。
4. **v0.3-c(LD 变体)**:缺失语义选项(对拍 allel 需含 -1 模式);整数精确 r² 思路(PLINK 式)。
5. **v0.3-d**:Fst·Hudson 有限样本校正变体(snputils 的 π·n/(n−1));WC84 多群体(r>2)泛化。
- **v0.4.0（v14，已落地）**：`reed stats --hist LIST` 把 QUAL/DP/DP-site/GQ/
  IDD/AF 直方图与任意数值 `INFO/`/`FMT/` 标签输出为 `PREFIX.<m>.hist` TSV
  （整数指标附精确分位数；对拍 `gold/popgen/run_hist.sh`）；逐等位
  singleton/doubleton 表 `--singletons`（vcftools 布局，对拍
  `gold/popgen/run_singletons.sh`）；`--pop` 分层 QC——`.afreq.pop`/
  `.vmiss.pop` 每群体一行 + `.singletons` POP 列（分层值对拍按群体切样本
  的 plink2 跑法，`gold/popgen/run_popqc.sh`）。
6. **v1.0 方向**:GRM/PCA(Halko 随机 + Patterson scaler)、KING、ROH(双家族)、EHH/iHS(相位层)、πN/πS(注释层)。
