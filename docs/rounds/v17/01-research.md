# doc_v17 / 01 调研台账：三家工具的个体×个体矩阵

> 三个子代理分别用 `codegraph` 索引 + 源码精读，逐条给出 file:line 与精确
> 公式。判决用语：**源码实证**（读到并引用）｜**数值实测**（构造数据跑参照
> 工具）。参照源码：
> - plink-ng 1.9 + 2.0：`/sri/home/sri2025201067/reed/plink-ng`（已索引）
> - vcftools 0.1.18：`/sri/home/sri2025201067/reed/vcftools`（已索引）
> - bcftools：`/sri/home/sri2025201067/reed/bcftools`（本轮新建索引，149 文件 / 4229 节点）

## 0. 矩阵总表

| 矩阵 | 参照选项 | 输出文件 | 公式核心 |
|---|---|---|---|
| 等位基因计数距离 | plink `--distance` | `.dist` | Σ\|d_i − d_j\|（两两非缺失） |
| IBS | plink `--distance ibs` | `.mibs` | 1 − dist·0.5/m |
| 1−IBS | plink `--distance 1-ibs` | `.mdist` | dist·0.5/m |
| GRM（关系矩阵） | plink `--make-rel`、`--make-grm-*` | `.rel`/`.grm*` | Σ x_i x_j / N_valid |
| 协方差矩阵 | plink `--make-rel cov` | `.rel` | 同上但只中心化 |
| KING 稳健亲缘 | plink2 `--make-king` | `.king` | 0.5 − (4·IBS0+H1H2+H2H1)/(4·(HETHET+min(H1H2,H2H1))) |
| PI_HAT / IBD | plink `--genome` | `.genome` | 矩估计 Z0/Z1/Z2，PI_HAT=Z1/2+Z2 |
| Manichaikul φ | vcftools `--relatedness2` | `.relatedness2` | (N_AaAa − 2·N_AAaa)/(N_Aa[i]+N_Aa[j]) |
| Yang Ajk | vcftools `--relatedness` | `.relatedness` | Σ(x_i−2p)(x_j−2p)/(2p(1−p)) / M_ij |
| 不一致度 | bcftools `gtcheck` | stdout `.gtcheck` | 失配数 或 Σ min_G[−log P_a(G) − log P_b(G)] |

---

## 1. plink `--distance`（等位基因计数距离）

**入口**：`1.9/plink_calc.c:7601 calc_distance()`；写出
`distance_d_write()`（`:3756`）→ `distance_d_write_tri_emitn()`（`:3356`）。

### 1.1 分子

基因型 2-bit 编码经 XOR `01` 变换后为 `00=缺失, 01=hom1, 10=hom2, 11=het`；
掩码 `((g|(g>>1))&01)*11` 在非缺失时为 `11`（注释 `:7850-7891`）。

累加 `incr_dists_i()`（`:973-1024`）：对每对样本做
`popcount((geno_i ^ geno_j) & mask_i & mask_j)`（`:1014`/`:1019`）。
2-bit 编码下这恰好是**剂量差** |d_i − d_j| ∈ {0,1,2}（het vs hom = 1 位，
相反纯合 = 2 位，相同 = 0）。任一方缺失的位点贡献 0。

```
A_ij = Σ_v [两方均非缺失] |dosage_i(v) − dosage_j(v)|
```

`--distance` 默认类型即 allele-ct（`:6186-6188`）。文本输出直接写该值，
不除以位点数（`:3368`/`:3375`）。

### 1.2 缺失校正

默认（非 `flat-missing`）按 HWE 权重校正。每个位点
`p = set_allele_freqs[v]`：

```
dxx = p(1−p)(p² − p + 1)                    // :7777-7787
dyy = Σ dxx
dyy = (2^32 − marker_ct) / dyy              // :7803-7807
uii = round(weight_v * dyy)                 // 每变异权重
marker_weight_sum = Σ uii                   // :7810
```

每对：

```
g_dists[i,j] = marker_weight_sum_d
             / (marker_weight_sum − sample_missing[i] − sample_missing[j]
                + missing_both_weight[i,j])
             * g_idists[i,j]                // :8183（无 --distance-wts）
```

分母 = 两方均非缺失位点的权重和。无缺失时退化为
`Σ|剂量差|`。

### 1.3 修饰符（`:6162-6188`）

| 修饰符 | 效果 | 输出 |
|---|---|---|
| 默认 / `allele-ct` | 上述距离 | `.dist` |
| `ibs` | `1 − dist·0.5/m`（IBS 比例） | `.mibs` |
| `1-ibs` | `dist·0.5/m`（遗传距离比例） | `.mdist` |
| `flat-missing` | 跳过 HWE 权重，用普通计数：`(m/(m−miss_i−miss_j+miss_both))·idists` | 同 `.dist` |

`--distance-matrix`（已弃用）= `1-ibs flat-missing square`。

### 1.4 输出布局

下三角、行主序：第 1 行只有 `d(1,2)`；第 2 行 `d(1,3) d(2,3)`；……
（帮助 `1.9/plink_help.c:708-711`；emitter `:3366-3384`）。

---

## 2. plink `--make-rel` / `--make-grm-*`（GRM）

**入口**：`2.0/plink2_matrix_calc.cc:5534 CalcGrm()`。

### 2.1 标准化（`:4179 ExpandCenteredVarmaj`）

```
alt_freq = 1 − ref_freq
variance = 2 · ref_freq · alt_freq            // :4183
inv_stdev = 1 / sqrt(variance)                // :4210
x_s(v) = (dosage_alt(s,v) − 2·alt_freq) · inv_stdev   // :4224
```

缺失（及无硬 call 的剂量）取 `missing_val = 0.0`，即中心化后的均值。
多等位：`variance = Σ p_a(1−p_a)`（`:4267`），`inv_stdev = (1/√2)/sqrt(variance)`
（`:4376`）。

GRM 累加为 BLAS 自叉积（`:4635`/`:4662`）：

```
G_ij = Σ_v x_i(v) · x_j(v)      （除以有效位点数之前）
```

### 2.2 `cov` 修饰符（`:5670`，flag `kfGrmCov`）

```
inv_stdev = is_haploid ? 0.5 : 1.0            // :4222
x_s(v) = dosage_alt − 2·alt_freq              // 仅中心化
```

### 2.3 缺失策略与分母

| 模式 | 分母 | 代码 |
|---|---|---|
| `meanimpute` | `variant_ct`（全体位点数） | `:5758-5763` |
| 默认（pairwise-complete） | `N_valid(i,j) = variant_ct − miss_i − miss_j + miss_dbl(i,j)`；对角 `variant_ct − miss_i` | `:5746-5755` |

`miss_dbl` 由 `CalcMissingMatrix()`（`:4744`）与
`CalcDblMissingThread`（`:4681-4742`）累计。

**没有 N−1 或 N 的除数**——除数是有效变异数（GCTA 约定）。

### 2.4 GRM 导出格式

- `.grm.bin`：单精度下三角，第 r 行 r+1 个浮点（`:5940-5948`）。
- `.grm.N.bin`：对应计数（`:5953-5999`）；无缺失时为常数 `variant_ct`。
- `.grm` 列表：每行 `row+1 col+1 obs_ct rel_coeff`（`:6059-6085`）。
- `--make-grm-sparse <thr>`：仅 `rel_coeff >= thr`，0-based（`:6041-6058`）。
- `.id` 追加（`:6105-6119`）。

---

## 3. plink2 `--make-king`（KING 稳健亲缘）

**入口**：`2.0/plink2_matrix_calc.cc:1932 CalcKing()`。

### 3.1 公式（`:1569-1576 ComputeKinship`）

```cpp
const uint32_t ibs0_ct    = king_counts_entry[kKingOffsetIbs0] + singleton_hom1_ct + singleton_hom2_ct;
const uint32_t hethet_ct  = king_counts_entry[kKingOffsetHethet];
const uint32_t het2hom1_ct= king_counts_entry[kKingOffsetHet2Hom1] + singleton_het2_ct;
const uint32_t het1hom2_ct= king_counts_entry[kKingOffsetHet1Hom2] + singleton_het1_ct;
const intptr_t smaller_het_ct = hethet_ct + MINV(het1hom2_ct, het2hom1_ct);
return 0.5 - ( (double)(4*ibs0_ct + het1hom2_ct + het2hom1_ct)
               / (double)(4 * smaller_het_ct) );
```

即

```
kinship_ij = 0.5 − (4·IBS0 + HET1_HOM2 + HET2_HOM1)
                  / (4·(HETHET + min(HET1_HOM2, HET2_HOM1)))
```

代数上等价于 Manichaikul 2010 的
`(HETHET − 2·IBS0)/(het_i + het_j)`（当 `HET1_HOM2 = HET2_HOM1` 时精确；
一般用 `min` 的稳健对称形式）。`--make-king-table` 内联同一式
（`:2592-2593`）。

### 3.2 计数累加（`:1349-1429 IncrKing`）

```
het1 = hom1 & ~ref2het1;   het2 = hom2 & ~ref2het2;
agg_ibs0     = (ref2het1 ^ ref2het2) & (hom1 & hom2);   // 相反纯合
agg_hethet   = het1 & het2;
agg_het2hom1 = hom1 & het2;
agg_het1hom2 = hom2 & het1;
```

缺失被 `hom`/`ref2het` 拆分自然屏蔽，只计两方均非缺失的位点。

### 3.3 约束

- **只用常染色体**：`CalcKing()` 调 `CountNonAutosomalVariants` +
  `ExcludeNonAutosomalVariants`（`:1981-1996`）；帮助 `2.0/plink2_help.cc:1561`。
- **无 within-family 分支**：帮助明说"pedigree information is currently
  ignored; the between-family estimator is used for all pairs"（`:1562-1563`）。

### 3.4 `--make-king-table` 列

默认 `maybefid,id,maybesid,nsnp,hethet,ibs0,kinship`
（`2.0/plink2_help.cc:1602`，`kfKingColDefault`）。`NSNP` =
`het1hom2 + het2hom1 + homhom + hethet`（`:2611`）。

---

## 4. plink `--genome`（IBD / PI_HAT）

**入口**：`1.9/plink_calc.c:4516 calc_genome()`；发射
`calc_genome_emitn()`（`:4232`）。

### 4.1 矩估计（`:4381-4385`）

`nn` = 两方均非缺失的位点数：

```
nn  = marker_ct − missing_i − missing_j + missing_dbl(i,j)      // :4381
oo  = nn − IBS1 − IBS0                                          // :4382 (IBS2)
Z0  = IBS0 / (e00 · nn)                                         // :4383
Z1  = (IBS1 − Z0·e01·nn) / (e11 · nn)                           // :4384
Z2  = (oo − nn·(Z0·e02 + Z1·e12)) / nn                          // :4385
```

`e00,e01,e02,e11,e12` 是期望 IBS|IBD 系数，按位点平均后归一化
（`:4848-4868`，`:4890-4895`）。

### 4.2 PI_HAT 与派生列

```
PI_HAT = Z1·0.5 + Z2                                            // :4418
DST    = 1 − (IBS1 + 2·IBS0) / (2·nn)                           // :4454
PPC    = normdist((HETHET/(HOMHOM+HETHET) − 2/3) / sqrt(2/9/(HOMHOM+HETHET)))  // :4450-4458
RATIO  = HETHET / HOMHOM（HOMHOM=0 时 NA）                       // :4460-4463
```

### 4.3 裁剪与修饰符

- 默认裁剪（`:4386-4416`）：任一 Z > 1 → 置 1、其余置 0；任一 Z < 0 →
  置 0、其余按比例缩放到和为 1。
- `unbounded`（flag `GENOME_IBD_UNBOUNDED`）：跳过整个裁剪块（`:4386`）。
- `nudge`（`:4423-4427`）：当 `PI_HAT² < Z2` 时
  `Z0=(1−PI_HAT)², Z1=2·PI_HAT·(1−PI_HAT), Z2=PI_HAT²`。

### 4.4 输出列（`:4993`/`:4995`）

```
FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO
[full: IBS0 IBS1 IBS2 HOMHOM HETHET]
```

只用常染色体（`:4669-4673`），拒绝单倍体基因组（`:4592-4595`）。

---

## 5. vcftools `--relatedness2`（Manichaikul φ）

**入口**：`src/cpp/variant_file_output.cpp:4646-4759
output_indv_relatedness_Manichaikul()`。

### 5.1 计数（`:4706-4741` 逐字）

```cpp
for ui:
  get_indv_GENOTYPE_ids(ui, geno_id);
  if (geno_id.first != geno_id.second && first>=0 && second>=0) N_Aa[ui]++;
  for uj:
    get_indv_GENOTYPE_ids(uj, geno_id2);
    if (ui is het) {
      if (uj is het) N_AaAa[ui][uj]++;
    }
    if (ui is hom) {
      if (uj is hom && geno_id.first != geno_id2.first) N_AAaa[ui][uj]++;
    }
```

- `N_AaAa`：**双方均为杂合**。
- `N_AAaa`：**双方均为纯合且基因型不同**（0/0 vs 1/1）。两个方向都累加，
  故对称。
- `N_Aa[ui]`：**边际**杂合计数，与外层 `uj` 无关。
- 缺失：跳过（`(-1,-1)` 既不满足 het 也不满足 hom 测试）。

### 5.2 公式（`:4752`）

```
phi[ui][uj] = (N_AaAa[ui][uj] − 2·N_AAaa[ui][uj]) / (N_Aa[ui] + N_Aa[uj])
```

### 5.3 输出与语义要点

- 头：`INDV1 INDV2 N_AaAa N_AAaa N1_Aa N2_Aa RELATEDNESS_PHI`（`:4668`）。
- **完整 N×N 有序矩阵**（`ui`、`uj` 都跑 0..N−1），**对角也写**，
  对角值 = `N_Aa/(2·N_Aa) = 0.5`（不是 1）。
- 每个无序对输出两次（对称）。
- **不是 pairwise-complete**：分子是成对的，分母是边际的
  `N_Aa[ui] + N_Aa[uj]`，两者位点集不同。
- 分母为 0 时不做保护（得 inf/nan）。
- 只计双等位（`:4693-4697`）与二倍体（`:4699-4704`）。
- **不检查 `include_genotype`**（与 `--relatedness` 不同），所以通过
  GQ/DP 过滤但非字面缺失的基因型仍被计入——这是与 `--relatedness` 的
  真实行为分歧。

---

## 6. vcftools `--relatedness`（Yang Ajk）

**入口**：`variant_file_output.cpp:4761-4869 output_indv_relatedness_Yang()`。

- 头：`INDV1 INDV2 RELATEDNESS_AJK`（`:4785`）。
- 过滤：双等位、二倍体，且要求 `0 < alt_freq < 1`（`:4821-4825`）。
- 剂量 `x[ui] = geno_id.first + geno_id.second`（缺失 −2，`:4827-4835`）。
- 累加（`:4837-4851`）：

```
div = 1/(2·p(1−p))
对角:  Ajk[i][i] += (x_i² − (1+2p)x_i + 2p²)·div;   N_sites[i][i]++
非对角: Ajk[i][j] += (x_i − 2p)(x_j − 2p)·div;      N_sites[i][j]++   (i<j)
```

- 收尾（`:4854-4867`）：

```
Ajk[i][i] = 1 + Ajk[i][i]/N_sites[i][i]
Ajk[i][j] = Ajk[i][j]/N_sites[i][j]        // 仅 i<j
```

即 `Ajk_ij = (1/M_ij) Σ (x_i−2p)(x_j−2p)/(2p(1−p))`，**pairwise-complete**。

与 relatedness2 的五点差异：加权 vs 计数；上三角 vs 全矩阵；pairwise
分母 vs 边际分母；尊重 `include_genotype` vs 不尊重；需要频率 vs 不需要。

---

## 7. bcftools `gtcheck`（不一致度）

**入口**：`vcfgtcheck.c:1238 main_vcfgtcheck()`；核心 `process_line`
（`:579-815`），报告 `report()`（`:872-1096`）。

> **重要更正**：旧版 gtcheck 的 `n_checked`/`n_mismatch`/`n_hethom` 等
> 数组**在本版本已不存在**；当前是 DCv2 实现，字段为
> `nused[2][2]`、`pdiff`、`ndiff`、`ncnt`、`nmatch`、`hwe_prob`
> （`:69-75`）。旧手册里的 `SM` 行也不存在（NEWS 记载 1.19 换成 DCv2）。

### 7.1 剂量掩码

- `gt_to_dsg`（`:481-486`）：非缺失时 `1 << (allele0 + allele1)` →
  `1=hom-ref, 2=het, 4=hom-alt`；缺失 = 0。
- `pl_to_dsg`（`:487-497`）：取 `min(PL0,PL1,PL2)`，**所有等于最小值**的
  索引都置位（可为 3、5、6、7）。

### 7.2 配对比较（`:758-770`，无 `-g` 时 `:795-812`）

```
match = dsg_i & dsg_j;
if (!match) ndiff[idx]++;
else if (calc_hwe_prob) { hwe_prob[idx] += hwe_dsg[match]; nmatch[idx]++; }
ncnt[idx]++;
```

- 无 `-g`：`for j < i` 的**严格下三角**，每无序对一次，无自配对。
- 有 `-g`：`nqry × ngt` 全矩形。

### 7.3 不一致度（两模式）

| 模式 | 值 | 打印 |
|---|---|---|
| `-E 0`（或 `--distinctive-sites` 强制） | `ndiff`（失配计数） | `%u` |
| 默认 `-E 40` | `pdiff`（概率分数） | `%e` |

**概率分数**：`pl2prob[i] = 10^(−0.1·i)`（`:428`）。PL 路径
（`pl_to_prob` `:509-525`）：

```
prob[k] = pl2prob[PL_k]（PL>255 饱和到 255）
prob[k] /= Σ prob
prob[k] = −log(prob[k])
```

GT 路径（`gt_to_prob` `:498-508`）用预计算表，`eprob = 10^(−0.1·gt_err)`：

```
hom-ref: {0,        −log eprob,   −2log eprob}
het:     {−log eprob, 0,          −log eprob}
hom-alt: {−2log eprob, −log eprob, 0}
```

每对：`pdiff += min_G(qry_prob[G] + gt_prob[G])`（`:716-721`/`:799-802`），
即 **min of sums**，不是单个联合概率的 −log。

### 7.4 HWE 列

每站点（`:613-659`）：

```
af = AC_alt / (AC_ref + AC_alt)（AC 总和为 0 时取 1e-6）
hwe[0] = −log((1−af)²)
hwe[1] = −log(2·af·(1−af))
hwe[2] = −log(af²)
hwe_dsg[mask] = min over set bits k of hwe[k]
```

匹配站点累加 `hwe_dsg[dsg_i & dsg_j]`，输出
`hwe_prob/nmatch`（平均 −log P(HWE)；越大越像同一人）。

### 7.5 输出

```
#DCv2  [2]Query Sample  [3]Genotyped Sample  [4]Discordance  [5]Average -log P(HWE)  [6]Number of sites compared  [7]Number of matching genotypes
DCv2 <qry> <gt> <disc> <avg_hwe> <ncnt> <nmatch>
```

前面有 `INFO` 计数器（`sites-compared`、`sites-skipped-*`、
`sites-used-*-vs-*`，`:875-886`）。

### 7.6 站点过滤（`is_input_okay` `:1098-1182`、`set_data` `:527-578`）

| 计数器 | 条件 |
|---|---|
| `multiallelic` | `n_allele > 2` |
| `monoallelic` | 全部为 REF 且**未给 `--keep-refs`** |
| `no-data` | GT 与 PL 都取不到 |
| `GT-not-diploid` / `PL-not-diploid` | 返回长度不等于 `2·n` / `3·n` |
| `filtering-expression` | `-i/-e` 失败 |
| `no-match` | 两文件记录不重叠（仅 `-g` 模式） |

---

## 8. 结论：本轮要实现什么

**做**（纯 VCF 下游、公式明确、可对拍）：

| 矩阵 | 对齐对象 | 备注 |
|---|---|---|
| `allele-ct` 距离 | plink `--distance` | 含 HWE 缺失校正 |
| `ibs` / `1-ibs` | plink `--distance ibs/1-ibs` | 由距离换算 |
| `grm` / `cov` | plink `--make-rel [cov]`、GCTA | 含 pairwise/meanimpute |
| `king` | plink2 `--make-king` | 常染色体、between-family |
| `pihat` | plink `--genome` | 矩估计 + 裁剪 + nudge |
| `phi`（Manichaikul） | vcftools `--relatedness2` | 边际分母 |
| `ajk`（Yang） | vcftools `--relatedness` | pairwise 分母 |
| `discordance` | bcftools `gtcheck` | 失配计数 + 概率分数 |

**不做**（本轮登记）：

| 项 | 理由 |
|---|---|
| `--genome` 的 PPC/IBD 片段 | 需 HMM 与 gap 过滤，另一轮 |
| `--distinctive-sites` | 是逐站点输出，不是矩阵 |
| `--make-grm-sparse` | 稀疏输出，待需求 |
| `--distance-wts exp=` | 加权距离，参数空间大 |
| vcftools `--geno-chisq` | 是 SNP×SNP 卡方，不是个体×个体 |
| bcftools 的 relatedness/kinship | 该版本不存在（已核实） |
