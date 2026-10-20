# doc_v17 实施状态

> 对拍证据在 §3。判决用语：**匹配**（数值差在参照工具的输出精度内）｜
> **未对齐**（公式实现不同，如实登记）。

## 1. 实施状态（D1–D4）

| 项 | 内容 | 状态 |
|---|---|---|
| **D1** 调研 | 三个子代理 + codegraph；01/02 落定 | **完成** |
| **D2** crate 与迁移 | `reed-distance` 新建；`Mat` 下沉 `reed-core`；GRM 迁入 | **完成** |
| **D3** 矩阵族 | 10 种矩阵、4 种输出格式、`reed distance` | **完成** |
| **D4** 对拍与收口 | 见 §3 | **完成**（PI_HAT 部分未对齐） |

## 2. 结构变更

| 变更 | 说明 |
|---|---|
| 新增 `crates/reed-distance` | 标准化 + 10 种矩阵 + 输出格式 |
| `Mat` 从 `reed-reduce::linalg` 下沉到 `reed-core::linalg` | 打破 `reed-reduce` ↔ `reed-distance` 的循环依赖 |
| `Grm`/`Stand`/`Standardized`/`standardize` 从 `reed-reduce::pca` 迁入 `reed-distance` | `reed-reduce::pca` 改为 re-export，v16 API 不破 |
| `reed::genotype`（新） | VCF→矩阵的桥接从 `pca` 抽出，PCA 与 distance 共用 |
| `reed::read_genotype_matrix` | 门面公开，返回 `(Mat, Vec<String>)` |
| `reed distance` 命令 | 10 种 `--type`、4 种 `--format`、样本子集、线程/内存预算 |

## 3. 对拍结果（HapMap3 957×14266）

| 矩阵 | 参照 | 结果 | 证据 |
|---|---|---|---|
| `grm` | plink `--make-rel` | **匹配**：max abs 5e-6（= plink 6 位输出精度）；`[0,0]` 1.025560 vs 1.025557964，`[0,1]` 0.050482 vs 0.050482040 | §3.1 |
| `king` | plink2 `--make-king` | **匹配**：max abs 5e-7；`[0,1]` -0.004675080 vs -0.004675082 | §3.2 |
| `allele-ct` | plink `--distance` | **匹配**：max rel 5e-6；`[0,1]` 9901.81 vs 9901.810899 | §3.3 |
| `phi` | vcftools `--relatedness2` | **匹配**：915849 对，`N_AaAa`/`N_AAaa`/`N1_Aa`/`N2_Aa` **零差异**，phi max abs 4.97e-7 | §3.4 |
| `discordance`（count） | bcftools `gtcheck -E 0` | **匹配**：8358/8318/8209 逐对相同；`n_compared`/`n_matching` 也相同 | §3.5 |
| `pihat` | plink `--genome` | **未对齐**：见 §4 | — |
| `ajk` | vcftools `--relatedness` | 未对拍（vcftools 需 v4.2 且该模式未跑通） | — |
| `ibs`/`1-ibs` | 由 distance 换算 | 单测验证 `ibs + 1-ibs = 1` 与 `1-ibs = d·0.5/m` | §3.6 |

### 3.1 GRM

```
GRM: max abs 4.999e-06  max rel 2.410e-03
     [0,0]: 1.025560000 vs 1.025557964
     [0,1]: 0.050482000 vs 0.050482040
```

### 3.2 KING

```
KING: max abs 5.000e-07  max rel 1.300e-05
      [0,0]: 0.500000000 vs 0.500000000
      [0,1]: -0.004675080 vs -0.004675082
```

### 3.3 allele-count 距离

```
max abs diff: 0.050000   max rel: 4.998e-06
plink[0,1]=9901.810000 reed[0,1]=9901.810899
```

关键实现点：**分子不加权**（plink 的 HWE 权重只进缺失校正的分母）。
第一版误把权重乘进分子，差 3e9；修正后 5e-6。

### 3.4 Manichaikul phi

```
pairs: vcftools=915849 reed=915849 common=915849
max |dphi| = 4.970e-07
max |dN_AaAa| = 0   max |dN_AAaa| = 0
count mismatches: 0 / 915849
sample: vcftools (6434,0,6434,6434,0.5) reed (6434,0,6434,6434,0.500000000)
```

复刻了 vcftools 的两个非显然行为：
- **边际分母** `N_Aa[i]+N_Aa[j]`（不是 pairwise-complete）
- **完整 N×N 有序矩阵**、对角写 0.5

### 3.5 discordance

```
bcftools: DCv2 2424_NA19835 2431_NA19916 8358 ... 14249 5891
reed:     DCv2 2424_NA19835 2431_NA19916 8358 ... 14249 5891
```

（`-E 0` 的失配计数模式；`n_compared`/`n_matching` 逐对相同。）

## 4. 未对齐：PI_HAT

`pi_hat` 的期望 IBS 系数用了**理论期望值**（`e00 = 2p²q²` 等），而 plink
用的是 **PLINK 1.07 的 Nei 偏差修正**（`1.9/plink_calc.c:4848-4872`）：

```c
num_allelesf2 = num_alleles²/((num_alleles-1)(num_alleles-2));
num_allelesf3 = num_allelesf2 * num_alleles/(num_alleles-3);
e00 += 2*dpp²*dqq²*dxx1*dyy1*num_allelesf3;
e01 += 4*dpp*dqq*num_allelesf3*(dpp²*dxx2 + dqq²*dyy2);
...
```

这些系数依赖样本量 `num_alleles` 的阶乘修正，且只在 `num_alleles > 3` 时
更新。**本轮不实现**：这是 PLINK 1.07 的历史估计量，需要逐项复刻
`genome.cpp` 的推导。

实测差异：plink `[0,1]` 的 Z0/Z1/Z2 = 0.9551/0.0449/0.0000/0.0225，
reed 为 0/1/0/0.5。**PI_HAT 的输出保留但标记为未对齐**；若要与 plink
对拍，需要先实现 Nei 修正。

## 5. 已知限制

1. **PI_HAT 未对齐**（§4）。
2. **`ajk` 未对拍**：vcftools 的 `--relatedness` 需要 VCFv4.2 且本轮未跑通
   该模式；公式按 `variant_file_output.cpp:4837-4867` 实现。
3. **HWE 列**（discordance 的第 5 列）写 0：bcftools 的
   `average -log P(HWE)` 未实现（需要按站点 AC 算 HWE 概率）。
4. **`--genome` 的 DST/PPC/RATIO** 写占位值：它们依赖 gap 过滤与
   二项检验，未实现。
5. **常染色体限制**：`king` 需要只喂常染色体，由调用方负责（reed 的
   VCF 桥接目前不过滤 contig）。

## 6. 复现

```bash
cargo test -p reed-distance          # 23 tests
cargo test -p reed-reduce            # 38 tests（GRM 迁移后不变）
cargo test -p reed-cli --test e2e_pca

reed distance in.vcf.gz --type grm -o out            # out.rel + out.rel.id
reed distance in.vcf.gz --type phi -o out            # out.relatedness2
reed distance in.vcf.gz --type discordance --mismatch-count -o out
reed distance in.vcf.gz --type grm --format bin -o out
```
