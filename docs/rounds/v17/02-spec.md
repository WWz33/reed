# doc_v17 / 02 规范：`reed-distance`

## 1. 目标与边界

一个 crate 承载**所有个体×个体矩阵**：标准化基因型 → 各种统计量 → 统一
输出。PCA 不再自己建 GRM，而是调用本 crate。

**不依赖** `reed-io`/`reed-engine`（纯计算，输入是 `Mat`）；VCF 桥接在
`reed` 门面层（复用 `reed pca` 已有的 `PcaJob`）。

## 2. 模块

```
reed-distance
├── mat        与 reed-reduce 共享的 Mat（或复用 reed-reduce::Mat）
├── standard   Stand、Standardized、standardize()
├── grm        Grm（从 reed-reduce::pca 迁移）、mean_impute/pairwise
├── matrices   allele_ct_distance、ibs、grm_rel、king、pi_hat、phi、ajk、discordance
├── output     MatrixFormat（triangle/square/square0/bin）、write_matrix()
└── lib
```

**迁移决定**：`reed-reduce::pca` 的 `standardize`/`Standardized`/`Grm`
移入 `reed-distance::standard`/`grm`；`reed-reduce` 改为依赖
`reed-distance` 并 re-export（保持 v16 的 API 不破）。

## 3. 类型

### 3.1 `Stand`

```rust
pub enum Stand { Binom2, Binom }
```

与 v16 相同（binom2 = sqrt(2p(1-p))，binom = sqrt(p(1-p))）。

### 3.2 `Standardized`

```rust
pub struct Standardized {
    pub w: Mat,        // n_samples x n_sites，缺失置 0
    pub mean: Vec<f64>,
    pub sd: Vec<f64>,
    pub called: Vec<bool>,  // n_samples * n_sites，是否非缺失
    pub n_sites: usize,
}
```

`called` 是新增的：距离与计数类矩阵需要"哪些位点两方都有"。

### 3.3 `Matrix`

```rust
pub struct Matrix {
    pub values: Vec<f64>,  // 完整 n x n（对称），行主序
    pub n: usize,
    pub kind: MatrixKind,
    pub ids: Vec<String>,
}

pub enum MatrixKind {
    AlleleCountDistance,
    Ibs,
    OneMinusIbs,
    Grm,
    Covariance,
    King,
    PiHat,
    Phi,
    Ajk,
    Discordance,
}
```

内部**总是存完整对称矩阵**，输出时按格式裁剪。理由：下三角输出只是
渲染问题，存完整矩阵让各算法实现简单且不易错。

### 3.4 `MatrixFormat`

```rust
pub enum MatrixFormat {
    Triangle,   // plink 默认：下三角逐行
    Square,     // 完整方阵
    Square0,    // 完整方阵，右上填 0
    Bin,        // 双精度二进制方阵（R 可读）
}
```

`Triangle` 的行布局严格按 plink：第 1 行只有 `d(1,2)`，第 2 行
`d(1,3) d(2,3)`，……（01 §1.4）。

## 4. 各矩阵的规范

### 4.1 `allele_ct_distance`

```
A_ij = Σ_v [两方 called] |dosage_i(v) − dosage_j(v)|
```

缺失校正（默认开启，对齐 plink）：

```
w_v = p(1−p)(p²−p+1),  归一化到 Σ round(w_v·dyy) = 2^32 − m
d_ij = A_ij · Σw / (Σw − miss_i − miss_j + miss_both_ij)
```

`flat_missing = true` 时用普通计数
`(m/(m−miss_i−miss_j+miss_both_ij))·A_ij`。

### 4.2 `ibs` / `one_minus_ibs`

```
ibs_ij     = 1 − d_ij · 0.5 / m
one_minus  = d_ij · 0.5 / m
```

其中 `m` 是**总位点数**（plink 语义），不是有效位点数。

### 4.3 `grm` / `covariance`

```
grm_ij = Σ_v w_i(v)·w_j(v) / N_valid(i,j)
N_valid(i,j) = m − miss_i − miss_j + miss_both_ij     （pairwise）
N_valid(i,j) = m                                       （meanimpute）
cov: w = dosage − 2p（不除 sd）
```

**契约**：与 plink2 `--make-rel` 的差 < 1e-9（同数据同集合）。

### 4.4 `king`

```
king_ij = 0.5 − (4·IBS0 + H1H2 + H2H1) / (4·(HETHET + min(H1H2, H2H1)))
```

其中（只计两方 called 的位点）：
- `IBS0` = 相反纯合的对数
- `HETHET` = 双方杂合
- `H1H2` = i 杂合 & j 纯合（alt）
- `H2H1` = i 纯合 & j 杂合

**常染色体限制**：本 crate 不知道 contig，由调用方过滤（CLI 层负责）。

### 4.5 `pi_hat`

按 plink 矩估计（01 §4）：

```
nn  = m − miss_i − miss_j + miss_both_ij
IBS1, IBS0 由位点计数
Z0 = IBS0/(e00·nn)
Z1 = (IBS1 − Z0·e01·nn)/(e11·nn)
Z2 = (oo − nn·(Z0·e02 + Z1·e12))/nn
PI_HAT = Z1/2 + Z2
```

`e00..e12` 按位点平均后归一化。默认裁剪到 [0,1]；`unbounded` 关闭；
`nudge` 在 `PI_HAT² < Z2` 时重设。

### 4.6 `phi`（Manichaikul，vcftools 口径）

```
N_AaAa[i][j] = # 双方均杂合
N_AAaa[i][j] = # 双方均纯合且基因型不同
N_Aa[i]      = 个体 i 的边际杂合数
phi = (N_AaAa − 2·N_AAaa) / (N_Aa[i] + N_Aa[j])
```

**必须复刻 vcftools 的边际分母**（不是 pairwise），否则对不上。分母为 0
时写 `NaN`（vcftools 写 inf/nan，reed 显式写 NaN 并在 stderr 警告）。

### 4.7 `ajk`（Yang，vcftools 口径）

```
Ajk_ij = (1/M_ij) Σ [ (x_i−2p)(x_j−2p) / (2p(1−p)) ]
Ajk_ii = 1 + (1/M_ii) Σ [ (x_i²−(1+2p)x_i+2p²) / (2p(1−p)) ]
```

`M_ij` = 两方 called 的位点数（pairwise）。要求 `0 < p < 1`。

### 4.8 `discordance`（bcftools gtcheck 口径）

两种模式：

```
mode=count:  discordance = # 失配（dsg_i & dsg_j == 0）
mode=prob:   discordance = Σ_sites min_G [ −log P_i(G) − log P_j(G) ]
```

其中 GT 的概率表按 `eprob = 10^(−0.1·E)` 生成（E 默认 40）。
另输出 `n_compared` 与 `n_matching`（bcftools 的第 6、7 列）。

## 5. 输出

### 5.1 文件名

| kind | 后缀 | ID 文件 |
|---|---|---|
| AlleleCountDistance | `.dist` | `.dist.id` |
| Ibs | `.mibs` | `.mibs.id` |
| OneMinusIbs | `.mdist` | `.mdist.id` |
| Grm | `.rel` | `.rel.id` |
| Covariance | `.rel`（同 rel） | `.rel.id` |
| King | `.king` | `.king.id` |
| PiHat | `.genome` | —（列内自带 ID） |
| Phi | `.relatedness2` | —（列内自带 ID） |
| Ajk | `.relatedness` | — |
| Discordance | `.gtcheck` | — |

### 5.2 文本格式

- 矩阵类（dist/ibs/mdist/rel/king）：下三角 tab 分隔，`{:.6}` 或更高精度
  （plink 用 6 位；reed 用 9 位以便对拍）。
- 表格类（genome/relatedness2/relatedness/gtcheck）：逐行带 ID 与表头，
  列名逐字对齐参照工具。

### 5.3 二进制

`Bin` = 双精度方阵，小端，`n*n*8` 字节，无头（R 的 `readBin` 直接读）。

## 6. 精度契约

| 矩阵 | 契约 | 验证 |
|---|---|---|
| `grm` | 与 plink2 `--make-rel` 相对差 < 1e-9 | 实测对拍 |
| `allele_ct_distance` | 与 plink `--distance` 相对差 < 1e-9 | 实测对拍 |
| `king` | 与 plink2 `--make-king` 差 < 1e-9 | 实测对拍 |
| `phi` | 与 vcftools `--relatedness2` 差 < 1e-9 | 实测对拍 |
| `ajk` | 与 vcftools `--relatedness` 差 < 1e-9 | 实测对拍 |
| `discordance`(count) | 与 bcftools `gtcheck -E 0` 整数相等 | 实测对拍 |
| `pi_hat` | 与 plink `--genome` 差 < 1e-6 | 实测对拍（矩估计的系数实现细节多） |

## 7. 不做

见 01 §8 的"不做"表。核心原则：**只做有明确公式且能对拍的**；需要 HMM、
需要新数据形态、或不是个体×个体矩阵的，一律登记。
