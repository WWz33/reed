# doc_v16 / 02 规范：`reed-reduce` 的 API 与语义

> 本文件定义 `crates/reed-reduce` 的模块边界、类型、语义与精度契约。
> 与既有 crate 的关系：**完全独立**，不依赖 `reed-core`/`reed-io`/
> `reed-engine`，也不被它们依赖。它是原型库，不是命令面的一部分。
>
> 依赖策略：默认启用 `faer`（`faer-backend` feature）用于稠密线性代数；
> `--no-default-features` 必须仍能构建（纯 `std` 回退实现），这是可审计性
> 与受限环境的双重要求。两套后端的结果必须逐项一致。

## 1. 模块边界

```
reed-reduce
├── linalg      Mat、分块 matmul、gram、距离/向量工具
├── eig         sym_eig（tred2 + tql2）、SymEig
├── pca         standardize_binom2、Grm、Pca、pca_dense、pca_dense_pairwise
├── lanczos     LanczosConfig、lanczos_top_k、PartialEig
├── graph       Embedding、pcoa、kernel_pca、knn_graph、
│               laplacian_eigenmaps、diffusion_map、isomap、Kernel
├── grm_reduce  grm_distances、grm_similarity、grm_pcoa、grm_kernel_pca、
│               grm_laplacian、grm_diffusion、grm_isomap、grm_kernel、
│               grm_knn_weights
├── umap        UmapConfig、FuzzyGraph、fuzzy_simplicial_set、
│               find_ab_params、umap
└── extra       TsneConfig、tsne、nmf、lle、ppca
```

依赖方向：`linalg ← eig ← {pca, lanczos, graph} ← {grm_reduce, umap, extra}`。
**无环、无外部 crate。**

## 2. 数据类型

### 2.1 `Mat`

行主序稠密矩阵：`rows`、`cols`、`data: Vec<f64>`。

| 方法 | 语义 |
|---|---|
| `zeros(r, c)` / `from_fn(r, c, f)` | 构造；`from_fn` 接受 `FnMut` |
| `get(i, j)` / `set(i, j, v)` | 逐元素 |
| `row(i)` / `row_mut(i)` | 行切片 |
| `matmul(&other)` | 分块（32×32）矩阵乘，单线程 |
| `transpose()` | 转置（复制） |
| `matvec(&v)` / `tmatvec(&v)` | 矩阵-向量 |
| `gram()` | `self' * self`，只算下三角再镜像 |

**不变量**：`data.len() == rows * cols`；所有索引方法在越界时 panic
（原型阶段不做边界容错）。

### 2.2 `SymEig`

`values: Vec<f64>`（**升序**）、`vectors: Mat`（n×n，特征向量按列）。

`top_k(k) -> (Vec<f64>, Mat)` 返回**降序**的前 k 个特征值与对应列。

**精度契约**：对 n ≤ 200 的随机对称矩阵，`||Av - λv||/|λ| < 1e-10`
（单测 `random_symmetric_residual_is_small` 断言）。

### 2.3 `Standardized`

`w: Mat`（n×m，binom2 标准化后，缺失置 0）、`mean: Vec<f64>`
（每 SNP 的平均剂量 `2p`）、`sd: Vec<f64>`（每 SNP 的
`sqrt(2p(1-p))`）。

### 2.4 `Grm`

`a: Mat`（n×n 关系矩阵）、`n_pairs: Option<Mat>`（仅 pairwise 时非空）、
`m: usize`（使用的 SNP 数）。

| 构造 | 语义 | 缺失策略 |
|---|---|---|
| `Grm::mean_impute(&w)` | `W W' / m` | 均值填补（缺失在 W 中已是 0） |
| `Grm::pairwise(&x)` | 逐对除以共同非缺失数 | pairwise-complete（GCTA `.grm.N` 语义） |

**不变量**：`a` 对称；pairwise 时 `n_pairs(i,j)` = 两人共同非缺失的 SNP
数；无缺失时两种构造给出相同 `a`（差 < 1e-10，单测
`pairwise_and_mean_impute_agree_without_missing_data`）。

### 2.5 `Pca`

`values: Vec<f64>`（降序，即 λ/m）、`vectors: Mat`（n×k，单位特征向量）、
`pve: Vec<f64>`（每 PC 占总谱的比例）、`all_values: Vec<f64>`（升序全谱）。

| 方法 | 语义 |
|---|---|
| `scores() -> Mat` | `U * sqrt(λ)`，样本坐标（n×k） |
| `loadings(&w) -> Mat` | `W' U / sqrt(λ)`，SNP 载荷（m×k） |

**注意**：`scores` **不**再乘 W —— GRM 的特征向量已在样本空间。此前原型
误写成 `W · V`，导致坐标全错（本轮修复，见 STATUS）。

### 2.6 `PartialEig`

`values`（降序）、`vectors`（n×k）、`iterations`（块迭代次数）、
`max_residual`（返回对的最大相对残差）。

**精度契约（诚实版）**：特征值与稠密解一致到打印精度；`max_residual`
在 HapMap3 上约 1e-2（固定维数子空间、无 thick restart）。**生产使用应
优先稠密路径，或按 `docs/pca-algorithms.md` §14.11 实现带 thick restart
的块 Lanczos。**

### 2.7 `Embedding`

`coords: Mat`（n×k）、`values: Vec<f64>`（该方法的特征值）。

### 2.8 `FuzzyGraph`

`weights: Mat`（n×n 对称隶属度）、`knn_indices`、`knn_dists`。

**不变量**：`weights` 对称、值域 [0,1]（单测
`fuzzy_graph_is_symmetric_and_bounded`）。

## 3. 语义规范

### 3.1 binom2 标准化（唯一口径）

```
p_j  = mean(x_ij over called i) / 2
W_ij = (x_ij - 2 p_j) / sqrt(2 p_j (1 - p_j))
```

- 缺失 → 0（即中心化后的均值）。
- `sd_j <= 1e-9`（单态列）→ `W_ij = mean_j`（flashpca 矩阵路径语义）。
  注：GCTA/plink2 在此处置 0；reed-reduce 选 flashpca 语义以保持列居中，
  **这是本轮明确记录的差异**，若后续要对齐 GCTA 需加开关。
- 非 {0,1,2} 输入不检查（调用方责任）。

### 3.2 GRM 与 PCA

`A = W W' / m`；`pca_dense` 对 A 做全谱分解，取前 k。

- `pve_i = λ_i / Σλ`（分母是**全谱**，不是前 k 之和）。
- `scores = U sqrt(λ)`；`loadings = W' U / sqrt(λ)`。

### 3.3 图方法

| 方法 | 算子 | 取哪些特征对 |
|---|---|---|
| `pcoa` | `B = -0.5 J D² J` | 最大 k 个 |
| `kernel_pca` | 中心化后的核矩阵 K | 最大 k 个 |
| `laplacian_eigenmaps` | `L_sym = I - D^{-1/2} W D^{-1/2}` | **最小** k 个（跳过平凡零特征值） |
| `diffusion_map` | `S = D^{-1/2} W D^{-1/2}` | 最大 k 个（跳过 λ=1） |
| `isomap` | 测地距离 → PCoA | 最大 k 个 |

**退化谱的处理**：三个对称簇的 diffusion 谱会简并（λ₂=λ₃），简并子空间
内的任意旋转都合法 —— 测试用**二维联合分离度**而非单看 PC1（单测
`laplacian_and_diffusion_separate_clusters`）。

### 3.4 UMAP

严格按 umap-learn 0.5.12：

| 步骤 | 规范 |
|---|---|
| kNN | 精确欧氏（原型不做 NN-descent） |
| `smooth_knn_dist` | 二分搜索 sigma 使 `Σ exp(-(d-rho)/sigma) = log2(k)`；64 次上限；`SMOOTH_K_TOLERANCE=1e-5` |
| rho | 第 `local_connectivity` 个非零距离（默认 1） |
| 隶属度 | `d-rho <= 0 → 1.0`，否则 `exp(-(d-rho)/sigma)` |
| 模糊并集 | `A + A' - A*A'` |
| `find_ab_params` | 网格搜索 + Gauss-Newton 拟合 `1/(1+a x^{2b})` |
| 初始化 | 归一化拉普拉斯的最小非平凡特征向量 |
| SGD | 吸引 `-2ab d^{2(b-1)}/(a d^{2b}+1)`；负采样 `2γb/((0.001+d²)(a d^{2b}+1))`；clip ±4；alpha 衰减到 1/4 |

**ab 参数契约**：`find_ab_params(1.0, 0.1) = (1.576943618, 0.895060710)`，
与 umap-learn 的 `(1.576943460, 0.895060878)` 差 < 1e-7（单测断言）。

### 3.5 t-SNE / NMF / LLE / PPCA

- **t-SNE**：perplexity 二分搜索、对称化、Student-t 核、动量 0.5→0.8、
  每轮零中心化。
- **NMF**：Lee-Seung 乘性更新，`X ≈ W H`，非负。
- **LLE**：局部 Gram 正则化求解重建权重，取 `(I-W)'(I-W)` 的**最小**
  非平凡特征向量。
- **PPCA**：EM，E 步跳过缺失项 → **缺失值原生**。

## 4. 与既有 crate 的关系

| 问题 | 决定 | 理由 |
|---|---|---|
| 是否依赖 `reed-core` 的 `Mat`？ | **否** | 原型要能独立验证；接入时再适配 |
| 是否接入 CLI？ | **本轮不接** | 见 06 §"本轮不做" |
| 是否进 workspace？ | **是**（`members = ["crates/*"]` 自动包含） | 便于 `cargo test --workspace` 覆盖 |
| 是否受 `unsafe_code = "deny"` 约束？ | **是**（workspace lint 生效） | 全库零 unsafe |
| 是否要求 `missing_docs`？ | 是（workspace lint `warn`） | 当前零警告 |

## 5. 精度契约汇总

| 组件 | 契约 | 验证方式 |
|---|---|---|
| `sym_eig` | 残差 < 1e-10（n ≤ 200） | 单测 |
| `pca_dense` | 与 numpy `eigh` 相对误差 < 1e-12 | 手工对拍（实测 4.5e-16） |
| `pca_dense` | 与 plink2 exact 打印精度一致 | HapMap3 对拍 |
| `lanczos_top_k` | 特征值与稠密解一致；残差 ~1e-2（**已知限制**） | 单测 + 实测 |
| `find_ab_params` | 与 umap-learn 差 < 1e-7 | 单测 |
| `fuzzy_simplicial_set` | 对称、值域 [0,1]、基数 = log2(k) ± 0.05 | 单测 |
| 图方法 | 三簇 fixture 上分离比 > 1.5 | 单测 |

## 6. 不在本轮范围

- **CLI 命令**（`reed pca` / `reed reduce`）：需要先定命令面、输出格式、
  与 `popg`/`convert` 的边界，属独立一轮。
- **VCF 直读**：`reed-reduce` 只吃 `Mat`；VCF→`Mat` 的桥接在接入时做。
- **外部 BLAS（OpenBLAS/MKL）**：不做；性能通过 faer（纯 Rust、多线程）获得。
- **NN-descent / Barnes-Hut / Dijkstra**：大规模近似，属性能轮。
- **PanSN / hap 级**：与降维正交，属泛基因组轮。
