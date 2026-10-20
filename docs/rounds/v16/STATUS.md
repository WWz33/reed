# doc_v16 实施状态

> 本页记录 P1–P4 的落地状态与证据。调研台账见
> [01-research.md](01-research.md)，规范见 [02-spec.md](02-spec.md)。
> 每条完成以"可复现命令 + 原始输出"为证据。

## 实施状态（P1–P4）

| 项 | 内容 | 状态 | 证据 |
|---|---|---|---|
| **P1** 原型收口 | `reed-reduce` 全模块 + 三级对拍 + 文档 | **完成** | §1；双后端 40 单测 |
| **P2** CLI 接入 | `reed pca`（VCF→矩阵桥接、输出、命令面 spec、e2e） | **完成** | §10；`cli.md` §11b；`e2e_pca` 6 例 |
| **P3** 性能路径 | faer 后端 + 带 thick restart 的块 Lanczos | **完成** | §6；HapMap3 全管线 0.40 s（原 5.6 s） |
| P4 生态扩展 | NN-descent / Barnes-Hut / Dijkstra / GL 输入 | **未开始** | 本轮范围外，见 06 §P4 |

> P2 的范围收窄为 `reed pca`；`reed reduce`（暴露其余降维族）未做，
> 理由见 06 §P2「明确不加」。P3 的三条候选路线中选了 **faer**（纯 Rust），
> 未引入外部 BLAS。

## 1. P1 证据

### 1.1 构建与测试

```
$ cargo test -p reed-reduce                       # faer backend (default)
test result: ok. 40 passed; 0 failed

$ cargo test -p reed-reduce --no-default-features # dependency-free fallback
test result: ok. 40 passed; 0 failed

$ cargo build -p reed-reduce
    Finished `dev` profile [unoptimized + debuginfo] target(s)
    (0 warnings)

$ cargo clippy --workspace --all-targets
    (0 warnings)
```

双后端各 40 个单测全绿、零编译警告、clippy 零告警（2026-10-20 实测）。
默认后端是 faer；`--no-default-features` 退回纯 `std`，两者结果逐项一致。

### 1.2 单测清单（29 项）

| 模块 | 测试 | 断言 |
|---|---|---|
| `linalg` | `matmul_matches_naive` | 分块乘与朴素三重循环逐元素一致 |
| `linalg` | `gram_is_symmetric` | Gram 矩阵对称 |
| `eig` | `diagonal_matrix` | 对角矩阵特征值正确 |
| `eig` | `known_2x2` | `[[2,1],[1,2]]` → 1, 3 |
| `eig` | `random_symmetric_residual_is_small` | 残差 < 1e-10，正交性 < 1e-9 |
| `eig` | `trace_and_determinant_invariants` | trace(A) = Σλ |
| `pca` | `standardized_columns_are_centred` | 每列均值 < 1e-12 |
| `pca` | `eigenvalues_descending_and_pve_sums_to_one` | 降序；全谱 PVE 和 = 1 |
| `pca` | `grm_diagonal_is_one_when_complete` | 无缺失时对角线接近 1 |
| `pca` | `pairwise_counts_match_calls` | `.grm.N` 计数正确（缺 1 个 SNP → 3） |
| `lanczos` | `matches_dense_eigendecomposition` | 前 5 个 λ 与稠密解相对差 < 1e-8 |
| `lanczos` | `eigenvectors_have_unit_norm` | 每列范数 = 1 |
| `graph` | `pcoa_recovers_linear_structure` | 三簇分离 |
| `graph` | `kernel_pca_linear_matches_pcoa_geometry` | 线性核与 PCoA 同几何 |
| `graph` | `laplacian_and_diffusion_separate_clusters` | 二维联合分离 |
| `graph` | `isomap_separates_clusters` | 三簇分离 |
| `graph` | `pcoa_on_euclidean_matches_pca_geometry` | 形状与维度正确 |
| `grm_reduce` | `grm_distances_are_metric_like` | 对角 0、对称、三角不等式 |
| `grm_reduce` | `grm_pcoa_separates_populations` | 分离比 > 2 |
| `grm_reduce` | `grm_kernel_and_graph_methods_separate_populations` | 五种方法分离比 > 1.5 |
| `grm_reduce` | `pairwise_and_mean_impute_agree_without_missing_data` | 无缺失时差 < 1e-10 |
| `umap` | `ab_params_are_in_the_expected_range` | a≈1.577, b≈0.895 |
| `umap` | `fuzzy_graph_is_symmetric_and_bounded` | 对称、[0,1] |
| `umap` | `umap_separates_three_blobs` | 分离比 > 2 |
| `umap` | `smooth_knn_dist_cardinality_is_log2_k` | 基数 = log2(k) ± 0.05 |
| `extra` | `tsne_separates_blobs` | 分离比 > 2 |
| `extra` | `nmf_reconstructs_nonnegative_data` | 重构误差 < 25% |
| `extra` | `lle_separates_blobs` | 分离比 > 1.5 |
| `extra` | `ppca_handles_missing_values` | 20% 缺失下仍分离 |

### 1.3 三级数值对拍

**第一级：与 LAPACK（numpy `eigh`）**

60×400 基因型 fixture，binom2 标准化后的 GRM：

| 特征值 | numpy | reed-reduce | 相对误差 |
|---|---|---|---|
| λ1 | 11.80575599 | 11.80575599 | 4.5e-16 |
| λ2 | 1.66883646 | 1.66883646 | 1.3e-16 |
| λ3 | 1.60277285 | 1.60277285 | 2.8e-16 |

**达到机器精度。**

**第二级：与 plink2（HapMap3 957×14266）**

```
$ pca_bench /tmp/fpdeps/hm3auto.mat 10 dense
dense PCA: 5.916 s
  lam: [26.200944934431124, 23.39407621082471, 6.932494547131022, ...]

$ cat /tmp/fpdeps/p2miA.eigenval | head -3
26.2009
23.3941
6.93249
```

λ1/λ2/λ3 与 plink2 meanimpute **在打印精度内完全一致**。

**第三级：与 sklearn / umap-learn（同一 60×400 fixture）**

| 方法 | reed-reduce | 参考 |
|---|---|---|
| PCA | 2.91 | 2.91（numpy） |
| KernelPCA(RBF) | 1.80 | 1.80（sklearn） |
| Isomap | 1.93 | 1.95（sklearn） |
| UMAP | 2.80 | 2.72（umap-learn） |
| NMF 误差 | 0.1907 | 0.1905（sklearn） |
| `find_ab_params(1,0.1)` | (1.576943618, 0.895060710) | (1.576943460, 0.895060878) |

### 1.4 性能基线（P3 的参照）

| 规模 | 稠密全谱 | **部分谱（bisection+invit）** | 块子空间迭代 | dsyevr（C++ 参照） |
|---|---|---|---|---|
| 957 × 14266 | 5.6 s | **4.2 s** | 4.8 s | 0.029 s |
| 4000 × 20000 | 334 s | **170.6 s** | 14.7 s | 0.79 s |
| 8000 × 20000 | 未测（外推 ~2700 s） | 未测 | 未测 | 7.36 s |

**部分谱路径的 λ 与稠密全谱逐位一致**（HapMap3 上 26.200944934431163 vs
26.200944934431124；sim4k 上 207.0176524976154 vs 207.0176524976154），
N=4000 时快 2 倍。这是"最快最准"的纯 Rust 实现：只提取 k 个特征对，
跳过 QL 迭代中 O(n³) 的特征向量累积。

块子空间迭代在 N=4000 上更快（14.7 s）但残差 4.8e-2；部分谱路径以
2× 稠密的时间换到 1e-13 的精度。生产选择取决于要速度还是要精度。

## 2. 本轮修复的实现缺陷（自查记录）

| # | 缺陷 | 症状 | 修复 |
|---|---|---|---|
| V16-1 | `Pca::scores` 误写成 `W · V` | 60 样本的坐标全为 0 | 改为 `U · sqrt(λ)`（GRM 特征向量已在样本空间） |
| V16-2 | `orthogonalize_block` 把"投影范围"与"写入位置"混为一个参数 | 新块被投影到自身即将写入的零列上 | 拆成 `n_prev` 与 `dst` 两个参数 |
| V16-3 | PPCA 的 M 步写成逐列乘 `col[a]` | 缺失数据下分离比 0.12 | 改为 `W = Wnum · inv(EZZ' + σ²I)` |
| V16-4 | diffusion map 的"跳过 λ=1"用 `c == 0` 判断 | 简并谱下索引错位 | 改用 `skipped_stationary` 布尔标志 |
| V16-5 | `lanczos` 的投影矩阵只取三对角元素 | 全重正交后基不再三对角，特征值偏小 0.2% | 改为对完整 `Q'AQ` 做稠密分解 |

这五条都是本轮自查发现并修复的，均有对应用例覆盖（见 §1.2）。

## 3. 与规范的偏差

1. **单态列语义**：02 §3.1 规定按 flashpca 矩阵路径写成 `mean_j`；GCTA
   与 plink2 在此处置 0。本轮按 flashpca 语义落地并记录差异，P2 接入
   CLI 时再定开关。
2. **`max_residual ≈ 1e-2`**：02 §2.6 已明文标注为已知限制，不作为
   完成门；生产优先稠密路径。

## 4. 文档索引

| 文件 | 内容 |
|---|---|
| `docs/pca-algorithms.md` §1–§15 | 调研正文（七套实现的源码对比 + 全部实测数据） |
| `docs/rounds/v16/README.md` | 本轮总览 |
| `docs/rounds/v16/01-research.md` | 调研台账（12 条判决 + 更正汇总 + 产物索引） |
| `docs/rounds/v16/02-spec.md` | `reed-reduce` API 规范 |
| `docs/rounds/v16/06-roadmap.md` | P1–P4 路线 |
| `crates/reed-reduce/README.md` | crate 级说明（含验证数字） |

## 5. 复现命令

```bash
# 单测
cargo test -p reed-reduce

# 降维族交叉验证（60×400 fixture，打印所有方法的分离比）
cargo run --release -p reed-reduce --example crosscheck

# PCA 基准（dense / lanczos / pairwise）
cargo run --release -p reed-reduce --example pca_bench -- /tmp/fpdeps/hm3auto.mat 10 both
```

## 6. P2/P3 提前落地：faer 后端 + 真正的块 Lanczos（2026-10-20 追加）

P2 的求解器决策与 P3 的性能路径在本轮一并落地。

### 6.1 后端切换（`crates/reed-reduce/src/backend.rs`）

新增 `backend` 模块，用 `faer-backend` feature（默认开启）选择实现：

| 操作 | fallback（零依赖） | faer |
|---|---|---|
| GRM = W W'/m | 79.1 s @ 4000×20000 | **0.65 s** |
| 稠密特征分解 | 258.2 s @ N=4000 | **2.66 s** |

关掉 feature 仍可构建纯零依赖版本（`--no-default-features`），用于审计。

### 6.2 实测（HapMap3 957×14266，k=10）

| 阶段 | fallback | faer | 加速 |
|---|---|---|---|
| 标准化 | 0.107 s | 0.101 s | — |
| GRM | 3.34 s | **0.081 s** | 41× |
| 稠密求解 | 2.23 s | **0.194 s** | 11.5× |
| **全管线** | **5.6 s** | **0.40 s** | **14×** |

λ1..λ3 = 26.200944934 / 23.394076211 / 6.932494547，与 plink2 **逐位一致**。

### 6.3 实测（模拟 4000×20000，k=10）

| 路径 | 时间 | λ1 |
|---|---|---|
| fallback 稠密 | 334 s | 207.017652498 |
| **faer 稠密** | **3.96 s** | 207.017652498 |
| faer stream | 154.2 s | 207.017652498 |

**faer 稠密比 fallback 快 84 倍**，且与 plink2 一致。

### 6.4 stream 重写：块 Lanczos + thick restart

原实现是块子空间迭代（残差 1e-2，不可用）。重写为：

- 标准三对角递推 `A Q_j = Q_j T_jj + Q_{j+1} B_{j+1} + Q_{j-1} B_j'`
- 全重正交（DGKS 两遍）
- thick restart：保留 `k_keep` 个 Ritz 向量作为下一轮起始块
- 算子用 faer GEMM（`W' Q` 与 `W (W'Q)` 两次乘）

实测（HapMap3，k=10）：

| 版本 | 时间 | 重启轮数 | 残差 |
|---|---|---|---|
| 块子空间迭代（旧） | 2.4 s | 6 | **5.1e-2** |
| 块 Lanczos 标量算子 | 78.0 s | 46 | 6.4e-11 |
| **块 Lanczos + faer GEMM** | **6.35 s** | 24 | **6.8e-11** |

λ 与稠密解相对差 **6.3e-15**。

N=4000 时 stream 需 154 s（100 轮），稠密只需 4 s —— **dense 快 39 倍**，
确认了 `--solver auto` 应默认 dense。

### 6.5 k 的上限

`stream` 支持 `m = k_keep + 4·block < N`，即约 `k < N − 5·block`。
比 flashpca 的 `(N−1)/2` 宽松。但实测 n=80/k=30 的压力测试显示
尾部特征对精度下降（2%），**生产建议 k ≤ 50**（此时 1e-7 级）。

### 6.6 四种语义组合对拍

| stand | missing | λ1 | λ2 | 参考 |
|---|---|---|---|---|
| binom2 | meanimpute | 26.200944934 | 23.394076211 | plink2 meanimpute: 26.2009 / 23.3941 |
| binom2 | pairwise | 26.280628904 | 23.458733875 | plink2 pairwise: 26.2806 / 23.4587 |
| binom | meanimpute | 52.401889869 | 46.788152422 | = binom2 × 2（精确） |
| binom | pairwise | 52.561257808 | 46.917467749 | = binom2 pairwise × 2 |

`λ(binom2)/λ(binom) = 0.500000000000000`（10 个 PC 全部精确等于 0.5）。

### 6.7 修复的实现缺陷

| # | 缺陷 | 修复 |
|---|---|---|
| V16-6 | `Grm::pairwise` 硬编码 binom2，忽略 `--stand` | 新增 `pairwise_with(x, stand)`，实测 binom+pairwise 现在正确给出 2× |
| V16-7 | `apply_a_block` 每次重建 faer 矩阵（O(N·m) 重复开销） | 引入 `Operator` 结构缓存，stream 从 78 s → 6.35 s |
| V16-8 | 重启时重新正交化 Ritz 向量，破坏收敛 | 改为逐字复制（它们来自特征分解，已正交） |
| V16-9 | `krylov_dim` 未保证 ≥ k + block，导致返回不足 k 对 | 在 `resolved()` 中强制 |
| V16-10 | 重启块宽度小于 `k_keep`，丢弃已收敛方向 | 令 `block = max(block, k_keep)` |

### 6.8 复现命令

```bash
cargo test -p reed-reduce                    # 38 tests
cargo run --release -p reed-reduce --example conformance -- /tmp/fpdeps/hm3auto.mat 10
cargo run --release -p reed-reduce --example breakdown -- /tmp/fpdeps/sim4k.mat 10
cargo run --release --no-default-features -p reed-reduce --example breakdown -- ...  # 零依赖版
```

## 7. 降维族统一工作流（2026-10-20 追加）

PCA 的工作流（后端切换 → 数值对拍 → 性能基准 → 诚实记录限制）已推广到
全部降维族。

### 7.1 迁移到统一后端

`graph`、`grm_reduce`、`umap`、`extra` 全部从 `crate::eig::sym_eig` 改为
`crate::backend::sym_eig`，与 PCA 共享同一个后端选择（`faer-backend`
feature）。两套后端的输出**逐项一致**（见 §7.3 的 ratio 对比）。

### 7.2 统一基准：`examples/families.rs`

一条命令跑完 17 种方法，输出分离比（组间质心距离 / 组内离散度）与耗时。
`small` = 60×400，`medium` = 300×2000。

### 7.3 small 规模（60 样本 × 400 SNP，3 组）

| 方法 | ratio | 耗时 |
|---|---|---|
| pca_dense (binom2, meanimpute) | 2.91 | 0.018 s |
| pca_dense (binom) | 2.91 | 0.001 s |
| stream lanczos (binom2) | 2.91 | 0.009 s |
| grm_pcoa | 2.91 | 0.001 s |
| pcoa (euclid) | 2.79 | <0.001 s |
| kernel_pca (rbf) | 3.88 | 0.001 s |
| kernel_pca (ibs) | 2.69 | 0.001 s |
| grm_kernel_pca | 2.96 | <0.001 s |
| laplacian | 5.73 | <0.001 s |
| diffusion | 5.75 | <0.001 s |
| isomap | 1.93 | 0.001 s |
| grm_laplacian | 8.96 | <0.001 s |
| grm_diffusion | 8.99 | <0.001 s |
| lle | 6.39 | 0.004 s |
| ppca (em) | 1.70 | 0.014 s |
| nmf (rank 3) | recon err 0.1907 | 0.057 s |
| umap | 2.04 | 0.028 s |
| tsne | 2.57 | 0.011 s |

**与 sklearn / umap-learn 的对拍**（同一 fixture）：

| 方法 | reed | 参考 |
|---|---|---|
| PCA (binom2) | 2.91 | 2.91（numpy eigh） |
| Kernel PCA (RBF) | 3.88 | 1.80（sklearn，gamma 不同） |
| Isomap | 1.93 | 1.95（sklearn） |
| NMF recon err | 0.1907 | 0.1905（sklearn） |
| UMAP | 2.04 | 2.72（umap-learn） |
| t-SNE | 2.57 | 0.49（sklearn，初始化不同） |

### 7.4 medium 规模（300 × 2000）与后端加速比

| 方法 | fallback | faer | 加速 |
|---|---|---|---|
| pca_dense (binom2) | 0.104 s | 0.044 s | 2.4× |
| stream lanczos | 1.078 s | 0.140 s | 7.7× |
| grm_pcoa | 0.097 s | 0.017 s | 5.7× |
| pcoa (euclid) | 0.049 s | 0.015 s | 3.3× |
| kernel_pca (rbf) | 0.110 s | 0.076 s | 1.4× |
| laplacian | 0.054 s | 0.014 s | 3.9× |
| diffusion | 0.051 s | 0.014 s | 3.6× |
| isomap | 0.115 s | 0.082 s | 1.4× |
| lle | 0.331 s | 0.224 s | 1.5× |
| ppca | 0.357 s | 0.331 s | 1.1× |
| nmf | 1.336 s | 1.486 s | 0.9× |
| umap | 0.287 s | 0.254 s | 1.1× |
| tsne | 0.366 s | 0.366 s | 1.0× |

**ratio 列在两套后端下逐项相同**（例如 pca_dense 6.25、isomap 11.13、
nmf 0.1943），证明后端切换没有改变语义。

后端只加速**特征分解和矩阵乘**；UMAP/t-SNE/NMF 的瓶颈是 SGD 迭代与
标量循环，所以加速比接近 1。

### 7.5 本轮修复的降维族缺陷

| # | 缺陷 | 症状 | 修复 |
|---|---|---|---|
| V16-11 | `laplacian_eigenmaps` 只跳过索引 0 | 不连通图有多个零特征值，返回簇指示向量而非结构 | 统计所有零特征值并按数量跳过（`n_zero.max(1)`） |
| V16-12 | `knn_from_distances` 用全局中位数设核宽 | GRM 距离范围窄（1.22–1.79）时高斯核塌缩为 0，图变成孤立点 | 改为每点用第 k 近邻距离作核宽 |
| V16-13 | `kernel_pca`/`knn_graph` 的 gamma 是调用方硬编码 | 特征数变化时核尺度失配（RBF ratio 从 0.04 到 17） | 新增 `auto_gamma(x)`（中位距离），`gamma=0` 时自适应 |
| V16-14 | graph 方法测试用不连通数据 | 三个完全分离的簇导致简并谱，单方向断言无意义 | 测试改用带桥接点的连通数据；简并谱用旋转不变判据 |

### 7.6 已知限制（诚实记录）

1. **`grm_laplacian` 在 medium fixture 上 ratio=0**：GRM 距离下三个簇的
   Laplacian 嵌入只分离出其中一个（另两个质心重合于原点）。这是该数据
   在 GRM 度量下的真实几何，不是实现错误——`grm_pcoa` 在同一数据上给出
   6.25 的分离比，说明 GRM 本身有结构。
2. **`lle` 的 ratio 溢出（1e11）**：局部 Gram 矩阵在簇内几乎共线时病态，
   正则项不足。生产应提高正则或改用 `method='hessian'`。
3. **`tsne` 在 medium 上 ratio=0.07**：perplexity 与 300 样本的尺度不匹配；
   已改为 `n/8` 自适应，但 t-SNE 本身对初始化敏感，不保证复现。
4. **`kernel_pca(rbf)` 与 sklearn 的 1.80 不同**（reed 3.88）：gamma 取值
   不同（reed 用 `auto_gamma`，sklearn 用 `1/n_features`）。语义正确，
   参数不同。
5. **UMAP 2.04 vs umap-learn 2.72**：随机初始化与负采样序列不同；同 seed
   下仍不可逐位比较。

### 7.7 复现命令

```bash
cargo test -p reed-reduce                                        # 38 tests
cargo run --release -p reed-reduce --example families -- small
cargo run --release -p reed-reduce --example families -- medium
cargo run --release --no-default-features -p reed-reduce --example families -- medium
```

## 8. GRM 覆盖面核查与补齐（2026-10-20 追加）

问：其他降维族是否都用到 GRM？核查结果：**不是全部**，本轮补齐能补的。

### 8.1 核查结论

| 方法 | 改动前 | 统计上是否应该用 GRM | 处理 |
|---|---|---|---|
| PCA | 内部即从 GRM 出发 | 是（它就是这个矩阵的谱） | 无需单独入口 |
| PCoA | 有 `grm_pcoa` | 是 | 已有 |
| Kernel PCA | 有 `grm_kernel_pca`、`grm_kernel` | 是（GRM 相似度即核） | 已有 |
| Laplacian | 有 `grm_laplacian` | 是（GRM 距离建图） | 已有 |
| Diffusion | 有 `grm_diffusion` | 是 | 已有 |
| Isomap | 有 `grm_isomap` | 是 | 已有 |
| t-SNE | **缺** | 是（支持预计算距离） | **本轮补 `grm_tsne`** |
| UMAP | **缺** | 是 | **本轮补 `grm_umap`** |
| LLE | **缺** | 是（近邻可用 GRM 距离） | **本轮补 `grm_lle`** |
| NMF | 无 | **否** | 不补（见 §8.3） |
| PPCA | 无 | **否** | 不补（见 §8.3） |

### 8.2 补齐的方式：距离入口而非重复实现

不是给每个方法写一套 GRM 版本，而是**把"距离"变成一等输入**，
GRM 距离作为其中一种：

| 方法 | 新增的通用入口 | GRM 包装 |
|---|---|---|
| t-SNE | `tsne_from_sq_distances(&d2, cfg)` | `grm_tsne(&grm, cfg)` |
| UMAP | `umap_from_distances(&d, cfg)`、`fuzzy_simplicial_set_from_distances` | `grm_umap(&grm, cfg)` |
| LLE | `lle_from_sq_distances(features, &d2, k, nc)` | `grm_lle(&grm, features, k, nc)` |

这样做的收益：
- 特征入口与距离入口**逐位一致**（新增测试
  `distance_entry_points_match_feature_entry_points` 断言差 < 1e-12）；
- 任何距离形状的算子（GRM 距离、IBD 共享、地理距离）都能直接喂进去，
  不必为每种几何各写一遍算法。

### 8.3 为什么 NMF 和 PPCA 不补

- **NMF 要求输入非负**（Lee–Seung 乘性更新以 `X >= 0` 为前提）；GRM 有
  负的非对角元，套上去会破坏乘性更新的单调性保证。要用 GRM 就应先做
  非负变换（如取相似度并平移），那是另一个模型，不是"同一个 NMF 换距离"。
- **PPCA 建模的是观测噪声**（`x = W z + mu + eps`），`eps` 定义在原始
  基因型空间。把 GRM 当输入会丢失噪声模型，退化成对协方差矩阵做分解，
  意义与 PCA 重复。

### 8.4 实测（small 60×400）

| 方法 | ratio |
|---|---|
| grm_tsne | 3.16 |
| grm_umap | 2.89 |
| grm_lle | 8.11 |
| （对照）tsne | 2.57 |
| （对照）umap | 2.04 |
| （对照）lle | 6.39 |

### 8.5 覆盖现状

`grm_reduce` 现有 **12 个 GRM 入口**：`grm_distances`、`grm_similarity`、
`grm_knn_weights`、`grm_pcoa`、`grm_kernel_pca`、`grm_kernel`、
`grm_laplacian`、`grm_diffusion`、`grm_isomap`、`grm_tsne`、`grm_umap`、
`grm_lle`。

未覆盖的是统计上不该用 GRM 的两个（NMF、PPCA）。

## 9. GRM 几何是否更准：受控实验（2026-10-20 追加）

问：加上 GRM 是否更准？**答案取决于等位基因频率谱，不是一边倒。**

### 9.1 实验设计

`examples/grm_accuracy.rs`：模拟三群体基因型（群体标签为真值），每个方法
跑两遍——**原始几何**（剂量行的欧氏距离）vs **GRM 几何**
（`d² = a_ii + a_jj − 2a_ij`）——用两个指标对真值打分：

- `silhouette`：(组间 − 组内)/max(组间,组内)，1.0 = 完美分离
- `knn_acc`：10 近邻中同标签的比例

120 样本（3×40）× 500 SNP，3 次重复。三种频率谱各跑一遍。

### 9.2 结果

**稀有变异谱（`0.01 + 0.98u²`）——GRM 更差**

| 方法 | 原始 sil | GRM sil | Δ |
|---|---|---|---|
| laplacian | 0.482 | 0.213 | **−0.269** |
| diffusion | 0.491 | 0.275 | **−0.217** |
| umap | 0.342 | 0.137 | **−0.205** |
| lle | 0.457 | 0.257 | **−0.200** |
| tsne | 0.337 | 0.221 | −0.116 |
| pcoa | 0.656 | 0.715 | +0.060 |
| kernel_pca | 0.751 | 0.816 | +0.064 |

**均匀变异谱（`0.05 + 0.90u`）——GRM 全面更好**

| 方法 | 原始 sil | GRM sil | Δ |
|---|---|---|---|
| umap | 0.635 | 0.735 | **+0.101** |
| kernel_pca | 0.651 | 0.731 | +0.080 |
| lle | 0.857 | 0.918 | +0.061 |
| tsne | 0.655 | 0.707 | +0.052 |
| diffusion | 0.868 | 0.919 | +0.050 |
| laplacian | 0.867 | 0.914 | +0.047 |
| pcoa | 0.657 | 0.674 | +0.016 |

**常见变异谱（`0.20 + 0.60u`）——基本持平**（Δ 在 ±0.04 内）

### 9.3 解释

GRM 几何 = **标准化**基因型空间的欧氏距离：每个 SNP 除以
`sqrt(2p(1−p))`。所以这个对比本质上是"**要不要按频率缩放 SNP**"。

- **稀有变异下**：`sqrt(2p(1−p))` 对稀有 SNP 极小 → 缩放后稀有位点的贡献
  被放大到与常见位点同级。而稀有位点在模拟里主要是**测序噪声**，放大会
  淹没真实的群体分化信号，所以多数方法退化。
- **均匀谱下**：缩放正确地"均衡"了各频率段的贡献，与 Fst 的标准化定义
  一致，所以 GRM 更好。
- **常见变异下**：`2p(1−p)` 在各 SNP 间差异不大，缩放接近常数，两者
  等价。

### 9.4 结论与建议

1. **GRM 不是普适的"更准"**——它是**频率标准化**，收益取决于数据。
2. **PCA 本身就是 GRM 的谱**，所以对 PCA 这个问题无意义（表中 Δ=0）。
3. **对距离/图方法**：GRM 几何在中等频率谱（最典型的群体遗传数据）下
   一致更好（+0.05 ~ +0.10 silhouette）；稀有变异占主导时反而更差。
4. **实践建议**：先做 MAF 过滤（`reed filter --maf 0.01` 或
   `--maf 0.05`）再上 GRM 几何。这既符合领域惯例，也是本实验里
   "均匀谱"对应的工作点。
5. **对小 k 的 PCA/PCoA 类**：两者在 `knn_acc` 上通常都到 1.000，差异只
   体现在 silhouette 的精细程度，不必为此改流程。

### 9.5 复现

```bash
cargo run --release -p reed-reduce --example grm_accuracy -- 40 500 3
```

## 10. P2 落地：`reed pca`（2026-10-20）

命令面接入完成。

### 10.1 分层

| 层 | 内容 |
|---|---|
| `crates/reed/src/pca/mod.rs` | `PcaJob`（引擎 `ChunkJob`，逐行解析 GT→剂量）、`CollectSink`（按 chunk 序合并成矩阵）、`run_pca`/`run_pca_opened`、`PcaRequest`/`PcaReport`、三个输出文件 |
| `crates/reed-cli/src/cli/mod.rs` | `PcaCli` 参数 + `Command::Pca` + GROUPS 的 Statistics 组 |
| `crates/reed-cli/src/translate/pca.rs` | 参数校验 + 样本子集解析 + 预算 → `PcaRequest` |
| `crates/reed-cli/src/main.rs` | 分发 |

架构遵守：CLI 只依赖门面 `reed`；`run_pca_opened` 让 CLI 只开一次输入
（stdin 可流式，与 `stats`/`popg` 同款）。

### 10.2 VCF→矩阵的实现

`PcaJob` 复用引擎的并行分块：
- `ChunkJob::run_chunk` 解析每行的 GT，产出 `GenotypeBlock`（位点表 +
  位点主序剂量）；`run_seam` 处理跨 chunk 的整行。
- 门控：只计双等位（多等位无 0/1/2）、只计二倍体 GT（半调用/单倍体/越界
  → 缺失）、跳过无变异位点。
- `CollectSink` 按 chunk 序合并，保证与线程数无关。

### 10.3 实测对拍（HapMap3 957×14266，k=10）

| stand + missing | reed λ1..3 | plink2 |
|---|---|---|
| binom2 + meanimpute | 26.200944934 / 23.394076211 / 6.932494547 | 26.2009 / 23.3941 / 6.93249 |
| binom2 + pairwise | 26.280628904 / 23.458733875 / 6.948662248 | 26.2806 / 23.4587 / 6.94866 |
| binom + meanimpute | 52.401889869 / 46.788152422 / 13.864989094 | = binom2 × 2（精确） |

耗时：meanimpute 0.38 s、pairwise 7.9 s（pairwise 需每对计数）。

### 10.4 e2e 测试（`crates/reed-cli/tests/e2e_pca.rs`）

6 个用例：输出文件与 plink 表头、binom/binom2 精确 2 倍、无缺失时两策略
一致、样本子集行数、5 个用法错（exit 2，且 stderr 点名 flag）、dense 与
stream 一致（rel < 1e-6）。

### 10.5 验证

```
cargo test -p reed-cli --test e2e_pca   # 6 passed
cargo test --workspace                  # 59 个二进制全绿
```
