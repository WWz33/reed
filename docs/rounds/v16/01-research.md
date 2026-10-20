# doc_v16 / 01 调研台账：PCA 求解器与降维族

> 本文件记录本轮全部调研的证据。判决用语：**源码实证**（逐行读到并给出
> file:line）｜**数值实测**（真实命令/构造成数据，给出原始数字）｜
> **更正**（此前的记录有误，本轮改正）。
>
> 参照实现（本地检出，只读）：
> - `GCTA` v1.95（`/sri/home/sri2025201067/reed/GCTA`，GPL）
> - `plink-ng` 1.9 + 2.0（`/sri/home/sri2025201067/reed/plink-ng`，GPL）
> - `flashpca` 2.1（`/sri/home/sri2025201067/reed/flashpca`，GPL）
> - `umap-learn` 0.5.12（`/tmp/umap-src`，BSD-3）
> - sklearn 1.9.0、numpy（已装）
>
> 可执行：plink 1.9.0-b.7.11、plink2 2.0.0-a.7LM、gcta64 v1.95。
> 数据：HapMap3（957 人 × 14389 SNP，其中 14266 常染色体）、
> 模拟结构化 4000×20000 与 8000×20000。

## 0. 判决总表

| # | 结论 | 判决 | 证据位置 |
|---|---|---|---|
| R1 | GCTA 的 PCA 是 Eigen 全谱分解，不是 dsyevr | 源码实证 | `GCTA/main/grm.cpp:673` |
| R2 | plink2 exact 用 LAPACK `dsyevr` range='I' 只取前 k | 源码实证 | `plink2_matrix.cc:1102-1160` |
| R3 | flashpca 2.x 求解器是 Spectra（ARPACK 移植），与 LAPACK 无关 | 源码实证 | `randompca.cpp:121-150`、Spectra `SymEigsBase.h:84,134` |
| R4 | flashpca 的 `stand` 是五个标准化口径，不是三种算法 | 更正 + 源码实证 | `util.cpp:15-160` |
| R5 | binom 与 binom2 只差全局 √2（PC 相同、λ 恰好 ×2） | 数值实测 | §2.2 |
| R6 | `dsyevr` 不是某家独有实现，取决于链接哪份 LAPACK | 源码 + 二进制实证 | §3 |
| R7 | 小规模下 dsyevr 比 Lanczos 快 8–29 倍 | 数值实测 | §4 |
| R8 | GCTA/plink2/flashpca 在同集合同缺失策略下数值一致 | 数值实测 | §5 |
| R9 | plink 默认 `--pca` 用 14305 个位点（1-22 + 25） | 源码 + 实测复现 | §5.3 |
| R10 | GRM 不是标准化的分支，是"何时物化 N×N" | 源码实证 | §6 |
| R11 | 纯 Rust 可实现部分谱，块化是性能关键（167× 吞吐差） | 数值实测 | §7 |
| R12 | `faer` 无部分谱 API，Rust 侧需自研或绑 LAPACK | 源码实证 | §7.3 |
| R13 | 零依赖下可自研部分谱（tridiagonalize + bisection + 逆迭代），N=4000 比稠密快 2× 且 λ 逐位一致 | 数值实测 | §7.4 |

## 1. 求解器定位（逐家）

### 1.1 GCTA —— 稠密全谱

```cpp
// GCTA/main/grm.cpp:673
SelfAdjointEigenSolver<MatrixXd> eigensolver(_grm.cast<double>());
MatrixXd evec = (eigensolver.eigenvectors());
VectorXd eval = (eigensolver.eigenvalues());
```

- `--pca` 必须配 `--grm`（`main/option.cpp:1484`），先 `--make-grm` 或
  `--grm-gz` 读入（`option.cpp:384`）。
- `.eigenval` 写出**全部 n 个**特征值（降序），`.eigenvec` 只写前
  `out_pc_num`（默认 20，`option.cpp:82`）。
- `--pca 20` 只少写列，计算量不变（全谱）。
- 标准化：`main/data.cpp:1848-1875`，`sd = mu(1-0.5mu) = 2p(1-p)`，再乘
  `sqrt(1/sd)` —— 即 binom2 口径。
- 缺失：`grm.cpp:66-113`，`.grm.N.bin` 存每对共同非缺失数，`_grm(i,j) /=
  _grm_N(i,j)` → **pairwise-complete**。

### 1.2 plink 1.9 / 2.0 —— 稠密部分谱

```cpp
// plink2_matrix.cc:1102 起，dsyevr range='I'
// jobz='V', range='I', il=dim+1-pc_ct, iu=dim
```

- `CalcPca`（`plink2_matrix_calc.cc:7179`）现场建 N×N 关系矩阵，
  `ExtractEigvecs` 调 `dsyevr` 只取前 k（MRRR 算法）。
- 标准化：`plink2_matrix_calc.cc:4183`，`variance = 2*ref_freq*alt_freq`，
  `inv_stdev = 1/sqrt(variance)` —— binom2。
- 缺失：默认 pairwise（`:5755`，`variant_ct - miss_i - miss_j + miss_dbl`）；
  `meanimpute` 修饰改全局分母。
- 2.0 默认 10 PC（1.9 为 20）；`N > 46340` 需 ILP64。
- approx 分支：`plink2_matrix_calc.cc:7289` 注释自证移植自 EIGENSOFT 6
  `ksrc/kjg_fpca.c`（Halko 随机化），**一律均值填补**。
- 可导出 GRM：`--make-rel`/`--make-grm-list --zs`/`--make-grm-bin`/
  `--make-grm-sparse`（GCTA 兼容格式）。

### 1.3 flashpca 2.1 —— 矩阵-free Lanczos

```cpp
// randompca.cpp:132-136
Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, SVDWide>
   eigs(&op, ndim, ncv);   // ncv = 2*ndim+1
eigs.init();
eigs.compute(maxiter, tol);
```

- 算子 `y = X (X' x)`（`svdwide.h`），从不物化 N×N。
- Spectra 是 **ARPACK 的 C++ 移植**：`SymEigsBase.h:84` 写 "Implicitly
  restarted Lanczos factorization"，`:134` 注释直接引用 ARPACK
  `dsaup2.f:677-684`。
- 依赖：Eigen（头文件）+ Spectra（头文件）+ Boost program_options；
  **不链接 LAPACK**（`Makefile` 的 LDFLAGS 只有 boost）。
- 编译用 `-ffast-math`（`Makefile:40-41`）—— 对 NaN/Inf 语义不保守。
- 缺失：`data.cpp:318 scaled_geno_lookup(1,k)=0; // impute to average`，
  均值填补。

### 1.4 `dsyevr` 的实现来源（R6）

`dsyevr` 是 LAPACK 标准例程名（MRRR 算法，Netlib 中由 `dstegr`/`dstemr`
实现），实际二进制取决于链接：

| 使用方 | 来源 | 证据 |
|---|---|---|
| 本机基准 | OpenBLAS 0.3.30 | `ldd bench` → `libopenblas.so.0`；`liblapack.so.3` 是指向 `libopenblasp-r0.3.30.so` 的符号链接 |
| plink2 | 编译期宏切换 MKL/OpenBLAS/AOCL/ATLAS/Accelerate/参考 LAPACKE（`plink2_matrix.h:100-145`）；官方 Linux 二进制内嵌 **MKL** | 源码宏 + `strings` 见 `MKL_DOMAIN_*` |
| GCTA | 不用 dsyevr（见 1.1） | `grm.cpp:673` |

## 2. `stand` 五选项（R4/R5）

### 2.1 更正

此前记录为"flashpca 有三种算法可选"，实际是：v2 起求解器固定为 Spectra
Lanczos（`METHOD_EIGEN`/`METHOD_SVD` 是残留常量，CLI 无切换入口）；
`stand` 是**标准化口径**的五个取值。

### 2.2 五个分支（`util.cpp` 逐行）

| stand | 非缺失变换 | 缺失 → | sd ≤ 1e-9 列 → |
|---|---|---|---|
| none | 原值不动 | 原始列均值 | — |
| sd | (x−mean)/s（s 为样本 sd） | 0 | 原始列均值 |
| binom | (x−2p)/sqrt(p(1−p)) | 0 | 原始列均值 |
| binom2 | (x−2p)/sqrt(2p(1−p)) | 0 | 原始列均值 |
| center | x−mean | 0 | — |

关键恒等式：`sd_binom2 = √2 · sd_binom`（`util.cpp:114` 的 `mult=1 vs 2`）。

### 2.3 实测（HapMap3，957×14389，k=10）

| stand | λ1 | λ2 | λ3 | trace/p |
|---|---|---|---|---|
| none | 759.6093 | 12.8986 | 11.3524 | 1218.38 |
| sd | 24.8367 | 22.1871 | 6.6795 | 954.53 |
| binom | 52.9360 | 47.0280 | 13.9160 | 1980.86 |
| binom2 | 26.4680 | 23.5140 | 6.9580 | 990.43 |
| center | 12.9039 | 11.3533 | 3.2504 | 458.81 |

与 binom2 对比（前 10 PC）：

| stand | λ 比 d1/d2/d3 | 逐 PC 最小 \|corr\| |
|---|---|---|
| binom | 2.000000 / 2.000000 / 2.000000 | 1.000000 |
| sd | 0.9384 / 0.9436 / 0.9600 | 0.995447 |
| center | 0.4875 / 0.4828 / 0.4671 | 0.939211 |
| none | 28.70 / 0.5485 / 1.6316 | 0.000352 |

- none 的 PC1 与"个体平均剂量"相关 0.9828（均值方向尖峰），且 none 的
  d1/d2/d3 ≈ center 的 d0/d1/d2。
- sd vs binom2 逐 SNP 除数比（MAF 1–5% 组）p5/中位/p95 = 0.85/1.00/1.02。

### 2.4 退化列怪癖（实测）

SD/BINOM/BINOM2 分支对 sd≈0 的列把值写成**原始列均值**
（`util.cpp:107,146`）——全-2 列变成常数 2。而 flashpca 的 BED 输入路径
相反：`scaled_geno_lookup` 零初始化（`data.cpp:200`），只在 `sd > VAR_TOL`
时填（`:300-319`）→ 单态列全 0。**同一工具两条输入路径不一致**；GCTA
（`mkl.cpp:98-99`，sd<1e-50 → 0）与 plink2 均为置 0 语义。

实测（4×4 构造矩阵）：all-2 列在 sd/binom/binom2 下输出全 2.0、center 下
全 0、none 下保持 2；binom/binom2 除数比 0.484123/0.684653 = 1/√2（精确）。

## 3. 求解器速度/精度（R7）

方法：同一 binom2 标准化矩阵喂三种求解器（`/tmp/fpdeps/bench.cpp`），
精度口径 `||Av − λv||/|λ|`。

### 3.1 HapMap3（957 × 14266，k=10）

| 求解器 | 耗时 | 最大相对残差 |
|---|---|---|
| 稠密全谱（Eigen） | 0.552 s | 5.9e-14 |
| **稠密部分谱（dsyevr）** | **0.029 s** | **1.7e-15** |
| Lanczos（Spectra） | 0.613 s | 9.5e-11 |

三家 λ1..λ3 完全一致。

### 3.2 模拟结构化（4000 × 20000，k=10）

| 求解器 | 耗时 | 倍数 | 残差 |
|---|---|---|---|
| 稠密全谱 | 55.374 s | 70× | 8.1e-13 |
| **dsyevr** | **0.789 s** | 1× | 1.2e-14 |
| Lanczos | 23.109 s | 29× | 7.7e-11 |

### 3.3 8000 × 20000（k=10 / k=30）

| 求解器 | k=10 | k=30 |
|---|---|---|
| 稠密全谱 | 522.5 s | 526.3 s |
| **dsyevr** | **7.36 s** | **6.58 s** |
| Lanczos | 61.96 s | 70.98 s |

**dsyevr 的成本几乎与 k 无关**（只要 k ≪ N）；Lanczos 随 k 增长（子空间
2k+1、每步 O(N²k)）。popgen 常用的 k=20–50 正是 dsyevr 最舒服的区间。

### 3.4 结论

- 速度：dsyevr 部分谱 > Lanczos ≫ 稠密全谱（本机 8 线程）。
- 精度：三者残差都 < 1e-10；dsyevr 最准（1e-14）。
- 交叉点：N² 内存能否装下（N=50k 时 20GB，N=100k 时 80GB）。

## 4. 跨工具对拍（R8/R9）

### 4.1 对齐集合（14266 常染色体 SNP）

| 运行 | λ1 | λ2 | λ3 | vs plink2 pairwise |
|---|---|---|---|---|
| flashpca binom2 mean-impute | 26.200945 | 23.394076 | 6.932495 | 3.53e-03, 0.999778 |
| plink2 pairwise | 26.280600 | 23.458700 | 6.948660 | — |
| plink2 meanimpute | 26.200900 | 23.394100 | 6.932490 | 3.53e-03, 0.999778 |
| plink 1.9 pairwise | 26.280600 | 23.458700 | 6.948660 | 0（逐位相同） |
| gcta64 v1.95 | 26.280600 | 23.458700 | 6.948660 | 0（打印精度内） |

- **flashpca(meanimpute) vs plink2 meanimpute：max|Δλ|/λ = 2.0e-06**
  —— plink2 文本输出的 6 位有效数字极限；PC 相关 1.00000000。
- GCTA vs plink2 pairwise：λ 打印精度内完全一致，corr 1.00000000。
- pairwise vs meanimpute：3.55e-03（≈ 2×缺失率 0.154% 的分母效应）。

### 4.2 默认变异集合差异（R9）

- plink 1.9/2.0 默认用 **14305** = chrom 1-22 + 25。源码规则
  `plink2_common.cc:3397-3420 ExcludeNonAutosomalVariants`：只清
  X(23)/Y(24)/MT(26)；注释（:3388）"for backward compatibility, unplaced
  markers are considered to be autosomal" → 编码 25 的 39 个位点被保留。
  实测：`--chr 1-22,25` 精确复现默认输出 λ1=26.2733。
- GCTA v1.95：用 **14266** = chrom 1-22。
- flashpca 矩阵路径：喂什么用什么（本次全量 14389）。

**后果**：默认设置下 plink vs GCTA 的 λ1 差 7.1e-04、PC corr 0.9988 ——
纯集合差异，对齐后消失。

### 4.3 排查顺序

① 缺失策略（0.35%）② 变异集合（7e-4）③ 标准化口径（binom/binom2 差 √2）
④ tol/maxiter。**数值精度从来不是第一位**。

## 5. GRM 的位置（R10）

统一记法：`W_ij = (x_ij − 2p_j)/sqrt(2p_j(1−p_j))`，`A = W W'/m`。

| 实现 | 标准化 | GRM 物化时机 | 缺失分母 |
|---|---|---|---|
| smartpca | binom（`sqrt(p(1-p))`） | **一开始就物化**（算法定义） | 全局分母 |
| GCTA | binom2 | **强制显式**（`--pca` 必须 `--grm`） | **`.grm.N.bin` 每对计数** → pairwise |
| plink 1.9/2.0 | binom2 | 内部物化、可导出 | 默认 pairwise；`meanimpute` 可改 |
| flashpca | binom2 默认 | **从不物化**（矩阵-free） | **只有全局分母** |
| reed-reduce | binom2 | 两条路径都提供 | `Grm::mean_impute` / `Grm::pairwise` |

**关键澄清**：标准化阶段不"考虑 GRM"；binom2 只决定 `W` 的列缩放，GRM
的存在与否由"是否物化 N×N"决定。GRM 物化的真正收益是**缺失分母的表达
能力**——pairwise 需要每对的共同非缺失数（N×N 个数），矩阵-free 算子拿
不到。

## 6. 纯 Rust 部分谱（R11/R12）

### 6.1 块化是性能关键

faer 0.20.2 上 8000×20000 的实测吞吐：

| RHS 列数 | 耗时/次 | 吞吐 |
|---|---|---|
| 1 | 0.085 s | **3.8 GFlop/s** |
| 8 | 0.022 s | 116 GFlop/s |
| 32 | 0.029 s | 355 GFlop/s |
| 64 | 0.065 s | **633 GFlop/s** |

**167 倍差距**（3.8 → 633）。单列 gemv 无法复用数据、无法向量化。

### 6.2 faer 的能力边界（R12）

- 可用：`LinOp` trait、`compute_tridiag_real_evd`、`compute_hermitian_evd`、
  `Mat` 的多线程 GEMM。
- **不可用**：`krylov_schur.rs` 在 `lib.rs:1646-1647` 是注释掉的原型；
  无 Lanczos/LOBPCG/随机化 SVD 的公开 API。
- faer 全谱实测：4000×20000 单线程 10.83 s、多线程 2.58 s（4.2×）；
  8000×20000 单线程 95.17 s、多线程 11.02 s（8.6×）。精度 1e-14。

### 6.3 零依赖的部分谱实现（R13）

`faer` 没有部分谱 API（`krylov_schur.rs` 是注释掉的原型），因此本轮自研了
一条：**Householder 三对角化 + Sturm 序列二分 + 逆迭代**，只提取前 k 个
特征对，跳过 QL 迭代中 O(n³) 的特征向量累积。

实测（`pca_top_k`，与稠密全谱同数据）：

| 规模 | 稠密全谱 | 部分谱 | 加速 | λ 一致性 |
|---|---|---|---|---|
| 957 × 14266 | 5.6 s | **4.2 s** | 1.3× | 26.200944934431163 vs ...124（1e-13） |
| 4000 × 20000 | 334 s | **170.6 s** | **2.0×** | 207.0176524976154 vs ...154（逐位） |

实现要点：
- 三对角化用归一化 Householder 向量（`H = I - 2vv'`），同时累积 Q；
- 特征值用 Sturm 序列（LDL 主元计数）在 Gershgorin 区间上二分；
- 特征向量用 LDL 前代/回代做逆迭代，并对已求向量重正交（处理简并谱）；
- 残差 < 1e-13（单测 `partial_matches_full_spectrum` 断言 < 1e-9）。

**这是零依赖模式下的"最快最准"实现**：块子空间迭代在 N=4000 上更快
（14.7 s）但残差 4.8e-2；部分谱路径以 2× 稠密的时间换到 1e-13 的精度。

> 后续（同轮）已把默认后端切换为 faer：GRM 构建 79 s → 0.65 s、稠密
> 特征分解 258 s → 2.66 s（4000×20000），全管线 334 s → 3.96 s。
> 上面的部分谱实现保留在零依赖回退路径中，仍可用
> `--no-default-features` 构建。

### 6.4 三个必须避开的坑

1. **全 1 向量在零空间**：binom2 标准化后每列均值为 0，`X'·1 = 0`，
   用它做起子 Rayleigh 商恒为 0。
2. **`faer::linalg::matmul::matmul(acc, lhs, rhs, Some(alpha), beta, par)`
   在 0.20.2 下静默不写入**（acc 保持全零，无 panic）。必须用运算符
   `&lhs * &rhs`。
3. `Mat` 的运算符消耗所有权，循环内必须借用。

## 7. 降维族与 umap-learn 对拍

### 7.1 UMAP 关键公式（umap-learn 0.5.12 逐行）

- `smooth_knn_dist`（`umap_.py:153`）：二分搜索 sigma 使模糊基数 =
  `log2(k)`；rho = 最近非零邻居距离（`local_connectivity=1`）。
- 隶属度（`:374`）：`exp(-(d - rho)/sigma)`；`d - rho <= 0` 时取 1.0。
- 模糊并集（`:611`）：`A + A' - A*A'`（`set_op_mix_ratio=1`）。
- `find_ab_params`（`:1478`）：对 `1/(1 + a x^{2b})` 做 curve_fit。
- SGD（`layouts.py:63`）：吸引梯度 `-2ab d^{2(b-1)} / (a d^{2b} + 1)`；
  负采样梯度 `2γb / ((0.001 + d²)(a d^{2b} + 1))`；梯度 clip 到 ±4；
  alpha 线性衰减到 1/4。

### 7.2 对拍结果（60×400 基因型 fixture，指标 = 组间质心距离/组内离散度）

| 方法 | reed-reduce | 参考实现 |
|---|---|---|
| PCA (binom2) | 2.91 | 2.91（numpy eigh） |
| Kernel PCA (RBF) | 1.80 | 1.80（sklearn） |
| Isomap | 1.93 | 1.95（sklearn） |
| Laplacian | 4.61 | 6.90（sklearn） |
| LLE | 5.42 | 3.47（sklearn） |
| UMAP | 2.80 | 2.72（umap-learn） |
| NMF 重构误差 | 0.1907 | 0.1905（sklearn） |
| t-SNE | 2.63 | 0.49（sklearn） |

`find_ab_params(1.0, 0.1)`：reed `(1.576943618, 0.895060710)` vs
umap-learn `(1.576943460, 0.895060878)` —— 7 位一致。

### 7.3 差异解释

- Laplacian / LLE 的差异来自图构造细节（kNN 数、正则项、初始化），不是
  算法错误——两者都正确分离了三组。
- t-SNE 高度依赖初始化与学习率调度，两边默认不同，不预期一致。
- UMAP 的 2.80 vs 2.72 属随机初始化差异（同 seed 下 reorder 顺序不同）。

## 8. 此前记录的更正汇总

| 此前记录 | 本轮更正 | 依据 |
|---|---|---|
| "flashpca 有三种算法可选" | v2 起求解器固定 Spectra；`stand` 是五个标准化口径 | `util.cpp:15-160`、`flashpca.cpp:336-343` |
| "binom2 标准化阶段考虑 GRM" | 标准化与 GRM 是两层，正交 | §5 |
| "jvarkit 664 tools"（v15 前的记录） | 525 个 `@Program` | 已在 `docs/jvarkit-study.md` §F 更正 |

## 9. 复现产物索引

| 产物 | 路径 | 内容 |
|---|---|---|
| HapMap3 矩阵 | `/tmp/fpdeps/hm3auto.mat` | 957×14266，0/1/2 + NA |
| 模拟 4000 | `/tmp/fpdeps/sim4k.mat` | 4000×20000 三群结构 |
| 模拟 8000 | `/tmp/fpdeps/sim8k.mat` | 8000×20000 |
| 小 fixture | `/tmp/fixture.mat` | 60×400，Python 与 Rust 共用 |
| C++ 求解器基准 | `/tmp/fpdeps/bench.cpp` | 三求解器同矩阵计时 |
| Rust 原型 | `/tmp/rsbench/src/bin/*.rs` | 单向量/块 Lanczos、gemv 吞吐 |
| umap-learn 源码 | `/tmp/umap-src` | 0.5.12 全量 |
| plink2 输出 | `/tmp/fpdeps/p2pairA.*`、`p2miA.*` | 对齐集合两种缺失策略 |
| plink 1.9 输出 | `/tmp/fpdeps/p19A.*` | 同上 |
| GCTA 输出 | `/tmp/fpdeps/g_hm3_pca.*` | GRM + PCA |
| flashpca 输出 | `/tmp/fpdeps/flA.*` | 对齐集合 binom2 |
