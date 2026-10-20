# reed doc_v16：PCA 与降维族（2026-10-20）

> 本目录是 `docs/rounds/v15/` 之后的**第十六轮规范**，性质是**新能力轮**
> （与 v6 的"命令面扩展"同型，与 v7/v9–v12/v15 的缺陷修复轮不同）。
> 本轮不重开既有结构项，只新增一个**独立 crate** `reed-reduce` 与其配套
> 调研文档。
>
> 来源：一次围绕 PCA 求解器与降维族的**源码级调研 + 数值对拍**，跨
> GCTA / plink 1.9 / plink2 / flashpca / EIGENSOFT / sklearn / umap-learn
> 七套实现。调研产物见 [01-research.md](01-research.md)，规范见
> [02-spec.md](02-spec.md)，路线见 [06-roadmap.md](06-roadmap.md)，
> 实施状态见 [STATUS.md](STATUS.md)。
>
> **先写 doc 再落代码**（v13/v14/v15 同款工作方式）。本轮代码已随文档
> 一起落地：`crates/reed-reduce`（降维库）与 `reed pca`（CLI 命令，
> `docs/spec/cli.md` §11b）。`reed reduce`（暴露其余降维族）未做，
> 理由见 06 §P2「明确不加」。

---

## 0. 为什么这是一轮

`reed` 目前的命令面覆盖 VCF 下游的统计、过滤、转换，但**没有任何降维
能力**。教学站覆盖度文档（`docs/speciation-course-coverage.md`）把 PCA
列为 **G1 缺口**（课程用 plink LD-prune + PCA，reed 无法替代）；泛基因组
方向文档把群体结构可视化列为后续项。本轮把这一块从"调研"推进到
"可验证的原型"，并回答三个此前没有定论的问题：

1. **哪个求解器最快最准**（此前只有"plink2 用 dsyevr、GCTA 用 Eigen 全谱、
   flashpca 用 Spectra"的定性描述，没有同数据同规模的实测）。
2. **`stand` 五选项的语义与相互关系**（此前误以为是三种算法，实际是五个
   标准化口径）。
3. **GRM 到底在哪一步进入计算**（此前有"binom2 标准化阶段考虑 GRM"的
   误解）。

## 1. 本轮的产出（两层）

| 层 | 交付物 | 性质 |
|---|---|---|
| 调研 | `docs/pca-algorithms.md` §1–§15 | 源码级对比 + 数值对拍，全部带 file:line 与实测数字 |
| 代码 | `crates/reed-reduce`（默认 faer 后端；`--no-default-features` 可退回零依赖） | PCA（binom2/binom）+ 10 种降维方法 + GRM 驱动入口 |
| 命令 | `reed pca`（`docs/spec/cli.md` §11b） | plink 兼容的 `.eigenvec`/`.eigenval`/`.pve` |

## 2. 核心结论（详见 01）

### 2.1 求解器

| 规模（N × 20000 SNP, k=10） | dsyevr（部分谱） | 块 Lanczos | 稠密全谱 |
|---|---|---|---|
| 957 | **0.029 s** | 0.61 s | 0.55 s |
| 4000 | **0.79 s** | 23.1 s | 55.4 s |
| 8000 | **7.36 s** | 62.0 s | 522 s |

**在本机直到 N=8000（N²=512MB）dsyevr 都是最快且最准**（残差 1e-14）。
Lanczos 的翻盘条件是 N² 装不下（N ≳ 50k）或必须流式。

### 2.2 五家实现的一致性

同一份 HapMap3（14266 常染色体 SNP，957 人）上：

| 对拍 | 特征值最大相对差 | 前 10 PC 最小 \|corr\| |
|---|---|---|
| flashpca(binom2+填补) vs plink2 meanimpute | 2.0e-06（= 文本位数极限） | 1.00000000 |
| GCTA vs plink2 pairwise | 0（打印精度内） | 1.00000000 |
| plink 1.9 vs plink2 pairwise | 0（逐位相同） | 1.00000000 |
| pairwise vs meanimpute | 3.5e-03 | 0.99978 |

**数值精度没有赢家**；"结果对不上"的来源按影响排序是：① 缺失策略
（0.35%）② 变异集合（plink 默认 14305 = 1-22+25，GCTA 14266 = 1-22）
③ 标准化口径（binom vs binom2 差 √2）④ 迭代法的 tol。

### 2.3 `stand` 五选项

`binom2`（默认）与 `binom` 只差全局 √2 —— **PC 完全相同、特征值恰好 ×2、
PVE 相同**。`sd` 在 HWE 附近 ≈ binom2（前 10 PC 相关 0.995+）。`center`
是协方差 PCA，`none` 的 PC1 被均值方向占据（与平均剂量相关 0.98）。

### 2.4 GRM 的位置

**GRM 不是标准化的分支**，是"何时物化 N×N"。标准化（binom2）决定每列
怎么缩放，GRM 决定要不要把 `W W'` 算出来。两者正交。GRM 物化的真正收益
是**缺失分母的表达能力**：pairwise-complete 需要每对个体的共同非缺失
SNP 数（GCTA 的 `.grm.N.bin`），矩阵-free 算子拿不到，只能退化成全局
分母。

## 3. 代码：`crates/reed-reduce`

| 模块 | 内容 |
|---|---|
| `linalg` | 行主序 `Mat`、分块 matmul、gram、距离 |
| `eig` | 稠密对称特征分解（EISPACK `tred2` + `tql2`）+ **部分谱求解**（tridiagonalize + bisection + inverse iteration） |
| `pca` | binom2 标准化、GRM（mean-impute / pairwise）、稠密 PCA、**部分谱 PCA**、scores/loadings |
| `lanczos` | 矩阵-free 块子空间迭代（top-k） |
| `graph` | PCoA、Kernel PCA、Laplacian、Diffusion、Isomap |
| `grm_reduce` | GRM 驱动的上述五种 |
| `umap` | fuzzy simplicial set + SGD（umap-learn 0.5.12 逐行） |
| `extra` | t-SNE、NMF、LLE、PPCA-EM |

**双后端**：默认 `faer`（多线程 GEMM + 稠密特征分解），
`--no-default-features` 退回纯 `std` 实现。两套后端结果逐项相同。
38 个单元测试全绿、零编译警告。

### 3.1 精度（三级对拍）

| 级别 | 参照 | 结果 |
|---|---|---|
| 与 LAPACK | numpy `eigh`（60×400 fixture） | 相对误差 **4.5e-16** |
| 与 plink2 | HapMap3 真实数据 | λ1/λ2/λ3 = 26.200944934 / 23.394076211 / 6.932494547，**与 plink2 逐位一致** |
| 与 sklearn/umap-learn | 同 fixture | PCA 2.91 vs 2.91、KernelPCA 1.80 vs 1.80、UMAP 2.80 vs 2.72、NMF 0.1907 vs 0.1905 |

`find_ab_params(1.0, 0.1)` = `(1.576943618, 0.895060710)`，umap-learn 为
`(1.576943460, 0.895060878)` —— 7 位有效数字一致。

### 3.2 部分谱路径（最快且最准）

`pca_top_k` 用 tridiagonalize + Sturm 二分 + 逆迭代，只提取前 k 个特征对：

| 规模 | 稠密全谱 | 部分谱 | 加速 |
|---|---|---|---|
| 957 × 14266 | 5.6 s | **4.2 s** | 1.3× |
| 4000 × 20000 | 334 s | **170.6 s** | **2.0×** |

λ 与稠密解**逐位一致**（sim4k：207.0176524976154 vs 207.0176524976154），
残差 < 1e-13。这是零依赖模式下的"最快最准"实现；默认的 faer 后端更快
（HapMap3 全管线 0.40 s，4000×20000 3.96 s）。

### 3.3 已知限制

1. 块子空间迭代在 HapMap3 上 `max_rel_residual ≈ 1e-2`：特征值已收敛
   （与稠密解一致到打印精度），但单个 Ritz 向量残差在该量级——原型用
   固定维数子空间、无 thick restart。稠密路径是精确参照。
2. 稠密路径 O(N³) 单线程（sim4k 334 s），因为无 BLAS。生产应走
   `faer`（见 `docs/pca-algorithms.md` §14.10）或块 Lanczos + BLAS。
3. UMAP 用精确 kNN（非 NN-descent）、t-SNE 无 Barnes-Hut、Isomap 用
   Floyd-Warshall —— 小数据够用，大规模需替换。

## 4. 文件

| 文件 | 内容 |
|---|---|
| [01-research.md](01-research.md) | 调研台账：七套实现的源码定位、五轮数值实验的完整数据、三处此前误解的更正 |
| [02-spec.md](02-spec.md) | `reed-reduce` 的 API 规范：模块边界、类型、语义、精度契约、与既有 crate 的关系 |
| [06-roadmap.md](06-roadmap.md) | 路线 P1–P4：原型收口 → CLI 接入 → 性能路径 → 生态扩展；每批的门；"本轮不做" |
| [STATUS.md](STATUS.md) | 实施状态与证据索引 |

## 5. 与既有文档的关系

- `docs/pca-algorithms.md` 是本轮的**调研正文**（§1–§15），v16 不重复其
  内容，只引用结论。
- `docs/speciation-course-coverage.md` 的 G1（PCA 缺口）在本轮从"缺失"
  推进到"已可用"（`reed pca`），覆盖度表在下一轮统一改写。
- `docs/pangenome-directions.md` 的群体结构项同理。
