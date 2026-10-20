# doc_v16 / 06 路线：P1–P4

> 顺序原则：**先收口原型**（P1，本轮已完成），再定命令面（P2，需独立
> 设计），再做性能（P3，有明确的实测基线），最后扩展生态（P4）。
> 每批完成即更新 STATUS.md 并跑本批门。

| 批次 | 内容 | 门 | 状态 |
|---|---|---|---|
| **P1** 原型收口 | `reed-reduce` 全模块落地；三级对拍；README/文档 | 40 单测绿（双后端各 40）、零警告；HapMap3 与 plink2 逐位一致 | **完成** |
| **P2** CLI 接入 | `reed pca`：VCF→`Mat` 桥接、输出格式、与 `popg`/`convert` 的边界 | 命令面 spec（cli.md §11b）+ e2e；输出与 plink2 对拍 | **完成** |
| **P3** 性能路径 | `faer` feature（稠密多线程）+ 带 thick restart 的块 Lanczos | N=4000 稠密 3.96 s（原 334 s）；Lanczos 残差 6.8e-11 | **完成** |
| **P4** 生态扩展 | 大规模近似（NN-descent/Barnes-Hut/Dijkstra）；PCAngsd 式 GL 输入；PanSN/hap 级 | 各自独立门 | 未开始 |

## P1：原型收口（本轮，已完成）

### 交付

| 项 | 证据 |
|---|---|
| 10 个源文件、4599 行 | `wc -l crates/reed-reduce/src/*.rs` |
| 双后端 | 默认 faer；`--no-default-features` 零依赖，两者 38 测试均全绿 |
| 40 个单测（双后端） | `cargo test -p reed-reduce [--no-default-features]` |
| 零编译警告 | `cargo build -p reed-reduce` |
| 与 numpy 逐位一致 | 相对误差 4.5e-16 |
| 与 plink2 逐位一致 | λ1/2/3 = 26.200944934 / 23.394076211 / 6.932494547 |
| 与 umap-learn 一致 | `find_ab_params` 7 位有效数字 |

### 门（已过）

```
cargo test -p reed-reduce          # 40 passed
cargo build -p reed-reduce         # 0 warnings
cargo tree -p reed-reduce          # 1 line
```

## P2：CLI 接入（下一轮）

### 必须先定的问题

1. **命令面**：新增 `reed pca` 还是并入 `reed convert --pca`？
   - 倾向新增 `reed pca`：PCA 有自己的参数空间（`--pcs`、`--stand`、
     `--missing`、`--solver`、`--ld-prune`），塞进 convert 会让其选项爆炸。
2. **VCF→`Mat` 桥接**：
   - 复用 `reed-io` 的解码 + `reed-filter` 的 LD 修剪？
   - 还是独立读一遍（简单，但重复 IO）？
   - 倾向复用：`reed pca` 内部调 `run_filter` 的 `--indep-pairwise` 产
     位点列表，再扫一遍取值。
3. **输出格式**：
   - `.eigenvec` / `.eigenval` 对齐 plink（课程 R 脚本直接消费）；
   - 附加 `.pve`（flashpca 风格）；
   - 可选 `--out-grm` 写 GCTA 兼容格式。
4. **缺失策略默认**：
   - `pairwise`（对齐 plink/GCTA 默认，课程友好）还是 `meanimpute`
     （对齐 flashpca）？
   - 倾向 `pairwise` + `--missing meanimpute` 开关。
5. **求解器选择**：
   - 默认稠密部分谱（需 `dsyevr`，即引入 LAPACK 依赖或自研）；
   - 还是先用 `reed-reduce` 的稠密全谱（faer 后端）？
   - 倾向：P2 先接 `reed-reduce` 的稠密路径（零新依赖），P3 再换。

### 门

- cli.md 新增 §"reed pca"（选项表、输出文件、与 plink 的差异）；
- e2e：HapMap3 上 `reed pca` 输出与 plink2 `--pca` 对拍（容差 1e-6）；
- `--out-grm` 与 `--make-grm-bin` 逐字节或数值一致。

## P3：性能路径

### 实测基线（当前）

| 规模 | 稠密（faer） | 目标 |
|---|---|---|
| 957 × 14266 | 5.9 s | ≤ 1 s（dsyevr 0.029 s） |
| 4000 × 20000 | 334 s | ≤ 15 s（dsyevr 0.79 s + 标准化） |
| 8000 × 20000 | ~2700 s（外推） | ≤ 60 s |

### 三条候选路线

| 路线 | 依赖 | 预期 | 代价 |
|---|---|---|---|
| **faer**（纯 Rust） | `faer` crate | 8000 规模 11 s（实测） | 引入一个依赖，但仍是纯 Rust |
| **ndarray-linalg + OpenBLAS** | `ndarray-linalg`、`openblas-src` | 与 C++ dsyevr 同速（7.4 s） | 系统库依赖 |
| **自研块 Lanczos + thick restart** | 无 | 按 633 GFlop/s 外推 ~1–2 s | 实现量大，需仔细处理收敛 |

### 门

- N=8000 稠密 ≤ 15 s（faer 路线）；
- 或 Lanczos 残差 ≤ 1e-8（自研路线）；
- 精度不回归（与 plink2 对拍仍 1e-6 内）。

## P4：生态扩展

| 项 | 说明 | 依赖 |
|---|---|---|
| NN-descent | UMAP 的 kNN 从 O(N²m) 降到近似 | 独立实现或绑 `pynndescent` 的算法 |
| Barnes-Hut t-SNE | O(N log N) | 四叉树 |
| Dijkstra/双向 | Isomap 的 O(N³) Floyd-Warshall | 优先队列 |
| PCAngsd 式 GL 输入 | 用 VCF 的 GL/PL 代替硬 call | 需 `reed-io` 暴露 GL |
| PanSN / hap 级 | 单倍型级降维 | 需泛基因组轮先落地 |

## 本轮不做

1. **CLI 接入**——需要先定命令面（P2 的五问），属独立一轮；本轮只交付
   可验证的原型库。
2. **引入外部 BLAS**——性能优化用 faer（纯 Rust），不引入 OpenBLAS/MKL。
3. **修改 `docs/speciation-course-coverage.md` 的 G1 状态**——PCA 尚未
   接入 CLI，覆盖度表暂不改写；P2 完成后再改。
4. **修改 `docs/pangenome-directions.md`**——同上。
5. **外部 BLAS**——不做（faer 已提供多线程）。
6. **大规模近似算法**——属 P4。
7. **EIGENSOFT smartpca 的对拍**——本地无源码（GPL，且未检出），
   其算法与 GCTA 同构（显式 GRM + 稠密全谱），已有 GCTA 对拍覆盖。

## 风险

| 风险 | 缓解 |
|---|---|
| `reed-reduce` 的 `Mat` 与 `reed-core` 的类型不一致，接入时需转换层 | 原型只吃 `Mat`，接入时写一个 `From<&reed_core::X>` 适配；或直接复用其矩阵类型（P2 决定） |
| 块子空间迭代残差 1e-2，若被误用于生产会给出不精确的 PC 坐标 | README 与 02 §2.6 明文标注"生产优先稠密路径"；P3 修复 |
| 零依赖下 O(N³) 稠密路径无法处理 N > 10k | P2 的 CLI 在 N 超阈值时明确报错并指向 P3 路线 |
| 单态列的 `mean_j` vs `0` 语义与 GCTA/plink2 不同 | 02 §3.1 已记录差异；P2 定 `--monomorphic` 开关 |

## P2 命令面决策（2026-10-20 定稿）

P2 的开放问题按三条判据收敛：**对齐领域标准**、**只做 VCF 下游**、
**不重复三款软件已有的**。定稿如下。

### 保留的维度

| 维度 | 保留值 | 默认 | 理由 |
|---|---|---|---|
| 标准化 | `binom2`、`binom` | **binom2** | binom2 是 GCTA/plink2/flashpca 共同口径；binom 与它只差全局 √2（PC 相同、λ ×2），一行代码支持老 EIGENSOFT 对拍 |
| 缺失 | `pairwise`、`meanimpute` | **pairwise** | pairwise 是 plink/GCTA 默认、课程管线友好；meanimpute 对齐 flashpca。两条均已实测逐位复现 plink2 对应模式 |
| 求解器 | `dense`、`stream` | **auto** | auto：`N²·8 < 内存预算` → dense，否则 stream |
| 模式 | `global`、`local` | global | local = Local PCA，三款软件都没有，纯 VCF 下游 |
| k 的选择 | `--tw` | 关闭 | Tracy–Widom 显著性检验，三款软件都缺 |

### 明确不加的参数

| 参数 | 决定 | 理由 |
|---|---|---|
| `--missing em` | **不加** | 换的是估计量（PPCA 概率模型），不是同一 PCA 的缺失策略；且没有"特征值"这一量，`.eigenval` 语义需重新定义。若要做应是独立的 `--method ppca` |
| `--use-gl` | **不加** | 需要 `reed-io` 暴露 GL/PL 读侧；PCAngsd 是独立算法（似然模型），不是"PCA 换个输入"，工程量远超一个 flag |
| `--indep-pairwise` | **不加** | 已在共享 `SiteFilterArgs` 里，但只有 `filter` 实现；`stats`/`popg` 均显式拒绝（`translate/common.rs:253`）。LD 修剪是上游预处理，与统计量无关；`reed pca` 接受它会在同一参数名下引入第二条实现路径（违背 v15 "多份实现各自漂移"的教训）。帮助里写明"先跑 `reed filter --indep-pairwise`" |
| `sd`/`center`/`none` 标准化 | 不加 | 面向非基因型矩阵；reed 只吃 0/1/2，`sd` 在 HWE 附近 ≈ binom2，无增量 |
| Robust/sparse/Bayesian PCA、flash | 不加 | 各带新超参或重实现；先看 Local PCA 的实际需求 |
| admixture 族、树/图 | 不加 | 不同统计问题，另立项 |
| UMAP/t-SNE 进 `reed pca` | 不加 | 可视化层，不是 PCA；原型已有，若要暴露走独立命令 |
| `reed reduce`（暴露其余降维族） | **本轮不做，登记为后续轮** | 见下方 §P2.1 |

### P2.1 为什么本轮不做 `reed reduce`

`reed-reduce` 里除 PCA 外的 10 种方法（PCoA、kernel PCA、Laplacian、
diffusion、Isomap、LLE、NMF、PPCA、UMAP、t-SNE）目前只有库入口，没有
CLI。本轮不暴露它们，理由：

1. **输出契约不同**。PCA 有领域标准输出（plink 的 `.eigenvec/.eigenval`），
   可以直接对拍；其余方法没有公认格式，每个都要单独设计并解释语义
   （尺度、符号、随机性），是一轮独立工作。
2. **随机方法的复现性问题**。UMAP/t-SNE 依赖初始化与 SGD 轨迹，跨后端
   不逐位一致（v16 STATUS §7.6 实测）；发布为命令前需要先定 seed/收敛
   契约。
3. **参数尺度敏感**。v16 STATUS §7.5 记录过 RBF 核宽度与 kNN 图尺度失配
   导致结果退化；CLI 必须先把自适应策略（`auto_gamma`）和默认值钉死。
4. **与 `popg`/`convert` 的边界未定**。例如 `popg --ld-decay` 已经做
   pairwise r²，Isomap/LLE 的邻域图与它是什么关系，需要先写清。

结论：**本轮只交付 `reed pca`**（有标准输出、可对拍、参数语义明确）；
其余方法留在库中，作为后续轮的输入。若后续要做，建议按
`reed reduce --method pcoa|tsne|umap|...` 一次设计统一输出契约。

### 命令面

```
reed pca in.vcf.gz --pcs 20 -o out
  --stand binom2|binom            默认 binom2
  --missing pairwise|meanimpute   默认 pairwise
  --solver auto|dense|stream      默认 auto
  --window 100k                   启用 Local PCA（每窗口一份输出）
  --tw                            输出 Tracy–Widom 检验
  --out-grm                       写 GCTA 兼容 .grm.bin/.grm.id/.grm.N.bin
  --threads N / --memory MB       复用现有预算机制
```

输出：
- `out.eigenvec` / `out.eigenval`（plink 兼容）
- `out.pve`
- `--window` 时 `out.window.<contig>_<start>_<end>.eigenvec`
- `--tw` 时 `out.tw`（每 PC 的 Tracy–Widom 检验 p 值）

### 求解器语义（重要）

| 路径 | 数据流 | 算法 | k 上限 | 对标 |
|---|---|---|---|---|
| `dense` | 物化 N×N | faer `compute_hermitian_evd`（全谱） | 无 | **数据流形态同 GCTA**；性能争取接近 plink2（faer 8.8s vs dsyevr 7.36s @N=8000） |
| `stream` | 矩阵-free | 块 Lanczos + thick restart | `m = b·⌈(2k+1)/b⌉ < N`，b=8 时约 N−8 | **算法族同 flashpca**（IRAM），但用块化提升 BLAS 效率 |

**注意**：`dense` 只是数据流形态像 GCTA，**算法上不采用 GCTA 的 Eigen 全谱路线**
（那是实测最慢的：N=8000 需 522 s）。我们用 faer 的多线程全谱（8.8 s）。

`stream` 当前原型是块子空间迭代（残差 1e-2，**不可用**），P3 必须重写为带
thick restart 的块 Lanczos。k 上限比 flashpca 的 `(N−1)/2` 更宽松
（因为用固定子空间维数 m 而非 ncv=2k+1），但实际 k 很少超过 50，
真正的价值是块化的速度。

### 对拍矩阵（P2 的门）

4 种语义组合 × 3 列参考：

| 我们的组合 | 参考 | 容差 |
|---|---|---|
| binom2 + pairwise | plink2 pairwise | 1e-6 |
| binom2 + meanimpute | plink2 meanimpute | 1e-6 |
| binom + pairwise | plink2 pairwise 的 λ × 0.5（或 flashpca --standx binom） | 1e-6 |
| binom + meanimpute | flashpca --standx binom | 1e-6 |
