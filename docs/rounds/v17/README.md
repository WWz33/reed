# reed doc_v17：个体间距离/关系矩阵族（2026-10-20）

> 本目录是 `docs/rounds/v16/` 之后的**第十七轮规范**，性质是**新能力轮**。
> 目标：把三家参照工具（plink 1.9/2.0、vcftools、bcftools）的**个体×个体
> 矩阵**能力收进 reed，并把 v16 里散在 `reed-reduce::pca` 的 GRM 构建
> 迁到独立模块，成为所有矩阵的统一底座。
>
> 先写 doc 再落代码（v13–v16 同款）。

---

## 0. 为什么这是一轮

v16 交付了 `reed pca`，但 GRM 构建写在 `reed-reduce::pca` 里，只服务于
PCA。实际上**同一份标准化基因型能派生出多种个体间矩阵**，三家参照工具
各有一套：

| 工具 | 矩阵 | 选项 |
|---|---|---|
| plink 1.9 | 基因组距离（allele-count 单位）、IBS 矩阵、1−IBS | `--distance [ibs] [1-ibs] [allele-ct] [flat-missing]` |
| plink 1.9 | IBD（PI_HAT） | `--genome [full] [unbounded] [nudge]` |
| plink 1.9/2.0 | GRM（variance-standardized 关系矩阵）、协方差矩阵 | `--make-rel [cov] [meanimpute]`、`--make-grm-gz/bin` |
| plink 2.0 | KING-robust kinship | `--make-king`、`--make-king-table` |
| vcftools | Manichaikul 2010 relatedness（φ） | `--relatedness2`（`.relatedness2`：`INDV1 INDV2 N_AaAa N_AAaa N1_Aa N2_Aa RELATEDNESS_PHI`） |
| bcftools | 基因型不一致（discordance） | `gtcheck`（`.gtcheck`：`DCv2` 行 + 每对 `n_mismatch`/`n_checked`） |

这些矩阵是同一批统计量的不同侧面，应当共享标准化与缺失处理，而不是各写
一遍（v15 的教训："多份实现各自漂移"）。

## 1. 本轮产出

| 层 | 交付物 |
|---|---|
| 库 | 新 crate `reed-distance`：标准化 + 全部矩阵族 + 统一输出格式 |
| 迁移 | `reed-reduce::pca` 的 `Grm`/`standardize` 移入 `reed-distance`，PCA 改为依赖它 |
| CLI | `reed distance --type ibs\|1-ibs\|grm\|king\|pihat\|kinship\|discordance` |
| 对拍 | plink `--distance`/`--make-rel`/`--make-king`、vcftools `--relatedness2`、bcftools `gtcheck` |

## 2. 关键决定（详见 02-spec）

1. **一个标准化底座，多个矩阵**。`standardize(stand)` 产出的 `W` 是所有
   矩阵的共同输入；缺失策略（meanimpute / pairwise）在矩阵层决定，不在
   每个矩阵里各写一遍。
2. **输出格式统一**：`triangle`（plink 默认，下三角逐行）、`square`、
   `square0`（右上填 0）、`bin`（R 可直接读的双精度方阵）。文件后缀
   沿用 plink（`.dist`/`.rel`/`.king`/`.genome` 等），`.id` 并列。
3. **公式逐条对照参照实现**，不做"看起来对"的近似（见 01 §2）。
4. **不做的**：`--genome` 的 PPC/IBD 片段推断（需 HMM，另一轮）、
   `--make-grm-sparse`（稀疏输出，待有需求）、`--distance-wts`（加权）。

## 3. 文件

| 文件 | 内容 |
|---|---|
| [01-research.md](01-research.md) | 三家工具的选项与公式核实（file:line），每个矩阵的数学定义 |
| [02-spec.md](02-spec.md) | `reed-distance` 的 API、输出格式、精度契约 |
| [06-roadmap.md](06-roadmap.md) | D1–D4 路线与门 |
| [STATUS.md](STATUS.md) | 实施状态与对拍证据 |
