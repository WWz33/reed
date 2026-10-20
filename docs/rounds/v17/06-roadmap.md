# doc_v17 / 06 路线：D1–D4

| 批次 | 内容 | 门 | 状态 |
|---|---|---|---|
| **D1** 调研与规范 | 三个子代理用 codegraph 探索 plink/vcftools/bcftools；01/02 落定 | 每个矩阵有 file:line 与公式 | **完成** |
| **D2** crate 与迁移 | 新建 `reed-distance`；`Mat` 下沉 `reed-core`；GRM 从 `reed-reduce::pca` 迁入 | 双 crate 测试全绿；PCA 数值不变 | **完成** |
| **D3** 矩阵族实现 | 10 种矩阵 + 4 种输出格式 + CLI `reed distance` | 与参照工具对拍 | **完成**（PI_HAT 部分） |
| **D4** 对拍与收口 | 逐矩阵数值对拍、文档、CHANGELOG | 见 STATUS 的对拍表 | **完成** |

## 本轮不做

见 `01-research.md` §8 的"不做"表：`--genome` 的 PPC/IBD 片段（需 HMM）、
`--distinctive-sites`（逐站点非矩阵）、`--make-grm-sparse`、`--distance-wts`、
vcftools `--geno-chisq`（SNP×SNP）。

## 风险与已发生的问题

| 风险 | 处理 |
|---|---|
| `reed-reduce` 与 `reed-distance` 循环依赖 | 把 `Mat` 下沉到 `reed-core`，两 crate 都依赖它 |
| 两个 crate 各有一套标准化会漂移 | 标准化只在 `reed-distance`，`reed-reduce::pca` re-export |
| plink 的 PI_HAT 用 Nei 偏差修正，与理论期望值不同 | **如实标记未对齐**（见 STATUS §4），不声称匹配 |
| plink 的 triangle 含对角线 | 实测发现后修正（原实现漏了第 0 行） |
