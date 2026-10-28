# 06 后续路线（v21 之后）

> v21 清掉了 must-fix M1–M11 与评审的 bug 类 should-fix（S11–S22、S29、S30、
> S32–S34）。以下余量按优先级记录；来源＝评审 2026-10-27 报告与 v20 06 的未完项。

## P0：会撒谎的门 / 静默错果

| 项 | 内容 | 备注 |
|---|---|---|
| CI-1 | data-gates 不构建被测二进制、`run_all.sh` 不调 `run_reed.sh`（旧 v20 CI-1 未修） | 与 v20 同一条 |
| M9 残余 | reed UPGMA 的**最近邻扫描 tie 序**与 R `hclust` 不同（槽位已对齐，扫描未对齐）；真实矩阵 clade 仍差 100/237 | 需按 R `hcass2` 复刻扫描顺序，或明确以槽位规则为准并文档声明 |
| S23 | `PopgenRequest.regions` 被库入口忽略（`run_popgen` 用 `input::open`，不读 `regions`） | 库 API 契约 |
| S24 | `convert` 的 SFS-only 跑仍逐样本词法化 GT（`has_row_targets` 只用于是否 sanitize） | 纯开销 |
| S25 | f-stats/LD 每位点分配 contig 名（`jackknife.rs`、`sink.rs`）；`--ld-pairs` 不可内存容纳 | 规模 |
| S28 | `N_MISSING`/`F_MISSING` 表达式只数全缺失，`-g miss` 数半call——同 run 口径分裂 | 语义 |
| S31 | 区域谓词每条记录线性扫全部区间（M4 修好 BED 后更吃紧） | 需按 contig 分组/排序/二分 |

## P1：验证面与文档

| 项 | 内容 |
|---|---|
| S26/S27 | sklearn/UMAP/t-SNE 数值差距与参考锚点缺失；PPCA 文档写 EM 实为 ALS；isomap NaN 排序点；断连 kNN 图静默填最大距离 |
| S35–S53 | nice-to-have：`--threads` 未约束 eigensolve、`grm` 重复标准化矩阵、`--bed-file` 每区深拷贝累加器、`decay_range` 短块提前返回、若干公共 API 潜在 bug、死代码、`--ld-pairs` 缓存按条目非字节、convert 默认 prefix 丢目录、MSRV 1.91 vs 钉 1.99 |

## P2：测试金字塔与 lean

- `QcReportAcc::merge` 缺定律测试；`empty_group_is_rejected_at_build_time` 伪造期望错；
- `reed-filter`/`reed-io`/`reed` 若干子模块（`--hwe` 规则、head/tail 单元、region 源）补单元；
- 多区域 `-r` 折叠成 min..max 全量扫描（间断 seek）；
- BCF `-r` POS 谓词 vs bcftools 重叠语义（改谓词或文档注明）。

## 明确不做

- 不引入新 CLI 形态、不改统计定义；
- 不为逐位复刻 R 的 tie 扫描顺序而重写 UPGMA（除非确有用户对拍需求）。
