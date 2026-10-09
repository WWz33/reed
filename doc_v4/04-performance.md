# 04 性能（2026-10-09 实测与下一步）

> 数字见 01-review.md §3；本文给出解释、测量纪律与下一步优化。性能已远超 v1 目标，本轮优化只针对三处"多付的工作"与 reader 吞吐；每项仍以"无收益则撤回"为准。

---

## 1. CPU 账（全量，56 线程）

| 路径 | user CPU | 构成（由对照实验与 `--timing` 推得） |
|---|---|---|
| filter gold `-Oz6` | 455 s | inflate ≈ 60、扫描+掩码+计数+编辑 ≈ 140、deflate 6 级 ≈ 255 |
| filter gold `-Ob` | 519 s | 同上，deflate 换成 BCF 原生编码 + 0 级外的 BGZF 压缩，+64 s |
| filter gold，BCF 输入 | 630 s | 比 VCF 输入多 175 s：BCF→文本渲染（每条约 11 KB 文本）+ 文本重扫 |
| stats 五统计（VCF） | 71 s | inflate ≈ 60 + GT 扫描与累加 ≈ 10 |
| stats 五统计（BCF） | 56 s | inflate ≈ 45（BCF 更大但更规则）+ GT 查表 ≈ 10 |
| stats 边过滤边统计 | 170 s | 比纯统计多 99 s：`FilteredJob` 对每条存活记录 `edit_vcf_line` + 内层重新扫描 |
| reports 两族 | 232 s | 两遍读取各约 110 s |
| index VCF / BCF | 77 / 46 s | inflate + 行/记录定界 |

deflate 仍是 filter 的最大项（56%），级别 6 保持（与 htslib 默认一致，doc_v3 O8）。

---

## 2. 下一步优化（按收益排序）

| # | 项 | 现状 → 预期 | 做法 | 门 |
|---|---|---|---|---|
| O-1 | 报告两族一次读取（01 F4-6） | 5.3 s → ≈3 s | summary 累加器并入 report partial | 两族输出 md5 不变 |
| O-2 | 边过滤边统计不重写文本（01 F4-7） | 170 CPU-s → ≈ 90 | `FilteredJob` 向内层传 `codes`/列 | stats 输出 md5 不变（现有 `stats_site_filters_match_the_filter_then_stats_composition` 测试） |
| O-3 | reader 吞吐（03 §3.4） | BCF stats 3.1 s → ≤ 2.7 s | fadvise + 大缓冲 → 整 chunk 读 → mmap | `--timing` reader 累计；BCF stats ≤ VCF stats（doc_v3 02 §4 的门） |
| O-4 | BCF 输入的 filter 原生解码（01 F4-5） | 12.1 s → ≈ 9.5 s | `RecordPipeline` 接收 `SiteView` + `ScanOutput`，BCF 解码器直接产出（GT 表 + DP/GQ 类型化解码 + INFO 按字典下标取值） | 正文 md5 不变；worker_busy 下降 ≥ 25% |
| O-5 | `--timing` 增加 reader/reducer/链等待的累计（03 A4-5） | — | 引擎 `StageMetrics` 扩展 | 每个基准行附带 |
| O-6 | 扫描第 2 级（AVX2 位图，doc_v3 O6） | 扫描部分 1.5–3× | 运行时分派 | criterion A/B，p < 0.05 |

不做：原生 BCF 手术（B3b，已记录理由）；投机单遍（O9-2，已记录理由）；mimalloc/缓冲池（V1 实测无收益）。

---

## 3. 目标（本机全量，56 线程，负载 < 1，热缓存）

| 命令 | 现在 | 目标 | 依赖 |
|---|---|---|---|
| filter gold `-Oz6` | 8.8 s | ≤ 8.8 s（不回退） | — |
| filter gold，BCF 输入 | 12.1 s | ≤ 9.5 s | O-4 |
| stats 五统计 VCF / BCF | 2.7 / 3.1 s | ≤ 2.7 / ≤ 2.7 s | O-3 |
| stats 边过滤 `--pi` | 3.4 s | ≤ 3.0 s | O-2 |
| reports 两族 | 5.3 s | ≤ 3.2 s | O-1 |
| index VCF / BCF | 2.4 / 2.2 s | 不回退 | — |

---

## 4. 测量纪律（补充 doc_v3 03 §3）

本轮在同一命令、同一二进制、同一参数下测到 56 线程 **29.3 s 与 5.5 s** 两种结果，差异来自机器上其他用户的负载（1 分钟负载 2.97 vs 0.32，13 个登录用户）。
`worker_busy` 两次都约 270 s，说明 CPU 工作量不变、是被抢占。因此：

1. 每次基准前后记录 `uptime` 的 1 分钟负载，负载 > 1 的样本不入 `results.tsv`；
2. 每个数字取 3 次运行的中位数；
3. 判定"是否 CPU 受限"用 `wall ≈ worker_busy ÷ workers`（误差 ≤ 15%），否则标注"受其他因素影响"并附 `--timing`；
4. 线程扩展曲线必须在一次低负载时段内连续完成；
5. `results.tsv` 新增列：`load_1m`、`median_of`、`timing_path`；V0–V6 的数字从 `doc_v3/STATUS.md` 补登（01 F4-11）。

---

## 5. 线程扩展的参考曲线（CPU 受限模型）

全量 `-Oz1` 全部存活：user ≈ 265 s，模型 `wall ≈ 265 / (threads − 2)`：4 → 131 s（实测 131）、8 → 44 s（实测 44.7）、16 → 19 s、32 → 8.8 s、56 → 4.9 s（实测 5.5）。
在 16 核以下的机器上 reed 仍然是 CPU 受限的，deflate 级别对总时长的影响与 56 核时相同比例；文档的"性能提示"应写明：中间文件用 `-Oz1`/`-Ou`，最终交付用默认级别。
