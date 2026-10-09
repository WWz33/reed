# reed doc_v3：第三轮复审 + BCF 编解码方案（2026-10-08）

> **后继修订**（见 `docs/README.md`）：doc_v4 修订参数方向（`--geno/--mind` → `--max-missing` 家族）
> 与 CLI/文档一致性；doc_v3 04 §C3 的"`TagTable` 未做"在 `docs/rounds/v3/STATUS.md` 更正为**已做**（doc_v4 01 F4-2）。
> V0–V6 的实施结果仍以 `docs/rounds/v3/STATUS.md` 为准。

> 延续 doc_v2。doc_v2 的规范（CLI、对齐原则、架构）继续有效；本目录给出 R0–R8 落地后的复审、R6 剩余的 BCF 方案、
> 性能与清理计划。不修改代码。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 第三轮复审：状态核实、全量性能实测（对照实验分离 CPU 去向）、发现清单（带 `file:line` 与证据） |
| [02-bcf.md](02-bcf.md) | **BCF 编解码层**：格式要点（BCFv2.2/htslib）、长度前缀链并行定界、格式无关的消费接口、GT 查表解码、四种输入输出组合、字典与 `IDX` 规则、测试门、分期 B0–B4 |
| [03-performance.md](03-performance.md) | 性能优化：消除 4 处串行段（O1–O5）、降低每条记录 CPU（O6–O12）、基准协议补充、目标表 |
| [04-cleanup.md](04-cleanup.md) | 语义修正（`--mind`、样本规则互斥）、结构清理（清单残留、reed-stats 公开面、reed-filter 纯化、QC 入口合并、类型化错误）、CLI 一致性清单 |
| [05-roadmap.md](05-roadmap.md) | V0–V6 的顺序与验收门 |

## 结论速览

| 维度 | 结论 |
|---|---|
| 完成度 | doc_v2 的 R0–R8 除 BCF 外基本完成；P0 缺陷全部修复并有回归测试；测试连续 3 次全绿 |
| 正确性 | 过滤与 gold 一致；π/dxy/θ/Tajima/Fst 与 pixy 逐位一致；**新发现两处语义缺陷**：`--mind` 的实现与声明的 plink2 语义不符；样本规则互斥导致 `--mind`/`--keep` 被静默忽略 |
| 性能 | 并行路径 CPU 受限（wall ≈ CPU ÷ 54）；deflate 6 级占 filter CPU 的 56%（默认值不改，与 htslib 一致）；**4 处串行段**：`--hardy` 13.5 s、全部存活的 `-Ov` 13.2 s、`reed index` 31–49 s、stdin 整体落盘 |
| 资源 | stdin 落盘文件从不删除（/tmp 已残留 5 个） |
| BCF | 方案：长度前缀链解决无换行符的并行定界；先抽象格式无关的消费接口（B0），再接入解码（B1）、渲染（B3a）、编码（B2，先用 noodles）、原生手术（B3b）；`clean.bcf` 与 `clean.vcf.gz` 是同一份数据，可做全量跨格式门 |
| 对齐残留 | `##reed_ascertainment` 仍写出；reed-stats 仍导出已撤下的方法；reed-filter 仍依赖 noodles；`reed qc` 与 `stats --summary` 两个入口；表达式缺"陈旧 INFO"警告 |

## 需要你确认

| 项 | 推荐 |
|---|---|
| `--mind` 改为 plink2 语义（C0） | 是：参数名来自 plink2，就按 plink2 计算；gold 继续用 `-S` 复现原始口径 |
| `reed qc` 子命令 | 并入 `reed stats --summary`，保留一个版本的隐藏别名 |
| BCF 编码先用 noodles（B2） | 是：先保证正确，profile 证明是瓶颈后再写原生编码器（B4） |
| reed-stats 中没有 CLI 出口的方法 | 改为 `pub(crate)`，接入 CLI 并通过对拍后再公开 |
