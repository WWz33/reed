# reed doc_v21：真实数据的复审 + 修复（2026-10-28）

> 触发：仓库外评审报告 `reed-rs-review-2026-10-27.md`（HEAD `fa8131f`）。
> 方法：**只读复审该报告**，但把每一处"复现"从合成小文件换成仓库里真实的
> `subvcf.vcf.gz`（348 样本 / 3000 记录，20 contig 头、只有 `Chr05` 有记录）
> 与 `clean.vcf.gz`（7.0 GB WGS 子集），逐条在本机构建 `reed 0.4.0` 重跑。
> 结论：报告的 10 条 must-fix **全部属实**；但其中两条结论有出入（见下），
> 且用真实数据一测，就翻出**一条报告明确判为"误报"的真 bug**。
>
> 本目录：README（本轮总览+决策）、[01-review.md](01-review.md)（对评审的核实
> 与更正、新发现的 M3）、[STATUS.md](STATUS.md)（M1–M11 + S 系列修复表）、
> [06-roadmap.md](06-roadmap.md)（未修条目）。

## 结论速览

| 维度 | 结论 |
|---|---|
| 评审可信度 | must-fix 10 条全部已在真实数据上复现；should-fix 抽查逐条属实 |
| 报告的出入（已更正） | ① #2 的"CSI is unaffected"不成立（`.csi` 同样失败）；② 报告判为"不可复现/误报"的 `input ends inside a record` 是真 bug（本轮 M3） |
| 真实数据揭示的新 bug | **M3**：区域读在 BGZF 块跨行时截断报错——真实 `subvcf` 上 40 个真实区间里 36 个失败（90%），`clean.vcf.gz` 上 0 个（该文件恰好块块行对齐，正是报告作者只测它才漏掉的原因） |
| 修复 | M1–M11（must-fix）11 条 + S 系列 bug 类 12 条（命令行静默 no-op/坏输出/口径），全部"先写失败测试再修" |
| 文档 | 更正 `docs/spec/cli.md`、`README.md` 的偏差（`--hwe midp`、`-W`、BED、`--max-mac 0`、`--indep-pairwise` r² 范围等），并补齐本轮文档 |
| 基线 | `cargo test --workspace` exit 0（全部目标 0 失败）；`cargo clippy --workspace --all-targets` 0 告警 |

## 关键更正

1. **M3 是评审漏判的真 bug。** 评审第 299 行称 `-r` 的
   `input ends inside a record` "65 个随机区域无法复现，Treat that reviewer's
   Must-fix 1 as a false positive"。在真实 `subvcf.vcf.gz` 上按真实记录区间取
   40 个，36 个失败。机理：`RegionChunkSource` 以"整块"为单位读取，块边界
   跨行时，最后一个块含半个记录，引擎的 `stitcher.finish()` 把它当截断。
   决定性实验：同一份文本、同一个 reed 自建索引，仅 BGZF 切块方式不同——
   跨行切必失败、行对齐切必通过。
2. **评审 #2 的"CSI unaffected"有误。** `.tbi`（1 ref，首现顺序）与
   `bcftools index -c` 重建的 `.csi`（同样只含出现的 contig）都以
   `invalid reference sequence ID` 失败；只有头与索引 ref 顺序一致时才通过。
   修法必须对 `.tbi` 和 `.csi` 一起用索引自身的名字表。

## 决策

| 决策 | 状态 |
|---|---|
| `--memory` 单位 | 已核实是 **MiB**（`m << 20`）；报告 #10 的数值复现需按 MiB 读 |
| M9 UPGMA 修法 | 已定：`(survivor, gone) = (older, newer)`，与 R `hclust` 一致；本机整数矩阵穷举 25/25 命中（旧规则 19/25） |
| M9 真实矩阵残余 | `reed` 与 R 在本机 348 样本矩阵上 clade 组成仍差 100/237——归因于 R 的最近邻扫描 tie 序与全对扫描不同，非槽位规则；**记录为残余，不在本轮强求逐位一致** |
| S16 `-W` 校验位置 | 已定：放在 **CLI 翻译层**，不进库——库默认 `WriteIndex::Csi`，库级拒绝会伤到未指定压缩的库调用方 |
| S29 `--snps-only` 语义 | 已定：**裸标志 = 任意 1-bp 替换（bcftools）**，只有 `=just-acgt` 限制 A/C/G/T；原 `default_missing_value="just-acgt"` 使两者等价，已改 |
| should-fix 范围 | 本轮清 **bug 类**（静默 no-op / 坏输出 / 口径分裂）；验证面/CI/lean 仍留 roadmap |

## 标记约定

【据】有实测；【更正】对评审结论的修正；**M#** 本轮 must-fix 编号；**S#** 沿用评审 should-fix 编号。
