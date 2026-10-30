---
title: 文档
template: splash
hero:
  title: reed
  tagline: 群体遗传学统计与硬过滤，纯 Rust 实现
  actions:
    - text: 快速开始
      link: quickstart/
      icon: right-arrow
      variant: primary
    - text: 命令一览
      link: commands/
      variant: minimal
---

reed 是用 Rust 写的群体遗传学软件，读 VCF/BCF，做硬过滤、QC 报告、窗口化群体统计、主成分与个体间关系矩阵、格式转换和索引。当前版本 0.4.0。

英文页是主文档，中文页是它的镜像；两版不一致时以英文页为准。

二进制叫 `reed`，是 `reed` crate 的一层封装；同一个库也可以从 Rust 代码直接调用，见[库调用](library/)。

## 命令

| 命令 | 别名 | 作用 |
|---|---|---|
| `filter` | | 硬过滤：位点与样本选择、基因型掩码、QC 阈值、LD 剪枝 |
| `stats` | | 位点级与样本级 QC 报告：等位频率、HWE、杂合度、缺失、深度、直方图 |
| `popg` | `popgen` | 窗口化群体统计：π、dxy、Fst、θ、Tajima's D、PBS、SFS、f 统计量、LD 衰减 |
| `pca` | | 样本关系矩阵（GRM）的主成分；`--dapc` 判别分析，`--xval` 交叉验证 |
| `distance` | | 个体间矩阵：等位差距离、IBS、GRM、KING、PI_HAT、相关度、不一致率；`--tree` 建树 |
| `convert` | | 系统发育矩阵、基因型矩阵、hapmap、参考 MSA、投影 SFS；hapmap 反向转 VCF |
| `rename` | `reheader` | 标识符编辑：样本名、contig 名、变异 ID、header 行 |
| `index` | | 建 CSI/TBI 索引 |

## 设计上的几条约定

1. 参数名和语义全库一致：样本子集（`-s/-S/--keep/--remove`）、区域选择（`-r/-R/-t/-T`）、输出（`-o/-O`）在每个命令里拼写相同。
2. 阈值参数都是保留条件：要保证的量写下限（`--maf/--mac/--hwe`），要抑制的量写上限（`--max-maf/--max-mac/--geno/--mind`），命令行上不做 `1−x` 换算。
3. 缺失阈值（`--geno`、`--mind`）的边界值按保留处理：恰好等于阈值时位点或样本留下。
4. 畸形记录跳过并计数，不崩溃；计数进日志或 JSON 报告。
5. 参数在处理第一条记录之前全部校验，出错退出码 2，不产生半截输出。

## 输入输出约定

输入：VCF、bgzip 压缩的 VCF、BCF，或 stdin（`-`）。区域查询（`-r/-R`）需要输入旁边的 `.csi`/`.tbi` 索引（`reed index` 生成）；按位置流式选择（`-t/-T`）不需要索引。

输出：`filter`、`rename`、`convert`、`index` 把 `-o` 当作输出文件名；`stats`、`popg`、`pca`、`distance` 把 `-o` 当作输出前缀，各自追加后缀（`.pi.txt`、`.eigenvec` 等）。零分母的估计量写 `NA`，不写 0。坐标一律 1-based 闭区间（BED 文件输入除外，按 BED 的 0-based 半开约定读）。

`--threads N` 控制工作线程（默认全部核），`--memory MB` 控制在途内存预算（MiB，默认物理内存一半）。长度参数接受 `100000`、`100k`、`100kb`、`1M`、`1Mb`，后缀大小写不敏感。
