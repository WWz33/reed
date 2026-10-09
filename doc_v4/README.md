# reed doc_v4：第四轮全项目复审 + 参数语义统一规范（2026-10-09）

> 延续 doc_v2（规范）与 doc_v3（BCF/性能/清理）。doc_v3 的实施结果记录在 `doc_v3/STATUS.md`，本目录不重复。
> 本目录给出 V0–V6 落地之后的**全项目重新审阅**，以及一条新的用户裁决：**参数语义以"保留条件"统一方向，不必全部照搬主流软件**。
> 本目录不修改代码；实现者按 06-roadmap.md 执行并把结果写进 `doc_v4/STATUS.md`。

## 用户裁决（2026-10-09）

1. **参数语义必须对齐一致**，优先于与某个主流软件的名字一致。
2. 缺失率阈值采用 vcftools `--max-missing 0.8` 的方向：**值表示必须存在（被调用）的比例**，而不是 plink2 `--geno 0.2` 的"允许缺失的比例"。
3. `--maf 0.05` 保持现状（keep if MAF ≥ 0.05）。

由此得到贯穿本目录的原则：**所有阈值参数都用"保留条件"表述**——下限写作 `≥`（`--maf`、`--mac`、`--hwe`、`--max-missing`），上限写作 `≤`（`--max-maf`、`--max-mac`、`-Q`）。
当前 CLI 中唯一违反这条的是 `--geno` 与 `--mind`（"丢弃条件"），02-cli-semantics.md 把它们改掉。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 全项目复审：状态核实（构建/测试/clippy）、正确性（gold、跨格式、三工具边界对拍）、全量性能实测与方法论修正、发现清单（带 `file:line`、证据、修法） |
| [02-cli-semantics.md](02-cli-semantics.md) | **参数语义统一规范**：方向规则、每个阈值的分子/分母/阶段/精确整数判定、`--geno/--mind` → `--max-missing/--max-missing-sample` 的切换、边界等价性证明与测试、CLI 全表、迁移 |
| [03-architecture.md](03-architecture.md) | 架构现状与残留：crate 图与纯度、引擎（票据/字节预算/长度前缀链）、输入输出层、文档与代码的漂移、命名残留 |
| [04-performance.md](04-performance.md) | 今日实测（56 线程、负载 0.3）、BCF 路径的 CPU 账、线程扩展、测量纪律（共享机器的污染样本）、下一步优化与目标 |
| [05-algorithms.md](05-algorithms.md) | 统计与过滤语义的对齐状态：已逐位对齐的、尚未接入 CLI 的、未做的（多等位逐等位重算、SFS 投影 CLI、LD、f 统计）、gold/popgen 的缺口 |
| [06-roadmap.md](06-roadmap.md) | W0–W6 顺序与验收门 |
| [STATUS.md](STATUS.md) | 实施状态（实现者填写，格式沿用 doc_v3/STATUS.md） |

## 结论速览

| 维度 | 结论 |
|---|---|
| 完成度 | doc_v3 的 V0–V5 全部落地，V6 的原生 BCF 手术按文档"可选"未做；`cargo test` 344 通过（连续 3 次）；release 构建 0 警告；clippy 0 error |
| 正确性 | gold 等价命令 3,315,604 → 1,768,759 与 gold 一致；`-Ob` 正文 ≡ `-Oz` 正文；BCF 输入正文 ≡ VCF 输入正文；BCF/VCF 输入的 stats 输出 md5 相同；**`--geno 0.1` 在掩码后与 vcftools `--max-missing 0.9`、plink2 `--geno 0.1` 的位点集合逐个相同（263,420 / 381,672）** |
| 参数语义 | 两处"丢弃方向"的阈值（`--geno`、`--mind`）与其余"保留方向"的阈值不一致——这是本轮的首要修正（02） |
| 性能（56 线程） | filter gold 规则 `-Oz6` 8.8 s、`-Ob` 9.9 s；BCF 输入 12.1 s（渲染成文本再扫描，CPU +38%）；stats 五统计 2.7 s（BCF 3.1 s，reader 受限）；`reed index` 2.4 s（VCF）/2.2 s（BCF） |
| 新发现 | `--verbosity` 在 stats 中被静默忽略；STATUS 对 C3 的记录与代码不符（reed-filter 已不依赖 noodles）；`reed index -` 对 stdin 直接报 "No such file"；`##reed_command` 写入绝对路径的 argv；BCF 输入的 filter 多付一遍渲染+重扫；`stats` 报告模式两族同时请求时读两遍 |
| 测量纪律 | 共享机器上同一命令 56 线程出现过 5.5 s 与 29 s 两种结果（负载 0.3 vs 3.0）；基准必须记录负载并取 3 次中位数 |

## 需要你确认

| 项 | 推荐 | 备选 |
|---|---|---|
| 位点缺失阈值的名字 | `--max-missing FLOAT`（vcftools 名字，值 = 必须被调用的比例） | `--min-call-rate FLOAT`（自解释，但不是你习惯的 "miss 0.8"） |
| 样本缺失阈值的名字 | `--max-missing-sample FLOAT`（同方向、同分母定义） | `--min-sample-call-rate FLOAT` |
| `--geno`/`--mind` 的去向 | **删除**，输入时报错并给出换算后的新参数（避免静默反向） | 保留为隐藏别名并自动换算（有静默反向风险） |
| 缺失率的分子 | **按等位基因计**（vcftools 语义；本数据没有半调用，与按基因型计的结果相同） | 按基因型计（plink2 语义） |
