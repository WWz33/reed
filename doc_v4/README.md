# reed doc_v4：第四轮全项目复审 + 参数语义统一规范（2026-10-09）

> 延续 doc_v2（规范）与 doc_v3（BCF/性能/清理）。doc_v3 的实施结果记录在 `doc_v3/STATUS.md`，本目录不重复。
> 本目录给出 V0–V6 落地之后的**全项目重新审阅**；参数部分的命名在 2026-10-09 晚间经过一次回退（见下第 3 条）。
> 本目录不修改代码；实现者按 06-roadmap.md 执行并把结果写进 `doc_v4/STATUS.md`。
>
> **后继修订**：2026-10-09 晚的第五轮复审在 `doc_v5/`。本目录的 02（CLI）、05 §2（统计出口）、06（路线 W5/W6）、07（目录树）已被 doc_v5 的 02 / 05 / 07 / 06 取代；W0–W4 的实施记录仍以本目录的 STATUS 为准。

## 用户裁决（2026-10-09）

1. **参数语义必须对齐一致**，优先于与某个主流软件的名字一致。
2. `--maf 0.05` 保持现状（keep if MAF ≥ 0.05）。
3. **缺失阈值只有位点与样本两级，`--geno`/`--mind` 保留原样**（plink2 短名字，值 = 允许的缺失比例）；撤回首版的 `--max-missing`、`--max-missing-sample`、`--max-missing-count`。理由：这两个名字短、与 plink2/plink 1.9 的习惯一致，而位点与样本两级本来就自洽（同为缺失率、同为"上限 `≤`"、同为等位计、同取掩码后）。

由此得到贯穿本目录的原则：**所有阈值参数都用"保留条件"表述**——要保证的量写下限（`--maf`、`--mac`、`--hwe` 用 `≥`），要抑制的量写上限（`--max-maf`、`--max-mac`、`-Q`，以及 **`--geno`、`--mind`** 用 `≤`）。首版把 `--geno/--mind` 误判成"方向例外"，才引出改名；实际上它们与 `--max-maf` 同属 `≤` 家族，方向本来一致。详见 02-cli-semantics.md §1。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 全项目复审：状态核实（构建/测试/clippy）、正确性（gold、跨格式、三工具边界对拍）、全量性能实测与方法论修正、发现清单（带 `file:line`、证据、修法） |
| [02-cli-semantics.md](02-cli-semantics.md) | **参数语义统一规范（命名回退版）**：方向规则、每个阈值的分子/分母/阶段/精确整数判定、`--geno`/`--mind` 的精确定义、撤回 `--max-missing` 家族、边界等价性证明与测试、CLI 全表、W0b 实施清单 |
| [03-architecture.md](03-architecture.md) | 架构现状与残留：crate 图与纯度、引擎（票据/字节预算/长度前缀链）、输入输出层、文档与代码的漂移、命名残留 |
| [04-performance.md](04-performance.md) | 今日实测（56 线程、负载 0.3）、BCF 路径的 CPU 账、线程扩展、测量纪律（共享机器的污染样本）、下一步优化与目标 |
| [05-algorithms.md](05-algorithms.md) | 统计与过滤语义的对齐状态：已逐位对齐的、尚未接入 CLI 的、未做的（多等位逐等位重算、SFS 投影 CLI、LD、f 统计）、gold/popgen 的缺口 |
| [06-roadmap.md](06-roadmap.md) | W0–W6 顺序与验收门（W0b = CLI 拼写回退，先于 W1） |
| [07-layout.md](07-layout.md) | **目录树规范**：现状（含 79 G 生成物、三个文档根、版本后缀 crate）、目标树（`reed-rs/` 6 个顶层目录 + `~/reed/` 四分法）、现状→目标映射与 `git mv` 命令、要同步改的引用、注意事项 |
| [STATUS.md](STATUS.md) | 实施状态（实现者填写，格式沿用 doc_v3/STATUS.md）；W0 表中涉及 CLI 拼写的行由 W0b 行取代 |

## 结论速览

| 维度 | 结论 |
|---|---|
| 完成度 | doc_v3 的 V0–V5 全部落地，V6 的原生 BCF 手术按文档"可选"未做；`cargo test` 344 通过（连续 3 次）；release 构建 0 警告；clippy 0 error |
| 正确性 | gold 等价命令 3,315,604 → 1,768,759 与 gold 一致；`-Ob` 正文 ≡ `-Oz` 正文；BCF 输入正文 ≡ VCF 输入正文；BCF/VCF 输入的 stats 输出 md5 相同；**`--geno 0.1` 在掩码后与 vcftools `--max-missing 0.9`、plink2 `--geno 0.1` 的位点集合逐个相同（263,420 / 381,672）** |
| 参数语义 | 首版判"`--geno`/`--mind` 方向不一致"并改名；2026-10-09 晚裁决**撤回改名**（这两个短名字保留，方向本就属 `≤` 家族），改为把语义写成条文：分子、分母、阶段、精确边界（02 §3/§4）。W0 的实现层修正全部保留 |
| 性能（56 线程） | filter gold 规则 `-Oz6` 8.8 s、`-Ob` 9.9 s；BCF 输入 12.1 s（渲染成文本再扫描，CPU +38%）；stats 五统计 2.7 s（BCF 3.1 s，reader 受限）；`reed index` 2.4 s（VCF）/2.2 s（BCF） |
| 新发现 | `--verbosity` 在 stats 中被静默忽略；STATUS 对 C3 的记录与代码不符（reed-filter 已不依赖 noodles）；`reed index -` 对 stdin 直接报 "No such file"；`##reed_command` 写入绝对路径的 argv；BCF 输入的 filter 多付一遍渲染+重扫；`stats` 报告模式两族同时请求时读两遍 |
| 测量纪律 | 共享机器上同一命令 56 线程出现过 5.5 s 与 29 s 两种结果（负载 0.3 vs 3.0）；基准必须记录负载并取 3 次中位数 |

## 需要你确认

| 项 | 推荐 | 备选 |
|---|---|---|
| 位点缺失阈值的名字 | **`--geno FLOAT`**（plink2 短名；值 = 允许的缺失比例，`--geno 0.2` ⟺ 调用率 ≥ 0.8） | ~~`--max-missing FLOAT`~~（首版方案，已撤回） |
| 样本缺失阈值的名字 | **`--mind FLOAT`**（同上；仅 `reed filter`） | ~~`--max-missing-sample FLOAT`~~（已撤回） |
| `--max-missing*` 的去向 | 解析期报错 + **精确十进制**换算提示（`--max-missing 0.9` → `--geno 0.1`），不提供静默别名 | 完全忽略（无提示；会让人以为参数不存在） |
| 计数形式 | **删除** `--max-missing-count`（本族定为位点/样本两个参数） | 将来确有需求时加 `--geno-count INT` |
| 缺失率的分子 | **按等位基因计**（vcftools 语义；本数据 0 个半调用，与按基因型计结果相同） | 按基因型计（plink2 语义） |
| 多群体 Fst 的出口（05 §2.1） | **已裁决**：`--fst [wc84\|hudson]` 不变，新增布尔开关 `--fst-global`（三个群体以上时算合并而非成对；默认成对；不能与 `hudson` 组合；两个群体时等于成对值，只提示）。`--fst-pops` 作废 | 备选拼法 `--fst-scope pairwise\|global` |
| 目录树（07） | 按 07 §2 的树 `git mv`；生成物进 `target/`；参照仓库与数据移到 `~/reed/{refs,data}` | 只做 07 §4 的仓库内部分，仓库外暂不动 |
