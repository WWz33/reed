# 05 统计与过滤语义：对齐状态

> 原则不变（doc_v2 04-algorithms.md）：经典方法 + pixy 的缺失数据处理，参数按科学方法命名，软件只在测试中作参照物。
> 本文是对齐状态的台账：哪些已逐位对齐、哪些有实现但没有 CLI 出口、哪些未做。

---

## 1. 已对齐（有参照物、有数字）

| 统计量 / 行为 | 参照物 | 对拍结果 | 来源 |
|---|---|---|---|
| π、dxy（含不变位点，默认） | pixy 2.2.4 | 全量 28,470 行，最大相对误差 4.8e-14 | doc_v3 STATUS R7 |
| π，`--variant-only` | vcftools `--window-pi` | 最大相对差 1.05e-5（vcftools 的 n/(n−1) 校正） | 同上 |
| Watterson θ、Tajima's D | pixy 2.2.4 | 3.1e-14 / 8.9e-11 | 同上 |
| Fst wc84 / hudson | pixy 2.2.4 | 最大绝对差 2.9e-15 / 2.2e-15 | 同上 |
| PBS（独立 Hudson 对） | scikit-allel | 逐列一致 | 同上 |
| `--freq/--hardy/--het/--missing` | plink2 | 逐行一致（键 = 染色体编号, 位置） | 同上 |
| `--depth` | vcftools | 逐行一致 | 同上 |
| `--summary` | bcftools stats | 九节逐行一致 | doc_v3 STATUS |
| 过滤：gold 规则 | bcftools（gold） | 1,768,759 位点，GT/AC/AN/NS 一致，`check_parity.py` PASS | 本轮复测 |
| 过滤：掩码后缺失阈值 | vcftools、plink2 | 263,420 / 381,672 位点集合相同 | 本轮（01 §2） |
| 过滤：`-i/-e` 表达式、`FMT/` + `--set-gts .`、`-q/-Q` | bcftools | 位点集合 / GT md5 一致 | doc_v3 STATUS |
| 过滤：`--mind`（plink2 口径） | plink2 | fixture 样本集合一致 | doc_v3 STATUS C0 |
| BCF 渲染 | `bcftools view -H` | 全量 3,315,604 条 sha256 一致 | doc_v3 STATUS |
| BCF 输出 | bcftools 读回 | 正文 sha256 == VCF 输出 | 本轮复测 |

---

## 2. 有实现、无 CLI 出口（`#[doc(hidden)]`，不在文档化 API 中）

| 累加器 | 位置 | 参照物（未建） | 接入条件 |
|---|---|---|---|
| 投影 SFS（`ProjectedSfsAcc`、`projection_weights`） | `reed-stats/src/sfs.rs` | easySFS（`--proj`） | `--sfs --proj N[,N…]`，dadi/fastsimcoal2 文件格式，`gold/popgen/run_easysfs.sh` |
| 联合 SFS（`JointSfsAcc`） | `reed-stats/src/individual.rs` | easySFS 二维谱 | 同上 |
| 合并形式（r 群体）WC84（`WcMultiAcc`） | `reed-stats/src/fst.rs` | vcftools `--weir-fst-pop` ×r（**合并**，weighted Fst）、hierfstat `wc()`、scikit-allel `weir_cockerham_fst(subpops)` | 见 §2.1；**待裁决**：`--fst wc84-global`（推荐）或 `--fst-pops A,B,C`（doc_v2 原案） |
| 成对矩阵（`FstMatrixAcc`） | `reed-stats/src/fst.rs` | pixy（全部群体对）、plink2 `--fst`（默认全部对） | **不需要新参数**：`--fst` 已经是"所有群体对"（`cli.rs:222`；`request.rs:40-49` 逐对建 `FstAcc`，群体来自 `--pop`）。矩阵累加器只是把"逐对累加"换成"一遍建成"，是实现替换 |
| LD（`LdDecayAcc`、`LdPruner`、`pair_sums`） | `reed-stats/src/ld.rs` | plink2 `--r2-unphased`、`--indep-pairwise`；vcftools `--geno-r2` | `--ld-decay`、`--indep-pairwise`；需要平面存储（整 contig） |
| Patterson f2/f3/f4、D、f4-ratio（`F2Acc/F3Acc/F4Acc`、`f4_ratio`） | `reed-stats/src/patterson.rs` | scikit-allel `admixture.py` | `--f2/--f3/--f4/--dstat/--f4-ratio`，block jackknife 输出 |
| Block jackknife（`BlockJackknife`） | `reed-stats/src/jackknife.rs` | scikit-allel `average_patterson_d` | 随 f 统计接入 |

接入顺序建议：SFS（需求最常见、参照物最简单）→ `--fst` 的合并形式（名字见上表，待裁决）→ f 统计 → LD。每项接入时把 `#[doc(hidden)]` 去掉，并在 `gold/popgen` 增加对应的 `run_<tool>.sh` 与 compare 脚本。

### 2.1 多群体 Fst 的两种形式（`--fst-pops` 到底是什么，供裁决）

**定义。** 给定 r ≥ 2 个群体，Fst 有两种算法，不是同一个数：

| 形式 | 算什么 | 估计量 | reed 现状 |
|---|---|---|---|
| **成对**（pairwise） | 每个群体对 (i, j) 各一个 Fst，r 个群体给 r(r−1)/2 个值 | WC84（r = 2）或 Hudson（Bhatia 2013 比值之和） | `--fst wc84\|hudson`，群体来自 `--pop`，逐对建 `FstAcc`（`request.rs:40-49`），输出 pixy `_fst.txt` 布局，已与 pixy 逐位对拍 |
| **合并**（r 群体，常称 global/overall Fst） | 把 r 个群体放进**同一个** WC84 方差分析，得到**一个** θ = Σa / (Σa + Σb + Σc)；a = 群体间、b = 群体内个体间、c = 个体内的方差分量，按 r 群体公式用 n̄、n_c、s²、h̄ 计算 | 只对 WC84 有定义；Hudson 本质是两群体量，没有合并形式 | `WcMultiAcc`（`fst.rs:249-343`），按 scikit-allel `weir_cockerham_fst` 的公式：n_i = an_i // 2，h̄ 由 n_het 得到；任一群体 an_i = 0、n̄ ≤ 1 或 n_c ≤ 0 的位点跳过并计数；`m7_fst.rs` 用独立参考实现验过 r = 3 的 a、b、c 与比值（误差 < 1e-12）。**无 CLI 出口** |

合并值**不是**成对值的平均：分量在 r 个群体上一起求，分母里 n_c、s² 都依赖整个群体集合。

**三家参照物对 ≥ 3 个群体的做法（源码核实）：**

| 工具 | 做法 | 证据 |
|---|---|---|
| vcftools `--weir-fst-pop` 给 r 个文件 | **合并**：每位点一个 `WEIR_AND_COCKERHAM_FST`；日志给 weighted Fst = Σa / Σ(a+b+c) 与 mean Fst | `variant_file_output.cpp:3466` 起：`N_pops` 循环、`r = N_pops`、`nbar`、`nc`，一列输出 |
| pixy `--stats fst` | **成对**：`combinations(pop_names, 2)` | `pixy/stats/summary.py:384` |
| plink2 `--fst` | **成对**："Estimate Wright's Fst between pairs of populations"；默认全部对，`base=` / `ids=` 限制对 | `plink2_help.cc:1061`；`plink2_misc.cc:13612-13665` 建 `pop_pairs` |

所以"`--fst` 与 pixy、vcftools 类似"只在**两个群体**时成立。三个以上群体：pixy、plink2 给成对值（reed 的 `--fst` 即此），vcftools 给的是合并值——`--fst-pops` 要补的正是 vcftools 这一种。

**什么时候要合并形式：**
- 要一个数回答"这 r 个群体之间的分化占总方差多少"（物种或区域水平的 global Fst，层级设计里的一层）；
- 全基因组扫描时每个窗口/位点一个"所有群体同时"的分化值，比 r(r−1)/2 列成对值更便于找 outlier（vcftools 的逐位点输出就是这种用法）；
- 复现 vcftools、hierfstat 的结果。
不适合回答"谁和谁分化"——那是成对矩阵的事。

**与 `--fst` 的关系与对拍：**
- 群体来源：成对用 `--pop` 的全部群体。合并形式要么也用 `--pop` 全部（= vcftools 行为），要么在参数里点名子集（= `--pbs A,B,C` 的样子），这是下面两案的分水岭。
- 输出：成对 = 每对一行（pixy 布局）；合并 = 每窗口一行，列 `pops`（逗号连接的群体名）、`no_snps`、`fst`，另附 a、b、c 分量和便于对拍。
- 对拍：vcftools r 个 `--weir-fst-pop` 的 weighted Fst 直接等于 `WcMultiAcc::value()`——无半调用、双等位时逐位一致（vcftools 的 n_i 用有基因型的个体数，reed 用 an_i // 2，无半调用时相同；vcftools 支持多等位，`WcMultiAcc` 是双等位，多等位位点先被默认门剔除）。对拍脚本尚未写（§4）。

**两种出口（已裁决：取 A 的"单独开关"变体，见下文定稿形态）：**

| | A：`--fst wc84-global`（推荐） | B：`--fst-pops A,B,C`（doc_v2 02 §5.2 原案） |
|---|---|---|
| Fst 参数个数 | 一个 | 两个 |
| 群体来源 | `--pop` 全部（与 vcftools 相同） | 参数里点名（与 `--pbs A,B,C` 相同） |
| 要子集 | 改 pop 文件；将来若需要，为**所有**统计量统一加一个群体子集选择器 | 直接写在参数里 |
| 与 `hudson` 的关系 | 方法槽里天然只有 wc84 有 global 形态，无需额外规则 | 必须规定 `--fst-pops` 只用 WC84，并定义与 `--fst hudson` 同时出现时的行为 |
| 与 `--pi/--dxy/--theta` 的一致性 | 一致：统计量参数不点名群体 | 多一种点名机制 |
| 一次运行同时要成对与合并 | 不能（一个方法槽），分两次跑；若常用再让 `--fst` 收多个值（`--fst wc84,wc84-global`） | 能 |

推荐 A：少一个参数，估计量差异进方法槽，不与 `--pop` 争"选群体"的角色；"Hudson 没有合并形式"也说明合并是 WC84 的一个形态，而不是独立的统计量。选 B 的理由只有一个：想在参数里点名子集并与成对矩阵同跑。

**定稿形态（用户裁决：用一个单独的开关选成对/全局，不塞进 `--fst` 的值里）：** `--fst [wc84|hudson]` 不变，新增布尔开关 **`--fst-global`**，与 `--pbs-normed` 同类（修饰某个统计量的布尔开关）。

- 默认（不给 `--fst-global`）：成对，全部群体对。现状不变。
- `--fst --fst-global`：合并（r 群体 WC84），群体为 `--pop` 的全部，一个值；写 `PREFIX_fst_global.txt`（每窗口一行：`pops`、`chromosome`、`window_pos_1`、`window_pos_2`、`avg_wc_fst`、`no_snps`）。成对与合并二选一，两种都要就跑两次。
- 只对三个及以上群体有意义：两个群体时合并值等于 `--fst wc84` 的成对值（WC84 r = 2 就是成对估计量，`m7_fst.rs` 已锁定），此时打印一条提示而不报错。
- `--fst-global` 没有 `--fst` → 解析期报错（clap `requires = "fst"`）；`--fst hudson --fst-global` → 报错 `hudson Fst has no r-population (global) form; use wc84`。
- `StatKey`：成对仍为 `fst.wc84.0.1`，合并为 `fst.wc84.global`。
- 备选拼法 `--fst-scope pairwise|global`（值开关，默认 pairwise）：更显式但更长；用户偏好短开关，取 `--fst-global`。

B 案（`--fst-pops`）随之不需要；`lib.rs:47` 的注释与 doc_v2 02 §5.2 的 `--fst-pops` 行在 W5 一并清掉。上面实现者注里的键名 `fst.wc84-global` 以本段的 `fst.wc84.global` 为准。

**实现者注（定了之后归 W5）：** `StatKind::FstGlobal { groups }` + `WcMultiAcc`；`StatKey` `fst.wc84-global`；每窗口一行的输出；`gold/popgen/run_vcftools_fst.sh`（3 群体 weighted Fst 对拍）；去掉 `#[doc(hidden)]`；更新 `reed-stats/src/lib.rs:47` 的注释（现在写的是 `--fst-pops`）。

---

## 3. 未做

| 项 | 现状 | 风险 | 处理 |
|---|---|---|---|
| 多等位位点的逐等位重算（AC/AF Number=A） | 默认双等位门剔除；显式允许多等位且需重算时**报错**（不再静默写错） | 多等位数据（indel 集）无法直接用 `--fill-tags` | 解码层已保留原始等位下标（`ScanOutput::alleles`）；`fill_tags` 按 Number=A 输出逐 ALT 向量，MAF 取第二常见等位（fill-tags 口径）；pixy 式多等位 π/dxy 推广已在累加器中 |
| `reed stats` 的样本缺失阈值 | 不提供（明确报错） | 用户需先 `reed filter` | 维持；02 的 `--mind` 同样只在 filter 中 |
| 报告模式的 `-t/-T` | 报错 | 可用性 | 01 F4-10 |
| `FMT/` 下标 `[0]/[*]`、通用正则 | 明确报 "not supported yet" | 少见 | 需求出现再做 |
| HWE mid-p 的 plink2 数值对拍 | `midp` 已实现 | 极端位点末位差异未量化 | 补对拍用例（plink2 `--hwe P midp`） |

---

## 4. gold/popgen 的缺口

现有：`run_pixy.sh`、`run_vcftools.sh`、`compare.py`、`compare_pbs.py`、`compare_reports.py`。
缺：easySFS、plink2 LD、scikit-allel f 统计、bcftools stats 的自动化比较脚本（现为手工）、**边界 fixture 生成器**（02 §7 的缺失阈值边界、半调用、`X·2N` 整数/非整数）。
建议把 fixture 生成器放在 `gold/popgen/make_fixtures.py`，输出固定种子的小 VCF，所有对拍脚本共用。
