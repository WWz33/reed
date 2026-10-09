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
| 多群体 WC84（`WcMultiAcc`、`FstMatrixAcc`） | `reed-stats/src/fst.rs` | vcftools 多个 `--weir-fst-pop` | `--fst-pops A,B,C` |
| LD（`LdDecayAcc`、`LdPruner`、`pair_sums`） | `reed-stats/src/ld.rs` | plink2 `--r2-unphased`、`--indep-pairwise`；vcftools `--geno-r2` | `--ld-decay`、`--indep-pairwise`；需要平面存储（整 contig） |
| Patterson f2/f3/f4、D、f4-ratio（`F2Acc/F3Acc/F4Acc`、`f4_ratio`） | `reed-stats/src/patterson.rs` | scikit-allel `admixture.py` | `--f2/--f3/--f4/--dstat/--f4-ratio`，block jackknife 输出 |
| Block jackknife（`BlockJackknife`） | `reed-stats/src/jackknife.rs` | scikit-allel `average_patterson_d` | 随 f 统计接入 |

接入顺序建议：SFS（需求最常见、参照物最简单）→ `--fst-pops` → f 统计 → LD。每项接入时把 `#[doc(hidden)]` 去掉，并在 `gold/popgen` 增加对应的 `run_<tool>.sh` 与 compare 脚本。

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
