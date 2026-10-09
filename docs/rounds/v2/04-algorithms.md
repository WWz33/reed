# 04 算法：经典方法 + pixy 的缺失数据处理

> 用户裁决（2026-10-06）：
> - reed **不做算法创新，只做架构创新**；
> - 统计参数**只报告科学方法（谁发明的方法），不报告软件**；
> - 缺失数据的处理以 **pixy（Korunes & Samuk 2021）为默认**。
>
> 本文给出：①各软件在经典统计量上到底哪里一致、哪里不一致（源码证据）；②reed 采用的唯一定义；
> ③需要撤下的 reed 自创语义；④对拍方案（软件只在测试中出现）。

---

## 1. 撤下的 reed 自创语义

| 自创语义 | 位置 | 处理 |
|---|---|---|
| 清单门：按 `SiteCompleteness`/`Ascertainment` 拒绝或标注统计量 | DESIGN §4、v1 §5.8/§11.7、`reed-core/src/manifest.rs`、v0.2 各估计器 `availability()` | **删除**。改为 pixy 的不变位点检查 + `--variant-only` 声明（02-cli.md §5.3）。`##reed_ascertainment` 不再写出 |
| SFS/θ/Tajima 的 complete-case | v0.2 `spectrum.rs`、M4 `reed-stats/src/{sfs,theta,tajima}.rs` | **删除**，θ/Tajima 用 §3 的定义；SFS 用投影（投影到全样本量时等价于 complete-case） |
| 投影作为 θ/Tajima 的默认 | v1 §9.4 | **撤回** |
| π/dxy 的 `PerPairComparable`/`AccessibleBp` 两套估计器 | v0.2、M4 | 合并为一个定义：pixy 的比值之和；`--variant-only` 只改变分母中不变位点的来源 |
| PBS 的三群体交集位点集合 | v0.2 `pbs.rs` | **删除**，三对 Hudson Fst 各自独立计算（Yi et al. 2010 的原始定义） |
| Nei86 Fst、Hudson 有限样本校正 | v1 §9.2 | 不实现 |
| 软件兼容模式 `--compat` | doc_v2 上一版 | **撤回**（用户裁决：参数只写方法，不写软件） |
| reed 自定义窗口网格 | M4 | 改为 §4 |
| `--mind` 在原始基因型上计算、表达式偏离 bcftools、`--min-dp` 的缺失规则 | doc_v2 早先版本 | 已在 02-cli.md 中改为 plink2/bcftools 语义 |

---

## 2. 各软件在经典统计量上的异同（源码核实）

结论：**公式完全相同**。差异只出现在下表三个维度上，其中缺失数据是主因。
无缺失数据、且分母约定相同时，vcftools、scikit-allel、pixy 的结果逐位一致。

| 维度 | vcftools 0.1.18 | scikit-allel 1.3 | pixy 2.2.4 |
|---|---|---|---|
| π 的分母 | 变异位点用实际调用的比较对数；**窗口中其余位置视为单态且 2N 条染色体全部调用**（`variant_file_output.cpp:4252-4272`） | 每位点平均差异之和 / 窗口碱基数（`diversity.py:203, 397`） | Σ差异 / Σ比较，**只统计文件中出现的位点**，因此要求不变位点在文件中（`calc.py:208`） |
| Tajima's D 的 n | 固定为 2 × 个体数（`variant_file_output.cpp:3954`） | 整条 contig 的最大调用等位数（`diversity.py:970`） | 窗口平均调用数，方差只算一次（`calc.py:636-740`，ThetaRecov） |
| Watterson θ 的 n | 不提供 | 整条 contig 的最大调用等位数（`diversity.py:674`） | 逐位点调用数的调和数（`calc.py:549`） |
| WC84 的样本量 | 基因型计数 N_hom + 0.5·N_het（`variant_file_output.cpp:3564-3620`） | `an // 2`（`fst.py:136-241`） | 调用 allel，同左 |
| Hudson Fst 的位点纳入 | 不提供 | 分母要求两侧有调用，分子要求两侧 n ≥ 2（`fst.py:247`） | 分子无定义时该位点从分子分母同时剔除（`calc.py:332-376`） |
| 窗口网格 | `[s·step+1, s·step+size]`，空窗不输出 | 从首个变异位点开始（`window.py:96`） | 从 1 开始，按窗口大小不重叠 |

---

## 3. reed 的定义（唯一，参数中只出现方法名）

记号：某窗口内文件中出现的位点集合为 S；位点 i 上群体 A 的已调用等位数为 n_i，其中 ALT 计数为 k_i（多等位时按等位分别计数）。

| 统计量 | 定义 | 文献 |
|---|---|---|
| π | `Σ_i diffs_i / Σ_i comps_i`，其中 `comps_i = n_i(n_i−1)/2`，`diffs_i` 为不同等位的配对数（多等位推广：`comps_i − Σ_a C(k_ia, 2)`）。`--variant-only` 时分母再加上 `(L − |S|) · C(2N, 2)`，L 为窗口长度（或 BED 窗口的长度），N 为群体个体数 | Nei & Li 1979；缺失数据：Korunes & Samuk 2021 |
| dxy | `Σ_i diffs^{AB}_i / Σ_i n^A_i n^B_i`；`--variant-only` 时分母加上 `(L − |S|) · 2N_A · 2N_B` | Nei 1987；Korunes & Samuk 2021 |
| Watterson θ | `θ_raw = Σ_{i∈S, 分离} (k_i^{等位数}−1) / a(n_i)`，`a(n) = Σ_{j=1}^{n−1} 1/j`；窗口值 = `θ_raw / 位点数`，`--variant-only` 时除以 L | Watterson 1975；逐位点 n：pixy |
| Tajima's D | `(θ_π − θ_W) / sqrt(e1·S + e2·S(S−1))`，常数 a1、a2、b1、b2、c1、c2、e1、e2 按 Tajima 1989；n 取窗口内（有调用的位点上）平均调用数并四舍五入，方差只计算一次；S 为突变数 Σ(等位数−1) | Tajima 1989；方差的 n 取法：pixy 2.x（ThetaRecov） |
| Fst `wc84` | Weir & Cockerham 1984 方差分量 a、b、c，每群体样本量 `n = an // 2`，窗口值 Σa / Σ(a+b+c)；r 个群体时按原文一般式 | Weir & Cockerham 1984 |
| Fst `hudson` | 每位点 `num = H_between − (h_A + h_B)/2`、`den = H_between`，h 为无偏期望杂合度；任一群体 n < 2 时该位点分子分母同时剔除；窗口值 Σnum / Σden | Hudson et al. 1992；Bhatia et al. 2013 |
| PBS | 三对 Hudson Fst（各自独立的位点集合），`T = −ln(1 − Fst)`，Fst 截断到 `[0, 0.99999]`（截断约定取自 scikit-allel，原文未规定）；`PBS = (T_AB + T_AC − T_BC)/2`；`--pbs-normed` 时除以 `1 + (T_AB + T_AC + T_BC)/2` | Yi et al. 2010；Malaspinas et al. 2016 |
| SFS | 调用数 `m ≥ N` 的位点按超几何分布投影到 N 个单倍型；`m < N` 的位点丢弃 | Marth et al. 2004；Gutenkunst et al. 2009 |
| HWE | 精确检验；`midp` 时为 mid-p | Wigginton et al. 2005；Graffelman & Moreno 2013 |
| 个体 F | `(O_hom − E_hom) / (N − E_hom)` | 矩估计（plink `--het` 的定义） |

说明：

- 以上就是 pixy 2.2.4 对 π、dxy、θ、Tajima's D、两种 Fst 的定义；reed 只额外提供 `--variant-only` 的分母补全。
  这个补全与 vcftools `--window-pi` 的做法相同（窗口内其余位置视为单态且全部调用），因此 `--variant-only` 下的 π 可以直接用 vcftools 对拍。
- Fst、PBS、Tajima's D 是比值统计量，单态位点对分子分母的贡献都是 0，不受 `--variant-only` 影响。
- pixy 只有不重叠窗口；reed 的 `--step` 产生的重叠窗口，等价于把这些窗口写成 BED 交给 pixy 的 `--bed_file`。

---

## 4. 窗口网格（唯一）

- 默认网格从 1 开始：`[1, w]`、`[1+s, w+s]`、…，s 默认等于 w（与 pixy 一致；s < w 时与 vcftools 的步长网格一致）。
- 末窗：截断到 contig 长度；header 缺少长度时报错（可用 `--contig-lengths` 补充），**不得编造**。
- 空窗：输出，数值写 `NA`（与 pixy 一致）。
- 分箱累加器（`gcd(w, s)` 分箱）作为架构层的实现方式保留，分箱边界从 1 开始对齐网格。

---

## 5. 对拍方案（软件只出现在测试中）

| 对象 | 参照物 | 判据 |
|---|---|---|
| π、dxy、θ、Tajima's D、Fst（wc84/hudson），含缺失数据、含不变位点 | pixy 2.2.4 | 逐列一致（浮点相对误差 1e-9） |
| `--variant-only` 的 π | vcftools `--window-pi`（逐群体 `--keep`） | 按 vcftools 的打印精度逐行一致 |
| 无缺失数据的 fixture | vcftools、scikit-allel、pixy 三家同时对拍 | 三家与 reed 互相一致，用来证明"经典统计量在各软件间统一"这一前提 |
| PBS | scikit-allel `pbs`（`normed=False` 对应原始定义，`normed=True` 对应 `--pbs-normed`），Fst 部分用 pixy 的 Hudson 位点纳入规则构造 fixture | 逐列一致 |
| SFS | easySFS（同样的 `--proj`） | 逐分箱一致 |
| freq、hardy、het、missing | plink2 | 按打印精度逐行一致 |
| depth | vcftools | 逐行一致 |
| summary | bcftools stats | 逐行一致 |

`gold/popgen/` 的布局：`fixtures/`（合成：含/不含不变位点 × 有/无缺失 × 不等样本量 × 多等位 × 3 群体）、`run_<tool>.sh`、`compare.py`。
全量数据上至少对拍一次 vcftools 的 `--window-pi`、`--TajimaD`（在无缺失数据的子集上）与 plink2 的 `--freq/--hardy/--het/--missing`。

`reed-core-v02-oracle` 中与上述定义相同的部分（Hudson/WC84 的分量计算、等位计数）继续作为差分基准；complete-case 与清单门相关的 oracle 用例随语义撤下而删除。

---

## 6. 输入契约（I/O 正确性，继续保留）

| 契约 | 处理 |
|---|---|
| 同一 contig 内位置升序 | 违反即报错（worker 检查 chunk 内部，reducer 检查 chunk 之间） |
| contig 不回访 | 报错 |
| 窗口需要 contig 长度而 header 缺失 | 报错（`--contig-lengths` 补充），不编造（01-review D7） |
| 相同 POS 的多条记录 | 每条作为独立位点 |
| 多等位位点 | 按 §3 的多等位推广逐等位计数（pixy 的做法），不折叠 |

---

## 7. 过滤侧的对齐要点

- `--fill-tags`：复现 `bcftools +fill-tags`（多等位 AC/AF 按 Number=A，MAF 取第二常见等位，退化位点 HWE = 1）。
  现实现把多等位的所有 ALT 合并计数，需要修正。
- `-i/-e`：bcftools 表达式语义。
- `--min-dp/--max-dp/--min-gq`、`--mind/--geno/--maf/--mac/--hwe`：plink2 语义（缺失 DP/GQ 不掩码；`--mind` 在掩码之后计算）。
- `--hwe` 的检验方法写明 Wigginton 2005 / mid-p；与 plink2 的数值对拍（plink2 在极端位点上用 log 域计算，末位可能不同，对拍时记录容差）。
- 求值顺序只影响漏斗报告，不影响结果（02-cli.md §6）。
