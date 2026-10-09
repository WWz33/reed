# 05 统计与过滤语义：对齐状态与 W5/W6 的设计

> 原则不变：经典方法 + pixy 的缺失数据处理；参数按科学方法命名；软件只在测试中作参照物；**只做架构创新，不做算法创新**。
> 本文三部分：①已逐位对齐的台账；②尚未接入 CLI 的累加器及其**出口设计**（doc_v4 STATUS 指出"输出形状未固定，属功能设计"，本文把形状定下来）；③多等位重算（W6）。

---

## 1. 已对齐（有参照物、有数字）

| 统计量 / 行为 | 参照物 | 结果 | 来源 |
|---|---|---|---|
| π、dxy（含不变位点，默认） | pixy 2.2.4 | 全量 28,470 行，最大相对误差 4.8e-14 | doc_v3 R7 |
| π，`--variant-only` | vcftools `--window-pi` | 最大相对差 1.05e-5（vcftools 的 n/(n−1) 校正） | 同上 |
| Watterson θ、Tajima's D | pixy | 3.1e-14 / 8.9e-11 | 同上 |
| Fst wc84 / hudson（成对） | pixy | 2.9e-15 / 2.2e-15 | 同上 |
| PBS | scikit-allel | 逐列一致 | 同上 |
| `--freq/--hardy/--het/--missing` | plink2 | 逐行一致 | 同上 |
| `--depth` | vcftools | 逐行一致 | 同上 |
| `--summary` | bcftools stats | 九节逐行一致 | doc_v3 |
| 过滤：gold 规则 | bcftools（gold） | 1,768,759 位点，`check_parity` PASS，正文 md5 `7aedc59ac798…` | W0b / W3 / W4 重跑 |
| 过滤：`--geno/--mind` 边界 | plink2、vcftools | N=10/7、含/不含半调用、重复小数边界逐位点一致；子集 263,420 / 381,672 | W0b |
| 过滤：`-i/-e`、`FMT/`+`--set-gts .`、`-q/-Q` | bcftools | 位点集合 / GT md5 一致 | doc_v3 |
| BCF 渲染 / 输出 / 原生过滤输入 | bcftools | 3,315,604 条 sha256 一致；`-Ob` 正文 ≡ `-Oz`；BCF 输入正文 ≡ VCF 输入 | doc_v3 / W4 |
| 跨 reader | mmap vs 流式 | 10 个 stats 输出、CSI sha256、gold 正文全部相同 | W3 |

---

## 2. 有实现、无 CLI 出口（`#[doc(hidden)]`）

| 累加器 | 位置 | 参照物 | 出口（本文 §3） |
|---|---|---|---|
| 投影 SFS（`ProjectedSfsAcc`、`projection_weights`、`ProjectionCache`） | `reed-stats/src/sfs.rs` | easySFS + dadi 2.4.4（已装；`sfs.vcf --proj 6` 的输出已固定） | `--sfs --proj`（§3.1，**在途**：`StatKind::ProjectedSfs` 已在工作区） |
| 联合 SFS（`JointSfsAcc`） | `reed-stats/src/individual.rs` | easySFS 二维谱 | `--sfs --proj` 两群体（§3.1） |
| 合并（r 群体）WC84（`WcMultiAcc`） | `reed-stats/src/fst.rs` | vcftools `--weir-fst-pop` ×r（合并，weighted Fst）、hierfstat `wc()` | `--fst-global`（§3.2，**已裁决**） |
| 成对矩阵（`FstMatrixAcc`） | 同上 | pixy、plink2 `--fst`（默认全部对） | 不需要新参数：`--fst` 已是所有群体对；实现替换 |
| Patterson f2/f3/f4、D、f4-ratio + block jackknife | `reed-stats/src/{patterson,jackknife}.rs` | scikit-allel `admixture.py`（1.3.13 已装） | §3.3 |
| LD（`LdDecayAcc`、`LdPruner`、`pair_sums`） | `reed-stats/src/ld.rs` | plink2 `--r2-unphased`、`--indep-pairwise` | §3.4 |

接入顺序：SFS → `--fst-global` → f 统计 → LD。每项"逐位/逐行对拍通过后才去掉 `#[doc(hidden)]`"，并在 `gold/popgen` 增加 `run_<tool>.sh` 与 compare 脚本（fixture 已由 `make_fixtures.py` 交付：`sfs.vcf`、`fstats.vcf`、`ld.vcf`、`multiallelic.vcf`、`missing_rate.vcf`）。

---

## 3. 出口设计【规】

### 3.1 投影 SFS：`--sfs [folded] --proj N[,N…]`

**语义。** 对 `--pop` 中的每个群体（或 `--sfs` 时未给 `--pop` 则全体样本为一个群体），把每个位点上该群体的已调用染色体 `m`、ALT 计数 `alt` 按超几何分布投影到 `n` 条染色体：`P(k) = C(alt,k)·C(m−alt,n−k)/C(m,n)`，`k ∈ 0..=n`（Marth 2004；dadi `project`）。`m < n` 的位点丢弃并计数（写入 `.log`）。`--proj` 给一个值则所有群体同用，给 r 个值按 `--pop` 中群体首次出现的顺序对应，个数不符报错。
- **非分离质量**（k = 0 与 k = n）保留在数组里（在途改动已把 `bins` 扩成 `n+1`），因为 dadi/fsc2 文件就是整个数组；π/S 等比值仍只用 `1..n`。
- **祖先态**：默认 unfolded，以 REF 为祖先态（与 easySFS 的 `--unfolded` 相同；文档写明这是约定）；`--sfs folded` 把 `k` 与 `n−k` 合并到 `min(k, n−k)`。
- **不变位点**：all-sites 输入时单态位点（`alt = 0`，`m ≥ n`）全部质量落在 k = 0；`--variant-only` 时 k = 0 只含文件里出现的位点。两种模式都写入 dadi 的 mask 行：k = 0 与 k = n 置 1（easySFS 默认），用户自行决定是否 unmask。
- **两群体**：`--pop` 有两个群体时额外写二维谱（`JointSfsAcc`，`(nA+1)×(nB+1)`）；三个以上群体只写各自一维谱（多维谱不做）。
- **不是窗口统计**：genome-wide 一个谱，不进 `_*.txt` 窗口行模型。需要新的 genome-wide 请求路径：`StatsRequest.sfs: Option<SfsRequest>`，每 chunk 的 partial 与窗口 partial 并列上交，sink 在收尾时合并写出（这是 doc_v4 STATUS 指出的"套不进现有 StatKind/窗口行模型"的解法；在途改动把它塞进 `StatKind::ProjectedSfs` 作为窗口统计，**与本节冲突**，应改为 genome-wide 路径，否则 `--levels genome` 之外的层级会产生无意义的分窗谱）。

**输出。** `PREFIX_sfs/dadi/<pop>-<n>.sfs`（三行：`<n+1> unfolded|folded`、数值行、mask 行）；`PREFIX_sfs/fastsimcoal2/<pop>_MAFpop0.obs`（`1 observations. No. of demes and sample sizes are on next line` / `1\t<n>` / 数据行）；二维：`<A>-<B>-<nA>-<nB>.sfs` 与 `jointMAFpop1_0.obs`。数值默认按 easySFS `--dtype int` 四舍五入？**否**：写浮点（`%.6f`），加 `--sfs-int` 可选取整——投影后的谱本来就是实数，取整是 easySFS 的默认而非方法要求。

**对拍。** `gold/popgen/run_easysfs.sh`：fixture `sfs.vcf --proj 6` → dadi 行 `1 3 3 2 1 1 0`、mask `1 0 0 0 0 0 1`；全量上以 `clean.3pop.txt` 的三个群体各投影到 2·min(n_i) 与 easySFS 逐分箱一致（容差 1e-6 相对，取整前）。门：fixture 逐值相等；全量每群体最大相对差 ≤ 1e-6。

### 3.2 合并 Fst：`--fst-global`（已裁决 2026-10-09）

布尔开关，修饰 `--fst`。无 `--fst` 报错；`--fst hudson --fst-global` 报错 `hudson Fst has no r-population (global) form; use wc84`；两个群体时值等于成对 wc84，提示不报错。算法 = `WcMultiAcc`：r 群体 WC84 方差分量 a、b、c（`n_i = an_i // 2`，`h̄` 由 `n_het`），窗口值 `Σa/(Σa+Σb+Σc)`，任一群体 `an_i = 0`、`n̄ ≤ 1`、`n_c ≤ 0` 的位点跳过并计数。输出 `PREFIX_fst_global.txt`：`pops chromosome window_pos_1 window_pos_2 avg_wc_fst no_snps`，`pops` 为逗号连接的群体名。`StatKey` `fst.wc84.global`。

对拍：`gold/popgen/run_vcftools_fst.sh`：vcftools `--weir-fst-pop` ×3 的日志 weighted Fst（= Σa/Σ(a+b+c)）与 reed 的 genome 层级值相等；逐位点 `WEIR_AND_COCKERHAM_FST` 与 reed 单位点窗口值相等（无半调用、双等位位点；vcftools 的 `n_i` 用有基因型的个体数，等价）。门：fixture 逐位点相对差 ≤ 1e-9；全量 genome 值相对差 ≤ 1e-9。

### 3.3 f 统计：`--f2 A,B`、`--f3 C;A,B`、`--f4 A,B;C,D`、`--dstat A,B;C,D`、`--f4-ratio A,B;C,D;X`

Patterson et al. 2012 的定义，逐位点用等位频率（含缺失时按已调用等位算频率，群体在该位点无调用则该位点跳过）：
- `f2(A,B) = Σ (pA − pB)²`（无偏：减去 `pA(1−pA)/(nA−1) + pB(1−pB)/(nB−1)`，scikit-allel `patterson_f2`）；
- `f3(C;A,B) = Σ (pC − pA)(pC − pB)`（减 `pC(1−pC)/(nC−1)`）；
- `f4(A,B;C,D) = Σ (pA − pB)(pC − pD)`；`D = Σ num / Σ den`，`num = (pA−pB)(pC−pD)`、`den = (pA+pB−2pApB)(pC+pD−2pCpD)`；
- `f4-ratio = f4(A,O;X,C) / f4(A,O;B,C)`（admixtools 的 α 估计，群体顺序写在参数里）。
分号分隔角色组，逗号分隔组内群体；可重复。标准误：block jackknife（`BlockJackknife`），块按 `--block-size`（默认 `5M` bp）切，`Z = est / SE`。

输出 `PREFIX.fstats`：`stat pops estimate se z n_blocks n_snps`，一行一个检验。对拍：`gold/popgen/compare_fstats.py` 对 scikit-allel `average_patterson_f3/d`（`blen` 对应块大小）与 `patterson_f2`；门：fixture `fstats.vcf` 逐值 1e-12，全量三群体 1e-9。

### 3.4 LD：`--ld-decay [MAX_BP]`、`--indep-pairwise W S R2`

- `--ld-decay`：对距离 ≤ MAX_BP（默认 `100k`）的位点对计算 r²（unphased，plink2 `--r2-unphased` 的基因型相关），按距离分箱（默认 1 kb）输出 `PREFIX.ld.decay`：`bin_start bin_end n_pairs mean_r2`；`--ld-pairs` 另写逐对 `PREFIX.ld`（plink2 `.vcor` 列：`CHROM_A POS_A ID_A CHROM_B POS_B ID_B UNPHASED_R2`）。
- `--indep-pairwise W S R2`：plink2 的变体修剪，W 窗口（位点数或 `kb`）、S 步长、R2 阈值；输出 `PREFIX.prune.in`、`PREFIX.prune.out`。
- **架构**：距离受限的逐对统计可流式：worker 对 chunk 内的对直接算，跨 chunk 的对由"尾窗携带"（每 chunk 末尾 MAX_BP 内的位点随 partial 上交，reducer 补算跨界对）；`--indep-pairwise` 有顺序依赖，按 contig 串行、contig 之间并行（`LdPruner`）。
- 对拍：`run_plink2_ld.sh` 对 fixture `ld.vcf`：r² 逐对 1e-9；prune 集合相同。

---

## 4. 多等位逐等位重算（W6）

现状：默认双等位门剔除；显式允许多等位且需 AC/MAF/HWE 重算时**报错**（不静默写错）。设计：
- `--fill-tags`：AC/AF 按 `Number=A` 逐 ALT 输出（`ScanOutput::alleles` 已保留原始等位下标）；MAF = 第二常见等位频率（`bcftools +fill-tags` 口径）；HWE 只在双等位位点算，多等位写 `.`；NS、AN、F_MISSING 与双等位相同。
- `--maf/--mac`：多等位时次等位 = 第二常见等位（同上口径）；`--geno` 不变。
- stats：π/dxy 的多等位推广已在累加器（`comps − Σ_a C(k_a,2)`）；θ 的突变数 `Σ(等位数−1)`。
- 对拍：`multiallelic.vcf` 上 `--fill-tags AC,AN,AF,MAF,HWE,NS` 与 `bcftools +fill-tags` 逐字节一致；全量 indel 集抽样 1,000 条逐字节一致。门：默认双等位门行为不变（gold md5 不变）。

---

## 5. gold/popgen 的缺口（接入时补齐）

现有：`run_pixy.sh`、`run_vcftools.sh`、`compare.py`、`compare_pbs.py`、`compare_reports.py`、`make_fixtures.py`。
缺：`run_easysfs.sh` + `compare_sfs.py`、`run_vcftools_fst.sh`、`compare_fstats.py`、`run_plink2_ld.sh` + `compare_ld.py`、`run_filltags.sh`（W6）。每个脚本的判据写在 §3/§4 的"门"里；通过后把结果数字登记到当轮 STATUS。
