# 01 复核清单：sgkit API 覆盖、科学符合性与 API 工程质量

> 本文件分两部分：**§1–§5** 是 sgkit 覆盖与科学符合性（V30-01…09）；**§6–§13** 是
> reed-rs 作为 API 库的工程质量（V30-Q01…Q36）。
>
> 立项 2026-11-11。对象：sgkit 公开 API（`/sri/home/sri2025201067/reed/sgkit/docs/api.rst`）
> 与 reed 当前工作树的对照。
>
> **方法**：6 个只读子代理按域审查（基础统计/QC、群体遗传、LD、亲缘/PCA、关联/系谱、
> IO/窗口），主代理逐条复核关键结论并亲自复现三个缺陷。
>
> **证据标记**：
> `repro` = 用当前源码构建的二进制实测复现；`code` = 确定性代码路径可证明；
> `formula` = 公式层核验（未跑 sgkit 数值夹具）；`doc` = 文档/契约层；`N/A` = 架构上不适用。
>
> **基线**：HEAD `e920663` + 未提交改动；复现二进制 `reed 0.4.0 (v0.3.0-177-ge920663-dirty
> 2026-11-10)`，构建于 `/tmp/reed-target/release/reed`（`.cargo/config.toml` 把 target-dir
> 重定向到本地盘）。
>
> **重要边界**：仓库内**没有 sgkit 数值夹具**。因此"与 sgkit 等价"绝大多数是**公式层**
> 结论；现有对拍基准是 pixy / scikit-allel / plink2 / vcftools / bcftools / ADMIXTOOLS。
> 任何"与 sgkit 数值相同"的说法在本轮都**未经数值验证**，只在 §2.1 明示处有对拍。

---

## §1 覆盖总表

口径：`docs/api.rst` 的 autosummary 共 **54** 项；`sgkit/__init__.py` 另导出 `pca` 与
`regenie_loco_regression` 未被 api.rst 收录，合计 **56** 项。`create_genotype_call_dataset` /
`create_genotype_dosage_dataset` 是内存数据集构造器，与流式 CLI 无可比面，不纳入 56。

| 判定 | 数量 |
|---|---|
| 完整或等价 | **12** |
| 部分 / 换算法实现 | **13** |
| 未实现 | **27** |
| 不适用 | **4** |

### 1.1 IO / Dataset（3）

| sgkit API | reed 对应 | 状态 | 证据 |
|---|---|---|---|
| `unpack_variables` | — | N/A | sgkit 仅还原旧 packed BGEN（`sgkit/io/bgen/utils.py:7-84`）；reed 无 BGEN/dosage 层 |
| `load_dataset` | — | N/A | reed 无 Zarr/xarray Dataset（全仓无 zarr 依赖） |
| `save_dataset` | — | N/A | 输出是 VCF/BCF/文本矩阵，非 Zarr |

### 1.2 基础统计（8）

| sgkit API | reed 对应 | 状态 | 证据 |
|---|---|---|---|
| `call_allele_frequencies` | — | 未实现 | 无 per-sample×allele 输出；`reed-core/src/count/site_counts.rs` 只有群体级计数 |
| `cohort_allele_frequencies` | `stats --pop --freq`（`.afreq.pop`） | 完整（二态） | `reed/src/stats/writers/reports.rs:1430` |
| `count_call_alleles` | — | 未实现 | 同 `call_allele_frequencies` |
| `count_cohort_alleles` | `stats --pop --freq` 路径 | 完整（二态） | `reports.rs:1425-1438` |
| `count_variant_alleles` | `filter --fill-tags AC,AN`；`stats --freq` | 完整（二态） | `reed-filter/src/eval/edit_plan.rs:39-47`；`reports.rs:843-857` |
| `count_variant_genotypes` | `stats --hardy`（hom_ref/het/hom_alt） | 完整（二态二倍体） | `reports.rs:859-897` |
| `variant_stats` | `.afreq`+`.vmiss`+`.hardy`+`.summary` 组合 | 部分（无单一入口；无 `n_non_ref` 列；半调用计入口径不同） | `reports.rs:843-909` |
| `sample_stats` | `--missing`(.smiss) + `--summary` PSC | 部分 | `reports.rs:1159`；`reed-stats/src/summary/report.rs:165-192` |

### 1.3 群体遗传（8）

| sgkit API | reed 对应 | 状态 | 证据 |
|---|---|---|---|
| `divergence` | `popg --dxy`；`DxyAcc` | 完整（口径变体） | `reed-stats/src/dxy.rs:11` |
| `diversity` | `popg --pi`；`PiAcc` | 完整（口径变体） | `reed-stats/src/pi.rs:10` |
| `Fst` | `--fst hudson\|wc84`、`--fst-global` | 部分（缺 Nei；Hudson+WC84 完整） | `reed-stats/src/fst.rs:14-19` |
| `Garud_H` | — | 未实现 | 全仓无 `Garud` 命中 |
| `individual_heterozygosity` | — | 未实现 | 无 per-call 输出；`reed-stats/src/het.rs` 是 plink `--het` 矩估计（不同统计量） |
| `observed_heterozygosity` | `stats --hardy` 的 Ho 列 | 部分（全样本逐位点，非按群体/窗口） | `reports.rs:859-875` |
| `pbs` | `--pbs`、`--pbs-normed`；`PbsAcc` | 完整 | `reed-stats/src/pbs.rs:15-20` |
| `Tajimas_D` | `--tajima-d`；`TajimaAcc` | 完整（n 用窗口均值，pixy 口径） | `reed-stats/src/tajima.rs:21` |

### 1.4 亲缘（10）

| sgkit API | reed 对应 | 状态 | 证据 |
|---|---|---|---|
| `Weir_Goudet_beta` | — | 未实现 | 无 β；`phi`/`ajk` 均非"相对平均 kinship 中心化" |
| `genomic_relationship` | `distance --type grm\|cov` | 部分（GCTA/plink per-marker 加权，非 sgkit 单分母 VanRaden；无 E-J 收缩、无 ancestral AF） | `reed-distance/src/matrices.rs:287-350` |
| `pc_relate` | — | 未实现 | 无 PC 回归；`king`/`phi` 不等价（结构群体下会被亲缘/结构混淆） |
| `pedigree_contribution` | — | 未实现 | 全仓无 pedigree |
| `pedigree_inbreeding` | — | 未实现 | 同上（`--het` 是基因组矩估计，输入与概念都不同） |
| `pedigree_inverse_kinship` | — | 未实现 | 同上 |
| `pedigree_kinship` | — | 未实现 | 同上 |
| `hybrid_inverse_relationship` | — | 未实现 | 同上 |
| `hybrid_relationship` | — | 未实现 | 同上 |
| `identity_by_state` | `distance --type ibs\|1-ibs` | 部分（定义不同：plink IBS het-het=1，sgkit 概率=0.5） | `reed-distance/src/matrices.rs:143-270` |

### 1.5 LD（2）

| sgkit API | reed 对应 | 状态 | 证据 |
|---|---|---|---|
| `ld_matrix` | `popg --ld-decay --ld-pairs` | 部分（无 r² 门限/`cmp`/对角线；窗口只能按 bp） | `reed-stats/src/matrix.rs:284-357` |
| `ld_prune` | `filter --indep-pairwise` | 部分（plink2 分窗贪心 + MAF tie-break，非 sgkit 全局 MIS；阈值 `>` vs `>=`） | `reed-stats/src/ld.rs:452-590` |

### 1.6 关联（4）

| sgkit API | reed 对应 | 状态 | 证据 |
|---|---|---|---|
| `gwas_linear_regression` | — | 未实现 | `rg -i "gwas\|linear_regression\|phenotype"` 0 命中 |
| `regenie` | — | 未实现 | 0 |
| `hardy_weinberg_test` | `stats --hardy` P 列；`filter --hwe`；`--fill-tags HWE` | 完整（退化位点返回 1.0 哨兵，非 sgkit 的 NaN） | `reed-core/src/math/hwe.rs:67-131`；`reports.rs:889-895` |
| `genee` | — | 未实现 | 0 |

### 1.7 Miscellaneous（2）

| sgkit API | reed 对应 | 状态 | 证据 |
|---|---|---|---|
| `maximal_independent_set` | 无独立 API（并入 `--indep-pairwise` 的窗口内贪心） | 部分 | `reed-stats/src/ld.rs:452-570` |
| `pairwise_distance` | `reed_reduce::graph::pairwise_distances`（仅欧氏） | 部分（无 correlation；无逐对缺失屏蔽） | `reed-reduce/src/graph.rs:29-40` |

### 1.8 Utilities（17）

| sgkit API | reed 对应 | 状态 | 证据 |
|---|---|---|---|
| `convert_call_to_index` | `convert --012` | 部分（二态二倍体；单倍体→−1、多等位跳过） | `reed/src/convert/job.rs:1221` |
| `convert_probability_to_call` | — | 未实现 | 无 GP 解码（全仓仅在 BCF 往返测试串出现 `GP` 标签） |
| `display_genotypes` | — | 未实现 | 无 dump/display 子命令 |
| `display_pedigree` | — | N/A | 无 pedigree 数据模型 |
| `filter_partial_calls` | `HalfCallPolicy::TreatAsMissing`（库内能力，CLI 未暴露） | 部分 | `reed-core/src/genotype/gt_code.rs:231-259` |
| `infer_call_ploidy` | — | 未实现 | 运行级 `PloidyPolicy`，无 per-call 推断 |
| `infer_sample_ploidy` | — | 未实现 | 同上 |
| `infer_variant_ploidy` | — | 未实现 | 同上 |
| `invert_relationship_matrix` | — | 未实现 | 无公开矩阵求逆 |
| `parent_indices` | — | 未实现 | 0 |
| `pedigree_sel` | — | 未实现 | 0 |
| `simulate_genedrop` | — | 未实现 | 0 |
| `simulate_genotype_call_dataset` | — | 未实现 | 仅 `crates/reed-reduce/examples` 内的私有模拟器 |
| `window_by_genome` | `popg --levels genome` | 完整 | `reed/src/stats/sink.rs:268-272` |
| `window_by_interval` | `--bed-file`（BED3） | 完整 | `reed-cli/src/translate/popg.rs:771-820` |
| `window_by_position` | `--window/--step` | 部分（无 offset；需 contig 长度） | `reed-core/src/grid.rs:118-146` |
| `window_by_variant` | — | 未实现 | 结构存在（`grid.rs:29`）但驱动明确拒绝：`reed/src/stats/driver.rs:233-235` |

### 1.9 `__init__` 额外导出（2）

| sgkit API | reed 对应 | 状态 | 证据 |
|---|---|---|---|
| `pca` | `reed pca`（binom/binom2、dense/stream、缺失 meanimpute/pairwise） | 完整（`--stand binom` 同 Patterson；`.eigenval` 标度不同，见 §2.2-9） | `reed/src/pca/mod.rs`；`reed-reduce/src/pca.rs` |
| `regenie_loco_regression` | — | 未实现 | 0 |

---

## §2 科学符合性

### 2.1 站得住（公式层等价 + 有对拍）`formula`

| 项 | 结论 | 依据 |
|---|---|---|
| π / dxy | 逐位点公式与 sgkit 代数等价（`n_diff/n_pairs`、`comps=C(n,2)`） | `reed-stats/src/pi.rs:122-125`；对拍 pixy 1e-9（`gold/popgen/README.md`） |
| θ_W | `a1(n)=H(n−1)` 与 sgkit `_harmonic_sum(n,1)` 逐项相同 | `reed-core/src/math/harmonic.rs:79-82` |
| Tajima's D | 常数块 `e1,e2,b1,b2,c1,c2` 与 sgkit 逐项相同 | `reed-stats/src/tajima.rs`；对拍 pixy |
| Fst·Hudson | ratio-of-sums（allel/pixy 口径），两群体均要求 `n≥2` | `reed-stats/src/fst.rs:209-232`；对拍 vcftools/allel |
| Fst·WC84 | `Σa/(Σa+Σb+Σc)` 标准式（sgkit 根本没有 WC84） | `fst.rs:299-306` |
| PBS | `PbsN1` 与 sgkit `_pbs_cohorts` 归一化式逐字一致 | `pbs.rs`；对拍 scikit-allel |
| HWE | Wigginton 2005 精确检验，log 域递推；与 sgkit 数值等价 | `reed-core/src/math/hwe.rs:67-119`；oracle 自 sgkit 可译移植 |
| LD r² | 整数精确 Pearson，与 Rogers–Huff 代数等价；成对剔除 | `reed-stats/src/ld.rs:98-115`；对拍 plink2 `--r2-unphased` rtol 2e-6 |
| f2/f3/f4/D/f4-ratio | Patterson 2012 无偏校正代数与 scikit-allel 一致；块刀切为 ADMIXTOOLS 语义 | `reed-stats/src/patterson.rs:21-87`；`jackknife.rs`；qpDstat 对拍 2e-12（`gold/popgen/run_admixtools_jackknife.sh`） |
| PCA（binom） | `W=(x−2p)/√(p(1−p))` 与 Patterson 2006 一致；PVE 与 sgkit `ratio` 等价 | `reed-reduce/src/pca.rs` |
| VCF/BCF 读写 | 与 bcftools fixture 逐字节一致、与 noodles 差分一致 | `reed-io/tests/bcf_render.rs`、`bcf_vs_noodles.rs` |

### 2.2 口径变体（科学上成立，但需在文档写清，否则会被当缺陷） `formula`/`doc`

1. **窗口聚合**：sgkit 是"逐位点值求和"，reed 是"比值之和"（pixy）。缺失率异质时不同。
   `pi.rs:122-125`；`docsite/src/content/docs/popg.md:169`（已写）。
2. **Tajima's D 的 n**：sgkit 用 `max(总等位计数)`，reed 用窗口均值（round half-to-even）。
   恒定 n、完整数据下重合。`reed-stats/src/tajima.rs`。
3. **θ_W 按每位点可用 n 分箱**（pixy/Korunes & Samuk），sgkit 用单一 n。
4. **f-stat 的 `estimate` 是 Σ per-site，不是位点均值**：`reed/src/stats/request.rs:700-733`
   注释明确对齐 ADMIXTOOLS 的 sum 尺度。z 与 D/f4-ratio 不受影响，但点估计与 SE 比文献
   大 ~n_snps 倍；`popg.md:236-253` 只给逐位点公式，**未显式说明"estimate=Σ"**。
5. **`f3*` 过度承诺**：`popg.md:242` 称"normalized variant also reports ΣT/ΣB"，但 CLI 无
   归一化选项（`value_normalized` 只存在于未接线的 `F3Acc`）。
6. **`--f4-ratio` 角色序**：`reed-cli/src/translate/popg.rs:671-681` 映射
   `A=parts[0], O=parts[1], C=parts[2], B_ref=parts[3], X=parts[4]`，即
   `f4(A,O;X,C)/f4(A,O;B,C)`；CLI/文档只写 `A,B;C,D;X`，按字面把 D 当 outgroup 会**静默**得到另一个 alpha。
7. **`observed_heterozygosity` 非按群体**：`--hardy` 是全样本逐位点（`reports.rs:859-875`），
   `--pop` 只分层 `.afreq`/`.vmiss`，无法用 CLI 复现 sgkit 的 per-cohort Ho。
8. **半调用口径**：`.vmiss` 的 `MISSING_CT` 只统计完全缺失（`reports.rs`），半调用算"已调用"，
   与 sgkit `variant_stats` 的 `n_called`（只计完整调用）不同 → call_rate 偏乐观。
9. **PCA 输出标度**：默认 `binom2` 的特征值与 sgkit `explained_variance` 差一个由 m、n
   决定的常数（PVE 不受影响）；`.eigenvec` 是单位 U，sgkit 默认给 projection `U√λ`；
   `Pca::loadings` 的缩放与自身注释不符（`reed-reduce/src/pca.rs:58-77`，CLI 未暴露）。
10. **`genomic_relationship`**：reed 是 GCTA/plink2 的 per-marker 加权 GRM，
    不是 sgkit 的单分母 VanRaden（仅当各标记 p 相等时两者等价）；无 Endelman-Jannink 收缩、
    无 ancestral-frequency 参数。
11. **`identity_by_state`**：reed `--type ibs` 是 plink 的 IBS（het-het=1），
    sgkit 是"两随机等位基因 identical-in-state 的概率"（het-het=0.5）；缺失时 reed 用全 m + HWE 权重，
    sgkit 用成对 called。**不能互相代入**。
12. **`convert --012`**：单倍体→−1（sgkit 给等位索引）；多等位跳过（sgkit 有 comb 索引）。
13. **MSA 全缺失→REF**：`reed/src/convert/job.rs:1246-1290`（vcf2msa step 4 原语义）。
    文档已写，但科学上等于把缺失伪造成参考纯合，做系统发育前必须知道。
14. **`window_by_position` 无 offset**；**`window_by_variant` 不可用**（见 §1.8）。

### 2.3 真缺陷

#### V30-01（P0，`repro`+`code`）PI_HAT 的期望 IBS 系数错位

`reed-distance/src/matrices.rs:570-592` 写入的系数是：

```
e00 = 2p²q²   e01 = 0   e02 = 0   e11 = 2pq(p²+q²)   e12 = 2pq
```

plink 1.9 的无限样本极限（按 `plink-ng/1.9/plink_calc.c:4863-4867` 逐项退化，
`dxx1/dxx2/dyy1/dyy2→1`、`num_allelesf2/f3→1`）应为：

```
e00 = 2p²q²   e01 = 4pq(p²+q²)   e02 = (p²+q²)²   e11 = 2pq   e12 = p²+q²
```

即 `e11` 拿到了 `e01` 的公式、`e12` 拿到了 `e11` 的公式，`e01/e02` 被清零。
求解器 `z0=ibs0/(e00·n)`、`z1=(ibs1−z0·e01·n)/(e11·n)`、`z2=oo/n−z0·e02−z1·e12`
（`matrices.rs:519-521`）与 plink 同构，系数一错结果全错。

**复现**（`reed distance subvcf.vcf.gz --type pihat`，348 样本 / 3000 位点 / 60378 对）：

```
Z0      Z1      Z2      PI_HAT   对数
0.0000  1.0000  0.0000  0.5000   60346
...（其余 32 对为近重复样本，值各异）
```

即**绝大多数无关个体对被报成 PI_HAT=0.5**。`docs/rounds/v17/STATUS.md:88-107` 把它归因于
"plink 用 PLINK 1.07 的 Nei 偏差修正"——该归因**不充分**：即使忽略有限样本修正，
系数排列本身也与 plink 不符。

**注**：`--type ibs`/`1-ibs`/`king`/`phi`/`grm` 不受影响（各自独立实现，且 v17 有对拍）。

#### V30-02（P1，`repro`+`code`）`filter -I` 文档与行为相反

- 文档页 `docsite/src/content/docs/filter.md:155-156`：
  `-I, --no-update  Leave INFO/AC and INFO/AN alone. By default they are recomputed whenever the sample set changes`。
  `--help`（`reed-cli/src/cli/mod.rs:365-366`）措辞更弱：`Do not update INFO/AC,INFO/AN when samples are dropped`。
- 实现：默认陈旧清单含 AC/AN（`reed-filter/src/spec/recompute.rs:40-54`），
  `-I` 只关闭**重算**（`reed/src/filter/mod.rs:520-538`），**不把 AC/AN 移出删除清单**。
  因此 `filter.md` 的"Leave alone"为**假**；`--help` 的"do not update"也**不完整**
  （实际是删除，且无样本变化时也删）。

**复现**（先 `--fill-tags AC,AN` 造出带 AC 的输入，原值 AC=331，再随机只留 1 个样本）：

| 参数 | 输出带 AC 的记录 | 首条 AC |
|---|---|---|
| `-s S`（默认重算） | 3000 | `AC=1`（正确：1 样本） |
| `-s S -I` | **0** | —（被删） |
| `-s S -I -x ''` | 3000 | `AC=331`（原值保留） |
| `-I`（无样本变化） | **0** | —（被删） |

即真正"保留 AC/AN"必须写 `-I -x ''`；且默认 `-x` 在**没有任何样本变化**时也删除这 8 个 tag，
比 bcftools 的默认（不删）破坏性大。bcftools `view -I` 的语义是保留原值（不删）。

#### V30-03（P1，`code`+`doc`）多等位三条路径不一致

| 命令 | 多等位行为 | 证据 |
|---|---|---|
| `stats`（QC 报告） | **静默折叠** `MultiallelicPolicy::CollapseAlt`；`.afreq` 写单一 `ALT_FREQS` 列 | `reed/src/stats/run.rs:298-304`；`reports.rs:1383` |
| `popg` | 站点门直接跳过（计入 `n_multiallelic_skipped`） | `reed/src/stats/run.rs:452` |
| `filter --fill-tags AC/AN/MAF/HWE` | 报 usage 错，拒绝 | `reed/src/filter/mod.rs:543-554` |

V27-39 已拍板"项目不做多等位"，但该拍板只覆盖了 `popg`；`stats` 的 QC 报告仍以
`CollapseAlt` 静默折叠 —— 对多等位位点，AC/AF/MAF/HWE 相对 sgkit 的 per-allele 定义不成立。
三者应统一：要么都拒绝（符合 V27-39），要么都明确文档化折叠语义。

#### 文档层（`doc`）

- **V30-04** f-stat `estimate=Σ per-site` 未在任何用户文档说明（§2.2-4）。
- **V30-05** `--f4-ratio` 角色序缺失（§2.2-6）。
- **V30-06** `f3*` 的 ΣT/ΣB 承诺不可达（§2.2-5）。
- **V30-07** MSA 全缺失→REF 未在 `convert.md` 显著警示（§2.2-13）。
- **V30-08** `identity_by_state` 与 reed `--type ibs` 定义不同，命名易误导（§2.2-11）。

---

## §3 实现但 sgkit 没有对应 API（用户点名）

以下能力在 reed 有实现与对拍，但 **sgkit 的公开 API 里没有对应函数**（因此
"覆盖 sgkit API"的分子/分母里都没有它们；它们是 reed 相对 sgkit 的**超出部分**）。

| # | 能力 | reed 位置 | sgkit 对应 |
|---|---|---|---|
| 1 | **WC84 Fst**（含 r-population 全局 `--fst-global`） | `reed-stats/src/fst.rs:14-19`；`WcMultiAcc` | `Fst` 只有 hudson/nei |
| 2 | **未归一化 PBS（Yi2010）** | `reed-stats/src/pbs.rs:15-20`；`--pbs` | `pbs` 只输出归一化 PBSn1 |
| 3 | **Watterson θ 作为一等输出** | `reed-stats/src/theta.rs:15`；`popg --theta` | 仅在 `Tajimas_D` 内部计算 |
| 4 | **SFS**（1D/2D 联合、folded/unfolded、hypergeometric 投影） | `reed-stats/src/sfs.rs:15`；`--sfs --proj` | 无 SFS API |
| 5 | **f2 / f3 / f4 / D / f4-ratio**（Patterson 无偏校正 + 块刀切） | `reed-stats/src/patterson.rs:21-87`、`jackknife.rs`；`--f2/--f3/--f4/--dstat/--f4-ratio/--block-size` | 无 f-statistics |
| 6 | **LD decay**（r²-距离分箱曲线） | `reed-stats/src/ld.rs`；`--ld-decay` | 只有 `ld_matrix`/`ld_prune` |
| 7 | **距离/亲缘矩阵族**：KING-robust、PI_HAT、phi(vcftools relatedness2)、Ajk(Yang)、discordance(bcftools gtcheck)、allele-ct、IBS/1-IBS、cov | `reed-cli/src/translate/distance.rs:30-49`（10 个 `--type`） | 只有 `pairwise_distance`(euclidean/correlation) 与 `identity_by_state` |
| 8 | **建树**：UPGMA / NJ + 位点 bootstrap clade support | `reed-graph/src/tree.rs:100,203`、`bootstrap.rs:110` | 无建树 |
| 9 | **降维族**：PCoA/kernel PCA/Laplacian/diffusion/Isomap/LLE/NMF/PPCA-EM/UMAP/t-SNE/DAPC/xval-DAPC | `reed-reduce/src/lib.rs:45-55` | 只有 `pca` |
| 10 | **QC 报告族**：`.afreq/.vmiss/.smiss/.hardy/.het/.depth/.summary`、直方图+精确分位数、singleton/doubleton 表、`--pop` 分层 | `reed/src/stats/writers/reports.rs` | `variant_stats`/`sample_stats`/`hardy_weinberg_test`/`individual_heterozygosity`，无报告布局 |
| 11 | **VCF/BCF 原生 IO**：plain/BGZF/gzip 嗅探、BCF 原生解码编码、CSI/TBI 构建、region 查询 | `reed-io/src/format.rs:31-73`、`write/index.rs:19-21` | **已移除** VCF/PLINK/BGEN 读写（交 bio2zarr/vcztools） |
| 12 | **硬过滤命令族**：掩码、QC 阈值、LD pruning、表达式语言、`--fill-tags`、样本/位点选择 | `reed-filter/**` | 只有 `filter_partial_calls`/`ld_prune` |
| 13 | **convert 族**：PHYLIP/FASTA/NEXUS/二元 NEXUS/012/hapmap/MSA/参考 MSA/SFS 导出、hapmap→VCF 反向 | `reed/src/convert/job.rs`、`convert/mod.rs:25+` | 无格式转换 |
| 14 | **rename/reheader 族**：header/样本名/contig/ID 编辑 | `reed/src/rename/mod.rs` | 无 |
| 15 | **index 族**：CSI/TBI | `reed-io/src/write/index.rs:19-21` | 无 |
| 16 | **pixy 口径的 accessible-bp 分母**（`--variant-only` / `--bypass-invariant-check`） | `reed-stats/src/pi.rs`；`--variant-only` | 无对应开关 |

> 第 11–15 项是"sgkit 主动移除、reed 自研"的能力，按用户口径仍属"实现但 sgkit 无 API"。

---

## §4 缺口清单（27 项未实现）

| 分组 | 项 | 是否属"统计内核" |
|---|---|---|
| 关联（4） | `gwas_linear_regression`、`regenie`、`regenie_loco_regression`、`genee` | 否（独立产品） |
| 系谱（8） | `pedigree_contribution`、`pedigree_inbreeding`、`pedigree_inverse_kinship`、`pedigree_kinship`、`hybrid_relationship`、`hybrid_inverse_relationship`、`parent_indices`、`pedigree_sel` | 否（独立产品） |
| 概率/倍性（5） | `convert_probability_to_call`、`infer_call_ploidy`、`infer_sample_ploidy`、`infer_variant_ploidy`、`call_allele_frequencies` | 半（需数据层决策） |
| 亲缘（4） | `pc_relate`、`Weir_Goudet_beta`、`invert_relationship_matrix`、`count_call_alleles` | **是**（`pc_relate` 在结构群体下不可替代） |
| 选择/杂合（2） | `Garud_H`、`individual_heterozygosity` | **是** |
| 窗口/展示（4） | `window_by_variant`、`display_genotypes`、`simulate_genedrop`、`simulate_genotype_call_dataset` | 否 |

剔除关联 4、系谱 8、概率/倍性 5、窗口/展示与模拟 4、Zarr/Dataset N/A 4（合计 25 项非内核），
剩 **31 项统计内核**：其中 **25 项有实现**（12 完整 + 13 部分）、**6 项未实现**
（`pc_relate`、`Weir_Goudet_beta`、`Garud_H`、`individual_heterozygosity`、
`invert_relationship_matrix`、`count_call_alleles`）→ **25/31 ≈ 81%**。

---

## §5 边界与未核实

1. **无 sgkit 数值夹具**：§2.1 的"等价"是公式层结论。要对外声称兼容 sgkit，需新增
   sgkit 夹具做直接对拍（见 02 V30-09）。
2. **未核实项**：sgkit `_Fst_Nei` 与 tskit/Nei GST 的精确关系（本组只确认 reed 未实现该式）；
   `Pca::loadings` 的缩放是否符合某外部工具的约定（CLI 未暴露）。
3. **本轮不裁决**功能取舍（是否要做 GWAS/系谱），只记录事实。
4. 本清单对 **HEAD+未提交改动** 有效。未提交改动覆盖 `reed-io`、`reed-filter`、
   `reed/src/stats/writers/reports.rs` 等；PI_HAT（`reed-distance`）、`stats_policies`
   与 `-I` 逻辑均已按当前工作树逐行复核，且 V30-01/02 用当前源码构建的二进制复现。

---

# 第二部分：API 工程质量（V30-Q01…V30-Q36）

> 2026-11-11 并入 doc_v30。对象：reed-rs
> 作为**可被外部使用的 Rust API 库**是否"工程上优雅"。功能覆盖与科学符合性见本文件 §1–§5。
>
> **方法**：6 个只读子代理按维度审查（公开面、类型/错误、文档、人体工学、契约/CI、
> 发布/依赖），主代理逐条复核关键结论，并实测三项：`origin/main` 树内容、
> `cargo check -p reed-reduce --no-default-features`、`cargo doc -D warnings`。
>
> **证据标记**：`repro` = 本次实测；`code` = 确定性代码路径/清单可证；`doc` = 文档与
> 代码不一致；`unverified` = 子代理结论，主代理未复跑。
>
> **基线**：HEAD `f190708`（分支 `dev`）；`origin/main` 是 `ci/publish.sh` 生成的公开子集。
> **编号**：V30-Q01…Q36（Q = quality），与第一部分的 V30-01…09 不冲突。

## §6 先定义"优雅"

一个 API 库的优雅可以拆成 8 个可检查的性质。本轮按这张表找缺口：

| 性质 | 含义 | 本仓现状 |
|---|---|---|
| 可发现 | 常用类型在根/`prelude` 可达，命名有规律 | 部分（无 prelude，`FilterSpec` 不在根） |
| 可命名 | 公开签名里不出现调用方无法命名的类型 | **不成立**（`StatOutput`/`BcfDict`） |
| 可构造 | 合法请求容易构造，非法状态不可表示 | **不成立**（全 pub 字段 + Default） |
| 可消费 | 结果是数据，写出是显式一步 | 部分（popgen/filter 是数据，stats 是文件） |
| 可诊断 | 错误可编程分类、带位置 | 部分（仅 Usage/Runtime，`Foreign` 丢类型） |
| 可演化 | semver、`#[non_exhaustive]`、公共 API 快照 | **不成立**（无快照/无 semver 政策） |
| 可验证 | doctest、契约测试、CI 门 | **不成立**（3 个可运行 doctest；main 只 build） |
| 可发布 | `cargo add` 可用、元数据完整、打包成功 | **不成立**（门面 `publish=false`） |

结论先行：**骨架是库的，脸面还不是。** 分层、sealed `Error`、零 `pub mod`、typed request、
602 个 `#[must_use]` 这些做得好；但"可命名/可构造/可演化/可验证/可发布"五项硬指标未过。

---

## §7 公开面与封装

### V30-Q01（P0，`code`）门面把整棵下层 crate 当公共 API

- **位置**：`crates/reed/src/lib.rs:138-176`（30 条 `pub use`，其中 11 条是整模块/整 crate 别名）。
- **现状**：`pub use reed_distance as distance`、`reed_graph as graph`、`reed_reduce as reduce`、
  `reed_io::{bgzf, chunk}`、`reed_io::write::atomic as atomic_file`、`reed_io::header as header`、
  `reed_engine::{budget, metrics}`、`reed_filter::{expr, funnel, plan, spec}`、7 个 `writers::*`
  子模块。展开后门面根约 82 个命名项 + 20 个模块别名 + **≈410 个下层 pub 项**可直达。
  其中 `reed::reduce` 全仓**零引用**（`rg -n "reed::reduce" crates` 无命中）。
- **为什么不是优雅 API**：下层 crate 的任何 minor 改动都成为门面的 breaking change；而门面
  自己 `publish=false`，既没有 semver 承诺，又拥有巨大表面。用户面对 `reed::reduce::extra::*`
  这种路径，无法判断哪些是稳定契约。
- **建议**：逐项再导出（`pub use reed_io::chunk::region::Region;`、
  `pub use reed_engine::budget::{ThreadBudget, MemBudget};` 等）；`distance/graph/reduce`
  先有驱动函数（V30-Q16）再决定是否导出；加守卫：禁止裸 `as <module>;` 别名。
- **附**：`lib.rs:74` 的 `#![allow(unreachable_pub)]` 关闭了整 crate 的可见性检查，公开面
  完全靠人工维护。

### V30-Q02（P0，`code`）公开签名里出现调用方无法命名的类型

- **位置**：`crates/reed/src/stats/request.rs:264`（`pub outputs: Vec<reed_stats::request::StatOutput>`）、
  `crates/reed/src/input.rs:58`（`bcf_dict() -> Option<&Arc<reed_io::bcf::BcfDict>>`）、
  `crates/reed/src/stats/writers/reports.rs:201`（`pub qc: Option<QcReportAcc>`，其 12 个字段
  类型在 reed-stats 深层且未导出）。
- **现状**：`StatOutput`、`BcfDict`、`QcReportAcc` 的字段类型都不在门面清单（`rg` 零命中）。
  本仓测试已经绕过门面：**13 个** `crates/reed/tests/*.rs` 直接
  `use reed_core/io/engine/stats/filter`，其中 `e2e_filter.rs` 79 处、`e2e_bcf.rs` 57 处。
- **为什么不是优雅 API**：API 的第一个外部消费者（集成测试）都不得不走旁路，说明门面
  作为 API 不成立；用户无法为 `match` 出的值写函数签名。
- **建议**：`pub use reed_stats::request::{StatAcc, StatOutput};`（或把 `accs` 换成稳定结构）；
  `bcf_dict` 删除或返回 opaque 包装；加一条文档门："公共签名中出现的类型必须可在 `reed::*` 命名"。

### V30-Q03（P1，`code`）`doc(hidden)` 只有 4 处，且被唯一前端依赖

- **位置**：`lib.rs:145`（`open_regions`）、`reed-filter/src/plan/strategy.rs:13`（`ExecStrategy`）、
  `reed/src/input.rs:64`（`source_mut`）、`reed/src/filter/mod.rs:67`（`FilterOptions::strategy`）。
- **现状**：`open_regions` 被 `reed-cli` 调用 **6 次**；`ExecStrategy` 经 `reed::plan` 整模块
  再导出，被隐藏 flag `--strategy` 使用。真正该隐藏的（`atomic_file` 整模块、
  `budget::BytePermits`、`header::MAX_HEADER_BYTES`、`expr::lexer/parser`、`spec::render`）没有隐藏。
- **为什么不是优雅 API**：rustdoc 上看不到的项恰恰是前端契约的一部分，文档与真实契约不一致；
  `doc(hidden)` 被当成"不好意思导出"的补丁，而不是稳定性工具。
- **建议**：给"前端专用"单独一个 `#[doc(hidden)] pub mod frontend`；其余要么正式文档化，
  要么移出公共面。

### V30-Q04（P1，`code`）根可达性不对称、同一类型多路径

- **位置**：`reed::spec::FilterSpec`（不在根）vs `reed::PopgenRequest`（在根）；
  `Mat` 同时 `reed::Mat`（`lib.rs:110`）与 `reed::reduce::Mat`（`reed-reduce/src/lib.rs:80`）；
  `Stand` 有 4 条路径（`reed-distance` / `reed-reduce::pca` / 门面别名）。
- **为什么不是优雅 API**：文档无法给出唯一写法；`use reed::*` 拿不到全部常用类型。
- **建议**：加 `pub mod prelude`；根上补 `pub use spec::FilterSpec`；每个类型只保留一条规范路径。

### V30-Q05（P1，`code`）区域类型三份，互斥字段只在运行期报错

- **位置**：`reed-filter/src/spec/region.rs:43`（`Region { contig, start, end }`）、
  `reed-io/src/chunk/region.rs:23`（`Region { contig, start1, end1 }`）、
  `reed/src/convert/mod.rs:109`（`MsaRegion { start1, end1 }`）；`PopgenRequest.regions` +
  `select` 并存（`stats/request.rs:359,364`）；`FilterOptions.regions` + `FilterSpec.targets.indexed`
  并存（两者都填时运行期报 usage 错）。
- **为什么不是优雅 API**：同名不同字段名、同义三份；"二选一"应由类型表达，而不是错误消息。
- **建议**：唯一 `Region` + `enum RegionSource { None, Seek(Vec<Region>), Predicate(SiteSelector) }`。

### V30-Q06（P2，`code`）没有 `prelude`

`rg "pub mod prelude"` 全仓 0。常用类型散在根/模块/下层别名中。建议加一个
只含请求、报告、错误、`ThreadBudget`、`WindowGrid`、`Mat` 的 prelude。

---

## §8 请求 / 类型 / 错误

### V30-Q07（P0，`code`）`PcaRequest`/`DapcRequest` 的 `threads`/`memory` 是死字段

- **位置**：`crates/reed/src/pca/mod.rs:272-273`、`crates/reed/src/dapc/mod.rs`（同形字段）；
  CLI 真的会写它们（`reed-cli/src/translate/pca.rs:86-87`）。
- **现状**：`rg "request\.(threads|memory)" crates/reed/src/pca crates/reed/src/dapc` **零命中**；
  真实并行度来自 `run_pca(..., &ThreadBudget)`。`PcaRequest::default()` 把 `prefix` 设为
  `"out"`（`pca/mod.rs:270`），库用户不改就会在当前目录写文件。
- **为什么不是优雅 API**："看起来能控制、实际无效"的字段比没有字段更糟；且引擎已有受约束
  newtype `Threads`/`MemMiB`，请求却用裸 `Option<usize>`/`Option<u64>`。
- **建议**：删除死字段（预算只走 `&ThreadBudget`），或让 `run_*` 真正读取并复用 newtype；
  `prefix` 改为必填或 `Option`（不默认写 CWD）。

### V30-Q08（P1，`code`）请求结构"全 pub 字段 + Default"，非法状态可表示

- **位置**：`FilterSpec`（`reed-filter/src/spec/mod.rs:30`）、`PopgenRequest`
  （`reed/src/stats/request.rs:316`）、`PcaRequest`、`DapcRequest`、`ConvertRequest`、
  `StatsReportRequest`（`reports.rs:40`）、报告类型（全部 `pub` 可变字段）。
- **现状**：
  - `PopgenRequest::default()` 的 `kinds=[]` → 驱动走 `GridMode::None`
    （`stats/driver.rs:204-209`）→ **静默成功、0 行**，不报错。
  - `FilterSpec.policy.multiallelic=Reject` 可与 `site.max_alt_alleles=Some(>1)` 自相矛盾
    （`multiallelics_allowed()` 反而为 true）。
  - `Levels { window:false, contig:false, genome:false }` 可直接构造，绕过 `Levels::parse`
    的拒绝。
  - `min_maf > max_maf`、`LdPruneSpec.r2` 越界等均可在构造期表示。
- **为什么不是优雅 API**：库用户"填 Default 再改几字段"是主流用法；非法组合应不可构造，
  而不是运行期才 Usage 错或静默产空结果。
- **建议**：请求字段私有 + `new()`/builder + `validate() -> Result<Self, Error>`；
  非法组合（maf 区间、policy/alt 数）在构造期拒绝；报告类型提供只读访问器。

### V30-Q09（P1，`code`）枚举用 `String` 承载

- **位置**：`PcaRequest.stand/missing/solver`（`pca/mod.rs:248-260`）、
  `DapcRequest.stand/missing/solver/result`、`XvalRequest.result`。
- **现状**：非法值在 `run_pca_opened`/`run_dapc_opened` 里用 `matches!` 运行期校验
  （`pca/mod.rs:326-338`）。同一门面的 `PopgenRequest` 却用 typed `StatKind/FstKind/PbsKind`。
- **为什么不是优雅 API**：拼错只有运行时才知道；同一 API 面两套哲学。
- **建议**：`enum Stand/Missing/Solver/XvalResult`，`FromStr` 仅供 CLI 边界使用。

### V30-Q10（P1，`code`）`FromStr` 几乎不存在，`parse` 返回四种错误类型

- **现状**：全仓 `impl FromStr` **2 个**；但有一批 inherent `parse(&str)`：
  `MissingRate::parse -> Result<_, FilterError>`、`Levels::parse -> Result<_, reed::Error>`、
  `Stand::parse -> Result<_, String>`、`WriteIndex::parse -> Result<_, IoError>`。
- **为什么不是优雅 API**：`"0.1".parse::<MissingRate>()` 不可用；`?` 不能直通 `reed::Error`，
  调用方每次都要 `map_err`。
- **建议**：核心类型统一 `impl FromStr<Err = reed::Error>`（或单一 `ParseError`）；
  inherent `parse` 保留为别名或删除。

### V30-Q11（P1，`code`）错误分类不可编程

- **位置**：`crates/reed/src/error.rs:43`（`ErrorKind` 非 `#[non_exhaustive]`）、
  `crates/reed-filter/src/error.rs:11`（`FilterError` 非 `#[non_exhaustive]`）、
  `error.rs:38`（`Foreign(Box<dyn Error>)` 丢类型）。
- **现状**：`Error` 本体 sealed 是**正确**的；但 `ErrorKind` 加第三类会破坏调用方穷尽 match；
  下层 `FilterError` 未 non_exhaustive 却参与门面 `kind()` 的穷尽匹配；`Foreign` 只剩字符串，
  无法 `downcast_ref`；结构化 `locus`（`FilterError::Decode { contig, pos }`）不可达；
  只有 Usage/Runtime 两类，对库消费者过粗。
- **建议**：`ErrorKind` 加 `#[non_exhaustive]`；`FilterError` 加 `#[non_exhaustive]`；
  提供 `Error::locus() -> Option<&Locus>` 或把位置类型 re-export；`Foreign` 保留 source 类型。

### V30-Q12（P2，`code`）请求结构未 `#[non_exhaustive]`

全仓仅 4 处 `#[non_exhaustive]`，都在错误枚举。请求/报告结构加字段即 breaking，
与 semver 演化直接冲突。

### V30-Q13（P2，`code`）`Mat` 暴露 `pub rows/cols/data`

`reed-core/src/linalg/mod.rs:9-15`：外部可构造 `rows*cols != data.len()` 的矩阵，
破坏所有下游 kernel 的不变量。建议字段私有 + `from_vec` 校验 + `as_slice()`。

### V30-Q14（P2，`code`）`FilterOptions` 混入测试钩子与重复格式字段

`filter/mod.rs:67` 的 `#[doc(hidden)] pub strategy` 是测试钩子，却在一个用户可见请求结构里；
`output_compressed: Option<bool>` 与 `output_bcf: Option<u8>` 两个字段编码同一"输出格式"决策。

### V30-Q15（P2，`code`）`Populations` 两套表示 / `Option` 语义混用

`PopgenRequest.populations: Populations`（非 Option，`All`=全体）vs
`StatsReportRequest.populations: Option<Populations>`（`None`=全体，`All` 也=全体）。
另有 `Option` 同时表示"默认/关闭/未知"（`layout: Option<Layout>` 实际是模式开关）。

---

## §9 入口 / 结果 / 易用性

### V30-Q16（P0，`code`）`run_distance` 缺失，前端承担后端工作

- **位置**：`crates/reed-cli/src/translate/distance.rs:56-472`（约 472 行）；
  `rg "pub fn run_distance|pub struct DistanceRequest"` **零命中**。
- **现状**：前端自己选矩阵函数、构 `BootstrapConfig`、跑 `graph::support`、写 `.nwk`/矩阵/ID
  文件，并匹配 `GraphError` 变体做 usage/runtime 分类（`distance.rs:22-28`）。
- **为什么不是优雅 API**：这是 `v5/08 §6.1`"前端只构造请求"唯一没落地的命令，也是 V30-Q01
  整 crate 别名的根因；库用户根本拿不到 distance 能力（只能自己去 `reed::distance` 拼）。
- **建议**：`DistanceRequest`/`DistanceReport` + `run_distance(input, &request, &budget)`；
  把 `GraphError` 分类搬进门面；随后删除 `reed::distance`/`graph`/`reduce` 别名。

### V30-Q17（P1，`code`）`run_*` 家族签名不统一

- **现状**：
  - `run_filter(input, &spec, &FilterOptions)`：无 `&ThreadBudget`（预算藏在 `FilterOptions`）。
  - `run_index(input: &Path, flavor, output, threads: usize)`（`index.rs:524`）：无 request
    结构、无 budget、不能 stdin。
  - `run_rename(input, spec, output: &RenameOutput, budget)`：比其他多一个 `output`。
  - `_opened` 变体：pca/dapc/rename/convert/popgen/stats 有，**filter 没有**。
  - "打开输入"两套：`crate::input::open` 与 `open_regions(path, &[])`（pca/mod.rs:309）。
- **建议**：统一 `run_<cmd>(input: Input<'_>, request: &XRequest, budget: &ThreadBudget)`；
  补齐 `run_filter_opened`；`open`/`open_regions` 合并。

### V30-Q18（P1，`code`）`run_stats` 是"文件进、文件出"

- **位置**：`StatsReportRequest.prefix: PathBuf`（`reports.rs:88`）、
  `StatsOutcome.files: Vec<PathBuf>`（`reports.rs:194-207`）。
- **现状**：`run_stats` 必然落盘，返回路径；用户想对 `.afreq` 编程处理只能自己 `fs::read`。
  同一门面的 `run_popgen` 返回 `PopgenReport.windows: Vec<WindowRow>`（数据），
  `run_filter` 返回 `FilterReport`（数据）——"结果是什么"取决于命令。
- **建议**：统一"数据 in → 数据 out，落盘显式 `write_*`"；报告族提供 `Vec<Row>` 迭代器。

### V30-Q19（P1，`code`）构造样板大，缺便捷构造器

- **位置**：`window_kinds(n_groups, pi, dxy, theta, tajima, fst, pbs)`
  （`stats/request.rs:94-112`，7 个位置参数、4 个 bool，已被
  `#[allow(clippy::fn_params_excessive_bools)]` 豁免）。
- **现状**：CLI 构造一个窗口 `PopgenRequest` 约 40 行/17 字段
  （`translate/popg.rs:582`）；无 `StatKind::pi(group)` / `PopgenRequest::pi_scan(...)`。
- **为什么不是优雅 API**：位置 bool 无类型安全；"一行做常见事"的入口缺失。
- **建议**：`StatKind` 加构造器；`PopgenRequest` 加 `pi_scan/dxy_scan/...` 预设；
  `window_kinds` 改为接收 `&[WindowStat]`。

### V30-Q20（P1，`code`）单位解析锁在前端，库只收裸数字

`--window 100k` 的 `parse_size` 是 `reed-cli/src/units.rs:12` 的 `pub(crate)`；
库调用者只能写 `WindowGrid::Bp { size: 100_000, step: 100_000 }`。建议在 `reed-core`
放 `Length`/`WindowSize`（含 `FromStr` 与 kb/Mb 后缀），CLI 复用之。

### V30-Q21（P2，`code`）`to_args` 四种签名，仅 4 个请求实现

`FilterSpec::to_args(&self) -> Vec<String>`、`WriteIndex::to_args(self) -> Vec<String>`、
`PopgenRequest::to_args(&self, group_names) -> Result<...>`、
`ConvertRequest::to_args(&self) -> Result<...>`。`PcaRequest/DapcRequest/RenameSpec/
StatsReportRequest` 无文本契约；`default_kinds` 产出的 kind 无 flag 拼写，`to_args` 直接
返回 Usage 错（`stats/render.rs:92`），"便捷默认集"不可往返。

### V30-Q22（P2，`code`）请求夹带前端/渲染关注点

`PopgenRequest.window_output`（含 prefix/command/version，是 CLI 流式写盘）、
`compact_window_accumulators`（CLI 性能开关）、`layout: Option<WindowLayoutOptions>`
（`None` = 另一条渲染路径）。库默认请求不应表达这些。附带：`PcaRequest::default()`
`prefix="out"` 会静默写 CWD。

### V30-Q23（P2，`code`）`FilterSpec` 表达式存原文 `Vec<String>`，每次运行重解析

`SiteRules.include/exclude`、`TailRules.include/exclude`、`InfoPredicate.keep_if` 等均为
`String`；已解析 AST 不缓存。建议请求存 AST、原文只用于 `to_args`。

---

## §10 文档 / 示例

### V30-Q24（P0，`code`）零用户级 doctest

- **现状**：10 个 lib crate 中 9 个 `lib.rs` **零 doctest**；全仓仅 **3 个可运行 doctest**
  （`reed` lib 2 + `reed-core::rng` 1）+ 3 个 `no_run`；`# Examples` 段 **2 个**，
  `# Arguments`/`# Returns` **0 个**。
- **examples/**：仅 `reed-reduce/examples/` 6 个文件，文件头自述为 benchmark/crosscheck/
  conformance，不是用户示例；其余 9 crate 无 examples。
- **为什么不是优雅 API**：rustdoc 目前是"内部架构索引"，不是用户文档；公开 API 无法自学。
- **建议**：每个 crate 至少 1 个 crate 级可运行示例 + 每个主要入口 1 个 `# Examples`；
  `README` 用 `doc = include_str!` 纳入编译门。

### V30-Q25（P1，`doc`）README 与 docsite 的库定位不一致

- 根 README 全是 `reed ...`/`cargo install`，无 `cargo add`、无库示例、无 feature 说明
  （MSRV 声明有：Rust 1.91）。
- `docsite/src/content/docs/library.md:16-19` 称"crate root exports **six** entry points"
  （run_filter/run_popgen/run_stats/run_rename/run_index/run_convert），实际根上还有
  `run_pca`、`run_dapc` 及各自 `_opened`（≥8）；同页称下层 crate 是"internal, may change"，
  但门面 `pub use reed_distance as distance` 等把整 crate 公开。示例块是普通 markdown
  fence，**不进 doctest**，会持续腐烂。

### V30-Q26（P1，`code`）公开 rustdoc 被内部票号污染，且无稳定性声明

`lib.rs`、`pipeline.rs` 等公开文档大量引用 `doc_vNN`/`D7-22`/`V26-08`/`A11-2`，外部读者
无法解析且无映射表。全仓 `# Stability` = 0；`docs/rounds` 与 `docsite` 对"稳定什么"
三种说法（CLI+输出 / 仅 crate root / 只承诺 CLI）。

### V30-Q27（P2，`code`）oracle crate 文档首行错误

`crates/oracle/src/lib.rs:1` 仍写 "`reed-core`: canonical population-genetics statistics"，
并指向不存在的 `reed/docs/rounds/v1/design.md`；Cargo.toml 自述是冻结 0.2 oracle。
该 crate `missing_docs=allow`，错误长期无人发现。

### V30-Q28（P2，`code`）文档 lint 无回压

`missing_docs = "warn"`（且 oracle allow）；`missing_errors_doc`/`missing_panics_doc = "allow"`。
Result 函数的 `# Errors` 覆盖不均，`# Panics` 少。建议门面与公开 crate 提级为 deny
（至少 `missing_docs`），并给 `run_*` 强制 `# Errors`。

---

## §11 版本 / 发布 / CI

### V30-Q29（P0，`code`）发布定位与"API 库"直接矛盾

- **现状**：门面 `reed` 与 `reed-core/io/engine/stats/filter/cli/oracle` **`publish=false`**；
  反而 `reed-distance`/`reed-reduce`/`reed-graph` **未写** `publish=false`（默认可发布），
  但 `reed-reduce` 依赖 `reed-distance = { path = ... }`、workspace 三者也只有 path 无
  `version` → `cargo package` 必失败。全仓 0 个 `repository/readme/keywords/categories/
  documentation/[package.metadata.docs.rs]`。没有任何 crates.io 发布通路（`release.yml` 只发
  GitHub tar；`publish.sh` 是 git force-push）。
- **为什么不是优雅 API**：用户不能 `cargo add reed`；"可发布"的三个 crate 恰是最没准备好的三个。
- **建议**：先定发布策略（发布门面 + 必需闭包，或明确"仅 git 依赖"并删除可发布暗示）；
  补元数据与 crate 级 LICENSE；加 `cargo package --dry-run` 门。

### V30-Q30（P1，`code`）版本 / edition / lint 三套标准

`reed-distance`/`reed-reduce`：`version = "0.1.0"`、`edition = "2021"`、无 `rust-version`、
无 `[lints] workspace = true`；`reed-graph` 无 `rust-version`；其余 8 个走
`version.workspace`(0.4.0)/edition 2024/1.91。workspace 注释自称 "lockstep versioned"，
与实际双轨不符。后果：distance/reduce 不受 `unsafe_code=deny`/`unwrap_used=deny`/
`missing_docs` 约束。

### V30-Q31（P1，`repro`）`--no-default-features` 承诺是假的

`cargo check -p reed-reduce --no-default-features` **失败**：

```
error[E0063]: missing field `parallel` in initializer of `Operator<'a>`
  --> crates/reed-reduce/src/lanczos.rs:146:13
```

另有 `backend.rs:155/206` 未使用参数警告。`reed-reduce` 的 lib 文档与 Cargo feature 承诺
"dependency-free fallback"，CI 无该矩阵。`reed` 门面也无条件依赖 reduce，下游无法关闭 faer。

### V30-Q32（P1，`code`）CI 门错位

`.github/workflows/ci.yml:12-26`：`build` job 条件 `(base_ref||ref_name)=='main'` **只
`cargo build`**；`gates` job 条件 `!= 'main'` 才跑 test/clippy/doctest/doc-drift。
因此 **push 到 main 与所有以 main 为目标的 PR 不跑测试**。另无 `cargo doc`、
`cargo fmt --check`、`cargo deny/audit`、`--no-default-features` 矩阵；
`ci/check-structure.sh` 无任何 workflow 调用（`rg` 仅命中脚本自身注释）。

### V30-Q33（P1，`repro`）公开树被主动剥离 tests/docs

`ci/publish.sh` 的 `PUB_PATHS` 不含 `docs/`/`bench/`，并显式
`':!crates/*/tests/'`。实测：`origin/main` 259 文件、**0** 个 `/tests/`、**0** 个 `docs/`、
无 `CHANGELOG.md`；`origin/dev` 579 文件、**65** 个测试、**151** 个 docs。
若 main 是访客默认分支，则 API 库的"门面页"没有测试与文档。
（这是发布策略，不是事故；但对"API 库"定位是反向的。）

### V30-Q34（P2，`repro`）忽略规则静默丢弃新文件（`crates/*/tests/`、`docs/`）

7 个 crate 的 `.gitignore` 内容为 `tests/`（`crates/reed` 为 `tests/*` + `!tests/support/`）。
实测 `git check-ignore -v crates/reed/tests/new_contract.rs` 命中：
新测试文件会被 `git add` 静默忽略，直到 CI/发布才发现缺失。

同一个坑还有根 `.gitignore:19` 的 `/docs/`（注释写 "committed on dev, never merged to
main"，实际靠 `git add -f`）：**本轮产物 `docs/rounds/v30/`
在 `git status` 里完全不可见**，`git ls-files docs/rounds/v30` 为 0；已跟踪的 151 个
docs 文件是历史遗留。`docs/README.md` 的索引行是已跟踪文件，所以只有它显示为 `M`。
后果：`git add -A`/`git commit -a` 不会带上新轮次文档，记录可能静默丢失。

### V30-Q35（P2，`repro`）rustdoc 已坏且无门

`RUSTDOCFLAGS='-D warnings' cargo doc --no-deps -p reed-distance -p reed-graph -p reed-io`
**exit 101**：

- `reed-graph`：public doc 链接私有 `MAX_DEPTH`；
- `reed-distance`：`unresolved link to covariance`、`unresolved link to 0,1`、
  `` `grm` is both a function and a module ``。

CI 无 `cargo doc` 门，所以这些长期存在。

### V30-Q36（P2，`code`）依赖卫生与许可证打包

`deny.toml` 存在但无任何 workflow 调用（`cargo deny/audit` 空转）；workspace 依赖全 caret
无策略；`Cargo.lock` 同名多版本（`thiserror` 1/2、`syn` 2/3、`getrandom` 0.3/0.4 等 7 组）。
11 个 crate 均无 crate 级 LICENSE 文件（仅根有 `LICENSE-MIT`/`LICENSE-APACHE`）。

---

## §12 正面结论（避免只记录坏的）

这些已经做得对，后续修复不应破坏：

- 门面 `crates/reed/src/lib.rs` **0 个 `pub mod`**，实现模块全私有（`mod`），符合 v5/08 §6.2。
- 后端 crate **0 处 `println!/eprintln!`**；错误/日志走 `tracing` 与 typed error。
- `Error` 本体 sealed（`Box<ErrorInner>`），`source()` 链完整；`ErrorKind` 的 Usage/Runtime
  两分类对 CLI 退出码是正确的。
- 请求/报告类型普遍 `Clone + Debug + Default`，`#[must_use]` **602 处**。
- L0–L5 单向分层清晰；`reed-stats`/`reed-filter` 保持纯策略（无 I/O、无线程）。
- `ci/check-doc-drift.sh` 98 条守卫 + `ci/check-structure.sh` 16 条，内部一致性纪律强。
- 602 个 `#[must_use]`、`unsafe_code = deny`、`unwrap_used = deny`（对 8 个 crate 生效）。

---

## §13 边界与未核实（第二部分）

1. **未复跑**：`cargo test --workspace`（全量测试耗时）；`cargo package --dry-run`（结论
   由 manifest 的 path 无 version 推出，标 `code` 而非 `repro`）。
2. **子代理结论、主代理抽查**：`FilterSpec` 的 policy/alt 矛盾细节、`QcReportAcc` 字段树、
   `run_stats` 的 `compact_window_accumulators` 路径、`deny.toml` 的具体 license 白名单。
   实施前应再核对应行。
3. **不裁决**：是否要做 crates.io 发布、是否公开 tests/docs、是否保留三个 0.1.0 crate 的
   独立版本——这些是产品决策，本清单只记录现状与代价。
4. 本轮**不改代码**。修复项见 [02-roadmap.md](02-roadmap.md)。


---

# 第三部分：根因审查（V30-R01…V30-R14）

> 2026-11-11 追加。方法来自用户建议的四项根因审查：**语义归属 / 能力接口 / 生命周期与所有权 /
> 变更传播**。目的不是再加症状，而是把"同一条规则反复在多模块修正、加一种输入要改好几层、
> 边界判断不断复制"这三类现象归到**同一项根因**，再决定是否需要调整数据模型/能力接口/状态
> 所有权。
>
> 写法：每条给 **证据 → 根因 → 建议归属 → 如何验证**。证据级：`repro`（主代理本轮实测，
> 命令与输出见 [STATUS.md](STATUS.md) §6）、`code`（确定性代码路径/计数）、`commit`（git 历史）、
> `hypothesis`（待验证）。**本轮只审查，不改代码**。
>
> 严重度（本部分专用）：**P0 = 会产出错误数值/错误结果或使结构约束失效**；P1 = 结构性根因，
> 造成重复修改与误用风险。共 **14 条：4 P0 / 10 P1**。

## §14.1 语义归属（缺失值/半缺失、坐标、等位基因由谁定义）

### V30-R01（P0，`repro`）GT/等位语义没有单一所有者：5 份 ALT 计数、4 个 GT 读取者、convert 自定"有数据"

- **证据**
  - `crates/reed-io/src/vcf/site.rs:22-34` 的注释声明 `count_alts` 是 "the single implementation
    of the rule"，实际至少 5 份：`crates/reed/src/stats/writers/summary.rs:86-92`、
    `crates/reed/src/stats/writers/reports.rs:634-639`、`crates/reed/src/convert/job.rs:824-829`、
    `crates/reed-filter/src/eval/head.rs:249-254`、`crates/reed/src/genotype.rs:263-269`
    （最后一份用 `split(b',').count()`，空 ALT 的语义与 `count_alts` 不同）。
  - GT 文本读取者至少 4 处：`crates/reed-io/src/vcf/genotype_scan.rs`、`crates/reed-io/src/bcf/gt.rs`、
    `crates/reed/src/genotype.rs:65-71`（`gt_subfield` + 自定 `dosage`）、
    `crates/reed/src/convert/job.rs:818-845`（自写逐样本 tokenize）；另有 5 条 BCF native 入口
    （`genotype.rs:112`、`stats/run.rs:833`、`filter/single_pass.rs:291`、`filter/pass0/matrix.rs:466`、`index.rs`）。
  - **分歧点**：`crates/reed/src/convert/job.rs:851` 用 `gt.present = !raw.starts_with(b".")`
    定义"样本有数据"，而 `crates/reed-core/src/genotype/gt_code.rs:146-157` 把 `x/.` 与 `./x`
    一律归为 `GtClass::HalfCall`（`crates/reed-core/src/count/site_counts.rs:163-164` 也把半调用计为 1 个已调用等位）。
- **实测**：8 样本两个位点（`./1`×4 + `./.`×4，与 `0/.`×4 + `./.`×4）：
  `reed convert half.vcf --phylip --min-samples-locus 4` 只保留 `0/.` 位点（输出 `8 1`）；
  `--012 --min-samples-locus 4` 的 `.012.pos` 只剩 `chr1 20`；
  而 `reed stats half.vcf --missing` 两个位点都是 `MISSING_CT=4 OBS_CT=8 F_MISS=0.5`。
- **根因**：等位/GT 规则的归属写在传输层注释里，没有类型；每个消费方按自己的局部需要重写一遍。
- **建议归属**：`reed-core`——`GtCode::has_called_allele()`（`class() != Missing`）作为唯一判定；
  新增 `AltCount` newtype 承载 ALT 计数；`reed-io` 只做借位解析，消费方不再二次计数。
- **如何验证**：convert 的 `present` 换成 `has_called_allele()` 后，断言 `--phylip/--012/--msa`
  与 `stats --missing` 在半缺失上一致；结构门让 `rg -n "fn count_alts|fn n_alt_of" crates` 收敛为 1。

### V30-R02（P0，`repro`）AC/AN 有第二份实现：`stats --summary` 手写 INFO 解析，且同一条命令内两个 AC 不同

- **证据**
  - 权威实现：`crates/reed-filter/src/semantics.rs:195-271`（`info_af_indexed`：`Number=A` 长度域、
    `AC` 溢出 checked、`ΣAC>AN` 报错）。
  - 第二份：`crates/reed/src/stats/writers/summary.rs:294-328`——`:323 let nac: u64 = ac.iter().sum();`
    无 checked 求和，`:324` 只做 `an >= nac`；`:299` 用 `info_fragments` 手扫，不区分 `AcLength/AcOverflow`。
  - 第三来源：同一命令的 `.afreq` 走 GT 计数（`crates/reed/src/stats/writers/reports.rs:843-857`）。
- **实测 1（溢出）**：`AC=18446744073709551615,1;AN=100` 的单记录 VCF：
  debug 二进制 `reed stats ovf.vcf --summary` → `attempt to add with overflow` →
  `reed: worker panicked on chunk 0`；release 二进制不报错，同一位点 `.stats` 的 AF 行给出 `0.000000`
  （求和回绕）。而 `reed-filter` 的单测 `info_af_rejects_ac_overflow_instead_of_wrapping`
  （`semantics.rs:376`）对同一输入是**拒绝**——两条路径结论相反。
- **实测 2（同命令两个 AC）**：`GT` 全 `0/0`、`INFO AC=5;AN=8` 的位点：
  `reed stats disagree.vcf --freq --summary` 产出 `.afreq` 的 `ALT_FREQS=0`，
  而 `.stats` 的 AF 段 `0.610000`。
- **根因**：allele count / AC-AN 没有拥有类型；"唯一实现"的声明放在 L1（io），消费方各自重写。
- **建议归属**：`reed-core::count::InfoCounts::from_info(index, info, n_alt) -> InfoAfOutcome`
  （含长度域与 checked 求和），filter 与 stats 共用；同时明确 `stats` 以 INFO 还是 GT 为准。
- **如何验证**：溢出夹具在两条路径返回同一个 typed 错误；`.afreq` 与 `.stats` 的 AF 同源，
  或在文档里显式声明差异。

### V30-R03（P0，`code`）BCF `pos0→pos1` 四处转换、三种哨兵语义；`coord` 类型零调用

- **证据**
  - `crates/reed-io/src/bcf/render.rs:76-89`：`pos0 == -1` → POS 0（htslib `vcf.c:1649-1650`），`<-1` 才报错。
  - `crates/reed/src/stats/run.rs:873-878`：任何负值 → `BCF POS is negative`（连哨兵也拒绝）。
  - `crates/reed/src/filter/record_pipeline.rs:440` 与 `crates/reed/src/genotype.rs:134`：
    `u64::try_from(header.pos0).unwrap_or(0) + 1`，`-1` 静默变成 POS 1。
  - `crates/reed-core/src/coord.rs:1-7` 要求"所有单位换算必须走 `From`"，但
    `to_pos0/to_pos1/to_interval0/to_interval1/Interval1` 在 `crates/` 内除定义处**零调用**，
    实际全是裸 `+1/-1`（如 `crates/reed/src/stats/sink.rs:853,896`、`convert/job.rs:1556,1564`）。
- **实测（子代理，主代理只复核代码路径）**：同一份 `pos0 = -1` 的 BCF 记录，
  `pca` 接受（`1 sites`）、`stats` 报 `BCF POS is negative`、`filter` 报 malformed —— 三种行为。
- **根因**：哨兵值由各消费方各自解释；`coord` 的转换类型存在但没被当作 API 使用。
- **建议归属**：`reed-io::bcf::record` 提供唯一的 `pos1()`（`-1` 显式命名为 telomere → POS 0，`<-1` 报错），
  内部统一 `reed_core::coord::{Pos0,Pos1}`。
- **如何验证**：`rg -n "try_from\(header\.pos0\)|pos0 == -1" crates` 收敛为 1；BCF 哨兵夹具下四命令行为一致。

### V30-R04（P1，`code`）区间规则三份类型、两份判定；convert `--msa --region` 用无 header 的启发式解析

- **证据**
  - `crates/reed-filter/src/spec/region.rs:113-125`（生产 `region_allows`）与 `:241-245`
    （`Region::contains`，仅测试 `:302-304` 调用）是同一 1-based closed 规则的两份实现。
  - 第二类型：`crates/reed-io/src/chunk/region.rs:23`（`start1/end1`）；第三类型：
    `crates/reed/src/convert/mod.rs:107-116` 的 `MsaRegion` + `crates/reed/src/convert/job.rs:127-131`
    的 0-based half-open `pos0 >= s && pos0 < e`。
  - 解析入口：`crates/reed-cli/src/translate/convert.rs:39-47` 调 header-less `Region::parse`，
    而该函数自己的文档说生产入口应带 contig 词表（`crates/reed-filter/src/spec/region.rs:190-203`）。
- **实测（子代理）**：参考 FASTA 的 contig 字面名是 `chr1:100` 时，
  `reed convert --msa --reference ref.fa --region 'chr1:100'` 报 "reference FASTA does not contain chr1"。
- **根因**：区间基址的所有权没有唯一类型；解析需要 header 词表，却由 L3/CLI 各自拼。
- **建议归属**：`reed-core::coord::Interval1` 为唯一区间；`Region` 移到 `reed-io`；
  解析统一 `parse_in(text, contig_names)`，`convert --msa` 传 `reference.keys()`。
- **如何验证**：`rg -n "fn contains|fn region_allows|fn allows" crates` 收敛为 1；
  补一个 `convert --msa --region` 用例（现有 `declared_contigs_disambiguate_colons` 只覆盖解析）。

## §14.2 能力接口（顺序读取/随机访问/可重读/可定位）

### V30-R05（P1，`repro`）"可重读/可寻址"不是类型能力，被压成 `path: Option<PathBuf>`；two-pass 策略在无 replay 能力时仍被选中

- **证据**
  - `crates/reed/src/input.rs:30-41`：`OpenedInput { source, path: Option<PathBuf> }`——
    stdin（一次性）与文件（可重开）在类型上完全同构；`:70-78` 的 "the caller must read it exactly once" 只是注释。
  - `crates/reed/src/filter/mod.rs:130-145`：`open_pass1` 的"重读"就是 `first.take()` 或按
    `replay_path` 再 `open_regions`。
  - `crates/reed/src/filter/mod.rs:345-374`：`None if is_bcf => ExecStrategy::TwoPassReparse`
    ——选择两趟策略时**不问**是否有 replay 能力（只有 `--mind` 分支才检查 `replay_path.is_some()`）。
- **实测**：`cat small.vcf | reed filter - --strategy two_pass_reparse --timing t.txt -o out`
  → 摘要打印 `strategy two_pass_reparse`，`--timing` 只有 `phase.pass1_ms`，exit 0（实际单趟）。
- **根因**：可重读能力是 driver 记住的路径，不是输入对象的属性，所以策略选择方无法询问"你能重读吗"。
- **建议归属**：input 层拆 `StreamedInput` / `ReplayableInput::replay(regions)`；
  `ExecStrategy` 的选择函数接收 replay 能力（或 `Option<&ReplayableInput>`）；
  `StdinSpill` 变成 `StreamedInput::into_replayable(tmp_dir)`。
- **如何验证**：`cat x.vcf | reed filter - --strategy two_pass_reparse` 必须报 usage（或如实降级并打印）；
  `strace` 统计 `-r` 时输入只 `open` 1 次（现状 2 次）。

### V30-R06（P1，`repro`）同一能力缺失有三种错误类别与文案：stdin exit 2 / 缺索引 exit 1 / plain+gzip 变成 "malformed BGZF" exit 1

- **证据**
  - `crates/reed/src/input.rs:100-104`（stdin+regions → `IoError::UnsupportedInput`）、
    `:141-153`（缺 `.csi/.tbi` 同一变体）、`crates/reed/src/error.rs:77-80`（`Io` → `ErrorKind::Runtime`）。
  - `crates/reed/src/input.rs:110-123`：region 分支检出 `Format` 后仍直接 `BgzfChunkSource::open`，
    不对"regions 只能用于 BGZF/BCF"做构造期检查。
  - stdin 的同类用法错误却走 `Error::usage`：`crates/reed/src/filter/mod.rs:267-271`、
    `crates/reed-cli/src/translate/common.rs:376-380`。
- **实测（同一能力，三种结果）**
  - `reed filter small.vcf -r chr1:10-20` → `malformed BGZF data: not a BGZF block header`，**exit 1**；
  - `reed filter small.vcf.gz -r chr1:10-20`（无 sidecar）→ `unsupported input: regions need an index next to the input`，**exit 1**；
  - `cat small.vcf | reed filter - -r chr1:10-20` → `regions need a seekable file, not stdin`，**exit 2**。
- **根因**：`UnsupportedInput` 同时装"数据格式不支持"与"命令行用法组合不支持"，
  而 region 能力检查落在读块失败处而非构造期。
- **建议归属**：`IoError::Capability { needed: Seek|Index|Bgzf|Replay, fix: &'static str }`，
  由 `Error::kind()` 映射为 Usage；`open_regions` 在构造 `RegionChunkSource` 前校验 `Format`。
- **如何验证**：上面三条命令的 exit code 与文案进 `crates/reed-cli/tests/cli_contract.rs`。

### V30-R07（P1，`code`）同一个 source 有三种获取方式，BCF/框架能力靠 `Option` 字段在运行期推断

- **证据**
  - `crates/reed/src/input.rs:62-67` 的 `source_mut` 文档称"仅本 crate 集成测试用、非稳定面"，
    生产代码在用：`crates/reed/src/genotype.rs:465`、`crates/reed/src/rename/mod.rs:751`。
  - 第二种是直接访问 `pub(crate)` 字段：`stats/driver.rs:449-463`、`stats/writers/reports.rs:1441-1475`、
    `convert/job.rs:1904-1905`、`filter/mod.rs:231`、`filter/pass0/matrix.rs:1096`、`index.rs:624-639`。
  - 第三种是按路径重开（R05）。BCF 能力用 `Option<Arc<BcfDict>>` 表达，`.bcf.is_some()` 共 **5 处**：
    `bcf_input.rs:279`、`index.rs:560,583`、`pca/mod.rs:169`、`filter/mod.rs:321`。
  - 第四种嗅探：`crates/reed/src/convert/hapmap.rs:507-546` 自写 stdin/文件/gzip 探测并重开。
- **根因**：`ChunkSource` 只暴露 `next_chunk`，格式/框架/原生解码能力留在 `OpenedInput` 的可选字段里。
- **建议归属**：`InputMeta { format, framing, bcf_dict }` 由 `ChunkSource`/`OpenedInput::into_parts`
  唯一带出；删除 `source_mut`。
- **如何验证**：`rg -n "bcf\.is_some\(\)" crates` 与 `rg -n "source_mut\(\)" crates` 收敛为 0/1；
  单测断言 BCF 夹具 `format == Format::Bcf`，而不是由字典是否存在推断。

## §14.3 生命周期与资源所有权（缓冲/临时文件/线程/预算/取消）

### V30-R08（P1，`code`）prune 线程在失败/取消路径被 detach：无 Drop、无取消点、与 stats 的 Abort 是两条通道

- **证据**
  - `crates/reed/src/filter/pass0/matrix.rs:1027`：`self.prune_jobs.push(std::thread::spawn(move || ...))`；
    句柄只存字段（`:773`），只有 `finish` 的 drain 分支 join（`:810-820`）；
    `rg "impl Drop for MatrixSink"` 无结果——`JoinHandle` 被 drop 即 detach，线程各自持有 `SiteMatrix`。
  - `crates/reed/src/filter/pass0/prune.rs:128` 的 `prune_one` 不接收取消信号，站点循环无检查。
  - 对照做对的两处：`crates/reed/src/stats/sink.rs:545`（`LdExecutor::drop` 先 trigger 再关池）、
    `crates/reed/src/rename/dedup.rs:107`（`Pool::drop` join）。
- **根因**：线程生命周期挂在 sink 字段上，字段 drop 只 detach；取消通道按 stage 各建一套。
- **建议归属**：`PrunePool`（仿 `dedup::Pool`）唯一拥有线程与取消标志，`Drop` 先置标志再 join；
  `prune_one` 接收 `&AbortSignal`。
- **如何验证**：故障注入（第 3 个 contig 后返回 Err）后断言 prune 线程计数归零；当前会失败。

### V30-R09（P1，`code`）内存预算无单一所有者：`--memory` 不参与策略选择，convert spill 用编译期常量

- **证据**
  - `crates/reed/src/filter/mod.rs:357-361` 用 `MemBudget::default_half_memory().bytes` 调
    `choose_for_mind`，**忽略** `opts.mem_bytes`；`:445-447`、`:450-455`、
    `filter/pass0/matrix.rs:1052-1057` 各自再推导一次；`default_half_memory` 在 `crates/` 内共 **8 处**调用。
  - `crates/reed/src/convert/job.rs:452` 比对常量 `crates/reed/src/convert/spill.rs:41
    SPILL_TOTAL_BUDGET: usize = 8 << 20`，而 `ConvertSink` 构造（`convert/job.rs:1885-1905`）不接收 budget。
  - LD 侧 `crates/reed/src/stats/sink.rs:167-173` 累积 `ContigMatrix`，只受 `max_ld_jobs`
    （线程数，`stats/driver.rs:406`）限制，无字节记账。
- **根因**：`MemBudget` 不是唯一所有者——各 stage 从宿主内存重算一次，convert 用常量替代预算。
- **建议归属**：`ThreadBudget`（含 `MemBudget`）单向注入每个 sink/决策点；convert 的 8 MiB 作为下限而非上限。
- **如何验证**：`--memory 64 --mind` 跑大输入，报告 strategy 必须反映 64 MiB；`rg -n "default_half_memory"
  crates/reed/src` 收敛为构造处 1 次。

### V30-R10（P1，`code`）`TempSpill` 只拥有路径：append/read 按名重开，崩溃后无清扫

- **证据**
  - `crates/reed/src/spill.rs:117-152`：`TempSpill { path: PathBuf }`，`create` 返回 `(File, Self)`
    把写句柄交给调用者，`Drop` 只 unlink 名字。
  - `crates/reed/src/convert/spill.rs:52-55`：后续 `OpenOptions::new().append(true).open(spill.path())`
    按名重开，无 `create_new`/`O_NOFOLLOW`/inode 校验；`:92`、`:140` 同样按名 `File::open`。
  - 对照 `crates/reed-io/src/write/atomic.rs:88-145`（存活表 + 自有 pid 扫除 + 其他进程残留告警）；
    spill 侧没有对应清扫，SIGKILL/OOM 后 `reed-stdin-<pid>-*.vcf`、`reed-convert-*.spill` 永久残留。
- **根因**：临时文件身份用可变路径表示，"创建时的独占"与"之后的读写"不是同一个对象。
- **建议归属**：`TempSpill` 唯一拥有 `File`（或 inode + 私有目录），读写经 `try_clone`；
  注册进与 `AtomicFile` 相同的存活表与启动扫除。
- **如何验证**：单测——`create` 后 `remove_file(path)` 并在同址 symlink `victim`，再 append，
  断言 `victim` 未被改写（当前会被改写）；子进程 `_exit(0)` 后孤儿 spill 应被扫除。

### V30-R11（P1，`code`）输出文件双 owner：`adopting` 可覆盖 `create` 自建的 `AtomicFile`，写入被静默丢弃

- **证据**
  - `crates/reed-io/src/write/ordered_output.rs:88-95`：`adopting(mut self, atomic) { self.atomic = Some(atomic); self }`
    ——无 `assert`；`:112-126`：`create` 内部已 `AtomicFile::create` 并让 `BufWriter` 绑定它的句柄。
  - 两种装配并存：`crates/reed/src/rename/mod.rs:698-705` 走 `create`；
    `crates/reed/src/filter/mod.rs:719-724`+`:180-186`、`crates/reed/src/convert/hapmap.rs:586-609` 手工配对。
- **根因**：输出 owner 既可能由 sink 自建、也可能由调用方交付，类型上不可区分，`adopting` 允许叠加。
- **建议归属**：拆 `create`（自建自持）与 `from_parts(AtomicFile, writer)`（只接外部 owner）；
  `adopting` 至少 `assert!(self.atomic.is_none())`。
- **如何验证**：单测 `OrderedOutput::create(p).adopting(a)` 后写 header，断言 `p` 含写入内容；
  当前 `p` 只含 `a` 的内容（写入随被 unlink 的 temp 丢失）。

## §14.4 变更传播（加一种输入/改一条规则要动几层）

**共变基线（主代理实测）**：`git log --no-merges` 中触碰 `crates/` 的提交 **211** 个，其中
fix/feat/refactor **162** 个，**132** 个同时改 ≥3 个 `.rs`。crate 对共现前五：
`reed×reed-cli` **57**、`reed×reed-filter` **47**、`reed×reed-io` **41**、`reed-cli×reed-filter` **31**、
`reed×reed-stats` **28**；文件对共现前三：`reed-filter/src/plan/mod.rs × reed/tests/e2e_filter.rs` **23**、
`reed-cli/src/translate/filter.rs × reed-filter/src/plan/mod.rs` **22**、
`reed-filter/src/eval/head.rs × reed-filter/src/plan/mod.rs` **20**（含回归测试，属预期；但规则本体跨 4–5 层不是）。

### V30-R12（P0，`commit`）分层倒置：L3 策略 crate 依赖 L1 传输 crate，让格式/INFO/区域语义可在两层任意落位

- **证据**：`crates/reed-filter/Cargo.toml:12` `reed-io.workspace = true`，6 个文件 **23** 处 `reed_io::` 引用；
  而 `docs/spec/architecture.md:1843` 写明"L3 策略 crate 只能依赖 reed-core、serde、thiserror、smallvec、bitflags（保持纯）"
  （`§14.6`，同页 `:621` 再次引用）。
- **影响**：这是 R01（ALT/GT）、R02（AC/AN）、R04（区域）在同一条规则上有 4–5 份实现的共同原因——
  L1 语义可以在 L1 与 L3 之间任意落位，于是"改一条规则要动几层"变成常态。
- **建议归属**：把"格式/标签/区域"的语义收回 L0/L1；`reed-filter` 只保留策略并删除 `reed-io` 依赖
  （或改成 L1 提供 trait、L3 只消费）。
- **如何验证**：架构门禁——`crates/reed-filter/src` 内 `reed_io::` 出现次数为 0；
  或 `cargo tree -p reed-filter` 不含 `reed-io`。

### V30-R13（P1，`commit`）INFO 标签没有注册表：一个 `--fill-tags` 标签有 6 份手工枚举

- **证据（git）**：`39ba997`（`--fill-tags` 校验）6 个文件；`d610f88`（GT-less 记录保留 AC/AN）6 个文件；
  `7dd8381`（AC 求和溢出 + INFO 约定）6 个文件。
- **现状 6 份枚举**：`crates/reed-filter/src/spec/recompute.rs:60 SUPPORTED_FILL_TAGS: [&str; 5]`、
  `crates/reed-filter/src/plan/compile/finish.rs:74` 的 `RecomputeFlags { wants("AC"), … }`、
  `crates/reed-filter/src/plan/mod.rs:293` 的 5 个 bool、
  `crates/reed-filter/src/eval/edit_plan.rs:29-56` 的 `fill_info_set` 逐 tag if、
  `crates/reed/src/filter/mod.rs:1029` 的 `recompute_definitions`（拼 `##INFO` 行）、
  `crates/reed-cli/src/cli/mod.rs:375-379` 的帮助文本。
  且 `docs/spec/cli.md:142` 已写 AF/F_MISSING"待实现"——下一次加标签是必然事件。
- **根因**：没有标签注册表；"名字 → header 定义 → 取值函数 → CLI 文案"四者之间没有编译期/测试期约束。
- **建议归属**：`reed-core` 一张 `RecomputeTag { name, ty, number, description, compute: fn(&SiteCounts) -> InfoValue }`
  静态表；`reed-filter` 只消费，`reed-io` 只按表渲染 `##INFO`，CLI 帮助由表生成；`RecomputeFlags` 退化为 `Vec<&'static RecomputeTag>`。
- **如何验证**："5 个标签 × 6 个枚举点"矩阵检查进 `ci/check-structure.sh`（只加一处枚举即 DRIFT）；
  `cli.md` 的标签集合 == 表。

### V30-R14（P1，`commit`）没有"幸存记录决策"值：每条过滤规则要在 3 个 `EditPlan` 组装点重复落地

- **证据（git）**：`d610f88` 的规则本体只有一行（`!has_gt` 时不更新 AC/AN），却要改
  `crates/reed/src/filter/mid_stage.rs:30` 的 `SiteEmit`、`filter/mod.rs:535` 的 `auto_stale_tags`、
  `filter/emit.rs:174` 与 `filter/single_pass.rs:361` 两处 `EditPlan` 组装、
  `crates/reed-filter/src/eval/edit_plan.rs:88` 的 `stale_for_record`、`crates/reed-filter/src/plan/mod.rs:360`。
- **现状**：`rg -n 'EditPlan \{' crates/reed/src` = **3** 处（`filter/emit.rs:174`、`filter/single_pass.rs:361`、
  `crates/reed/src/filtered.rs:166`）；两条执行路径（SinglePass / SketchThenEmit）各自把 counts+mask+stale 拼一遍。
- **根因**：没有单一的"幸存记录决策"值，规则必须在每条执行路径上各写一次。
- **建议归属**：`reed-filter` 产出纯值 `SurvivorEdit { info_set, stale, mask_gt, keep }`，
  两个 job 只负责交给 `edit_vcf_line`/原生编码器。
- **如何验证**：`EditPlan {` 组装点收敛为 1；已有的
  `single_pass_and_two_pass_reparse_produce_identical_output`（`crates/reed/tests/e2e_filter.rs:213`）
  固定为契约，再补 `FilteredInput` 角色的三角测试。

## §15 根因与已记症状的关系

- **根因解释症状**：R04 ↔ V30-Q05（区域类型三份、互斥字段运行期报错）；R07 ↔ V30-Q02（不可命名类型
  `BcfDict` 进入公开签名）；R12 解释 V30-Q01/Q02 为什么"门面一开就是 410 个下层项"——
  下层 crate 的边界本身没有被守住。
- **根因不在 Q 系列**：R08/R09/R10/R11 属运行期/资源维度（线程、预算、临时文件、输出 owner），
  Q 系列只审了 API 形状，没有覆盖；R13/R14 属"规则没有归属"的传播维度。
- **两条建议归属落在同一处**：R01/R02/R03/R04 都指向"把语义收回 `reed-core`（类型 + 唯一函数）"，
  与 Q-A1（公开面收口）共用同一批改动。

## §16 本部分边界

1. **未装工具**：`rust-analyzer`（本机只有 rustup 代理，无组件）、`cargo-modules`、`cargo-public-api`、
   `cargo-semver-checks`、`cargo-check-external-types`、`cargo-hack`、`cargo-fuzz`、`cargo-deny`
   均未安装，本轮用 `rg` + `codegraph`（索引 274 文件 / 6187 节点）+ `cargo metadata` 替代；
   **公共 API 快照与 semver 基线仍未生成**。
2. **CodeScene change-coupling 未使用**：共变统计由主代理用 `git log --no-merges --name-only` 手工统计，
   已排除合并提交，但**未排除批量改名/格式化提交**（可能高估）。
3. **端到端未复跑**：R03 的 BCF 哨兵夹具、R04 的 FASTA `chr1:100` 用例由子代理实测，
   主代理只复核了代码路径（标 `code`）；R08/R10/R11 的故障注入测试尚未编写（结论是"当前会失败"，
   为 `code` 推理）。
4. **不裁决**：每条"建议归属"是选项，不是决定；是否重构 `reed-filter` 的依赖方向（R12）需要单独拍板。


---

# 第四部分：工具化全库审计（V30-T01…V30-T20）

> 2026-11-11 追加。**六个子代理**按维度覆盖全部 11 个 workspace crate，主代理复核关键数字后再记录。
> 与前三部分的区别：这一部分以**工具输出**为主证据，不再只靠 `rg`/阅读；同时修正前文若干计数。

## §17 工具矩阵与全库基线

| 工具 | 版本 | 实际产出 | 结论 |
|---|---|---|---|
| rust-analyzer | 1.99.1 | `analysis-stats .`：109,819 LOC、`panics: 0`、MIR failed 7 (0%)；`diagnostics .`：**620 条 WeakWarning**（618 `inactive_code`） | 可用；CLI 无 find-references，引用分析用 `rg`/codegraph 补 |
| cargo-modules | 0.28.0 | 9 个 crate 的 structure/dependencies/orphans；**orphans 全为 0**；模块级环 **2 处** | 可用；`reed-cli`/`oracle` 未跑 structure |
| cargo-public-api | 0.52.0 | 9 个 crate 快照；`reed -s` **1120 行**、32 条 `pub use`、**241 条声明项 + 8 个 `pub mod`** | 可用；**不展开整模块再导出的内容** |
| cargo-semver-checks | 0.51.0 | `--workspace` 只检查 **3/11**，其余打印 `Skipping` | 可用；覆盖是假绿（T01） |
| cargo-check-external-types | 0.5.0 | 默认 nightly 报 JSON v61≠v57；固定 `nightly-2026-03-20` 后 9 个 crate 全跑通 | 可用，**必须固定工具链**（T06） |
| cargo-hack | 0.6.46 | `--each-feature` 14 plans / 1 失败；`--all-features` 11/11 通过；`--no-default-features` 唯一失败 | 可用（T07） |

补充命令：`cargo tree -d`（宿主 **7 组/14 条目**；`--target all` **10 组/20 条目**）、
`cargo package`（3 个"可发布" crate 全失败）、`cargo doc -D warnings`（**exit 101，5 crate/19 处**）、
`cargo test --doc --workspace`（**6 个 doctest，3 个实际执行**）、`cargo clippy --workspace --all-targets`
（**exit 0 / 0 warning**，但只覆盖继承 workspace lints 的 8 个 crate）。

**逐 crate 基线**（`cargo public-api -s` 的声明项 = `pub fn/struct/enum/trait/const/type`，含 trait-impl 方法；edition/rust-version/lints 来自 `cargo metadata`）：

| crate | 版本 | edition/MSRV | lints | publish | 声明项 | doctest | rustdoc -D |
|---|---|---|---|---|---|---|---|
| reed（门面） | 0.4.0 | 2024/1.91 | 继承 | false | 241 | 4 | pass |
| reed-cli | 0.4.0 | 2024/1.91 | 继承 | false | 无 lib | — | **FAIL** |
| reed-core | 0.4.0 | 2024/1.91 | 继承 | false | 969 | 1 | pass |
| reed-io | 0.4.0 | 2024/1.91 | 继承 | false | 657 | 1 | **FAIL** |
| reed-engine | 0.4.0 | 2024/1.91 | 继承 | false | 165 | 0 | pass |
| reed-stats | 0.4.0 | 2024/1.91 | 继承 | false | 782 | 0 | pass |
| reed-filter | 0.4.0 | 2024/1.91 | 继承 | false | 494 | 0 | **FAIL** |
| reed-graph | 0.4.0 | 2024/1.91 | 继承 | **可发布** | 60 | 0 | **FAIL** |
| reed-distance | 0.1.0 | 2021/**无** | **无** | **可发布** | 170 | 0 | **FAIL** |
| reed-reduce | 0.1.0 | 2021/**无** | **无** | **可发布** | 283 | 0 | pass |
| reed-oracle | 0.2.0 | 2021/无 | 仅 `[lints.rust]` | false | 313 | 0 | pass |

## §18 发现

### 18.1 公共 API 与 semver

#### V30-T01（P0，`repro`）semver 门只覆盖 3/11，且恰好是最没准备好的三个

- **证据**：`cargo semver-checks check-release --workspace --baseline-rev HEAD -v` 输出
  `Skipping reed / reed-core / reed-engine / reed-filter / reed-io / reed-oracle / reed-stats`
  （7 个 `publish = false`），`reed-cli` 无 lib target；实际只检查 `reed-distance`/`reed-graph`/`reed-reduce`
  —— 三个 0.1.0/2021/无 lints 的 crate。显式 `-p reed --baseline-rev HEAD` 则可跑通（202 checks 全过）。
  更早基线 `v0.3.0`（4491b21）因当时 `crates/` 布局尚不存在而 `failed to retrieve local crate data`，**不可比**。
- **问题**：门的"绿"是幸存者偏差：真正的 API 面（门面 + L0–L5 内核）因 `publish=false` 整片漏检；
  发布前没有可回退的契约基线。
- **建议**：CI 用逐 crate 的 `-p` 清单跑 semver-checks（或先裁决发布策略）；发布点打 tag 作为基线。
- **验证**：CI 日志中 `Skipping` 为 0；每个公开 crate 至少 202 checks 有结果。

#### V30-T02（P0，`repro`）公开签名仍含不可命名类型，且新增 `EngineError`

- **证据**：
  - `crates/reed/src/stats/request.rs:264` `pub outputs: Vec<reed_stats::request::StatOutput>`；
    `cargo public-api -p reed -s` 中无 `reed::StatOutput`。
  - `crates/reed/src/input.rs:58` `bcf_dict() -> Option<&Arc<reed_io::bcf::dict::BcfDict>>`；门面不导出 `reed_io::bcf`。
  - **新**：4 个公开写出函数返回 `reed_engine::error::EngineError`：
    `crates/reed/src/stats/writers/fstats.rs:22`、`ld.rs:38`、`sfs.rs:110/176`；门面未导出 `EngineError`。
  - `QcReportAcc` 已补导出（`crates/reed/src/lib.rs:134`），但其字段类型仍未导出——半修复。
  - 现有 `crates/reed/tests/api_surface.rs` 只校验手挑的名字，测不到这些泄漏；**12 个**门面测试
    仍直接 `use reed_core/io/engine/stats/filter::`（`e2e_filter.rs` 79 处、`e2e_bcf.rs` 57 处）。
- **问题**：调用方读得到字段/返回值，却写不出签名；API 的第一个消费者（自家测试）在绕门面。
- **建议**：导出 `StatOutput`、`EngineError`（或统一 `reed::Error`）；`bcf_dict` 删除/opaque 化；
  把"公共签名中的类型必须可在 `reed::*` 命名"做成编译门（遍历 rustdoc JSON 或 trybuild）。
- **验证**：`cargo public-api` 快照中不存在无法从 `reed::*` 命名的路径；12 个旁路测试降为 0。

#### V30-T03（P0，`repro`）门面整树再导出的规模被低估：20 条模块再导出，不是 11 条

- **证据**：`crates/reed/src/lib.rs` 有 **30** 条源 `pub use`，展开为 **32** 条 `pub use reed::X`；
  其中整 crate/整模块别名 **13** 条（`distance`、`graph`、`reduce`、`bgzf`、`chunk`、`atomic_file`、
  `header`、`budget`、`metrics`、`expr`、`funnel`、`plan`、`spec`）+ 7 个 `writers::*` 子模块 =
  **20 条模块再导出**。`reed::reduce` 全仓 **0** 引用；`reed::bgzf`、`reed::header` 也 0 直接引用
  （内部走 `reed_io::bgzf` 21 次 / `reed_io::header` 14 次）。`cargo public-api` 只列 `pub use` 语句、
  不展开被再导出 crate 的内容，故快照无法反映真实可达面；`ci/` 下没有 `api/` 快照。
- **问题**：前文 Q01 的"11 条"与"≈410 个下层项"都不准确；真实可达面是 8 个 crate 的公开项总和
  （8 个被再导出 crate 的声明项合计 **≈3,580**），且没有机器门约束。
- **建议**：逐项再导出；删除 `distance/graph/reduce` 别名或先补 `run_distance`；生成并入库
  `ci/api/*.txt`，CI diff；用 `cargo public-api --diff` 检查跨提交变化。
- **验证**：`grep -E '^pub use .* as [a-z_]+;' crates/reed/src/lib.rs` 为 0；公共 API 快照入库。

### 18.2 发布与外部契约

#### V30-T04（P1，`repro`）三个"可发布" crate 全部 `cargo package` 失败，元数据 11/11 缺

- **证据**（`cargo package --manifest-path crates/<p>/Cargo.toml --no-verify --allow-dirty`）：
  - `reed-distance` / `reed-graph`：`no matching package named 'reed-core' found`（`reed-core` 是
    `publish=false`，registry 中不存在；依赖本身带 `version="0.4.0"`）。
  - `reed-reduce`：`all dependencies must have a version requirement specified ... dependency 'reed-distance' does not specify a version`。
  - `cargo metadata --no-deps`：11 个 crate 的 `repository/documentation/keywords/categories/homepage`
    **全缺**；`readme` 仅 `reed-reduce` 有。
- **问题**：三个默认可发布的 crate 一个也发不了；元数据不满足 crates.io 的展示/检索要求。
- **建议**：先裁决发布策略（发布闭包 vs 仅 git 依赖）；补齐元数据与 crate 级 README/LICENSE；
  修 `reed-reduce → reed-distance` 的 `version`；CI 加 `cargo package --dry-run`。
- **验证**：三个 crate 的 `cargo package --dry-run` 全绿。

#### V30-T05（P1，`repro`）第三方类型进入公共契约：noodles 10 项、serde 37 类型、bitflags 6、memmap2 1、crossbeam 2

- **证据**（`cargo +nightly-2026-03-20 check-external-types --output-format markdown-table`，9 份表）：
  - `reed-io`：`FileIndex::Csi(Box<noodles_csi::Index>)`、`FileIndex::Tbi(...)`（`chunk/region.rs:383/385`）；
    `Indexer` 的 5 个字段直接是 noodles 的深层泛型（`write/index.rs:95-110`）；
    `ChunkPayload::Mapped { map: Arc<memmap2::Mmap> }`（`chunk/mod.rs:37`）；
    `impl From<libdeflater::DecompressionError> for IoError`（`error.rs:70`）。
  - `reed-core`：`bitflags` 生成的 `SiteFlags::iter/iter_names/FLAGS/IntoIter` 与 2 个 trait impl（`site.rs:72`）；
    `Provenance::with_spec<S: Serialize>` 返回 `Result<Self, serde_json::Error>`（`provenance.rs:48`）。
  - `reed-engine`：`AbortSignal::raced_by_abort(&Receiver<()>)` / `receiver()`（`abort.rs:75/84`）。
  - `reed-filter/reed-stats`：27 + 9 个公开类型实现 `Serialize`（其中 22 个也实现 `Deserialize`）。
  - **门面可达**：`pub use reed_io::{bgzf, chunk}` 使 `Box<noodles_csi::Index>` 与
    `Arc<memmap2::Mmap>` 经 `reed::chunk` 可达，但门面不 re-export `noodles`/`memmap2`——下游拿到类型却无法命名。
- **问题**：noodles 0.119（0.x）与 memmap2/libdeflater 的内部结构被钉进 semver；serde derive 是 37 个类型的
  公开合约；`serde_json::Error` 出现在返回类型。
- **建议**：`FileIndex`/`Indexer` 改 opaque 包装；`ChunkPayload` 的 map 私有 + 访问器；
  若要保留 noodles 作为 public dependency，则显式 `pub use noodles;` 并写进发布说明；
  serde 放 feature 或明确 allow-list（注意真实路径是 `serde_core::*`）。
- **验证**：`cargo check-external-types` 的 allow-list 只剩显式批准项；门面编译探针能在不添加
  `noodles`/`memmap2` 依赖的情况下命名所有公开类型。

#### V30-T06（P2，`repro`）external-types 工具需要固定 `nightly-2026-03-20`

- **证据**：默认 nightly（1.99，rustdoc JSON **v61**）→ 工具报 `requires format version 57`；
  `rustup toolchain install nightly-2026-03-20` 后 9 个 crate 全部跑通（`/tmp/audit-a3/et-*.md`）。
  工具 0.5.0 已是最新（`cargo search --registry crates-io`）。
- **问题**：不固定工具链，CI 会静默失去这道门；工具按 crate 分析，门面的 module re-export 传递项
  需要"逐依赖 crate 跑 + 编译探针"才能覆盖。
- **建议**：CI 固定 `nightly-2026-03-20`（或升级到支持 v61 的构建）；门面额外做 union/探针。
- **验证**：CI 中该命令 exit 0，且 9 份表入库。

#### V30-T07（P1，`repro`）feature 矩阵名存实亡：只有 `reed-reduce` 有 feature，唯一失败点是 E0063

- **证据**：`cargo metadata --no-deps` 显示仅 `reed-reduce` 有 `default=["faer-backend"]`，其余 10 个
  crate `features=[]`。`cargo hack check --workspace --each-feature` 14 plans / **1 fail**：
  `reed-reduce --no-default-features` → `lanczos.rs:146 E0063 missing field 'parallel'` +
  `backend.rs:155/206 unused 'par'`。`--all-features` 11/11 通过；`--feature-powerset --depth 2` 同样
  只在该点失败。门面无条件依赖 `reed-reduce`（默认 feature），下游无法关闭 faer。
- **问题**：`reed-reduce` 文档承诺的 "dependency-free fallback" 不可编译；其余 crate 的 feature
  矩阵只是重复默认构建。
- **建议**：修 `Operator` 初始化与未用参数；CI 加 `cargo hack check --workspace --each-feature`
  与 `--no-default-features`；若要让 serde/noodles 可选，再引入 feature（否则删掉暗示）。
- **验证**：上述两条命令 exit 0。

#### V30-T08（P1，`code`）未使用的直接依赖扩大构建图与 MSRV 暴露面

- **证据**（逐 crate `rg` 源码 0 命中）：`reed` 的 `bstr` 完全未用、`sha2` 仅测试用（应移 dev-deps）；
  `reed-io` 的 `bitflags`、`reed-engine` 的 `tracing`、`reed-filter` 的 `sha2`、
  `reed-graph`/`reed-stats` 的 `thiserror` 均无源码引用；dev-deps：`reed-core` 的 `reed-oracle`、
  `reed-io` 的 `crossbeam-channel` 未用。
- **问题**：库用户被迫编译未用依赖；MSRV 与供应链面被无谓放大；与 `deny.toml` 的
  `multiple-versions = "warn"` 互相放大噪声。
- **建议**：删除或移到 dev-deps；用 `cargo machete`/`cargo udeps` 复核后加 CI 检查。
- **验证**：`cargo machete` 0 命中；`cargo tree -e normal -p reed` 不含 sha2/bstr。

### 18.3 结构与分层

#### V30-T09（P1，`commit`）分层违规复核成立：`reed-filter → reed-io` 是唯一真实越层，另有 2 处模块环

- **证据**：`reed-filter/Cargo.toml:12` 依赖 `reed-io`；源码 6 文件 **23** 处 `reed_io::`
  （`error.rs:6` 公开 `impl From<IoError> for FilterError`、`eval/head.rs:6-8`、`expr/mod.rs:23,74,99,119`、
  `eval/edit_plan.rs:11`、`plan/mod.rs:321`、`semantics.rs:167,196`），与第三部分 R12 完全一致。
  `architecture.md:1843`（§14.6）白名单只允许 `reed-core/serde/thiserror/smallvec/bitflags`，
  `memchr/serde_json/sha2/tracing` 也超白名单。
  模块级 2 段聚合后全仓只有 2 个环：`reed-engine: abort ↔ budget`、`reed-filter: expr ↔ semantics`
  （源码双向 `use` 已核）。`reed-engine → reed-io` 属 DAG 规定方向，`rayon` 只受 L3 约束，均**不**违规。
- **问题**：L3 无法脱离 L1 单测；reed-io 的接口成为 filter 的准 ABI；两处同层纠缠使模块不可独立替换。
- **建议**：把 filter 需要的最小传输视图（列偏移/INFO 片段/`IoError` 映射）抽象到 `reed-core`
  trait/值类型；`AbortSignal` 与 `CmpOp` 下沉到叶子模块；删除 `sha2`/`serde_json` 的非测试用途。
- **验证**：`cargo tree -p reed-filter` 不含 `reed-io`；模块级 SCC = 0。

#### V30-T10（P1，`commit`）分层规格与实际 crate 清单脱节：`reed-qc` 不存在，4 个 crate 未分配层级

- **证据**：`cargo metadata` 的 11 个成员中，`architecture.md` §3.1/§4/§20 D2 只定义 8 个
  （core/io/engine/stats/**qc**/filter/reed/cli）；实际无 `reed-qc`（QC 在 `reed-stats/src/summary/*`
  与 `reed/src/stats/*`），且 `reed-distance`/`reed-reduce`/`reed-graph`/`reed-oracle` 未被 §4 收录、
  未分配 Lx。它们的依赖边（如 `reed-reduce → reed-distance`、`reed → distance/reduce/graph`）
  既非合规也非违规，无法被 §3.2/§14.6 检查。
- **问题**：规格漂移会让人误以为"L3 纯净"仍成立（实际已破）；4/11 的架构位置没有裁决依据。
- **建议**：补齐 §3/§4 的层级归属（distance/reduce/graph 建议 L3，oracle 标注"迁移期测试 oracle"），
  记录 `reed-qc` 合并入 `reed-stats` 的决定。
- **验证**：`cargo metadata` 的每个成员都能在 architecture.md 找到 Lx 归属。

#### V30-T11（P1，`code`）结构卫生：`convert/job.rs` 2105 行神模块，28 个 src 文件 >800 行

- **证据**：`crates/reed/src/convert/job.rs` **2105** 行、无 `#[cfg(test)]`，同时容纳
  `ConvertJob/ConvertSink/ConvertPartial`、`run_convert/_opened` 与 5 种输出格式；
  `crates/*/src` 中 >800 行的文件 **28** 个（次大：`reed-io/src/bcf/encode_native.rs` 1851、
  `reed/src/stats/writers/reports.rs` 1526、`reed/src/rename/mod.rs` 1458）。
  正面：9 个 crate `cargo modules orphans` 全 0，`find` 交叉验证 180 个 src 文件无孤儿，
  `foo.rs` + `foo/` 无混用，模块树最大深度 3。
- **问题**：单文件承载多概念，改动面/合并冲突/审查成本高（与 R14 变更传播同源）。
- **建议**：按格式拆成 `convert/{phylo,nexus,count012,hapmap,msa,job}.rs`；其余超限文件逐步治理。
- **验证**：`wc -l crates/reed/src/convert/job.rs` < 800；结构守卫加"单文件上限"。

### 18.4 文档与可验证性

#### V30-T12（P1，`repro`）rustdoc 门红得比前文更广：5 个 crate / 19 处

- **证据**：`RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --workspace` **exit 101**；
  失败 crate：`reed-cli`（9 处）、`reed-io`（5）、`reed-distance`（3）、`reed-filter`（1）、`reed-graph`（1）。
  类型：4 处链接私有项（`MAX_DEPTH`、`format_fields`、`admit_tag`×2）、9 处未解析链接
  （`SiteSelector`×2、`0,1`×2、`covariance`、`InfoTagOverflow`、`IoError::MalformedRecord`、
  `FilterSpec::multiallelics_allowed` 等）、6 处未闭合 HTML（`<TAG>`×2、`<TAB>`、`<suffix>`、`<m>`、`<contig>`）。
  另有 `output filename collision`（lib `reed` 与 bin `reed` 同名）。
- **问题**：库文档无法用 `-D warnings` 构建；CI 无 doc 门。
- **建议**：逐条修复 19 处；CI 加 `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --workspace`；
  消除 `reed` 同名 target 冲突（bin 改名或 `doc = false`）。
- **验证**：该命令 exit 0。

#### V30-T13（P1，`repro`）文档示例近零，`missing_docs` 是假绿：强制后 252 项未文档化

- **证据**：`cargo test --doc --workspace` 共 **6** 个 doctest（reed 4：2 执行 + 2 `no_run`；
  reed-core 1；reed-io 1 compile-only），其余 7 个 lib crate 为 0；`# Examples` 全仓 **2** 个。
  普通 `cargo doc` 的 `missing documentation` 为 0，但强制 `-W missing_docs` 后：
  `reed-distance` **3**、`reed-reduce` **32**、`reed-oracle` **217**——因为 distance/reduce
  根本没有 `[lints] workspace = true`，oracle 显式 allow。
- **问题**：`missing_docs = "warn"` 在最需要约束的公开 crate 上不生效；用户的 rustdoc 自学路径几乎不存在。
- **建议**：distance/reduce 继承 workspace lints；公开 crate `missing_docs = "deny"`；
  每个 lib crate 至少 1 个 crate 级可运行示例；`cargo test --doc --workspace` 进门。
- **验证**：`-W missing_docs` 为 0；每 crate doctest ≥1。

#### V30-T14（P1，`doc`）library.md 的 6 个入口 vs 真实 14 个驱动；examples/ 全是内部基准

- **证据**：`docsite/src/content/docs/library.md:16-19`（及 zh-cn）写 "six entry points"，
  实际裸驱动 **8** 个 + `_opened` 6 个 = **14**；`:22` 称签名统一，反例是 `run_filter`（无 budget）、
  `run_index`（`&Path`/`WriteIndex`/`usize`）、`run_rename`（多 `output`）；`:73-80` 的表漏
  `Pca/Dapc` 的 Request/Report。`examples/` 仅 6 个文件，全在 `reed-reduce/`，自述为
  timing/conformance/crosscheck/benchmark，且 0 个使用门面 `reed`。README 无库用法
  （`rg 'run_|cargo add' README.md` = 0）。
- **问题**：README 承诺"Rust API"手册，但用户按文档只会看到 6 个入口、错误的统一签名和不完整的表；
  `examples/` 暗示有用户示例，实际是内部 PCA 基准。
- **建议**：由 `lib.rs` 的导出表生成/校对入口表；README 加最小 `run_popgen` 示例；
  `crates/reed/examples/` 增加只依赖门面的示例；把 6 个 harness 移到 `bench/` 或改名；
  用 `cargo test --doc` 或 mdbook test 把 library.md 片段纳入门。
- **验证**：文档入口数 == `pub fn run_*` 数；README 示例可编译。

#### V30-T15（P2，`code`）公开文档被内部票号污染，0 个 `# Stability`，三处口径矛盾

- **证据**：`crates/*/src` 中 `doc_vNN` 出现 **1,235** 次（更宽的 V/A/D 票号合计约 1,400）；
  渲染后的公开 rustdoc 约 1,366 次，门面 landing page 35 次。`rg '# Stability' cryptes` = 0。
  稳定性口径三处不同：`docsite/.../library.md:12`（只承诺 crate root）、
  `docs/spec/cli.md:26`（只承诺 CLI 与输出格式）、`docs/spec/naming.md:9-12`（同 cli）。
  同页称下层 crate 是 internal，但门面用 13 条整 crate/模块别名把它们暴露；
  `docs/spec/cli.md:26` 称所有 `reed-*` 都 `publish=false`（实际 distance/graph/reduce 未设）；
  `crates/oracle/src/lib.rs:1-3` 仍写 "`reed-core`: canonical ..." 并指向不存在的路径。
- **问题**：外部读者无法解析内部票号；没有机器可读的稳定性声明；三份文档对"稳定什么"给出三种答案。
- **建议**：票号移入 `# Implementation notes`；门面根加 `# Stability` 并统一三处口径；
  修正 oracle 文档与 publish 表述。
- **验证**：公开 rustdoc 的 `doc_vNN` 为 0；`# Stability` 存在且三处一致。

### 18.5 门与运行时

#### V30-T16（P1，`repro`）CI 门错位，5 类检查未进门；release.yml 无测试

- **证据**：`.github/workflows/ci.yml:11-24`：`build` 条件 `(base_ref||ref_name)=='main'`
  只跑 `cargo build --release`；`gates` 条件 `!= 'main'` 才跑 test/clippy/doctest/doc-drift。
  因此 push main 与所有 PR→main **不跑测试**。`gates` 的 doctest 只跑 `-p reed`；
  `cargo fmt --check`、`cargo doc`、`cargo deny`、`cargo hack --no-default-features`、
  `ci/check-structure.sh` 均无 workflow 调用（`check-structure.sh` 16 条守卫本身 exit 0，但没人跑）。
  `data-gates` 需 PR 标签 + self-hosted；`gold/popgen/` 的 21 个对拍脚本在 workflow 中零引用。
- **问题**：门与风险成反比；Q35/T12 的 rustdoc 破坏能长期存在正是因为无 doc 门。
- **建议**：合并 main/non-main 条件（main 跑全量）；`gates` 增加 fmt/doc/deny/no-default/structure；
  `cargo test --doc --workspace` 替代 `-p reed`；把 gold popgen 对拍接入（或夜间）。
- **验证**：main 的 push/PR 也出现 test/clippy/doc job；`rg 'check-structure' .github` 命中。

#### V30-T17（P1，`repro`）lints 只覆盖 8/11：distance/reduce/oracle 逃逸，强制后 28 条会失败

- **证据**：只有 8 个 crate 有 `[lints] workspace = true`；`reed-distance`/`reed-reduce` 无任何
  `[lints]`，`reed-oracle` 只有 `[lints.rust] missing_docs = "allow"`。因此
  `unsafe_code=deny`/`unwrap_used=deny`/`missing_docs` 对它们无效。
  `cargo clippy -p reed-reduce --lib -- -W clippy::all -W clippy::pedantic` 报 **201** 条，
  扣除 workspace allow 后仍有 **28** 条会触发 CI 的 `-D warnings`（manual_midpoint 12、
  explicit_iter_loop 9、if_not_else 3 等）。生产路径有 1 个 `unwrap()`（`umap.rs:110`）；
  `expect` 的 "invariant:" 约定无机器检查，已有 1 处违约（`backend.rs:251`）。
  unsafe 实为 **5 处操作**（1 个 `unsafe fn` + 4 个 block），全在 `reed-io/src/chunk/bgzf.rs`，
  SAFETY 注释齐全；其中 3 个 block 在测试内，`MmapChunkSource::open` 无生产调用者。
- **问题**：门本身有效（对 8 个 crate clippy exit 0），但覆盖面与策略文档不一致，
  "workspace 全绿"会被误读为"全部受约束"。
- **建议**：distance/reduce 补 `[lints] workspace = true` 并先修 28 条；oracle 显式写
  `unsafe_code = "deny"` 或在结构守卫中列入白名单；把 `expect("invariant: …")` 做成 grep 守卫。
- **验证**：`cargo clippy -p reed-reduce --lib -- -D warnings` exit 0。

#### V30-T18（P2，`repro`）rust-analyzer 620 条 WeakWarning：618 条 `inactive_code`，fallback 从未被编译

- **证据**：`rust-analyzer diagnostics .` → 620 条，全部 WeakWarning：`inactive_code` 618
  （其中 src 内 144 条 "test is disabled"、5 条 `feature = "faer-backend" is enabled`）、
  `remove-unnecessary-else` 1（`hapmap.rs:393`）；0 error。被标 inactive 的正是
  `lanczos.rs:142-146/195-218`、`backend.rs:165-174/215-223`、`pca.rs:112-116` 等 fallback 分支。
- **问题**：这些代码既无编译门也无 RA 覆盖（"写了但永不验证"），与 T07 的 E0063 同源。
- **建议**：修复 fallback 并加 `--no-default-features` 门；其余 inactive 分支要么删除，
  要么用 feature/测试覆盖并加注释。
- **验证**：`--no-default-features` 编译通过后，`inactive_code` 中 faer 相关项减少到 0。

#### V30-T19（P2，`code`）R08 补充：batch prune 是第二处取消缺口

- **证据**：`crates/reed/src/filter/pass0/prune.rs:191 prune_matrices(...)` 同样只收
  `&ThreadBudget`、无取消参数；调用点 `crates/reed/src/filter/mod.rs:603`
  （`--mind` + `--indep-pairwise` 的 staged 路径）。流式版 `prune_one`（`prune.rs:128`）已由 R08 记录。
  反证可注入：`matrix.rs:1060` 已创建 `abort` 并传给 `run(...)`，只是 sink 字面量没有 abort 字段。
- **问题**：R08 的"prune 线程 detach"范围应扩大为**两个 prune 引擎**；批量路径同样无法停止。
- **建议**：两个入口都接 `&AbortSignal`，由 `PrunePool::drop` 先置标志再 join。
- **验证**：故障注入后，两条路径的线程计数都归零。

#### V30-T20（P2，`code`）R10/R11 补充：`AtomicFile::commit` 自身按名重开；`reed-distance` 公开 writer 绕过原子发布

- **证据**：
  - `crates/reed-io/src/write/atomic.rs:289`：`commit` 里 `OpenOptions::new().write(true).open(&self.temp)`
    ——参考实现本身也按名重开；create 与 commit 之间若名字被换，`sync_all` 作用在受害者上。
  - `crates/reed-distance/src/output.rs:55-66 write_matrix` 与 `:151-157 write_ids` 是 `pub`，
    经 `pub use reed_distance as distance` 暴露，内部直接 `File::create(path)`（非原子、失败即截断）；
    CLI 正确走 `write_matrix_to`，这两个函数目前只有 crate 内测试调用——"死但公开"的危险 API。
  - `adopting` 当前没有 `create + adopting` 组合（`filter/mod.rs:186`、`convert/hapmap.rs:609` 都是
    `new(...).adopting`；`rename/mod.rs:698` 用 `create`），故 R11 是**潜伏**缺陷，不是现行 bug。
- **问题**：R10 推荐照抄的"正确实现"自身也有按名重开；门面暴露了绕过原子发布的路径式 writer。
- **建议**：`AtomicFile` 持有所属 `File`（或 commit 前 fstat 校验 ino/dev）；
  `reed-distance` 的两个 pub writer 删除、`#[deprecated]` 或改为调用 `write_matrix_to`。
- **验证**：`rg 'OpenOptions.*open\(&self\.temp\)' crates` 为 0；`rg 'File::create' crates/reed-distance/src` 为 0。

## §19 对前三部分记录的修正

| 既有条目 | 原记录 | 本轮实测/工具结果 | 处理 |
|---|---|---|---|
| V30-Q01 | 11 条整模块别名；≈410 个下层 pub 项 | **20 条模块再导出**（13 crate/模块 + 7 writers）；"≈410"不可复现且偏低 | 以 T03 为准 |
| V30-Q02 | `StatOutput`/`BcfDict`/`QcReportAcc` 不可命名；13 个测试绕过门面 | `QcReportAcc` 已导出（半修复）；**新增 `EngineError` 出现在 4 个公开返回值**；旁路测试 **12** 个、`e2e_filter` **79** 处 | 以 T02 为准 |
| V30-Q12 | 4 处 `#[non_exhaustive]`，都在错误枚举 | 4 处，但 1 处是 `site::VariantClass`（非错误） | 修正表述 |
| V30-Q24 | 3 个可运行 doctest | **6** 个 doctest（3 执行 + 3 compile-only）；`# Examples` 2 | 以 T13 为准 |
| V30-Q29 | distance/reduce/graph 因 "path 依赖缺 version" 打包失败 | distance/graph 是 **registry 无 `reed-core`**（它 `publish=false`）；只有 reduce 是缺 version；readme 仅 reduce 有 | 以 T04 为准 |
| V30-Q30 | `reed-graph` 无 `rust-version` | **有**（`crates/reed-graph/Cargo.toml:6`，1.91）；distance/reduce 才是没有 | **删除该错误** |
| V30-Q35 | 3 个 crate / 9 处 rustdoc 错误 | **5 个 crate / 19 处**（新增 reed-cli 9 处、reed-io 5 处等） | 以 T12 为准 |
| V30-Q36 | 同名多版本 7 组 | 宿主图 **7 组/14 条目**；`--target all` **10 组/20 条目**（+r-efi、thiserror、thiserror-impl） | 修正计数口径 |
| V30-R12 | `reed-filter → reed-io` 6 文件 23 处 | 复核一致；另发现 `memchr/serde_json/sha2/tracing` 超 L3 白名单 | 以 T09 为准 |
| V30-R08 | prune 流式线程 detach | 另有 **batch `prune_matrices`** 第二处取消缺口 | 以 T19 为准 |
| V30-R10/R11 | TempSpill 只拥有路径；adopting 可覆盖 | 补充：`AtomicFile::commit` 自身按名重开；`reed-distance` 公开 writer 绕过原子发布；`create+adopting` 组合**当前不存在**（潜伏） | 以 T20 为准 |

## §20 第四部分边界

1. **主代理亲自复核**：T01 的 `Skipping` 输出、T02/T03 的 `pub use`/返回值/导出行、T04 的 package 报错、
   T07 的 hack 失败、T09 的 6 文件 23 处、T11 的 2105 行/28 个超限文件、T12 的 19 处/5 crate、
   T13 的 6 个 doctest、T17 的 5 处 unsafe、T20 的 commit 按名重开。
2. **子代理实测、主代理未逐条复跑**：T05 的 9 份 external-types 表（工具链为 `nightly-2026-03-20`）、
   T08 的未用依赖清单、T10 的 architecture 对照、T14 的文档入口表、T15 的票号渲染计数、
   T16 的 CI 条件、T18 的 RA 620 条诊断、T19 的批量 prune 调用点。
3. **工具边界**：`cargo public-api` 不展开整模块再导出，故 T03 的"可达面"是下界；
   `cargo semver-checks` 对 `publish=false` 的 crate 需显式 `-p`；
   `cargo check-external-types` 需要固定旧 nightly；
   `rust-analyzer diagnostics` 只给 WeakWarning，不能替代 clippy/rustc 门。
4. **未做**：没有修改任何代码或 CI；没有生成 `ci/api/` 快照入库；没有在 CI 上实跑新门；
   `reed-cli`/`oracle` 的 `cargo modules structure` 未跑（仅有独立孤儿核对）。
