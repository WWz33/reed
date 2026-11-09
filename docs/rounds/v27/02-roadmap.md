# 02 v27 路线与验收门

> **2026-11-09：第 0–11 步已落地。** 第 1 步的 V27-05 与方法论同源的 V27-20 已在第三批补齐；
> 逐条证据见 [STATUS.md](STATUS.md)。

## 1. 顺序

| 步 | 条目 | 前置 | 验收 |
|---|---|---|---|
| 0 | **V27-01** 裸路由 `-o` 不得等于输入 ✅ 已落地 | 无 | `-o <输入>` 时 exit 2 且输入 md5 不变；窗口路由行为不变；帮助文本区分「文件路径」与「前缀」 |
| 1 | **V27-02 + V27-05** 输出会话（唯一设计，两条一起做）——V27-02 / V27-05 ✅ 已落地 | 无 | 引入显式 `commit()` 阶段：输入完整性校验通过后才发布。截断 BCF + 已存在好 `.ld` → 失败后旧文件 md5 不变；零配对/空结果 → 本轮承诺的输出全部被替换（表头/零谱/NA 行），旧前缀文件不残留。覆盖 `.ld`、`.ld.decay`、`_sfs/`、`.fstats` |
| 2 | **V27-08 + V27-09** SFS 文件名 ✅ 已落地 | 无 | 组名中的路径成分被**消毒**（`[A-Za-z0-9_.-]` 白名单，`/`、`..` 等映射为 `_`；同 stem 追加 `-2`/`-3`），不拒绝——因为群体名是合法输入，只是不能当路径用。popg 与 `convert --sfs` 同一行为；一维谱与联合谱不再互相覆盖（本次运行内重复路径报 `Job` 错） |
| 3 | **V27-14 + V27-15** 区域谓词 AND ✅ 已落地 | 无 | `--sites-file` 与 `-t/-T` 同时给出 = 交集；`-r/-R` 与 `-t/-T` = 交集，与 `reed filter`/`bcftools view` 同结果；popg 与 stats 共用同一实现；只有单一谓词时行为逐字节不变 |
| 4 | **V27-06 + V27-07** 坐标游标与可选窗口计划 ✅ 已落地 | 1 | `--step > --window` 跨 chunk 不再误报；有序位点落在窗口空隙合法；仅 `--sfs`/`--ld-decay`/`--f2` 时无 contig 长度也能跑；输入顺序校验保留（乱序仍报错） |
| 5 | **V27-04** LD 单元素尾块 ✅ 已落地 | 1 | 4095/4096/4097/4098/8193 位点的配对数正确；`.ld` 行数 == `.ld.decay` 的 `n_pairs` 之和 |
| 6 | **V27-03** Hudson/PBS 分母 ✅ 已落地 | 无 | 双位点窗口 `Some(1/3)`；低覆盖单点 `None`；`fst_multi` 测试同步改写；低覆盖夹具与 pixy 对齐到 1e-9；`no_snps` 口径一并确定（V27-27） |
| 7 | **V27-10 + V27-12 + V27-13** variant-only 与 `no_sites`（V27-10/12/13 ✅ 已落地） | 6 | 空窗口在 variant-only 下为 0（非 NA）；π/dxy/θ 的 `no_sites` 在同一窗口一致；与 vcftools 或文档公式的对齐方向在 spec 中写明；`gold/popgen/run_pixy_bed.sh` 的 `no_sites` 列通过 |
| 8 | **V27-11** contig/genome `bp` 语义 ✅ 已落地 | 7 | 同一数据只改 `--step` 时 genome `bp` 与 π 不变；`--bed-file` 超出 contig 长度时被裁剪或拒绝；窗口和 = contig 行的不变量在 `--step > --window`、`--skip-empty`、BED 下也成立 |
| 9 | **V27-16 + V27-17 + V27-37** LD 执行所有权 ✅ 已落地 | 无 | 线程峰值 ≤ `--threads`（8/16/56 三档实测）；`--threads 1` 只 1 线程；库调用返回 `Err` 后 500 ms 内线程归零；线程创建失败返回 typed error 而非 panic |
| 10 | **V27-18 + V27-19** 内存预算可见性 ✅ 已落地 | 9 | `--window 1` 的峰值 RSS 不再随窗口数线性增长；`--bed-file` 的内存不随 chunk 段数放大；`--memory` 的保证范围写进 spec，或明确它不覆盖这些状态 |
| 11 | **V27-21…V27-36** 其余 P2（V27-20/21/22/23/24/25/26/28/29/30/31/32/33/34/36 已落地；V27-35/38 第十批落地） | 0–10 | 每条按 01-review 的「最小处理方向」逐项验收；V27-34 只作为分配优化，不作为运行时缺陷声称；V27-36 为库 API 防御性修复，不声称 CLI 可触发 |

## 2. 实施纪律

1. **先写能失败的测试，再改共享底层**。v20/v26 的做法（R1–R14、V26-01…）已证明有效；不要在多个调用点重复打 guard。
2. **V27-02 与 V27-05 必须一起设计**。只修 `.ld` 会在 `_sfs/`、`.fstats` 留下同类缺陷；只修「零配对」会漏掉「残缺输入」这条更严重的路径。输出会话的语义应是：运行开始按请求确定输出集合 → 全部校验通过 → 统一发布；而不是每个 writer 各自决定。
3. **V27-14/V27-15 改一处，popg 与 stats 同时受益**，验收必须两个命令都跑。
4. **V27-08/V27-09 的修复点在 `reed/src/stats/writers/sfs.rs`**，不是 CLI——`convert --sfs` 共用同一 writer，实测同缺陷。
5. **数值类条目（V27-03/10/12/13/27）必须带手工可算的最小夹具**（2–4 个样本、1–2 个位点），并与参照实现（pixy 2.2.4 / scikit-allel 1.3.13 / vcftools）比对到 `gold` 门要求的精度。
6. **性能类条目（V27-16…V27-19、V27-30…V27-35）先测基线再改**，报告里区分「诊断值」与「优化前后对照」。v27 立项时的数字都是诊断值。
7. **改既有测试要写清理由**。V27-03 会推翻 D15-11 的部分结论，V27-05 已把 `streamed_pairs_file_is_not_created_without_pairs` 改为 `streamed_pairs_file_is_header_only_without_pairs`，V27-25 会推翻 `cli_contract.rs:719` 的 exit code 断言——每处都要在提交信息里引用本文件对应 ID。

## 3. 本轮不作为缺陷处理

- **V27-36**（`ThetaAcc` merge 形状）：常规 CLI 未复现 panic，不把此项升级为端到端正确性缺陷；但合法的零尾短表示不应被当作不同统计请求，已于第九批完成库 API 防御性修复。只豁免 θ 的零尾，不放宽全局 `merge_zip_eq`；不同 group/n_chrom、compact 和内部非法形状仍拒绝。
- **V27-34**（jackknife 每 push 分配）：分配斜率真实，但端到端 <2%。只作为分配优化，不在正确性条目里混级。
- **`--threads 1` 多一个线程（V27-37）**：字面成立，但不是预算突破，也不是平方问题。与 V27-16 一起顺手修，不单独列为资源 bug。
- **多等位位点的 π 公式**：`04-algorithms.md:51` 写了多等位推广 `comps_i − Σ_a C(k_ia, 2)`，但驱动在 `n_alt > 1` 时整条跳过（v5 裁决「不做多等位」）。这是**文档与实现的口径差**，不是数值错误；若要修，应与「是否支持多等位」的功能决策一起做。
- **`reed target/` 产物陈旧**：仓库 `.cargo/config.toml` 指定 `target-dir=/tmp/reed-target`，仓库内 `target/release/reed`（2026-11-03）早于源码（2026-11-06）。这是环境事实，不是代码缺陷；但复现时若不注意会得到旧行为。

## 4. 需要用户拍板的选择（含建议，不预先选定）

### 4.1 V27-03：Hudson/PBS 的位点剔除规则

**冲突**：`docs/rounds/v2/04-algorithms.md:56` 与 pixy `calc.py:374` 都要求「任一群体 n < 2 时分子分母**同时**剔除」；但现有实现（`fst.rs:213`）是「n ≥ 1 进分母、n ≥ 2 才进分子」，且测试 `fst_multi::hudson_zero_numerator_is_zero_not_na` 把后者钉成预期。D15-11 修的是「单点全 NaN 输出 0.0 而非 NA」——那个结论对单点成立，但顺带把窗口级语义固化了。

- **A（建议）**：按文档与 pixy 修——`n < 2` 的位点 num/den 同时剔除。代价：改写 D15-11 的测试（拆成「单点 0.0」+「双位点 1/3」两条）；低覆盖数据上的 Hudson Fst/PBS 数值会变（方向：不再被稀释）。
- **B**：保留现状，改文档——把 04-algorithms.md:56 改成「n ≥ 1 进分母」，并删除与 pixy 对齐的承诺。代价：与 pixy 的偏差保留（低覆盖夹具 −39.9%）。

### 4.2 V27-12：`--variant-only` 的 `|S|` 定义

**冲突**：`04-algorithms.md:51` 定义 `|S|` = 「窗口内文件中出现的位点集合」（现状符合），`cli.md:274` 承诺 `--variant-only` **== vcftools `--window-pi`**（vcftools 用「有调用的位点数」）。两者不可同时成立。

- **A（曾被建议，已被 4bis.6 取代）**：保 vcftools 对拍——`absent` 只扣除「该群体至少一条被调用等位」的位置；多等位 skipped 仍单独计 present。代价：改 `PiAcc`/`DxyAcc` 的审计语义与 `no_sites` 列（V27-13），并加一个 `n_occupied` 计数器。
- **B**：保现有公式——删掉 `cli.md:274` 的「== vcftools」承诺，在 `04-algorithms.md` 写明与 vcftools 的差异及偏差方向。

### 4.3 V27-11：contig/genome 行的 `bp`

**现状**：`bp` = 「从 1 到最后一个位点所在窗口末尾」的覆盖跨度，而它被当作 variant-only 的全基因组分母 `L`。所以 `--step` 越大、`L` 越小、π 越高（实测 +33%）；`--bed-file` 超出 contig 长度时反向稀释。

- **A（建议）**：`bp` 改用声明的 contig 长度；窗口裁剪只影响窗口行。代价：改变现有 genome 行的 `window_pos_2`，若有下游依赖需同步。
- **B**：保留「覆盖跨度」语义，但改列名（如 `covered_bp`）并在日志/文档声明它不是全基因组分母；variant-only 分母另用 contig 长度。

### 4.4 V27-01 的 `-o` 语义（次要，但需确认）

裸路由的 `-o` 现在是**文件路径**（写一个 TSV），窗口路由的 `-o` 是 **PREFIX**。这个不对称本身是既成设计，问题只在于裸路由没做「不得等于输入」的检查。

- **A（建议）**：保留不对称，只加检查 + 帮助文本说明。
- **B**：统一为 PREFIX 语义（裸路由写 `PREFIX.tsv`），但这会破坏现有脚本，属破坏性变更。

## 4bis. 已拍板（2026-11-08）

**用户裁决：pixy 是主旋律。** 凡 `docs/spec/cli.md` 与 pixy 不一致的地方，以 pixy 为准；`cli.md` 的既有承诺（如「已对拍 pixy」）是本次裁决的依据而非待改对象。

裁决同时确定了下述三项，并附带一条新发现。

### 4bis.1 V27-03 定案：跟 pixy，改实现

**关键事实（本轮实测，推翻了 01-review 的「文档 vs 实现」二分）**：`n=1` 位点上 **scikit-allel 与 pixy 本身互相矛盾**。

| 参照 | 同一夹具的窗口 Hudson Fst | 规则 |
|---|---:|---|
| scikit-allel `hudson_fst` | 0.22222222222222224 | `n=1` → num=NaN、den 保留；`nansum` 得 0.2222 |
| pixy `_hudson_fst` | 0.33333333333333337 | `den = where(isnan(num), nan, h_between)` → 分母剔除 |
| reed 现状 | 0.22222222222222224 | 与 allel 一致 |

夹具（2 位点，A 组 2 个体、B 组 2 个体，site1 的 B 组仅 `0/.` 一条调用等位）：

```
chr1 10 GT 0/1 0/0 0/. ./.
chr1 20 GT 0/0 0/0 0/1 0/1
pixy: 0.33333333333333   allel/reed: 0.22222222222222224
```

**定案**：按 pixy 修 `reed-stats/src/fst.rs:213`——`mpd_ab` 门槛改为 `n_a >= 2 && n_b >= 2`，num/den 同时进出。`pbs.rs:119` 同构代码一并改（保持「同一统计量一套实现」）。

**连带后果（必须一并处理，否则会造出两条互相冲突的门）**：
- `gold/popgen/compare_pbs.py` 对拍的是 **allelogram `allel.pbs`**（`compare_pbs.py:4-6` 明写）。改完规则后，PBS 在含 `n<2` 位点的夹具上与 allel 不再一致。
- 现有 PBS 夹具（`e2e_pixy.rs:286` 的 `PBS_FIXTURE`）实测每个位点每组都是 `n=8`，**不含 `n<2`**，所以该测试不会红。
- 结论：PBS 的 gold 门需改为「在不含 `n<2` 的夹具上对拍 allel」（现状夹具即满足），或改对拍 `pixy` 的三对 `fst`。**建议前者**：PBS 是 Yi et al. 2010 的定义，其参照实现就是 `allel.pbs`；只要夹具不含 `n<2`，两个规则给出相同结果，门仍然有效。

### 4bis.2 V27-27 定案：`no_snps` 跟 pixy（只数变异位点）

实测（1 条变异 + 2 条不变，`--fst hudson --bypass-invariant-check --window 100`）：

```
pixy: -0.16666666666667  no_snps=1
reed: -0.16666666666666652  no_snps=3
```

数值一致、只有 `no_snps` 不同。按 pixy 修：`FstAcc`/`PbsAcc` 只对 `all_variant` 递增 `n_sites`。

### 4bis.3 V27-12/V27-13 定案：`|S|` 与 `no_sites` 都跟 pixy（有调用的位点）

pixy 的 `no_sites` 是「该群体有等位调用的位点数」（`pixy/stats/summary.py:208`：`np.count_nonzero(np.sum(gt_pop.count_alleles(), 1))`）。实测（1 条全缺失 + 1 条变异，`--pi --bypass-invariant-check`）：

```
pixy:  avg_pi=0.5  no_sites=1  count_diffs=3  count_comparisons=6  count_missing=6
reed:  avg_pi=0.5  no_sites=2  count_diffs=3  count_comparisons=6  count_missing=6   ← 只有 no_sites 不同
```

`avg_pi` 相同是因为 `--bypass-invariant-check` 不使用 `window_bp` 补分母；**分歧只在 `--variant-only` 的 `absent` 计算**。实测全缺失记录被算作 present：`4/594`（reed）vs `4/600`（vcftools/pixy 口径）。

**定案**：
- `no_sites` 改为「该群体已调用等位数 > 0 的记录数」（pixy 口径），π/dxy/θ 三个 writer 共用同一计数器。
- `--variant-only` 的 `absent = window_bp − n_occupied`，其中 `n_occupied` 同上。
  **4bis.6 已取代这段的旧 D7-11 兼容方案**：本项目不做多等位，过滤期即剔除，不再保留 skipped-present 计数。
- 连带：`docs/rounds/v2/04-algorithms.md:51` 的 `|S|` 定义需改写为「窗口内该群体有调用的位点集合 + 多等位 present 位点」；`cli.md:274` 的「== vcftools」保留（pixy 与 vcftools 在此口径上一致）。
- 连带：`e2e_stats.rs::trailing_skipped_site_advances_the_windows` 的两条断言 `3/594`、`13/2994` 需按新分母重算（该夹具无全缺失位点，但尾窗的 `absent` 由 `n_present_file` 改为 `n_occupied`；需重新推导后更新）。

### 4bis.4 新发现：pixy 的 dxy 是**成对**的 `no_sites`，π 是逐群体的

实测（A 组 2 条全缺失、B 组全调用）：

```
pixy pi:  A ... no_sites=1 ; B ... no_sites=2
pixy dxy: A B ... no_sites=1
reed pi: A ... no_sites=2 ; B ... no_sites=2
reed dxy: A B ... no_sites=2
```

即 pixy 的 dxy `no_sites` 是 `pop1_sites AND pop2_sites`（`summary.py:316`），**不是**任一群体。这与 `|S|` 在 π 里是逐群体、在 dxy 里是成对交集的结构一致。V27-13 的修复必须区分这两者：`PiAcc` 用「该群体有调用」，`DxyAcc` 用「**两个群体都有调用**」。

### 4bis.5 V27-11 定案：接受 C（新增分母字段，输出列不变）

`bp` 保持「覆盖跨度」语义（`window_pos_2` 不变，下游脚本不受影响），另加一个**只供 variant-only 分母使用**的 `denom_bp`：contig/genome 行取声明的 contig 长度；无声明长度时退回 `span_end` 并给 warning。`--bed-file` 取「BED 区间并入后裁剪到 contig 长度」的长度。`.log` 写明分母来源。

### 4bis.6 多等位：项目不做（2026-11-08 用户裁决），popg 改为**过滤期剔除**

**用户裁决**：「本项目不考虑多等位」。据此 popg 不再让多等位记录进入统计阶段。

**现状（实测）**：`-M` 在 popg 里**是生效的**，但默认不套用文档声称的 `-M 2` 门：

```
multiallelic.vcf（1 条多等位 + 1 条双等位），--variant-only --window 100
  不写 -M   avg_pi = 3/594 = 0.0050505   ← 多等位通过过滤，在统计阶段被 run.rs:433 跳过
  -M 2      avg_pi = 3/600 = 0.005       ← 过滤阶段剔除
  -M 3      avg_pi = 3/594 = 0.0050505   ← 与默认相同
```

根因：`MultiallelicPolicy::Reject` 是 `reed-filter` 的默认值（`spec/mod.rs:369`），但 popg 的统计通路走的是 `StatsJob` 的 `stats_policies()`（`run.rs:248-253`），那里是 `CollapseAlt`——即在 genotype 层把多等位压到 ALT 位，随后由 `run.rs:433` 的 `n_alt > 1` 分支整条跳过。所以多等位**既不是被过滤掉的，也不是被统计的**，而是「进来再丢掉」。

**定案 A（用户选定）**：

1. popg 的路由在构造请求时把 `MultiallelicPolicy::Reject` 接入过滤链（与 `filter` 的默认门一致），使多等位记录在**过滤阶段**被剔除，永不进入 `StatsJob`。
2. `-M > 2` 报 usage 错（exit 2），消息说明「本项目不支持多等位」；`-M 2` 与不写等价。
3. `n_multiallelic_skipped` 从「统计期跳过计数」改为「过滤期剔除计数」，并在**窗口路由的 stderr 汇总与实际 `.log`** 中都输出（裸路由已有，窗口路由当前完全没有，见 V27-39）。

**连带删除（这是 A 的主要收益）**：多等位不再进入统计后，「skipped present」这一整套机制失去存在理由，可整体删除：

| 位置 | 删除内容 |
|---|---|
| `reed-stats/src/pi.rs` | `n_skipped_present` 字段、`note_skipped_site()`、`n_present_file()` 里的 `+ n_skipped_present` |
| `reed-stats/src/dxy.rs` | 同上（8 处） |
| `reed-stats/src/request.rs` | `StatAcc::note_skipped()` 分派 |
| `reed-stats/src/window.rs` | `ContigBins::push_skipped()`、`RegionAccs::push_skipped()` |
| `reed/src/stats/run.rs` | `ContigAcc::push_skipped()`、`SegmentPartial` 的多等位分支里的 `segment.acc.push_skipped(pos0)` |
| `reed-stats/src/pi.rs:195` | 测试 `skipped_multiallelic_is_present_for_variant_only` |

删除后 `n_present_file()` ≡ `n_occupied`（pixy 口径），V27-12 的实现变成单一计数器，不再有「两套 present 定义」的张力。**D7-11 随之作废**（它的前提是多等位会进入统计）。

**注意保留**：`run.rs` 多等位分支里的 `segment.first_pos0/last_pos0` 更新属于 D15-20（尾窗必须被跳过位点推进）。A 落地后多等位在过滤期就被剔除、不再产生 segment，这段自然失效；但**验收时必须重新验证 D15-20 的场景仍然通过**（`e2e_stats::trailing_skipped_site_advances_the_windows` 用的夹具就是一条多等位记录——它需要改成用「被 `-t` 排除的位点」或「无调用位点」来构造尾窗场景，否则会随多等位一起失效）。

**行为变更（需在 CHANGELOG 与 `cli.md` 写明）**：
- `-M 3` 及以上：从「静默跳过」变为 **usage 错**。
- 含多等位的输入：`--variant-only` 的 π/dxy 分母不再把该位点算作 present（分母变大、数值变小，向 pixy/vcftools 靠拢）。

### 4bis.6 V27-01 定案：A（保留不对称 + 加检查）

保留「裸路由 `-o` = 文件路径、窗口路由 = PREFIX」的既成不对称，在 `prelude` 之前调用 `refuse_output_eq_input`，并在帮助文本写明裸路由的 `-o` 是文件路径。

## 5. 与既有决策的关系

| 既有 ID | 本轮影响 | 说明 |
|---|---|---|
| D15-11 | **V27-03 会推翻其部分结论** | D15-11 的「单点 0.0 非 NA」保留；「分母累加而分子不累加」按 4bis.1 改为分母也剔除。 |
| D7-11 | **4bis.6 / V27-39 使其作废** | 多等位在过滤期剔除，不再进入统计；`n_present_file` ≡ `n_occupied`。4bis.3 早先的兼容文字已被 4bis.6 明确取代。 |
| D7-23 | V27-38 需重新审视 | D7-23 的前提「CLI 只产 `per_bp:false`」在裸路由不成立。 |
| V23-14 | V27-05/V27-20 已修复 | 保留流式输出与逐文件原子；改为运行开始创建临时 pair 表、成功后统一发布，并由 sink → driver → CLI 回报实际路径，不回退到内存累积。 |
| V23-61 | **V27-15 曾违反，第二批已恢复** | `SiteSelector` 将 `-r/-R` 与 `-t/-T` 分成独立族并逐族 AND；`cli.md` §7.3 已补齐统计命令的组合语义。 |
| V27-14/15 实施决议（2026-11-08） | 新增 `reed_filter::spec::region::SiteSelector`（`indexed`/`include`/`exclude`/`sites` 四族 + `allows` 逐族 AND）；`PopgenRequest`/`StatsReportRequest` 的三个扁平字段收敛为单个 `select`。 | 不做「不支持组合就报 usage 错」的 D7-6 备选——三族语义都已定义，交集才是 `cli.md` 的承诺；结构上收敛使再次被折叠成 OR 时编译期就能看见。 |
| V26-08 | **V27-02 曾违反，第一批已恢复** | pair 临时文件即使在运行开始创建，也要等完整性检查通过后才 commit；失败路径删除临时文件，旧 `.ld` 不变。 |
| D15-20 | V27-11 是其延伸 | D15-20 断言「窗口和 = contig 行」，但只在默认参数下测过；`--step > --window`、BED 下不成立。 |
| **`cli.md` §7.1 原则** | **4bis 的直接依据** | 「缺失数据统一采用 Korunes & Samuk 2021（pixy）的逐位点计数法」——pixy 主旋律是既有的规范原则，不是新引入的偏好。 |
| `gold/compare_pbs.py` | 4bis.1 需同步 | 该门对拍 `allel.pbs`；改规则后须确保夹具不含 `n<2`（现夹具 `n=8`，已满足）或改对拍 pixy 三对 Fst。 |

## 6. 实施顺序记录

**第十批（2026-11-10，已落地）**：V27-35 与 V27-38。以主流统计定义而非最小改动为准：
f-statistics 使用 ADMIXTOOLS 观察位点动态块和有效位点加权，统一共享网格，确保线程/切分不变；
真实 qpDstat 的自然估计与 SE 全精度一致。裸路由 variant-only 与 pair/bp 分派修正，
`theta_w.total=θ_raw`、`theta_w.bp=θ_raw/L`，窗口 avg 保持 pixy 的 `raw/no_sites`。
全 workspace 1126 项全绿，可用 gold 全通过，文档漂移 98 ok / 0 DRIFT。

**第九批（2026-11-09，已落地）**：V27-36。`ThetaAcc::merge` 在相同 group/n_chrom 且未 compact 时，允许省略零尾的直方图互合并；保留语义和内部形状断言，拒绝异请求与 compact 终态。失败先行夹具 `(3 vs 7)` 正反方向均与整体累加一致；正常 CLI 输出与 θ 公式不变。

**第一批（2026-11-08，已落地）**：V27-01、V27-02、V27-08、V27-09。理由是「数据安全 + 失败不覆盖 + 路径可控」，
不涉及数值口径、不需拍板，且改动面各自集中在一处（`popg.rs` 裸路由、sink 的 commit 时序、`sfs.rs` 的命名函数）。
V27-05 当时**有意推迟**：它与 V27-02 共用「运行开始确定输出集合 + 统一发布」的设计，单独修 `.ld`
会在 `_sfs/`、`.fstats` 留下同类缺陷。

**第二批（2026-11-08，已落地）**：V27-14 + V27-15。理由见 §2 纪律 3：改一处、popg 与 stats 同时受益；
可用 `bcftools view` 直接对拍。实施细节与真实数据三方对拍见 [STATUS.md](STATUS.md) 的「第二批实施记录」。

**第八批（2026-11-09，已落地）**：V27-10、V27-21…V27-34（除 V27-35/36 另列）。
补齐 levels/群体/BED/长度/退出码/日志契约、多等位升序、临时文件告警与性能项；
V27-35 的 jackknife 网格算法与 V27-38 的裸路由列语义已按科学定义在第十批实施。

**第七批（2026-11-09，已落地）**：V27-18、V27-19。流式窗口输出会话、标量快照与 raw rows；
BED 按 contig 共享累加器，跨 worker/seam 只合并一次。`--memory` 的输入在途块保证与全进程 RSS 区分写明。

**第六批（2026-11-09，已落地）**：V27-16、V27-17、V27-37。每 run 一个 LD 执行池，
pipeline 与 LD 共享总线程预算；失败时 abort 并 join，成功不触发 abort。V27-18/19 内存项另列。

**第五批（2026-11-09，已落地）**：V27-03、V27-11、V27-12、V27-13、V27-27、V27-39。
Hudson/PBS 与 occupied-site 计数跟 pixy；多等位在统计 gate 前剔除，`-M > 2` 报 usage；
`denom_bp` 将 contig/genome 分母与公开 span 解耦。实际 CLI 与 pixy 小夹具逐列对拍通过。

**第四批（2026-11-09，已落地）**：V27-04。`decay_range` 的单元素尾块不再提前返回；
以 4095/4096/4097/4098/8193 边界与 CLI `.ld` == decay `n_pairs` 总和锁定。

**第三批（2026-11-09，已落地）**：V27-05、V27-06、V27-07 与 V27-20。先分离输入游标和窗口 frontier，
再让无窗口统计走 `GridMode::None`，最后按请求初始化空结果输出；流式 pair 文件仍由 driver 在完整性校验后 commit。

**收尾顺序（源自立项时的建议，保留作记录）**：原计划「第 0 步 + 第 1 步 + 第 2 步」一次做完，
实际拆成了三批——V27-05 先与 V27-02 共用输出会话设计，再与第 4 步（V27-06/V27-07）同批落地。