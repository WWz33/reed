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

### 3.4 LD：`--indep-pairwise W[kb] [S] R2`（过滤规则，裁决 2026-10-10）与 `--ld-decay [MAX_BP]`（stats）

#### 3.4.1 `--indep-pairwise`：一条命令出剪枝后的 VCF

用户裁决：不照搬 plink 的 `.prune.in / .prune.out` 两段式，但"第一步选位点"要的就是 plink2 的那个算法。因此：

- **归属**：`SiteFilterArgs` 的位点规则（filter 与 stats 共用）。`reed filter … --indep-pairwise 50 5 0.2 -Oz -o pruned.vcf.gz` 直接输出剪枝后的文件；`reed stats … --indep-pairwise …` 只在剪枝后的位点上统计。
- **清单**：`--prune-list PREFIX` 可选，写 `PREFIX.prune.in` / `PREFIX.prune.out`（plink2 格式：ID 一行一个；ID 为 `.` 用 `chrom_pos`）。
- **阶段**：02 §6 第 6b 步——`--mind`、`--geno/--maf/--mac/--hwe` 与依赖基因型的 `-i/-e` 之后、INFO 更新之前；只看存活位点、掩码后的基因型（plink2 同序：导入期 `--vcf-min-gq` 等先掩码，`--maf/--geno` 先于剪枝）。
- **漏斗 id**：`after_mask.ld_prune`。

**算法 = plink2 `--indep-pairwise`（`plink2_ld.cc` `IndepPairwiseThread`、`LdPruneNextWindow`，2026-10-10 源码核实）：**

1. 参数：`W` 个位点，或 `W kb` 按 bp；`S` 步长（位点数）。plink2 2.0：`S` 缺省 1，`kb` 窗口时**必须**为 1（`plink2_help.cc:1204`）。reed 相同：`--indep-pairwise W[kb] [S] R2`，`S` 缺省 1，kb 窗口给 `S ≠ 1` 报错。按 contig 独立（plink2 的 subcontig）；reed 无性别信息，X/Y 不特殊处理，文档写明。
2. 每个位点进入窗口时：基因型按**主等位**计数编码 −1/0/+1（`FillVaggs`：`sum = ref2 − alt2`，`ssq = hom_ct`，`nm = 已调用数`）；窗口内**单态**位点（`FillVaggs` 返回真）立即剪掉。
3. 窗口内每一对（plink2 默认从窗口尾向前扫，"nearest new pairs first"；PLINK 1.x 顺序不做）：只在**两位点都已调用**的样本上求 `n`、`Σx`、`Σy`、`Σx²`、`Σy²`、`Σxy`（`ComputeIndepPairwiseR2Components` 的三种情形），`cov = n·Σxy − Σx·Σy`，`var_i = n·Σx_i² − (Σx_i)²`，判定 **`cov² > R2 · var1 · var2`**（严格大于：零共同观测的对不剪）。
4. **剪谁**：比较两位点的主等位频率 `maj`（全样本、导入后，plink2 `--freq` 口径；`--pop` 不参与）：`maj_first ≤ maj_second · (1 + ε)` 时剪**后者**（位置靠后的），否则剪前者。等价说法：**次等位频率更低的被剪，平局剪靠后的**（ε = plink2 的 `kSmallEpsilon`）。`--indep-preferred` 不做。
5. 推进：位点窗 `start += S`；bp 窗"start 至少前进到 end 能前进一个位点"；被剪位点不再进入后续窗口；每窗完成后清当前窗的剪除标记。
6. 产物：被剪集合 → reed pass0 的"保留位图"。

**架构（reed）**：有顺序依赖 → 两遍，复用 `--mind` 的 pass0 / sketch 机制。pass0：workers 并行解码每个位点的 2 位基因型向量与主等位频率（走 cohort 通道，不渲染文本），reducer 按 chunk 顺序喂给**每 contig 一个** `LdPruner`（contig 之间并行、contig 内串行）；窗口只需 W 个位点的向量常驻（348 样本 → 每位点 2 × 6 个 u64），不整 contig 常驻。pass1 = 现有发射路径按位图筛。与 `--mind` 同时给出时同一 pass0：先样本位图，再剪枝。`LdPruner`（纯数值）留在 `reed-stats/src/ld.rs`，pass0 驱动在 `reed/src/filter/ld_prune.rs`（门面有权同时接 filter 与 stats 层）。

**代码示意**（非最终代码；类型名取现有的 `ChunkJob`/`OrderedSink`（`reed-engine/src/pipeline.rs:49,100`）、`LdPruner::push(site, pos0, planes)`（`reed-stats/src/ld.rs:227`）、`SitePlanes` 的 4 个位平面；`MissingJob`/`SumSink`（`reed/src/filter/reduce.rs`）是同一模子的样本级版本，照抄其骨架）：

```rust
// reed/src/filter/ld_prune.rs —— pass0：只产出"保留位图"，不写任何记录。

/// worker 交给 reducer 的最小单元：一个通过了 02 §6 第 1–6 步的位点。
pub(super) struct PruneSite {
    ordinal_in_chunk: u32,          // chunk 内的存活序号；reducer 加基址得全局序号
    contig: ContigId,
    pos0: u64,
    planes: [Box<[u64]>; 4],        // 掩码后的位平面（与 sketch 的布局相同）
    maj_freq: f64,                  // 主等位频率：由平面 popcount 得 alt/called，全样本口径（plink2 --freq）
}

#[derive(Default)]
pub(super) struct LdPrunePartial {
    sites: Vec<PruneSite>,
    n_survivors: u32,               // 本 chunk 存活位点数 = 下一 chunk 的基址增量
}

impl ChunkJob for LdPruneJob {
    type Partial = LdPrunePartial;
    type Scratch = PipelineScratch;
    fn run_chunk(&self, chunk: &ChunkBytes, scratch: &mut Self::Scratch, out: &mut Self::Partial)
        -> Result<(), EngineError>
    {
        for record in self.framing.records(chunk) {
            // 与 pass1 完全相同的 head → mask → count → tail 判定，保证两遍看到同一存活序列
            let Some(scan) = self.pipeline.survivor_planes(record, scratch)? else { continue };
            out.sites.push(PruneSite {
                ordinal_in_chunk: out.n_survivors,
                contig: scan.contig, pos0: scan.pos0,
                planes: scan.planes_after_mask(),
                maj_freq: scan.major_allele_freq(),
            });
            out.n_survivors += 1;
        }
        Ok(())
    }
    // run_seam：接缝记录与 pass1 同样处理（同一函数），否则两遍的序号错位
}

/// reducer：按 chunk 顺序接收；每个 contig 一条串行流，交给各自的线程跑 LdPruner（contig 之间并行）。
pub(super) struct LdPruneSink {
    base: u64,                                   // 已接收的全局存活序号基址
    current: Option<(ContigId, Sender<PruneSite>)>,
    workers: Vec<JoinHandle<(ContigId, Vec<u64>)>>,   // 每 contig 一个：返回该 contig 的保留序号列表
    cfg: PruneConfig,                            // threshold, LdWindow::{Sites(w)|Bp(w)}, step
}

impl OrderedSink<LdPrunePartial> for LdPruneSink {
    fn accept(&mut self, _id: ChunkId, p: &mut LdPrunePartial) -> Result<(), EngineError> {
        for mut s in p.sites.drain(..) {
            let ordinal = self.base + u64::from(s.ordinal_in_chunk);
            if self.current.as_ref().map(|(c, _)| *c) != Some(s.contig) {
                self.spawn_contig(s.contig);       // 换 contig：关闭旧 sender（线程自然收尾），开新线程
            }
            s.ordinal_in_chunk = 0;                // 线程端只认全局序号
            self.current.as_ref().unwrap().1.send((ordinal, s)).map_err(|_| EngineError::Internal)?;
        }
        self.base += u64::from(p.n_survivors);
        Ok(())
    }
    fn finish(mut self) -> Result<u64, EngineError> {
        self.current.take();                        // 关最后一个 sender
        let mut keep = KeepBitmap::new(self.base);  // 全局位点序号 → 保留
        for w in self.workers { let (_, kept) = w.join()?; for o in kept { keep.set(o); } }
        self.report.n_pruned = self.base - keep.count();   // 漏斗 after_mask.ld_prune
        Ok(keep.count())
    }
}

/// 每 contig 线程：严格按到达顺序喂 LdPruner；窗口语义、剪谁、单态剔除全在 LdPruner 里（照 plink2，见上文算法 1–5）。
fn contig_pruner(rx: Receiver<(u64, PruneSite)>, cfg: PruneConfig, n_samples: usize) -> Vec<u64> {
    let mut pruner = LdPruner::new(cfg.threshold, cfg.window, cfg.step, cfg.w, n_samples);
    for (ordinal, s) in rx {                       // rx 关闭 = contig 结束
        let planes = [&s.planes[0][..], &s.planes[1][..], &s.planes[2][..], &s.planes[3][..]];
        pruner.push(ordinal as usize, s.pos0, planes);   // 现有 API；需补 maj_freq 入参以实现"剪 MAF 低者"
    }
    pruner.finish()                                // 需新增：推进到 contig 末尾，返回最终保留的全局序号
        .into_iter().map(|o| o as u64).collect()
}
```

要点：① 两遍必须用**同一个** `RecordPipeline` 判定函数产生存活序列，位图按"存活序号"索引，不按文件偏移；② `LdPruner` 现有签名缺主等位频率与 `finish()`，补这两处即可，其余照 plink2 规则实现在它内部；③ pass1 = 现有 `sketch_then_emit`/`two_pass_reparse` 的发射路径多一条 tail 规则 `after_mask.ld_prune`：`keep.get(ordinal)`；④ `--prune-list` 在 pass1 发射时顺手写两个清单（存活且保留 → `.prune.in`，存活但被剪 → `.prune.out`），不需要第三遍。

**对拍**：`gold/popgen/run_plink2_ld.sh`：fixture `ld.vcf` 与全量子集 Chr01 上 `plink2 --vcf … --indep-pairwise 50 5 0.2 --out x` 的 `x.prune.in` 与 reed `--prune-list` 的 `.prune.in` **逐 ID 相同**；再以 `--maf 0.05 --geno 0.1 --indep-pairwise 50 5 0.2`（plink2 同参数）验证顺序；两种窗口单位（`50` 与 `50kb`）各一例。门：集合差异为 0。剪枝后 VCF 的正文 = 用该清单 `--extract` 得到的正文（同一命令内自洽）。

#### 3.4.2 `--ld-decay [MAX_BP]`、`--ld-pairs`（stats）

对距离 ≤ MAX_BP（默认 `100k`）的位点对计算 r²（unphased，plink2 `--r2-unphased` 的基因型相关，缺失按两位点都已调用的样本），按距离分箱（默认 1 kb）输出 `PREFIX.ld.decay`：`bin_start bin_end n_pairs mean_r2`；`--ld-pairs` 另写逐对 `PREFIX.ld`（plink2 `.vcor` 列：`CHROM_A POS_A ID_A CHROM_B POS_B ID_B UNPHASED_R2`）。架构：距离受限的逐对统计可流式——worker 对 chunk 内的对直接算，跨 chunk 的对由"尾窗携带"（每 chunk 末尾 MAX_BP 内的位点随 partial 上交，reducer 补算跨界对）。对拍：plink2 `--r2-unphased --ld-window-kb` 在 `ld.vcf` 上逐对 1e-9。

**代码示意**（非最终代码；`pair_sums(x, y, n_samples).r2()` 与 `LdDecayAcc::{push, merge}` 是 `reed-stats/src/ld.rs` 的现有 API）：

```rust
/// 携带的位点：只在 chunk 边界附近才需要，数量 ≤ MAX_BP 内的位点数。
#[derive(Clone)]
struct SiteVec { contig: ContigId, pos0: u64, planes: [Box<[u64]>; 4] }

#[derive(Default)]
pub struct LdDecayPartial {
    acc: LdDecayAcc,          // chunk 内部、距离 ≤ MAX_BP 的对
    head: Vec<SiteVec>,       // 距本 chunk 首个位点 ≤ MAX_BP 的位点：供上一 chunk 的尾窗配对
    tail: Vec<SiteVec>,       // 距本 chunk 末位点 ≤ MAX_BP 的位点：留给下一 chunk
    pairs: Vec<PairRow>,      // --ld-pairs 时的逐对行（chunk 内）；跨界对由 reducer 追加
}

// worker：一个按距离滑动的 ring buffer，和 pass1 完全相同的存活判定。
fn run_chunk(&self, chunk: &ChunkBytes, scratch: &mut Scratch, out: &mut LdDecayPartial) -> Result<(), EngineError> {
    let mut ring: VecDeque<SiteVec> = VecDeque::new();
    let mut first_pos0: Option<(ContigId, u64)> = None;
    for record in self.framing.records(chunk) {
        let Some(site) = self.pipeline.survivor_planes(record, scratch)?.map(SiteVec::from) else { continue };
        while ring.front().is_some_and(|f| f.contig != site.contig || site.pos0 - f.pos0 > self.max_bp) {
            ring.pop_front();
        }
        for prev in &ring {
            if let Some(r2) = pair_sums(prev.planes_ref(), site.planes_ref(), self.n_samples).r2() {
                out.acc.push(site.pos0 - prev.pos0, r2);
                if self.want_pairs { out.pairs.push(PairRow::new(prev, &site, r2)); }
            }
        }
        let first = *first_pos0.get_or_insert((site.contig, site.pos0));
        if first.0 == site.contig && site.pos0 - first.1 <= self.max_bp { out.head.push(site.clone()); }
        ring.push_back(site);
    }
    out.tail = ring.into_iter().collect();        // 末尾 MAX_BP 内的位点，按位置有序
    Ok(())
}

// reducer：跨界对 = 上一 chunk 的 tail × 本 chunk 的 head（同 contig、距离 ≤ MAX_BP），然后换届。
impl OrderedSink<LdDecayPartial> for LdDecaySink {
    fn accept(&mut self, _id: ChunkId, p: &mut LdDecayPartial) -> Result<(), EngineError> {
        for a in &self.carry {
            for b in &p.head {
                if a.contig != b.contig { continue; }
                let d = b.pos0 - a.pos0;                       // head 在 tail 之后，d ≥ 0
                if d > self.max_bp { break; }                  // head 按位置有序，可提前退出
                if let Some(r2) = pair_sums(a.planes_ref(), b.planes_ref(), self.n_samples).r2() {
                    self.total.push(d, r2);
                    if self.want_pairs { self.pairs_out.write(PairRow::new(a, b, r2))?; }
                }
            }
        }
        self.total.merge(&p.acc);
        self.pairs_out.write_all(p.pairs.drain(..))?;          // chunk 内的对按序追加
        self.carry = std::mem::take(&mut p.tail);
        Ok(())
    }
    fn finish(self) -> Result<u64, EngineError> { self.writers.write_decay(&self.total) }
}
```

要点：① head/tail 的"≤ MAX_BP"都以**本 chunk 的首/末存活位点**为基准，这已足够覆盖所有跨界对（离本 chunk 首位点超过 MAX_BP 的位点离上一 chunk 任何位点都更远）；② 接缝记录归属哪个 chunk由 `run_seam` 决定，head/tail 按归属后的序列构造，否则会重复或漏掉一对；③ `--ld-pairs` 的逐对行要按位置有序写出，所以 chunk 内的行在 reducer 里追加、跨界对插在前面，顺序即"上一 chunk 的所有行 → 跨界行 → 本 chunk 的行"；④ 内存：carry 只有 MAX_BP 内的位点，与 contig 长度无关（04 §3 P5-5 的门）。

#### 3.4.3 对在途 X4.4 代码的处置

2026-10-10 00:59 之后的工作树把 `--ld-decay`、`--ld-pairs`、`--indep-pairwise` 都放在 `reed stats`（`cli.rs` stats 块）并以清单为剪枝产物。处置：`--ld-decay / --ld-pairs` 保留不动；`--indep-pairwise` 移到共用的 `SiteFilterArgs`，按 3.4.1 实现 pass0 位图 + pass1 发射，`--prune-list PREFIX` 写清单；stats 路径 = 位图筛位点后统计。`LdPruner` 的数值部分可原样复用。

---

## 4. 多等位逐等位重算（W6 / X5）——**不做**（裁决 2026-10-10；以下为存档设计）

用户裁决：**本项目目前不做多等位**。现状即终态：默认双等位门（`-M 2`）剔除多等位位点；显式放开多等位又要求 AC/MAF/HWE 重算时**报错**（不静默写错）；stats 侧 π/dxy 的多等位推广已在累加器里，不再扩展。下面的设计与代码示意**存档**，将来若放开时直接用。设计：
- `--fill-tags`：AC/AF 按 `Number=A` 逐 ALT 输出（`ScanOutput::alleles` 已保留原始等位下标）；MAF = 第二常见等位频率（`bcftools +fill-tags` 口径）；HWE 只在双等位位点算，多等位写 `.`；NS、AN、F_MISSING 与双等位相同。
- `--maf/--mac`：多等位时次等位 = 第二常见等位（同上口径）；`--geno` 不变。
- stats：π/dxy 的多等位推广已在累加器（`comps − Σ_a C(k_a,2)`）；θ 的突变数 `Σ(等位数−1)`。
- 对拍：`multiallelic.vcf` 上 `--fill-tags AC,AN,AF,MAF,HWE,NS` 与 `bcftools +fill-tags` 逐字节一致；全量 indel 集抽样 1,000 条逐字节一致。门：默认双等位门行为不变（gold md5 不变）。

**代码示意**（非最终代码；`ScanOutput::alleles` 是解码层已保留的原始等位下标，`SampleMask` 是现有的保留样本位图；落点 `reed-filter/src/spec/recompute.rs` 与 `eval/edit_plan.rs`）：

```rust
/// 多等位位点的逐等位计数（INFO Number=A 的输入）。
pub struct AlleleTally {
    pub an: u32,            // 已调用等位数（半调用计 1）
    pub ac: Vec<u32>,       // 每个 ALT 的计数，长度 = n_alt
    pub ns: u32,            // 至少一条染色体已调用的样本数
    pub n_samples: u32,     // 保留样本数（F_MISSING 的分母）
}

pub fn tally_alleles(scan: &ScanOutput, n_alt: usize, keep: &SampleMask) -> AlleleTally {
    let mut t = AlleleTally { an: 0, ac: vec![0; n_alt], ns: 0, n_samples: keep.count() };
    for (sample, gt) in scan.alleles.iter().enumerate() {      // gt: [Option<u8>; 2]，掩码后为 [None, None]
        if !keep.contains(sample) { continue; }
        let mut called = false;
        for a in gt.iter().flatten() {
            called = true;
            t.an += 1;
            if *a > 0 { t.ac[usize::from(*a) - 1] += 1; }      // 0 = REF 不进 ac
        }
        if called { t.ns += 1; }
    }
    t
}

/// bcftools +fill-tags 口径的输出值；双等位时与现有路径逐字节相同（这是门）。
pub fn fill_tag_values(t: &AlleleTally, n_alt: usize, het_hom: Option<(u32, u32, u32)>) -> Vec<(&'static str, InfoValue)> {
    let af: Vec<f32> = t.ac.iter().map(|&c| if t.an == 0 { 0.0 } else { c as f32 / t.an as f32 }).collect();
    // MAF：把 REF 也放进来排序，取第二常见等位的频率（bcftools 的定义）
    let mut freqs: Vec<f32> = std::iter::once((t.an - t.ac.iter().sum::<u32>()) as f32 / t.an.max(1) as f32)
        .chain(af.iter().copied()).collect();
    freqs.sort_by(|a, b| b.total_cmp(a));
    let maf = freqs.get(1).copied().unwrap_or(0.0);
    let mut out = vec![
        ("AC", InfoValue::Ints(t.ac.iter().map(|&c| i64::from(c)).collect())),   // Number=A
        ("AN", InfoValue::Int(i64::from(t.an))),
        ("AF", InfoValue::Floats(af)),                                              // Number=A
        ("NS", InfoValue::Int(i64::from(t.ns))),
        ("MAF", InfoValue::Float(maf)),
        ("F_MISSING", InfoValue::Float(1.0 - t.ns as f32 / t.n_samples as f32)),
    ];
    // HWE 只对双等位位点定义；多等位写 "."（bcftools 同）
    out.push(("HWE", match (n_alt, het_hom) {
        (1, Some((hom_ref, het, hom_alt))) => InfoValue::Float(hwe_exact_p(hom_ref, het, hom_alt)),
        _ => InfoValue::Missing,
    }));
    out
}

// --maf / --mac 在多等位位点上的"次等位"同样取第二常见等位（含 REF）：
// minor_count = 按计数降序排列的 [an−Σac, ac…] 的第 2 项；双等位时退化为 min(ac, an−ac)，与现有 tail 规则一致。
```

要点：① 双等位位点必须与现有路径逐字节相同（浮点格式、字段顺序、`Number=A` 的单元素向量仍写成标量形式与 bcftools 一致）；② `tally_alleles` 走的是 `ScanOutput` 已有的下标，不需要第二次解析；③ 默认双等位门（`-M 2`）不变，本改动只在用户显式放开多等位时生效。

---

## 5. gold/popgen 的缺口（接入时补齐）

现有：`run_pixy.sh`、`run_vcftools.sh`、`compare.py`、`compare_pbs.py`、`compare_reports.py`、`make_fixtures.py`。
缺：`run_easysfs.sh` + `compare_sfs.py`、`run_vcftools_fst.sh`、`compare_fstats.py`、`run_plink2_ld.sh` + `compare_ld.py`、`run_filltags.sh`（W6）。每个脚本的判据写在 §3/§4 的"门"里；通过后把结果数字登记到当轮 STATUS。
