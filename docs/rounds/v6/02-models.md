# 02 模型：数据模型、执行模型、统计模型

> "模型"在本项目里有三层：**数据模型**（基因型与计数在内存里长什么样）、**执行模型**（一条分析如何在 56 个核上展开）、**统计模型**（估计量的定义与数值行为）。本文逐层给出现状、与参照实现的比较、以及 v6 的最优方案。所有现状陈述带 `file:line`。

---

## 1. 数据模型

### 1.1 基因型的三种表示

| 表示 | 定义 | 大小（348 样本） | 生产者 | 消费者 |
|---|---|---|---|---|
| `GtCode` | 1 字节/样本：低 4 位 = 两个 slot 的状态（REF/ALT/MISSING），bit4 相位，bit5 单倍体（`reed-core/src/genotype/gt_code.rs:3-31`） | 348 B/位点 | 文本扫描器（`reed-io/src/vcf/genotype_scan.rs:146-188`）、BCF 64K 查表（`reed-io/src/bcf/gt.rs`） | filter 的逐样本计数循环（`reed/src/filter/record_pipeline.rs:519-551`）、stats 的 `add_code`（`reed/src/stats/job.rs:930-937`）、编辑器的掩码位 |
| `SitePlanes` | 4 个位平面 `slot0_alt / slot0_miss / slot1_alt / slot1_miss`，每平面 `W = ⌈N/64⌉` 字；可选相位平面（`reed-core/src/genotype/planes.rs:1-37`） | 4 × 6 × 8 = 192 B/位点 | `PlaneBlock::push_codes`（由 GtCode 转置） | sketch pass0（`reed/src/filter/sketch.rs:45-57`）、LD（`reed-stats/src/ld.rs:33-72`：成对和 = 平面 AND + popcount）、`--het` |
| `SiteGroupCounts` | 6 × u32：hom_ref/het/hom_alt/n_missing… 派生 an/ac（`reed-core/src/count/mod.rs:23-56`） | 24 B/群体 | `add_code` 或 popcount 内核 | 全部累加器（`reed-stats/src/request.rs` 的 `push(row)`） |

**与 plink2 的 2 位 `nyp` 比较**：plink2 每基因型 2 位（hom-ref/het/hom-alt/missing），相位另存，88 B/位点；类计数 = 3 次 popcount。reed 的 4 平面是它的 2.2 倍，但能表达 **半调用**（`0/.`：一个 slot 缺失）——这是规范里"缺失按等位计"（`--geno/--mind`、`--sfs` 的 m）的前提，2 位编码表达不了。结论（与引擎审计一致）：**不把 2 位编码定为规范表示**；`GtCode` 继续做解码层的唯一交换格式，`SitePlanes` 做分析层的紧凑格式。

**计数内核的现状**（核实）：`reed-core/src/count/kernel.rs` 的 AVX-512 `VPOPCNTDQ` 路径只在 8 字块上向量化（`kernel.rs:121-122` `vector_words = n_words − n_words % 8`），348 样本时 W = 6，向量循环**一次都不执行**；且工作区没有任何生产调用者（`count_site_group`/`CountTable` 只在测试里）。filter 与 stats 的计数都是逐样本分支循环：3.3M × 348 ≈ 1.15G 次 `add_code`，约 2–3 CPU-s，占 stats 的 3%。所以计数不是性能杠杆；内核的去留见 03 §5（soundness）。

### 1.2 记录模型（filter）

现状是**文本优先**：判定用 `GtCode` + 类型化 FORMAT 列（文本扫描或 BCF 原生解码），编辑在文本上做（`line_edit.rs:77-116`，逐列拷贝、只改 INFO 与被掩码样本的 GT），输出 `-Oz` 直接压文本；`-Ob` 把编辑后的文本再交给原生编码器**重新解析**（`single_pass.rs:338-345` → `native_encode.rs:131-139`）。BCF 输入时存活记录先 `render_record` 成文本再编辑（`single_pass.rs:305-326`），所以 BCF→BCF 是 解码→渲染→解析→编辑→解析→编码。代价实测：`-Ob` 比 `-Oz6` 多 61 CPU-s（04 §1）。

最优形态：一个**格式中立的记录模型**只用于需要重编码的路径（`-Ob`、BCF 输入），文本输出保持字节级编辑（它本来就是最省的）。见 03 §8。

### 1.3 计数与累加器

- 累加器是**可加 partial**：整数精确，浮点用 `NeumaierSum`（`reed-core/src/math/sum.rs`），`merge` 逐分量相加；跨线程确定性由 `determinism.rs` 锁定。
- `StatKind`/`StatAcc` 是**封闭枚举**（`reed-stats/src/request.rs`），`push(row)` 一次 `match` 分派（`request.rs` 的 `push`）。请求数 ≤ 10 时这比 trait 对象快且可内联；不需要注册表。
- 窗口分箱：`ContigBins` 用 `BTreeMap<bin, Vec<StatAcc>>`，每位点一次 `entry()`（`reed-stats/src/window.rs:73-80`），bin 首次触碰时克隆 seed；窗口收尾按 bin 顺序合并。O(log bins) 每位点，忽略不计。
- genome-wide 量（SFS、f 统计、LD）作为与窗口 partial **并列**的 partial 上交（`reed/src/stats/job.rs:1296-1340` 的 `sfs/fstats/ld_partial`），sink 收尾合并——这是 v5 定下的正确形态。

### 1.4 位点身份

`--indep-pairwise` 的 pass0 → pass1 通过 `(contig, pos0)` 的有序表二分查找传递保留集合（`reed/src/filter/ld_prune.rs:233-247`），而不是序号位图。正确性依赖 pass0 与 pass1 看到同一存活序列（同一 `RecordPipeline` 判定），这一点代码满足（`ld_prune.rs:100-144` 与 `single_pass.rs` 共用 `process`）。开销 1.77M × log₂(235k) 次比较，可忽略。

### 1.5 标识类型

`reed_core::ContigId` 根导出的是 `pub type ContigId = u32`（`reed-core/src/contig.rs:9`，`lib.rs:24`），而 `id.rs:35` 还有一个 newtype `ContigId(u32)`（`SiteKey` 用它）。两个同名类型、一个裸整数——违反 `id.rs:5-6` 自己写的"公开 API 不接受裸整数"。统一为 newtype。

---

## 2. 执行模型

### 2.1 三种执行形态

| 形态 | 适用 | 现状 | 效率（04 §3） |
|---|---|---|---|
| **流式一遍**（BAPP：块对齐 chunk → worker 融合 inflate+计算+deflate → 有序归约） | filter 单遍、stats 窗口统计、reports、index、SFS/f 统计 | `reed-engine/src/pipeline.rs:129-479`；worker 只交换 chunk 与小 partial | filter 96%，stats 76%（0.4 s 串行） |
| **两遍**（全局样本规则） | `--mind` | pass0 = sketch（每头部存活位点 `(4+R)×W` 字）或 `MissingJob` 全遍，pass1 发射（`reed/src/filter/mod.rs:235-300`）；BCF 输入强制 two-pass（估计按行数，`mod.rs:209-211`） | 两遍各 ≈ 一遍 |
| **顺序依赖 / 全局分析** | `--indep-pairwise`、`--ld-decay`（将来 PCA/GRM/IBS） | 被塞进有序归约的 sink 里**串行**（`ld_prune.rs:318-369`；`stats/job.rs` 的 `ld_tail` 跨 chunk 配对） | 剪枝 66 s、衰减 169 s（04 §6） |

前两种形态是对的，第三种是错的：有序归约是一个**单线程串行点**，只该做 O(chunks) 的拼接，不该做 O(sites × window) 的计算。

### 2.2 全局分析的矩阵化模型【规】

plink2 的答案是先把基因型矩阵放进内存再并行（它的 `--indep-pairwise` 从 pgen 只要 0.66 s）。reed 的对应方案：

1. **pass0 并行解码到按 contig 的位平面矩阵**：worker 对每个 chunk 产出一个 `PlaneBlock`（存活位点的 4 平面）+ `pos0` 数组 + 主等位频率数组；reducer 只把 `Arc<PlaneBlock>` 按顺序推进 `Vec`（O(chunks)）。全量 3.3M 位点 × 192 B = **634 MB**，按 contig 分段每段几十 MB；超过 `--memory` 预算时按 contig 逐个处理（contig 内仍是矩阵）。
2. **分析阶段并行**：LD 衰减按位置区块并行（每区块回看 ≤ MAX_BP 的前驱，只读，无依赖）；剪枝按 contig 并行（plink2 的 subcontig 并行，`plink2_ld.cc` 的 `subcontig_thread_assignments`），contig 内保持 plink2 的顺序语义；r² 用平面 AND + popcount（现有 `pair_sums`）。
3. **pass1** 按保留位图发射（现有路径）。

目标：剪枝 66 s → ≤ 3 s（CPU 已够：plink2 1.5 CPU-s），衰减 169 s → ≤ 10 s（343 CPU-s ÷ 54 = 6.4 s 下限）。

**代码示意**（非最终代码；类型名取现有的 `PlaneBlock`、`pair_sums`、`LdPruner`、`ChunkJob`/`OrderedSink`）：

```rust
// reed-core/src/genotype/matrix.rs —— 按 contig 的位平面矩阵（只读分析输入）
pub struct ContigMatrix {
    pub contig: ContigId,
    pub blocks: Vec<Arc<PlaneBlock>>,   // 每 chunk 一块，site-major 4 平面
    pub block_base: Vec<u32>,           // 前缀和：块首位点在 contig 内的序号
    pub pos0: Vec<u64>,                 // 每位点位置（contig 内升序）
    pub maj_freq: Vec<f32>,             // 剪枝用主等位频率（plink2 --freq 口径）
    pub w: usize,                       // 每平面字数
}
impl ContigMatrix {
    pub fn n_sites(&self) -> usize { self.pos0.len() }
    /// 第 i 个位点的 4 平面切片（块号二分或直接查前缀和）。
    pub fn planes(&self, i: usize) -> [&[u64]; 4] { /* block_base 查块，PlaneBlock::site_planes */ }
}

// reed/src/filter/ld_prune.rs —— pass0：worker 填块，reducer 只拼接
#[derive(Default)]
pub(super) struct MatrixPartial { blocks: Vec<(ContigId, Arc<PlaneBlock>, Vec<u64>, Vec<f32>)> }
impl ChunkJob for MatrixJob {
    type Partial = MatrixPartial;
    fn run_chunk(&self, chunk: &ChunkBytes, body: &[u8], scratch: &mut Scratch) -> Result<MatrixPartial, EngineError> {
        // 与 pass1 相同的 RecordPipeline 判定；存活位点写入 PlaneBlock（按 contig 切块）
    }
}
impl OrderedSink<MatrixPartial> for MatrixSink {
    fn accept(&mut self, _id: ChunkId, p: &mut MatrixPartial) -> Result<(), EngineError> {
        for (contig, block, pos0, maj) in p.blocks.drain(..) { self.matrices.entry(contig).or_default().push(block, pos0, maj); }
        Ok(())                              // O(blocks)，没有逐位点工作
    }
}

// 分析阶段：contig 之间并行（rayon），contig 内按 plink2 顺序
pub fn prune_all(matrices: &[ContigMatrix], cfg: &PruneConfig) -> Vec<KeepBitmap> {
    matrices.par_iter().map(|m| {
        let mut pruner = LdPruner::new(cfg.threshold, cfg.window, cfg.step, m.w, m.n_samples);
        for i in 0..m.n_sites() { pruner.push(i, m.pos0[i], m.planes(i), m.maj_freq[i]); }
        pruner.finish_contig(); pruner.keep_bitmap()
    }).collect()
}

/// LD 衰减：位置区块并行，每区块只回看 ≤ max_bp 的前驱（只读）。
pub fn ld_decay(m: &ContigMatrix, max_bp: u64, edges: &[u64], n: usize) -> LdDecayAcc {
    (0..m.n_sites()).into_par_iter().by_uniform_blocks(4096)
        .map(|range| {
            let mut acc = LdDecayAcc::new(edges.to_vec());
            for i in range {
                let mut j = i;                                  // 向前回看
                while j > 0 && m.pos0[i] - m.pos0[j - 1] <= max_bp {
                    j -= 1;
                    if let Some(r2) = pair_sums(m.planes(j), m.planes(i), n).r2() { acc.push(m.pos0[i] - m.pos0[j], r2); }
                }
            }
            acc
        })
        .reduce(|| LdDecayAcc::new(edges.to_vec()), |mut a, b| { a.merge(&b); a })
}
```

要点：① reducer 只做 O(chunks) 的拼接，串行点消失；② 并行度来自 contig 数（剪枝）与位点区块（衰减），与线程数无关的确定性由"每 contig 一个顺序 pruner"与"衰减的可交换合并"保证；③ 内存按 contig 封顶，超预算时逐 contig 流式处理并提示；④ 同一矩阵以后直接给 PCA/GRM/IBS 用（plink2 的全部矩阵类分析都是这个模型）。

### 2.3 两遍策略的空白

- BCF 输入的 `--mind` 一律 two-pass（`filter/mod.rs:209-211`），因为 sketch 的记录数估计按行；BCF 可以用 chunk 的 `blocks` 表 + 首块记录数估计，或直接用矩阵模型（2.2）替代 sketch——sketch 本质就是一个带掩码位图的平面矩阵。**建议统一**：`--mind` 与 `--indep-pairwise` 共用一个 pass0（填矩阵 + 每样本已调用计数），不再维护 `MissingJob`/`SketchJob`/`LdPruneJob` 三套 pass0。
- 顺序模式（`threads ≤ 2`）规范有、代码无（05 §4）。

### 2.4 内存模型

在途 = tickets（2 × workers）× ~1 MiB × `MEM_FACTOR 8` ≈ 864 MB（引擎审计 §2.2），与实测在途峰值 912 MB 一致；字节许可只在小内存机器上先于票据生效。矩阵模型新增 ≤ 634 MB（全基因组一次性）或每 contig ≤ ~60 MB。规则：矩阵模式的估计（位点数 × 192 B）必须进 `--memory` 预算判断，超了就按 contig 分段。

---

## 3. 统计模型

### 3.1 定义与对拍状态（全部已逐位对拍，详见 docs/rounds/v5/05 §1 与 STATUS）

| 统计量 | 定义要点 | 参照物 / 结果 |
|---|---|---|
| π、dxy | Σ差异/Σ比较（pixy 逐位点计数法）；`--variant-only` 补单态分母 | pixy 4.8e-14；vcftools `--window-pi` 1.05e-5 |
| θ_W、Tajima's D | 逐位点 n 的调和数；D 的方差按窗口平均调用数算一次 | pixy 3.1e-14 / 8.9e-11 |
| Fst wc84 / hudson（成对） | WC84 分量 Σa/Σ(a+b+c)，n_i = an_i // 2；Hudson 比值之和 | pixy 2.9e-15 / 2.2e-15 |
| Fst global | r 群体 WC84（`WcMultiAcc`） | vcftools weighted Fst 2.2e-16；逐位点 2e-6（vcftools 打印 6 位） |
| PBS | 三对独立 Hudson，T 截断 [0, 0.99999] | scikit-allel 逐列 |
| 投影 SFS | 超几何投影到 n，全支撑 k ∈ 0..=n | easySFS 1e-15（夹具）、≤ 5.7e-11（子集） |
| f2/f3/f4/D | Patterson 2012，bp 分块 jackknife | allel 估计量 1.3e-16；SE 不可比（块定义不同） |
| LD r² | unphased，两位点共同调用样本 | plink2 2e-6（打印精度）；剪枝集合全量 235,403 == plink2 |
| freq/hardy/het/missing/depth/summary | plink2 / vcftools / bcftools 布局 | 逐行一致 |

### 3.2 已登记的偏差与局限（模型层面）

| 项 | 性质 | 处理 |
|---|---|---|
| f 统计 SE：reed 按 bp 分块，allel 按变异数分块 | 定义不同，非错误 | 文档写明（已） |
| `--indep-pairwise Wkb`：区域内 reed 保留 131 vs plink2 158 | **窗口推进语义未对齐**（plink2 `LdPruneNextWindow`：bp 窗口下 start 至少前进到 end 能前进一个位点） | M6-4 |
| `--ld-decay` 需要 `--pop` | LD 与群体无关，却被当作"窗口统计"要求 `--pop`（实测报错 `a windowed statistic needs --pop FILE`） | M6-3 |
| `--sfs-int` 半整数取整与 easySFS 差 1 | 浮点噪声方向不同 | 对拍以取整前为准（已） |
| 联合 SFS、多等位 | 未接 / 裁决不做 | — |
| Tajima's D 的 n 取窗口平均 | pixy 2.x 口径，与 vcftools（2N 固定）不同 | 规范已选 pixy |

### 3.3 数值模型

f64 累加 + `NeumaierSum`；INFO 浮点以 f32 读写（与 htslib 相同）；比值统计在收尾处一次相除；分母为 0 写 `NA`。确定性：partial 按 chunk 顺序合并，浮点加法顺序固定，故与线程数无关（`determinism.rs` 在 1/2/3/7/56 线程上 sha256 相同）。

### 3.4 模型层建议（M6）

| # | 建议 | 落点 |
|---|---|---|
| M6-1 | **位平面矩阵成为全局分析的规范输入**（§2.2），`--mind`/`--indep-pairwise`/`--ld-decay` 共用一个 pass0 | 03 §7、06 Y1 |
| M6-2 | `ContigId` 统一为 newtype | 06 Y3 |
| M6-3 | `--ld-decay` 不要求 `--pop`（LD 与群体无关）；窗口统计与 genome-wide 量的 `--pop` 要求分开判定 | 06 Y2 |
| M6-4 | LD kb 窗口按 plink2 `LdPruneNextWindow` 对齐并对拍 | 06 Y2 |
| M6-5 | 表达式的 INFO 访问复用 head 阶段的片段索引（一记录一次解析），不要每个 tag 重扫（`expr/mod.rs` `info_value` vs `eval/head.rs:190-214` 的 `HeadScratch.fragments`） | 03 §4 |
| M6-6 | `GtCode`/`SitePlanes` 维持；popcount 内核修 soundness 或删除（03 §5） | 06 Y3 |
| M6-7 | 不做 2 位规范表示；若将来群体数 > 8 的计数成为瓶颈，再加派生的"2 位类平面" | 不排期 |
