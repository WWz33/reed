# 04 性能（2026-10-10 v6 实测）

> 全量 `clean.vcf.gz`（7.0 GB，3,315,604 条，348 样本）/ `clean.bcf`（8.1 GB）；二进制从 `8709f1f`（工作树干净）于 10:22 构建；56 线程；热缓存；批次开始 1 分钟负载 **1.28**（docs/rounds/v5 04 §4 的纪律：批次开始采样，逐次采样只作参考）；每个数字取 3 次中位数，线程扫描每点 1 次；计时用 bash `time`（`%R %U %S`）。参照物：plink2 v2.0.0-a.7LM（AVX2）、bcftools 1.24。

---

## 1. 全量基线（与 v5 对照）

| 命令 | v5（10-09）wall / user | **v6 wall / user** | 3 次 wall | 说明 |
|---|---|---|---|---|
| filter gold `-Oz6`，VCF 输入 | 8.79 / 461 s | **8.85 / 460 s** | 9.03 8.85 8.82 | 不回退 |
| filter gold `-Oz6`，BCF 输入 | 8.87 / 465 s | 9.00 / 464 s | 9.13 9.00 8.96 | BCF ≈ VCF（W4） |
| filter gold `-Ob` | 9.94 / 523 s | 9.98 / 521 s | 9.98 9.98 | 文本→BCF 再解析 +61 CPU-s |
| filter gold **`-Oz1`** | — | **5.64 / 289 s** | 5.65 5.63 | 6 级比 1 级多 **171 CPU-s**（37%） |
| filter 无位点规则（全部存活）`-Oz6` | — | 11.0 / 572 s | 10.98 11.04 | 输出记录数 ×1.9，deflate 随之 |
| stats 五统计（π/dxy/θ/D/Fst），VCF | 1.74 / 71 s | **1.78 / 72 s** | 1.78 1.81 1.80；复测 1.75 1.75 1.76 | |
| stats 五统计，BCF | 1.51 / 55 s | 1.65 / 56 s | 1.56 1.65 1.68 | |
| reports 两族 + summary | 3.76 / 186 s | **3.0–5.4 s（双峰）/ 147 s** | 3.13 5.02 4.55；复测 3.04 5.32 5.41 | X3 省了 39 CPU-s；wall 双峰见 §7 |
| stats 边过滤 `--pi` | 3.04 / 152 s | 3.09 / 153.6 s | 3.10 3.09 3.06 | |
| `reed index` VCF.gz | 1.60 / 75 s | 1.63 / 75 s | 1.62 1.63 | |
| `reed index` BCF | 1.05 / 44 s | 1.12 / 44.5 s | 1.11 1.13 | |
| **新**：`--sfs --proj 214` 三群体 | — | 1.60 / 73 s | 1.62 1.60 | +1–2 CPU-s |
| **新**：`--fst --fst-global` | — | 1.52 / 71 s | 1.49 1.56 | ≈ 0 |
| **新**：`--f3` + `--dstat` | — | 1.56 / 71.5 s | 1.57 1.56 | ≈ 0 |
| **新**：`--ld-decay 100k`（需 `--pop`） | — | **169 / 343 s** | 169.2 168.5 | **串行**（§6） |
| **新**：`filter --indep-pairwise 50 5 0.2 -Oz1` | — | **66.6 / 289 s** | 66.3 66.9 | pass0 64.4 s 全在 reducer（§6） |

跨格式一致：五统计 VCF 与 BCF 输入的 `_pi.txt` md5 `45c2f2e1632076fa601eb5575d6288f3`、`_fst.txt` `fa144045…`，与 v4/v5 记录相同。

## 2. gold 对拍（本轮重跑）

`gold/scripts/run_reed.sh`（8 线程，76.4 s）→ 正文 md5 **`7aedc59ac798729a2674e6e7ef19adfc`**（`bcftools view -H` 与 `zcat | grep -v '^#'` 两种算法相同）；`check_parity.py` 对 `bcftools.full.vcf.gz`：1,768,759 / 1,768,759 位点，Jaccard 1.000000，348 样本匹配，基因型不一致 mean=0 max=0，**PASS**。X4/X6 之后正文不变。

## 3. 线程扩展（并行效率）

模型：wall ≈ user ÷ (threads − 2)（reader 与 reducer 各占一线程）。

| 线程 | filter gold `-Oz6` wall | 理想 | 效率 | stats 五统计 wall | 理想 | 效率 |
|---|---|---|---|---|---|---|
| 8 | 76.3 s | 76.2 | 100% | 12.0 s | 11.7 | 98% |
| 16 | 32.8 s | 32.6 | 99% | 5.33 s | 5.0 | 94% |
| 32 | 15.7 s | 15.2 | 97% | 2.70 s | 2.36 | 87% |
| 56 | 8.80 s | 8.46 | 96% | 1.75 s | 1.33 | 76% |

user 在四个点上恒定（filter 457–459 s，stats 70.3–71.7 s）：没有随线程增长的同步开销。
结论：**filter 在 Amdahl 极限上**，唯一杠杆是每条记录的 CPU 工作量；**stats 在 56 线程有 ≈ 0.4 s 的串行部分**（reader 0.35 s 的 mmap 块头扫描 + 收尾合并 + 启动），再往下要缩的是串行段而不是 worker。

`--timing`（filter VCF）：worker_busy 465.4 s、reader 0.39 s、reduce 0.45 s、6,647 chunk、在途峰值 912 MB——reducer 只占 0.45 s，输出追加不是瓶颈。

## 4. CPU 账

### 4.1 deflate 级别（子集 Chr01-03，1.03 GB 输入，gold 规则，56 线程）

| 级别 | wall | user | 输出大小 | 相对 6 级 |
|---|---|---|---|---|
| 1 | 0.93 s | 42.1 s | 705 MB | CPU −37%，体积 +25% |
| 4 | 1.15 s | 53.6 s | 604 MB | CPU −20%，体积 +7% |
| 6（默认） | 1.40 s | 66.8 s | 564 MB | — |

全量：6 级比 1 级多 171 CPU-s。htslib 把 `-l6` 映射到 libdeflate 7 级（`bgzf.c` `lvl_map`），reed 直接传 6，已经比 bcftools 默认略快略松；输出兼容性与级别无关（读者不看压缩字节，对拍门是解压后 md5）。

### 4.2 构成（全量 filter gold，460 CPU-s）

| 项 | CPU-s | 依据 |
|---|---|---|
| deflate 6 级 | ≈ 230 | 6 级 − 1 级 = 171，1 级本身 ≈ 60 |
| 扫描 + 掩码 + 计数 + 编辑 | ≈ 140 | 对照实验（doc_v4 04 §1） |
| inflate | ≈ 60 | 同上 |
| 其中可去的"纯浪费" | ≈ 30–40 | 每块 65 KiB 的 `resize(…, 0)` memset（deflate 与 inflate 各一次，约 23 + 15 GB memset）、每块一次 `Vec` 分配、每条存活记录的 `stale` Vec / `edit_info` 两个 Vec / `format_g6` 的 String、每 worker 建一张 128 KiB GT 表、`--fill-tags HWE` 不走 `HweCache`（03 §4） |
| `-Ob` 额外 | +61 | 编辑后的文本被原生编码器**再解析一遍** |

### 4.3 stats（72 CPU-s）

inflate ≈ 60 + 每样本 `add_code` 计数 ≈ 2–3 + 分箱累加 ≈ 5；popcount 内核在生产路径没有被调用（03 §4）。新接入的四个统计量各 +0–2 CPU-s。

## 5. 与参照工具

| 任务 | reed（从 VCF.gz） | plink2 | bcftools |
|---|---|---|---|
| LD 剪枝 `--indep-pairwise 50 5 0.2` | **66.6 s**（235,403 保留） | 从 pgen **0.66 s / 1.5 CPU-s**（235,403 保留，与 reed 完全相同）；含 VCF 导入 26.6 s | — |
| VCF 导入 | （流式，无导入） | `--make-pgen` 27.8 s / 105 CPU-s | — |
| freq/hardy/het/missing | 3.0 s（双峰到 5.4）/ 147 CPU-s | 从 pgen 1.48 s / 2.75 CPU-s；从 VCF 需 +27.8 s | — |
| gold 过滤管线 | 8.85 s | — | 三段管线约 40 min（docs/rounds/v5/STATUS X2 备注，单线程） |
| 建索引 | 1.63 s | — | 142 s（docs/rounds/v3 记录） |

两点结论：① 对"流式一遍"的任务（过滤、窗口统计、报告、索引）reed 领先一到两个数量级；② 对"全局、有顺序依赖"的任务（LD 剪枝、LD 衰减）reed 现在**慢 plink2 两个数量级**——plink2 把 2 位基因型矩阵放进内存后，剪枝只要 0.66 s。这不是常数因子，是执行模型的差别（§6、03 §7）。

## 6. LD 的执行模型问题（本轮最大的性能发现）

`--indep-pairwise` 的 `--timing`：`phase.pass0_ld_ms 64,436`，其中 **`reduce_ms 63,967`**；worker_busy 229 s 均摊到 54 个 worker 只有 4.2 s。即 pass0 的全部剪枝在 reducer（sink）线程里**串行**执行：`crates/reed/src/filter/ld_prune.rs:318-369` 的 `LdPruneSink` 只有一个 `Option<LdPruner>`，`accept` 逐位点 `pruner.push`；模块注释（`ld_prune.rs:7-9`）写的"每 contig 自己的线程、并行"没有实现。

`--ld-decay 100k`：169 s wall、343 CPU-s；chunk 平均 ≈ 500 条记录，基因组跨度小于 100 kb 窗口，所以**绝大多数位点对是跨 chunk 对**，而跨 chunk 对由 sink 在 `accept` 里串行补算（docs/rounds/v5 05 §3.4.2 的"尾窗携带"设计在窗口 ≫ chunk 跨度时退化为串行）。

根因相同：**顺序依赖/全局的分析被塞进了"有序归约"这个串行点**。plink2 的做法是先把基因型矩阵载入内存（2 位/基因型，3.3M × 348 = 290 MB），再按 contig/区块并行计算。reed 的对应设计（03 §7）：pass0 并行解码到**按 contig 的位平面矩阵**（4 平面 × 6 字 × 8 B = 192 B/位点，全量 634 MB；可按 contig 分段），然后 LD 衰减按位置区块并行（只读，无依赖），剪枝按 contig 并行（plink2 的 subcontig 并行），pass1 按位图发射。目标：剪枝 66 s → ≤ 3 s，衰减 169 s → ≤ 10 s（CPU 343 s ÷ 54 ≈ 6.4 s 是下限）。

## 7. reports 的双峰

两组各 3 次：3.13 / 5.02 / 4.55 与 3.04 / 5.32 / 5.41，CPU 恒为 147 s。不是负载（stats 五统计同一时段稳定在 1.75 s），而是某个串行环节时快时慢。本轮的诊断命令被权限分类器拦下，未能定位；留给实现者：用 `--verbosity 2` 记录每次的 `reader_ms / reduce_ms`，并把输出前缀换到 tmpfs 对照一次——嫌疑是 reducer 侧写出约 600 MB 逐位点报告时的页缓存回写，或 summary 收尾。门：3 次极差 ≤ 3%。

## 8. 下一步优化（按收益排序，细则与门在 03/06）

| # | 项 | 现状 → 目标 | 做法 |
|---|---|---|---|
| P6-1 | **LD 矩阵化执行模型** | 剪枝 66.6 s → ≤ 3 s；衰减 169 s → ≤ 10 s | pass0 产出按 contig 的位平面矩阵；衰减按区块并行、剪枝按 contig 并行；pass1 按位图发射（03 §7） |
| P6-2 | deflate 浪费 + 级别策略 | −30–40 CPU-s（memset/alloc）；级别 4 再 −90 CPU-s 换 +7% 体积 | worker 常驻压缩缓冲（`deflate.rs:134`、`block_writer.rs:183`）；默认级别保持 6（与 htslib 行为对齐），文档给出 4 级的数字供用户选 |
| P6-3 | inflate 浪费 | −5–8 CPU-s | 常驻 64 KiB 缓冲，不每块 `resize(…,0)`；不重复解析块头（`inflate.rs:36,49`） |
| P6-4 | `-Ob` 不再文本往返 | −61 CPU-s | 记录模型：解码一次 → 改模型 → 直接编码（03 §3） |
| P6-5 | 扫描 SIMD | 140 CPU-s 中的 −30–50 | 对样本区一次 AVX2 扫 `\t` 与 `:`，替代每样本逐字节 `split_needed`；`memrchr` 替代 1 MiB chunk 的标量 `rposition`（`chunk/mod.rs:231-246`） |
| P6-6 | 每记录分配清零 | −5–10 CPU-s | `stale` Vec、`edit_info` 两个 Vec、`format_g6` String、`Locus` String、每 worker 128 KiB GT 表改 `Arc` |
| P6-7 | `--fill-tags HWE` 走缓存 | −2–5 CPU-s | `SiteCounts::hwe_p` 改用 `HweCache`（同 `(n_called, mac)` 大量重复） |
| P6-8 | stats 串行段 | 1.75 → ≤ 1.45 s | 先补 `reader_wait / reducer_wait` 计量（引擎审计 §2.5），再决定是 mmap 缺页还是收尾合并 |
| P6-9 | reports 双峰 | 稳定 ≤ 3.2 s | §7 的诊断 |
| P6-10 | 索引写出 | 只影响 `-W` | reducer 侧每记录的 contig `String` 与 noodles 对象改整数路径 |

不做：改并行形态（工作窃取、两段流水、NUMA）——引擎审计与 §3 的数字都说明并行结构已在极限上；`target-cpu=native`（规范禁止，可移植性）。

## 9. 测量纪律（沿用 v5 04 §4，本轮新增两条）

- 线程扫描的 user 必须在各点上恒定（±1%），否则说明有随线程增长的同步开销，数字作废。
- 出现双峰（同 CPU、wall 差 > 30%）时不取中位数，按 §7 单独诊断并登记。
