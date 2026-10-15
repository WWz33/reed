# reed 目标架构与工程规范 Spec（v1，2026-10-05）

> **定位**：本文是 reed 的**目标态**（end-state）总规范：架构分层、crate/文件布局、数据模型、
> I/O 与并行执行算法、统计与过滤算法、命名规则（crate/文件/模块/类型/函数/变量/常量/配置/CLI/
> 输出列/产物）、代码规范、测试与性能工程、迁移路线。它允许并建议重构；不对现有代码做任何修改。
>
> **与现有文档的关系**：
> - `docs/rounds/v1/design.md`：统计语义（清单/估计器封死定义/计数优先/拓扑正交）**全部沿用**，本文只重做其实现架构；
>   与本文冲突的工程裁决（如"不写解析代码"）在 §20 决策登记表中逐条列出，等你拍板。
> - `docs/spec/filter.md`：语义（求值序、保真策略、gold 五级验收）沿用；§3 分层、§7 执行引擎
>   由本文 §6–§8、§11 取代。
> - `docs/spec/naming.md`：大部分规则沿用；本文 §13 是其**超集与终版**，冲突处以本文为准（逐条标注）。
> - `docs/spec/statistics-map.md`：作为统计正典来源继续有效。
>
> **标记**：【规】= 规范，必须遵守；【荐】= 推荐默认，可在 ADR 中改；【决】= 需要你拍板（§20）；
> 【据】= 有实测/源码证据。代码引用格式 `crate/path.rs:line`（以 2026-10-05 工作区为准）。

---

## 目录

0. 一页结论
1. 现状审计（保留什么、必须改什么）
2. 设计原则（十二条铁律）
3. 目标架构总览
4. Workspace 与 crate 布局（完整文件树）
5. 数据模型（reed-core）
6. I/O 层（reed-io）：块对齐并行解析、字段计划、快速样本列扫描、并行写出
7. 执行引擎（reed-engine）：融合式 chunk 流水线、确定性、背压、线程预算
8. 统计层（reed-stats）：计数表内核 + 可加累加器 + 窗口分箱
9. 统计算法细则（逐族）
10. QC 统计层（reed-qc）：对标 bcftools stats / vcftools 报告，带 INFO
11. 过滤层（reed-filter）：通用类型化谓词、计划编译器、三种执行策略、字节级编辑
12. 输出格式与产物
13. 命名规范（终版）
14. 代码规范
15. 错误、日志、可观测性
16. 测试与验证体系
17. 性能工程
18. CLI 设计
19. 迁移路线（M0–M8，每步带验收门）
20. 决策登记表（需要你拍板）
附录 A 旧→新文件映射；附录 B 术语表；附录 C 本机数据画像

---

## 0. 一页结论

**你的方向是对的**：bcftools/vcftools 单核，plink2 多核但丢 INFO/FORMAT、改 contig 名，三者都没有
"多核 + 全字段保真 + INFO/FORMAT 过滤 + 正确的群体遗传统计语义"。reed 的统计语义设计（清单、封死
估计器、计数优先、gold 对拍）已经是同类工具里最严谨的，**应当全部保留**。

**但当前实现离"多核快"还差三个结构性台阶**：

1. **统计路径完全没有并行**。`reed-io::for_each_window`（`reed-io/src/lib.rs:168`）是单线程回调流；
   多核只存在于 reed-filter。你最想加速的 "stats vcf" 目前一个核都没多用。
2. **并行切错了位置**。reed-filter 让单个 reader 线程完成"inflate + 全量物化 RecordBuf"，再把对象
   扔给 worker，所以 worker 大半时间在等（文档 §7.2 已诊断）。正确的切法是：**reader 只切压缩字节，
   inflate、解析、计算、序列化、压缩全部在 worker 上按 chunk 融合完成**（§6.2、§7.2）。
   这样串行部分只剩切块和按序 append，开销接近 memcpy。
3. **统计内核重复劳动**。π、dxy、Fst（三种）、PBS、SFS、θ、Tajima、f 统计、HWE 都只依赖
   "每位点 × 每群体的 (ref, alt, het, hom_alt, called) 计数"，但每个估计器各自对 `Vec<i8>` 逐样本
   重扫一遍。目标是每个位点只解码一次，压成位平面（4 bit/样本），用 popcount 一次性算出所有群体的
   计数表；之后所有估计器都是 O(位点 × 群体) 的标量运算，并且全部写成**可加累加器**。窗口、
   滑窗、跨 chunk 合并、block jackknife、多分辨率汇总都由此直接得到（§8）。

**算法层面最重要的三处修正**：

- SFS、θ_W、Tajima's D 当前用 complete-case 语义。本机数据零缺失位点不到 1%（profile：
  site_missing p1 = 0.0086），所以这三个统计量在真实数据上会丢掉 99% 以上的位点。需要引入**超几何
  投影（projection）估计器**作为默认，complete-case 保留为可选的封死定义（§9.4）。
- 过滤器的 INFO 规则写死成 QD/FS/MQ/SOR（`reed-filter/src/spec.rs` `InfoRules`）。"按 INFO 过滤"
  是 reed 相对 plink2 的核心卖点，需要做成**通用、类型化、按 header 校验**的 INFO/FORMAT 谓词（§11.2）。
- 带 cohort 规则的过滤目前要完整解析两遍。可以改成 **sketch-then-emit**：第一遍把掩码后的基因型
  存成位平面（本机约 290 MB）；第二遍只对存活者做字节级编辑，非存活者连字段都不解析（§11.4）。

**性能目标（本机：clean.vcf.gz 6.6 GB，348 样本 × 3.3M 位点，56 核）**：过滤全流程 ≤ 60 s
（现状约 9.6 min，bcftools 约 40 min，plink2 约 2.5 min 但丢字段）；窗口统计（3 群体，π/θ/Tajima/
dxy/Fst，100 kb 窗）≤ 30 s；QC 全报告 ≤ 30 s。估算依据见 §17.1。

---

## 1. 现状审计

### 1.1 必须保留的资产（不许在重构中丢失）

| 资产 | 现位置 | 为什么保留 |
|---|---|---|
| 清单 `SiteCompleteness`/`Ascertainment` + `availability()` 门 | reed-core `manifest.rs`、各估计器 | 同类工具没有的正确性护栏 |
| 估计器 = 封死定义（公式 + 分母 + 缺失语义 + 引用 + `id()`） | reed-core `estimators/*` | 用户选不到"错分母配错公式" |
| 计数优先：结果携带整数分子/分母/缺失 | `result.rs`、`PairCounts` | 可审计，也是并行合并的数学基础 |
| 空分母 → `None`，拒绝 → `Err(Unavailable)` | 全部估计器 | 语义清楚 |
| gold/ 三方金标准 + 五级验收（L_site/L_gt/L_tags/L_manifest/L_props） | `gold/` | 重构期间的回归锚点 |
| 跨线程数字节一致性测试 | `reed-filter/src/exec.rs` tests | 确定性是硬约束 |
| ErrorSlot 首错即停、BatchID 重排、漏斗整数部分和 | `exec.rs` | 机制正确，迁到引擎层复用 |
| 不变量校验（contig 长度必需、升序、ContigRevisit、PositionBeyondContigLength） | `reed-io/src/lib.rs`、`window.rs` | 不静默丢数据 |
| 许可证白名单政策 | `docs/rounds/v1/design.md §8` | 法律约束 |

### 1.2 结构性问题清单（按严重度排序，每条带证据）

| # | 严重度 | 问题 | 证据 | 目标态 |
|---|---|---|---|---|
| A1 | 致命（性能） | 统计路径零并行 | `reed-io/src/lib.rs:168` `for_each_window` 串行回调；无任何线程 | §7 融合式 chunk 并行 + §8 可加累加器 |
| A2 | 致命（性能） | 过滤传输为全物化 RecordBuf，reader 串行 inflate + parse | `reed-io/src/record/mod.rs` `for_each_record` → `record_bufs`；文档 §7.2 实测 283 s/遍 | §6.2 块对齐并行解析 |
| A3 | 高（性能） | emit 线程串行把 RecordBuf 序列化成文本；MT BGZF writer 自建线程池，与计算池叠加导致线程超订 | `exec.rs:653` `emitter.write_record`；`sink.rs` `MultithreadedWriter::with_worker_count(threads)` + `compute_threads` | worker 内完成序列化 + BGZF 压缩，写出阶段只做有序 append（§6.6） |
| A4 | 高（算法） | 每个估计器重扫原始基因型；N 个统计 × M 对群体 = N×M 遍 | `pi.rs`/`dxy.rs`/`fst.rs`/`pbs.rs`/`spectrum.rs` 均逐位点调 `window.allele_counts(site, samples)` | §8 计数表一次构建，所有估计器共享 |
| A5 | 高（算法） | SFS/θ/Tajima 用 complete-case，真实数据下几乎丢光位点 | `spectrum.rs:44` `complete_counts`；profile site_missing p1=0.0086 | §9.4 投影估计器 |
| A6 | 高（语义） | 两套 GT 解码语义不一致：窗口路径拒绝多等位/非二倍体；记录路径静默取前两个等位、单倍体当半调用 | `reed-io/src/lib.rs:545,567` vs `record/vcf.rs:244-260` `take(2)` | §6.4 单一 GT 解码器 + 显式 `PloidyPolicy`/`MultiallelicPolicy` |
| A7 | 高（语义） | `biallelic_snps=false` 时，多等位位点的 AC 把所有 ALT 合并，却按 Number=1 写回 | `site_counts.rs` `count_site`（`a1 > 0` 即 ALT）+ `record/vcf.rs` `apply` | 重算按 Number=A 逐 ALT 计数，或在 plan 期拒绝 |
| A8 | 高（产品） | INFO 过滤写死 QD/FS/MQ/SOR | `spec.rs` `InfoRules`、`record/mod.rs` `RecordView.info_qd..` | §11.2 通用类型化谓词 |
| A9 | 中（架构） | 策略泄漏进传输层：`GtMask`、`Recompute`、`STALE_INFO`、`apply` 在 reed-io | `reed-io/src/record/mod.rs`、`record/vcf.rs` | 传输层只提供"字段访问 + 字节编辑原语"，策略回到 reed-filter |
| A10 | 中（架构） | 依赖方向倒置：reed-io（传输）依赖 reed-core（含估计器）；reed-filter 从 reed-core 估计器模块拿 `SiteCounts` | `reed-io/Cargo.toml`；`exec.rs` `use reed_core::estimators::site_counts` | §3.2 依赖 DAG：模型/内核在最底层，估计器不被传输依赖 |
| A11 | 中（正确性） | `ContigId` 按"首次出现序"分配，`header_info().contigs` 按 header 序给出，调用方无法可靠地从 id 映射回名字 | `reed-io/src/lib.rs:345` `next_contig_id += 1` | `ContigId` = header 序下标，并导出 `ContigTable` |
| A12 | 中（正确性） | contig 长度 `as u32` 静默截断 | `reed-io/src/lib.rs:249,294` | 位置统一 `u64`，或用 checked 转换显式报错 |
| A13 | 中（正确性） | 掩码短于样本数时默认 keep（静默） | `record/mod.rs:146`、`record/vcf.rs:122,205`、`site_counts.rs:89` `unwrap_or(true)` | `SampleMask` 构造时长度必须等于 `n_samples`，不提供隐式默认 |
| A14 | 中（性能/内存） | 重排缓冲 `pending: BTreeMap` 无界；pass0 部分和通道 `unbounded` | `exec.rs:631`、`exec.rs:246` | §7.4 票据窗口（ticket window）背压 |
| A15 | 中（性能） | 滑窗把每个位点的基因型复制 size/step 次；空窗逐个分配 | `lib.rs:419-468` `build_next_window` + `pending.retain` | §8.5 分箱累加：每位点只贡献一次 |
| A16 | 中（性能） | LD 逐样本 f64 循环，O(k²·S)，且只在窗口内配对 | `ld.rs:63-69` `dosage`/`pair_r` | §9.7 位平面 popcount 整数 r²，距离窗口跨窗口边界 |
| A17 | 低（性能） | HWE 精确检验每位点分配 `Vec`，无记忆化 | `hwe.rs:31` | §9.6 worker 内 scratch 复用 + (n, mac, het) 记忆化 |
| A18 | 中（语义） | 文本输出把 VCFv4.4+ header 静默降级为 4.3 | `sink.rs:179` `text_header` | 字节直通后不再重渲染 GT，降级理由消失；需要时显式记录 |
| A19 | 中（可复现） | provenance 哈希基于 `Debug` 输出，跨版本不稳定，且没有写进 header | `plan.rs:48`、`exec.rs:473` `let _plan = ...` | 规范 JSON + 版本号 + sha256，写入 header |
| A20 | 低（工程） | 库代码里有 `eprintln!` 临时探针；CLI 参数解析手写；workspace 没有统一 lints/profile/依赖版本；没有全局分配器 | `exec.rs:214,297,318,480,500,619,633,671,683`；`main.rs` | §14、§15：`tracing` + `clap` + `[workspace.lints]` + release profile + mimalloc |
| A21 | 低（API） | `PopulationView` 要求每个样本都必须归组，没有群体名，组不能重叠 | `population.rs` `validate_against` 用组大小之和比较样本数 | §5.5 `GroupTable`：命名、掩码表示、允许未分组和重叠 |
| A22 | 低（流程） | 性能记录不一致：`gold/out/clean.filtered.reed.log` 记 2422.8 s，文档写约 9.6 min | 日志末行 | §17.3 基准日志协议：每个数字带 commit、命令、缓存状态 |
| A23 | 低（API） | 偏差说明是 `&'static str`，机器不可读 | `manifest.rs` `Availability::Biased(&str)` | `BiasCode` 枚举 + 文案 |
| A24 | 低（API） | 漏斗规则名是字符串 `Vec<(String, u64)>` | `exec.rs` `RunReport.funnel` | `Vec<FunnelStage { rule: RuleId, n_passed }>` |

### 1.3 对现有性能诊断的补充判断

`docs/spec/filter.md §7.3` 的"惰性批传输"方向正确，但它只解决了"物化太重"的问题，
**reader 仍然是单线程（inflate + 切行）**。P3 计划用 noodles `MultithreadedReader` 做并行 inflate，
但它只是把 inflate 并行化，解压后的字节仍然通过一个 `BufRead` 串行流出，切行、建 `Record` 还在一个线程上。
本机解压后约 35 GB 文本（附录 C），单线程 memchr 切行加记录构建要 20–40 s，会成为新的地板。
§6.2 的块对齐并行解析把切行也分给 worker，串行段只剩"按 BGZF 块头切压缩字节"（每块读 18 字节头），
这是能做到的最低串行开销。

---

## 2. 设计原则（十二条铁律）【规】

1. **库优先，组合权在调用方**（沿用 DESIGN §1）。CLI 是库的薄壳，不允许 CLI 有库做不到的能力。
2. **语义与执行分离**。`*Spec`（声明：做什么）→ `*Plan`（编译：怎么做，含资源与策略）→ 引擎执行。
   用户只写 Spec；策略由 plan 编译器按成本模型选，并写进 provenance。
3. **每个字节只解析一次，每个位点只解码一次**。解码结果是位平面和计数表，之后的所有统计都从这里取。
4. **串行段只允许是 O(压缩字节) 的 memcpy 级工作**：读压缩块、按序 append 已压缩的输出块、
   按序合并小型累加器。inflate、解析、计算、序列化、deflate 一律在 worker 上并行。
5. **chunk 是并行单位，chunk 边界只由输入决定**（固定压缩字节目标，并对齐 BGZF 块边界），与线程数无关。
   归约严格按 chunk 序进行。由此输出对线程数逐字节一致（含浮点）。
6. **所有统计量都写成可加累加器**（`push_site` / `merge` / `finish`），`merge` 满足结合律。
   窗口、滑窗、跨 chunk、block jackknife、多分辨率汇总都只是不同的合并方式。
7. **保真靠构造**：未被编辑的记录字节原样直通；被编辑的记录只在字节层面做最小手术（GT 置 `./.`、
   删样本列、替换 INFO 片段）。不经过"解析 → 对象 → 重新渲染"的往返。
8. **传输层不懂策略，策略层不碰字节格式**。reed-io 提供字段访问和编辑原语；"DP<5 就掩码"
   "删掉陈旧 AC" 这类决定属于 reed-filter。
9. **显式胜于默认**：不允许静默截断、静默取前两个等位、静默 keep、静默降级 header。
   遇到不支持的输入要么显式报错，要么按 Spec 里声明的 Policy 处理，并把处理结果计入报告。
10. **确定性 + 可复现**：同一输入、同一 Spec、同一版本 ⇒ 逐字节相同的输出，与线程数和机器无关。
    provenance（版本、规范化 Spec、sha256）写进输出 header 或 sidecar。
11. **内存有界且可预测**：在途内存 = 票据数 × chunk 字节 × 放大系数。超过 `--mem` 预算时，
    plan 改选更省内存的策略，而不是 OOM。
12. **性能主张必须有实测**（沿用 DESIGN §16.4）：每个优化都要有前后对照的基准记录（§17.3），没有提升就撤回。

---

## 3. 目标架构总览

### 3.1 分层

```
L5  reed-cli        二进制 `reed`：filter / stats / qc / … 子命令（clap），只做参数 → Spec 的翻译
L4  reed            门面库：把 L3 的纯策略接到 L2 引擎上的 driver（run_filter / run_stats / run_qc），
                    以及对外 re-export。下游 Rust 用户只依赖这一个 crate。
L3  reed-stats      群体遗传统计：封死估计器、可加累加器、窗口网格        （纯：无 I/O、无线程）
    reed-qc         QC 汇总统计：对标 bcftools stats / vcftools 报告      （纯）
    reed-filter     过滤：Spec、计划编译器、谓词求值、编辑计划、漏斗     （纯）
L2  reed-engine     并行执行底座：chunk 调度、融合流水线、票据背压、
                    首错即停、确定性归约、线程预算、区域分片、阶段指标
L1  reed-io         传输：格式探测、BGZF 并行 inflate/deflate、块对齐切分、
                    VCF/BCF 惰性记录视图、字段计划、快速样本列扫描、字节编辑原语、
                    写出器、CSI/TBI 索引读写
L0  reed-core       模型 + 内核：ID/坐标/样本掩码/群体表/位平面/计数表/清单/窗口网格；
                    popcount 计数内核、HWE 精确检验、调和数表等纯数学     （无依赖）
```

### 3.2 依赖 DAG【规】

```
reed-core ◄── reed-io ◄── reed-engine ◄──┐
    ▲   ▲                                │
    │   └── reed-stats ◄─────────────────┤
    │   └── reed-qc    ◄─────────────────┼── reed ◄── reed-cli
    └────── reed-filter ◄────────────────┘
```

- L3 三个策略 crate **只依赖 reed-core**，单测不需要任何文件和线程，可以直接用内存数据构造。
- 只有 L4 `reed` 同时看得到引擎和策略。"怎么并行"和"算什么"在类型层面就隔开了。
- 禁止反向依赖和同层互相依赖。stats 和 filter 共享的东西（如位点计数）下沉到 reed-core。
- 第三方依赖白名单见 §14.6。

### 3.3 三条主数据流

**(1) 统计 `reed stats`（单遍）**

```
[reader 1 线程]  读 BGZF 块头 → 按固定压缩字节切 chunk（对齐块边界）→ chunk 描述符（Arc<[u8]> 压缩字节 + 续块）
        │ 票据窗口（有界）
[worker × N]     inflate(libdeflate) → 认领本 chunk 起始的行 → 快速扫描 CHROM/POS/REF/ALT + 仅 GT 列
                 → 位平面 GenotypePlanes → popcount 计数内核 → CountTable（位点 × 群体）
                 → 按窗口分箱推入累加器 → ChunkPartial { 分箱部分累加器, 统计计数, 不变量检查状态 }
        │
[reducer 1 线程] 按 chunk 序合并分箱部分 → 水位线推进 → 完成的窗口 finish() → 写出 TSV/JSON
```

**(2) 过滤 `reed filter`（按 Spec 自动选三种策略之一，§11.4）**

```
SinglePass（无 cohort 规则）:
  worker: inflate → 头部谓词（位点字段）→ 需要时解 GT/DP/GQ… → 掩码 → 计数 → 尾部谓词
          → 存活者：LineEditor 字节编辑 → 序列化到 worker 本地缓冲 → BGZF deflate → 压缩块
  writer: 按 chunk 序 append 压缩块 + 同步构建 .tbi/.csi

SketchThenEmit（有 cohort 规则，sketch 能放进内存预算）:
  pass0 worker: 头部谓词位图 + 原始缺失计数 + 掩码后位平面（每 chunk 一份）
  中段:         keep = f(全局缺失)；对位平面按 keep 做 popcount → 尾部谓词 → 存活位图（并行）
  pass1 worker: 非存活行只做 memchr 跳过；存活行 LineEditor → deflate
  writer:       同 SinglePass

TwoPassReparse（sketch 放不下时的兜底）:
  pass0 只算样本缺失归约；pass1 = SinglePass
```

**(3) QC `reed qc`（单遍）**：与 (1) 同构，按需解码的字段换成 QUAL/INFO/FORMAT 数值，
累加器换成直方图和逐样本向量（§10）。

**组合**：`reed stats --filter spec.toml` = 过滤的 pass0（或 SinglePass 的判定部分）+ 统计累加。
不落中间文件，Ascertainment v2 由 plan 自动推导（§11.7）。这就是 filter-architecture.md §9 预留的
`WindowFeeder` 的最终形态。

---

## 4. Workspace 与 crate 布局（完整文件树）【规】

```
reed-rs/
├── Cargo.toml                    # [workspace] members / [workspace.package] / [workspace.dependencies]
│                                 # [workspace.lints] / [profile.*]（§14.3）
├── rust-toolchain.toml           # channel = "1.99"（与本机一致；MSRV 由 rust-version 声明）
├── rustfmt.toml  clippy.toml  deny.toml
├── .cargo/config.toml            # 机器相关：target-dir=/tmp/reed-target、镜像源（不入库，已有）
├── README.md  README.zh-CN.md  CHANGELOG.md
├── docs/rounds/v1/design.md                     # 统计语义决策记录（冻结为历史，新增决策写 ADR）
├── docs/
│   ├── architecture-spec.md      # 本文
│   ├── adr/NNNN-<kebab-title>.md # 架构决策记录（§16.7）
│   ├── statistics-reference-map.md
│   ├── filter-architecture.md    # 语义部分继续有效
│   └── naming-conventions.md     # 被本文 §13 取代后，改为指向本文的短页
├── crates/
│   ├── reed-core/                # L0
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── id.rs             # ContigId SampleId GroupId ChunkId RecordOrdinal（newtype）
│   │   │   ├── coord.rs          # Pos0 Pos1 Interval0 Interval1 + 唯一的 pos0()/pos1() 转换
│   │   │   ├── contig.rs         # ContigTable（名字 ↔ id ↔ 长度）
│   │   │   ├── sample.rs         # SampleTable（名字 ↔ SampleId）、SampleMask
│   │   │   ├── group.rs          # GroupTable（命名群体、掩码、允许重叠/未分组）
│   │   │   ├── bitset.rs         # 定长位集（替代 bits.rs；word 级操作）
│   │   │   ├── genotype/
│   │   │   │   ├── mod.rs        # GenotypePlanes（4 位平面/样本/位点）、Ploidy、Allele 编码
│   │   │   │   ├── planes.rs     # 位平面块 PlaneBlock（多位点连续存放）
│   │   │   │   └── policy.rs     # PloidyPolicy、MultiallelicPolicy、HalfCallPolicy
│   │   │   ├── count/
│   │   │   │   ├── mod.rs        # SiteGroupCounts、CountTable（SoA）
│   │   │   │   └── kernel.rs     # popcount 计数内核（标量 + 运行时分派 SIMD）
│   │   │   ├── site.rs           # SiteKey（contig_id, pos0, ref/alt 类别）、VariantClass
│   │   │   ├── manifest.rs       # SiteCompleteness / Ascertainment v2 / Availability / BiasCode
│   │   │   ├── grid.rs           # WindowGrid（Bp / Sites / Regions）、BinId、窗口 ↔ 分箱映射
│   │   │   ├── math/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── harmonic.rs   # a1(n)/a2(n) 表、Tajima 常数块
│   │   │   │   ├── hypergeom.rs  # 投影用超几何权重（log 域、缓存）
│   │   │   │   ├── hwe.rs        # Wigginton 精确检验 + mid-p + 记忆化
│   │   │   │   └── sum.rs        # Neumaier 补偿求和 NeumaierSum
│   │   │   └── error.rs          # CoreError
│   │   ├── benches/kernel.rs
│   │   └── tests/{genotype_planes,count_kernel,grid}.rs
│   ├── reed-io/                  # L1
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── format.rs         # Format::{Vcf, VcfBgzf, VcfGzip, Bcf} + detect(魔数优先，扩展名兜底)
│   │   │   ├── header.rs         # 读 header；HeaderFacts（样本、contig、INFO/FORMAT 定义与类型）
│   │   │   ├── bgzf/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── block.rs      # 块头解析、BlockRef、虚拟偏移
│   │   │   │   ├── inflate.rs    # libdeflate 解压（worker 本地 Decompressor 复用）
│   │   │   │   └── deflate.rs    # 压缩成 BGZF 块（worker 本地 Compressor 复用）+ EOF 块
│   │   ├── chunk/
│   │   │   │   ├── mod.rs        # ChunkPlan、ChunkSource trait
│   │   │   │   ├── bgzf.rs       # 按块边界切压缩字节 + 续块
│   │   │   │   ├── plain.rs      # 未压缩 VCF：mmap 字节区间
│   │   │   │   ├── stream.rs     # stdin / 普通 gzip：串行 inflate 后切
│   │   │   │   └── index.rs      # 有 CSI/TBI 时按区域切（区域并行、BCF 无索引兜底）
│   │   │   ├── vcf/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── line.rs       # LineRef：一行的字节 + 8 个固定列 + FORMAT + 样本区的偏移
│   │   │   │   ├── fixed.rs      # CHROM/POS/ID/REF/ALT/QUAL/FILTER 快速解析
│   │   │   │   ├── info.rs       # INFO 片段迭代、按 TagId 取类型化值
│   │   │   │   ├── sample_scan.rs# 单遍多字段样本列扫描（GT→位平面，DP/GQ→u16 列）
│   │   │   │   └── edit.rs       # LineEditor：字节级手术（删样本列、GT 置 ./.、INFO 替换）
│   │   │   ├── bcf/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── record.rs     # BcfRecordRef（shared/indiv 两段偏移）
│   │   │   │   ├── decode.rs     # GT int8/int16 向量 → 位平面；类型化 FORMAT 数组
│   │   │   │   └── edit.rs       # BCF 记录编辑（必要时经 noodles RecordBuf 兜底）
│   │   │   ├── field_plan.rs     # FieldPlan：本 pass 要解码哪些列（由各 plan 求并集）
│   │   │   ├── write/
│   │   │   │   ├── mod.rs        # OutputFormat、OutputSpec
│   │   │   │   ├── bgzf_writer.rs# 有序 append 压缩块 + 虚拟偏移记账
│   │   │   │   ├── index_writer.rs # 边写边建 .tbi/.csi
│   │   │   │   └── header_writer.rs# header 渲染 + provenance 行注入
│   │   │   └── error.rs          # IoError（含 Locus 定位）
│   │   ├── benches/{inflate,sample_scan}.rs
│   │   ├── fuzz/                 # cargo-fuzz：sample_scan / line / info 解析
│   │   └── tests/{format_detect,chunk_seams,sample_scan_vs_noodles,line_editor,bcf_decode,roundtrip}.rs
│   ├── reed-engine/              # L2
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── budget.rs         # ThreadBudget、MemBudget
│   │   │   ├── ticket.rs         # TicketWindow：在途 chunk 上限（背压）
│   │   │   ├── pipeline.rs       # run_chunks：reader → worker 池 → 有序 sink 的通用骨架
│   │   │   ├── reorder.rs        # ReorderBuffer（有界，按 ChunkId）
│   │   │   ├── abort.rs          # AbortSignal + FirstError（原 ErrorSlot）
│   │   │   ├── shard.rs          # 区域分片（索引可用时）
│   │   │   ├── metrics.rs        # StageMetrics（每阶段 busy/idle/bytes），tracing 输出
│   │   │   └── error.rs          # EngineError
│   │   └── tests/{determinism,backpressure,abort}.rs
│   ├── reed-stats/               # L3
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── statistic.rs      # Statistic/Accumulator trait、StatSpec 封闭枚举、StatOutput
│   │   │   ├── request.rs        # StatsRequest（值类型 builder）+ 构建期校验
│   │   │   ├── window.rs         # 分箱累加、窗口完成、水位线
│   │   │   ├── jackknife.rs      # BlockJackknife<A>
│   │   │   ├── pi.rs  dxy.rs  fst.rs  pbs.rs
│   │   │   ├── sfs.rs            # 1D/2D SFS（complete-case / projected）
│   │   │   ├── theta.rs          # Watterson θ
│   │   │   ├── tajima.rs         # Tajima's D
│   │   │   ├── fstats.rs         # f2/f3/f4/D/f4-ratio（含 Patterson 校正）
│   │   │   ├── hwe.rs  heterozygosity.rs  inbreeding.rs
│   │   │   ├── ld.rs             # 位平面 r²、LD 衰减分箱
│   │   │   ├── allele_freq.rs    # 原 allele_counts.rs
│   │   │   └── error.rs          # Unavailable、StatsError
│   │   └── tests/{contracts,parity_allel,unequal_groups,accumulator_laws,window_bins}.rs
│   ├── reed-qc/                  # L3
│   │   ├── src/{lib.rs, report.rs, site.rs, sample.rs, info_dist.rs, histogram.rs, sketch.rs, tstv.rs, error.rs}
│   │   └── tests/{parity_bcftools_stats,merge_laws}.rs
│   ├── reed-filter/              # L3
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── spec/
│   │   │   │   ├── mod.rs        # FilterSpec（serde）+ builder + 构建期校验
│   │   │   │   ├── predicate.rs  # SitePredicate / FieldPredicate / Cmp / OnMissing
│   │   │   │   ├── cohort.rs     # CohortRule（MaxMissing / KeepList / RemoveList）
│   │   │   │   ├── region.rs     # RegionRule、Region
│   │   │   │   ├── mask.rs       # GenotypeMaskRule（FORMAT 谓词 → 置 ./.）
│   │   │   │   ├── recompute.rs  # RecomputeTags、StaleTagPolicy
│   │   │   │   └── compat.rs     # from_gold_conf（KEY=VALUE 兼容读入，仅为 gold 对拍）
│   │   │   ├── plan/
│   │   │   │   ├── mod.rs        # FilterPlan + compile(spec, header_facts, budget)
│   │   │   │   ├── strategy.rs   # ExecStrategy::{SinglePass, SketchThenEmit, TwoPassReparse} + 成本模型
│   │   │   │   ├── field_set.rs  # 由谓词推导 FieldPlan
│   │   │   │   ├── ascertainment.rs # Ascertainment v2 降级推导
│   │   │   │   └── provenance.rs # 规范 JSON + sha256 + header 行
│   │   │   ├── eval/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── head.rs       # 位点级谓词求值（只看固定列 + INFO）
│   │   │   │   ├── mask.rs       # 基因型掩码求值（FORMAT 列 → 掩码位图）
│   │   │   │   ├── tail.rs       # 掩码后计数谓词（缺失/MAF/HWE/AC…）
│   │   │   │   └── edit_plan.rs  # 存活者 → EditPlan（交给 reed-io 编辑原语执行）
│   │   │   ├── funnel.rs         # RuleId、FunnelCounter（可加）、FunnelReport
│   │   │   ├── sketch.rs         # pass0 sketch 数据结构（每 chunk 位平面 + 头部位图）
│   │   │   └── error.rs          # FilterError
│   │   └── tests/{spec_validation,predicate_semantics,plan_strategy,funnel}.rs
│   ├── reed/                     # L4 门面
│   │   ├── src/{lib.rs, filter.rs, stats.rs, qc.rs, prelude.rs}
│   │   └── tests/{e2e_filter,e2e_stats,e2e_qc,determinism}.rs
│   └── reed-cli/                 # L5
│       ├── src/{main.rs, cli.rs, cmd/{filter.rs,stats.rs,qc.rs}, output.rs}
│       └── tests/cli.rs          # assert_cmd 黑盒测试
├── bench/                        # 端到端基准与历史后端对照（独立 workspace，保留）
│   ├── e2e/                      # hyperfine 脚本 + 结果日志（§17.3 格式）
│   └── hts/ noo/                 # 既有后端对照（冻结）
├── gold/                         # 外部三工具金标准（规则沿用，不引入 reed 代码）
├── xtask/                        # 【荐】cargo xtask：gold 对拍、基准、发布（替代散落脚本）
└── ci/publish.sh            # 现有发布流程（保留，最终迁入 xtask）
```

**布局规则**【规】：

- 一个 crate 就是一个依赖边界。只有在需要隔离依赖方向或独立复用时才新建 crate。
- 目录模块统一使用 `mod.rs` 风格（与现有代码一致），不要同时存在 `foo.rs` 和 `foo/`。
- 每个文件只承载一个概念。文件名需要用"和"才能描述时，就拆成两个文件。
- 单个文件超过约 800 行（不含测试）就应拆分。`exec.rs` 995 行是反例。
- 单元测试写在文件内的 `#[cfg(test)] mod tests`；跨模块或需要 fixture 的测试放 `tests/`。

---

## 5. 数据模型（reed-core）

### 5.1 标识符（newtype，全部 `Copy + Ord + Hash`）【规】

| 类型 | 底层 | 语义 | 取值来源 |
|---|---|---|---|
| `ContigId(u32)` | header `##contig` 序下标 | 稳定，可经 `ContigTable` 反查名字/长度 | header；不在 header 的 contig → 报错（或 `--allow-undeclared-contig` 时追加） |
| `SampleId(u32)` | header 样本列序下标 | 稳定 | header |
| `GroupId(u16)` | 群体表序下标 | 稳定 | 群体文件出现序 |
| `TagId(u32)` | header INFO/FORMAT 定义序下标 | 字段计划与编辑使用 | header |
| `ChunkId(u64)` | chunk 输入序 | 归约顺序的唯一依据 | reader |
| `RecordOrdinal { chunk: ChunkId, local: u32 }` | 记录在 chunk 内的序号 | sketch 与存活位图对齐 | worker |

- 公共 API 里禁止用裸 `u32`/`usize` 表示上述概念。内部热循环可以 `.get()` 拆箱。
- 现有 `ContigId` 按"首次出现序"分配（A11）。目标态改为 header 序，并导出 `ContigTable`。
  所有结果都携带 `ContigId`，输出层再通过 `ContigTable` 渲染名字。

### 5.2 坐标【规】

```rust
pub struct Pos0(pub u64);              // 0-based
pub struct Pos1(pub u64);              // 1-based（= VCF POS）
pub struct Interval0 { pub start: Pos0, pub end: Pos0 }   // 半开 [start, end)
pub struct Interval1 { pub start: Pos1, pub end: Pos1 }   // 闭区间 [start, end]
impl From<Pos1> for Pos0 { … }   // 全库只有这一处做 -1
impl From<Pos0> for Pos1 { … }
```

- 位置一律用 **u64**（htslib 已支持超过 2^31 的大坐标，见 `htslib/README.large_positions.md`；
  肺鱼、葱属等物种的染色体超过 4 Gb）。热数据（`SiteKey`）里位置只占 8 字节，相对基因型数据可以忽略。
- 内部计算统一用 0-based 半开区间；VCF 文本、用户输入的 region、错误信息、TSV 输出用 1-based（§12.2）。
- 坐标换算只能通过 `From` 实现进行。调用点写 `pos - 1` 这类裸算术，code review 时直接拒绝。

### 5.3 位集与样本掩码

```rust
pub struct Bitset { words: Box<[u64]>, len: usize }           // 定长；尾部多余 bit 恒为 0（不变量）
pub struct SampleMask(Bitset);                                  // len == n_samples（构造时校验）
impl SampleMask {
    pub fn all(n_samples: usize) -> Self;
    pub fn from_ids(n_samples: usize, ids: impl IntoIterator<Item = SampleId>) -> Result<Self, CoreError>;
    pub fn contains(&self, id: SampleId) -> bool;
    pub fn n_selected(&self) -> u32;
    pub fn and(&self, other: &SampleMask) -> SampleMask;        // 长度不同 → 编译期类型相同但运行期 Err
    pub fn words(&self) -> &[u64];                              // 内核用
}
```

- 不再提供"短掩码默认 keep"（A13）。长度不一致是构造错误。
- "保留"统一用 `keep` 一词（cohort 规则），"置缺失"统一用 `mask` 一词（基因型规则），见 §13.4。

### 5.4 基因型位平面（核心表示）

**编码（每样本每位点 4 bit，按"单倍型槽位"拆分）**：

| 平面 | 含义 | 置 1 条件 |
|---|---|---|
| `slot0_alt` | 第 1 个等位基因为 ALT | 槽位 0 的等位基因 ≥ 1（多等位见策略） |
| `slot0_miss` | 第 1 个等位基因缺失 | 槽位 0 为 `.` |
| `slot1_alt` | 第 2 个等位基因为 ALT | 槽位 1 的等位基因 ≥ 1 |
| `slot1_miss` | 第 2 个等位基因缺失 | 槽位 1 为 `.`，或单倍体且策略为半调用 |
| `phased`（可选） | 该基因型是否定相 | 分隔符为 `\|`；只有 `FieldPlan.phase` 开启时才写 |

不变量：`slotX_alt & slotX_miss == 0`，尾部 padding 位全为 0。

- 槽位拆分的理由：能无损表示半调用（`0/.`、`./1`），这是逐等位基因删除语义需要的；
  定相数据可以直接作为单倍型矩阵供未来的 `HaplotypeStream`（EHH/iHS）使用，不必另行解码；
  所有计数都能用 AND/XOR + popcount 得到（§5.6）。
- 内存：每位点 4 × ⌈S/64⌉ 个 u64。348 样本时为 24 words = 192 B/位点；
  现有 `Vec<i8>×2 + missing Bits + phased Bits` 约为 700 B/位点。
- 布局 `PlaneBlock`：一个 chunk 内所有位点按"位点主序、位点内平面连续"存放
  （`words[(site * N_PLANES + plane) * W + w]`），一个位点的 4 个平面落在同一段连续缓存里。

**策略（显式，写进 Spec/Request，默认值保守）**【规】：

```rust
pub enum PloidyPolicy      { RequireDiploid, HaploidAsHalfCall, HaploidAsHomozygous }  // 默认 RequireDiploid
pub enum MultiallelicPolicy { Reject, CollapseAlt, SplitPerAlt }                      // 默认 Reject
pub enum HalfCallPolicy    { AlleleLevel, TreatAsMissing }                              // 默认 AlleleLevel
```

- 每个策略触发的次数都计入报告（如 `n_haploid_as_half_call`）。
- 统计层、过滤层、QC 层共用同一个解码器（修复 A6）。

### 5.5 群体表

```rust
pub struct GroupTable {
    groups: Vec<Group>,            // GroupId = 下标
    n_samples: usize,
}
pub struct Group { pub id: GroupId, pub name: String, pub mask: SampleMask, pub n_samples: u32 }
impl GroupTable {
    pub fn from_tsv(path: &Path, samples: &SampleTable, opts: GroupFileOptions) -> Result<Self, CoreError>;
    pub fn from_assignment(samples: &SampleTable, pairs: &[(&str, &str)]) -> Result<Self, CoreError>;
    pub fn restrict(&self, keep: &SampleMask) -> GroupTable;   // cohort 过滤后组成员自动收缩
    pub fn get(&self, id: GroupId) -> Option<&Group>;
    pub fn by_name(&self, name: &str) -> Option<GroupId>;
}
```

- 允许有样本不属于任何组（被忽略），也允许组之间重叠（例如自动生成的 `all` 组）。
  这修正了 A21 的"必须完整划分"约束。
- 群体文件格式沿用 `sample<TAB>population`（与 `clean.3pop.txt`、pixy、vcftools 一致）。
  出现 VCF 里不存在的样本名时，默认报错；`--ignore-missing-samples` 时降级为警告并计数。

### 5.6 计数内核与计数表

对群体掩码 `G`，单个位点上（`pc` = popcount，均限定在 `& G` 内）：

```
n_alt_alleles      = pc(slot0_alt) + pc(slot1_alt)
n_missing_alleles  = pc(slot0_miss) + pc(slot1_miss)
n_called_alleles   = 2·|G| − n_missing_alleles
n_ref_alleles      = n_called_alleles − n_alt_alleles
full               = ¬slot0_miss ∧ ¬slot1_miss
n_full_called      = pc(full)
n_het              = pc((slot0_alt ⊕ slot1_alt) ∧ full)
n_hom_alt          = pc(slot0_alt ∧ slot1_alt)                 // alt 位隐含已调用
n_hom_ref          = n_full_called − n_het − n_hom_alt
n_both_missing     = pc(slot0_miss ∧ slot1_miss)
n_half_called      = |G| − n_full_called − n_both_missing
```

```rust
#[derive(Clone, Copy, Default)]
pub struct SiteGroupCounts {      // 24 B，按位点主序存放
    pub n_ref_alleles: u32, pub n_alt_alleles: u32,
    pub n_hom_ref: u32, pub n_het: u32, pub n_hom_alt: u32,
    pub n_half_called: u32,
}
pub struct CountTable { n_groups: u16, counts: Vec<SiteGroupCounts> }   // counts[site * n_groups + group]
```

- 现有 `SiteCounts`（AC/AN/NS/…，VCF 标签口径）由 `SiteGroupCounts` 派生，作为 reed-filter 的
  视图函数，不再是 reed-core 估计器模块里的类型（修复 A10）。
- 内核实现：先写标量版作为正典，再加运行时分派的 SIMD 版（AVX2 Harley-Seal、AVX-512 VPOPCNTDQ；
  本机 Xeon Gold 5520+ 支持）。两者用 proptest 做逐位对拍。348 样本时每位点每群体约 40 次
  word 操作，计数阶段相对解析可以忽略；样本数上万时 SIMD 才有意义。

### 5.7 位点键与 chunk 解码块

```rust
pub struct SiteKey {
    pub contig: ContigId, pub pos: Pos0,
    pub class: VariantClass,          // Snv / Mnv / Indel / Symbolic / Ref（单态）/ Other
    pub n_alt: u8,                    // ALT 个数（策略处理前）
    pub flags: SiteFlags,             // is_pass、is_segregating（样本内是否分离，沿用现语义）等
}
pub struct SiteBlock {                // 一个 chunk 解码后的列式结果（worker 本地，可复用）
    pub keys: Vec<SiteKey>,
    pub planes: Option<PlaneBlock>,   // 只有需要个体级数据的统计（LD、个体 F、QC 逐样本）才保留
    pub counts: CountTable,           // 统计层的主输入
    pub columns: FormatColumns,       // 按 FieldPlan 解码的 DP/GQ 等（过滤、QC 用）
}
```

### 5.8 清单 v2（Ascertainment）

采纳 filter-architecture.md §6 的四维结构，本文给出终版：

```rust
pub struct Ascertainment {
    pub sites: SiteSet,            // Complete | Filtered { plan_hash } | Unknown
    pub invariants: InvariantSet,  // Present | PartiallyRemoved { plan_hash } | Absent
    pub genotypes: GenotypeState,  // Original | Masked { plan_hash }
    pub cohort: CohortState,       // Full | Subset { keep_hash, n_samples }
}
pub enum Availability { Available, Biased(BiasCode), Unavailable(RefusalCode) }
#[non_exhaustive] pub enum BiasCode { VariantOnlyDenominator, FilteredNumerator, MaskedGenotypes, … }
```

- 现有 `SiteCompleteness` 通过 `From` 映射到 v2，作为兼容层保留一个大版本。
- `BiasCode`/`RefusalCode` 各有 `fn message(&self) -> &'static str`。机器可读的 code 写进输出
  JSON sidecar（修复 A23）。

### 5.9 窗口网格

```rust
pub enum WindowGrid {
    Bp { size: u64, step: u64 },                 // 沿用现语义：从 contig 首记录所在窗到末记录所在窗；空窗照发
    Sites { n_sites: u32, step: u32 },           // 按位点数分窗（每 contig 重新计数）
    Regions(Arc<RegionSet>),                     // BED 给定（可重叠）
}
```

- `Bp`：分箱宽度 `w = gcd(size, step)`，每个位点只落入一个分箱；窗口 k = 分箱
  `[k·step/w, (k·step + size)/w)` 的合并（§8.5）。
- `Sites`：worker 不知道全局位点序号，所以 worker 输出逐位点的最小充分统计量，由 reducer
  分箱。每位点只有十几个字节，reducer 开销 O(位点数)。
- `Regions`：worker 用区间树把位点分到一个或多个区域；重叠区域各自累加。
- 窗口永不跨 contig（沿用）。空窗发射改为由网格惰性推导，不再逐个分配（修复 A15）。
- contig 末窗截断到 contig 长度（沿用）。另提供 `AccessibilityMask`（BED），使 AccessibleBp 估计器的
  分母等于窗口内可及碱基数，而不是窗口跨度（DESIGN §12.1 留的 v0.3 `is_accessible`）。

---

## 6. I/O 层（reed-io）

### 6.1 格式探测【规】

按魔数判断，扩展名只用于命名输出：

| 首字节 | 进一步判断 | 结果 |
|---|---|---|
| `1f 8b` 且 FEXTRA 含 `BC` 子字段 | 解压首块，以 `BCF\x02` 开头 | `Bcf` |
| 同上 | 以 `##fileformat=VCF` 开头 | `VcfBgzf`（可并行） |
| `1f 8b` 但不是 BGZF | — | `VcfGzip`（只能串行 inflate；发出警告并建议 `bgzip`） |
| `##fileformat=VCF` | — | `Vcf`（未压缩，mmap 并行） |
| 其他 | — | `IoError::UnsupportedInput` |

现有三份扩展名匹配器（`for_each_window`、`for_each_record`、`record_header`）以及 `.gz` 与
`.vcf.gz` 判断不一致的问题，全部由 `Format::detect` 替代。

### 6.2 块对齐并行解析（BAPP：Block-Aligned Parallel Parsing）【规】

这是 reed 在 VCF 吞吐上超过 bcftools/plink2 的核心算法。

**步骤**：

1. **header**：主线程串行 inflate 开头若干块，交给 noodles 解析 header，直到 `#CHROM` 行结束。
   记下第一条记录的虚拟偏移 `(block_offset, in_block_offset)`。
2. **切块（reader 线程，串行，O(压缩字节)）**：顺序读压缩流，每个 BGZF 块只解析 18 字节块头拿到
   `BSIZE`，不做解压。把连续的块累积到压缩字节数 ≥ `CHUNK_TARGET_COMPRESSED`（常量，默认 1 MiB），
   在块边界处切出一个 chunk，分配递增的 `ChunkId`。产出 `ChunkBytes { id, compressed: Arc<[u8]>,
   first_in_block_offset }`。chunk 边界只由输入字节决定，与线程数无关（原则 5）。
3. **worker（并行）**：
   - 用 worker 本地复用的 libdeflate `Decompressor` 把 chunk 的所有块解压到 worker 本地缓冲；
   - 行归属规则：chunk k 的**头片段**是从开头到第一个 `\n`（含）的字节，**尾片段**是最后一个 `\n` 之后
     的字节；中间的完整行归 chunk k 自己处理；
   - 对自有行执行融合流水（§7.2），产出 `ChunkOutput { head_fragment, tail_fragment, partial, out_blocks }`。
4. **接缝缝合（reducer，按 ChunkId 序）**：`tail(k−1) ++ head(k)` 恰好是一条完整记录。
   chunk k−1 以 `\n` 结尾时 tail 为空，head(k) 本身就是完整一行，规则仍然成立。
   reducer 在线处理这一条"接缝行"，并把结果插在 chunk k−1 与 chunk k 的结果之间。
   超长行跨越多个 chunk 时，片段链式拼接，同样成立。
   接缝行约每 500 条记录出现一条，占总工作量约 0.2%。
5. **首个 chunk** 从 header 结束的偏移开始，没有头片段。

**为什么不用 noodles `MultithreadedReader`**：它只并行了 inflate，解压后的字节仍然通过一个
`BufRead` 串行流出，切行和建记录还在单线程上（§1.3）。BAPP 让 worker 自己 inflate 自己的 chunk，
解压后的 35 GB 文本从不经过任何共享线程。

**其他输入的切分方式**：

| 输入 | 切分方式 | 串行段 |
|---|---|---|
| VCF.gz（BGZF） | 上述 BAPP | 读块头 |
| VCF（未压缩） | mmap，按固定字节区间切（对齐 4 KiB），同样做接缝缝合 | 无 |
| BCF + CSI 索引 | 用索引 bin 的虚拟偏移作为记录起点切分（BCF 是二进制，无法靠 `\n` 定界） | 无 |
| BCF 无索引 | 并行 inflate，串行定界（每条记录只读 8 字节长度前缀后跳过），再按批分发 | 定界（廉价） |
| 普通 gzip / stdin | 串行 inflate 线程 → 按字节切 → 同 BAPP | inflate（警告） |

**块头校验**【规】（doc_v10 N4）：块头解析必须拒绝 `BSIZE + 1 < 18` 的块
（htslib `bgzf.c:1219-1225` 同款检查）；畸形块走 `MalformedBgzf`，不得让
下游对短块做切片或长度下溢。

**多成员 gzip**【规】（doc_v10 N6）：普通 gzip 按多成员流解码——第一个成员结束后
继续读下一个成员，等价 htslib `inflate_gzip_block` 的 `Z_STREAM_END` + `inflateReset`
（`bgzf.c:832-908`）。`cat a.vcf.gz b.vcf.gz` 必须两个成员全读。BGZF 的空成员
（EOF marker、串联接缝）同样只跳过、不结束流。

**常量**【规】：`CHUNK_TARGET_COMPRESSED = 1 MiB`（本机约 5 MiB 文本、约 500 条记录），作为引擎常量
而不是用户参数，因为它参与确定性定义（原则 5）。修改它要升版本并记 ADR。

### 6.3 字段计划（FieldPlan）

```rust
pub struct FieldPlan {
    pub fixed: FixedColumns,          // bitflags: CHROM|POS|ID|REF|ALT|QUAL|FILTER
    pub info: SmallVec<[TagId; 8]>,   // 需要的 INFO 标签
    pub format: SmallVec<[TagId; 4]>, // 需要的 FORMAT 键（GT 单独标记）
    pub gt: GtNeed,                   // None | Planes | PlanesWithPhase
}
impl FieldPlan { pub fn union(&self, other: &FieldPlan) -> FieldPlan; }
```

- 每个 plan（过滤、统计、QC）各自导出 FieldPlan，L4 取并集后交给解码器。没进计划的列一个字节都不解析。
- 过滤器的头部谓词只用固定列和 INFO，因此不会触碰样本区（约占每行 99% 的字节）。
  被头部谓词淘汰的记录完全跳过样本区。

### 6.4 VCF 快速解析

**固定列**：用 memchr 找前 9 个 `\t`，得到 `LineRef` 的列偏移。POS 用手写十进制解析；
QUAL 用 `fast-float`，`.` 解析为 None；FILTER 用字节比较判断是否等于 `PASS`。

**INFO**：对 `;` 分隔的片段，用 header 建立的 `TagDict`（完美哈希或排序后二分）把键映射成 `TagId`；
只解析计划内的标签。类型和 Number 按 header 定义：Number=A 在双等位位点取 1 个值；
Flag 判断是否存在；`.` 为 None。

**样本列（热点）：单遍多字段扫描** `sample_scan.rs`【决 D1】：

1. FORMAT 串 → 计划内各键在冒号列表中的下标。按 FORMAT 字节缓存（小型 LRU，同一文件的 FORMAT 通常只有一两种）。
2. 对每个样本：从当前位置开始，只走到最后一个需要的子字段为止，然后用 memchr 跳到下一个 `\t`。
3. **GT 快路径**：如果字节是 `[a, sep, b]`，其中 a、b ∈ {`0`,`1`,`.`}，sep ∈ {`/`,`|`}，且后跟
   `:`/`\t`/`\n`，就用 3 字节查表直接得到 4 个平面位和相位位。否则进入慢路径通用解析：
   多位数等位下标、单倍体、多倍体、整列为 `.`、BCF 风格的 end-of-vector，统一经过 §5.4 的策略处理。
4. DP/GQ 等整数：手写无符号十进制解析，写入 worker 本地复用的 `Vec<u32>` 列，缺失用 `u32::MAX` 表示。
5. 正确性保证：
   - 差分测试 `tests/sample_scan_vs_noodles.rs`：对 fixture 和随机生成的行，快扫描结果必须与基于
     noodles 的参考解码器逐位一致，或双方都报错；
   - cargo-fuzz：随机字节不 panic，结果与参考一致；
   - gold 全量对拍（L_gt/L_tags）。

**依据**【据】：DESIGN §16.2 实测 noodles 每多解一列 FORMAT 多付 28.3%，因为 `series.iter` 对每个字段
都要重扫一遍样本区。单遍扫描对 K 个字段只扫一次。这是本文唯一需要推翻原【定】"不写解析代码"的地方，
理由和适用范围见 §20 D1。header、BCF 结构、输出渲染仍然交给 noodles。

### 6.5 BCF 解码

- 记录头的 `l_shared`/`l_indiv` 给出两段偏移；位点字段从 shared 段惰性读取。
- BCF 的 FORMAT 是按字段连续存放所有样本（字段主序）。GT 是 int8/int16 向量，编码为 `(allele+1)<<1 | phased`，
  可以直接逐字节（或 SIMD）转成 4 个平面。DP/GQ 是类型化整数数组，直接拷贝。
- 不支持的情况（如 GT 为 int32 编码的超大等位下标）走 noodles `RecordBuf` 兜底，并记数。

**描述符与字典的合法性**【规】（doc_v10 N2/N3/N8）：

- `(NULL, count > 0)` 是非法描述符（NULL 无载荷），在唯一的描述符读取关卡拒绝，
  与 htslib `bcf_record_check` 一致（`vcf.c:2231-2236`）；`(NULL, 0)`（INFO flag）合法。
- 头部 `IDX=` 是外部输入：解析时拒绝 `idx >= INT_MAX - 1`（htslib `vcf.c:877-891`），
  字典以索引为键稀疏存储，任何 IDX 都不得驱动按值大小的分配（htslib 的
  `bcf_hdr_set_idx` 会 resize，但以可报告失败收尾，`vcf.c:800-830`）。
- 记录字段字节数 = `width × count × n_sample`，一律用 checked 乘法计算
  （htslib 的 FORMAT 2 GiB 上限同源，`vcf.c:2236-2241`）；溢出即畸形记录。
- **GT 位宽按编码值判**：`(allele+1)<<1 | phased` 的值域决定 INT8/INT16/INT32
  （等位 ≤ 62 → INT8；≤ 16382 → INT16；否则 INT32），与 htslib `bcf_enc_vint`
  按编码值扫 min/max 等价（`vcf.c:2955-2993`）。按等位域判宽会在
  `as i16` 处静默截断（doc_v10 D10-1）。

### 6.6 写出：worker 内压缩 + 有序 append【规】

- 每个 worker 把本 chunk 的输出（原样直通行或编辑后的行）追加到本地字节缓冲，
  再按每块 ≤ 65 280 字节未压缩数据（与 bgzip 一致）切成 BGZF 块，用本地 libdeflate `Compressor` 压缩，
  同时记录每条输出记录在本地缓冲中的未压缩偏移。
- writer（reducer 线程）按 ChunkId 序 append 压缩块，由块长度累计出每条记录的虚拟偏移，
  并同步喂给 noodles 的 CSI/tabix indexer。这样 `--write-index` 不需要额外再读一遍文件。
- 块不跨 chunk，压缩率损失可以忽略（每个 chunk 末尾多一个不满的块）。
- 未压缩 VCF 输出（stdout 默认）同样由 worker 序列化，writer 只做按序 `write_all`。
- header 由 noodles writer 渲染到缓冲区，注入 provenance 行（§11.8），作为文件的第一个（或前几个）BGZF 块。
- 文件尾写 BGZF EOF 块。
- 压缩级别：`--compress-level`，默认 6（与 bgzip/bcftools 默认一致）；`-O u` 输出未压缩 BCF。
- 不再使用 noodles `MultithreadedWriter`：它自带线程池，会和计算池叠加造成超订（A3）。

### 6.7 字节编辑原语（LineEditor）

```rust
pub struct EditPlan<'a> {
    pub keep_samples: Option<&'a SampleMask>,  // None = 不删列
    pub mask_gt: Option<&'a Bitset>,           // 置 ./. 的样本（按输入样本序）
    pub info_remove: &'a TagSet,               // 陈旧标签
    pub info_set: &'a [(TagId, InfoValueOut)], // 重算标签（已格式化或待格式化）
}
pub fn edit_vcf_line(line: &LineRef, plan: &EditPlan, header: &HeaderFacts, out: &mut Vec<u8>) -> Result<(), IoError>;
```

- 列 1–7 和 FORMAT 原样拷贝。
- INFO：遍历原始片段，跳过 `info_remove` 中的标签；被覆盖的标签保留在原位置并替换值；
  新标签按 `info_set` 给定的规范顺序追加在末尾。具体顺序以与 gold（bcftools fill-tags）逐字节一致为准，
  由测试锁定。
- 样本列：被删的列跳过；被掩码的样本只替换 GT 子字段（二倍体为 `./.`，单倍体为 `.`），
  其余子字段（AD/DP/GQ/PL…）原样保留。
- 浮点格式化：`FloatFormat::HtslibCompat`，与 htslib 渲染 float32 的方式一致，保证 L_tags 可以精确对拍。
- 没有任何编辑的存活者直接整行拷贝（零重序列化）。
- 输出 BCF 时：文本 → BCF 需要完整编码。先走 noodles `RecordBuf`（正确但较慢），后续再写原生 BCF 编辑器（M8）。
  文档要说明：VCF.gz 是最快的输出格式，BCF 是体积更小的格式。

### 6.8 索引与区域

- 读：noodles-csi / noodles-tabix。`-r/--regions` 要求有索引，按索引跳读；`-t/--targets` 是流式谓词，
  不需要索引（沿用 bcftools 语义区分）。
- 写：§6.6 的边写边建。
- 区域分片（§7.10）用索引把基因组切成压缩字节数大致相等的分片。

---

## 7. 执行引擎（reed-engine）

### 7.1 通用骨架

```rust
pub trait ChunkJob: Sync {
    type Partial: Send;
    type Scratch: Send + Default;                      // worker 本地复用缓冲（解压、列、平面、输出）
    fn run_chunk(&self, chunk: &ChunkBytes, scratch: &mut Self::Scratch)
        -> Result<ChunkOutput<Self::Partial>, EngineError>;
    fn run_seam(&self, line: &[u8], scratch: &mut Self::Scratch)
        -> Result<ChunkOutput<Self::Partial>, EngineError>;   // 接缝行，reducer 上执行
}
pub trait OrderedSink<P>: Send {
    fn accept(&mut self, id: ChunkId, output: ChunkOutput<P>) -> Result<(), EngineError>;
    fn finish(self) -> Result<SinkReport, EngineError>;
}
pub fn run<J: ChunkJob, S: OrderedSink<J::Partial>>(
    source: ChunkSource, job: &J, sink: S, budget: &ThreadBudget,
) -> Result<RunReport, EngineError>;
```

reed-filter、reed-stats、reed-qc 都不直接碰线程。L4 的 driver 把它们的纯逻辑包装成 `ChunkJob` 和 `OrderedSink`。

### 7.2 融合式 chunk 处理【规】

worker 对一个 chunk 一口气完成：**inflate → 认领行 → 解析（按 FieldPlan）→ 解码到位平面/列 →
计数 → 谓词/累加 → 编辑/序列化 → deflate**。

- 中间结果全部留在 worker 本地缓存和本地内存里，不在阶段之间跨线程传递对象。
- 跨线程传递的只有：压缩输入（`Arc<[u8]>`，约 1 MiB）、压缩输出块、小型部分累加器。
- 不在线程之间传递 per-record 堆对象，所以 `filter-architecture.md §7.2` 诊断的远程 free 问题从结构上消失。

### 7.3 线程预算【规】

`--threads N` 是唯一的计算资源旋钮，表示 reed 使用的线程总数（默认为可用核数）：

| 角色 | 线程数 | 说明 |
|---|---|---|
| reader | 1 | 读块头、切 chunk；大部分时间阻塞在 I/O 或票据上 |
| reducer/writer | 1 | 按序合并、接缝缝合、append、建索引 |
| worker | max(1, N − 2) | 全部 CPU 重活 |

- 进程里不允许出现第二个线程池（不用 noodles MT reader/writer；rayon 只在 pass 之间的批量运算里使用，
  例如 sketch 的尾部判定，并且复用同一个预算）。**落地形态（doc_v8 D8-1/D8-6）**：rayon 阶段不使用
  全局池，经 `ThreadBudget::auxiliary()`（= `workers()`，顺序预算 = 1）建每调用一次的临时池；
  串行预算下闭包收到 `parallel == false`，必须走普通迭代器（`par_iter` 在串行分支仍会惰性
  拉起全局池，D8-6 的教训）；per-contig 线程上限（filter pass0 的 `max_prune`、stats LD sink
  的 in-flight cap）同取此值。
- N ≤ 2 时退化为 reader 和 worker 合一的顺序模式。结果必须与并行模式逐字节相同（测试锁定）。
  辅助 rayon 阶段与 per-contig 线程同样走串行/受限路径（`auxiliary() == 1`，doc_v8 D8-1）。

### 7.4 背压：票据窗口【规】

- 容量 `T = TICKETS_PER_WORKER(=2) × n_workers`。reader 每切出一个 chunk 先取一张票；reducer 按序消费掉
  chunk k 后归还它的票。
- 在途 chunk 数 ≤ T，因此重排缓冲可以用定长环形数组（下标 `id % T`），不需要 `BTreeMap`（修复 A14）。
- 内存上界 ≈ T × (压缩输入 + 输出块 + 部分累加器) + n_workers × scratch。
  本机为 108 × (1 MiB + 约 3 MiB) + 54 × 约 20 MiB ≈ 1.5 GB，可预测。
- 所有通道都有界。现有 pass0 的 `unbounded` 通道取消。

### 7.5 中止与错误

- `AbortSignal`：一个只用来关闭的 crossbeam 通道。所有阻塞点都用 `select!` 同时等待数据和 abort，
  中止即时生效，不再用 20 ms 轮询。
- `FirstError`：首个错误胜出（沿用 ErrorSlot 语义）。worker panic 时由 `catch_unwind` 转成
  `EngineError::WorkerPanicked { chunk }`；reducer panic 转成 `EngineError::ReducerPanicked { chunk }`。
  **任何线程（reader / worker / reducer）的错误与 panic 都必须先经 `FirstError` 记录
  （记录即触发 abort），不允许绕过它直接返回**——否则其余阻塞点失去唯一的无条件唤醒来源，
  reducer 侧出错会变成整个 run 的挂死（doc_v9 §2 规约 W2）。
- 错误必须带定位：`Locus { contig: String, pos1: u64 }` 加 `ChunkId`。

### 7.6 确定性契约【规】

1. chunk 边界只由输入决定（§6.2）。
2. worker 内部按记录序处理，部分累加器按记录序累加。
3. reducer 严格按 `ChunkId` 序合并，接缝行插在两个 chunk 之间。
4. 浮点累加使用 `NeumaierSum`（补偿求和）。即使将来 chunk 常量改变，结果差异也只在最后一两位。
5. 由 1–4 可得：同一版本、同一输入、同一 Spec ⇒ 输出对线程数逐字节一致（含浮点和压缩块）。
6. 测试：`reed-engine/tests/determinism.rs` 与 `reed/tests/determinism.rs`，覆盖 threads ∈ {1, 2, 3, 7, 56}，
   比较输出文件的 sha256。

### 7.7 内存预算与策略选择接口

`MemBudget`（`--mem`，默认为物理内存的 50%）由 plan 编译器读取，用来在策略之间选择（§11.4）。
引擎负责运行时监控（在途**估算足迹**计数，= 压缩字节 × `MEM_FACTOR`）。预算以**背压**方式执行：
超出时 reader 阻塞等待许可归还，**不报错**；`EngineError::MemBudgetExceeded` 只保留给"单个 chunk
的足迹就超过全部预算"这种无法继续的情况（依据 `docs/rounds/v2/03-architecture.md` §4.2；
本段由 doc_v9 D9-7 回改，不改代码）。不会默默 OOM。

### 7.8 NUMA 与分配器

- 解压缓冲、列缓冲、输出缓冲都由 worker 自己分配（first-touch 落在本地 NUMA 节点）。reader 只分配压缩字节。
- 二进制使用 mimalloc 作为全局分配器【决 D7】。库 crate 不设全局分配器。
- `--pin-threads`（可选）：按 NUMA 节点交错绑核。默认关闭，只有实测有收益才默认打开。

### 7.9 区域分片（有索引时）

- 用 CSI/TBI 把基因组切成 S 个压缩字节大致相等的分片，每个分片是一个独立的 chunk 源，
  使用独立文件句柄，可以并发读（对 /sri 这类并行文件系统有利）。
- ChunkId 为 `(shard, local)` 字典序，确定性保持。
- 用于：`-r` 区域查询；远程/慢速存储；以及无索引 BCF 之外的另一条并行路径。

### 7.10 阶段指标

每个阶段记录忙碌时间、空闲时间、字节数、记录数（`StageMetrics`），通过 `tracing` 输出；
`--profile run.json` 把指标落盘，用于基准日志（§17.3）。现有 `TEMP instrumentation` 的 `eprintln!`
全部由此替代（A20）。

---

## 8. 统计层（reed-stats）

### 8.1 核心抽象【规】

```rust
/// 一个位点交给累加器的全部输入（只读视图）。
pub struct SiteInput<'a> {
    pub key: &'a SiteKey,
    pub counts: &'a [SiteGroupCounts],       // 长度 = n_groups，下标 = GroupId
    pub planes: Option<SitePlanes<'a>>,      // 仅个体级统计需要
    pub derived: &'a SiteDerived,            // 每位点每群体只算一次的派生量（p_alt、组内 mpd…）
}

pub trait Accumulator: Clone + Send {
    fn push_site(&mut self, site: &SiteInput<'_>);
    /// 必须满足结合律；整数分量必须精确，浮点分量用 NeumaierSum。
    fn merge(&mut self, other: &Self);
}

pub trait Statistic: Clone + fmt::Debug + Send + Sync {
    type Acc: Accumulator;
    type Output;
    fn id(&self) -> StatId;                                   // 稳定标识，如 "fst.hudson"
    fn citation(&self) -> &'static str;
    fn availability(&self, asc: &Ascertainment) -> Availability;
    fn needs(&self) -> InputNeeds;                            // Counts | Planes | Positions | Columns
    fn new_acc(&self, ctx: &StatContext) -> Self::Acc;        // ctx：群体大小、网格、可及性掩码
    fn finish(&self, acc: &Self::Acc, frame: &WindowFrame) -> Self::Output;
}
```

- 估计器仍然是**封死定义**：每个 `Statistic` 实现就是 `(统计量, 估计器, 群体组合)` 三元组
  （如 `HudsonFst { a, b }`），公式、分母、缺失语义、引用和 `availability` 绑在一起。
- 异构请求集用**封闭枚举**静态分发（`StatSpec`/`StatAcc`/`StatOutput`），不用 `dyn`（沿用 DESIGN §16.3 的裁决）。
- `SiteDerived` 是性能关键：π、dxy、Fst、PBS 共享"每群体组内平均成对差异"和"每群体等位频率"，
  每个位点每个群体只算一次。群体对之间的量在需要时才计算。
- 保留一次性自由函数（`reed_stats::pi(&block, group, PiEstimator::AccessibleBp)`），供库用户和测试使用。
  它们内部就是 `new_acc → push_site* → finish`。

### 8.2 请求（StatsRequest）

```rust
let request = StatsRequest::builder(&groups, WindowGrid::Bp { size: 100_000, step: 50_000 })
    .ascertainment(asc)                         // 或从输入 header 的 ##reed_ascertainment 读取（§11.7）
    .add(Pi::per_pair_comparable(pop1))
    .add(Fst::hudson(pop1, pop2))
    .add(Fst::matrix(FstEstimator::Hudson))     // 所有群体对，共享逐群体中间量
    .add(Sfs::projected(pop1, 600))
    .add(TajimaD::projected(pop1, 600).min_sites(3))
    .build()?;                                   // 构建期：组存在、互不相同、availability 拒绝即 Err
```

- 请求是值类型（`Clone + Debug + serde`），可以从 TOML 读入（CLI 用）。
- `build()` 时就按清单完成全部可用性判定：被拒绝的组合在构建期报错，而不是跑到一半才失败。

### 8.3 计算流（L4 driver 内）

1. worker：解码 chunk → `CountTable` + `SiteKey`（+ 需要时的 `PlaneBlock`）→ 对每个位点算
   `SiteDerived` → 确定分箱 → 推入该分箱的 `StatAcc[]`。
2. worker 输出 `ChunkPartial { bins: Vec<(BinKey, Vec<StatAcc>)>, n_sites, policy_counters, invariants }`，
   只包含被触及的分箱（稀疏）。
3. reducer：按 ChunkId 序合并分箱部分；用水位线判定哪些窗口已经完整（窗口最后一个分箱 < 当前 chunk 的
   第一个分箱，或 contig 已切换）；完整的窗口 `finish` 后按序写出。
4. 不变量检查（位置升序、ContigRevisit、超出 contig 长度）跨 chunk 进行：worker 检查 chunk 内部，
   reducer 检查 chunk 之间的衔接（前一 chunk 末位点 vs 后一 chunk 首位点）。

### 8.4 累加器定律（测试锁定）【规】

- 单位元：`Acc::default()`（或 `new_acc`）与任何 `a` 合并都等于 `a`。
- 结合律：`(a⊕b)⊕c == a⊕(b⊕c)`。整数分量要求精确相等，浮点分量要求 ≤ 4 ulp。
- 切分不变：把位点序列任意切成若干段分别累加再合并，结果等于整体顺序累加（同样的精度要求）。
- `tests/accumulator_laws.rs` 用 proptest 对每个 `Statistic` 自动验证以上三条。

### 8.5 分箱窗口与滑窗【规】

- `Bp { size, step }`：分箱宽度 `w = gcd(size, step)`。每个位点只推入一个分箱；
  窗口 = 连续 `size / w` 个分箱的合并。每个窗口的合并成本是 O(size/w)。典型 `size/w ≤ 20`。
- 不使用"前缀和相减"：浮点减法会破坏确定性和精度。
- 与现状对比：现在每个位点被复制 size/step 次（A15）；分箱后每个位点只解码、只累加一次。
- 空窗：网格推导出来、没有位点的窗口照常输出（`value = None`，计数为 0），可用 `--skip-empty` 关闭。

### 8.6 Block jackknife（通用）

- 对任意 `Statistic`，`BlockJackknife<S>` 以分箱（或独立的块网格，如 5 Mb / 500 位点）为块。
  用前缀合并数组和后缀合并数组得到 leave-one-out 累加器，总成本 O(B)，并且全部是确定性的合并操作。
- 块大小不等时用加权 jackknife（Busing et al. 1999，ADMIXTOOLS 语义；snputils BSD 可移植）。
- 输出 `est / se / z / p`。f 统计、Fst、D 统计都可以直接使用。

### 8.7 多分辨率

全基因组值 = 所有分箱合并；逐 contig 值 = 该 contig 分箱合并。不需要再扫一遍数据，一次运行即可同时输出
窗口、contig、全基因组三个层级（`--levels window,contig,genome`）。这也修复了 DESIGN §16.1-11 指出的
PBS 无法重聚合问题：PBS 的累加器保存三对的 num/den，而不是 clip 之后的值。

### 8.8 逐位点输出

`--per-site` 用于位点级统计（逐位点 Fst、HWE、Hobs/Hexp、等位频率）。由 worker 直接序列化成文本片段，
按 chunk 序 append，不经过分箱。

---

## 9. 统计算法细则（逐族）

> 每一族都写清：输入、公式要点、缺失语义、变体、对拍正典、相对现状的变化。正典来源沿用
> `docs/spec/statistics-map.md`。

### 9.1 π / dxy

- 输入：`n_ref_alleles`、`n_alt_alleles`（每群体）。
  π：`diffs = ref·alt`，`comps = n(n−1)/2`；dxy：`diffs = refA·altB + altA·refB`，`comps = nA·nB`。
- 估计器沿用 `PerPairComparable` 与 `AccessibleBp`。`AccessibleBp` 的分母改为
  "窗口内可及碱基数"（有 `AccessibilityMask` 时）或窗口跨度（没有时，与现状一致）。
- 整数溢出分析：单位点 comps ≤ (2S)²/2。窗口内求和用 u64；全基因组层级用 u128。debug 构建下 checked_add。
- 与现状逐位一致（整数计数完全相同，浮点 ≤ 1e-12）是 M5 的验收门。

### 9.2 Fst

- **Hudson**：保持 allel 逐位点归一化语义（num 需双侧 n≥2，den 需双侧 n≥1）。
- **WC84**：推广到 r 个群体（去掉 `fst.rs` 里写死的 `(2.0 − 1.0)`），n_i 采用 allel 的 `an // 2` 地板语义，
  h̄ 来自计数表的 `n_het`，不再单独调 `genotype_counts`。
- **Nei86**：推广到多群体。
- **新增**：`HudsonFiniteSample`（snputils 的 `π·n/(n−1)` 校正）作为第四个估计器。
- **矩阵请求**：`Fst::matrix` 一次计算所有群体对。每位点每群体的组内项只算一次，群体对的项是
  O(G²) 次标量运算。
- 审计线（`PairCounts`）改为按需：结果里保留整数分量，但只有请求了 `--audit` 时才输出全部审计列。

### 9.3 PBS

- 复用三对 Hudson Fst 的累加器（消除 `pbs.rs` 对 `fst.rs` 的代码复制）。
- 位点集合变体：`SiteSet::Intersection`（现语义）与 `SiteSet::PerPair`（allel 语义）作为两个封死定义
  （解决 DESIGN §11 的待定项），默认值【决 D5】。
- clip 区间保持 `[0, 0.99999]`（allel）。

### 9.4 SFS / θ_W / Tajima's D：缺失数据（重点修正）

**问题**【据】：本机 profile 的位点缺失分位数为 p0 = 0、p1 = 0.0086，即零缺失位点不到 1%（这还是掩码前）。
现在的 complete-case 语义（`spectrum.rs:44` `complete_counts`）会丢弃 99% 以上的位点。结果仍然"正确"，
但几乎没有可用信息，而且被保留的位点本身有偏（高覆盖区域）。

**三种封死定义**：

| 估计器 | 定义 | 适用 |
|---|---|---|
| `CompleteCase` | 现语义：组内任一等位缺失即跳过 | 无缺失或缺失极少的数据；与 allel 对拍 |
| `Projected { n }`（**新默认**【决 D4】） | 超几何投影（Marth et al. 2004；dadi `project`；easySFS）：调用数 `m ≥ n` 的位点按 `P(k \| n, m, alt) = C(alt,k)·C(m−alt, n−k) / C(m, n)` 把一个位点分摊到 SFS 的各个分箱；`m < n` 的位点跳过并计数 | 真实数据默认。θ_π、θ_W、Tajima's D 都在投影后的 SFS 上计算，n 固定，方差公式成立 |
| `PerSiteN` | θ_W = Σ_分离位点 1/a_{n_i}，π 用逐位点 mpd（allel 语义） | 只需要 θ、不需要 D 时；D 的可变 n 方差（Ferretti et al. 2012）列为远期 |

- 默认 `n`：取 `⌊(1 − q) · 2·|G|⌋`，其中 q 为该群体逐位点缺失等位比例的 90 分位（由 pass 内统计得出，
  会使用 2 遍）。也可以由用户显式给出 `n`。CLI 推荐显式给出，以保证可复现。
- 实现：对每个 `(m, alt)` 用比值递推 `P(k+1)/P(k) = (alt−k)(n−k) / ((k+1)(m−alt−n+k+1))` 生成权重，
  并在相对权重低于 1e−15 时截断，实际成本约 O(√n)/位点。worker 内用 LRU 缓存 `(m, alt) → 权重切片`。
- θ_π = Σ_k k(n−k)·ξ_k / C(n,2)；S = Σ_{k=1}^{n−1} ξ_k；Tajima 常数块使用 `math/harmonic.rs` 的预计算表。
- 对拍：projected SFS 与 easySFS/dadi 的投影结果对拍（新增 gold 项，§16.4）。

### 9.5 f 统计

- f2/f3/f4/D 加 Patterson 无偏校正（allel `admixture.py` MIT / snputils BSD），以 `Corrected`/`Naive`
  两个封死定义提供。默认 `Corrected`。
- 新增 f4-ratio。
- 标准误统一走 §8.6 的 `BlockJackknife`。

### 9.6 HWE 精确检验

- 算法沿用 Wigginton 2005（与 vcftools/sgkit/bcftools 同源），并新增 mid-p 变体（Graffelman & Moreno 2013）。
- **记忆化**：键为 `(n_full_called, n_minor_alleles)`，缓存该 `(n, mac)` 下全部 het 取值的 p 值数组
  （已归一化并按概率排序后累加）。之后任意 het 的查询都是 O(1)。worker 内使用有界 LRU（默认 8192 项）。
  过滤器的 HWE 规则和 QC 的 HWE 分布都受益。
- scratch 复用，消除逐位点 `Vec` 分配（A17）。
- 退化位点哨兵（fill-tags 口径 p = 1）保留在 reed-filter 的视图函数中，不放进数学内核。

### 9.7 LD

- **位平面整数 r²**（Rogers-Huff 未定相基因型相关）：对两个位点 x、y 和有效集合 `V = full_x ∧ full_y ∧ G`，
  记 `a = slot0_alt`、`b = slot1_alt`：
  ```
  n    = pc(V)
  Σx   = pc(a_x∧V) + pc(b_x∧V)
  Σx²  = Σx + 2·pc(a_x∧b_x∧V)                 // a² = a
  Σxy  = pc(a_x∧a_y∧V) + pc(a_x∧b_y∧V) + pc(b_x∧a_y∧V) + pc(b_x∧b_y∧V)
  ```
  协方差和方差用 i64 精确计算，最后一步才转 f64。每对位点约 12 次 popcount × W；
  现状是逐样本的 f64 分支循环（`ld.rs:63-69`）。整数化同时解决确定性问题。
- 定相时提供 Hill-Robertson r²、D、D′（直接用单倍型平面）。
- **拓扑**：LD 的配对距离（常见 100 kb–1 Mb）远大于一个 chunk 的跨度（约 150 kb），不适合在 chunk 内做。
  采用"两段式"：
  1. 并行把整个 contig 解码成平面存储（本机 192 B × 3.3M ≈ 630 MB；超出预算时按 contig 溢写）；
  2. 把位点切成块并行计算，每块向后回看 `max_bp`/`max_sites`。
  这同时是 DESIGN §2 所说 `MatrixAccumulator` 拓扑的第一个实例。
- 输出：逐对（`--ld-pairs`，按距离和 r² 阈值裁剪）、LD 衰减（按距离分箱的 Σr² 与计数，可加累加器）、
  窗口汇总（现有 `LdSummary` 语义）。
- LD pruning（ordered-carry，plink2 语义）复用平面存储：每个 contig 内顺序贪心，contig 之间并行
  （以及在大间隔处切成子 contig 并行，与 plink2 `LdPruneSubcontigSplitAll` 同构）。

### 9.8 杂合度与近交

- 逐位点 Hobs/Hexp/Fis 沿用（p = alt/(ref+alt) 的修复保留）。
- 逐个体 F（plink `--het` 矩估计：(O_hom − E_hom)/(N − E_hom)）是逐样本可加向量。
  从位平面累加逐样本计数时使用 positional popcount（位切片进位保存加法器，每 64 个位点做一次进位折叠），
  避免逐样本逐位点循环。

### 9.9 等位频率与联合 SFS

- 逐位点、逐群体 AF/AC（对标 vcftools `--freq`/`--counts`、plink `--freq`）。
- 2D 联合 SFS（allel `sf.py` MIT），支持投影，直接输出 dadi/fastsimcoal 可读格式。

### 9.10 数值规范【规】

- 整数分量：窗口层级 u64，全基因组层级 u128，加法溢出在 debug 下 panic、release 下由测试覆盖。
- 浮点分量：`NeumaierSum`。`value()` 最后一步才做除法。
- 输出格式：`ryu` 最短往返表示，避免 `{:.6}` 截断。

---

## 10. QC 统计层（reed-qc）：`reed qc`

**动机**：你提到的 "stats vcf" 除了群体遗传统计，还包括 bcftools stats / vcftools 风格的质控汇总。
这类汇总全部是可加的，最容易做到线性加速；而且只有 reed 能在多核下同时给出 INFO/FORMAT 的分布（plink2 做不到）。

### 10.1 报告目录（对标表）

| 节 | 内容 | 对标 |
|---|---|---|
| `summary` | 记录数、样本数、SNV/MNV/indel/symbolic/多等位/无 ALT 计数 | bcftools stats SN |
| `tstv` | 总体 Ts/Tv；按 AF 分箱的 Ts/Tv | bcftools TSTV / AF；vcftools `--TsTv-summary` |
| `substitution` | 12 种替换类型计数 | bcftools ST |
| `indel_length` | indel 长度分布 | bcftools IDD |
| `af_spectrum` | ALT 频率分箱计数（全体/每群体） | bcftools AF |
| `qual` | QUAL 直方图 + 分位数 | bcftools QUAL；vcftools `--site-quality` |
| `info.<TAG>` | **任意数值 INFO 标签**（QD/FS/MQ/SOR/ReadPosRankSum/MQRankSum/DP/ExcessHet/…）的直方图 + 分位数 + 缺失数 | **独有**：plink2 不能，bcftools 需要逐个 query |
| `format.<KEY>` | 基因型级 DP/GQ 等直方图（精确整数直方图） | bcftools DP |
| `sample` | 逐样本：hom_ref/het/hom_alt/missing/singleton/ts/tv/平均 DP/平均 GQ/het 比/个体 F | bcftools PSC；vcftools `--missing-indv`、`--het`、`--depth`、`--singletons` |
| `site`（流式，可选） | 逐位点：缺失率、MAF、HWE p、平均 DP | vcftools `--missing-site`、`--freq`、`--hardy`、`--site-mean-depth` |
| `density` | 每窗口位点数 | vcftools `--SNPdensity` |
| `group`（给出群体时） | 每群体私有等位、固定差异、分离位点数 | — |
| `mask_impact`（给出过滤 Spec 时） | 掩码前后缺失率、掩码基因型数按规则归因 | 独有（profile_input.py 的产品化） |

### 10.2 可合并的分布结构【规】

- 整数字段（DP、GQ、长度）：稠密精确直方图 `0..=MAX` 加溢出桶。分位数精确。
- 浮点字段：固定边界直方图（默认边界随标签内置，可配置）用于绘图；**DDSketch**（相对误差 1%，可合并、确定性）用于分位数。
  报告中写明分位数是近似值及其误差界。
- 逐样本向量：`Vec<u64>` 按样本累加，使用 positional popcount。
- 以上全部满足 §8.4 的累加器定律，因此 QC 与统计层共用同一套 chunk 并行驱动。

### 10.3 输出

- `<prefix>.qc.json`：全部节的结构化结果，带 provenance。
- `<prefix>.<section>.tsv`：每节一个 TSV，供 R/pandas 直接读入。
- `--format bcftools-stats`【荐】：额外输出 bcftools stats 兼容文本，可以直接喂给 `plot-vcfstats` 和 MultiQC，降低用户迁移成本。

---

## 11. 过滤层（reed-filter）

### 11.1 Spec v2（类型化、可序列化、按 header 校验）

不引入字符串表达式语言（沿用 DESIGN §16 的裁决）。Spec 是结构化数据，可以从 TOML/JSON 读入，也可以用 builder 构造：

```toml
# filter.toml
[cohort]
max_missing = 0.20                 # 或 keep_list = "samples.keep.txt" / remove_list = "…"

[targets]                          # 流式谓词（bcftools -t 语义）；按索引跳读用 CLI -r
include = ["Chr01:1-2000000"]
exclude = []

[site]
variant_class = ["snv"]            # snv | mnv | indel | symbolic | ref
max_alt_alleles = 1
require_filter = ["PASS"]          # FILTER ∈ 集合；"." 不算通过
min_qual = 30.0

[[site.info]]                      # 任意数值 INFO；按 header 校验存在性、Type、Number
tag = "QD"
keep_if = ">= 4.0"
on_missing = "pass"                # pass | fail；默认 pass（bcftools 缺失比较为假的语义）

[[site.info]]
tag = "FS"
keep_if = "<= 50"

[[site.info]]
tag = "ReadPosRankSum"             # 现状无法表达
keep_if = ">= -8.0"

[[genotype.mask_if]]               # 任一为真 → GT 置 ./.（其余 FORMAT 保留）
field = "DP"
when = "< 5"
on_missing = "keep"                # keep = 缺失值不触发掩码（bcftools 语义）

[[genotype.mask_if]]
field = "DP"
when = "> 100"

[[genotype.mask_if]]
field = "GQ"
when = "< 20"

[site_after_mask]                  # 掩码后由 GT 重算
max_missing = 0.10
min_maf = 0.05
min_mac = 0                        # 0 = 关闭
min_hwe_p = 0.0                    # 0 = 关闭（自交物种默认关闭）
hwe_mode = "exact"                 # exact | midp

[recompute]
tags = ["AC", "AN", "NS", "MAF", "HWE"]
drop_stale = ["AC", "AN", "AF", "MLEAC", "MLEAF", "ExcessHet", "InbreedingCoeff", "DP"]

[policy]
ploidy = "require_diploid"
multiallelic = "reject"
```

- 比较式 `keep_if`/`when` 解析为 `Cmp { op: Ge | Gt | Le | Lt | Eq | Ne, value: Typed }`。
  这只是一个固定格式的二元比较，不是表达式语言。
- Number=A 且多等位时必须声明量词 `quantifier = "any" | "all"`，否则在 plan 编译期报错。
- **校验在 plan 编译期完成**（需要输入 header）：标签已声明、类型是数值、Number ∈ {1, A}；
  阈值范围合法（MAF ∈ [0, 0.5] 等，沿用现有 `build()` 的校验）。
- `compat::from_gold_conf` 保留现有 `filters.conf` 的 KEY=VALUE 读入，只服务 gold 对拍，并保证与 TOML 等价（测试锁定）。
- 群体维度不进入过滤器（沿用 filter-architecture.md §8）。

### 11.2 谓词分阶段与求值序【规】

| 阶段 | 输入 | 规则 | 求值序 |
|---|---|---|---|
| Cohort | 全基因组原始 GT 缺失 | `max_missing` / keep/remove 列表 | 全局归约，先于一切 |
| Head | 固定列 + INFO | targets → variant_class/max_alt → FILTER → QUAL → INFO 规则（按声明序） | 声明序，首个失败即短路 |
| Mask | FORMAT 列 | `mask_if`（任一为真） | 逐样本 |
| Tail | 掩码后计数（只计 keep 样本） | missing → MAF → MAC → HWE | 固定序（HWE 最贵，放最后） |

- 漏斗语义：第 i 条规则的计数 = "通过前 i−1 条且通过第 i 条"的记录数。因此求值序是语义的一部分，引擎不得重排。
- 这一顺序与 gold（bcftools 分支）逐位对齐，沿用 filter-architecture.md §4。

### 11.3 计划编译器

```rust
pub fn compile(spec: &FilterSpec, header: &HeaderFacts, budget: &PlanBudget) -> Result<FilterPlan, FilterError>;
pub struct FilterPlan {
    pub head: Vec<HeadRule>,  pub mask: Option<MaskRule>,  pub tail: Vec<TailRule>,
    pub cohort: Option<CohortRule>,
    pub field_plans: PassFieldPlans,      // 每个 pass 需要解码的列
    pub strategy: ExecStrategy,           // §11.4
    pub edit_template: EditTemplate,      // 陈旧标签集合、重算标签与格式、header 变更
    pub ascertainment: Ascertainment,     // §11.7
    pub provenance: Provenance,           // §11.8
}
```

现有 `exec.rs:473` 的 `let _plan = plan::compile(spec);`（计算了却不用）在目标态中不存在：
plan 是执行的唯一输入。

### 11.4 三种执行策略与成本模型【规】

| 策略 | 触发条件 | 解析次数 | 说明 |
|---|---|---|---|
| `SinglePass` | 无 cohort 规则，或 cohort 为显式 keep/remove 列表 | 1 次完整 | 头部 → 掩码 → 尾部 → 编辑 → 压缩，全部在 worker 内 |
| `SketchThenEmit` | 有 `max_missing` 规则，**且无尾部表达式/`-q` 的 INFO 依赖**（中段没有记录字节；doc_v12 N1），且 sketch 估计大小 ≤ 内存预算 | 1 次完整 + 1 次"只 inflate、只解析存活者" | 见下 |
| `TwoPassReparse` | sketch 放不下 | 2 次 | pass0 只解 GT 做缺失归约；pass1 = SinglePass |

**SketchThenEmit 细节**：

- **pass0（worker）**：对所有记录解 GT（原始缺失按样本累加，即 cohort 的口径，覆盖全部记录）；
  对通过头部规则的记录再解掩码字段，应用掩码，保存**掩码后位平面**和**掩码位图**（每位点 W 个 word），
  以及头部规则通过位图。sketch 按 chunk 存放。
- **中段（并行，rayon，同一预算）**：`keep = 原始缺失率 ≤ 阈值`；对每个 chunk 的 sketch 用
  `keep ∧ 平面` 做 popcount 得到尾部计数 → 尾部规则 → 存活位图，以及存活者的重算值（AC/AN/NS/MAF/HWE，很小）。
- **pass1（worker）**：再次 inflate 自己的 chunk；非存活记录只用 memchr 跳过换行（不解析任何字段）；
  存活记录按 sketch 中的掩码位图和重算值生成 `EditPlan` → LineEditor → deflate。
- **sketch 大小**：每条头部存活记录 `(4 + 1) × W × 8` 字节。本机：5 × 6 × 8 B × 3.3M ≈ 790 MB（含掩码位图），可以接受。
  记录数可以从 CSI 伪 bin 的 `n_mapped` 估计，没有索引时用前几个 chunk 的 "记录/压缩字节" 外推。
  运行中超出预算时，按 chunk 溢写到 `--tmp-dir`，不会失败。
- **本数据的特例**：cohort 最终保留了 348/348 个样本时，中段的计数等于 pass0 已经算好的计数，可以直接复用。

**等价性**【规】：三种策略对同一输入必须产出逐字节相同的输出与漏斗（`tests/plan_strategy.rs` + gold），
策略选择只影响速度。隐藏参数 `--strategy` 仅供测试使用，所选策略写进 provenance。

### 11.5 漏斗与报告

```rust
pub struct FunnelStage { pub rule: RuleId, pub n_passed: u64 }
pub struct FilterReport {
    pub n_input_records: u64, pub n_output_records: u64,
    pub n_samples_in: u32, pub n_samples_out: u32,
    pub funnel: Vec<FunnelStage>,
    pub n_masked_gt: u64, pub masked_by_rule: Vec<(RuleId, u64)>,    // 掩码归因（新增）
    pub policy_counters: PolicyCounters,                            // 单倍体/多等位等策略触发次数
    pub strategy: ExecStrategy, pub elapsed: Duration, pub metrics: RunMetrics,
}
```

- `RuleId` 是稳定标识（如 `site.info.QD`、`after_mask.maf`），JSON 中按字符串序列化（修复 A24）。
- 报告写到 stderr（人读表格）；`--report report.json` 输出 JSON。

### 11.6 重算语义修正

- 允许多等位输入（`multiallelic != reject` 或 `max_alt_alleles > 1`）时，AC 按 Number=A 逐 ALT 计数，
  MAF 按"第二常见等位"口径计算（修复 A7）。
- 补充的 header 定义使用规范 Number（AC/AF 为 `A`），不再是 `Count(1)`。
- 退化位点的 HWE 哨兵沿用 fill-tags 口径（p = 1）。

**PASS 声明**【规】（doc_v10 N7）：输出 header 必须声明
`##FILTER=<ID=PASS,Description="All filters passed">`，输入未声明时补上
（htslib 读写两侧都注入，`vcf.c:1485-1490`/`1724-1726`）。否则 `FILTER=PASS`
的记录在 BCF 编码时字典查不到，且写出的 BCF 头缺少 htslib 假设的第 0 位。

### 11.7 Ascertainment v2 自动推导与回写

| Spec 中出现 | 推导 |
|---|---|
| 任一位点规则（head 或 tail） | `sites = Filtered { plan_hash }` |
| `variant_class` 排除 `ref`，或 `min_maf > 0` / `min_mac > 0` | `invariants = Absent`（或在输入含不变位点时为 `PartiallyRemoved`） |
| `mask_if` 非空 | `genotypes = Masked { plan_hash }` |
| cohort 规则实际删除了样本 | `cohort = Subset { keep_hash, n_samples }` |

- 推导结果写进输出 header：`##reed_ascertainment=<紧凑 JSON>`。
- `reed stats` 读到这一行时自动采用（用户显式声明的清单优先，两者冲突时报错）。
  这样"过滤 → 统计"之间的清单声明与数据事实同源，不再依赖人记忆（DESIGN §11 L1.5 的初衷在此闭环）。

### 11.8 Provenance

```
##reed_version=0.3.0
##reed_command=reed filter --spec filter.toml -o out.vcf.gz --threads 56 clean.vcf.gz
##reed_filter_spec=<规范 JSON：键排序，浮点用最短往返表示>
##reed_filter_sha256=<sha256(规范 JSON || reed 版本)>
##reed_ascertainment=<JSON>
```

替换基于 `Debug` 输出的哈希（`plan.rs:48`，A19）。

---

## 12. 输出格式与产物

### 12.1 统计输出（长表，tidy）【规】

`<prefix>.stats.tsv`，每行对应一个 (窗口, 统计量, 群体组合)：

| 列 | 含义 |
|---|---|
| `chrom` | contig 名 |
| `start` / `end` | 窗口 **1-based 闭区间**（与 pixy/vcftools 一致）【决 D9】 |
| `level` | `window` / `contig` / `genome` |
| `stat` | 稳定 id，如 `pi`、`dxy`、`fst`、`pbs`、`theta_w`、`tajima_d` |
| `estimator` | 如 `per_pair_comparable`、`hudson`、`projected_n600` |
| `group_a` / `group_b` / `group_c` | 群体名（不适用时为空） |
| `value` | 估计值；分母为空时为 `NA` |
| `numerator` / `denominator` | 计数优先的分子/分母（有定义时） |
| `n_sites` / `n_sites_used` / `n_sites_skipped` | 位点计数 |
| `bias` | `BiasCode`（无偏时为空） |
| `se` / `z`（可选） | jackknife |

- 群体名永远不编码进列名（不出现 `pi_pop1` 这类列）。需要宽表时由 `--wide` 在输出层透视。
- `<prefix>.stats.json` sidecar：provenance、请求、清单、坐标约定（`"coordinates": "1-based-closed"`）、
  `BiasCode` 文案表、运行指标。

### 12.2 坐标约定汇总【规】

| 场合 | 约定 |
|---|---|
| 内部计算、`reed-core` 类型 | 0-based 半开（`Pos0`/`Interval0`） |
| VCF/BCF 记录、错误信息、CLI region 参数 | 1-based（VCF POS / `chr:beg-end` 闭区间） |
| TSV 输出 | 1-based 闭区间，列名 `start`/`end`，sidecar 声明 |
| BED 输入/输出 | 0-based 半开（BED 规范） |

### 12.3 过滤输出

VCF/VCF.gz/BCF + 可选索引；`--report` JSON；stderr 漏斗表。header 带 §11.8 的 provenance 行。

### 12.4 QC 输出

见 §10.3。

### 12.5 产物命名

沿用 `docs/spec/naming.md §6.3` 的语法 `<base>.<stage>[.<qualifier>][.<executor>].<format>`
以及 §6.4 的冻结布局，全部有效。新增：`<base>.stats.tsv`、`<base>.stats.json`、`<base>.qc.json`、
`<base>.qc.<section>.tsv`、`<base>.filter.report.json`、`<base>.run.profile.json`。

---

## 13. 命名规范（终版）

> 本节是 `docs/spec/naming.md` 的超集。相同的规则不再重复解释，冲突处以本节为准并标注 ⚠。

### 13.1 Crate 与包【规】

- 包名 `reed-<layer>`（kebab-case），库名自动为 `reed_<layer>`。门面包名 `reed`，二进制名 `reed`。
- 层名只用以下集合：`core`、`io`、`engine`、`stats`、`qc`、`filter`、`cli`。新增层要先写 ADR。
- 全部 crate 版本同步（lockstep）【决 D11】，由 `[workspace.package] version` 统一给出。

### 13.2 文件与目录【规】

| 对象 | 规则 | 例 |
|---|---|---|
| Rust 源文件 | `snake_case.rs`，单数名词，按**模块概念**命名而不是按类型名 | `plan/strategy.rs` 里放 `ExecStrategy` |
| 目录模块 | `<concept>/mod.rs` + 每个后端/变体一个文件 | `vcf/{mod,line,fixed,info,sample_scan,edit}.rs` |
| 估计器 | 一族一个文件，文件名 = 统计量短名 | `fst.rs`、`tajima.rs`（⚠ 现 `spectrum.rs` 拆为 `sfs.rs`/`theta.rs`/`tajima.rs`） |
| 集成测试 | `tests/<行为或契约>.rs` | `tests/chunk_seams.rs`、`tests/accumulator_laws.rs` |
| 基准 | `benches/<被测对象>.rs` | `benches/sample_scan.rs` |
| fuzz 目标 | `fuzz/fuzz_targets/<对象>.rs` | `fuzz_targets/sample_scan.rs` |
| fixture | `tests/fixtures/<场景>.<格式>`，场景用 snake_case | `half_calls.vcf`、`contig_revisit.vcf`、`multiallelic.bcf` |
| 文档 | `docs/<kebab-topic>.md`；ADR `docs/adr/NNNN-<kebab-title>.md` | `docs/adr/0001-hand-written-sample-scanner.md` |
| 脚本 | `<verb>_<object>.<sh\|py>`（沿用 §6.2） | `check_parity.py` |
| 配置样例 | `examples/<用途>.toml` | `examples/filter.gatk_hard.toml` |

禁止的文件名：`utils.rs`、`helpers.rs`、`common.rs`、`misc.rs`、`types.rs`（没有概念的"杂物箱"）。

### 13.3 模块【规】

- 小写单数。函数名不重复模块名（`fst::fst` 这种口吃，⚠ 在 reed-stats 中一次性消除：
  调用形式为 `reed_stats::fst(...)` 或 `Fst::hudson(a, b)`）。
- `pub use` 只在 crate 根和 `prelude` 里做，避免同一类型有多条公开路径。

### 13.4 领域词汇表（一词一义，一义一词）【规】

| 概念 | 唯一用词 | 禁用的同义词 |
|---|---|---|
| 样本级保留（cohort） | `keep`（`SampleMask`） | retain、select、include（作样本义时） |
| 基因型置缺失 | `mask` | filter、hide、blank |
| 输出中删除样本列 | `prune` | drop、remove（作列义时） |
| 位点被规则淘汰 | `reject` / 通过为 `pass` | drop、filter_out |
| 删除 INFO 标签 | `remove`（标签）/ `stale`（陈旧标签集合） | strip、clean |
| VCF 一行 | `record` | line（仅指字节层的 `LineRef`）、row、entry |
| 基因组上的一个位置（统计意义） | `site` | locus（`Locus` 只用于错误定位）、snp、variant |
| 个体 | `sample` | individual、indiv（BCF 规范术语只出现在 BCF 解码内部） |
| 群体 | `group`（类型）/ population（文档用语） | pop（标识符中禁用）、cohort（专指过滤后的样本全集） |
| 二倍体的两个等位位置 | `slot`（`slot0`/`slot1`） | hap、strand |
| 统计窗口 | `window` | bin（专指分箱）、region（专指用户区间） |
| 窗口分箱 | `bin` | bucket（专指直方图） |
| 并行处理单位 | `chunk` | batch（⚠ 退役）、block（专指 BGZF 块）、shard（专指区域分片） |
| 直方图的格 | `bucket` | bin |
| 并行部分结果 | `partial` | sub_result、piece |
| 累加器 | `acc` / `Accumulator` | state、agg |
| 合并两个累加器 | `merge` | combine、join、reduce（专指全局归约 pass） |
| 全局归约 pass | `reduce` | pass0（标识符中禁用，仅文档映射 gold） |
| 字节 → 字段 | `parse`（固定列/INFO）、`scan`（样本区单遍扫描） | read（专指 I/O） |
| 字段 → 位平面/列 | `decode` | convert |
| 声明 / 编译产物 / 运行结果 / 输入旋钮 | `Spec` / `Plan` / `Report` / `Options` | Config（只给 CLI 文件层）、Summary、Settings |

### 13.5 类型命名【规】

| 后缀/形式 | 含义 | 例 |
|---|---|---|
| `*Spec` | 用户声明（值类型，serde） | `FilterSpec` |
| `*Request` | 统计请求集合 | `StatsRequest` |
| `*Plan` | 编译后的执行计划 | `FilterPlan`、`FieldPlan`（⚠ 取代 `ExecutionPlan`/`Plan`） |
| `*Options` | 资源/行为旋钮，无副作用 | `ReadOptions` |
| `*Report` | 一次运行的结果与指标 | `FilterReport`、`RunReport` |
| `*Error` | 每个 crate 一个错误枚举 | `IoError`、`EngineError` |
| `*Id` | 稠密下标 newtype | `ContigId`、`TagId` |
| `*Mask` | 按样本的位集合 | `SampleMask` |
| `*Table` | id ↔ 名字 ↔ 属性的查找表 | `ContigTable`、`GroupTable` |
| `*Counts` | 整数计数结构 | `SiteGroupCounts` |
| `*Acc` | 累加器 | `HudsonFstAcc` |
| `*Output` | 统计量 `finish` 的结果 | `FstOutput`（⚠ 新代码不用 `*Result`，避免与 `std::result::Result` 混淆；旧类型在兼容层保留） |
| `*Ref<'a>` | 借用的零拷贝视图 | `LineRef`、`BcfRecordRef` |
| `*Block` | 一个 chunk 的列式数据 | `SiteBlock`、`PlaneBlock` |
| `*Reader` / `*Writer` | I/O 类型（`Source` 禁用，沿用） | `BgzfWriter` |
| `*Policy` | 对不支持输入的显式处理方式 | `PloidyPolicy` |
| `*Rule` | 过滤规则 | `HeadRule`、`MaskRule` |
| `*Estimator` | 封死定义的变体枚举（沿用） | `PiEstimator` |
| trait | 角色名词或能力（`Accumulator`、`Statistic`、`ChunkJob`、`OrderedSink`） | 不加 `Trait`/`I` 前缀 |
| 缩写 | 在类型名中按单词驼峰 | `VcfBgzf`、`HweAcc`、`LdPair`、`Fst`（不写 `FST`、`HWE`） |

### 13.6 函数与方法【规】

| 形式 | 用途 | 例 |
|---|---|---|
| `new` / `with_*` / `from_*` / `try_from_*` | 构造；会失败的返回 `Result` | `SampleMask::from_ids` |
| 无 `get_` 前缀的访问器 | 字段访问 | `group.name()` |
| `n_*()` | 计数访问器 | `mask.n_selected()` |
| `is_*` / `has_*` | 布尔 | `key.is_snv()` |
| `as_*`（廉价借用）/ `to_*`（分配）/ `into_*`（消耗） | 转换（Rust API Guidelines） | `pos1.to_pos0()` |
| `parse_*` / `scan_*` / `decode_*` | 见 §13.4 | `scan_samples` |
| `for_each_*` | 推式流（回调） | `for_each_record` |
| `iter*` | 拉式迭代器 | `table.iter()` |
| `push_*` / `merge` / `finish` | 累加器生命周期 | `acc.push_site(&input)` |
| `compile` | Spec → Plan | `filter::compile` |
| `run` / `run_*` | 驱动一次执行（L4/L2） | `reed::run_filter` |
| `write_*` / `emit` | 序列化 / sink 角色 | `write_header` |
| `try_*` | 与同名不失败版本并存的失败版本 | — |

### 13.7 变量命名总表【规】

**通用规则**：

- snake_case；名字要说明**单位或基准**（`_bp`、`_bytes`、`pos0`/`pos1`），而不是类型。
- 单字母只允许：≤ 10 行循环里的 `i/j/k`；数学公式里与论文符号一一对应的量（`n`、`p`、`s2`、`a`/`b`/`c`），
  且必须在函数文档中给出符号表；闭包参数 `|e|`（错误）。
- 公共 API 和字段中禁用缩写（沿用 §2.10 的禁用表 `buf/rec/iv/idx/src/tmp/cnt`），例外如下：
  ⚠ `*_idx` 作为后缀允许出现在局部变量和私有字段里（如 `site_idx`）；裸 `idx` 仍然禁用。
- VCF 标签镜像结构（字段与 VCF 标签一一对应，如 `InfoTags { ac, an, ns, maf, hwe }`）使用标签小写名；
  分析代码使用描述性名字（`n_alt_alleles`，而不是 `ac`）。

**领域变量**：

| 概念 | 规范名 | 反例（禁用） |
|---|---|---|
| 样本数 / 位点数 / 记录数 / 群体数 | `n_samples` / `n_sites` / `n_records` / `n_groups` | `nsamp`、`num_sites`、`site_cnt`、`npop` |
| 组内等位基因拷贝上限（2·\|G\|） | `n_allele_slots` | `n_chroms`、`n_chrom`、`chrom_a`（⚠ 旧名退役） |
| 已调用的等位基因拷贝数 | `n_called_alleles` | `n`、`an`（标签镜像除外）、`n_a` |
| REF/ALT 等位基因计数 | `n_ref_alleles` / `n_alt_alleles` | `ref_count`、`alt`、`ac` |
| 基因型类计数 | `n_hom_ref`、`n_het`、`n_hom_alt`、`n_full_called`、`n_half_called`、`n_missing_gt` | `het`（裸）、`hom1`/`hom2`、`obs_homr` |
| 频率 | `p_alt`、`p_ref`、`p_mean`、`maf` | `freq`、`af`（标签镜像除外）、`p1`/`p2`（公式内且有符号表时允许） |
| 平均成对差异 | `mpd_within_a`、`mpd_between` | `pa`、`pab`、`wa` |
| p 值 | `p_value`、`hwe_p_value` | 裸 `p`（与频率冲突） |
| 比例（0–1） | `*_rate` / `*_frac` | `*_pct`（只能用于 0–100） |
| 位置 | `pos0`、`pos1` | 未注明基准的 `pos`、`p` |
| 区间端点 | `start0`/`end0`（半开）、`start1`/`end1`（闭） | `beg`、`s`/`e` |
| 跨度 | `span_bp`、`window_bp` | `len`（与集合长度冲突） |
| 下标（局部） | `site_idx`、`sample_idx`、`group_idx`、`word_idx`、`bin_idx` | `s`、`ix`、`index`（作为下标变量） |
| 强类型 id | `contig_id`、`sample_id`、`group_id`、`chunk_id`、`tag_id` | `contig`（当它是 id 时） |
| 名字 | `contig_name`、`sample_name`、`group_name` | `name`（跨多种实体的函数里） |
| 样本保留掩码 | `keep: &SampleMask` | `keep: &[bool]`、`retain` |
| 待掩码基因型 | `mask_gt: Bitset` | `masked: Vec<bool>` |
| 群体掩码 | `group_mask` | `g`（公式外） |
| 位平面 | `slot0_alt`、`slot0_miss`、`slot1_alt`、`slot1_miss`、`phased` | `a0`、`m1`（公式文档内允许） |
| word 数 | `n_words` | `w`（内核里有符号表时允许） |
| 累加器 | `<stat>_acc`、`accs` | `state`、`agg` |
| 部分结果 | `partial`、`chunk_partial` | `sub`、`piece` |
| 统计输出 | `<stat>_output` | `res`、`r` |
| 时间 | `started_at: Instant`、`elapsed: Duration`、`<stage>_busy` | `t0`、`t_recv` |
| 字节量 | `n_bytes`、`compressed_bytes`、`<purpose>_bytes` | `size`（含义不明） |
| 通道 | `<payload>_tx` / `<payload>_rx` | 单个通道以外的裸 `tx`/`rx` |
| worker | `n_workers`、`worker_idx` | `w`、`t` |
| 临时缓冲 | `<purpose>_scratch`（如 `inflate_scratch`） | `buf`、`tmp` |
| 阈值 | `min_<x>`、`max_<x>`（含义为闭区间：`keep if x >= min_x`） | `threshold`、`t`、`m`、`cut` |
| 布尔 | `is_*`、`has_*`、`require_*`、`should_*` | `flag`、`b`、`ok`（`ok` 只允许在 ≤ 5 行的作用域） |
| 错误 | `err`（所有权变量）、`e`（闭包） | `error1`、`ex` |

### 13.8 常量与静态【规】

- `SCREAMING_SNAKE_CASE`，名字包含单位：`CHUNK_TARGET_COMPRESSED`（字节）、`BGZF_MAX_BLOCK_DATA = 65_280`、
  `TICKETS_PER_WORKER`、`HWE_MEMO_CAPACITY`。
- 用户可覆盖的默认值统一叫 `DEFAULT_<X>`（如 `DEFAULT_WINDOW_BP`）。
- 参与确定性定义的常量（chunk 大小等）在文档注释中标注 `// determinism-relevant`。

### 13.9 泛型与生命周期【规】

- 类型参数用有意义的单字母或短名：`R: Read`、`W: Write`、`S: Statistic`、`A: Accumulator`、`J: ChunkJob`。
  超过两个时使用完整名（`Job`、`Sink`）。
- 生命周期：单个时用 `'a`；多个时按语义命名：`'hdr`（header）、`'line`（行字节）、`'blk`（chunk 块）。

### 13.10 Feature flags【规】

kebab-case，名词或能力，不用 `use-`/`enable-` 前缀：`libdeflate`（默认开）、`simd`（默认开，运行时探测）、
`serde`、`mimalloc`（仅 reed-cli）、`python`（远期）、`parquet`（远期）。

### 13.11 错误变体【规】

- 按**失败模式**命名（沿用 §2.7）：`UnsupportedInput`、`HeaderMismatch`、`PositionNotAscending`、
  `ContigRevisit`、`UndeclaredTag`、`MemBudgetExceeded`、`WorkerPanicked`。
- 能定位的都携带 `Locus { contig_name, pos1 }` 或 `ChunkId`。
- ⚠ `ReadError::User`/`Core`、`FilterError::Format` 按旧台账更名，并在迁移中消失。

### 13.12 测试命名【规】

- 测试函数名是行为句：`<subject>_<expected_behaviour>[_when_<condition>]`，例如
  `sample_scan_matches_noodles_on_half_calls`、`sketch_then_emit_equals_two_pass_reparse`、
  `fst_matrix_equals_pairwise_calls`。
- 每个测试模块最多一个 `fixture()`；参数化构造器为 `fixture_*`/`make_*`；文件路径辅助函数为 `fixture_path(name)`（沿用 §2.11）。
- proptest 策略函数：`arb_<type>()`（如 `arb_genotype_planes()`）。

### 13.13 基准命名【规】

criterion 的 id 采用 `<crate>/<subject>/<parameter>`，例如 `core/group_counts/s348_g3`、
`io/sample_scan/gt_dp_gq`、`io/inflate/libdeflate_1mib`。

### 13.14 CLI 参数【规】

- 长参数 kebab-case；短参数**只用于与 bcftools 语义相同的字母**，避免熟悉 bcftools 的用户被误导：
  `-o/--output`、`-O/--output-type`（`v|z|b|u`）、`-r/--regions`（需要索引）、`-R/--regions-file`、
  `-t/--targets`（流式）、`-T/--targets-file`、`-s/--samples`、`-S/--samples-file`、`-W/--write-index`。
- ⚠ `--threads` 没有短参数（bcftools 的 `-t` 是 targets）。
- 数值带单位：`--window 100kb`、`--mem 16G`（解析器统一支持 `bp/kb/Mb/Gb` 与 `K/M/G`）。
- 子命令为单个动词或名词：`filter`、`stats`、`qc`、`index`。

### 13.15 配置键（TOML）【规】

snake_case；阈值 `min_*`/`max_*`；列表用复数（`tags`、`include`）；带单位的键名带后缀（`window_bp`、`max_ld_bp`）；
枚举值为 snake_case 字符串（`"require_diploid"`）。

### 13.16 输出列【规】

snake_case，见 §12.1；`n_*` 表示计数；坐标列固定为 `chrom/start/end`；群体列固定为 `group/group_a/group_b/group_c`。

### 13.17 漏斗规则 id【规】

点分小写层级 + 原始标签名：`site.variant_class`、`site.filter`、`site.qual`、
`site.info.QD`、`genotype.mask.DP`、`after_mask.missing`、
`after_mask.maf`、`after_mask.hwe`；样本级规则用 `sample.mind`
（doc_v4 02-cli-semantics.md §5.3 定名；首版 `after_mask.max_missing`/`sample.max_missing` 已撤回）。

### 13.18 Git【规】

- 分支：`dev`（沿用）、`feat/<crate>-<topic>`、`fix/…`、`perf/…`、`refactor/…`。
- 提交：Conventional Commits，scope 用 crate 短名：`perf(io): block-aligned parallel parsing`、
  `fix(stats): projected sfs weight underflow`。
- 标签：`v<semver>`；CHANGELOG 采用 Keep a Changelog 格式。

### 13.19 重命名台账（主要项；完整映射见附录 A）

| 现名 | 目标名 | 理由 |
|---|---|---|
| `WindowBuffer` | 由 `SiteBlock` + 分箱累加取代；保留为兼容适配器一个大版本 | 表示与算法变更（§5、§8） |
| `PopulationView` | `GroupTable` | 命名群体、掩码表示 |
| `Bits` | `Bitset` | 与 `SampleMask` 等组合时含义更清楚 |
| `SiteCounts`（reed-core estimators） | `SiteGroupCounts`（core::count）+ reed-filter 的 `InfoTags` 视图 | 依赖方向（A10） |
| `RecordView` | 退役：位点字段用 `LineRef`/`BcfRecordRef` 惰性访问，样本数据进 `SiteBlock` | 不再物化 per-record 对象 |
| `GtMask` | `GenotypeMaskRule`（reed-filter） | 策略回到策略层（A9） |
| `Recompute` | `RecomputeTags`（reed-filter） | 同上 |
| `STALE_INFO` | `DEFAULT_STALE_TAGS`（reed-filter） | 同上；常量命名规则 |
| `record::apply` | `reed_io::vcf::edit_vcf_line` + `EditPlan` | 字节级编辑 |
| `VariantSink` / `sink::open` | `OutputSpec` + `BgzfWriter`/`PlainWriter`（由 L4 选择） | 写出改为有序 append |
| `ExecOptions.threads` | `ThreadBudget` | 线程预算统一 |
| `RunReport.funnel: Vec<(String,u64)>` | `Vec<FunnelStage>` | A24 |
| `ReadError` | `IoError` | 每个 crate 一个 `*Error`，以 crate 命名 |
| `spectrum.rs` | `sfs.rs` / `theta.rs` / `tajima.rs` | 一族一个文件 |
| `allele_counts.rs` / `AlleleCountsResult` | `allele_freq.rs` / `AlleleFreqOutput` | 与 `SiteGroupCounts` 拉开距离 |
| `BATCH_SIZE`（4096 条） | `CHUNK_TARGET_COMPRESSED`（字节） | 按字节而非条数切分 |

---

## 14. 代码规范

### 14.1 语言版本

- `edition = "2024"`，`rust-version = "1.91"`（noodles 0.117 的 MSRV；本机 1.99）【决 D10】。
- `rust-toolchain.toml` 固定为本机工具链，CI 额外跑一次 MSRV 构建。

### 14.2 Workspace 清单（目标形态）

```toml
[workspace]
resolver = "3"
members = ["crates/*", "xtask"]

[workspace.package]
version = "0.3.0"
edition = "2024"
rust-version = "1.91"
license = "MIT OR Apache-2.0"

[workspace.dependencies]
noodles = { version = "0.117", features = ["bcf", "bgzf", "core", "csi", "tabix", "vcf"] }
libdeflater = "1"
memchr = "2"
memmap2 = "0.9"
crossbeam-channel = "0.5"
rayon = "1"
smallvec = "1"
bitflags = "2"
thiserror = "2"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
toml = "0.8"
sha2 = "0.10"
ryu = "1"
tracing = "0.1"
# 只给 cli：
clap = { version = "4", features = ["derive"] }
tracing-subscriber = { version = "0.3", features = ["env-filter", "json"] }
mimalloc = "0.1"
anyhow = "1"
# dev：
proptest = "1"
criterion = "0.5"
assert_cmd = "2"
```

浮点解析直接用标准库 `str::parse::<f64>`（Eisel-Lemire 算法，已经足够快），不引入 fast-float。

### 14.3 Lints 与 profile

```toml
[workspace.lints.rust]
unsafe_code = "deny"            # 例外模块显式 #![allow(unsafe_code)]，每个 unsafe 块写 // SAFETY:
missing_docs = "warn"
unreachable_pub = "warn"

[workspace.lints.clippy]
all = { level = "deny", priority = -1 }
pedantic = { level = "warn", priority = -1 }
unwrap_used = "deny"            # clippy.toml: allow-unwrap-in-tests = true
expect_used = "warn"
panic = "warn"
cast_possible_truncation = "warn"   # 抓 `as u32` 截断（A12）
cast_sign_loss = "warn"
module_name_repetitions = "allow"
must_use_candidate = "allow"
missing_errors_doc = "allow"

[profile.release]
lto = "fat"
codegen-units = 1
debug = "line-tables-only"      # perf/flamegraph 可读，体积影响很小
panic = "unwind"                # 引擎需要 catch_unwind

[profile.dev.package."*"]
opt-level = 2                   # 依赖在 dev 下也优化，测试大文件时不至于太慢

[profile.profiling]
inherits = "release"
debug = true
```

- 发布构建不用 `target-cpu=native`；SIMD 内核做运行时分派（`is_x86_feature_detected!`）。
  本地实验可以在 `.cargo/config.toml` 打开 native，但基准日志必须注明。

### 14.4 unsafe 政策

- 只允许出现在：`memmap2` 调用点（`reed-io/src/chunk/bgzf.rs`，doc_v6 06 Y0 起是工作区唯一一处）。
  doc_v6 Y0 删除了 `reed-core::count::kernel` 的 AVX-512 内核（release 下长度检查不成立 = 可达 UB，
  且 348 样本下永不向量化）；若将来重新引入 SIMD，必须先恢复本文档与 `ci/check-doc-drift.sh` 的计数守卫。
- 每个 unsafe 块都要写 `// SAFETY:` 说明前置条件；标量路径跑 Miri；解析器跑 fuzz。

### 14.5 API 规则

- 带不变量的类型不公开字段（现 `WindowBuffer` 字段全 pub，可以构造出非法状态）：私有字段 + 访问器 + 校验构造器。
- 可能扩展的公开枚举和结构体加 `#[non_exhaustive]`。
- 超过 3 个参数的构造使用 builder 或参数结构体（取代 `#[allow(clippy::too_many_arguments)]`）。
- 库不 panic：输入数据导致的任何问题都返回 `Err`。`expect` 只用于内部不变量，消息以 `"invariant: "` 开头。
- 库不打印：不出现 `println!`/`eprintln!`；诊断信息走 `tracing`。
- 热路径禁止：逐记录 `String` 分配（contig 名先查成 `ContigId`，带上次命中缓存）、逐记录 `Arc` clone、
  逐记录跨线程发送、`HashMap` 查找（用数组下标）、`Vec<bool>`（用 `Bitset`）、能用 popcount 的地方写逐样本 f64 循环。
- 每 worker 的计数器结构体按 128 B 对齐填充，避免伪共享。

### 14.6 依赖政策

- 新增依赖需要在 PR 中说明理由；`cargo-deny` 白名单许可证为 MIT/Apache-2.0/BSD/ISC/Zlib/Unicode；
  GPL/LGPL 禁止（与 DESIGN §8 的翻译政策一致：GPL 项目只读，不翻译、不链接）。
- L3 策略 crate 只能依赖 reed-core、serde、thiserror、smallvec、bitflags（保持纯）。

### 14.7 文档与注释

- 代码注释和 rustdoc 用英文（沿用现状）；设计文档用中文（沿用现状）。
- 每个模块的文档说明：用途、不变量、引用（论文与参考实现的 `file:line`）。
- 估计器文档包含符号表（单字母变量与论文符号的对照）以及公式。
- 新的工程决策写 ADR（§16.7），不再追加进 docs/rounds/v1/design.md。

---

## 15. 错误、日志与可观测性

- 每个库 crate 一个 `thiserror` 错误枚举，`#[non_exhaustive]`；`source()` 链完整。
- 只有 reed-cli 使用 `anyhow`。退出码：0 成功；1 运行期失败；2 用法错误；3 输入与 Spec/header 不一致
  （例如过滤谓词引用了未声明的 INFO 标签）。
- `tracing`：`error`（中止原因）、`warn`（策略降级，例如普通 gzip 只能串行解压）、`info`（阶段开始与结束、漏斗）、
  `debug`（chunk 级）、`trace`（行级，默认编译进去但不开启）。
- 进度条只在 stderr 是 TTY 时显示，按压缩字节进度计算。
- `--profile run.json`：输出各阶段 `StageMetrics`、内存峰值、所选策略，用于 §17.3 的基准日志。

---

## 16. 测试与验证体系

### 16.1 测试金字塔【规】

| 层 | 工具 | 对象 | 判据 |
|---|---|---|---|
| 单元 | `#[test]` | 每个函数的边界情况 | 每个错误变体至少有一个触发它的测试 |
| 性质 | proptest | SIMD 与标量内核；累加器三定律；LineEditor（编辑后再解析 = 预期）；分箱 = 直接窗口 | 逐位一致 / ≤ 4 ulp |
| 差分 | 自写 harness | 快速样本扫描 vs noodles；新统计 vs **冻结的 reed-core 0.2 实现** | 逐位一致（整数），≤ 1e−12（浮点） |
| 正典对拍 | Python 生成的 fixture | allel/sgkit/pixy/easySFS（沿用现有 families/unequal_groups） | 同上 |
| 三方金标准 | `gold/` | 过滤：bcftools/vcftools/plink2；QC：bcftools stats；LD：plink2 `--r2-unphased`（只比行为，不读代码） | L_site/L_gt/L_tags 精确 |
| 确定性 | e2e | threads ∈ {1,2,3,7,56} | 输出 sha256 相同 |
| 策略等价 | e2e | SinglePass / SketchThenEmit / TwoPassReparse | 输出 sha256 相同 |
| 模糊 | cargo-fuzz | sample_scan、line、info、BGZF 块头 | 不 panic，结果与参考一致 |
| CLI | assert_cmd | 参数、退出码、stdout/stderr | 黑盒快照 |
| 性能 | criterion + hyperfine | 内核、扫描、端到端 | §17.4 门 |

### 16.2 冻结旧实现作为 oracle【荐】

M1 开始前，把当前的 `reed-core 0.2.0` 原样复制为 `crates/reed-core-v02-oracle`（`publish = false`，只作为
dev-dependency）。新架构的统计结果在 fixture 和全量数据上与它对拍，直到 M7 有意改变语义为止
（例如投影 SFS 是新估计器，complete-case 估计器仍然要与 oracle 一致）。

### 16.3 Fixture 规则

- 小 fixture 入库：`tests/fixtures/<场景>.<格式>`。每个场景对应一个 DESIGN/本文中的规则：
  半调用、`.` 样本列、整列缺失、单倍体、多等位、ContigRevisit、乱序、超长行（跨 chunk）、
  恰好以 `\n` 结尾的 chunk 边界、VCFv4.4 header、无 GT 列、INFO 缺失标签、Number=A 多值等。
- 接缝测试必须把 `CHUNK_TARGET_COMPRESSED` 调到很小（测试专用构造参数），制造大量接缝。
- 大数据不入库；`xtask fixtures` 负责生成合成数据（随机缺失、不等覆盖、多群体）。

### 16.4 gold 扩展

- 在 `compare_tools.py` 中把 reed 作为第四个执行器（候选者），按 L_site/L_gt/L_tags 对拍。gold 本身仍然只由外部工具生成。
- 新增：QC 与 bcftools stats 的重叠节对拍；投影 SFS 与 easySFS 对拍；Fst/π 与 allel 窗口函数在全量 Chr01 上对拍。

### 16.5 CI（本地脚本或 xtask 也可以）

`fmt --check` → `clippy -D warnings` → `test`（debug 与 release 各一次）→ `miri`（core 标量）→ `cargo deny` →
`doc` → fuzz 冒烟（每个目标 60 s）→ 小规模 e2e 确定性测试。

### 16.6 覆盖要求

每个公开函数至少一个测试；每个错误变体至少一个测试；每个 `Policy` 的每个取值至少一个测试；
每个统计量至少一个正典对拍和一个累加器定律测试。

### 16.7 ADR

`docs/adr/NNNN-<kebab-title>.md`，固定小节：背景 / 决策 / 证据【据】/ 后果 / 被取代的决策。
§20 每个【决】项拍板后各写一篇 ADR。

---

## 17. 性能工程

### 17.1 本机预算（clean.vcf.gz：6.6 GB 压缩，约 35 GB 解压，348 样本 × 3.3M 记录，56 核 2×NUMA）

| 环节 | 单核量级 | 54 worker 理论 | 现实估计 | 依据 |
|---|---|---|---|---|
| 读压缩字节 | — | — | 热缓存约 1–2 s；/sri 冷读 3–7 s | 6.6 GB / 1–5 GB/s |
| inflate（libdeflate） | 约 1 GB/s（输出） | 0.7 s | 3–5 s（受内存带宽限制） | 35 GB |
| 固定列 + GT 扫描 | 约 1 GB/s | 0.7 s | 3–6 s | 35 GB |
| GT + DP + GQ 单遍扫描 | 约 0.5 GB/s | 1.3 s | 5–10 s | 同上 |
| 计数内核 + 累加 | 可忽略 | — | < 1 s | 每位点每群体约 40 次 word 运算 |
| 存活者字节编辑（53%） | 约 1 GB/s | 0.35 s | 1–3 s | 约 18.6 GB 输出文本 |
| deflate level 6（libdeflate） | 约 0.1–0.15 GB/s | 2.3–3.4 s | 4–8 s | 约 18.6 GB |
| **统计（单遍，GT-only）** | | | **≤ 30 s（目标）** | |
| **过滤 SketchThenEmit** | | | **≤ 60 s（目标）** | pass0 约 15–20 s + pass1 约 15–25 s |
| **QC 全报告** | | | **≤ 30 s（目标）** | 与过滤 pass0 同量级 |

对照：bcftools 约 40 min；vcftools 约 35 min；plink2 约 2.5 min（丢 FORMAT，输出只有 0.26 GB）；reed v0.2 约 9.6 min
（文档）/ 2422.8 s（日志，A22）。这些目标是基于上表的预算，不是承诺；每一步都要用 §17.3 的协议实测。

### 17.2 工具

`perf record -g` + `cargo flamegraph`、samply、heaptrack（分配热点）、`perf stat -e cache-misses,…`、
`numastat`、`--profile` JSON。每个 perf 类 PR 附上火焰图前后对比。

### 17.3 基准日志协议【规】

追加写入 `bench/suites/e2e/results.tsv`（只追加不修改），每行包含：`date`、`commit`、`dirty`、`command`、`input`、
`input_sha256_prefix`、`cache_state`（hot/cold，以及是否 drop_caches）、`threads`、`wall_s`、`user_s`、`sys_s`、
`max_rss_mb`、`strategy`、`stage_metrics_path`、`note`。
文档里引用任何性能数字都必须能在这个文件里找到对应的行。这解决 A22 中文档与日志不一致的问题。

### 17.4 回归门

- criterion 微基准：相对上一基线变慢超过 5% 时失败（噪声大的基准标注例外）。
- 子集端到端（`subset.Chr01-03`，48 万记录）：每个 perf/refactor 类 PR 都要跑，结果记入日志。
- 全量端到端：每个里程碑跑一次，作为验收门（§19）。

### 17.5 已知陷阱清单

逐记录的 String/Arc/通道/HashMap；`Vec<bool>`；跨 NUMA 节点的大缓冲由 reader 分配；第二个线程池造成超订；
无界队列；用 `{:?}` 做哈希；浮点归约顺序依赖线程调度；`as` 截断；debug 构建测性能。

---

## 18. CLI 设计（`reed`）

```
reed filter [OPTIONS] <INPUT>
    --spec <filter.toml> | --gold-conf <filters.conf>
    -o, --output <PATH>            默认 stdout 文本 VCF（沿用 Unix 过滤器语义）
    -O, --output-type <v|z|b|u>    默认按扩展名推断
        --compress-level <0-9>     默认 6
    -W, --write-index[=csi|tbi]
    -r/-R, -t/-T, -s/-S            bcftools 同义
        --threads <N>   --mem <SIZE>   --tmp-dir <DIR>
        --report <report.json>   --profile <run.json>   -v/-q

reed stats [OPTIONS] <INPUT>
    --pops <pops.tsv>
    --request <stats.toml> | --stats <list>   如 pi,theta_w,tajima_d,dxy,fst:hudson,pbs:pop1:pop2:pop3
    --window <100kb> [--step <50kb>] | --window-sites <500> [--step-sites <250>] | --windows-bed <BED>
    --accessible <mask.bed>
    --manifest <with-invariants|variant-only:complete|variant-only:filtered|variant-only:unknown>
    --filter <filter.toml>         边过滤边统计，不落中间文件
    --sfs-n <N>                    投影样本量（每群体可单独设置）
    --levels <window,contig,genome>   --per-site   --jackknife-block <5Mb>
    -o, --output-prefix <PREFIX>   --format <tsv|json>   --wide
    --threads / --mem / --profile

reed qc [OPTIONS] <INPUT>
    --pops <pops.tsv>   --sections <list|all>   --info-tags <QD,FS,…|all>
    --format <reed|bcftools-stats>   -o, --output-prefix <PREFIX>
    --filter <filter.toml>           同时报告掩码影响（mask_impact 节）

reed index <INPUT>   [-c|--csi | -t|--tbi]      （薄壳，调用 noodles indexer）
```

- CLI 不允许出现库做不到的能力（原则 1）。CLI 参数先翻译成 `FilterSpec`/`StatsRequest`/`QcRequest`，再调用 L4。
- 现在的 `reed-filter` 二进制在一个版本内保留为 `reed filter` 的别名，然后移除【决 D14】。

---

## 19. 迁移路线（每步都有验收门）

> 原则：先打地基（模型、传输、引擎），再迁移消费者（统计、过滤、QC），最后做算法升级。
> 每一步结束时：全部测试通过、gold 对拍不变、性能数字按 §17.3 入日志。任何一步都可以单独回滚。

| 步 | 内容 | 验收门 |
|---|---|---|
| **M0 卫生（行为不变）** | workspace 统一 package/deps/lints/profile；`tracing` 替换 `TEMP instrumentation`；CLI 改用 clap；按 §17.3 重新测量当前基线（解决 A22） | 测试全绿；gold 对拍不变；基线入日志 |
| **M1 core 模型** | `id/coord/contig/sample/group/bitset/genotype/count/manifest v2/grid/math`；冻结 `reed-core-v02-oracle`；旧 `WindowBuffer` 改为适配器 | proptest：位平面 + 计数内核 ≡ 旧 `allele_counts`/`genotype_counts`；SIMD ≡ 标量 |
| **M2 io 传输** | `Format::detect`；BGZF 切块 + libdeflate；`LineRef`/固定列/INFO；`sample_scan`；`FieldPlan`；BCF 解码 | 差分 ≡ noodles；fuzz 冒烟；**全量 GT-only 扫描 ≤ 20 s（56 线程）** |
| **M3 engine** | `run`、票据窗口、环形重排、`select!` 中止、确定性；第一个消费者：最小 QC（summary 节） | 确定性测试；内存上界测试（在途字节 ≤ 预算）；chunk 接缝测试 |
| **M4 stats v2** | 全部现有估计器改写为累加器；分箱窗口；L4 driver；`reed stats` | 与 oracle 对拍（整数精确、浮点 ≤ 1e−12）；累加器定律；**全量 3 群体统计 ≤ 30 s** |
| **M5 filter v2** | Spec v2 + compat；plan 编译器；SinglePass/TwoPassReparse；LineEditor；worker 内 deflate + 有序写出 + 同步建索引；provenance/清单回写 | gold L_site/L_gt/L_tags；确定性；**≤ 120 s** |
| **M5b SketchThenEmit** | sketch、中段并行判定、pass1 跳过非存活者 | 三策略输出 sha256 相同；**≤ 60 s** |
| **M6 qc** | reed-qc 全部节；bcftools-stats 兼容输出 | 重叠节与 bcftools stats 对拍；**≤ 30 s** |
| **M7 算法升级** | 投影 SFS（默认切换在 easySFS 对拍之后）；f 统计校正 + jackknife；HWE mid-p + 记忆化；位平面 LD + 平面存储 + 衰减 + pruning；WC84/Nei 多群体；Fst 矩阵；联合 SFS；个体 F | 每项都有正典对拍；性能前后对照 |
| **M8 扩展** | 原生 BCF 编辑器；区域分片；HaplotypeStream（EHH/iHS）；Python 绑定（pyo3，可选） | 各自 ADR |

版本：M4 发布 0.3.0（stats）、M5b 发布 0.4.0（filter）；兼容层（`WindowBuffer`、`SiteCompleteness`、`reed-filter` 二进制）在 0.5.0 移除。

---

## 20. 决策登记表（需要你拍板）

| # | 决策 | 推荐 | 备选 | 理由 / 证据 |
|---|---|---|---|---|
| D1 | VCF 样本列/固定列/INFO 用自写快速扫描（推翻 DESIGN【定】"不写解析代码"的这一部分） | **是**，限定在热路径；header、BCF 结构、输出渲染仍交给 noodles；差分测试 + fuzz 兜底 | 继续只用 noodles 惰性 API | noodles 每多解一列 FORMAT +28.3%（§16.2【据】）；单遍多字段扫描是吞吐上限的关键 |
| D2 | crate 拆分为 core/io/engine/stats/qc/filter/reed/cli 共 8 个 | **是** | 合并为 5 个（engine 并入 io，qc 并入 stats，reed 与 cli 合并） | 依赖方向靠 crate 边界强制；L3 保持纯净可单测 |
| D3 | 读写都绕过 noodles bgzf，自建 BGZF 切块 + libdeflate | **是** | 用 noodles `MultithreadedReader/Writer` | MT reader 解压后仍串行；MT writer 造成线程超订（§1.3、A3） |
| D4 | SFS/θ/Tajima 默认使用投影估计器 | **是**（CLI 要求显式给 `--sfs-n` 或使用确定性的默认规则） | 默认 complete-case | 本机零缺失位点 < 1%（profile【据】） |
| D5 | PBS 默认位点集合 | **PerPair**（与 allel 一致，信息利用更多）；Intersection 保留为封死定义 | 维持 Intersection 默认 | DESIGN §11 待定项 |
| D6 | Spec/Request 使用 TOML/JSON（serde），KEY=VALUE 只为 gold 兼容 | **是** | 继续使用 shell conf | 通用 INFO 谓词无法用 KEY=VALUE 表达 |
| D7 | 二进制使用 mimalloc | **是**（先测，有收益才保留） | glibc malloc | 现诊断里的跨线程释放问题在新架构中已经结构性消失，所以这只是保险 |
| D8 | 位置统一使用 u64 | **是** | u32 + 显式报错 | 大染色体；htslib large positions |
| D9 | TSV 输出坐标 1-based 闭区间 | **是**（与 pixy/vcftools 一致） | BED 风格 0-based 半开 | 用户习惯 |
| D10 | edition 2024，MSRV 1.91 | **是** | edition 2021 | 本机 1.99，noodles MSRV 1.91 |
| D11 | 全部 crate 版本同步（lockstep） | **是** | 各自独立版本 | 单人维护，接口同步演进 |
| D12 | Python 绑定 | **M8 之后再做** | 尽早 | 先把 Rust API 稳定下来 |
| D13 | 文本输出不再把 VCFv4.4 降级为 4.3 | **是**（字节直通后已不再重渲染 GT） | 维持降级 | A18 |
| D14 | 单一二进制 `reed` + 子命令 | **是** | 保留独立的 `reed-filter` 二进制 | 统一参数与输出约定 |

---

## 附录 A：旧 → 新文件映射

| 现文件 | 去向 |
|---|---|
| `reed-core/src/bits.rs` | `reed-core/src/bitset.rs` |
| `reed-core/src/window.rs` | `reed-core/src/genotype/*` + `count/*` + `grid.rs`；`WindowBuffer` 改为兼容适配器 |
| `reed-core/src/population.rs` | `reed-core/src/group.rs` |
| `reed-core/src/manifest.rs` | `reed-core/src/manifest.rs`（v2 + 兼容映射） |
| `reed-core/src/result.rs` | 拆入 `reed-stats` 各估计器文件的 `*Output` |
| `reed-core/src/estimators/{pi,dxy,fst,pbs,fstats,hwe,ld}.rs` | `reed-stats/src/{pi,dxy,fst,pbs,fstats,hwe,ld}.rs`（累加器化） |
| `reed-core/src/estimators/spectrum.rs` | `reed-stats/src/{sfs,theta,tajima}.rs` |
| `reed-core/src/estimators/allele_counts.rs` | `reed-stats/src/allele_freq.rs` |
| `reed-core/src/estimators/site_counts.rs` | 计数部分 → `reed-core/src/count/`；MAF/HWE 哨兵视图 → `reed-filter/src/eval/tail.rs` |
| `reed-core/src/estimators/mod.rs`（`samples_of`、bias 帮助函数） | `reed-stats/src/statistic.rs` |
| `reed-io/src/lib.rs`（窗口装配、`decode_gt`） | `reed-io/src/vcf/sample_scan.rs` + `bcf/decode.rs`；装配逻辑由 `reed-stats/src/window.rs` 的分箱取代 |
| `reed-io/src/record/mod.rs` | 传输部分 → `reed-io/src/vcf/line.rs`、`field_plan.rs`；`GtMask`/`Recompute` → `reed-filter/src/spec/{mask,recompute}.rs` |
| `reed-io/src/record/vcf.rs` | `decode` → `sample_scan.rs`/`info.rs`；`apply` → `vcf/edit.rs`；`filter_header` → `header.rs`；`STALE_INFO` → reed-filter |
| `reed-io/src/record/bcf.rs` | `reed-io/src/bcf/{record,decode,edit}.rs` |
| `reed-io/src/sink.rs` | `reed-io/src/write/*` |
| `reed-filter/src/spec.rs` | `reed-filter/src/spec/*` |
| `reed-filter/src/plan.rs` | `reed-filter/src/plan/*` |
| `reed-filter/src/exec.rs` | 引擎骨架 → `reed-engine`；过滤逻辑 → `reed-filter/src/eval/*`；driver → `reed/src/filter.rs` |
| `reed-filter/src/emitter.rs` | 删除（由 L4 选择写出器） |
| `reed-filter/src/main.rs` | `reed-cli/src/cmd/filter.rs` |
| `reed-io/examples/scan.rs` | `reed/examples/window_scan.rs` |

## 附录 B：术语表

| 术语 | 定义 |
|---|---|
| chunk | 并行处理单位：连续若干 BGZF 块，压缩字节约 1 MiB，边界只由输入决定 |
| 接缝行（seam） | 跨 chunk 边界的那一条记录，由 reducer 拼接后处理 |
| 位平面（plane） | 每样本一位的位集；一个位点有 4 个（加可选的相位平面） |
| 计数表（CountTable） | 位点 × 群体的 `SiteGroupCounts` |
| 分箱（bin） | 宽度为 gcd(size, step) 的窗口网格单元；窗口由连续分箱合并而成 |
| 累加器（accumulator） | 满足单位元、结合律、切分不变的统计中间状态 |
| sketch | 过滤 pass0 保存的、按 chunk 存放的掩码后位平面 + 头部位图 |
| 票据（ticket） | 在途 chunk 名额；用于背压与内存上界 |
| 投影（projection） | 用超几何分布把调用数不同的位点映射到固定样本量 n 的 SFS |

## 附录 C：本机数据画像（2026-10-05 实测）

- `clean.vcf.gz`：6.6 GB 压缩；抽样 40 MB 解压文本含 3727 条记录 → 约 10.7 KB/记录 → 全量约 35 GB 解压文本。
- 348 样本，3,315,604 位点，20 条染色体；FORMAT 字段 `GT:AD:DP:GQ:PGT:PID:PL:PS`。
- Chr01–03 画像（`gold/out/profile.md`）：位点缺失 p0 = 0、p1 = 0.0086、p50 = 0.089；GQ<20 占 11.8% 基因型；
  HWE p<1e−6 占 91.5%（自交物种）。
- 机器：2 × Xeon Gold 5520+（56 核，2 个 NUMA 节点，105 MiB L3），503 GB 内存，rustc 1.99。
- 当前过滤结果：1,768,759 / 3,315,604 位点存活（53.3%），348/348 样本保留，32,565,162 个基因型被掩码。
