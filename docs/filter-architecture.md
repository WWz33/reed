# reed-filter 目标架构（v2，2026-10-05 重写：v0.2 落地 + 传输层诊断后修订）

> 状态：v0.2 批式流水线已实现（spec/plan/exec/emitter + CLI 薄壳），2026-10-05
> instrumentation 判定**全物化 RecordBuf 批传输是性能瓶颈**（§7.2），本版把执行模型从
> "物化记录批传输"修订为"**惰性记录批传输 + 存活者才物化**"（§7.3–7.4），落地顺序见 §12。
> 沿用 2026-10-04 三裁决（无群体维度 §8 / stdout-VCF 默认 §9 / noodles 单后端 §3）与
> §16 判决体系。事实基础：gold/ 三方金标准（全量 1,768,759 位点三方逐位一致，
> 见 `gold/out/full/concordance.md`）。
> v1 存档要点：能力矩阵（§1）、多遍数据流模型（§2）、FilterSpec（§4）、计划编译器（§5）
> 全部按原设计实现并经测试，不重述其论证过程，只记录与现状的差异。

---

## 1. 动机：外部工具能力矩阵（全部实测，2026-10-04 本机；reed 列 2026-10-05 更新）

| 能力 | bcftools 1.24 | vcftools 0.1.18 | plink2 a.7LM | reed 目标 |
|---|---|---|---|---|
| 并行 | ✗（`--threads` 只管压缩） | ✗（单线程） | ✓（24 线程） | ✓（批式流水 + MT BGZF 两端） |
| FORMAT 字段保留 | ✓（11 个全保留） | ✓（11 个） | **✗（只剩 GT）** | ✓（全保留，掩码只改 GT） |
| 按 INFO 过滤（QD/FS/MQ/SOR） | ✓ | ✗ | ✗（只能透传） | ✓ |
| 全量耗时（3.3M×348，本机） | ~40 min | ~35 min | ~2.5 min | v0.2 实测 **~9.6 min**（传输瓶颈，§7.2）→ 修订目标 **1–2 min**（§7.7） |
| 染色体名/坐标保真 | ✓ | ✓ | **✗（Chr01→1，需还原）** | ✓ |
| 掩码后重算 AC/AN/MAF/HWE | ✓（需手工剔除陈旧 INFO） | ✓ | ✓ | ✓（天生按 GT 重算） |
| ascertainment 声明 | ✗ | ✗ | ✗ | **✓（自动降级 + provenance；接线缺口见 §5）** |
| 每条规则的漏斗归属 | 部分（逐条日志） | 部分（末尾汇总） | 部分 | **✓（全量 per-rule 计数，已实现）** |

**缺口判定**（不变）：并行 + 全字段保真 + INFO/FORMAT 过滤，三者无一工具同时满足；
ascertainment 与统计语义的打通更是三方都没有——这是 reed 作为库（而非 CLI）的独有落点。

**本数据集的实测陷阱**（外部工具各自的坑，reed 设计必须绕开）：
1. 输入 INFO/AC,AN,AF,MLEAC,MLEAF **陈旧**（348 样本是大队列子集、未重算）→ 任何依赖 INFO 频率的过滤都用错值；
2. bcftools 表达式 FMT 字段上 `||` 有 bug（须用 `|`）；
3. plink2 改染色体名（Chr01→1）；
4. plink2 丢 FORMAT（DP/GQ 掩码后无法溯源）；
5. HWE 对自交物种（大豆）语义失效（p<1e-6 占 91.5%），参数必须画像驱动。

---

## 2. 核心建模转变：过滤不是谓词链，是多遍数据流问题（v1 沿用，已被实现验证）

1. **样本级缺失率（mind）是全局归约**——必须扫全基因组才能决定保留谁，然后才轮到位点过滤；
2. **MAF/HWE/位点缺失依赖"掩码后"的基因型**——DP/GQ 置缺失 → 重算 AC/AN → 才能求 MAF/HWE。谓词求值序 = 数据依赖序，不能由用户在 shell 里手工排；
3. **LD pruning 是带携带状态的顺序扫描**（plink2 `LdPruneSubcontigSplitAll` 语义，DESIGN §16.1-2），窗口化会改变数学结果——ordered-carry pass 槽位保留（§12-P5）。

reed 的抽象（不变）：

> **FilterSpec（声明式、类型化）→ 计划编译器（依赖解析 + pass 划分 + ascertainment 降级 + provenance）→ 并行执行引擎（reduce / scan / ordered-carry 三类 pass）**

---

## 3. 分层与 crate 结构（2026-10-05 按建成现状修订）

```
L0 传输   reed-io（noodles 单后端【定】，htslib 全线移除）
          · 惰性记录批传输：vcf::Record / bcf::Record（owned、Send、字段零拷贝视图）
          · FieldSet 驱动惰性解码（只解计划声明的列；series.iter 列主序一趟）
          · TransparentRecord（= noodles RecordBuf）只在"存活且需编辑"时物化
          · 线程化 BGZF：MultithreadedReader（输入，P3 接入）/ MultithreadedWriter（输出，已用）
L1 物化   记录流：未触碰存活者原始字节直通（文本输出零重序列化）；
          被编辑者物化 TransparentRecord 手术式修改（掩码只写 GT、重算只换五个 INFO tag）
          ——保真 by construction。SoA/位打包留给窗口物化层（reed-core::Bits）
L2 计划   reed-filter::spec + plan：FilterSpec（Result 构造）→ Plan（固定调度展开，
          非查询优化器；provenance 哈希 + spec 摘要）
L3 执行   reed-filter::exec：reduce（pass0 全局归约）/ scan（pass1 过滤+编辑+emit）
          两 pass 批式流水，BatchID 重排保序，ErrorSlot 首错即停
```

建成现状（与 v1 §3 布局图的差异：record 层落在 reed-io 而非 reed-filter，exec 是单文件
未拆子模块，provenance.rs 未独立成文件）：

```
crates/reed-io/src/record/mod.rs   记录传输：for_each_record（RecordBuf 逐条回调，现状）
                                   + record_header；§7.3 新增 for_each_batch（惰性批）
crates/reed-io/src/record/vcf.rs   VCF 路径：decode/count_missing/apply/filter_header
crates/reed-io/src/record/bcf.rs   BCF ↔ TransparentRecord 转换
crates/reed-io/src/sink.rs         Emitter 实现栈：Stdout / PlainFile / VcfGz / Bcf
                                   （MultithreadedWriter 压缩池）；§7.6 新增 raw-line 快路
crates/reed-filter/src/spec.rs     FilterSpec + from_conf（gold filters.conf KEY=VALUE）
crates/reed-filter/src/plan.rs     固定调度展开 + provenance 哈希 + spec 摘要
crates/reed-filter/src/exec.rs     reduce/scan 两 pass 批式流水（本次重构对象）
crates/reed-filter/src/emitter.rs  reed-io sink 的薄 re-export（库不自带 I/O 策略）
crates/reed-filter/src/main.rs     CLI 薄壳：-c / -o / --threads，零策略旗标
```

关键裁决（不变）：**过滤器不进 `ReadOptions`**。`ReadOptions` 保持 `Copy`，Reader 只认识
FieldSet（解哪些列），plan 才认识过滤器。两个变化域解耦。

---

## 4. FilterSpec：类型化声明（已实现，与 v1 设计一致）

不引入 bcftools 式表达式解析器（DESIGN §11 裁决）；builder/构造期校验/serde 按 v1 落地。
求值序（与 gold 管线 bcftools 分支逐位对齐，exec.rs 模块文档）：
**RegionRule → biallelic SNP → FILTER==PASS → QUAL → INFO（原值）→ [GT 掩码] →
位点缺失（掩码后）→ MAF（掩码后重算）→ HWE（掩码后重算）**。

**RegionRule 是位点谓词，不是执行策略**（2026-10-04 三轮裁决，不变）：spec 里只有语义
（include 且不在 exclude，bcftools `-t` 流式语义）；索引跳过/区域并行是引擎层自选，
当前实现为逐条流式判定（无索引依赖），索引类执行留 §12-P5 槽位。

`from_conf()` 直接解析 `gold/config/filters.conf`：一份 spec、四个执行器对拍，参数漂移
在构造上不可能发生。

---

## 5. 计划编译器（已实现子集 + 接线缺口诚实记）

`plan::compile` 现做四件事中的两件半：**pass 激活判定**（cohort 规则 ⇒ reduce pass）、
**provenance 哈希**（canonical spec 渲染 → sha256）、**spec 摘要**（人读一行）。
谓词代价排序以"头部规则按 spec 声明序、尾部规则固定 missing→MAF→HWE"落地（HWE 精确
检验天然在最后，满足短路初衷）。

**接线缺口（2026-10-05 盘点）**：
1. `exec.rs` 里 `let _plan = plan::compile(spec);`——provenance 哈希与摘要**未写入输出
   header**（`##reed_filter=` / `##reed_filter_spec=` 未接线），随 §12-P2 一并接上；
2. FieldSet 求解在 v0.2 是名义的（RecordBuf 全物化无视字段需求）——§7.3 惰性传输落地后
   才成为真实优化（pass0 只解 GT、头部规则只解位点字段）；
3. ascertainment 自动降级推导仍是 spec 摘要级，结构体化迁移见 §6。

编译产物（gold 默认 spec 的实例，不变）：

```
Plan {
  passes: [ Reduce { rule: sample_missing<=0.20, fields: FMT/GT-validity },
            Scan   { stages: [mask(DP∈[5,100],GQ≥20)] → [recompute AC,AN,NS,MAF,HWE]
                            → [biallelic, PASS, QUAL≥30, QD≥4, FS≤50, MQ≥40, SOR≤3,
                               missing≤0.10, MAF≥0.05, HWE≥1e-6] } ],
  output: { emitter: StdoutVcf, fidelity: full-INFO-FORMAT, overwrite: [AC,AN,NS,MAF,HWE] },
  plan_hash: "sha256:…"
}
```

**输出保真策略**（不变）：未被重算覆盖的 INFO/FORMAT 原样透传；掩码只改 GT 为 `./.`，
DP/GQ 原值保留（可溯源，胜过 plink2）。§7.4 把"透传"强化到字节级。

---

## 6. Ascertainment v2：manifest 结构体化（方向不变，迁移待做）

`SiteCompleteness` 升维结构体的设计沿用 v1：

```rust
pub struct Ascertainment {
    pub sites: SiteSet,            // Complete | Filtered{plan_hash} | Unknown
    pub invariants: InvariantSet,  // AllPresent | PartiallyRemoved{plan_hash} | Absent
    pub genotypes: GenotypeState,  // Original | Masked{rules}
    pub cohort: CohortState,       // Full | Subset{keep_hash, n}
}
```

降级规则由计划编译器自动推导；估计器 availability 门从 match 枚举改为 match 结构体
（机械迁移）。**当前状态：未迁移**（plan.rs 只产出摘要字符串）。排在传输重构之后
（§12-P5），不与性能工作抢同一次 breaking change。

---

## 7. 执行引擎（2026-10-05 全节重写：v0.2 实测诊断 → 传输架构修订）

### 7.1 v0.2 建成形态（批式 RecordBuf 流水线）

两 pass 串行（样本归约必须先完成）。每 pass 三段流水：reader 线程拥有 noodles reader
（`for_each_record` 逐条产出 **RecordBuf**）→ 4096 条/批进 crossbeam 有界分发通道
（容量 4）→ `opts.threads` 个 compute worker 跑记录级管线（pass0：GT 缺失计数；
pass1：decode → 头部规则 → 掩码 → 计数 → 尾部规则 → 存活者 apply）→ emit 线程按
BatchID 重排（BTreeMap + next_batch_id），严格按输入序写 sink，漏斗部分和按批序整数
求和。确定性：跨线程数输出字节一致（有测试）。内存：通道容量界定在途批数。

### 7.2 实测诊断（2026-10-05 instrumentation，clean.vcf.gz 6.6GB 压缩 / 348 样本 × 3.3M 位点 / 56 核 2×28 Xeon Gold 5520+ 双 NUMA / 503GB RAM）

**观测**（exec.rs 临时探针输出）：

| 阶段 | 观测值 | 解读 |
|---|---|---|
| pass0 reader | wall=**283s**，callback(push+send)=**0.04s** | 283s 全在 `for_each_record` 内部（inflate + RecordBuf 物化）；发送从未阻塞 |
| pass0 workers | 近乎全闲 | 阻塞 recv，park 不占 CPU |
| pass1 reader | wall=**292s** | 同上 |
| pass1 compute | 224s/worker，全隐藏 | worker 已 ~77% 占用——reader 修好即轮到 compute 当墙 |
| pass1 emit | 全隐藏 | 写出/压缩不是瓶颈 |
| **总时长** | ≈ 2 × 串行物化 parse（~9.6 min） | 对照：v0.1 惰性 GT-only 单线程**全程** 187.3s（DESIGN §2 实测，同文件含 inflate 与 GT 解码） |

**根因排序**（机制归属；证据见各项）：

1. **传输物是全物化 RecordBuf（主因）**。`record_bufs()` = 行读取 + `parse_record_buf`
   全量物化（noodles-vcf `io/reader.rs:215-231` → `variant/record_buf/convert.rs`
   `try_clone_from_variant_record`）：每条记录把全部 INFO tag 物化为 owned IndexMap、
   8 个 FORMAT 字段 × 348 样本物化为 owned `Value`——**每条记录上千次小堆分配，且两
   pass 各付一遍**。而 pass0 只需 GT 有效性、pass1 头部规则只需位点级字段：全量物化是
   纯浪费。v0.1 走 `reader.records()` 惰性记录（`variant/record/fields.rs:12-15`：
   `Fields{ buf: String, bounds }`，一条记录一次分配、字段零拷贝视图）+ GT 列
   `series.iter` 一趟——同文件全程 187.3s 还包含了解码计算。
2. **跨线程分配/释放 + 系统分配器（流水线特有放大器）**。所有 RecordBuf 在 reader 线程
   分配、在 worker/emit 线程释放（pass0 整批 drop；pass1 非存活者 worker drop、存活者
   emit drop）。Cargo.toml 无自定义分配器 → glibc malloc：远程 free 必须锁拥有者
   （reader 的）arena，且 free 不回流 reader 的 tcache——reader 失去 tcache 快路径，
   每次 malloc 与数十个 worker 的远程 free 竞争同一把 arena 锁；双 NUMA 节点进一步
   放大。malloc 密集型 parse 在此模式下慢数倍是典型表现，且时间全部呈现为"parse
   inside for_each_record"——与观测精确吻合。
3. **在途足迹（第二放大器）**。4096 条/批 ×（4 分发 + 8 产出 + worker 在手）批的
   **全物化**记录，宽队列下在途数十 GB；reader parse 时持续首次触页（minor fault 记在
   reader 头上），而 v0.1 足迹≈一条行缓冲、全程热页复用。（503GB RAM 无 swap；页缓存
   122GB 使 pass1 重读输入近似免费。）
4. **单线程 BGZF inflate（独立地板）**。6.6GB 压缩输入单线程 inflate 地板 33–66s/遍
   （100–200MB/s）。v0.2 代码注释已自标 "the next performance lever"。

**排除项**（instrumentation 反证）：通道发送从未阻塞（callback 0.04s ⇒ 20ms
`send_timeout` 轮询未触发）；worker park 不占 CPU；`Instant::now()×2`/记录与
ErrorSlot 原子读为纳秒级。**调度/通道设计无罪——问题纯粹是"通道里运的是什么"。**

**探针语义**（`probe_reader_variants`，exec.rs 测试，`REED_PROBE_INPUT=<vcf.gz>
cargo test --release -- --ignored`）：V1 主线程顺序 = 基准；V2 spawn 线程顺序 = 隔离
"线程本身"（预期 ≈V1）；V3 批式流水线 = V3−V2 即批传输+通道+远程 free+足迹的净效应；
V1 的绝对水平对照 187.3s 即"物化传输 vs 惰性传输"的净效应。判定树：V1≈283s ⇒ 传输
物化是唯一故事；V1≈187s 而 V3≈283s ⇒ 远程 free/足迹叠加成立，两者都要修。

**为什么不是"换分配器 / 回收 RecordBuf"了事**：RecordBuf 在 noodles 0.117 无
clear/reuse API（源码核查 `variant/record_buf.rs`），回收不可行；mimalloc/jemalloc 能
吸收远程 free 惩罚，但每遍 3.3M × 全字段物化本身是约两个数量级的多余工作——分配器是
保险不是修复（保留为 §12-P4 实测项）。

### 7.3 修订：惰性记录批传输（lazy-batch transport）

**设计使能事实**（2026-10-05 vendored noodles 0.117 源码核查）：
- 惰性 `vcf::Record` = `Fields{ buf: String, bounds: Bounds }`——**拥有自身行缓冲**，
  自动 `Send`；所有字段访问是 `buf` 上的零拷贝视图；每条记录一次分配（可池化）。
- `bcf::Record` 同样自有二进制缓冲——同一批传输形状覆盖两种输入。
- `bgzf::io::MultithreadedReader` 存在且实现 `Seek`
  （`multithreaded_reader.rs:323` `seek_to_virtual_position`）。

```
reed-io::record 新增：
  for_each_batch(path, opts, |batch_id, RecordBatch| -> Result<...>)
  enum RecordBatch { Vcf(Vec<vcf::Record>), Bcf(Vec<bcf::Record>) }   // 惰性、owned、Send
  SiteView<'r>   // 借用惰性记录：contig/pos/qual/is_pass/INFO{QD,FS,MQ,SOR}，零拷贝
  gt_series / u32_series（GT/DP/GQ 列，series.iter 列主序一趟——v0.1 已验证姿势）
  count_missing_lazy（pass0 用：全 pass 零物化）
```

- **批次按未压缩字节数切**（目标 16–32MB/批），不按记录条数：宽队列下 4096 条物化批
  可达数 GB，同条数惰性行批只有原始文本大小；在途内存 ≈ 通道容量 × 批字节，可预测。
- **reader 线程成本回到"读"本身**：inflate + 行填充。解码、规则、物化全部下沉到
  worker（并行），reader 不再碰任何记录字段。
- **跨线程流量降三个数量级**：每条记录一次 String 分配（vs RecordBuf 上千次），且
  惰性记录由解析它的 worker 自行释放——远程 free 面只剩"存活者物化记录"（见 §7.4）。
- FieldSet 求解由此成为真实优化（§5-2 的落地）：pass0 只解 GT 列，头部规则只解位点
  字段，每少解一列少付一列的钱。

### 7.4 编辑路径：存活者分流 + 原始字节直通

存活者按是否需编辑二分：

| 存活者类型 | 条件 | 路径 |
|---|---|---|
| **Verbatim** | 无需剪列 && 未触发掩码 && recompute 全关 | 原始行字节直通文本 sink（**零重序列化**，保真强化到字节级） |
| **Edited** | 剪列 \|\| 掩码触发 \|\| recompute 激活 | worker 上物化 TransparentRecord → 现有 `apply` 手术（不变）→ sink |

- 默认 spec（recompute 全开）下所有存活者走 Edited——物化量仍从"每遍 × 每条输入记录"
  降到"仅存活者"（漏斗末端），且物化发生在 **worker（并行）而非 reader（串行）**。
- outcome 携带枚举 `{Verbatim(行字节), Edited(TransparentRecord)}`；emit 阶段重排逻辑
  不变，按变体调 sink 对应入口。
- BCF 输出无直通（二进制必须重序列化）：一切存活者物化。BCF 输入 → 文本输出时
  Verbatim 不适用（无原始文本行），存活者物化（二进制解码本就便宜，§16.2）。

### 7.5 索引策略与区域并行（沿用 2026-10-04 定稿，槽位不变）

索引是按需加速器，不是输入前提；RegionRule（targets 语义）永远不需要索引。区域并行
（引擎自选，需 .csi/.tbi）留 §12-P5：pass0/pass1 都按区域劈开（归约各区域部分和再合并），
区域并行时各区域用普通 `Reader`（MT reader 各带线程池代价重）。stdin/管道永远是
流式流水线形态。

### 7.6 输出汇编（沿用 + raw-line 快路）

单一共享 `MultithreadedWriter` + 按 BatchID 重排序提交（不变；实测 348 样本负载
5.2–5.5 GB/s，输出与单线程字节级同尺寸，经 bcftools view + `bgzip -t` 双验证）。
文本 sink 新增 **raw-line 快路**：header 仍经 `vcf::io::Writer` 写一次（VCF≥4.4 → 4.3
版本钳制逻辑保留），Verbatim 存活者的行字节直写下层 writer。`bcf::io::Reader::new(r)`
会二次包 BGZF 的坑（必须 `Reader::from(r)`）在输入侧 MT 化时同样适用，已记录。

### 7.7 性能模型（修订后诚实预算，全量 3.3M×348）

- **现状（v0.2 RecordBuf 批传输）**：pass0 283s + pass1 292s ≈ 9.6 min，瓶颈 = 两遍
  串行物化 parse（§7.2）。
- **地板**：单线程 inflate 33–66s/遍；v0.1 惰性 GT-only 单线程全程 187.3s（含 inflate
  与 GT 解码）是"单线程惰性扫描"的实证锚点。
- **修订后**：P1+P2 后 reader = inflate+行填充（串行 ~40–60s/遍），解码摊到 28–56 个
  worker（v0.1 的 187s 单线程 decode 工作量 ÷ 并行度）；P3 MT inflate（8–16 线程）
  后 reader 地板 ~5–10s。**总时长目标 ≈ 1–2 min**（两遍串行 + emit/压缩重叠），对齐
  v1 的流水线预估 1.5–2.5 min，对比 bcftools ~40 min 快一个数量级，与 plink2 2.5 min
  同量级但保真（plink2 丢 FORMAT、输出仅 0.26 GB）。
- **pass1 compute 预警**：v0.2 实测 worker 已 77% 占用；reader 修好后 compute 即墙——
  §7.3 同时削减 compute（非存活者零物化零 RecordBuf drop），§12-P4 留 scratch 池化
  与批次调参的实测槽。

---

## 8. 群体维度的裁决：过滤器层不引入（2026-10-04 定稿，不变）

PopScope / PerPop / MinOverPops 全部移出过滤器层：群体结构是统计层概念（reed-core 的
`PopulationView`）；逐群体过滤 = keep-list 组合（`CohortRule::KeepList` 已覆盖）；过滤后
分组零成本（统计层想怎么分就怎么分）。若未来成为真实需求，正确位置是统计层的
pre-processing 组合器，不是 FilterSpec 的字段。

## 9. Emitter：stdout-VCF 默认，CLI 级再加落盘格式（2026-10-04 裁决，已实现）

无 `-o` → stdout 文本 VCF（Unix 过滤器语义）；`-o out.vcf.gz` / `-o out.bcf` 按扩展名
走 MultithreadedWriter 池；`--threads` 是唯一资源旋钮（零策略旗标）。实现形态：reed-io
`sink::open` 返回 `Box<dyn VariantSink>`（`write_header`/`write_record`/`finish`），
reed-filter `emitter` 仅薄 re-export。§7.4 的 raw-line 快路为文本 sink 增加第三入口
（trait 加方法或内部特判，P2 定）。`WindowFeeder`（过滤直喂窗口 Assembler，零中间落盘）
保留为 §12-P5 槽位——§11"消除文件侧效应"在此闭环。

## 10. 验证架构（与 gold/ 合流，oracle 独立性原则；增加性能门）

**gold/ 永远由外部三工具生成，不引入 reed 代码**。五级验收不变：

| 级 | 内容 | 判据 |
|---|---|---|
| L_site | 位点集合 | 与 gold.vcf.gz 完全一致（(CHROM,POS,REF,ALT)） |
| L_gt | 基因型 | 同名样本逐位一致（掩码后的 GT） |
| L_tags | 重算标签 | AC/AN/NS 整数精确一致；MAF/HWE 相对容差 1e-9 |
| L_manifest | ascertainment | manifest 正确描述输出（自洽性） |
| L_props | 性质测试 | 幂等、阈值单调、cohort 子集不变性、跨线程数字节一致 |

**新增性能门（2026-10-05）**：性能回归与对拍同级——每次传输层改动记录 pass0/pass1
reader wall、worker work、总时长进 bench 日志；验收线：修订后总时长 ≤ 2 min
（clean.vcf.gz，本机）。check_parity.py 扩 `--tags` 级；CI 用合成 fixture，全量对拍按需。

## 11. 范围纪律（§16 判决沿用表，不变）

| §16 判决 | 本架构 |
|---|---|
| 不做 FilterPipeline<WindowBuffer>（整窗谓词） | 沿用——过滤在记录/批次层，不进窗口 |
| 谓词型/变换型 stage 分开 | 沿用——Scan 内 stages 有序组合，类型不同 |
| LD pruning 不做 per-window | 沿用——独立 ordered-carry pass 槽位（plink2 同构） |
| 不做 EstimatorRegistry / MultiScaleView | 沿用——估计器仍是自由函数 |
| 不做字符串表达式语言 | 沿用——类型化 builder + serde |
| 样本级 keep/remove 是真需求 | 吸收——CohortRule + 全局归约 pass |
| 过滤器必须携带清单降级规则 | 升级——Ascertainment v2 结构体（§6，待迁移） |
| GQ/INFO 过滤立项前先解后端 | 定稿推翻——noodles 单后端；§7.3 惰性解码使"多解一列多付一列"成为 FieldSet 语义而非风险 |

不做的（本版明确排除）：GL/PL 似然过滤、VQSR、SnpGap/IndelGap（数据无 indel；架构留
windowed-neighborhood 槽位）、CRAM、群体维度（§8）。

---

## 12. 重构落地顺序（2026-10-05 修订：P0–P5，每步带验收门）

> 原则：先确认机制（P0），再动传输（P1–P2），再加并行解压（P3），打磨只看 profile（P4）。
> 每步结束跑确定性测试 + gold 全量对拍，PASS 才算完；性能数字回填 §7.2/§7.7。

**P0 基线与机制确认（不改代码，半天）**
跑 `probe_reader_variants`（`REED_PROBE_INPUT=clean.vcf.gz`，`--release`）取 V1/V2/V3；
需要时补 heaptrack/perf 取 malloc 画像。判定树见 §7.2：V1≈283s ⇒ 传输物化是唯一故事
（调度无罪，P1–P2 即全部修复）；V1≈187s 且 V3≈283s ⇒ 远程 free/足迹叠加成立
（P1–P2 仍为主修，P4 分配器升优先级）。
**门**：机制归属与基线数字写入本文档 §7.2。

**P1 reed-io 惰性批次传输（VCF 文本先行，BCF 跟随）**
`for_each_batch` + `RecordBatch` + `SiteView` 惰性解码 + GT/DP/GQ series 解码 +
`count_missing_lazy`；批次按未压缩字节定界（target_bytes 参数，默认 16–32MB）。
RecordBuf 路径保留（apply/BCF 输出/header 操作仍需要）。
**门**：reed-io 单测——同一 fixture 上惰性解码 ≡ RecordBuf 解码（覆盖：GT 三形态
0/1、1|1、./.；半 call 0/.、./1；`.` 样本列；GT 字符串形；无 GT 列；INFO 缺失 tag）；
clippy 0 警告。

**P2 reed-filter exec 切换到惰性批（本计划的核心收益）**
- reduce pass：worker = GT-only 惰性 count_missing（全 pass 零物化）；
- scan pass：worker = SiteView 头部规则 → series 解码 → 掩码/计数/尾部规则 → 存活者
  分流 {Verbatim, Edited}（§7.4）；Edited 在 worker 上物化 + apply；
- sink：文本路径 raw-line 快路（§7.6）；
- 批次字节定界替换 4096 条定界；通道容量按字节预算重标定；
- 顺手接线：provenance 哈希/摘要/漏斗写入输出 header（§5 缺口）。
- **保留不动**：BatchID 重排 + ErrorSlot 首错即停 + 漏斗整数部分和（确定性三件套）。
**门**：跨线程数输出字节一致（现有确定性测试）；gold 全量对拍 L_site/L_gt/L_tags；
性能门：总时长 ≤ 3 min（未上 MT inflate 的中间态）。

**P3 MT BGZF 输入**
`.gz/.bcf` 走 `MultithreadedReader`（worker 数随 `--threads` 或独立旋钮）；
`bcf::io::Reader::from(r)` 姿势（`new` 二次解压坑）。
**门**：对拍不变；reader wall ≈ inflate 地板 × 安全系数；性能门：总时长 ≤ 2 min。

**P4 实测驱动的打磨（仅当 profile 说话）**
行缓冲池/批次 slab 回收（reader↔worker 返回通道）；mimalloc（若远程 free 仍可见）；
worker decode scratch 池化（gt/dp/gq Vec 复用）；target_bytes 调参；NUMA interleave。
**门**：每项独立前后对照，无提升即撤回。

**P5 路线图槽位（不变，不在本轮）**
BCF 输入全惰性链对拍；区域并行（§7.5）；ordered_carry（LD pruning，plink2 同构）；
Ascertainment v2 结构体迁移（§6，breaking，与 reed-core 0.3 合并）；WindowFeeder（§9）。

**风险登记**：
1. 惰性解码语义对等——半 call / `.` 样本列 / GT 两形态 / 无 GT 列：复用现有
   `gt_alleles` 函数族，P1 测试矩阵全覆盖；
2. INFO 标量加宽规则（Float f32→f64、Integer→f64、tag 缺失=None）与现有 `info_f64`
   逐行对齐；
3. raw-line 快路的 header 一致性：版本钳制只发生在 header 写一次处，行字节不触碰——
   与现实现逐字节对拍（P2 门）；
4. 字节定界批次的背压：批数上界仍由通道容量控，内存上界 ≈ 容量 × target_bytes ×
   (1 + 存活者物化系数)，P2 实测标定；
5. pass1 compute 成为新墙（§7.7 预警）：P2 后第一件事是重测 worker work，不为过。

---

## 13. v0.2 实施记录（2026-10-05 盘点，供后续轮次对照）

已建成并经测试：FilterSpec（builder/from_conf）、plan（固定调度 + provenance 哈希）、
exec（reduce/scan 两 pass 批式流水、BatchID 重排、ErrorSlot、RunReport 漏斗）、emitter
（stdout 默认 + VcfGz/Bcf MT 压缩池）、CLI 薄壳；跨线程数确定性测试（10,000 记录合成
fixture，threads=0/1/3/7 字节一致）；reed-io 记录层（for_each_record 三格式分派、
decode/count_missing/apply/filter_header、STALE_INFO）与 sink 三格式回环测试。
已知缺口：provenance 未写 header（§5）；Ascertainment v2 未迁移（§6）；RegionRule 仅
流式谓词（索引执行未做，§7.5）；BCF 输入路径未接入 reed-filter 端到端对拍。
实测性能瓶颈与修订方案：§7.2 → §7.3/§7.4，落地顺序 §12。
