# 02 BCF 编解码层（收尾 R6：BCF 输入与 `-Ob/-Ou` 输出）

> 目标：BCF 与 BGZF VCF 享有完全相同的保证——多核并行、字段保真、对线程数确定、与 bcftools 结果一致。
> 对齐原则不变：BCF 的一切行为以 **hts-specs BCFv2.2 与 htslib `vcf.c`** 为准，reed 不发明新的编码规则。
> 本文给出格式要点、解码框架（含"长度前缀链"并行定界算法）、输出编码、四种输入输出组合的实现路径、测试门与分期。

---

## 0. 现状与约束【据】

- 现在的 `reed::input::open` 能探测出 BCF，但直接报 "not supported yet"（`reed/src/input.rs`）；`-O b/u` 被 CLI 拒绝。
- 引擎是"文本行"导向的：chunk 内的行归属靠 `\n` 定界（`reed-io/src/chunk/mod.rs::inflate_chunk`），
  所有消费者都接收 `LineRef`/`SiteView`（文本）。BCF 没有换行符，不能直接套用。
- 本机数据 `clean.bcf`（8.06 GB，带 `.csi`）第一条记录的实测布局：`l_shared=134`、`l_indiv=11856`、
  8 个 FORMAT 字段：GT int8×2、AD int16×2、DP int16×1、GQ int8×1、PGT char×4、PID char×14、PL int16×3、PS int8×1，
  每样本 34 字节，与 VCF 文本的每行体积（约 10.7 KB）相当。
  **但 stats/filter 只需要 GT+DP+GQ（每样本 5 字节，约 1.7 KB/记录），其余字段按描述符长度直接跳过**——这是 BCF 输入相对文本的根本优势。

---

## 1. 格式要点（BCFv2.2，htslib 实现）

| 项 | 规则 | 出处 |
|---|---|---|
| 文件 | BGZF 压缩；解压后以 `BCF\x02\x02` 开头，接 `l_text:u32` 与以 NUL 结尾的 VCF header 文本，之后是记录 | hts-specs BCFv2.2 §6.1 |
| 记录 | `l_shared:u32`、`l_indiv:u32`，随后 shared 块（l_shared 字节）与 indiv 块（l_indiv 字节） | §6.3 |
| shared 块 | `CHROM:i32`（contig 字典下标）、`POS:i32`（0-based）、`rlen:i32`、`QUAL:f32`、`n_allele<<16 \| n_info:u32`、`n_fmt<<24 \| n_sample:u32`、ID（typed string）、n_allele 个等位（typed string）、FILTER（typed int 向量，字典下标）、n_info 个 INFO（typed int 键 + typed 值） | §6.3.1 |
| indiv 块 | n_fmt 个字段，**字段主序**：typed int 键、typed 描述符（类型 + 每样本元素数），随后 n_sample × 元素数个值 | §6.3.2 |
| typed 描述符 | 1 字节：低 4 位类型（0 MISSING，1 int8，2 int16，3 int32，5 float32，7 char），高 4 位个数；个数为 15 时后接一个 typed int 给出真实个数 | §6.3.3 |
| 整数哨兵 | int8：缺失 `0x80`（−128）、向量结束 `0x81`（−127），−128…−121 保留，**合法范围 −120…127**；int16/int32 同理（`BCF_MIN_BT_INT16 = −32760`、`BCF_MIN_BT_INT32 = −2147483640`） | `htslib/vcf.h:1517-1543` |
| 浮点哨兵 | 缺失 `0x7F800001`，向量结束 `0x7F800002`（按位比较，不是数值比较） | `vcf.h` `bcf_float_missing/vector_end` |
| 字符串 | char 向量，不足部分以 NUL 填充；FORMAT 字符串按样本等宽 | §6.3.3 |
| GT | `(allele+1)<<1 \| phased`；缺失等位 = 0；倍性低于向量宽度的样本用向量结束标记补齐 | `vcf.h:1033-1035` |
| Flag | INFO 值编码为类型 NULL、个数 0 | `vcf.c:3954` |
| 字典 | 字符串字典（FILTER/INFO/FORMAT 共用）与 contig 字典分开；`PASS` 固定为 0；header 行可带 `IDX=` 显式指定下标，否则按出现顺序 | §6.2.1 |
| 整数宽度选择 | 编码时取能容纳该向量所有值的最小类型（`bcf_enc_vint`，`vcf.c:2909-2960`） | htslib |

---

## 2. 输入：解码

### 2.1 探测与 header

- `Format::detect` 已能根据首块解压内容识别 `BCF\x02`。接受 v2.1 与 v2.2；更低版本报错。
- header：读 `l_text`，用 noodles 解析 header 文本（得到 `vcf::Header` 与 `StringMaps`，后者已支持 `IDX=`，
  `noodles-vcf/src/header/string_maps.rs`），并构建 reed 的查找表：

```rust
pub struct BcfDictionaries {
    pub contigs: Vec<ContigId>,            // BCF contig 下标 → ContigId（header 序）
    pub strings: Vec<Option<StringEntry>>, // 字符串字典下标 → { kind: Info|Format|Filter, tag: TagId }
    pub gt_key: Option<u32>, // 预先解析 GT/DP/GQ 等计划字段的字典下标，热路径只比较整数
}
```

- 记录区起点 = `5 + 4 + l_text`（解压后偏移），作为第一个 chunk 的 `first_record_offset`。

### 2.2 并行定界：长度前缀链【规】

BCF 记录没有同步标记，chunk（BGZF 块边界）内部无法独立找到第一条记录的起点。
但定界只需要"起点偏移"这一个整数沿 chunk 序传递，而沿一个 chunk 走完长度前缀只需 O(记录数) 次跳转（约 500 次、微秒级）。由此：

1. reader 不变：仍然只读 BGZF 块头、按压缩字节切 chunk（BAPP）。
2. worker 拿到 chunk k：**先 inflate**（昂贵、并行，不依赖任何人）。
3. 等待 `start_k`（chunk k 内第一条完整记录的偏移）。chunk 0 的 `start_0` 来自 header 解析。
4. 从 `start_k` 开始按 `8 + l_shared + l_indiv` 跳转，直到越过 chunk 末尾，得到 `start_{k+1} = 越界记录的结束偏移 − len(chunk_k)`，**立即发布**（先发布、后解码）。
5. 头片段 = `[0, start_k)`，体 = 完整记录，尾片段 = 越界记录在本 chunk 内的部分。之后与 VCF 完全一样：接缝记录由 reducer 拼接
   `tail(k−1) ++ head(k)`（现有 `stitch_seam`），超长记录跨多个 chunk 时 `start` 可以大于整个 chunk 长度，链照样成立（该 chunk 没有记录起点，整体是头片段）。

**无死锁论证**：chunk 通过 FIFO 通道按 id 分发，任意时刻 id 最小的未完成 chunk，其前驱都已完成并发布了起点，所以它不会无限等待；
其他 chunk 至多等待前驱的"走链"时间（微秒级），而前驱的 inflate 与自己的 inflate 是并行发生的。

**确定性**：chunk 边界仍只由输入决定；起点由数据唯一决定；与线程数无关。

**实现位置**：

```rust
// reed-io/src/chunk/framing.rs
pub enum RecordFraming { Newline, LengthPrefixed }
// reed-engine/src/frame_chain.rs：容量 = 票据数的环形数组，每格 (chunk_id, start) + Condvar；
// worker 调用 wait_start(k) / publish_start(k+1, offset)。中止信号到来时所有等待立即返回。
```

`ChunkJob` 不变，引擎在调用 `run_chunk` 前完成定界。`RecordFraming` 由 `ChunkSource` 声明（BGZF VCF / plain VCF 为 `Newline`，BCF 为 `LengthPrefixed`）。

**有索引时的快路径（可选）**：`-r/-R` 区域查询时，CSI 给出的虚拟偏移本身就是记录起点（低 16 位 = 块内偏移），
区域内第一个 chunk 的起点直接已知；后续 chunk 仍走链。现有 `RegionChunkSource` 只需把 `first_in_block_offset` 解释为记录起点。

### 2.3 记录视图与格式无关的消费接口【规】

现在所有消费者都拿文本 `LineRef`/`SiteView`。BCF 接入前先把消费接口抽象成格式无关的两部分：

```rust
/// 位点级字段（头部规则、表达式、报告用）。VCF 与 BCF 各实现一份，零分配、惰性。
pub trait SiteFields {
    fn contig(&self) -> ContigId;
    fn pos1(&self) -> u64;
    fn ref_allele(&self) -> &[u8];
    fn n_alt(&self) -> usize;
    fn alt_allele(&self, i: usize) -> &[u8];
    fn variant_class(&self) -> VariantClass;      // 两种格式用同一个分类函数
    fn qual(&self) -> Option<f32>;
    fn filter_ids(&self) -> FilterIds<'_>;        // 文本：按名字；BCF：按字典下标（统一成 TagId）
    fn info(&self, tag: TagId) -> Option<InfoValue<'_>>;
    fn id(&self) -> &[u8];
}

/// 基因型与计划内 FORMAT 列：输出与现在的 ScanOutput 完全相同（codes + int/float 列 + 可选等位下标）。
pub trait GenotypeSource { fn scan(&mut self, counters: &mut PolicyCounters) -> Result<ScanOutput<'_>, IoError>; }

pub enum RecordRef<'a> { Vcf(LineRef<'a>), Bcf(BcfRecordRef<'a>) }
```

- 表达式求值器的 `SiteValues`、过滤的头部/尾部规则、stats 的计数、reports、summary 都改为面向这两个 trait，不再直接接触 `LineRef`。
- **差异只允许出现在 trait 实现内部**；上层逻辑一份。这是 BCF 工作的前置步骤（05-roadmap.md V1）。

### 2.4 GT 解码

- int8、2 元素（最常见）：用 `(b0, b1)` 两字节查一张 65,536 项的 `GtCode` 表（64 KiB，进 L2），**每样本一次查表**，没有分支。
- 表的内容由 `GtCode::from_alleles` 生成，策略（ploidy / multiallelic / half-call）与 VCF 路径**同一个函数**，保证两种格式语义一致。
- 向量结束标记（倍性 < 2）按单倍体处理，走同一策略；int16/int32 编码的 GT（等位下标 > 62）走慢路径逐元素解码。
- 需要原始等位下标时（多等位的逐等位计数），从同一字节直接得到 `(v >> 1) − 1`。

### 2.5 FORMAT 数值列与跳过

- 计划内的 DP/GQ 等：按描述符类型逐元素拓宽到 `u32`/`f32`，缺失与向量结束映射到 `INT_MISSING`/NaN（与文本扫描相同语义）。
- 计划外字段：`跳过字节 = 类型宽度 × 元素数 × n_sample`，一次加法。
- 表达式中的 `FMT/` 字段同样按计划解码。

### 2.6 错误

字典下标越界、记录长度越过数据末尾、描述符类型非法、哨兵出现在不允许的位置 → `IoError::MalformedRecord { locus }`；
contig 下标不在字典中 → 与 VCF 一致的"未声明 contig"错误。

---

## 3. 输出：编码

### 3.1 输入 × 输出组合

| 输入 → 输出 | 未编辑的存活记录 | 被编辑的存活记录 | 分期 |
|---|---|---|---|
| VCF → VCF | 原样拷贝 | `LineEditor` 字节手术（现有） | 已有 |
| VCF → BCF | 文本 → BCF 编码 | `LineEditor` 后 → BCF 编码 | **B2** |
| BCF → VCF | BCF → 文本渲染 | 渲染后交给 `LineEditor`（复用） | **B3a** |
| BCF → BCF | **记录字节原样拷贝**（需满足 §3.3 的字典条件） | 阶段一：渲染 → `LineEditor` → 编码；阶段二：原生 BCF 手术 | B3a / **B3b** |

原则：文本 VCF 是编辑的规范中间形态（`LineEditor` 已经与 bcftools 逐字节对齐），BCF 首先通过"渲染/编码"接入，
性能需要时再做原生手术。这样 B2/B3a 完成后，四种组合都已正确，B3b 只是加速。

### 3.2 写出容器

- `-Ob`：BGZF，默认级别 6；`-Ou`：**BGZF 0 级块**（htslib 的 "uncompressed BCF" 就是 0 级 BGZF，不是裸流）。
  两者都走现有 `BlockWriter`（切块规则不变，记录可以跨块）。
- 文件头：`BCF\x02\x02` + `l_text` + header 文本（NUL 结尾）。
- 索引：BCF 只支持 CSI（`-W=tbi` 与 `-Ob` 同时出现时报错，与 bcftools 一致）；BCF 的 CSI **没有** tabix 辅助数据
  （VCF 的 CSI 有），`Indexer` 需要区分这两种 flavor；坐标用 0-based `POS` 与 `rlen` 计算区间。

### 3.3 字典与 `IDX`（BCF→BCF 正确性的关键）【规】

- BCF 记录里的 INFO/FORMAT/FILTER 都是**字典下标**。只要输出 header 的字典与输入不同，原样拷贝的记录就会被解读成错误的标签。
- 规则：
  1. **BCF 输出的 header 一律写显式 `IDX=`**（htslib 写 BCF header 时就是这样做的，`bcf_hdr_fmt_text(..., is_bcf=1)`）；
     noodles 能读 `IDX`，但写 header 时不会自动补上，所以由 reed 渲染 header 文本时补齐。
  2. **下标永不重排**：`-x` 删除的标签，其字典项保留在 header 中（或保留 `IDX` 空位），绝不能让后续标签前移；
     `--fill-tags` 新增的标签取 `max(IDX)+1` 起的新下标。
  3. 剔除样本（`-s/-S/--keep/--remove/--mind`）时所有记录的 indiv 块都要重建（没有原样拷贝路径）；contig 字典不变。
- 测试：构造"删除中间一个 INFO 标签 + 新增两个标签"的 BCF→BCF 用例，`bcftools view` 读回后与 VCF 路径逐字节一致。

### 3.4 B2：文本 → BCF 编码（先用 noodles，正确性优先）

- 每个 worker 持有一个 `noodles::bcf::io::Writer::from(Vec<u8>)`：启动时调用一次 `write_header(&out_header)`
  以建立 `StringMaps`，然后清空内部 `Vec`；之后每条存活记录：
  `LineEditor` 产出文本 → 构造惰性 `vcf::Record`（`vcf::io::Reader::new(bytes).read_record`）→
  `write_variant_record(&out_header, &record)` → 取出字节交给 `BlockWriter`。
  noodles 的编码器（`noodles-bcf/src/io/writer/record.rs`）对任何 `variant::Record` 都可用。
- 整数宽度选择、字符串填充、GT 编码由 noodles 完成；与 htslib 的差异（如果有）只能是"bcftools 读回后语义相同"范围内的，
  由 §4 的门验证。
- 预估成本：noodles 逐字段编码约 50–100 µs/记录（348 样本 × 8 字段）；本数据 1.77M 存活记录 ≈ 90–180 CPU-s，
  54 个 worker 分摊约 2–3 s。可以接受；超出时进入 B4 原生编码器。
- **header 中未声明的标签**：BCF 必须用字典下标，无法编码未声明的标签。htslib 在 VCF→BCF 时会自动补一个定义并警告；
  reed 是并行写出（header 已先写出），无法事后补充，因此**报错**，并提示用 `bcftools annotate -h` 补 header 行。

### 3.5 B3a：BCF → 文本渲染（与 `bcftools view` 逐字节一致）

- 原生实现（不用 noodles 的文本渲染器，它的浮点格式与 htslib 不同）：
  - 整数 → 十进制；浮点 → htslib 的 `%g` 规则（`line_edit.rs` 中现有的 `format_g6` 正是这一规则，复用）；
  - 缺失 → `.`；向量结束 → 截断（不输出）；Flag → 只输出键；
  - GT → `0/1`、`0|1`、`./.`、单倍体 `0`；FILTER 下标 → 名字（空向量 → `.`）；QUAL 缺失 → `.`。
- 门：`clean.bcf` 的渲染正文与 `bcftools view -H clean.bcf` 逐字节一致（只排除 reed 自己加的 header 行）。

### 3.6 B3b：原生 BCF 手术（性能阶段，可选）

- 未编辑：整条记录字节原样拷贝。
- GT 掩码：在 GT 数组中把该样本的两个字节写成 `0x00, 0x00`（`bcf_gt_missing`，非定相 `./.`，与 `bcftools filter -S .` 一致）；单倍体只写第一个。
- INFO 重写：拷贝未删除的键值对；删除的跳过；重算的标签按 htslib 的最小宽度规则编码（`bcf_enc_vint`）追加或替换；
  更新 `n_info` 与 `l_shared`。
- 剔除样本：对每个 FORMAT 字段按"每样本步长 = 类型宽度 × 元素数"收集保留样本的切片，更新 `n_sample` 与 `l_indiv`。
- 门：与 B3a 路径的输出经 `bcftools view` 读回后逐字节一致。

### 3.7 B4：原生文本 → BCF 编码器（只在 B2 成为瓶颈时做）

按 htslib `vcf_parse` 的规则把文本直接编码成 BCF（字段类型来自 header，整数宽度用 `bcf_enc_vint` 规则），
跳过 noodles 的中间对象。门同 B2，另加"与 B2 输出字节一致"（同样的宽度选择规则时应当完全相同）。

---

## 4. 测试与验收门【规】

| 层 | 内容 |
|---|---|
| 单元 | typed 值编解码往返（proptest：所有类型 × 缺失 × 向量结束 × 长度 15 以上的描述符）；65,536 项 GT 表与 `GtCode::from_alleles` 逐项一致；长度前缀链：在每个可能的块边界切开同一份数据，定界结果相同 |
| 差分 | `BcfRecordRef` 的 `SiteFields`/`GenotypeSource` 输出 ≡ noodles `bcf::Record` 解码（`clean.bcf` 前 10 万条，全部字段） |
| 跨格式一致 | 同一数据的 VCF.gz 与 BCF（由 `bcftools view -Ob` 生成）：`reed stats` 输出逐字节相同；`reed filter -Ov` 输出正文逐字节相同（B3a 后）；`--summary` 相同 |
| BCF 输出 | `bcftools view out.bcf` 正文 ≡ `bcftools view out.vcf.gz` 正文；`bcftools index -s out.bcf` 正常；CSI 区域查询与流式 `-t` 一致；`-Ou` 输出可被 bcftools 读取 |
| 字典 | §3.3 的删除/新增标签用例 |
| 确定性 | threads ∈ {1, 2, 3, 7, 56}：BCF 输入的 stats/filter 输出 sha256 相同；BCF 输出文件 sha256 相同 |
| 性能 | `clean.bcf` 上 `reed stats --pi` ≤ VCF.gz 上的同一命令（预期更快：只解 GT 字节）；`filter -Ob` 相对 `-Oz` 的额外 wall ≤ 30%（B2 阶段） |

注意：先确认 `clean.bcf` 与 `clean.vcf.gz` 是同一份数据（`bcftools view -H` 两者 md5 比较）；如果不是，用 `bcftools view -Ob` 从 VCF 子集生成配对数据作为跨格式门。

---

## 5. 分期

| 期 | 内容 | 依赖 | 交付后可用的功能 |
|---|---|---|---|
| **B0** | `SiteFields`/`GenotypeSource`/`RecordRef` 抽象，VCF 实现迁移；类型化 `EngineError`（BCF 新增大量错误来源） | — | 无功能变化，门 = 现有输出逐字节不变 |
| **B1** | BCF 输入：header、字典、长度前缀链（`RecordFraming`、`FrameChain`）、`BcfRecordRef`、GT 表、数值列 | B0 | `reed stats`、`reed index` 支持 BCF 输入；`reed filter` 对 BCF 输入暂时明确报错（写出需要 B3a 或 B2） |
| **B2** | `-Ob/-Ou`：BCF 头写出（显式 IDX）、noodles 记录编码、BGZF 0 级、BCF CSI | B0 | VCF → BCF |
| **B3a** | BCF → 文本渲染（htslib 格式） | B1 | BCF → VCF；BCF → BCF（经渲染 + 编码） |
| **B3b** | 原生 BCF 手术（原样拷贝、GT 掩码、INFO 重写、剔除样本） | B1、B2 | BCF → BCF 的高速路径 |
| **B4** | 原生文本 → BCF 编码器 | B2 | 只在 profile 显示 B2 是瓶颈时做 |

建议顺序：B0 → B1 → B3a → B2 → B3b →（B4）。B1+B3a 完成后，"BCF 输入"的所有命令都可用；B2 完成后"BCF 输出"可用。

---

## 6. 风险

1. **字典错位**（§3.3）：最危险，必须有专门测试。
2. **浮点渲染**：BCF→VCF 的浮点文本必须与 htslib 一致，否则对拍会出现大量末位差异。复用 `format_g6` 并对 `clean.bcf` 全量对拍。
3. **VCF 4.4/4.5 特性**：如 `Number=P` 等新类型、超长字符串、本地等位（LAA）字段，BCFv2.2 能否表达以 htslib 行为为准；无法表达时报错，而不是静默改写。
4. **noodles 编码器与 htslib 的宽度差异**：只要 bcftools 读回语义相同即可接受；若需要字节级一致（例如与 bcftools 输出比 sha），在 B4 中按 `bcf_enc_vint` 规则实现。
5. **长度前缀链的延迟**：每个链节点要等前驱 chunk 解压完并走完链。所有 chunk 的解压是同时进行的，正常情况下等待只有微秒级；
   如果 `--timing` 显示 worker 在 `wait_start` 上的时间明显，就增大 `TICKETS_PER_WORKER`，让更多 chunk 同时在途来隐藏等待。
