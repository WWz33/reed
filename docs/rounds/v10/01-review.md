# 01 BGZF/BCF/VCF 格式兼容性复审

> 复审对象：commit `850408a`（doc_v9 D9-1…D9-8 全部落地后的当前实现），工作区干净。
> 方法：reed 侧八文件通读 + htslib 引用抽查 + 每条缺陷的触发路径走查。
> 基准：`../../htslib`（bgzf.c、vcf.c、vcf.h）与 `../../bcftools` 的对应处理；
> 用户复审已给出的 htslib 行号中，三处承重引用经本轮独立复核属实。
> 本轮不修改任何代码；所有修法以规范形式给出（02），实施顺序见 06。
> 标记：【据】源码证据；【实】可复现触发路径。

---

## 0. 总评

v3–v8 把 BCF 的读、渲染、编码主路径建到了"与 `bcftools view` 字节兼容"
（`render.rs` 模块注释），v9 把引擎错误语义收干净。本轮以工作区 htslib
源码为基准逐条对照，发现的 8 处偏差有一个共同模式：**reed 在"自己生成的
数据"上全部正确，在"别人生成的边界数据"上出错**——大等位编号、畸形描述符、
大 IDX、小 BSIZE、跨界记录、串联 gzip、无 PASS 声明的头，全是 htslib
路径日常能遇到、reed 路径今天会坏的情形。

分档：

| 档 | 缺陷 | 后果形态 |
|---|---|---|
| P0 | D10-1 GT 位宽阈值错一倍 | **静默数据损坏**（写出的 BCF 基因型错误，任何读者都会读到坏数据） |
| P0 | D10-2 NULL 型描述符 | **进程 panic**（畸形 BCF 即崩，而非报错） |
| P0 | D10-5 index 跨 chunk 记录 | **功能必失败**（大样本文件几乎每个 chunk 边界都切在行内） |
| P1 | D10-3 IDX 无界分配 | **OOM abort**（恶意/损坏头） |
| P1 | D10-4 BGZF 最小块长缺失 | **进程 panic**（畸形块） |
| P2 | D10-6 多成员 gzip | 静默截断（串联 gzip 只读第一个成员） |
| P2 | D10-7 PASS 声明缺失 | 兼容输入报错 + 输出头缺行 |
| P3 | D10-8 计数乘法加固 | 64 位上无活 bug，平台无关化【荐】 |

已核减的说法（不计入缺陷）：

- **审核意见第 13 条"计数字段乘法会 wrap/越界"**：64 位目标上不成立。
  reed 的 `count` 经 `typed_int` 上界 `i32::MAX`（`record.rs:271-277` +
  `:284-294`），`n_sample` 是 24 位字段（`render.rs:71`），`type_width ≤ 4`，
  三处乘法（`render.rs:160`/`:207`、`native_encode.rs:773`/`:823`）的乘积
  ≤ 2^57，无溢出；越界由 `Cursor::take`（`record.rs:195-206`）兜底报
  `MalformedRecord`，与 htslib `vcf.c` 的 `end - ptr < bytes` 检查效果等价。
  降级为加固项 D10-8。

---

## 1. 缺陷清单

### D10-1 GT 位宽 INT16 阈值错一倍：大等位编号被静默截断【据】【实】

- **严重度：P0**。后果形态：写出的 BCF 数据损坏（无报错、无 panic）。
- **reed 现状**：
  - `native_encode.rs:522-528`（`encode_gt_column`，文本源）与
    `:905-911`（`encode_gt_column_bcf`，BCF 源）两处选宽：
    `max_allele <= 62 → INT8`，`max_allele <= 32766 → INT16`，否则 INT32。
  - 写出走 `push_int`（`native_encode.rs:1099-1105`）：width=2 时
    `value as i16`——Rust 的 `as` 截断**不 panic**。
- **错误论证**：BCF 基因型编码值是 `(allele+1)<<1 | phased`
  （`gt_codes_into`，`native_encode.rs:1165-1203`）。等位编号 `a` 的编码值
  上界是 `2a+3`。INT16 合法值域是 `[-32760, 32767]`
  （`record.rs:55-57`，htslib `BCF_MIN/MAX_BT_INT16`）：
  `2a+3 ≤ 32767 ⟺ a ≤ 16382`。reed 的阈值 32766 是值域上限减一，
  把"编码值的上界"误当成了"等位编号的上界"（INT8 阈值 62 恰好是对的：
  `2·62+3 = 127`，说明两处阈值一个对了一个错了）。
  等位编号 16383…=32766 的位点（即 ALT 数 > 16383）会走进 INT16 分支，
  编码值 32768…65535 被 `as i16` 截断：如等位 20000 → 编码值 40003 →
  截断为 -25533，解码端（reed 自己或 htslib）读到的是完全不同的基因型
  或哨兵值。
- **htslib 对照**：`bcf_enc_vint`（`vcf.c:2955-2993`）按**编码后的值**
  扫 min/max 判宽：`max ≤ BCF_MAX_BT_INT16(0x7fff)` 才进 INT16
  （`vcf.h:1539`）；单值走 `bcf_enc_int1`（`vcf.h:1649-1655`），同规则。
  等位 16382 → 编码值 32766/32767 → INT16；等位 16383 → 32768 → INT32。
- **附带问题**（同函数族，一并修）：`push_gt_code`（`native_encode.rs:1207-1216`）
  的合法上界 `(0..=i64::from(i32::MAX >> 1))`（`:1212`）多放一个等位：
  `a = 1073741823` 时 `(a+1)<<1 = 2^31`，`as i32` 同样回绕到 `i32::MIN`。
  上界应为 `a ≤ (i32::MAX - 1) >> 1` 即编码值 `(a+1)<<1 ≤ i32::MAX`。
- **正确先例（本仓库）**：非 GT 整数向量走 `int_width`
  （`native_encode.rs:1065-1073`），按值域判宽，与 htslib 一致——
  修法只需把两处 GT 站点对齐到同一规则。
- **修法**：阈值改为按编码值（或等价地等位阈值 62 / 16382）；
  `push_gt_code` 上界收紧一格。见 02 §1。
- **验收**：边界编码测试（等位 16382→INT16、16383→INT32、
  32766→INT32）；红绿对照（放回原阈值必红）。
  **必读性后果**：静默损坏是最坏的一类缺陷——下游统计在坏数据上
  照样跑出"看起来合理"的数。

### D10-2 NULL 型描述符（type=0, count>0）一路 panic【据】【实】

- **严重度：P0**。后果形态：进程 panic（畸形输入）。
- **reed 现状**：`Cursor::typed`（`record.rs:251-264`）只拒绝
  "未知类型"（`type_width(ty)==0 && ty != NULL`，`:260-262`），
  `(NULL, count>0)` 原样通过。此后：
  - INFO 路径：`render.rs:160` `c.take(type_width(0) * count)` = `take(0)`
    得空 slice，`:161` `fmt_array` → 非 CHAR 分支 width=0，
    `int_element(ty=0, &[])`（`render.rs:521-527` 的 `_` 分支）对空 slice
    取 `elem[0..4]` → **越界 panic**。
  - FORMAT 路径：`render.rs:207` 同理 `take(0)`；样本循环 `:235-243`
    对空 data 调 `fmt_gt`/`fmt_array1`，同样落入 `int_element` → panic。
  - FILTER 路径：`render.rs:130` `int_value(f_ty=NULL)` 返回
    `Err`（`record.rs:292`）——这条路径恰好是安全的。
  - 原生解码路径安全：`decode_format_columns` 有 `width == 0 → continue`
    （`render.rs:319`）；`gt.rs` 慢路径有 `matches!(ty, INT8|INT16|INT32)`
    检查（`gt.rs:166-171`）。**所以缺口只在渲染/文本导出路径。**
- **htslib 对照**：三层防护——`bcf_type_shift[BCF_BT_NULL]=0`
  （`vcf.c:92`，字节数恒 0）；`bcf_record_check` 明确把
  `type == BCF_BT_NULL && num > 0` 判 `BCF_ERR_TAG_INVALID`
  （`vcf.c:2207-2214`/`:2231-2236`），报 "Invalid INFO type 0 (null)"；
  格式化侧 `case BCF_BT_NULL` 只输出 `.` 不读数据（`vcf.c:6489`）。
- **修法**：在 `typed()` 拒绝 `(NULL, count>0)`（`MalformedRecord`）；
  `(NULL, 0)` 是合法 flag 写法（reed 编码侧 `encode_size(s, 0, NULL)`
  即此用法，`native_encode.rs:346`），保留。渲染侧无需改动。
  见 02 §2。
- **验收**：畸形 BCF 单测——手工拼 `(NULL, count=1)` 的 INFO 字段与
  FORMAT 字段各一，断言返回 `MalformedRecord` 而非 panic
  （`catch_unwind` 断言无 panic）。

### D10-5 `reed index` 文本路径把跨 chunk 记录残片当完整记录解析【据】【实】

- **严重度：P0**。后果形态：真实负载下功能必失败。
- **reed 现状**：文本 `IndexJob::run_chunk` 在 worker 端急切解析
  - `view.head`（`index.rs:165-170`）：chunk 首行的残片（上一 chunk
    断在行内时）被 `parse_record` 解析——残片不足 9 列时
    `LineRef::parse` 失败 → `EngineError::job("record line has fewer
    than 9 columns")` → 整个任务失败；
  - `view.tail`（`index.rs:185-193`）：chunk 末行的残片同理。
  归约器侧确实有 `prev_continues` 守门（`index.rs:431-434`），
  到那时才丢弃残片对应的 entry——**但 worker 已经先死了**。
  `run_seam` 收到引擎拼接好的完整行却忽略它（`index.rs:198-211`，
  `_line` 未用），注释说"the ordered sink already has it"——
  "it" 指的是被急切解析的 tail 残片，正是会失败的东西。
- **触发规模**：chunk 目标 1 MiB 压缩（`chunk/mod.rs:24`），边界按
  块切、与行无关。常规 VCF 行短、命中概率低；**大样本文件**
  （行宽数百 KB–数 MB，如 100k 样本）几乎每个边界都切在行内，
  → 索引必然失败。既有 e2e 夹具（`e2e_index.rs` 的 `rows()`，
  几十行 × 短行）总共几百字节、单 chunk，测不到；
  `chunks_cover_all_records_with_many_seams`（`chunk/bgzf.rs:567-634`）
  只数记录数，不经 `IndexJob`。
- **htslib 对照**：htslib 读取是纯流式，`bgzf_read` 在块边界内部
  自动补读，记录永远不跨块暴露给索引器；`bcf_index_build`/
  `tbx_index_build` 只消费流式记录 + `bgzf_tell`。不存在"残片"概念。
- **正确先例（本仓库）**：同文件 `BcfIndexJob`（`index.rs:216-315`）
  已是对的——`RecordFraming::Raw` 下 job 自己走长度前缀
  （`walk_bcf_records`，`:319-325`），跨界记录的前 24 字节 + 起始
  voff 经 `FrameChain` 传给下一 chunk（`:244-285`），到齐后
  `entry_at` 解析一次（`:349-378`）。
- **修法**：文本路径对齐 BCF 路径——worker 只**承载**残片
  （字节 + 文本偏移 + 预映射的 voff），**不解析**；解析推迟到
  记录完整之后（引擎拼好的 seam 行在 `run_seam` 收到，或归约器确认
  `!prev_continues` 后对 head 残片解析）。承载方式见 02 §5 的【决】。
- **验收**：多 chunk 夹具（小块强制出缝，仿
  `chunks_cover_all_records_with_many_seams` 的 32/40 字节切法）
  跑 `run_index` 成功；且与单线程/单 chunk 索引**字节一致**
  （复用 `index_bytes_are_thread_count_independent` 的对拍模式）。
  seam 记录若畸形，必须走 v9 W2 出口报错返回（不挂死、不静默跳过）。

### D10-3 头部 `IDX=` 直接驱动 Vec 分配【据】

- **严重度：P1**。后果形态：OOM abort（恶意或损坏的 BCF 头）。
- **reed 现状**：`BcfDict::from_header`（`dict.rs:26-73`）解析头部
  行里的 `IDX=`（`:58-65`），随后 `dict.strings.resize(idx + 1, None)`
  （`:66-68`）。idx 无上限：`IDX=4000000000` → Vec 分配约
  `4e9 × 24 B ≈ 96 GB` → 分配失败 abort（Rust 分配失败即 abort，
  不可恢复）。负值/非数字回退自动编号（`:61-65`），这条安全。
- **htslib 对照**：IDX 解析做 `strtol` + 范围检查
  （`idx < 0 || idx >= INT_MAX-1 → warning "Error parsing the IDX tag,
  skipping"`，`vcf.c:877-891`/`:925-935`）。**修正**：用户复审称
  "htslib 不用 IDX 做数组下标"——本轮核对 `bcf_hdr_set_idx`
  （`vcf.c:800-830`）后更正：htslib 同样把字典稠密数组 `hts_resize`
  到 `idx + 1`（`new_n = idinfo->id >= hdr->n ? idinfo->id+1 : ...`，
  `vcf.c:1041` 调用）。它的防护是两条：解析时拒绝
  `idx >= INT_MAX - 1`；`hts_resize` 失败返回 -1 由调用方报错
  （不是 abort）。reed 的缺口是：① 无范围检查；② `Vec::resize`
  分配失败即 abort。
- **修法**：① 解析时拒绝 `idx >= i32::MAX - 1`（丢弃该 IDX 并 warn，
  回退自动编号，语义同 htslib 的 skipping）；② 字典改稀疏存储
  （`HashMap<usize, _>`）使 IDX 不再驱动按值大小的分配——推荐，
  顺带消掉 `next` 游标状态。`honors_explicit_idx`（IDX=7）必须原样
  通过，所以界不能收到"已见定义数"。判据与取舍见 02 §3。
- **验收**：`IDX=10^9` / `IDX=4000000000` 的头解析不再触发大分配
  （行为测试 + 现有 `honors_explicit_idx` 回归）；损坏 IDX 回退
  自动编号有 warn。

### D10-4 BGZF 缺最小块长校验，畸形 BSIZE 下游 panic【据】【实】

- **严重度：P1**。后果形态：进程 panic（畸形块）。
- **reed 现状**：`parse_block_size`（`block.rs:55-67`）校验 magic /
  FEXTRA / BC 子字段后返回 `Some(bsize + 1)`，**没有**最小块长检查。
  `total < 18` 时两条下游路径 panic：
  - 流式 `read_block`（`chunk/bgzf.rs:202-242`）：`block.resize(total)` 后
    `self.inner.read_exact(&mut block[BGZF_HEADER_LEN..])`（`:233`）——
    `total < 18` 时切片起点越界 → panic；
  - `fill_to_target`（`chunk/bgzf.rs:128-150`）与 mmap 源
    （`:423-449`）都会拿 `total` 大小的块调
    `block_inflated_size`（`block.rs:45-48`）：`total < 4` 时 `n-4`
    下溢 → panic；`4 ≤ total < 8` 时不 panic 但读出垃圾 ISIZE，
    之后按错误的 `pos += total` 错位，下一次 `parse_block_size`
    报 "not a BGZF block header"（错误但非崩）。
  inflate 路径（`inflate.rs:42-47`）有 `total != block.len()` 检查，
  安全——缺口只在 chunk 走查路径。
- **htslib 对照**：`bgzf_read_block` 读到 BSIZE 立刻检查
  `block_length < BLOCK_HEADER_LENGTH(18)` →
  "Invalid BGZF block length" 报错返回（`bgzf.c:1219-1225`，本轮抽查）。
  注意 htslib 并不校验"块长 ≥ 头+尾"以外的上界（上界由 64 KiB 解压
  缓冲兜住，reed 的 `BGZF_INFLATE_MAX` 同，`inflate.rs:49-55`）。
- **修法**：`parse_block_size` 加 `bsize + 1 < BGZF_HEADER_LEN → None`。
  一处 if，三处下游全部自然落到已有的 `MalformedBgzf` 错误分支。
  见 02 §4。
- **验收**：`BSIZE ∈ {0, 3, 17}` 构造的流，在流式 scan_header、
  流式 next_chunk、mmap 三条入口都报 `MalformedBgzf` 而非 panic。

### D10-6 多成员 gzip 只读第一个成员【据】

- **严重度：P2**。后果形态：静默截断。
- **reed 现状**：`PlainChunkSource::open_gzip`（`plain.rs:45-48`）用
  `flate2::read::GzDecoder`——单成员解码器，第一个成员结束后返回
  EOF。`cat a.vcf.gz b.vcf.gz`（htslib 下两文件全读）在 reed 只处理 a，
  无警告。输入分流见 `input.rs:235-244`（`Format::VcfGzip`）。
- **htslib 对照**：普通 gzip 走 `inflate_gzip_block`（`bgzf.c:832-908`）：
  `Z_STREAM_END` 后 peek 是否还有数据，有则 `inflateReset` 续读下一
  成员；BGZF 路径遇空块只 continue（`bgzf.c:1239-1241`）。
- **边界说明**：BGZF 路径的串联有测试保护（`empty_members_anywhere_
  never_end_the_stream`、`concatenated_files_read_through_every_join`，
  `chunk/bgzf.rs:697-746`）；缺口**只在** plain-gzip 路径。
- **修法**：`GzDecoder` 换 `flate2::read::MultiGzDecoder`（等价于
  htslib 的 reset 循环）。见 02 §6。
- **验收**：两个 gzip VCF 拼接的 e2e（filter 记录数 = a + b）；
  单成员 gzip 回归（现有测试即覆盖，防换错类型）。

### D10-7 输出头不补 PASS 声明【据】

- **严重度：P2**。后果形态：兼容输入报错 + 输出头缺行。
- **reed 现状**：filter 输出头组装（`filter/mod.rs:469-487`）=
  输入头剪样本 + `ensure_info_defs`（只加 INFO 定义，`:830-890`），
  全仓无 PASS 注入（PASS 仅出现在 `dict.rs` 测试夹具）。输入头
  无 `##FILTER=<ID=PASS>` 声明时（VCF 规范允许不声明；htslib 读写
  两侧都自动补）：
  - 写 BCF：`FILTER=PASS` 的记录撞 `Self::undeclared("FILTER", "PASS")`
    （`native_encode.rs:291-296`，`:230-235`）→ 任务失败；
  - 写出的 BCF 头文本缺 PASS 行 → 下游 htslib 读到的是"不完整的头"
    （它能自补，但 `bcf_header_text` 写出的显式 `IDX=` 编号与
    htslib 的"PASS 占第 0 位"假设不一致，组合行为不可依赖）。
- **htslib 对照**：读头 `bcf_hdr_parse` 无条件先把 PASS 塞进字典
  （`vcf.c:1485-1490`，注释 "The filter PASS must appear first in the
  dictionary"）；新建写模式头同样先 append（`vcf.c:1724-1726`）。
- **修法**：输出头组装时若无 PASS 声明，补
  `##FILTER=<ID=PASS,Description="All filters passed">` 于 FILTER
  定义之首（与 htslib 的"PASS 在字典第 0 位"对齐——reed 的
  `bcf_header_text` 按 noodles StringMaps 顺序写显式 IDX，PASS 放
  首位即得 idx 0）。见 02 §7。
- **验收**：无 PASS 声明的输入跑 `-Ob` 不再报错；写出的 BCF 头含
  PASS 行且 PASS 的 IDX=0；htslib 侧用工作区 `bcftools view` 可读
  （真值门）。读侧 parity（`BcfDict` 不注入 PASS）列为【荐】，
  因为 htslib 写的 BCF 恒含 PASS，缺口只在理论输入上。

### D10-8 计数字段乘法：加固而非修复【据】

- **严重度：P3（加固）**。64 位目标上**无活 bug**，本节记录核减理由
  与加固规范，防止日后有人把它当 bug 重报或当已修忽略。
- **论证**：`count` 上界 `i32::MAX`（扩展计数走 `typed_int`，
  `record.rs:271-294`）；`n_sample` 上界 `2^24-1`（`render.rs:71`）；
  `type_width ≤ 4`（`record.rs:30-37`）。乘积上界
  `4 × 2^31 × 2^24 = 2^57 < 2^64`，不溢出。所有消费端经
  `Cursor::take`（`record.rs:195-206`）越界报错，行为等价 htslib
  `end - ptr < bytes`（`vcf.c`）与 FORMAT 的 2 GiB 上限（收紧但不矛盾）。
- **htslib 对照**：`vcf.c:2236-2241` 的 `uint64_t ndata = num * n_sample;
  ndata > (INT_MAX >> shift) → too_many_indiv` 是 32 位 `int` 尺寸
  世界的溢出防护；reed 的 usize 世界不需要同样的上限，但需要
  乘法不依赖"64 位"这个隐含前提。
- **修法**：`render.rs:160`/`:207`、`native_encode.rs:773`/`:823` 四处
  `type_width(ty) * count [* n_sample]` 改 `checked_mul` 链，溢出即
  `MalformedRecord`。先例 `gt.rs:287-291`。见 02 §8。
- **验收**：32 位目标（若 CI 覆盖）或注入极端 count 的畸形记录测试
  不 panic；64 位行为不变（现有测试全绿）。

---

## 2. 已核减 / 已消解的说法

| 说法 | 处置 | 依据 |
|---|---|---|
| 审核 #13"计数字段乘法 wrap/越界"是 bug | **核减**为加固项 D10-8 | §1 D10-8 的界论证；`take()` 兜底 |
| 审核 #1"reed index 多 chunk 失败"是格式处理差异 | **修正理解**：这是 reed 并行分块架构的结构性问题（htslib 纯流式，无从差异），修法是把 `BcfIndexJob` 已有模式对齐到文本路径 | §1 D10-5；`index.rs:216-315` 对照 `:163-211` |
| 审核 #14 只是"阈值口径"兼容性差异 | **上调**为 P0 数据损坏：INT16 分支内 `as i16` 静默截断，任何下游读者读到坏基因型 | §1 D10-1 |
| 审核 #16"渲染侧对 NULL 直接跳过" | **微调**：渲染侧无需改——`(NULL,>0)` 在 `typed()` 拒绝后到不了渲染；`(NULL,0)` 在 INFO 路径已被 `count == 0 → continue`（`render.rs:156-158`）正确处理 | §1 D10-2 |

---

## 3. 复核门（本轮新增，修复落地后由 `ci/check-doc-drift.sh` 登记）

1. **D10-1 守卫**：`native_encode.rs` 两处 GT 阈值出现 `16382`
   （放回原值 32766 即 DRIFT）。
2. **D10-2 守卫**：`record.rs` 的 `typed()` 含 `NULL` 且 `count` 的
   组合拒绝分支。
3. **D10-4 守卫**：`block.rs` 的 `parse_block_size` 含
   `BGZF_HEADER_LEN` 下界检查。
4. **D10-6 守卫**：`plain.rs` 使用 `MultiGzDecoder`。
5. **D10-7 守卫**：`filter/mod.rs`（或其实现落点）含 PASS 注入。
6. **D10-5 守卫**：`index.rs` 文本 job 的 `run_seam` 实际使用其
   `_line` 参数（残片不再在 worker 端解析）。
