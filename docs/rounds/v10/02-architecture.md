# 02 修订规范：BGZF/BCF/VCF 格式兼容性

> 目标：**reed 在"别人生成的边界数据"上与 htslib 行为一致**——
> 畸形输入报错（不 panic、不 OOM）、大等位编号不损坏、串联 gzip 全读、
> 输出头自带 PASS、跨界记录可索引。
> 本轮是缺陷修复轮，不重开 v6 的 Y 系列、v7 的 Z 系列、v8 的 P 系列、
> v9 的 Q 系列结构项；引擎终止语义（v9 W1/W2）本轮直接复用。
>
> 标记：【规】必须；【荐】推荐；【决】需拍板；【据】有源码证据。

---

## 0. 现状结构问题

| # | 问题 | 证据 |
|---|---|---|
| S1 | GT 选宽用"等位编号阈值"，与全仓库其它选宽（`int_width` 按值）和 htslib（`bcf_enc_vint` 按编码值）都不一致 | 01-review §1 D10-1 |
| S2 | `typed()` 是描述符合法的唯一关卡，但不拒绝 `(NULL, >0)`，把非法输入放行到多个消费端 panic | D10-2 |
| S3 | 字典用稠密 `Vec<Option<_>>` 以 IDX 为下标，IDX 来自不可信头部文本却无界 | D10-3 |
| S4 | `parse_block_size` 是块走查的唯一关卡，但不校验最小块长，三个调用点各自暴露panic/错位 | D10-4 |
| S5 | 文本 `IndexJob` 在 worker 端解析残片，与"残片归约侧裁决"的注释意图相反；`run_seam` 收完整行却不用 | D10-5 |
| S6 | plain-gzip 用单成员解码器 | D10-6 |
| S7 | 输出头组装只加 INFO 定义，FILTER 的 PASS 没有 htslib 那样的强制存在 | D10-7 |
| S8 | 四处计数乘法依赖"64 位 usize"隐含前提 | D10-8 |

---

## 1. GT 选宽按编码值【规】（修 D10-1）

> **规约 N1**：GT 列的整数宽度由**编码值的值域**决定，不由等位编号的
> 值域决定。编码值 `v = (allele+1)<<1 | phased`，值域
> `[-32760, 32767]`（INT16）/`[-120, 127]`（INT8）。

具体条款【规】：

1. `encode_gt_column`（`native_encode.rs:522-528`）与
   `encode_gt_column_bcf`（`:905-911`）的阈值改为
   `max_allele <= 62 → INT8`、`max_allele <= 16382 → INT16`、
   否则 INT32。推导：`2a+3 ≤ 127 ⟺ a ≤ 62`；`2a+3 ≤ 32767 ⟺ a ≤ 16382`。
2. **更稳的等价实现（推荐）**：像 `int_width` 那样扫编码值的
   min/max 按 `INT8_MIN_VALUE..INT8_MAX_VALUE` /
   `INT16_MIN_VALUE..INT16_MAX_VALUE` 判宽。两种实现选一即可，
   实现者按局部可读性决定，但阈值必须满足 N1 的等价关系，
   且**只改一处helper**供两个调用点使用（消除"一个对一个错"
   的漂移土壤）。
3. `push_gt_code`（`:1207-1216`）的合法上界收紧为
   `0..=(i32::MAX >> 1) - 1`（= 1073741822），保证
   `(a+1)<<1 | bit ≤ i32::MAX`（`a = 1073741823` 是首个越界者）；
   越界等位仍按现语义落为缺失（bare phasing bit），不 panic。
4. 解码侧（`gt.rs`/`render.rs`）**不改**：坏数据从来只可能来自
   本缺陷的写出侧；解码按渲染值走，天然免疫。

验收：等位 62/63、16382/16383、32766/32767 六组边界的编码测试
（断言描述符类型与字节）；修复前阈值做红绿对照（放回 32766 必红）。

## 2. NULL 型描述符非法【规】（修 D10-2）

> **规约 N2**：`(type == NULL, count > 0)` 是非法描述符，
> 必须在读取关卡 `Cursor::typed` 拒绝，返回 `MalformedRecord`；
> `(NULL, 0)` 合法（flag 的编码形态），放行。

条款【规】：

1. `typed()`（`record.rs:251-264`）在现有未知类型检查旁加
   `ty == value_type::NULL && count > 0 → Err(malformed("BCF null
   type with a nonzero count"))`。
2. 渲染侧（`render.rs` 的 `fmt_array`/`fmt_array1`/`fmt_gt`）与
   原生解码侧**不加 NULL 特判**——N2 保证非法输入到不了它们；
   `(NULL, 0)` 在 INFO 路径已被 `count == 0 → continue`
   （`render.rs:156-158`）正确处理，FORMAT 样本循环对 count=0
   输出空（现状行为，不改）。
3. 错误文案对齐 htslib 语义（"Invalid INFO type 0 (null)" 一族），
   便于对日志。

验收：手工拼 `(NULL, count=1)` 的 INFO 与 FORMAT 记录各一条，
`catch_unwind` 断言 `render_record` 返回 `Err` 且无 panic；
`(NULL, 0)` 的 flag 记录渲染回归（现有编码侧 `encode_size(s,0,NULL)`
   的 fixture 即覆盖）。

## 3. IDX 限界【规】（修 D10-3）

> **规约 N3**：头部文本里的 `IDX=` 必须过与 htslib 等价的范围检查，
> 且任何 IDX 驱动的分配失败必须是可报告的 `Header` 错误，不是 abort。

**修正说明【据】**：用户复审称 htslib "不用 IDX 做数组下标"，本轮核对
`bcf_hdr_set_idx`（`vcf.c:800-830`）后更正：htslib 同样把字典稠密数组
`hts_resize` 到 `idx + 1`（`new_n = idinfo->id >= hdr->n ? idinfo->id+1 : ...`）。
它的防护是两条：**`idx >= INT_MAX - 1` 在解析时即被拒**
（`vcf.c:877-891`，warn "Error parsing the IDX tag, skipping"），以及
`hts_resize` 失败时返回 -1（调用方报错）。reed 的缺口不在"用了数组下标"，
而在：① 无 `INT_MAX - 1` 量级的范围检查；② `Vec::resize` 失败即 abort。
`honors_explicit_idx`（`dict.rs:180-188`，IDX=7）证明稀疏小 IDX 必须保留，
因此"限到已见定义数"的表述作废——**界取 `idx >= INT_MAX - 1` 拒绝**。

条款【规】：

1. `BcfDict::from_header` 解析 IDX 时：`idx >= i32::MAX as usize - 1`
   即视为损坏——该行的 IDX **丢弃并 warn**（回退自动编号），语义对齐
   htslib 的 "Error parsing the IDX tag, skipping"。
   `honors_explicit_idx`（IDX=7）与 `builds_the_dictionary_in_header_order`
   必须原样通过。
2. **分配失败可报告**：字典 Vec 的 `resize` 不得以裸 `vec![None; n]` /
   `resize` 触发 abort。落法二选一：
   - **推荐**：`strings` 改 `HashMap<usize, Vec<u8>>`（或
     `BTreeMap<usize, _>`），`tag()` 改查表——IDX 不再驱动任何按值大小的
     分配，N3 的两条同时消解，且 `key_of` 的线性扫描一并消失（每 run 一次
     的冷路径，行为不变）；
   - 备选：保持 Vec 但在 resize 前用 `try_reserve`，失败返回
     `IoError::Header`（大 IDX 仍会分配，只是不再 abort）。
   **拍板推荐前者**：它把"IDX 是稀疏的"这一事实编码进类型，删掉
   `next` 游标与 resize 两条状态。
3. 现有行为保留：负/非数字 IDX 回退自动编号（htslib 对解析失败同样
   skip+warn；reed 现状静默，属可接受口径差，不在本轮扩大）。

验收：`IDX=10^9` 的头解析返回错误或回退（视 `idx >= INT_MAX-1` 判定）；
`IDX=4000000000` 同；`honors_explicit_idx` 回归绿；损坏 IDX 行产生 warn。

## 4. BGZF 最小块长【规】（修 D10-4）

> **规约 N4**：`parse_block_size` 必须拒绝总块长小于
> `BGZF_HEADER_LEN`（18）的块，与 htslib
> `block_length < BLOCK_HEADER_LENGTH`（`bgzf.c:1219-1225`）对齐。

条款【规】：

1. `parse_block_size`（`block.rs:55-67`）在返回前加
   `if bsize + 1 < BGZF_HEADER_LEN { return None; }`。
   调用点已有统一的 "not a BGZF block header" /
   `MalformedBgzf` 处理，三处入口（流式 scan、流式 chunk、
   mmap）无需各自加判断。
2. 上界不查：htslib 同样不查块长上界，解压侧由
   `BGZF_INFLATE_MAX`（`inflate.rs:49-55`）与 libdeflate 的
   "insufficient space" 错误兜住，reed 现状与此一致，不改。

验收：`BSIZE ∈ {0, 3, 17}`（总块长 1/4/18）构造的字节流，
三条入口各跑一遍：流式 `scan_header`、流式 `next_chunk`、
mmap `next_chunk`，全部返回 `Err`（`MalformedBgzf`）且无 panic；
`BSIZE=17`（total=18，恰过下界但必然在后续 read_exact/ISIZE
校验失败）断言同样无 panic。

## 5. index 残片承载：worker 只搬运，解析一次【规】（修 D10-5）

> **规约 N5**：跨 chunk 记录的所有权归引擎/归约侧。worker 对
> head/tail 残片**只承载（字节 + 起始 voff），不解析**；一条记录
> 在其字节完整之后被解析**恰好一次**。

背景【据】：引擎已保证 seam 行拼接（`NewlineStitcher`，
`pipeline.rs` v9 落地）并把完整行交给 `run_seam`；归约器持有
`prev_continues`（`index.rs:442`）。缺的只是"把 voff 也传过去"
和"worker 别先解析"。BCF 侧的 `FrameChain` 模式
（`index.rs:244-285`）证明本仓库已有此范式。

条款【规】：

1. `IndexPartial.head`/`tail` 的类型从 `Option<RecordStart>` 改为
   承载体，如
   `struct CarriedLine { voff: u64, fragment: Vec<u8> }`
   （fragment 为残片字节；完整行时即行字节）。
   voff 在 worker 端用现成的 `VirtualMap`（`index.rs:82-102`）
   预映射——worker 有 `chunk.blocks`，归约器没有，这是 voff 必须
   在 worker 端算好的原因。
2. 归约器（`IndexSink::accept`）裁决：
   - `prev_continues == false` 时本 chunk 的 head 是记录起点——
     在归约侧 `parse_record`（解析错误即 v9 W2 出口，fatal）；
   - `prev_continues == true` 时 head 是 seam 记录的尾部残片，
     丢弃（该记录已由 seam 路径收录）；
   - tail 残片暂存为 `pending_fragment`；引擎拼好 seam 行后
     `run_seam` 用**暂存 fragment 的 voff** 解析该完整行一次
     （替代现在"忽略 `_line`、靠 tail 预解析"的逻辑）。
3. `run_seam` 的错误走 v9 的 reducer 错误出口（fatal，与 body
   一致）——正好满足 v9 的 E1"seam 与 chunk 体错误判定一致"。
4. **【决】承载位置**：上面把解析放归约侧（`accept`/`run_seam`
   都在 reducer 线程跑）是本推荐，理由是残片归属裁决本来就在
   归约侧（`prev_continues`）。备选：worker 解析但只限"本 chunk
   内完整的行"，残片仍走承载——实现量相近，但保留"worker 解析、
   归约侧也解析"两条解析路径的漂移面。**拍板推荐前者。**
5. `BcfIndexJob` 不动（它已合规，是本条的模板）。

验收：

- 多 chunk 夹具（64 字节小块 + 40 字节 chunk 目标，仿
  `chunks_cover_all_records_with_many_seams`）跑 `run_index`
  成功，`n_records` 与串行扫描一致；
- 与单 chunk（同输入一整块）索引**字节一致**——对拍模式同
  `index_bytes_are_thread_count_independent`；
- seam 行畸形（把一条记录从中截断拼成坏行）时任务返回错误
  而非挂死/静默跳过（v9 W2 已保证出口，本条验证语义）。

## 6. 多成员 gzip【规】（修 D10-6）

> **规约 N6**：plain-gzip 输入按多成员流解码，与 htslib
> `inflate_gzip_block` 的 `inflateReset` 循环（`bgzf.c:832-908`）
> 等价。

条款【规】：

1. `plain.rs:46` `GzDecoder` → `flate2::read::MultiGzDecoder`。
   flate2 的 Multi 语义即"成员结束后自动续下一成员"，与 htslib
   的 peek+reset 循环等价（含空成员跳过）。
2. `open_gzip` 的 rustdoc 更新（"serial inflate" 描述不变，
   但不再暗示单成员）。

验收：两个 gzip VCF（各含头 + 记录）`cat` 拼接后跑 filter，
输出记录数 = a + b；单成员 gzip 回归绿（既有测试）。

## 7. PASS 注入【规】（修 D10-7）

> **规约 N7**：输出头必须声明 PASS，且位于 FILTER 定义之首
> （htslib 语义：PASS 在字典第 0 位，`vcf.c:1485-1490`）。

条款【规】：

1. filter 输出头组装（`filter/mod.rs:469-487`，`out_header`
   定型处）加一步：若 `out_header.filters()` 无 `PASS`，插入
   `##FILTER=<ID=PASS,Description="All filters passed">`。
   实现走 noodles 的 header map 插入，与 `ensure_info_defs`
   同风格同位置。
2. 放首位的原因：reed 的 BCF 头写出（`bcf_header_text`，
   `encode.rs:183-222`）给每行写显式 `IDX=`，编号来自
   StringMaps 顺序；PASS 在首位即 `IDX=0`，与 htslib 读者假设
   （`vcf.c:1485` 的 "must appear first"）一致。
3. 文本输出（`-Ov/-Oz`）同样受益：`bcftools view` 对无 PASS
   声明的头会自补，reed 主动补使两边输出头一致。
4. **【荐】读侧 parity**：`BcfDict::from_header` 读头时若无
   PASS 亦注入（记录引用 idx 0 = PASS 的畸形-but-htslib-可读
   BCF 才不报错）。非本轮【规】：htslib 写的 BCF 恒含 PASS，
   触发面是理论输入；列为后续项。

验收：无 PASS 声明的输入（头部只声明 LowQual 之类）跑
`-Ob`：不再报 `undeclared("FILTER","PASS")`；写出的 BCF 头
含 PASS 行且 `IDX=0`；工作区 `bcftools view` 可读该 BCF
（真值门，写进出索引/查询亦正常）。

## 8. 计数乘法平台无关化【荐】（修 D10-8）

> **规约 N8（【荐】）**：记录字段的字节数计算不得依赖
> "usize 是 64 位"的隐含前提；溢出即报 `MalformedRecord`。

条款【荐】：

1. `render.rs:160`/`:207`、`native_encode.rs:773`/`:823` 四处的
   `type_width(ty) * count [* n_sample]` 改 `checked_mul` 链，
   `None → Err(malformed("BCF field size overflows"))`。
   先例：`gt.rs:287-291`。
2. 不改行为：64 位上这些乘法本就不会溢出（01-review §1 D10-8
   的界论证），本条只为 32 位目标与"读代码者的安心"。

验收：现有测试全绿（行为不变）；若 CI 有 32 位 target，跑通。

---

## 9. spec 回改清单

| 文档 | 改动 |
|---|---|
| `docs/spec/architecture.md` §6（BGZF） | 补 N4（最小块长）与 N6（多成员 gzip）两条款 |
| `docs/spec/architecture.md` §6.3 / `docs/rounds/v3/02-bcf.md` | 补 N1（GT 按编码值）、N2（NULL 非法）、N3（IDX 限界）、N8（乘法） |
| `docs/README.md` | v10 行已由本轮登记 |
| `ci/check-doc-drift.sh` | 按 01-review §3 的六条守卫登记 |

## 10. 验收纪律（沿用 v9 §8）

1. 双 profile 全绿：`cargo test --workspace`（debug）与
   `--release` 都过；clippy 不新增警告。
2. 每个 P0/P1 先红后绿；红测试必须**在修复前验证过会失败/会
   panic**（D10-1 放回原阈值、D10-5 放回急切解析）。
3. htslib 真值门：凡输出 BCF 的改动，用工作区
   `../../bcftools/bcftools view` 读回验证（D10-1/7 必做）。
4. 确定性：threads ∈ {1, 4} 输出 sha256 相同。
5. 性能不回退 10%；`bench/suites/e2e/results.tsv` 追加一行。
6. 同一 PR 删除被替代的旧实现（v2 06 §1 纪律 7）。
