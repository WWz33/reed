# reed doc_v10：BGZF/BCF/VCF 格式兼容性复审 + 修复规范（2026-10-14）

> 本目录是 `docs/rounds/v9/` 之后的**第十轮规范**。它由两部分组成：
> 对当前实现（v9 修复全部落地后）与 htslib/bcftools **格式处理行为**的
> 逐条对照复审（BGZF 块、gzip 成员、BCF 记录校验、字典、GT 编码、头部声明），
> 以及在复审基础上的**修复规范**。
> 本目录不修改任何代码，所有改动以规范形式给出，由实现者执行。
> 历史轮次 `docs/rounds/vN/` 只增不改；当前规范索引见 `docs/README.md`。

## 背景

v7–v9 三轮把 CLI 语义、线程预算、引擎终止语义修到了"错误必返回、返回必一致"。
本轮换一个基准：**工作区内的 htslib/bcftools 源码**（`../../htslib`、
`../../bcftools`，与用户复审所引行号同一棵树）。对照中发现 reed 的
格式处理有 8 处偏差，按后果分三档：

- **静默数据损坏 / 进程 panic / 功能必失败**（D10-1/2/5，P0）：
  GT 位宽阈值错一倍导致大等位编号基因型被 `as i16` 截断写坏；
  NULL 型描述符一路 panic 而不是报错；
  `reed index` 文本路径把跨 chunk 记录残片当完整记录解析，必然失败。
- **畸形输入可打崩 / 打爆内存**（D10-3/4，P1）：
  头部 `IDX=` 直接驱动 Vec 分配（恶意头 OOM abort）；
  BGZF 缺最小块长校验（BSIZE=0 即 panic）。
- **兼容性缺口 / 加固**（D10-6/7/8，P2/P3）：
  多成员 gzip 只读第一个成员；输出头不补 PASS 声明；
  计数字段乘法补 checked_mul 做平台无关化加固（64 位上不是活 bug）。

好消息：**同仓库里已有正确范式可抄**——`BcfIndexJob` 的跨界记录承载
（D10-5 的修法模板）、`int_width` 的按值判宽（D10-1 的非 GT 先例）、
`gt.rs` 的 checked_mul（D10-8 的先例）、BGZF 串联测试（证明 D10-6
只影响 plain-gzip 路径）。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 格式兼容性复审：逐条缺陷清单 D10-1…D10-8（reed `file:line` + htslib 对照行号、严重度、后果、修法）、已核减的说法、复核门 |
| [02-architecture.md](02-architecture.md) | 修订规范：GT 宽度按编码值（N1）、NULL 描述符非法（N2）、IDX 限界（N3）、BGZF 最小块长（N4）、index 残片承载（N5）、gzip 多成员（N6）、PASS 注入（N7）、checked_mul（N8）；每条带验收标准 |
| [06-roadmap.md](06-roadmap.md) | 修复路线 R0–R3：红测试先行、单行级修复、结构对齐、兼容性收尾；通用验收门 |
| [STATUS.md](STATUS.md) | 实施状态（起始全部未修）+ 复审基线（代码证据摘录） |

## 复审方法【据】

1. **双向定位**：从 reed 侧 8 个被点名的文件（`plain.rs`、`filter/mod.rs`、
   `record.rs`、`render.rs`、`native_encode.rs`、`dict.rs`、`block.rs`、
   `index.rs`）通读实际实现，逐条核对审核意见对 reed 现状的描述是否属实；
2. **htslib 抽查**：对承担结论重量的三处 htslib 引用独立复核
   （`bgzf.c:1219-1225` 最小块长检查、`vcf.c:2238-2245` FORMAT 乘法上限、
   `vcf.c:2955-2993` `bcf_enc_vint` 按值判宽），行号与代码均属实；
   其余 htslib 引用（`vcf.c:877-935`、`vcf.c:1485/1724`、`bgzf.c:832-908` 等）
   与用户复审所引一致，本轮未逐行重读；
3. **现状证据**：每条缺陷在 reed 源码中的精确位置、触发路径与下游后果
   均已在复审中走查（01-review.md §1 总表、§3 逐条）；
4. **核减**：审核意见第 13 条（计数字段乘法）在 64 位目标上不构成活 bug
   （`count ≤ i32::MAX`、`n_sample ≤ 2^24-1`、`width ≤ 4`，乘积 ≤ 2^57，
   越界由 `take()` 兜底报错），降级为加固项 D10-8。

复审基线：commit `850408a`（doc_v9 D9-1…D9-8 全部落地），工作区干净。

## 结论速览

| 维度 | 结论 |
|---|---|
| GT 位宽（D10-1，P0） | `native_encode.rs:524`/`:907` 的 INT16 阈值 `max_allele <= 32766` 错一倍：等位 16383…32766 的编码值 32768…65535 越过 INT16 合法上限 32767，`push_int` 的 `as i16` 静默截断 → **写出的 BCF 基因型损坏**。htslib 等价规则：编码值 ≤32767 ⟺ 等位 ≤16382。INT8 阈值 62 恰好正确，非 GT 整数走 `int_width`（按值）也正确，只有两处 GT 站点错 |
| NULL 描述符（D10-2，P0） | `typed()`（`record.rs:251-264`）接受 `(NULL, count>0)`；INFO 路径（`render.rs:160`）与 FORMAT 路径（`render.rs:235-241`）的 `fmt_array`/`fmt_array1`/`fmt_gt` 落到 `int_element` 对空 slice 取 `elem[0..4]` → **panic**。htslib `vcf.c:2231-2236` 判 `BCF_ERR_TAG_INVALID` 报错；`(NULL, 0)` 是合法 flag，不能一刀切 |
| index 跨 chunk 记录（D10-5，P0） | 文本 `IndexJob` 在 worker 端急切解析 `view.head`/`view.tail` 残片（`index.rs:165-170`/`:188-193`），跨 chunk 记录行必然解析失败 → **整个任务报错**；e2e 夹具单 chunk 测不到。同文件 `BcfIndexJob`（`:216-315`）已实现正确架构（FrameChain 携带跨界记录前缀 + voff，到齐后解析一次），文本路径对齐它即可 |
| IDX 分配（D10-3，P1） | `dict.rs:66-69` `strings.resize(idx+1, None)` 无上限：`IDX=4000000000` 的头 → 约 96 GB 分配 → abort。htslib 同样把字典稠密数组 resize 到 idx+1（`bcf_hdr_set_idx`，`vcf.c:800-830`），但解析时拒绝 `idx >= INT_MAX-1`（`vcf.c:877-891`）且分配失败返回错误；reed 两条都缺 |
| BGZF 最小块长（D10-4，P1） | `parse_block_size`（`block.rs:55-67`）缺 `bsize+1 < 18` 检查（htslib `bgzf.c:1219-1225` 有）：`total<18` 时 `read_block` 的 `read_exact(&mut block[18..])`（`chunk/bgzf.rs:233`）切片越界 panic；`total<4` 时 `block_inflated_size`（`block.rs:45-48`）`n-4` 下溢 panic（经 `chunk/bgzf.rs:146` 与 mmap `:448`） |
| 多成员 gzip（D10-6，P2） | `plain.rs:46` 用单成员 `GzDecoder`：`cat a.vcf.gz b.vcf.gz` 只处理 a（静默截断）。htslib `inflate_gzip_block` 在 `Z_STREAM_END` 后 peek 并 `inflateReset` 续读（`bgzf.c:832-908`）。BGZF 路径有串联测试保护（`chunk/bgzf.rs:697-746`），缺口只在 plain-gzip |
| PASS 声明（D10-7，P2） | 输出头组装（`filter/mod.rs:469-487`）不补 PASS；输入头无 PASS 声明（htslib 读写两侧自动补，`vcf.c:1485-1490`/`:1724-1726`）时，`-Ob` 遇 `FILTER=PASS` 记录撞 `undeclared("FILTER","PASS")` 报错（`native_encode.rs:295`），写出的 BCF 头也缺这行 |
| 计数乘法（D10-8，P3） | 64 位上不会溢出，`take()` 越界报错已等价 htslib 的 `end-ptr` 检查；`gt.rs:288` 已有 checked_mul 先例。补 `render.rs:160`/`:207`、`native_encode.rs:773`/`:823` 属平台无关化加固 |
