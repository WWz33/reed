# doc_v10 实施状态

> 本页记录 `01-review.md` §1 缺陷清单与 `06-roadmap.md` R0–R3 的实施状态。
> 每条附证据（测试名 / 守卫 / 变更点）。复审基线（修复前的代码证据摘录）
> 见 §2。

## 实施状态（R0–R3）

| 项 | 严重度 | 状态 | 证据 |
|---|---|---|---|
| D10-1 GT 位宽 INT16 阈值错一倍 | **P0** | **已修** | 两处选宽合并为 `gt_width`（`native_encode.rs`）：等位 ≤ 62 → INT8，≤ 16382 → INT16，否则 INT32（按编码值域 `2a+3`）；`push_gt_code` 上界收紧为 1073741822。红测试 `gt_width_follows_the_encoded_value_domain`（修复前 allele 16383 得 INT16，修复后 INT32）；守卫 "GT width follows the encoded value domain" |
| D10-2 NULL 型描述符 panic | **P0** | **已修** | `Cursor::typed` 拒绝 `(NULL, count>0)`（`record.rs`），`(NULL,0)` 保留。红测试 `null_descriptor_with_nonzero_count_errors_without_panicking`（修复前 render.rs:525 越界 panic）；守卫 "null descriptors with a nonzero count are rejected" |
| D10-5 index 跨 chunk 记录失败 | **P0** | **已修** | worker 只承载（`CarriedLine { voff, bytes }` / `HeadEntry`），归约器 `IndexSink::feed` 拼接并**解码一次**；BCF 路径（`HeadEntry::Decoded`，FrameChain）不变。红测试 `index_handles_records_spanning_chunk_seams`（40 字节 chunk × 32 字节块 × 60 记录；修复前 "record line has fewer than 9 columns"；并行与顺序索引字节一致）；守卫 "the text index decodes carried fragments on the reducer" |
| D10-3 IDX 无界分配 | P1 | **已修** | `BcfDict.strings` 改 `BTreeMap<u32, Vec<u8>>`（稀疏）；`IDX_MAX = INT_MAX - 1` 之外丢弃并 warn、回退自动编号（htslib 语义）。测试 `huge_idx_does_not_allocate_and_falls_back`（`IDX=4e9` 回退、`IDX=1e6` 稀疏保留）+ `honors_explicit_idx` 回归；守卫 "the BCF dictionary is sparse and IDX-bounded" |
| D10-4 BGZF 最小块长缺失 | P1 | **已修** | `parse_block_size` 拒绝 `bsize + 1 < BGZF_HEADER_LEN`。红测试 `undersized_blocks_error_instead_of_panicking`（BSIZE ∈ {0,3,17} × 三条入口：流式 scan、流式 chunk 走查、mmap；修复前 `chunk/bgzf.rs:233` 切片越界 panic）；守卫 "undersized BGZF blocks are rejected" |
| D10-6 多成员 gzip | P2 | **已修** | `plain.rs` `GzDecoder` → `MultiGzDecoder`。测试 `multi_member_gzip_reads_every_member`（两成员拼接，两成员记录都在）；守卫 "plain gzip is decoded multi-member" |
| D10-7 PASS 声明缺失 | P2 | **已修** | `filter/mod.rs` 新增 `ensure_pass_filter`（输出头无 PASS 即注入 `Map::<Filter>::pass()`）。e2e `pass_is_declared_even_when_the_input_lacks_it`（无声明输入 `-Ob` 成功、头含 PASS、`dict.tag(0)==PASS`、记录渲染含 `FILTER=PASS`）；**htslib 真值门**：bcftools 1.24 读回 reed 写出的 BCF（3 记录 + 头声明 PASS）；守卫 "the output header declares PASS" |
| D10-8 计数乘法加固 | P3 | **已修** | 新增 `record.rs::field_bytes`（checked_mul 链），替换 6 处裸乘法（`render.rs` ×4、`native_encode.rs` ×2）。行为不变（全测试绿）；守卫 "BCF field sizes are computed with checked multiplication" |

## R0–R3 进度

| 步骤 | 状态 | commit |
|---|---|---|
| R0 红测试（a GT 边界 / b NULL 描述符 / c 畸形块 / d 多 chunk 索引） | 已完成（四红均验证在修复前红） | `5895fae` |
| R1 单行级修复（D10-1/2/3/4 + 守卫 1–4） | 已完成 | `4564c72` |
| R2 残片承载（D10-5 + 守卫 5） | 已完成 | `7bd8831` |
| R3 收尾（D10-6/7/8 + spec 回改 + 守卫 6–8 + STATUS） | 已完成 | 本提交 |

## 本轮未做（有意）

见 `06-roadmap.md` §4：`BcfIndexJob` 不动（已是 N5 模板）；读侧 PASS 注入
（`BcfDict`，【荐】项，htslib 写的 BCF 恒含 PASS）；`bcf_record_check` 不全量移植；
变异类型分类不复开。

## 复审基线（2026-10-14，**修复前**）

复审 commit：`850408a`（doc_v9 全部落地），工作区干净。
以下为修复前源码证据摘录，说明每条门修的是什么。

| 缺陷 | 基线证据 |
|---|---|
| D10-1 | `encode_gt_column`/`encode_gt_column_bcf` 阈值 62/32766；等位 16383 的编码值 32768 > 32767 → `as i16` 截断为 -32768。非 GT 的 `int_width` 按值判宽、正确——漂移面只在两个 GT 站点 |
| D10-2 | `(NULL,1)` INFO 字段：`render_site_part` `take(0)` → `fmt_array` → `int_element(0, &[])` 越界（实测 `render.rs:525` panic）。`(NULL,0)` 合法（flag 编码） |
| D10-5 | chunk 目标 1 MiB；文本 `IndexJob` 对 `view.head`/`view.tail` 残片 `parse_record` → `LineRef::parse` 失败 → `EngineError::job`（实测 40 字节 chunk 夹具）。对照：`BcfIndexJob` 的 `FrameChain` 承载已合规 |
| D10-3 | `IDX=4000000000` → `Vec::resize(4e9+1)` ≈ 96 GB → abort。htslib：`vcf.c:877-935` 解析时拒绝 `idx >= INT_MAX-1`；`bcf_hdr_set_idx`（`vcf.c:800-830`）同样 resize 到 idx+1，但 `hts_resize` 失败返回错误而非 abort |
| D10-4 | `BSIZE=0`（total=1）→ `read_block` `block[18..]` 越界 panic（实测）；total<4 → `block_inflated_size` `n-4` 下溢 panic。htslib：`bgzf.c:1219-1225` 下界检查 |
| D10-6 | `GzDecoder` 首成员后 EOF；htslib `bgzf.c:832-908` `Z_STREAM_END` 后 peek + `inflateReset`。BGZF 串联有测试保护，缺口仅 plain-gzip |
| D10-7 | 全仓无 PASS 注入（仅 `dict.rs` 测试夹具）；htslib 读写两侧强制（`vcf.c:1485-1490`/`:1724-1726`） |
| D10-8 | 乘积上界 2^57（`count ≤ 2^31-1`、`n_sample ≤ 2^24-1`、`width ≤ 4`），64 位不溢出；`take()` 兜底等价 htslib `end-ptr` 检查 |

## 测试与守卫基线

| 命令 | 结果 |
|---|---|
| `cargo test --workspace`（debug） | 全绿（含 v10 新增 5 个测试） |
| `cargo test --release --workspace` | 全绿 |
| `cargo clippy --release --workspace --all-targets -- -D warnings` | 0 警告 |
| `bash ci/check-doc-drift.sh` | 46 条全绿（含 8 条 v10 新增守卫） |
| htslib 真值门（bcftools 1.24） | `reed filter -Ob` 无 PASS 声明的输入 → bcftools `view` 读回 3/3 记录，头含 PASS |
