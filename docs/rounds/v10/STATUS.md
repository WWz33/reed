# doc_v10 实施状态

> 本页记录 `01-review.md` §1 缺陷清单与 `06-roadmap.md` R0–R3 的实施状态。
> 每条附证据（测试名 / 守卫 / 变更点）。起始状态：全部未修。
> 复审基线（修复前的代码证据摘录）见 §2，修复时逐条转为"已修"并附
> 红绿对照证据。

## 实施状态（R0–R3）

| 项 | 严重度 | 状态 | 证据 |
|---|---|---|---|
| D10-1 GT 位宽 INT16 阈值错一倍 | **P0** | **未修** | 现状：`native_encode.rs:524`/`:907` 阈值 `max_allele <= 32766`；`push_int` 的 `as i16`（`:1102`）静默截断；`push_gt_code` 上界（`:1212`）多放一个 |
| D10-2 NULL 型描述符 panic | **P0** | **未修** | 现状：`typed()`（`record.rs:251-264`）放行 `(NULL,>0)`；`fmt_array`/`fmt_gt` 经 `int_element`（`render.rs:521-527`）空 slice 越界 |
| D10-5 index 跨 chunk 记录失败 | **P0** | **未修** | 现状：`index.rs:165-170`/`:188-193` 急切解析残片；`run_seam` 忽略 `_line`（`:198-211`）；e2e 夹具单 chunk 测不到 |
| D10-3 IDX 无界分配 | P1 | **未修** | 现状：`dict.rs:66-69` `strings.resize(idx+1, None)`，idx 无上限 |
| D10-4 BGZF 最小块长缺失 | P1 | **未修** | 现状：`block.rs:55-67` 无下界检查；下游 panic 点 `chunk/bgzf.rs:233`（read_exact 切片）、`block.rs:45-48`（`n-4` 下溢，经 `chunk/bgzf.rs:146` / mmap `:448`） |
| D10-6 多成员 gzip | P2 | **未修** | 现状：`plain.rs:46` 单成员 `GzDecoder` |
| D10-7 PASS 声明缺失 | P2 | **未修** | 现状：输出头组装（`filter/mod.rs:469-487`）无 PASS 注入；`native_encode.rs:295` 遇 `FILTER=PASS` 记录报 undeclared |
| D10-8 计数乘法加固 | P3 | **未修** | 现状：`render.rs:160`/`:207`、`native_encode.rs:773`/`:823` 裸乘法（64 位无活 bug，见 01-review §1 D10-8 界论证） |

## R0–R3 进度

| 步骤 | 状态 |
|---|---|
| R0 红测试（a GT 边界 / b NULL 描述符 / c 畸形块 / d 多 chunk 索引） | 未开始 |
| R1 单行级修复（D10-1/2/3/4） | 未开始 |
| R2 残片承载（D10-5） | 未开始 |
| R3 收尾（D10-6/7/8 + spec 回改 + 守卫登记） | 未开始 |

## 本轮未做（有意）

见 `06-roadmap.md` §4（`BcfIndexJob` 不动、读侧 PASS 注入【荐】、
字典稀疏化【决】、变异类型不复开、`bcf_record_check` 不全量移植）。

## 复审基线（2026-10-14，**修复前**）

复审 commit：`850408a`（doc_v9 全部落地），工作区干净。
以下为修复前源码证据摘录，说明每条门修的是什么。

| 缺陷 | 基线证据 |
|---|---|
| D10-1 | `encode_gt_column`/`encode_gt_column_bcf` 阈值 62/32766；等位 16383 的编码值 32768 > 32767 → `as i16` 截断为 -32768。非 GT 的 `int_width`（`:1065-1073`）按值判宽、正确——漂移面只在两个 GT 站点 |
| D10-2 | `(NULL,1)` INFO 字段：`render_site_part` `take(0)` → `fmt_array` → `int_element(0, &[])` 越界。`(NULL,0)` 合法（`native_encode.rs:346` flag 编码） |
| D10-5 | chunk 目标 1 MiB（`chunk/mod.rs:24`）；文本 `IndexJob` 对 `view.head`/`view.tail` 残片 `parse_record` → `LineRef::parse` 失败 → `EngineError::job`。对照：`BcfIndexJob` 的 `FrameChain` 承载（`index.rs:244-285`）已合规 |
| D10-3 | `IDX=4000000000` → `Vec::resize(4e9+1)` ≈ 96 GB → abort。htslib：`vcf.c:877-935` 解析时拒绝 `idx >= INT_MAX-1`；`bcf_hdr_set_idx`（`vcf.c:800-830`）同样 resize 到 idx+1，但 `hts_resize` 失败返回错误而非 abort |
| D10-4 | `BSIZE=0`（total=1）→ `read_block` `block[18..]` 越界 panic；total<4 → `block_inflated_size` `n-4` 下溢 panic。htslib：`bgzf.c:1219-1225` 下界检查 |
| D10-6 | `GzDecoder` 首成员后 EOF；htslib `bgzf.c:832-908` `Z_STREAM_END` 后 peek + `inflateReset`。BGZF 串联有测试保护（`chunk/bgzf.rs:697-746`），缺口仅 plain-gzip |
| D10-7 | 全仓无 PASS 注入（仅 `dict.rs` 测试夹具）；htslib 读写两侧强制（`vcf.c:1485-1490`/`:1724-1726`） |
| D10-8 | 乘积上界 2^57（`count ≤ 2^31-1`、`n_sample ≤ 2^24-1`、`width ≤ 4`），64 位不溢出；`take()` 兜底等价 htslib `end-ptr` 检查 |

## 测试与守卫基线

| 命令 | 结果（基线，修复前） |
|---|---|
| `cargo test --workspace`（debug） | v9 记录的 52 个 test target 全绿（D10 各缺陷均无测试覆盖，除 `chunks_cover_all_records_with_many_seams` 间接证明 chunk 层可跨缝） |
| `cargo test --release --workspace` | 全绿 |
| `bash ci/check-doc-drift.sh` | v9 登记的 38 条全绿；本轮 6 条新守卫（01-review §3）待 R1–R3 落地后登记 |
