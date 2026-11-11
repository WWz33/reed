# v28 STATUS

> 更新 2026-11-10。P0/P1/P2 全部实施（V28-01…20）；V28-21…27 的下沉方向保留为设计项。
> “已修”只在有回归测试或对应验收证据时使用。

## 已实施

| ID | 状态 | 改动 | 验收 |
|---|---|---|---|
| V28-01 | **已修** | `bcf_header_text` 给 `##contig` 也写实际字典 `IDX=`；contig 用 `StringMaps::contigs()`，tag 继续用 `strings()`。 | `contig_idx_is_written_and_roundtrips`；`bcf_output_preserves_explicit_contig_idx`（VCF→BCF 回读 CHROM 正确，含 `-W` 索引后的 region 读）。 |
| V28-02 | **已修** | BCF 的 CSI 参考序列表按字典 tid 顺序构造；记录继续使用原生 BCF tid；文本输入仍用 header 顺序。 | `bcf_dictionary_tid_is_mapped_to_the_csi_header_order`：`idx.bcf` 的 chr1 region 只出两条 chr1；bcftools 同索引也能正确查询。 |
| V28-03 | **已修** | `RegionChunkSource::read_block` 区分 `filled == 0` 的合法 EOF 与 1–17 字节的短 header。 | `clean_eof_is_end_of_range_but_a_short_header_is_malformed`。 |
| V28-04 | **已修** | 新增 `record::gt_element`：vector-end 返回 `None`，其余负数报 `MalformedRecord`；render、table、native encoder 共用。 | `negative_gt_element_is_rejected_not_rendered_as_ref`；全 65,536 byte pair 的 table/text parity 继续通过。 |
| V28-05 | **已修** | integer-backed float FORMAT 列同时把 missing 和 vector-end 转 `NaN`。 | `char_typed_float_column_element_is_missing_not_a_panic` 扩展 INT8/16/32 vector-end 断言。 |
| V28-06 | **已修** | `OrderedOutput::finish_pending` / `PendingOrderedOutput`、`IndexSink` 暂存 index、`ReportSink` 暂存 `AtomicSink`；均在 `job.finish`/`FrameChain` 成功后发布。 | 新增截断 BCF 的 filter/index/reports 三个 sentinel 回归；普通成功路径的 e2e 全绿。 |
| V28-07 | **已修** | `GtTable` 同时构建 4.3 与 4.4 表；`GtQuery` 增加 `vcf44`，slow path 同享。 | `vcf44_first_allele_phasing_uses_the_header_aware_entries` 覆盖 `|0/1` 与 `|.`；现有 `/0|1` 等价测试保持。 |
| V28-08 | **已修** | `(NULL,count=0)` GT：render 输出 `.`，decode 返回 missing 且 `has_gt=true`，encoder 原样写 NULL 描述符。 | `null_count_zero_gt_decodes_as_missing`；renderer 单测。 |
| V28-09 | **已修** | 3 字节 GT 快路径在 `TreatAsMissing` 且恰一半缺失时交给 `from_alleles`。 | `fast_half_calls_respect_treat_as_missing`。 |
| V28-10 | **已修** | `encode_bcf` 在解包 shared header 后立即拒绝 `n_allele == 0`。 | `encode_bcf_rejects_a_zero_allele_record_before_writing`。 |
| V28-11 | **已修** | `LineRef::{terminator,body}`；无编辑整行复制，局部编辑复制原行尾。 | `crlf_terminator_is_preserved`（默认计划与 GT mask）。 |
| V28-12 | **已修** | plain/BGZF header scan 的上限改为内部可注入参数，生产仍传 `MAX_HEADER_BYTES`。 | 两个 cap 测试现在真正断言 `header exceeds 1048576`。 |
| V28-13 | **已修** | `InfoIndex::build` 与 `validate_info` 共用 `admit_tag`；新增 `classify_line_with_index`，分类器构建的索引由 head/表达式/tail 的 AC/AN 直接复用，每记录只走一次 INFO；`FilterPlan::reads_info` 决定是否构建。 | `a_prebuilt_index_gives_the_same_verdicts`（预建索引与自建索引的 verdict/funnel 完全一致）；`cargo test -p reed-filter`、`e2e_info`、`e2e_filter` 全绿。 |
| V28-14 | **已修** | `IndexQuery::query_ranges` + `merge_spans`；`RegionChunkSource` 按 span 队列 seek，跨 span 不混 chunk，每个 span 起点重置 `first_in_block_offset`。 | `v28_14_distant_regions_do_not_read_the_blocks_between_them`（两个远端 region 之间的块记录不出现）；`e2e_index` 全绿。 |
| V28-15 | **已修** | `PcaJob::run_bcf_records` 原生 hook：`gt_field` + `gt_raw_alleles`（无策略、policy-free），dosage 与文本路径同为 `NaN` 规则；`collect_via_engine` 在原生 BCF 输入时挂 hook。 | `m1_genotype_matrix_is_format_independent`（VCF/BCF 矩阵一致）；`cargo test -p reed` 的 PCA/距离 e2e 全绿。 |
| V28-16 | **已修** | 新增可复用 `FormatTable`；`decode_format_columns_with_table` 与 `decode_record_gt_with_format` 共用一次描述符遍历，`NativeScratch`/`PcaScratch` 按记录复用。 | 既有 table/text parity 测试与 `e2e_bcf` 全绿；`cargo clippy -D warnings` 通过。 |
| V28-17 | **已修** | `LineRef::sample_spans` 一趟完成样本列数校验与边界切分；`GenotypeScanner::scan` 用它替代 `validate_sample_cols` + 逐字段 memchr。 | `sample_spans_matches_the_validator_and_the_columns`（proptest 对照旧校验器与切片内容）；GT 扫描相关 e2e 全绿。 |
| V28-18 | **已修** | `BcfEncoder::new` 改为不可失败并实现 `Default`；删除 `record_ref` 中不可达的二次截断检查与过期注释。 | 编译期证明；`e2e_bcf`/`e2e_convert` 全绿。 |
| V28-19 | **已修** | plain chunk 不再复制进 worker scratch：`InflatedChunk` 直接借用 `ChunkBytes` 的 payload，生命周期同时绑定 chunk 与 scratch。 | 既有无压缩输入 e2e 全绿（chunk 语义由 `chunk/mod.rs` 单测覆盖）。 |
| V28-20 | **已修** | `NativeEncoder` 增加 `floats: Vec<f32>` scratch，Float INFO 分支复用同一缓冲，不再每条记录分配。 | `e2e_convert`/`e2e_bcf` 全绿。 |
| V28-28 | **已修（新发现，P0）** | `reed filter -Ob -W` 的 CSI 参考序列表按 `##contig` 行序建，而记录 tid 是 BCF 字典序（含显式 `IDX=`）；bcftools 因此对两个 contig 都返回 **0 条**记录（reed 自读自洽所以 e2e 未发现）。修法：`VcfHeader::contig_names_in_dict_order` 成为唯一字典序来源，filter 与 `reed index` 都按它建表。同时 `Indexer` 改为每 contig 一个 noodles indexer 后合并，解除 noodles “tid 必须递增”限制——htslib 只要求同一染色体连续，仓内 `idx.bcf` 正是 tid 递减文件，旧 `reed index` 会直接报错。 | `bcf_index_reference_table_follows_the_record_tid_space`、`index_accepts_descending_record_tids_like_htslib`（均在旧代码上失败、新代码上通过）；另用 `bcftools -r chr1/chr2` 对外验证 filter 输出与 `reed index` 输出；`contig_names_in_dict_order_honours_explicit_idx`、`contig_names_in_dict_order_fills_sparse_and_identity_cases` 覆盖 `IDX` 反序与稀疏洞。 |

## 下沉（V28-21…27，全部落地）

| ID | 状态 | 落点 | 验收 |
|---|---|---|---|
| V28-21 | **已实施** | `bcf::SharedHeader::parse` 成为 24 字节 shared 头唯一解码点（`site_header`/renderer/`encode_bcf` 共用，零等位拒绝随之统一）；`record::walk_format_fields` 成为唯一 FORMAT 描述符遍历（`format_fields` 与 `FormatTable::rebuild` 共用）。 | 既有 BCF 解析/渲染/编码测试全绿；`e2e_bcf` 17 项通过。 |
| V28-22 | **已实施** | `VcfHeader::contig_names_in_dict_order`（显式 `IDX=`、空槽保留）；filter 的 `contig_names`（记录 tid + 索引表）与 `reed index` 的 `index_names` 都改用它。 | V28-28 的两个回归。 |
| V28-23 | **已实施** | `reed_io::write::atomic::AtomicWriter`（拥有 `BufWriter`，`with_capacity`/`create_with_parents`/`finish`/`commit`）。迁移：convert、reports（含删掉“重开临时文件只为关句柄”的 `mem::replace`）、session、stats pairs、SFS、PCA、DAPC、CLI distance。 | `e2e_convert`、`e2e_reports`、`e2e_stats`、PCA/DAPC/distance e2e 全绿；截断 BCF 的三个“旧输出不变”sentinel 测试仍通过。 |
| V28-24 | **已实施** | `LineRef::{field_span,field_at,sample_span,sample_at,sample_ranges,sample_prefix_len}`；`format_layout::{subfield,subfield_span,subfield_spans}`。迁移：genotype（body+seam）、convert、reports hist、line_edit、`GenotypeScanner::split_needed`。 | `borrowed_span_accessors_match_the_validated_walk`（与 `sample_spans` 及固定列逐一对照，覆盖空区域/尾随空列/越界）。 |
| V28-25 | **已实施** | `bcf::{FrameStep, frame_chunk, record_boundary}`：跨 chunk 长度前缀与完整记录边界的唯一实现；filter 的 BCF 输入与 index job 各自保留 `Carried`、上限与解码职责。 | `e2e_bcf` 的 chunk 边界/线程无关测试全绿；`reed-engine` 24/23 全绿。 |
| V28-26 | **已实施** | `build_source(format, reader)` 统一四个格式分支（file/stdin 各一份 → 一份）；`open_file` 单次打开 + seek(0)，`sidecar_index` 唯一决定 `.csi`/`.tbi`；region 路径从 3 次打开降为 1 次（header scan 用 `try_clone`）。 | 全部 e2e（filter/index/stats/convert/region e2e）通过。 |
| V28-27 | **已实施** | `bgzf::address::{VirtualOffset, BlockTable}`；删除 `VirtualMap`/`VoMap` 的重复 walk（前者保留为薄别名），region 的 `span_end`、`ordered_output` 的 `last_block_end`/记录起点、`reed index` 的 EOF offset 都改用它。 | `bgzf::address` 三个单测（打包往返、表内行走/越界钳位、空成员保持压缩偏移连续）；region/index/write 全绿。 |

## 验收记录（2026-11-10）

- `cargo test -p reed-io -p reed-filter -p reed-engine -p reed -p reed-cli`：**41 个测试目标、732 项全部通过**（`reed-io` 179 unit，`reed` 67 unit，`reed-filter` 68，`reed-engine` 24+23）。
  - `reed` 集成：`e2e_bcf` 17、`e2e_convert` 25、`e2e_filter` 37、`e2e_index` 11、`e2e_info` 19、`e2e_reports` 4，其余 e2e 全绿。
- 新增回归：
  - `truncated_bcf_keeps_the_previous_filter_output`
  - `truncated_bcf_keeps_the_previous_index`
  - `truncated_bcf_keeps_the_previous_report`
  - `bcf_output_preserves_explicit_contig_idx`
  - `bcf_dictionary_tid_is_mapped_to_the_csi_header_order`
  - `a_prebuilt_index_gives_the_same_verdicts`（V28-13）
  - `v28_14_distant_regions_do_not_read_the_blocks_between_them`（V28-14）
  - `sample_spans_matches_the_validator_and_the_columns`（V28-17）
  - `bcf_index_reference_table_follows_the_record_tid_space`（V28-28，旧代码失败）
  - `index_accepts_descending_record_tids_like_htslib`（V28-28）
  - `contig_names_in_dict_order_honours_explicit_idx`、`contig_names_in_dict_order_fills_sparse_and_identity_cases`（V28-22）
  - `borrowed_span_accessors_match_the_validated_walk`（V28-24）
  - `bgzf::address` 三个单测（V28-27）
- `cargo clippy -p reed-io -p reed-filter -p reed-engine -p reed -p reed-cli --all-targets -- -D warnings`：通过。
- `cargo fmt --all -- --check`：通过；`git diff --check`：通过。
- 构建使用 `CARGO_TARGET_DIR=/tmp/reed-v28-target`，未污染仓库 target。
- 外部交叉验证（bcftools 1.24 / htslib 1.24）：filter 输出 BCF 的 `-r chr1`、`-r chr2` 各只返回自身记录；`reed index` 对 tid 递减的 `idx.bcf` 生成的 CSI 可被 bcftools 查询。

## 备注

- V28-02 的修法经过一次纠正：不能把 BCF 记录 tid 映射成 header 行序，否则 `idx.bcf` 的字典顺序 0(chr2)、1(chr1) 会被当成乱序。正确做法是让 BCF 的 CSI 表按字典顺序建，记录 tid 保持不变。
- V28-28 是同一类错误的第二处：V28-02 只修了 `reed index`，filter 的 `-W` 写入端仍在用行序。之前 `bcf_output_preserves_explicit_contig_idx` 的回读用 `-Oz`（文本输出丢弃 `IDX=`，行序恰好正确），因此没暴露。现在改为直接读 CSI 的参考表与 bcftools 交叉验证。
- V28-24 的一次修正：`sample_region_end` 不能复用冒号子字段 helper（样本区是 tab 分隔），保留独立的 tab 行走。
- 普通文本 VCF、无显式 `IDX` 的 BCF、LF 输入的既有输出保持不变。
