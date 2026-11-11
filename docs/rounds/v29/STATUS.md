# v29 STATUS

> 更新 2026-11-11。**8 项拍板；Batch 0/1/2 已实施，均有回归测试。**
> `V29-11…14` 是下沉设计项：V29-12 仅部分落地（`ConvertRequest::validate` 与
> 输出路径计划原语），V29-11/13/14 未实施；Batch 4 性能项未开始。

## 拍板状态

| 决策 | 关联 | 结论 | 状态 |
|---|---|---|---|
| D29-01 | W29-01（MSA 长 REF 契约） | **B：默认参考坐标对齐，重开 D15-26** | **已实施** |
| D29-02 | V29-04（HapMap 半缺失） | A：保留 `i/.` / `./i`；依据 bcftools 读入 half-call 原样保留 | **已实施** |
| D29-03 | V29-02（SFS `-s` 子集） | A：只校验/解码参与样本 | **已实施** |
| D29-04 | V29-03、V29-07（输出冲突） | A：运行前统一计划，冲突 exit 2 | **已实施** |
| D29-05 | V29-05（header-only SFS） | A：请求即写零谱 | **已实施** |
| D29-06 | W29-03（多倍体策略） | A：显式警告 + 计数，暂不支持 | **已实施** |
| D29-07 | V29-09（malformed 计数） | A：记录级一次 | **已实施** |
| D29-08 | V29-10（库校验 / `to_args`） | A：库内 validate + `to_args` target-only | **已实施** |

## 实施摘要

- **D29-01 B（MSA 参考坐标）**：记录占用 `[POS0, POS0+len(REF))`，游标推进
  `POS0+len(REF)`；块宽 `max(len(REF), len(ALTs))`，右补 `-`；插入只增加
  `W-L` 个 gap 列，不推进参考游标；同 POS 同 REF span 允许 last-wins，其余重叠
  报 `unsupported overlap`；REF 与 FASTA 大小写不敏感比较，越界/不一致/重复
  FASTA contig 报错。MSA spill frame 增加 `ref_len` 字段（内部格式）。
- **V29-01（GT 子字段 span）**：`GtSlots.raw` 记录 GT 子字段真实范围，
  binary NEXUS 直接切该范围，不再从样本列起点切。
- **V29-02（参与样本解码）**：新增 `GenotypeScanner::with_active_samples` 与
  `reed_io::bcf::decode_record_gt_masked`；convert 的 `sfs_accepts_gt` 用
  `StatsJob::active_samples()` 只校验参与样本。`-s` 等价于物理预裁剪。
- **V29-03（输出计划）**：新增 `convert/plan.rs` 的路径身份/别名检查；驱动在任何
  writer 打开前计算 `.phy/.fasta/.nexus/.bin.nexus/.used_sites.tsv/.012*/.hmp.txt/
  _sfs/*/_<contig>.fasta`，拒绝重复路径与 input/reference/populations 别名；
  hapmap→VCF 的 `-o` 先与输入比较。CLI 另外检查 `PREFIX.log` 别名。
- **V29-04**：HapMap 双字母格按两个 allele slot 处理，`AN→0/.`、`NG→1/.`、
  `NN→./.`，不再纯合化。
- **V29-05/V29-06**：请求的 SFS 即使没有 chunk 也写零谱；`merge_sfs` 合并
  `n_records`，多 chunk 计数正确。
- **V29-07/V29-08**：MSA contig 文件名不安全字节百分号编码（`a/b → a%2Fb`）；
  重复 FASTA contig、REF/FASTA 不一致、REF 越界、REF span 重叠均拒绝。
- **V29-09**：每条记录最多计一次 `malformed`。
- **V29-10**：`ConvertRequest::validate()` 覆盖空/重复 target、互斥、压缩 level
  1–9、SFS projection、MSA reference/region/collapse、VCF target 与 VCF 样本选项；
  `to_args` 文档明确为 target-only。
- **D29-06**：`ConvertCounters.polyploid` 与 `ConvertReport.polyploid` 新增；
  `--sfs` 遇到保留样本的多倍体 GT 时警告一次、计数并跳过；`--012` 保持
  vcftools 的 `-1`，同时警告一次并计数。

## 验证（2026-11-11）

- `cargo test -p reed -p reed-cli -p reed-io -p reed-stats`：**全部通过**。
- `cargo test -p reed --test e2e_convert`：**41 passed / 0 failed**。
- `cargo clippy -p reed -p reed-cli -p reed-io -p reed-stats --all-targets -- -D warnings`：
  通过。
- `cargo fmt --all -- --check`、`bash ci/check-doc-drift.sh`：通过。
- 回归均包含反向验证形态：旧实现会在对应断言处失败（`DP:GT` binary、SFS `-s`、
  半缺失、header-only SFS、多 chunk 计数、MSA 10 bp/300 bp/插入/重叠/REF 不一致、
  输出别名、malformed 多 target、polyploid 计数）。

新增测试锚点：

| ID | 测试 |
|---|---|
| V29-01 | `v29_01_binary_nexus_uses_the_gt_subfield_span` |
| V29-02 | `v29_02_sfs_sample_subset_ignores_the_excluded_genotypes`、`bcf::gt::tests::decode_record_gt_masked_skips_inactive_samples` |
| V29-03 | `v29_03_output_aliasing_the_input_is_refused_before_writing`、`v29_03_duplicate_msa_regions_are_refused_before_writing`、`v29_03_hapmap_to_vcf_refuses_output_eq_input` |
| V29-04 | `v29_04_hapmap_two_letter_half_calls_keep_a_missing_slot` |
| V29-05 | `v29_05_header_only_sfs_writes_the_zero_spectrum` |
| V29-06 | `v29_06_multi_chunk_sfs_counts_every_record` |
| V29-07 | `v29_07_msa_contig_slash_is_percent_encoded` |
| V29-08 | `v29_01_msa_reference_mismatch_bounds_and_duplicates` |
| V29-09 | `v29_09_malformed_counts_once_per_record_across_targets` |
| V29-10 | `v29_10_library_validation_and_target_only_to_args` |
| D29-01 | `v29_01_msa_reference_coordinates_for_deletions`、`v29_01_msa_insertion_adds_columns_without_advancing_the_cursor`、`v29_01_msa_overlap_rules`、`msa_long_deletion_is_not_truncated` |
| D29-06 | `d29_06_polyploid_records_are_counted_not_silent` |

## 行为/API 变化

- **MSA 输出宽度是显式契约升级**：长 REF/indel 的 FASTA 与旧版不同（300 bp 删除
  夹具 1299 → 1000）。D15-26 已在 v15 文档中标注“v29 已取代”。
- 新增公开 API：`ConvertRequest::validate`、`GenotypeScanner::with_active_samples`、
  `reed_io::bcf::decode_record_gt_masked`、`OpenedInput::path`；`ConvertCounters`
  与 `ConvertReport` 新增 `polyploid` 字段。
- `--012` 的多倍体格继续写 `-1`（vcftools 兼容），但不再是静默结果：运行会
  警告一次并计入 `polyploid`；`--sfs` 的同样记录从“静默零谱”改为“计数 + 一次性警告”。

## 逐条状态

| ID | 状态 | 证据 | 备注 |
|---|---|---|---|
| V29-01 | 已实施 | confirmed | GT 子字段 span + binary NEXUS 回归 |
| V29-02 | 已实施 | confirmed | text scanner + BCF native 都按参与样本解码 |
| V29-03 | 已实施 | code / CLI 场景 | 输出路径计划与别名检查 |
| V29-04 | 已实施 | confirmed | 双字母半缺失保留缺失 slot |
| V29-05 | 已实施 | confirmed | plain/BGZF/BCF header-only 都写零谱 |
| V29-06 | 已实施 | code | 多 chunk `n_records` 合并 |
| V29-07 | 已实施 | code | contig 文件名百分号编码 |
| V29-08 | 已实施 | code | REF/FASTA 一致性、越界、重复 contig、重叠 |
| V29-09 | 已实施 | code | 记录级 malformed 一次 |
| V29-10 | 已实施 | code | `validate` + target-only `to_args` |
| V29-11 | 未实施 | design | 共享 GT 词法/span 原语（V29-01 已直接用现有 `subfield_span`） |
| V29-12 | 部分落地 | design | `validate` 与输出路径原语已落地；`ConvertPlan` 结构未抽出 |
| V29-13 | 未实施 | design | `StatsJob::run_lines` + 统一空报告出口 |
| V29-14 | 未实施 | design | HapMap chunk/gzip、spill 清单 |
| W29-01 | 已撤回 | withdrawn | D29-01 B 已将其显式升级为参考坐标对齐 |
| W29-02 | 已撤回 | withdrawn | folded SFS 端点与 v5/05 §3.1 一致 |
| W29-03 | 已降级 | 策略空白 | 按 D29-06 A 处理：警告 + 计数，暂不支持 |

## 备注

- 当前工作树包含 v28 及更早的未提交改动；行号以最新文件为准。
- 本轮没有新增多文件事务：输出冲突在写前拒绝，逐文件仍保持原子提交。
- `V29-12` 的完整 `ConvertPlan` 与 CLI/target 解析去重、`V29-13` 的空报告出口、
  `V29-14` 的 chunk/gzip 复用仍是后续工作；Batch 4 性能/资源项需要先建基准。
