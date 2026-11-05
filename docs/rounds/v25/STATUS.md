# doc_v25 状态（2026-11-04）

> 起点：`6fd061d`（v24 收尾：945 passed / 0 failed，clippy 0 告警）。
> 组织原则见 [00-patch-sinking.md](00-patch-sinking.md)：每项动手前填
> "retired patch" 与 "sink target"，完成后按 00 §6 四问验收。
> 判定标准见 [README.md](README.md)；缺陷清单见 [01-review.md](01-review.md) §2；
> 各域规范见 02–06；顺序与验收门见 [07-roadmap.md](07-roadmap.md)；
> 大文件结构因果与拆分顺序见 [08-large-files.md](08-large-files.md)。

## 状态总表

| 编号 | 来源 | 状态 | 关键证据 / 规范 |
|---|---|---|---|
| V25-01 | F1 | **已完成** | 四条路径：scanner `parse_allele`（`n_alt` 上界）/ BCF `decode_pair`（解出 `(v>>1)-1` 后对照 `n_allele`）+ `decode_slow` / encoder text+BCF / convert `parse_gts` / PCA `dosage`；failing-first 证据见下 |
| V25-02 | F2 | **已完成** | GT 位置 **6→1**：`bcf/gt_decode.rs` 纯搬家（`3c377b9`）+ `format_layout::{key_index, gt_index, gt_subfield}` 单点；encoder/filter(3 处)/convert/PCA(body+seam) 全部收口；`grep gt_subfield_index` = 0 |
| V25-03 | F3 | **已完成** | `reed_io::vcf::line::{RecordClass, MalformedLine, classify_line, body_lines}` 单点；convert row/SFS、stats(popgen+reports)、rename、filter 五入口策略表收口；`line_is_sane` 删除；seam 空行 e2e（LF/CRLF） |
| V25-04 | F4 + R3 + R2 | **已完成** | `semantics::{compare_missing, info_af}` 单点；`Expr::Info` 进路由；`eval_head` 取 `plan.info_trusted`；R2 闭区间回放；三策略/子集 e2e |
| V25-05 | F5 | **已完成** | `reed_stats::merge_zip_eq`（release 断言 + ctx）；23 处 zip 收口；`ThetaAcc` resize 删除 |
| V25-06 | F6 + R5/R6/R7 | **已完成** | `Threads`/`MemMiB` + `ThreadBudget::from_cli`（7→1）；`params.rs` 统一 `-O`；`--ld-decay` value_parser；`ci/check-structure.sh` |
| V25-07 | F7 + R4 | **已完成** | `GtCode::class` + `SiteCounts::from_codes` 单点；`KeepMask`/`restrict` 可失败；幽灵 counter 删除；`PolicyCounters::FIELDS` 单表 |
| V25-08 | F8 + R15 + R11/R14 | **已完成** | facade 回填/错误分类/R11/R14（`b9353c7`）+ StageMetrics 渲染单点/五命令键集 e2e + provenance 提升 core/argv 引号（`c5c7829`） |
| V25-09 | F9 | **已完成** | `spill.rs` TempSpill（O_EXCL + Drop）；stdin/convert 收口；攻击测试随迁（`db5fff5`） |
| V25-10 | F10 + R1 + R12 | **已完成** | 残差混合尺度 + report 字段/告警；`fix_column_signs` 共用；`STREAM_DEBUG` 删除；`linalg.rs` 删除（`fd591f2`） |

## V25-01 实施记录（2026-11-04）

- **failing-first**：实现前 7 个新断言按预期失败——scanner 的 `0/9`/`1|2`/
  单倍体 `9`、BCF `decode_record_gt`、encoder text/BCF、PCA `run_chunk`/
  `run_seam`（`dosage` 把越界折成 NaN，位点被静默跳过）。
- **域规则**：合法 allele 是 `0..=n_alt`（`ALT=.` → 0）；BCF 侧存储值是
  `(allele+1)<<1|phased`，先解 `allele=(v>>1)-1` 再对照 `n_allele`；`ALT=.`
  时文本快路径的 allele 1 也拒绝。
- **首版实现的一个真 bug（被测试抓住）**：`decode_pair` 最初把编码字节当
  allele 索引，`0/1` 的首字节 2 被误报 "allele 2"——`table_matches_the_text_path`
  立即失败；改为解码后再判域，124 个 reed-io 单测全绿。
- **旧 fixture 修复（保留原意图）**：`scan_edge.vcf` 300/500 补足 ALT；
  scanner 的 `10/11` 测试补 12 个 ALT；fuzz reference 建模 `n_alt`；
  determinism 的 `ALT=.` 记录只出 `0/0`/`./.`；relaxed-policy e2e 打开
  `MultiallelicPolicy::Allow`（否则 ALT=G,T 记录被双等位门挡掉）。
- **验收**：`ci/test-parallel.sh` **955 passed / 0 failed**（945 基线 + 10
  新测试）；`cargo +stable clippy` 0；`ci/check-doc-drift.sh` 97 ok。
- **边界**：`convert` 保持 B13-9 政策（越界 GT 计 malformed 并跳过，exit 0），
  不新增第二套错误政策；`-Ob` 部分文件残留属 R13（v26 原子性）。

## V25-02 实施记录（2026-11-04）

- **3c377b9（纯搬家，行为不变）**：`gt_shape` / `strip_gt_prefix` /
  `gt_codes_into` / `push_gt_code` / `parse_i64` 从 `encode_native.rs` 移到
  `bcf/gt_decode.rs`，两个测试随迁；`encode_native` 只留 typed 编码器。
- **共享 API**：`format_layout::{key_index, gt_index, gt_subfield}`；`Layout`
  的 intern 也走 `key_index`。`gt_subfield` 对短列返回 `None`（GT 缺失），
  绝不回退到首个子字段。
- **6 个消费点收口**：scanner（本就走 `Layout.gt_index`）、encoder
  （`keys.position(b"GT")` → `gt_index`）、filter（`gt_subfield_index` 删除，
  `emit/single_pass/filtered` 三处 → `gt_index`）、convert（`memchr(b':')`
  → `gt_index` + `gt_subfield`）、PCA `run_chunk`/`run_seam`（两份首字段
  启发式 → 同一个 `dosage_from_column` 助手，顺带消掉重复的域错误映射）。
- **测试（5 个新断言）**：scanner `DP:GT`/末位 GT；encoder `DP:GT` 与
  noodles 逐字节一致；convert 012 `DP:GT` == `GT:DP`；filter 文本与 `-Ob`
  都解出 0/1；pca `DP:GT` 与 `GT:DP` 的 `.eigenvec` 相同。
- **验收 grep**：`gt_subfield_index` 0；`position(|k| *k == b"GT")` 0；
  `genotype.rs` 的 `position(|&c| c == b':')` 0；`parse_gts` 的
  `memchr(b':')` 0。
- **全量**：**960 passed / 0 failed**（955 + 5）；clippy 0；doc-drift 97。

## V25-03 实施记录（2026-11-04）

- **failing-first**：实现前 2 个新 e2e 按预期失败，错误都是
  `Error(Engine(Job("record line has fewer than 8 columns")))`——
  `convert --sfs` 的 seam 空行与 CRLF body 空行；seam fixture 用压缩块
  精确构造 1 MiB chunk 边界（不是靠 `--memory` 碰运气，原文修正见 02 §3）。
- **分类单点**：`reed_io::vcf::line` 新增 `RecordClass` / `MalformedLine` /
  `classify_line` / `body_lines`；`MalformedLine` 携带失败项，消费点保留
  原诊断（"fewer than 8 columns" 等）；POS 0 归 `Position`。
- **五入口**：
  - convert row：`line_is_sane` 删除；blank 跳过（不计、不占 draw slot），
    malformed 计数继续（原位保留 vcf2phylip 的 `snp_num` 口径）；
  - convert SFS：纯 SFS 也 sanitize（干净 body 零拷贝返回）；seam 只把
    `Record` 交给 stats —— 修掉了"mixed 跳过、pure 硬错"的第二处分叉；
  - stats（`StatsJob` + `ReportJob`）：blank 跳过，malformed 仍致命；
  - rename：blank 在 body 与 seam 都跳过（原先 body 丢 LF、seam 多出空行、
    CRLF 原样写出）；非空 malformed 仍逐字节透传；
  - filter（`RecordPipeline::process`/`process_scanned`）：blank/malformed
    都报错；body 行循环改用 `body_lines`，只丢终止符后的空片。
- **测试（9 个新断言）**：分类器 4 个单测（terminator 形式、失败项、
  `body_lines`、label）、stats seam/body 单测、convert seam LF/CRLF、
  convert CRLF body、filter blank 拒绝、rename body+seam。
- **结构轨道②**：`drive_pass1` / `output_sink` 合一两条 pass1 收尾段；
  `reed_io::bcf::bcf_file_header` 下沉 BCF magic；`filter` 无
  `BCF\x02\x02` 字面量。纯搬家，输出字节不变。
- **提交**：`0cef99f`（语义）、`9b61a37`（结构）；**970 passed / 0 failed**；
  clippy 0；doc-drift 97。

## V25-04 实施记录（2026-11-04）

- **路由（F4）**：`references_genotypes` 补 `Expr::Info`；`expr_reads_gt_fields`
  → `Expr::reads_gt_fields`；`INFO/AC` 与裸 `AC` 同路由（`[site]` 编译期拒绝，
  `[site_after_mask]` 接受）；`eval_head` 的硬编码 `true` 换成
  `plan.info_trusted`，注释删除。
- **semantics 模块**：`compare_missing`（6 op × 3 ctx 真值表）收口
  `any_pair`/`missing_value_passes`/列式 mask 三处；`info_af` 收口 head/tail
  两份 AC/AN 读取与畸形防护。畸形 `AC=200;AN=100` 现在 head/tail 同为
  "无可用 INFO"，tail 回退 genotype counts（有意行为变化，单测钉住）。
- **R2**：`when` 回放改用解析后的 `Cmp`；`DP <= 3` → `--min-dp 4`、
  `DP >= 20` → `--max-dp 19`；不可精确表达的形态省略（`to_args` 不可失败，
  与规范"usage 报错"有偏差，已记录）。DP=3 边界 render→parse mask 集逐值相同。
- **R3**：缺失 DP 在表达式语境 mask、在 `mask_if` 默认 `Keep` 不 mask——
  差异是契约，e2e 钉住。
- **R10**：分叉**不成立**（两策略共用同一 `mask_rule`），按 B 落地为等价性
  测试（重叠 field/expr 规则的 funnel/masked_by_rule 相同）；归因顺序确定为
  首个执行序命中。
- **三策略**：INFO built-ins 在 SinglePass/TwoPassReparse 逐字节一致；六个缺失
  比较在三策略逐字节一致；sketch 对"读记录"的 tail 表达式仍按既有契约拒绝
  （mid stage 无记录字节），测试断言拒绝而非换源算错。
- **验收**：`ci/test-parallel.sh` **984 passed / 0 failed**；clippy 0；
  doc-drift 97（N12 守卫改为 `fn info_af(info` 单点）；提交 `00e3a23`。

## V25-05 实施记录（2026-11-04）

- `merge_zip_eq` 落在 `reed-stats/src/merge.rs`，`assert_eq!` 而非
  `debug_assert!`；单测含 `catch_unwind` 消息断言。
- 23 处 component-wise merge（规范 15 + report 6 + histogram 1 + window 嵌套）
  全部改走；`theta.rs` 的 resize-then-zip 删除（形状一致成为合并前提）。
- `grep 'iter_mut().zip' crates/reed-stats/src/` 仅剩 `merge.rs` 实现一处。
- 985 passed / 0 failed；clippy 0；提交 `de05e12`。

## V25-06 实施记录（2026-11-04）

- `Threads`/`MemMiB`（reed-engine::budget）带 clap value_parser；8 个 CLI 字段
  改类型；`ThreadBudget::from_cli` 取代 7 处构造；CLI 的
  `available_parallelism().expect` 清零。
- `reed-cli/src/params.rs`：`-O` 的 letter/level 单点（v/z/b/u，1-9），
  `split_output_type` 与 R4 的三份拷贝删除；rename/hapmap 用文本子集。
- `--ld-decay` 由 value_parser 产出 `Option<u64>`（裸=100 kb，0 拒绝），translate
  二次解析与死分支删除。
- usage 矩阵 17 例全 exit 2；`ci/check-structure.sh` 四条门全绿。
- 偏差：distance 自建 pool 未换 `run_auxiliary`（【荐】项，闭包边界需重排，
  记 v26）。
- 988 passed / 0 failed；clippy 0；doc-drift 97；提交 `3c2deee`。

## V25-07 实施记录（2026-11-04）

- **计数 5→1**：`GtCode::class()`/`n_ref`/`n_alt` 单点；`from_codes(_kept)`
  收口 filter 计数、`add_code`、reports per-group、stats；`count_site` 死路径与
  短 mask 默认删除；kernel proptest 锁等价。
- **mask**：`restrict -> Result`；`KeepMask::try_new` 单点长度校验（四个驱动）；
  `pruned_header`/`from_assignments`/计数核心不再默认短 mask。
- **counter**：幽灵 `n_haploid_as_missing` 删除；`PolicyCounters::FIELDS` 唯一
  表驱动 io bit packing 与 filter 报告；reports 手抄 policy 改调用。
- 991 passed / 0 failed；clippy 0；structure 7/7；doc-drift 97；提交 `11b9b18`。

## V25-08/09/10 实施记录（2026-11-04）

- **V25-08a**：facade re-export + CLI 三依赖删除；`Error::Foreign`
  （runtime）；graph/LDA 按 variant 分流；rename 输出创建 runtime；R11 校验
  前移；R14 piHat 两列写样本名。e2e：退出码矩阵 + R11 三组合 + piHat。
- **V25-08b**：`StageMetrics::render_tsv` 唯一渲染器；五命令键集 e2e。
- **V25-08c**：`Provenance` 提升 `reed-core::provenance`；argv shell 引号
  单点（`params::command_line`）。
- **V25-09**：`spill.rs::TempSpill`（`create_new`、重试、Drop unlink）；stdin
  与 convert 共用；SEC-1 攻击测试随迁（含 convert tag）。
- **V25-10**：Lanczos 残差混合尺度 + `stream_max_residual` 报告/告警；
  `fix_column_signs` 共用（R1）；`STREAM_DEBUG` 删除；R12 死文件删除。
- **测试**：998 passed / 0 failed；clippy 0；structure 11 条门全绿；
  doc-drift 97。

## 本轮增补（2026-11-04 复核后）

| 增补 | 落点 | 状态 |
|---|---|---|
| PCA `PcaJob::run_chunk`/`run_seam` 两份 GT 位置并入 6→1 | 02 §2（V25-02） | 已写入规范 |
| PCA `dosage` 越界 allele 从 NaN 改为 `MalformedRecord` | 02 §1（V25-01） | 已写入规范 |
| R2（`when` 的 `<=`/`>=` 回放）并入 V25-04 | 03 §4 | 已写入规范 |
| R11（stats 退出码）/ R14（piHat FID/IID）并入 V25-08 | 06 §1 | 已写入规范 |
| R12（`reed-reduce/src/linalg.rs`）删除并入 V25-10 | 04 §5 | 已写入规范 |
| counter 三处镜像收成 core 单一定义 | 04 §2/§4（V25-07） | 已写入规范 |
| `ci/check-structure.sh` 升为【规】 | 07 §3，步 11 | 已写入规范 |

## 结构轨道（18 文件审查，08）

| 拆分 | 绑定 | 状态 |
|---|---|---|
| `encode_native.rs` → `bcf/gt_decode.rs` + `bcf/typed.rs` | V25-01/02 | **已完成**（`3c377b9`；`typed.rs` 未拆——`encode_native` 现在只剩 typed 编码器） |
| `filter/mod.rs` pass1 合一 + BCF magic 下沉 | V25-02/03 | **已完成**（`9b61a37`：`drive_pass1`/`output_sink`/`bcf_file_header`） |
| `convert/job.rs` spill → `convert/spill.rs` | V25-09 | **已完成**（`0373ed5`，纯搬家） |
| `reports.rs` 计数/policy 收口 | V25-07 | **已完成**（`11b9b18`：`code.class()` 单点、`stats_policies()` 单点、`classify_line` 收口；`count_site` 仅剩 oracle 的独立参照） |
| `plan/mod.rs compile` 纯搬家 | V25-04 后 | **已完成**（`7d235ab`：`plan/compile/{head,mask,tail,finish}.rs`，入口只做接线） |
| `request.rs` / `matrix.rs` / `popg.rs` | v26 | 不阻塞主线 |

- 判定：18 组判定（20 个 ≥670 行文件），9 组乱 / 1 组中度 / 8 组不乱；
  8/9 乱组与补丁热点重合。
- 因果：混杂 → 语义多入口 → 补丁只能打调用点 → bug；拆分只作为下沉载体，
  先纯搬家再改语义（08 §2–§3）。

## 待拍板

- ~~幽灵计数器 A（删，荐）/ B（补 policy 变体）~~（**已按 A 落地**：`n_haploid_as_missing`/`TRIG_HAPLOID_MISSING` 在 `11b9b18` 删除，常驻门已锁；04 §4）；
- ~~R10 mask 归因顺序~~（已决：分叉不成立，按 B 落地等价性测试，03 §5）；
- 07 §2 的 v26 候选九组（含 CI-1 钉子户，建议最先）；
- R8 popg 裸路由是否保留。

## 验收基线（2026-11-04，起点）

| 检查 | 结果 |
|---|---|
| `ci/test-parallel.sh` | **1005 passed / 0 failed**（V25-01…V25-10 累计 +60，945 基线） |
| clippy / doc-drift | 0 告警 / 97 项 ok（V25-04 后，N12 守卫指向 `info_af`） |

## 2026-11-05 全库 review 复核修复（三代理）

范围 `6fd061d..7d235ab`（17 commit），三域并行 review 后逐条用真实二进制 /
bcftools 1.24 / 新旧二进制差分复现，修 10 条：

| 修复 | 根因 | 回归测试 |
|---|---|---|
| pca 空行/seam 失败 | `run_chunk`/`run_seam` 未接 `classify_line`（V25-03 漏收口） | `genotype::tests::v25_03_*` ×2 |
| `check-structure.sh` GT 门恒过 | `! grep -q ... \| grep -v ...`：`-q` 无输出 | 注入违规手工验证（门现失败） |
| `convert --sfs` exit 0/1 分叉 | SFS 把越界/非数字/polyploid GT 交给严格 stats 扫描器 | `v25_01_convert_sfs_skips_...` |
| `-q/-Q` 对 `ΣAC > AN` 静默回退 | 无显式错误路径（bcftools 实测 abort，exit 1） | `ac_sum_above_an_aborts_the_run` |
| hapmap→VCF 输出创建失败 exit 2 | `File::create` 映射到 `Error::usage` | `output_create_failure_is_runtime_...`（扩展 hapmap 分支） |
| `--xval` 的 LDA 数据失败 exit 2 | 整个 `XvalError` 统一 `usage` | `xval_lda_failures_keep_the_lda_classification` |
| `convert --sfs` 全站跳过时 panic | SFS 累积器惰性初始化 + `expect`；V25-01 严格 GT 拒绝扩大了触发面 | `convert_sfs_writes_an_empty_spectrum_when_no_site_is_projectable` |
| `-q/-Q` 在 `-i/-e` 之前判 AF | 规则顺序与 bcftools 相反，`ΣAC>AN` 硬错误放大了影响 | `an_expression_excluded_record_does_not_reach_the_af_abort` |
| 前导零等位基因被全目标丢弃 | SFS 扫描器的 5 位上限泄漏进 row targets | `zero_padded_allele_stays_in_the_row_targets_but_not_the_sfs` |
| `-s` 子集下 `INFO/AC|AF` 丢记录 | V25-04 把 `INFO/` 与裸名路由合一；bcftools 在 subset 前求值 | `info_spelling_reads_the_pre_subset_info_like_bcftools`、`three_strategies_...` |

- AC 长度不匹配保持 D12-13 的基因型回退；`semantics.rs`/`head.rs` 注释改为
  明写与 bcftools `-q` 的分歧（bcftools 丢记录，reed 按 D12-13 回退）。
- `INFO/` 前缀表达式恢复为纯 INFO 读取（head 阶段），裸名走基因型路径（tail）：
  bcftools 在 `bcf_subset` 之前跑 `filter_test`，两者本就不同队列。v25 曾把
  两种拼写路由合一，是回归；`03-filter-semantics.md` 的验收门现已真正满足。
- 验收：`ci/test-parallel.sh` **1013 passed / 0 failed**（1005 + 8）；
  clippy 0；`check-structure.sh` 16 项 ok；`check-doc-drift.sh` 97 项 ok。
