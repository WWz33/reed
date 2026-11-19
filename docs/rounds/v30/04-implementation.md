# v30 实施记录（P0–P2 第一批）

> 2026-11-12。依据 [03-synthesis.md](03-synthesis.md) 的 P0→P1→P2 顺序落地。
> 本文件只记录**已改动的代码、验证命令与结果**；未完成的条目留在 §4，不勾销。
> 其中 CI/发布描述是 2026-11-12 的历史记录；当前 main/dev 边界和剩余问题以 [06-final-audit.md](06-final-audit.md) 为准。
> 审计结论仍以 [01-review.md](01-review.md) 与 [04-recheck.md](04-recheck.md) 为准。

## 1. P0 契约裁决

| 产出 | 内容 |
|---|---|
| `docs/spec/api-contract.md` | 发布面/稳定面/外部依赖/MSRV-feature/兼容承诺 5 问全部裁决，冻结至 2027-05-12 |
| 发布面 | 只发布 `reed` 门面；`reed-distance/graph/reduce` 从默认发布改为 `publish = false`（与其余 8 个内部 crate 一致） |
| `ci/api/reed.txt` | `cargo public-api -p reed -s` 快照，1976 行（`prelude` 的逐项再导出会展开；其他模块再导出仍不展开，是下界） |
| `ci/check-api.sh` | 重新生成并与快照 diff；漂移即 exit 1，要求人工评审 |
| `.github/workflows/ci.yml` | 删除"main 只 build"的分叉：main/PR 都跑 fmt、release build、debug/release test、doctest、clippy、`--no-default-features`、doc-drift、structure、API 快照；MSRV 1.91.1 对全部分支 |
| P3 接入 | `ci/check-semver.sh`（对最新可达 tag 跑 `cargo semver-checks -p reed`，CI 必过）、`ci/check-package.sh`（`reed` 为发布目标 + 每个 path 依赖带 version + `cargo package --list`）；见 ADR §6 |

验证（本机）：

```
bash ci/check-api.sh                  # exit 0
bash ci/check-doc-drift.sh            # exit 0
bash ci/check-structure.sh            # exit 0
cargo check -p reed-reduce --no-default-features   # exit 0
```

## 2. P1 语义内核

### V30-R01 GT 半缺失判定

- `reed-core::genotype::has_called_allele(&[u8])`：文本级唯一规则——只有 `.` 是缺失，
  `./1` 与 `0/.` 都算 present，`./.` 不算；畸形 token 仍算 present，让它走 malformed
  路径而不是被误当缺失。单测 `called_allele_rule_handles_half_calls` 覆盖 11 种拼写。
- `convert/job.rs` 的 `gt.present` 由 `!raw.starts_with(b".")` 改为该函数，删除
  vcf2phylip 的前缀规则；`--phylip/--012/--msa` 与 `stats --missing` 的半缺失口径一致。
- 未新增重复谓词：`GtCode::is_called_any` 仍是解码侧的同一规则。

### V30-R04 `--msa --region` 的 contig 词表

- `ConvertTarget::Msa.regions` 从 `Vec<MsaRegion>` 改为 `Vec<String>`（原始拼写）；
  `MsaRegion` 类型删除（未发布 API，仅测试使用）。
- `build_msa` 用参考 FASTA 的 contig 名调 `Region::parse_in`；`chr1:100` 作为字面
  contig 不再被拆成 `chr1` + 位置 100。0 边界与 start>end 的校验一并移到这里。
- 回归用例 `v30_r04_msa_region_uses_reference_contig_names`：参考 `>chr1:100`，
  `--region chr1:100` 产出 `out_chr1:100_1-10.fasta`。

### V30-R02 AC/AN 单一读取器

- 新 `reed-core::count::InfoCounts::from_info(ac, an, n_alt)`（含
  `InfoCountsOutcome`）：长度域、AN>0、checked 求和、`ΣAC>AN`、溢出五类结果；
  解码式 REF = `AN − ΣAC`。
- `reed-filter::semantics` 的 `InfoAfCounts/InfoAfOutcome` 改为 core 类型的 `pub use`
  别名；`info_af_indexed` 只负责从 `InfoIndex` 取两个值并保留 warn-once。
- `stats --summary` 删除手写 INFO/AC+AN 块（原 `ac.iter().sum()` 会 debug panic /
  release 回绕），改走同一函数；AC/AN 的 INFO 查找也改用 `InfoIndex`（首匹配、
  重复 tag 拒绝），与 filter 同一约定。
- 测试：core 5 个（含 `u64::MAX,1` 溢出）；`summary_uses_the_shared_ac_an_reader`
  比较"溢出 / 无 INFO"渲染逐字节相同，且"合法 INFO"与 GT 回退落在不同 AF 桶。

### V30-R03 BCF `pos0` 哨兵（本批补齐）

- `reed-io::bcf::record::pos1_from_pos0(i32) -> Result<Pos1, IoError>` 是唯一读取规则：
  `-1` = htslib 的 telomere 哨兵 → `POS 0`（`vcf.c:1649-1650`），`< -1` 报
  `MalformedRecord`；`0 → POS 1` 等普通值走 `Pos0::to_pos1`。
- 逆方向 `pos0_from_pos1(u64) -> Option<i32>` 拥有同一哨兵：`POS 0 → pos0 = -1`；
  text→BCF 以前把它写成 `0`（`POS 1`），现在保留 telomere。
- `SharedHeader::pos1()` 与 `SiteHeader::pos1()` 都委托该函数。
- 原有裸转换/分支全部收敛：`render_site_part`、`render_sample_part`、`stats/run.rs`、
  `filter/record_pipeline.rs`、`genotype.rs`、`index.rs` 两处、`encode_native` 的
  BCF→BCF 校验与 text→BCF 编码。
- 行为变化：`pos0 = -1` 在 filter/genotype/index 以前静默变成 `POS 1`、在 stats 以前
  报 `BCF POS is negative`；现在所有消费者都输出 `POS 0`，`< -1` 统一报
  "below the telomere sentinel (-1)"。
- 测试：`telomere_sentinel_is_the_only_pos0_conversion`、
  `pos0_from_pos1_is_the_inverse_with_the_same_sentinel`、
  `shared_header_pos1_applies_the_sentinel`；`e2e_bcf` 17 项、`e2e_index` 11 项全绿。
- 收敛门：`rg -n "try_from\((header\.)?pos0\)" crates --glob '*.rs'` 只剩定义体内 1 处。

### V30-R12 `reed-filter` 去 `reed-io`（语义类型下沉 L0）

- 用 `git mv` 把四组记录文本语义从 `reed-io/src/vcf/` 搬到
  `reed-core/src/vcf/`：`line`（`LineRef`/`RecordClass`/`classify_line`）、
  `info`（`InfoIndex`/`info_fragments`/`validate_info`）、
  `format_layout`（`IntSpec`/`FormatLayoutCache`）、
  `line_edit`（`EditPlan`/`InfoValueOut`/`edit_vcf_line`）；文件内只把
  `reed_core::` 改成 `crate::`。
- 新增 `reed-core::vcf::fixed`：`QualValue`+`parse_qual`（纯语法）与
  `parse_pos1_number`（POS 语法）。`reed-io::vcf::fixed` 的
  `parse_pos1`/`parse_qual` 保留为“core 语法 + 带 locus 的 `IoError`”包装，
  错误文本不变；`record.rs` 的 `classify_line_with_index` 改用
  `parse_pos1_number`，POS 校验不再经过 IO 错误类型。
- `count_alts` 移到 `reed-core::site`（`N_ALT` 与 filter 规则共用）；
  `reed-io::vcf::site` 只保留传输侧的 `ContigTable`/`SiteView`/
  `indexed_span_len` 并 `pub use` 该函数。
- `reed-io/src/vcf/mod.rs` 变为 `pub use reed_core::vcf::{...}` + 传输模块
  （`fixed`/`site`/`genotype_scan`），旧路径 `reed_io::vcf::line::...` 仍可用。
- `reed-filter` 的三处 `reed_io` 引用清零：`expr/mod.rs` 与 `eval/head.rs`
  的 QUAL 改读 `QualValue`（`Missing|Malformed` 都不满足 `MinQual`，与旧
  `.ok().flatten()` 一致）；删除 `From<reed_io::IoError> for FilterError`；
  `Cargo.toml` 删除 `reed-io.workspace = true`，描述改为 “policy layer over
  reed-core”。
- 门面侧 `PendingOrderedOutput::commit()` 的 `IoError` 改走 `reed/src/filter`
  已有的 `io_err`；`std::io::Error` 仍走 `From<std::io::Error>`。
- 结构门新增两条：`reed-filter/Cargo.toml` 与源码不得出现 `reed-io`；
  `LineRef`/`InfoIndex` 必须住在 `reed-core` 且 `reed-io` 旧文件已删。
- 验证：`cargo tree -p reed-filter -e normal` 输出不含 `reed-io`；
  `pos_grammar_is_shared`/`qual_grammar_is_shared` 单测通过。

### V30-R05–R07 输入能力成类型

- `crates/reed/src/input.rs` 新增四个类型：`InputMeta { format, framing,
  bcf_dict }`（`Framing::{Bgzf,Gzip,Plain}`）、`StreamedInput`（一次性流，
  `OpenedInput` 保留为它的公开别名）、`ReplayableInput`（文件或 stdin
  spill；`from_stdin(tmp_dir)` 在**首次打开前**落盘，`replay(regions)` 每次
  产出一条新流）、`SeekableInput`（已校验 BGZF/BCF 与 sidecar 的 replay 能力），
  以及 `InputCapability::{OneShot,Replayable,Seekable}`。
- `StreamedInput::into_parts()` 是取 source 的唯一出口；删除 `source_mut()`。
  `meta()`/`is_bcf()`/`bcf_dict()` 取代 `.bcf.is_some()` 的运行期推断，
  `InputJob` 构造时把 `RecordFraming` 记成字段（`self.bcf.is_some()` 归零）。
- `PendingInput`：一次 `open` 同时服务 header 探测与 source 构造——文件句柄
  用 `try_clone` 扫 header 后 `seek(0)`，region 分支复用同一句柄；
  filter 的 `-r` 规格解析与 CLI `prelude` 都改走它，删掉“先探测、再按路径
  重开”的两段式。
- `open_regions` 在构造 `RegionChunkSource` 前校验 framing 与 sidecar：
  plain VCF/普通 gzip 报 `IoError::Capability { needed: Bgzf }`，缺 sidecar
  报 `Needed::Index`，stdin+regions 报 `Needed::Seek`，一次性输入要求
  replay 时报 `Needed::Replay`。
- 退出码归类：`Error::kind()` 把 `IoError::Capability` 与 `FilterError::Usage`
  归为 Usage（exit 2）；filter 的 `io_err` 保留该分类，不再统一压成
  `FilterError::Io`。
- filter 策略选择改为询问能力：`open_pass1` 接收 `&InputCapability`；显式
  `--strategy two_pass_reparse`（或任何需要第二趟的 shape）遇到 OneShot 输入
  直接报 usage，不再“打印 two_pass_reparse 却只跑一趟”；BCF 的 sketch 例外
  只在确实选中 sketch 时改走 two-pass。
- 验证：`rg "bcf\.is_some\(\)"`/`rg "source_mut\(\)"` 全库 0 次；
  `strace -e openat` 下 `filter in.vcf.gz -r ... --strategy single_pass`
  对输入文件只有 1 次 `openat`；`cli_contract::input_capability_gaps_are_usage_errors`
  覆盖 plain+gzip、缺索引、stdin+`-r`、stdin+two-pass 四种 exit 2。

### V30-R09 预算单一所有者

- 新 `crates/reed/src/context.rs`：`RunContext { budget: &ThreadBudget,
  abort: &AbortSignal }` 是门面唯一的运行上下文。八个 driver（filter、popgen、
  stats、convert、rename、pca、dapc、index）都接收 `ctx: &RunContext<'_>`，
  不再各自构造预算或 abort 通道；`run_budget.rs` 与 `run_budget::for_run` 删除。
- 预算所有权：调用方构造 `ThreadBudget`（CLI 用 `ThreadBudget::from_cli`，
  测试显式构造）并借给 context；`FilterOptions` 的 `threads`/`mem_bytes`
  第二份预算已删除，策略选择、pass0、prune 池与 pipeline 都读 `ctx.budget`。
- 中止所有权：调用方持有 `Abort`（`abort_for(&budget)`），只把
  `&AbortSignal` 放进 context。信号仍携带触发句柄（`trigger_handle`），
  所以引擎内部的 `FirstError` 能唤醒所有等待者。`AbortSignal::idle()` 给
  没有取消源的库调用/测试使用；pipeline 检测到无句柄信号时会为本次运行
  换入一条本地通道，首错仍能唤醒全部 stage（回归测试
  `idle_signal_still_wakes_waiters_on_a_sink_error`）。
- convert：`ConvertSink` 的 `spill_budget = max(8 MiB, mem/4)` 来自 `ctx.budget`；
  `enforce_spill_budget` 用它替代固定常量，8 MiB 成为下限而非上限。
- LD：`ContigMatrix::size_bytes()` 给出驻留字节估计；`LdExecutor` 记录
  未合并矩阵的字节数，`finalize_ld` 在“数量 cap 或字节 cap”任一触发时先 join；
  单个大矩阵仍可入队（否则无法处理），但不会叠加。
- 验证：`ci/check-structure.sh` 的三条新门（RunContext 存在、八个 driver
  都接收它、`FilterOptions` 无第二预算）；`mind_strategy_reads_the_run_budget`
  用 12k 位点 + 200 条 mask 规则，经 `RunContext` 传入 16 MiB 预算时选
  `two_pass_reparse`、64 MiB 时选 `sketch_then_emit`。
- 顺带删除 `PcaRequest`/`DapcRequest` 的 `threads`/`memory` 死字段（Q07）：
  driver 只读 context，请求不再携带一份被忽略的副本；CLI 翻译层同步删除赋值。

## 3. P2 资源所有权

### V30-R11 `OrderedOutput::adopting`

- `adopting` 现在 `assert!(self.atomic.is_none())`：`create(p).adopting(a)` 会丢弃
  `p` 的临时文件，属于 API 误用；回归测试 `adopting_a_second_file_is_rejected`。
- 两处正常调用（filter 输出、hapmap→VCF）都是 `new(...).adopting(...)`，不受影响。

### V30-R10 `TempSpill` 持有 inode 并共享存活表

- `TempSpill` 增加 `file: File` 字段，`create` 时 `try_clone` 出所有权句柄；
  新增 `try_clone()`。convert 的追加写从"按路径 reopen + append"改为
  "clone 已持有 inode + seek end"，不再信任名字。
- **共享存活表**（本批补齐）：`reed-io::write::atomic` 暴露
  `register_live_temp` / `unregister_live_temp` / `live_temp_paths`
  （`#[doc(hidden)]`，仍是同一张 `LIVE_TEMPORARIES`）。`TempSpill::create`
  创建后立即注册，`Drop` 先注销再 unlink；`try_clone` 失败时手动注销并删除。
  因此输出临时文件与 spill 的任一次扫除都能看到对方仍在写的文件。
- **启动扫除**（本批补齐）：`spill::sweep_stale_spills` 在每个 spill 目录
  **首次创建前**运行一次（`SWEPT_DIRS` 记录已扫目录），只删除同时满足
  “mtime ≥ 60 s”且“pid 是本进程（未注册）或 `/proc/<pid>` 已不存在”的
  `reed-<tag>-<pid>-<counter>.<ext>`；另一个仍存活进程的 spill 只告警不删除
  （pid 复用无法从名字证明）。名字锚定完整形状，`reed-notes.txt`、缺 counter、
  非数字 pid 等用户文件永不触碰；`remove_file` 仍是 unlink，不跟随符号链接。
- 测试：`a_spill_is_registered_in_the_shared_live_table`、
  `a_live_spill_survives_the_sweep`、`an_unregistered_stale_spill_is_swept`、
  `the_startup_sweep_runs_once_per_directory`、`a_dead_owners_spill_is_swept`、
  `a_live_owners_spill_is_kept`、`spill_like_user_files_are_never_swept`。

### V30-R08 prune 线程 owner

- 新 `PrunePool`（`filter/pass0/prune.rs`）：`spawn` 在 cap 处 join 最旧；
  `finish` 按顺序 join 并返回全部 outcome；`Drop` 先置 cancel 再 join 所有句柄。
- `prune_one` 新增 `Option<&AtomicBool>`，每 1024 个位点轮询一次；取消时返回空
  outcome（该路径只在错误退出时触发，结果被丢弃）。
- `MatrixSink` 的 `prune_jobs/prune_done` 合并为一个 `PrunePool`。
- 测试：`the_pool_joins_at_the_cap_and_on_finish`、`a_cancelled_prune_returns_without_scanning`。
  worker panic 的 join 错误现在带 panic message，不再只有 "panicked"。

### V30-R13 `--fill-tags` 静态表

- `reed-core/src/count/recompute.rs`：`RECOMPUTE_TAGS: [RecomputeTag; 5]` 一行一个标签，
  字段为 `name/ty/number/description/compute`（AC/AN/NS/MAF/HWE）；`find`/`canonical`
  做大小写不敏感解析（`strcasecmp` 语义）。每个标签的 `compute` 直接从 `SiteCounts`
  （HWE 经 `HweCache`）取值，与先前 `fill_info_set` 的内联分支逐字节等价。
- 下游全部由表派生：
  - filter spec：`supported_fill_tags()`、`RecomputeSpec::all()`、`canonical_fill_tag`；
  - plan：`RecomputeFlags.tags: Vec<&'static RecomputeTag>` + `wants`/`enable`；
  - fill：`fill_info_set` 按表顺序写值；
  - header：`recompute_definitions` 用 `tag.ty`/`tag.number`/`tag.description` 生成 `##INFO`；
  - CLI：`--fill-tags` 帮助由 `fill_tags_help_names_every_registry_tag` 对着表断言，
    帮助与接受集不再能各自漂移。
- 结构门（`ci/check-structure.sh`）：`RECOMPUTE_TAGS` 全仓只定义一次；
  `crates/reed-filter/src/spec/recompute.rs` 必须引用该表；`recompute.ac/an/ns/maf/hwe`
  布尔字段与 `id: "AC"` 字面量不得回流。

### V30-R14 `SurvivorEdit`

- `reed-filter/src/eval/edit_plan.rs` 新增纯值
  `SurvivorEdit<'a> { keep_samples, mask_gt, info_remove, info_set }` 与唯一
  `edit_plan() -> EditPlan<'a>`。
- 三个发射路径（`single_pass.rs`、`filter/emit.rs`、`filtered.rs`）从各自拼
  `EditPlan { .. }` 改为构造该值；新增字段时不会漏掉某一路。
- 结构门：`crates/reed/src` 中 `EditPlan {` 字面量为 0；单元测试
  `survivor_edit_maps_every_plane_to_the_edit_plan` 锁四个平面的映射；
  契约测试 `single_pass_and_two_pass_reparse_produce_identical_output` 保持通过。

## 3.5 P3 门面与机器门（2026-11-12）

### 门面收口

- 删除 13 条整 crate/模块别名（`atomic_file`/`bgzf`/`budget`/`chunk`/`distance`/
  `expr`/`funnel`/`graph`/`header`/`metrics`/`plan`/`reduce`/`spec`）与 7 个 writers
  整模块；根模块只保留请求/报告/错误、`run_*`、`spec`（请求词表）、签名所需类型
  （`ThreadBudget`/`MemBudget`/`IoError`/`EngineError`/`BcfDict`/`Format`/`Region`/
  `ChunkSource`/`FilterSpec`/`SiteSelector`/`StatOutput`/`StageMetrics` 等）与
  `render_tsv`（报告示例的消费者）。
- 前端专用项统一收进 `#[doc(hidden)] pub mod frontend`（ADR §2 允许的唯一逃生门）：
  writers 七模块、`expr`/`plan`/`funnel`、`distance`/`graph`/`reduce`、
  `bgzf`/`chunk`/`header`/`tag_table`、`budget`/`metrics`、`PendingInput`/
  `open_regions`、`write_atomic`/`atomic_file`。`cargo public-api` 不列出
  `doc(hidden)`，快照从 1194 行降到约 1030 行；P3 的 `pub mod prelude` 再把逐项再导出展开，当前快照 1976 行。
- CLI/集成测试改为 `reed::frontend::…`（25 个文件）；公开签名类型全部根导出。
- 新增 `pub mod prelude`（P3 §1）：`use reed::prelude::*;` 一次导入 14 个
  `run_*`、请求/报告/错误/预算类型、`RunContext` 与 `OpenedInput`。实现只有
  一条 `pub use crate::{…}`，结构门禁止它直接再导出下层 crate；附可运行 doctest。

### 机器门

| 门 | 脚本 | 结果 |
|---|---|---|
| public-api 快照 | `ci/check-api.sh` | exit 0（快照 1976 行，含 `prelude` 展开） |
| 可命名性 | `ci/check-nameable.sh` | 22 个签名类型全部可在 `reed::*` 命名；探针 `crates/reed/tests/nameability.rs` |
| 外部类型 | `ci/check-external-types.sh` + `ci/external-types-allowlist.toml` | 0 errors（仅 workspace 模型 crate 在册） |
| feature | `ci/check-features.sh`（`cargo hack --each-feature`） | exit 0 |
| semver | `ci/check-semver.sh` | **CI 必过**；workspace 已 bump 到 0.5.0，对 v0.3.0 判定为 major（0.x 的破坏版本） |
| 结构 | `ci/check-structure.sh` | 新增 `prelude`/hidden `frontend` 模块门、RunContext 门与 R10 存活表/扫除门；31/31 |

`ci/check-doc-drift.sh` 的 "facade has no pub mod" 守卫改为 "唯一 pub mod 是
`#[doc(hidden)] frontend`"。CI 增加 nameability/feature 步与独立 external-types job。

## 3.6 P4 派生文档与结构（2026-11-12）

### 文档从代码派生

- 每个库 crate 的 `lib.rs` 现在有 ≥1 个可运行 doctest（此前 3 个；
  新增 reed-distance/engine/filter/graph/reduce/stats/oracle 各 1）。
  `cargo test --doc --workspace`：13 passed。
- `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --workspace` exit 0
  （修掉 private-intra-doc-link、失效链接、HTML 标签与 `line`/`grm`
  的 module/function 歧义）。
- 门面加 `# Stability`（根条目 + `run_*` 承诺；`frontend` 排除），
  library.md（en/zh）同步；入口表改为 8 个 `run_*` + 请求/报告两列。
- `ci/check-docs.sh`：rustdoc、doctest、`cargo modules orphans`（9 个库 crate
  全部 0）与入口表对 `pub use … run_*` 的一致性。

### 结构

- `convert/job.rs`（2142 行）拆为目录模块：`job/mod.rs`（类型/partials/
  OrderedSink/驱动，1044 行）、`job/convert_job.rs`（`ChunkJob` +
  `ConvertJob` impl，818 行）、`job/msa.rs`（等位表与 MSA 助手，293 行）。
  子模块可访问父模块私有类型，只有 8 个跨模块项加 `pub(super)`。
- 结构门新增：convert job 必须是目录模块；`crates/*/src/**/*.rs` 单文件
  ≤2000 行（当前最大 1850）。
- architecture.md 的 distance/reduce/graph/oracle 层级说明保留原样
  （doc_v25 已记录）；`reed-qc` 去向见 doc_v3 C4 守卫。

## 4. 未完成（按 03 §5 顺序）

| 条目 | 状态 | 原因/下一步 |
|---|---|---|
| R03 BCF `pos0` 哨兵 | **已完成（本批补齐）** | 原 §4 清单遗漏；实现见 §2 "V30-R03"，收敛门与测试见 §5 |
| R12 `reed-filter` 去 `reed-io` | **已完成（本批补齐）** | `LineRef`/`InfoIndex`/`info_fragments`/`parse_qual`/`count_alts`/`IntSpec`/`InfoValueOut` 已下沉 L0；见 §2 "V30-R12" 与 §5 的 `cargo tree`/结构门 |
| R05/R06/R07 能力接口 | **已完成（本批补齐）** | `StreamedInput/ReplayableInput/SeekableInput/InputMeta` + `IoError::Capability`；删除 `source_mut` 与按路径重开；见 §2 "V30-R05–R07" 与 §5 |
| R09 预算注入 | **已完成（本批补齐）** | `RunContext { budget, abort }` 统一传入八个 driver；删除 `run_budget.rs`、`FilterOptions` 第二预算与 Pca/Dapc 死字段；见 §2 "V30-R09" 与 §5 |
| R13/R14 | **已完成（本批补齐）** | `RECOMPUTE_TAGS` 驱动 spec/plan/fill/header/CLI help；`SurvivorEdit` 收敛三个发射路径；见 §3 "V30-R13/R14" 与 §5 |
| P3 门面收口 | **已完成（本批补齐）** | 13 别名 + 7 writers 移入 `#[doc(hidden)] frontend`；可命名性/外部类型/feature/semver/package 五门接入 CI；见 §3.5 |
| P4 文档/结构 | **已完成（本批补齐）** | 每 crate doctest、rustdoc `-D warnings`、`# Stability`、入口表门、`convert/job` 拆三文件；见 §3.6 |

## 5. 验证汇总（本批）

| 命令 | 结果 |
|---|---|
| `cargo fmt --all` | 已格式化 |
| `cargo clippy --workspace --all-targets -- -D warnings` | exit 0 |
| `cargo test --workspace` | exit 0（69 个 test-result 行，1206 passed，0 failed） |
| `cargo test -p reed-cli --bins` | 31 passed（含 `fill_tags_help_names_every_registry_tag`） |
| `cargo test -p reed --lib spill` | 16 passed（R10：共享存活表、启动扫除、死/活 pid、用户文件） |
| `cargo tree -p reed-filter -e normal` | 无 `reed-io` 节点（依赖只剩 `reed-core` 与第三方） |
| `rg "bcf\\.is_some\\(\\)"` / `rg "source_mut\\(\\)"` | 0 / 0 |
| `rg "default_half_memory" crates/reed/src` | 0 处（预算构造在调用方，structure 门 31/31） |
| `cargo test -p reed --test e2e_filter mind_strategy_reads_the_run_budget` | 1 passed（RunContext 的 16 MiB → two_pass_reparse，64 MiB → sketch_then_emit） |
| `strace -e openat`（`filter in.vcf.gz -r ... --strategy single_pass`） | 输入 `in.vcf.gz` 只 1 次 `openat` |
| `cargo test -p reed --lib` | 69 passed |
| `cargo test -p reed --test e2e_convert` | 42 passed |
| `cargo test -p reed --test e2e_filter` | 37 passed |
| `cargo test -p reed-io --lib bcf::record` | 9 passed（含 3 个 R03 用例） |
| `cargo test -p reed --test e2e_bcf --test e2e_index` | 17 + 11 passed |
| `rg -n "try_from\\((header\\.)?pos0\\)" crates --glob '*.rs'` | 只剩 `bcf/record.rs` 定义体内 1 处（另 1 条是文档说明） |
| `cargo check -p reed-reduce --no-default-features` | exit 0 |
| `bash ci/check-api.sh` / `check-doc-drift.sh` / `check-structure.sh` | 全部 exit 0（structure 30/30） |
| `bash ci/check-nameable.sh` / `check-external-types.sh` / `check-features.sh` | 全部 exit 0（21 个签名类型可命名；外部类型 0 error；feature 矩阵全绿） |
| `bash ci/check-docs.sh` | exit 0（rustdoc 0 warning；13 doctests；9 crate 无 orphan；入口表 8 项一致） |

> `04-recheck.md` 的 R01 行"convert `present` = vcf2phylip 口径"描述的是改动前状态；
> 本批按 03 的语义统一要求改成了 `has_called_allele`，半缺失现在与 `stats` 一致。
