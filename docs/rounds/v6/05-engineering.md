# 05 工程：构建、测试、守卫、CI、依赖、文档（2026-10-10 实测）

> 范围：仓库作为一个工程产品的健康度。数字全部来自本轮实测（`/tmp/v6_health.log`）与仓库文件；每条判断给出 `file:line`。

---

## 1. 健康状态

| 检查 | 结果 | 备注 |
|---|---|---|
| `cargo build --release --workspace --all-targets` | **0 警告 / 0 错误** | 工作树干净（`8709f1f`） |
| `cargo test --release --workspace` | **382 通过 / 0 失败**（27 个测试目标；`determinism` 单套 86 s） | |
| `cargo clippy --release --workspace --all-targets -- -D warnings` | **失败** | reed-filter lib 2 处 `doc_lazy_continuation`（`crates/reed-filter/src/eval/tail.rs:75-76`），在 `-D clippy::all` 下是 **error**；此外测试目标的 pedantic 警告约 60 条（`manual_assert_eq` 22、`cast_possible_truncation` 22、`ignored_unit_patterns` 6、`naive_bytecount` 5…）在 `-D warnings` 下同样变成 error。`docs/rounds/v5/STATUS.md` 写的 "clippy 0 error" 用的不是 CI 的命令 |
| `bash ci/check-doc-drift.sh` | 23/23 ok | 含两处 unsafe、前后端三条、目录搬迁、X3 等 |
| gold 全量对拍 | PASS，正文 md5 `7aedc59ac798729a2674e6e7ef19adfc` | 04 §2 |
| 跨格式 / 跨线程 | stats 五统计 VCF≡BCF；filter 8/16/32/56 线程 user 恒定 | 04 §1、§3 |

结论：**产品路径是绿的，CI 的 lint 门是红的**。`.github/workflows/ci.yml:26-27` 正是上面失败的那条命令，所以 GitHub 上的 CI 当前必然失败（远程 `origin` 指向 `github.com/WWz33/reed`）。

## 2. CI 与发布

`.github/workflows/ci.yml`：build / test / doctest / clippy / drift 守卫（`:20-29`），数据门（子集 263,420 / 381,672、gold 对拍）条件化到 `reed-data` self-hosted runner（`:34-47`）。缺：

- `cargo fmt --check`（`rustfmt.toml` 存在但无人校验）；
- `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps`（已知两处坏链：`reed-core/src/bitset.rs:3` 的 `crate::bits::Bits`、`reed-core/src/count/site_counts.rs:64-65` 的 `super::hwe::hardy_weinberg_p_value`）；
- MSRV 任务（`rust-version = "1.91"`，工具链固定 1.99；clippy 规则随版本变，lint 门必须在固定工具链上跑，这一点 ci.yml 通过 `rust-toolchain.toml` 做到了，但 MSRV 本身从未构建过）；
- `cargo deny check`（`deny.toml` 只有许可证白名单与 `multiple-versions = "warn"`，CI 没跑它）；
- 基准回归（`bench/suites/*` 与 criterion 基准没有任何自动化）。

发布：`ci/publish.sh` 把**剔除 `docs/` 与 `bench/`** 的子树强推到 `origin/main`（`publish.sh:1-9`）。后果：公开仓库没有规范、没有基准台账，只有代码与 README；CI 在公开树上跑时 `ci/check-doc-drift.sh` 里引用 `docs/…` 的守卫会失效或误判。这是一个有意的策略（内部文档不公开），但要在 README 与 CI 里写明并做成两套守卫（公开树一套）。

## 3. 测试形态

| 层 | 内容 | 位置 |
|---|---|---|
| 单元 | 模块内 `#[cfg(test)]`，核心 69、引擎 9、io/filter/stats 若干 | 各 `src/` |
| 性质（proptest） | 累加器合并律、分箱网格、位平面、扫描器 fuzz（2,000 用例 + 突变） | `reed-stats/tests/accumulator_laws.rs`、`reed-core/tests/{grid,genotype_planes}.rs`、`reed-io/tests/genotype_scan_fuzz.rs` |
| 差分 oracle | 冻结的 v0.2 实现（`crates/oracle`）；noodles 作扫描器与 BCF 渲染/编码的参照 | `reed-stats/tests/oracle_parity.rs`、`reed-io/tests/{genotype_scan_vs_noodles,bcf_vs_noodles,bcf_encode}.rs` |
| 端到端 | 13 个文件：跨线程确定性、跨格式 md5、缺失阈值边界（plink2/vcftools）、`-q :type`（bcftools）、SFS/f 统计/fst-global/LD（easySFS/allel/vcftools/plink2 脚本）、索引、stdin | `crates/reed/tests/*.rs`、`gold/popgen/*` |
| CLI 契约 | 帮助文本、参数往返（仅 filter）、verbosity、stdin | `crates/reed-cli/tests/{cli_contract,stdin_streaming}.rs` |
| 基准 | criterion：扫描器、行编辑（各 1 条记录）；外部套件：BCF 读通量 vs htslib / noodles | `crates/reed-io/benches/*`、`bench/suites/{vs-htslib,vs-noodles}` |

缺口（按风险）：

1. 引擎边界：一条记录跨 **>2 个 chunk**、空输入、`RecordFraming::Raw` 走引擎的测试、`FrameChain` 的毒化/多等待者/截断检查（引擎审计 §5.4）；文本路径的接缝片段无上限（`pipeline.rs:533-574`）没有测试。
2. 两种 chunk 源的等价性：`MmapChunkSource` 没有任何直接测试；流式源在丢弃空成员后 `compressed_start` 漂移（`reed-io/src/chunk/bgzf.rs:131-135,310`），没有"流式 vs mmap 的 chunk 元数据相同"的测试（03 §5）。
3. 索引内容：只有索引种类解析的测试，没有"块边界上的记录 / 跨块记录用 `bcftools view -r` 能查到"的内容测试。
4. 往返：stats 与 index 的 `to_args()` 往返测试缺（08 §6.4.4 遗留）。
5. fuzz：只有扫描器有 proptest；`LineRef`/INFO 解析器/表达式词法器/BCF 记录解析器没有 fuzz 目标。
6. 基准不进 CI；两个 criterion 基准只覆盖 2 个内核，deflate/inflate/渲染/编码/索引无基准。

## 4. 文档与代码的漂移（本轮核实）

| 位置 | 说的 | 实际 |
|---|---|---|
| `crates/reed/src/filter/ld_prune.rs:7-9` | "contig 各自线程、并行剪枝" | 一个 `Option<LdPruner>` 在 reducer 线程串行（`:318-369`；04 §6） |
| `docs/rounds/v5/STATUS.md` 最终验证 | clippy 0 error | CI 命令下 2 error（§1） |
| `crates/reed-io/src/vcf/line_edit.rs:8-10` | 无编辑的记录不该走 `edit_vcf_line` | 两处调用点无条件调用（`single_pass.rs:326,388`） |
| `crates/reed-io/src/error.rs:61-63` | 引用 `crate::ReadError` | 类型已不存在 |
| `crates/reed-core/src/tag_table.rs:91-105` | "按声明顺序" | 实际按名字排序 |
| `crates/reed-stats/src/individual.rs:44-47` | "按字 popcount、无逐样本循环" | 逐样本循环 |
| `crates/reed-core/src/count/kernel.rs:99-104` SAFETY | "dispatcher 断言长度" | release 下只有 `debug_assert`（03 §5） |
| `docs/spec/architecture.md:804`、`reed-engine/src/budget.rs:46-50` | `threads ≤ 2` 走顺序模式 | `run()` 无条件 `run_parallel`，`is_sequential` 无调用者 |

机制上的教训：守卫脚本只能抓"存在/不存在"的事实；"并行/串行"这类行为陈述要靠 `--timing` 的数字门（04 §9）或测试断言，不能靠注释。

## 5. 依赖与策略

| 项 | 现状 | 判断 |
|---|---|---|
| noodles 0.117 | 只在 reed-io（header 解析、索引读写、BCF header 文本）与 reed（header 类型）；热路径无 noodles；测试作 oracle | 合理；但 `HeaderFacts.header`/`FileIndex`/`Indexer` 把 noodles 类型放进公开签名，升级 noodles 即破坏内部 API |
| libdeflater 1.26.1 | 只用 `CompressionLvl`；每块 `resize(…,0)` | 04 §4.2 的浪费 |
| memchr 2.8 | 用于 tab/`:`/`memmem` | 但 1 MiB chunk 的尾部定位用标量 `rposition`（`chunk/mod.rs:231-246`） |
| rayon 1.12 | 只在 `filter/mid_stage.rs:78` | 保留 |
| crossbeam-channel、bitflags、smallvec、bstr、serde/serde_json、sha2、tracing、clap | 用途明确 | — |
| `publish = false` | 8 个 crate 全部 | 与 08 §6 一致 |
| `unsafe_code = deny` | 两处例外（AVX-512 内核、mmap），守卫锁定 | 内核例外的 SAFETY 不成立（03 §5），且该路径在 348 样本下永不执行（W=6 < 8） |
| edition 2024 / rust-version 1.91 / 工具链 1.99 | — | MSRV 从未验证 |

## 6. 工程建议（E6-1…E6-9，门写在 06）

| # | 建议 | 做法 |
|---|---|---|
| E6-1 | **把 lint 门变绿** | 修 `tail.rs:75-76` 的两处 doc lint；测试目标的 pedantic 要么修（`assert!(a==b)`→`assert_eq!`、`()` 模式、test 内 cast 用 `try_from`），要么把 `-D warnings` 限定到 lib/bin 目标而测试目标用 `-W`；二者选一并写进 ci.yml |
| E6-2 | CI 补四条 | `cargo fmt --check`、`RUSTDOCFLAGS=-D warnings cargo doc --workspace --no-deps`、`cargo deny check`、MSRV 构建（`1.91`） |
| E6-3 | 基准进 CI | `bench/suites/e2e/run_baseline.sh`（04 §9）在 `reed-data` runner 上跑并把 wall/user 与上一次比，回退 >5% 阻断 |
| E6-4 | 引擎与 chunk 源的测试 | 跨 >2 chunk 的记录、空输入、Raw 框架、FrameChain 毒化；`BgzfChunkReader` vs `MmapChunkSource` 元数据等价；索引内容测试（`bcftools view -r` 查块边界记录） |
| E6-5 | fuzz 目标 | `cargo-fuzz`：`LineRef::parse`、INFO 片段、表达式词法/语法、`record_ref`（BCF） |
| E6-6 | 往返测试补齐 | stats / index 的 `to_args()` 往返进 `cli_contract` |
| E6-7 | 文档漂移清单归零 | §4 的 8 条；`cargo doc -D warnings` 守住坏链 |
| E6-8 | 公开树的守卫 | `publish.sh` 剔除 docs 后，公开树的 ci.yml 不能引用 `docs/`；做两套守卫或公开 `docs/spec/` |
| E6-9 | 死代码清理 | `reed-io/src/field_plan.rs`（无消费者）、`reed-core` 的 `CountTable`/`count_site_group`（无生产调用者，见 03 §4 的取舍）、`math::hypergeom`（与 `reed-stats::sfs::projection_weights` 重复且规则不同）、`ChunkBytesExport`、`MemBudget::from_half_memory_for_tests`、`is_sequential`（要么实现顺序模式，要么删） |

---

## 7. 目录树、文件名、标识符（2026-10-10 实测）

### 7.1 目录树

根目录：`bench ci crates docs gold` + `.cargo .github` + `README.md README.zh-CN.md` + `Cargo.toml Cargo.lock rustfmt.toml clippy.toml deny.toml rust-toolchain.toml`；二级：`crates/{reed-core,reed-io,reed-engine,reed-stats,reed-filter,reed,reed-cli,oracle}`、`bench/suites`、`docs/{spec,rounds}`、`gold/{config,fixtures,popgen,scripts}`；生成物全在 `target/`。与 v5 06 的目标树一致，**合格**。唯一的旁枝是 `bench/.cargo`（bench 是独立 workspace，自带配置，可接受）。

### 7.2 文件与模块名

| 项 | 证据 | 判断 |
|---|---|---|
| 超长文件 | `reed/src/stats/job.rs` **1,936**；`reed-cli/src/translate/stats.rs` 842；`reed/src/stats/writers/reports.rs` 823；`reed/src/filter/mod.rs` 766；`reed-filter/src/plan/mod.rs` 766；`reed-io/src/chunk/bgzf.rs` 711；`reed-io/src/bcf/native_encode.rs` 683；`reed/src/filter/record_pipeline.rs` 610 | 8 个文件 > 600 行；可交付项目的常规上限约 500。`job.rs` 必拆（Y6 已列），其余按责任拆 |
| 模块风格 | 18 个 `mod.rs`；1 个 sibling 风格（`reed-cli/src/cli.rs` + `cli/`） | 一个离群点，统一成 `cli/mod.rs` |
| 测试文件名 | `reed-stats/tests/m7_fst.rs`（里程碑编号命名）、`qc_merge_laws.rs`（`qc` 是已退役词汇）、`individual_joint.rs`（随 C-5 的模块改名） | 改 `fst_multi.rs`、`summary_merge_laws.rs`、`het.rs` + `sfs_joint.rs` |
| 厂商词汇进生产代码 | `reed/src/stats/writers/pixy.rs`、门面 `pub mod pixy`、`pub fn pixy_kinds`、`pub struct PixyOptions` | 违反 `docs/spec/naming.md` §2.5"公开 API 不出现厂商词汇"与用户裁决"软件名只在测试里"。pixy 只是**列布局**的来源，应是 `writers/windows.rs`、`WindowLayoutOptions`、`window_kinds`，布局来源写在文档。测试里的 `e2e_pixy.rs`、`*_vs_noodles.rs`、`bench/suites/vs-htslib` 合规 |
| 标记词 | `# ponytail` ×3（`reed-core/src/math/hwe.rs:134`、`reed-stats/src/sfs.rs:153`、`reed-stats/src/jackknife.rs:8`） | 非标准标记，工具不识别；改 `TODO(#issue)` 或删 |

### 7.3 标识符（对照 `docs/spec/naming.md` §2）

| 规则 | 结果 |
|---|---|
| §2.10 禁用缩写 `buf/rec/idx/src/tmp/cnt` 不进公开项与字段 | 公开类型/函数/字段 **0 处** ✓；公开函数**参数名** 4 处（`reed-io/src/bcf/gt.rs:46` `rec`、`reed-io/src/bgzf/deflate.rs:41` `src`、`reed-stats/src/summary/histogram.rs:144,184` `idx`），参数名出现在 rustdoc 里，顺手改 |
| §2.1 计数 `n_*`、阈值 `min_/max_`、布尔 `is_/has_/require_` | 领域类型合规；`reed-cli/src/cli.rs` 的布尔字段按旗标名（`freq/het/csi/force/no_update`）是 clap 派生的可接受例外 |
| §2.6 坐标基 | `reed-core/src/coord.rs:21-32` 用 `Pos0`/`Pos1` 类型 ✓（好做法）；`reed-filter/src/expr/mod.rs:63` `pub pos: u64`、`spec/region.rs:34-36` `start/end: Option<u64>` 名字不带基，依赖文档注释（未逐条核对注释是否写明） |
| §2.7 每 crate 一个 `*Error` | `CoreError`/`IoError`/`EngineError`/`FilterError` ✓；reed-stats 无自有错误类型（用 `Option`/`CoreError`，可接受）；门面无统一 `reed::Error`（Y6） |
| §2.9 模块/函数 stutter | `hwe::{hwe_exact_p, hwe_midp_p, hwe_distribution}`（规则自注"0.2 暂容忍"） |
| plink 式 `_ct` 后缀、`#[allow(dead_code)]` | 0 ✓ |
| `cargo fmt --all --check` | 0 diff ✓ |

### 7.4 像不像可交付的 Rust 项目

| 项 | 现状 | 判断 |
|---|---|---|
| 许可证 | `Cargo.toml` 声明 `license = "MIT OR Apache-2.0"`，仓库里**没有** `LICENSE*`/`COPYING*` 文件 | **缺陷**：声明了许可证却无文本，`cargo deny`/GitHub 都识别不了；加 `LICENSE-MIT` + `LICENSE-APACHE` |
| 版本与发布 | 版本 0.3.0 lockstep；`git tag` 为空（`git describe` 只给哈希）；无 CHANGELOG；`reed --version` 只打 `reed 0.3.0`，无 commit/日期（无 `build.rs`） | 打 `v0.3.0`；Keep-a-Changelog；`--version` 带 `commit, date, rustc`（`build.rs` 注入 `git describe`） |
| 安装 | README 只有"Build and test"（`cargo build`），无 Installation（`cargo install --path crates/reed-cli`）、无预编译发布件、无 conda 配方 | 加安装节；release workflow 产出 `reed-x86_64-unknown-linux-gnu.tar.gz` |
| README 内容 | Commands / Verification / Performance / Repository layout / Build and test / License 齐全；缺 Installation、Citation、输入要求（header 需 `##contig` 长度、索引要求）、退出码 | 补四节 |
| 元数据 | 8 个 crate `publish = false`；`repository/readme/keywords` 缺 | 不发布可不填，`repository` 建议填 |
| 构建配置 | `[profile.release]`：`lto = "fat"`、`codegen-units = 1`、`debug = "line-tables-only"`、`panic = "unwind"`（`catch_unwind` 所需）；`profile.profiling` 继承 release 带符号；`rust-toolchain.toml` 固定 1.99；`Cargo.lock` 入库 | 合格；MSRV 1.91 未验证（§2） |
| 二进制 | `[[bin]] name = "reed"`（`reed-cli/Cargo.toml:16-18`），`cargo install --path crates/reed-cli` 可装 | ✓ |
| 治理文件 | 无 CONTRIBUTING / CODEOWNERS / issue 与 PR 模板；`.github/` 只有 `ci.yml` | 可选；CONTRIBUTING 至少写清 gold 数据门与 STATUS 纪律 |
| 公开树 | `publish.sh` 剔除 `docs/`、`bench/`（§2） | 见 E6-8 |

**结论**：代码层面像一个可交付的 Rust 项目（fmt 干净、lint 与 unsafe 策略、测试形态、release profile、分层 workspace、守卫脚本）；**发布层面不像**：没有许可证文本、没有 tag/CHANGELOG/安装说明、版本号不带构建信息。这些加起来是一天的工作，全部放进 Y0。

## 8. 工程建议补充（E6-10…E6-15）

| # | 建议 | 门 |
|---|---|---|
| E6-10 | `LICENSE-MIT` + `LICENSE-APACHE`；`cargo deny check licenses` 进 CI | `ls LICENSE-MIT LICENSE-APACHE` |
| E6-11 | `v0.3.0` tag；`CHANGELOG.md`（Keep a Changelog，从 v2 轮次的 STATUS 回填）；`build.rs` 注入 `git describe --always --dirty` 与日期，`reed --version` 打 `reed 0.3.0 (8709f1f 2026-10-10, rustc 1.99)` | `reed --version` 含 commit |
| E6-12 | README 加 Installation / Citation / Input requirements / Exit codes 四节（中英两份同步） | 人工 |
| E6-13 | `pixy` 词汇出生产代码：`writers/pixy.rs` → `writers/windows.rs`，`PixyOptions` → `WindowLayoutOptions`，`pixy_kinds` → `window_kinds`；门面 re-export 同步 | 守卫：`! grep -rniE '\bpixy' crates/*/src --include=*.rs`（测试目录除外）；输出 md5 不变 |
| E6-14 | 测试文件改名（`m7_fst`→`fst_multi`、`qc_merge_laws`→`summary_merge_laws`、`individual_joint`→`het`+`sfs_joint`）；`cli.rs` → `cli/mod.rs`；`ponytail` → `TODO(#…)`；4 个公开参数名去缩写 | `cargo test` 全绿；`! grep -rn ponytail crates` |
| E6-15 | 超长文件拆分：`stats/job.rs`（Y6 的 `request/run/sink`）、`translate/stats.rs`（按子命令段）、`writers/reports.rs`（按报告族）；目标无文件 > 600 行 | `find crates -name '*.rs' -path '*/src/*' \| xargs wc -l \| awk '$1>600'` 为空 |

### 7.5 `crates/` 内部树与源码形态（对照 noodles 与 rust-htslib）

参照物：`~/reed/noodles`（18 个 crate 的 workspace，`noodles-vcf/src` 全部 sibling 风格 `foo.rs` + `foo/`，0 个 `mod.rs`，最大文件 702 行，`lib.rs` = `pub mod` 每概念一个 + `pub use` 三个中心类型）；`~/reed/rust-htslib`（全部 `mod.rs` 风格，最大文件 9,672 行，带 `LICENSE.md`、`CHANGELOG.md`）。两种模块风格社区都接受，**混用**才是问题。

**合格的部分**

| 项 | 证据 |
|---|---|
| 分层与依赖方向 | `core ← io ← engine`、`core ← stats`、`core,io ← filter`、全部 ← `reed` ← `reed-cli`（只依赖 `reed`）；无环、无越层（各 `Cargo.toml`） |
| 纯度 | reed-stats 不碰 I/O、线程（grep 为空）；reed-core 唯一的文件读取是声明过的 `GroupTable::from_tsv`（`group.rs:97-103`）；reed-filter 只用 reed-io 的 4 个项（`LineRef`、`info_fragments`、`parse_qual`、`line_edit`），接近纯策略层 |
| 外部依赖按层收口 | noodles/libdeflater/memmap2 只在 reed-io；crossbeam 只在 engine；clap 只在 reed-cli；serde 在 core/stats/filter；rayon 只在 reed |
| 文档 | **131/131** 个 `src/*.rs` 都有 `//!` 模块文档；每个公开项有 `# Errors`/`# Panics` |
| 测试形态 | 56/131 文件带内联单元测试；每 crate 有 `tests/`；proptest、oracle crate、差分 oracle、e2e；`e2e_filter.rs` 21 个用例 |
| 清单卫生 | 7 个产品 crate 都有 `description`、`rust-version.workspace`、`[lints] workspace = true`（oracle 缺后两项，冻结 crate 可接受） |
| 门面形态 | `reed/src/lib.rs` 六个私有 `mod` + 纯 `pub use` 清单；`reed` 内 `pub(super)` 93 / `pub(crate)` 44，`reed-cli` 裸 `pub` 仅 1 处——可见性纪律比多数开源 Rust 项目严 |
| 文件尺寸 | 除 8 个离群点外 100–500 行，与 noodles 同量级，远好于 rust-htslib |

**不合格 / 不一致的部分**

| # | 项 | 证据 | 处理 |
|---|---|---|---|
| T-1 | 模块风格混用 | 18 个 `mod.rs` + 1 个 sibling（`reed-cli/src/cli.rs` + `cli/`）；noodles 全 sibling、rust-htslib 全 `mod.rs`，没有谁混用 | 统一；最省是把离群点改成 `cli/mod.rs`，并在 naming.md §2.9 写明"本仓库用 `mod.rs` 风格" |
| T-2 | crate 根的三种风格 | reed-core / reed-stats：`pub mod` 全部 + 主类型扁平 `pub use`（同一类型两条路径）；reed-io / reed-engine：只 `pub mod`（+ 错误类型）；reed-filter：`pub mod` + 把 `FilterError` **直接定义在 `lib.rs:28`**（其余三个 crate 都有 `error.rs`）；门面：隐藏模块只 `pub use`。noodles 的做法是 `pub mod` 每概念 + `pub use` 少数中心类型 | 统一为"`pub mod` 每概念一个 + `pub use` 2–4 个中心类型 + `error.rs`"；`FilterError` 搬到 `reed-filter/src/error.rs` |
| T-3 | 过期的 `#[doc(hidden)]` | `reed-stats/src/lib.rs:44-60` 注释写"CLI 尚未接入（`--fst-global`、`--sfs --proj`、LD/f 统计）"并把 `WcMultiAcc`、`ProjectedSfsAcc`、`LdPruner`、`F2Acc…`、`BpBlockJackknife` 全部藏起来；X4 已把它们接进 CLI 并对拍 | 去掉 `#[doc(hidden)]`，改注释；未接入的只剩 `JointSfsAcc` |
| T-4 | 尺寸离群点 | `reed/src/stats/job.rs` **2,039**、`writers/reports.rs` 853、`reed-cli/translate/stats.rs` 870、`filter/mod.rs` 824、`plan/mod.rs` 821、`chunk/bgzf.rs` 734、`native_encode.rs` 714、`record_pipeline.rs` 623；测试 `e2e_filter.rs` 1,188 | E6-15；测试按特性拆 |
| T-5 | `reed/src/filter/` 的概念缺目录 | 9 个文件平铺：三套 pass0（`reduce.rs`、`sketch.rs`、`ld_prune.rs`）与 `mid_stage.rs`、`emit.rs`、`single_pass.rs` 同级，看不出"pass0 / pass1"的结构 | 03 §7 的 `pass0/` 与 `emit/` 目录随 Y1 落地 |
| T-6 | `reed-stats` 一文件一统计量的例外 | `individual.rs` 装两个统计量（C-5） | Y6 |
| T-7 | 测试目录命名不统一 | 只有 `crates/reed/tests` 用 `e2e_` 前缀；`m7_fst.rs`、`qc_merge_laws.rs` 历史名 | E6-14 |
| T-8 | 基准只在 reed-io | `benches/` 只有 reed-io 两条；引擎、统计、过滤热路径无微基准 | E6-3 随 Y2/Y3 补 |

**判定**：`crates/` 的树与源码**合格偏上**——分层、纯度、文档、测试、清单、可见性都达到或超过 noodles 的水准；拉开差距的是三处一致性债（T-1 模块风格、T-2 根风格与 `FilterError` 位置、T-3 过期的 `doc(hidden)`）和 8 个尺寸离群点。前三项是半天的机械活，并入 Y0；尺寸问题随 Y1/Y6 的拆分解决。

## 8 补充：E6-16、E6-17

| # | 建议 | 门 |
|---|---|---|
| E6-16 | 统一模块与 crate 根风格：`cli.rs` → `cli/mod.rs`；每 crate 根 = `pub mod` 每概念 + `pub use` 中心类型 + `error.rs`；`FilterError` 搬入 `reed-filter/src/error.rs`；naming.md §2.9 写明约定 | `find crates -name '*.rs' -path '*/src/*' -mindepth 3` 无 sibling 对；`ls crates/*/src/error.rs` 四个产品 crate 都有 |
| E6-17 | 去掉 `reed-stats/src/lib.rs:44-60` 已接入项的 `#[doc(hidden)]`，注释改为只列 `JointSfsAcc` | `cargo doc -p reed-stats` 的索引含 `WcMultiAcc`、`ProjectedSfsAcc`、`LdPruner`、`F2Acc` |
