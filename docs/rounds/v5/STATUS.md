# doc_v5 实施状态

> 对照 `07-roadmap.md` 的 X0–X7。每项注明 commit 与验证方法（测试名或可复现命令）；偏差单独登记。格式沿用 `docs/rounds/v4/STATUS.md`。

## 已完成

### X0 解除工作树阻塞（07 X0；01 F5-1 / 03 A5-4）

复审当时工作树有两处编译错误（`reed-cli` bin `E0382`、`reed-stats` 测试 `E0004`）。处置：在途 W5 代码是**窗口式** `StatKind::ProjectedSfs`，与 05 §3.1 的 genome-wide 路径冲突（01 F5-3），且其 CLI 默认（folded）、输出路径（`PREFIX.sfs`）、文件格式均与 05 §3.1 不符；按 07 X0 认可的"`git stash` 保主干干净"一支，把在途代码整体 stash，留待 X4 按设计重做。

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 在途 W5 代码 stash | `git stash push -u -m "w5-inflight-window-sfs ..."`（stash@{0}）；备份 `/tmp/v5_w5_inflight_tracked.patch`、`/tmp/v5_w5_writers_sfs.rs.bak` | `d735802` | `git status`：`crates/`、`Cargo.lock` 干净 |
| 主干可测 | 主干回到 `c5b6bfe` 的代码 | `d735802` | `cargo test --release --workspace` 364 passed / 0 failed；`cargo build --release --workspace --all-targets` 0 警告 |
| 守卫基线 | — | `d735802` | `bash ci/check-doc-drift.sh` 8/8 ok |

### X1a crate 与文档卫生（C-1/2/3/8/9、A5-2/3）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| C-1 删 M3 `qc.rs` | 删 `writers/qc.rs`（"M3 first consumer"，生产路径无调用者）与其测试/示例消费者 `tests/e2e_qc.rs`、`examples/qc_summary.rs`；跨线程确定性门改挂到 reports+summary 路径 | `42b4787` | `tests/determinism.rs::reports_and_summary_are_byte_identical_across_thread_counts`（1/2/3/7/56 线程，全部报告文件 + 渲染的 `.stats` 同一 sha256） |
| C-2 删 `reed qc` | 删 `Command::Qc`、`QcCli`、`translate/qc.rs`，以及只服务于它的 JSON 驱动（`QcOptions`/`QcOutcome`/`QcSink`/`run_qc_report*`，`QcJob` 的 `ChunkJob` impl）；守卫 `! grep Command::Qc` | `42b4787` | `crates/reed-cli/tests/stdin_streaming.rs::single_pass_stdin_streams_and_leaves_nothing`（改走 `stats --summary`）；`--help` 不再有 `qc` |
| C-3 `qc_report` → `summary` | `writers/qc_report.rs` → `writers/summary.rs`；`PREFIX.stats` 头部改为 "produced by reed stats"；`.qc.json` 随 C-2 一并消失（无写出者） | `42b4787` | `crates/reed-cli/tests/cli_contract.rs::report_mode_accepts_targets` 等；`reed::summary::render_bcftools_stats` |
| C-8 reed-cli 依赖 | 删 `reed-core`、`reed-stats`（`reed-io` 仍被源码使用，保留，见偏差） | `42b4787` | `cargo check -p reed-cli` |
| C-9 oracle 改名 | `crates/reed-core-v02-oracle` → `crates/oracle`，包名 `reed-oracle`、lib `reed_oracle`、`publish = false`；`reed-core`/`reed-stats` 依赖与全部测试引用同步 | `42b4787` | `grep -rn reed_core_v02_oracle crates` 空；测试全绿 |
| A5-2 unsafe | `bgzf.rs` 注释改为"两处审计例外"（另一处 AVX-512 内核）；`docs/rounds/v4/STATUS.md` 两行同步 | `42b4787` | 守卫 `allow(unsafe_code)` 恰 2 处 |
| A5-3 注释清理 | `funnel.rs` 的 `sample.max_missing` → `sample.mind`；`reed-stats/src/lib.rs` 的 `--fst-pops` → `--fst-global`；stats `--tmp-dir` 帮助删 "read twice" | `42b4787` | 守卫 `no sample.max_missing`、`no --fst-pops` |
| 守卫扩展 | 新增 4 条（两处 unsafe、无 `sample.max_missing`、无 `--fst-pops`、无 `Command::Qc`） | `42b4787` | `bash ci/check-doc-drift.sh` 12/12 ok |

### X1b 前后端边界（08 §6）+ L 系列（08 §4）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 门面收口（L-2 / 08 §6.2） | `reed/src/lib.rs` 改为**纯 `pub use` 清单**（0 条 `pub mod`）；`input`/`bcf_input`/`filtered`/`stats`/`filter`/`index` 均私有；`#![allow(unreachable_pub)]`（门面按构造如此）；删除 `qc`/`qc_report` 历史 re-export；重导出请求类型、`run_*`、writers、`spec/expr/plan/funnel`、`budget/metrics`、`bgzf/chunk/header` | `86a4350` | 守卫 `grep -c '^pub mod' lib.rs == 0`；`cargo doc --no-deps -p reed` 公开项即该清单 |
| 前端只依赖门面（08 §6.4.1） | `reed-cli/Cargo.toml` 只留 `reed`（+ clap/tracing/tracing-subscriber）；CLI 源码与测试的 `reed_filter::`/`reed_engine::`/`reed_io::` 全部改 `reed::` | `86a4350` | 守卫 `! grep -E 'reed-(core\|io\|engine\|stats\|filter)' crates/reed-cli/Cargo.toml` |
| 后端不打 stderr（08 §6.4.2） | 5 处 `eprintln!` → `tracing::warn!`（`input.rs`、`filter/mod.rs`、`bgzf.rs`×2、`plan/mod.rs`）；`reed-io` 增 `tracing` 依赖；CLI 订阅器在 verbosity ≥ 1 装 `warn` 过滤（3 为 `RUST_LOG`/debug） | `86a4350` | 守卫 `! grep -rqE 'eprintln!\|println!' <backend>/src`；`cli/verbosity.rs::levels_map_and_above_three_is_rejected` |
| L-1 README 重写 | 两份 README 改写为 v0.3 现实：三条命令、对拍表、性能表、仓库结构；删 v0.2 的 `WindowBuffer`/`exec::execute`/`--groups`/`reed qc`/legacy binary 示例 | `86a4350` | 守卫 `! grep -qE 'WindowBuffer\|exec::execute\|--groups' README.md` |
| L-3 去 noodles re-export | 删 `reed-io` 的 `pub use noodles::vcf as noodles_vcf`；内部使用者改 `use noodles::vcf`；`reed-io` 测试用 noodles（dev-dep） | `86a4350` | 守卫 `! grep -rq noodles_vcf crates` |
| L-4 发布元数据 | 8 个 crate 全部 `publish = false`（oracle 原有，其余 7 个新增）；`reed-core` description 更新 | `86a4350` | `grep -c 'publish = false' crates/*/Cargo.toml` 各为 1 |
| L-5 稳定性文档 | `docs/spec/naming.md` 的 "published surfaces" 改为"内部 API，仅对 CLI 与输出格式承诺稳定"；02-cli §0 增第 9 条（稳定性承诺、前后端规则） | `86a4350` | 人工 |
| L-6 examples | 删 `crates/reed-io/examples/{gt_scan,bcf_dump}.rs`（M2 验收驱动，无引用） | `86a4350` | `cargo build --examples` 无目标 |
| L-7 doctest | 门面 3 个 doctest（`FilterSpec` 往返、`StatsRequest` 构造、`run_filter` 驱动） | `86a4350` | `cargo test --doc -p reed` 3 passed |

X1 门：`cargo test --release --workspace` 全绿；`cargo build --release --workspace --all-targets` 0 警告；`bash ci/check-doc-drift.sh` 17/17 ok；`cargo test --doc -p reed` ≥ 3。

### X2 目录搬迁（06-layout.md §4/§5）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 文档树 | `docs/{architecture-spec,filter-architecture,naming-conventions,statistics-reference-map}.md` → `docs/spec/{architecture,filter,naming,statistics-map}.md`；`docs/INDEX.md` → `docs/README.md`；`DESIGN.md` → `docs/rounds/v1/design.md`；`doc_v{2,3,4,5}/` → `docs/rounds/v{2,3,4,5}/`；`docs/spec/cli.md` = 合并版 CLI 规范（`docs/rounds/v5/02-cli.md` 的当前副本） | `2e11c01` | `ls docs`；守卫 `! ls -d doc_v[0-9]` |
| 度量/验证树 | `bench/{e2e,hts,noo}` → `bench/suites/{e2e,vs-htslib,vs-noodles}`；`bench/Cargo.toml` members 同步；`bench/target` 删除；`scripts/` → `ci/`（`check-doc-drift.sh`、`publish.sh`；`publish.sh` 的 `PUB_PATHS` scripts→ci） | `2e11c01` | `cargo metadata --manifest-path bench/Cargo.toml --no-deps` OK；守卫 `[ ! -d scripts ] && [ -f ci/check-doc-drift.sh ]` |
| 生成物出仓库 | `gold/work`（54 G）、`gold/out`（25 G）→ `target/gold/{work,out}`（同文件系统 `mv`，0.12 s）；新增 `gold/scripts/move_scratch.sh`（可重复、`--dry-run`）；`lib.sh` 用 `REED_SCRATCH`（默认 `target/gold`），`gold/out/work` 套娃并入 `target/gold/work`；`.gitignore` 只忽略 `gold/{work,out}/`，`gold/` 其余随仓库跟踪 | `2e11c01` | 守卫 `[ ! -d gold/work ] && [ ! -d gold/out ]`；`INPUT=/dev/null source gold/scripts/lib.sh` 解析出 `target/gold/{out,work}` |
| 引用一次 sed | `*.md/*.rs/*.sh/*.py/*.toml` 内 `doc_vN/` → `docs/rounds/vN/`、四个 spec 文件名、`scripts/` → `ci/`、bench 子目录；`docs/README.md` 重写为"当前规范/历史"两表（含 v5 与 spec/cli.md） | `2e11c01` | 守卫 `! grep -rq "doc_v[0-9]/" docs/README.md`；`bench/suites/e2e/results.tsv` 的引号内历史注记除外 |
| CI（X7） | `.github/workflows/ci.yml`：build/test/doctest/clippy/drift 守卫（可移植）；子集数字门与 gold 对拍作为数据门，条件化到带 `reed-data` 的 self-hosted runner | `2e11c01` | workflow 语法；本地五个可移植门等价命令全绿 |

X2 门：`cargo test --release --workspace` 全绿；`bash ci/check-doc-drift.sh` 21/21 ok；新路径 `run_all.sh` 烟测（`REGION=Chr01:1-2000000`）三方各 4,051 位点、Jaccard 1.0、同名样本基因型零不一致，168 s。**偏差**：全量 `run_all.sh` 的正文 md5 门未在本轮重跑（bcftools 分支约 40 min）；以同路径的小区域烟测证明脚本在新路径可跑通。

## 未完成 / 偏差

| 项 | 状态 | 说明 |
|---|---|---|
| C-8（部分） | 偏差 | `reed-io` 被 reed-cli 源码直接使用（`Region`、`tag_table`、`bgzf`），未删；C-9 的 7 项里仅 `reed-io` 保留。已由 08 §6 的门面 re-export 覆盖（CLI 现经 `reed::` 使用） |
| 08 §6.3 前端含后端逻辑 | 未做 | `translate/stats.rs` 仍调 `reed::plan::compile` 并自行 `assign_groups`；`translate/filter.rs` 的 `site_spec` 编译同理。移入 `run_filter`/`run_stats` 会改 `run_*` 签名（`StatsRequest` 需携带 `FilterSpec`/pop 文件），留待后续 PR |
| 08 §6.2 单一 `reed::Error` + `kind()` | 未做 | 现有 `FilterError`/`EngineError`/`IoError`/`CoreError` 未合并为 `reed::Error`；退出码仍由前端各自映射。留待后续 PR |
| 08 §6.4.4 cli_contract 往返 | 未做 | 现有往返测试只覆盖 filter；stats/index 的往返测试未补 |
| L-2 完全精选 | 部分 | 门面以 `pub use` 重导出模块（`spec/expr/plan/funnel`、`budget/metrics`、`bgzf/chunk/header`）与类型；`unreachable_pub` 为 crate 级 allow。未逐项把非精选项改 `pub(crate)`（量大，留待后续） |
| X2–X7 | 待做 | 见 `07-roadmap.md` |

## 基准（2026-10-09 v5 复审基线，56 线程，热缓存；负载与中位数见 `04-performance.md` §1）

（复审时的数字由 `04-performance.md` 给出；实现者在每步完成后按同一命令重测并补登 `bench/suites/e2e/results.tsv`。X0/X1 未触及计算路径，未重跑全量基准。）
