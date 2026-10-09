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

### X3 reports 门（03 §5）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| summary 消费共享扫描 | reports job 在请求含 `--summary` 时用 `.tracking_alleles()` 扫描器并解码 DP（`depth \|\| summary`）；`QcAccumulator::push_scan(&ScanOutput, &LineRef, &SiteView)` 取代自建扫描——删 `QcJob` 的 scanner/`QcScratch`/`push_line` 与 `GenotypeScanner` 依赖；`reports.rs` 把同一 `out` 交给 summary | `23ea380` | 守卫 `! grep -qE 'GenotypeScanner\|QcScratch' summary.rs` |
| 输出不变 + 门 | 全量 `clean.vcf.gz --freq --hardy --het --missing --depth --summary`（56 线程）：重构前后 8 个报告文件 md5 逐一相同（afreq `b2884092…`、hardy `fd8fbf27…`、het `790293d9…`、smiss `5f4fe4f8…`、vmiss `7c325d64…`、idepth `a4b6dd3a…`、ldepth.mean `475ee083…`、stats `189e025b…`）；wall 3.75 → **3.037 s**、CPU 184 → **152.0 s**（门 ≤ 3.2 s / ≤ 155 s；3 次 3.037/3.037/3.044 s、151.9/152.1/152.0 s，极差 ≤ 0.3%；批次开始负载 4.63 为自抬升，按 04 §4 纪律标注） | `23ea380` | 上述 md5 逐文件对比；`e2e_reports`、`determinism`、`e2e_bcf` 全绿 |

### X4.1 投影 SFS（05 §3.1）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| genome-wide 请求路径 | `StatsRequest.sfs: Option<SfsRequest { projections, folded }>`；job 在 `push_row` 里喂每个群体的 `ProjectedSfsAcc`（`ProjectionCache` 记忆权重），partial 与窗口 partial 并列上交，sink 收尾合并；`StatsReport.sfs: Option<SfsReport>`。不是窗口统计，不进 `_*.txt` 行模型 | `a49d318` | `crates/reed/tests/e2e_sfs.rs`（3 项）；buf 与 BCF 原生输入共用 `push_row`，`clean.bcf --sfs --proj 428` 三群体 1.1 s 跑通 |
| 累加器全支撑 | `ProjectedSfsAcc::bins` 扩到 `n+1`（folded: `n/2+1`），新增 `folded()`/`spectrum()`；`pi()`/`segregating()` 仍只用 `1..n`；`ProjectionCache` 加 `Clone` | `a49d318` | `cargo test -p reed-stats --test projected_sfs`（7 项，含改写后的全支撑断言） |
| 输出格式 | `PREFIX_sfs/dadi/<pop>-<n>.sfs`（`<n+1> unfolded\|folded` / 谱 / mask；folded 保持 n+1 长度、上半补零，mask 标 k=0 与结构零区）；`PREFIX_sfs/fastsimcoal2/<pop>_{DAF\|MAF}pop0.obs`（`1 observation` / `d0_0…d0_n` / 谱）；`--sfs-int` 取整 | `a49d318` | 与 easySFS 的实际文件逐字段对齐（`gold/popgen/compare_sfs.py`） |
| CLI | `--sfs [folded] --proj N[,N…]`（默认 unfolded；单值广播到各群体）、`--sfs-int`；无 `--pop` 时全体样本为一个群体（名为 `all`）；裸 `--sfs`（无窗口统计）不触发 all-sites 检查，levels 自动 genome-only | `a49d318` | `e2e_sfs`；手工：`reed stats gold/fixtures/sfs.vcf --pop … --sfs --proj 6` |
| 对拍 | `gold/popgen/run_easysfs.sh` + `compare_sfs.py`（rtol 1e-6，取整前） | `a49d318` | fixture：unfolded 7 bins 最大相对差 **1.04e-15**、folded **4.36e-16**；子集全量级（Chr01-03 × clean.3pop.txt 三群体，各投影到 2·min(n_i)=214）见下 |

**偏差**：`--sfs-int` 与 easySFS `--dtype int` 在“恰为半整数”的分箱可差 1（easySFS 的浮点噪声把精确 1.5 落成 1.499999999999999 → 取整 1；reed 的补偿求和得 1.5 → 取整 2）。规范里的 `1 3 3 2 1 1 0` 即此例；对拍以“取整前”的浮点 rtol 1e-6 为准（05 §3.1 的措辞“取整前”）。两群体（二维联合谱）与 `--sfs-int` 的半整数取整规则留待后续。

### X4.2 合并（r 群体）Fst：`--fst-global`（05 §3.2）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| `StatKind::FstGlobal` | 复用已实现的 `WcMultiAcc`（r 群体 WC84 方差分量，`m7_fst.rs` 已用独立参考验过 a/b/c）；`StatAcc::fst_global_columns` 给 `(value, no_snps)` | `50c1698` | `reed-stats/tests/accumulator_laws.rs`（proptest 加 `FstGlobal` 臂）；`crates/reed/tests/e2e_fst_global.rs`（3 项） |
| 输出 | `PREFIX_fst_global.txt`：`pops chromosome window_pos_1 window_pos_2 avg_wc_fst no_snps`，`pops` 为逗号连接的群体名；给 `--levels genome` 时多一行 `genome` | `50c1698` | 夹具与子集输出 |
| CLI | `--fst-global`（clap `requires = "fst"`）；`--fst hudson --fst-global` 报错 `hudson Fst has no r-population (global) form; use wc84`；两群体时提示不报错；`--levels` 在 pixy 布局下改为**透传用户值**（默认仍 window，行为不变） | `50c1698` | `cli_contract::help_lists_the_new_statistics_flags` |
| 对拍 | `gold/popgen/run_vcftools_fst.sh`：vcftools `--weir-fst-pop` ×3 的逐位点列与日志 weighted Fst | `50c1698` | **夹具**：逐位点最大相对差 1.0e-6（vcftools 只印 6 位有效数字，容差 2e-6）、genome 相对差 **2.2e-16**；**子集**（Chr01-03 × clean.3pop.txt 三群体，480,952 位点）：genome reed `0.0017395526794523797` vs vcftools `0.0017396` |

**偏差**：逐位点门取 2e-6 而非 05 §3.2 写的 1e-9——vcftools 的 `weir.fst` 以 `%g`（6 位有效数字）打印，参考文件本身的分辨率就是 5e-7；genome 级门 1e-9 达标（2.2e-16）。

### X4.3 Patterson f 统计（05 §3.3）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| genome-wide 请求路径 | `StatsRequest.fstats: Option<FStatsRequest { tests, block_bp }>`；`FStatTest::{F2,F3,F4,Dstat,F4Ratio}`；job 在 `push_row` 逐位点喂 `(num, den)`，partial 与窗口 partial 并列上交，sink 合并 | `2ca548f` | `crates/reed/tests/e2e_fstats.rs`（3 项） |
| 分块 jackknife | 新增 `BpBlockJackknife`：块键为 `(contig, pos0 / block_bp)`，与 chunk 划分和线程数无关（`--block-size`，默认 5M bp；Patterson/admixtools 口径）。原 `BlockJackknife`（按位点序数分块）保留 | `2ca548f` | 跨线程确定性测试；`jackknife.rs` 单测 |
| 公式复用 | `patterson.rs` 抽出 `f2_site`/`f3_site`/`f4_site` 三个逐位点函数，累加器改为调用它们（公式不再重复两份） | `2ca548f` | `patterson_fstats` 10 项全绿（抽取后未变） |
| CLI | `--f2 A,B`、`--f3 C;A,B`、`--f4 A,B;C,D`、`--dstat A,B;C,D`、`--f4-ratio A,B;C,D;X`、`--block-size SIZE`（可重复；分号分角色组，逗号分群体）；输出 `PREFIX.fstats`：`stat pops estimate se z n_blocks n_snps` | `2ca548f` | 夹具与子集输出；`cli_contract` 新旗标测试 |
| 对拍 | `gold/popgen/run_allel_fstats.sh` + `compare_fstats.py`（scikit-allel `patterson_f2/f3/patterson_d`） | `2ca548f` | **夹具**：f2 `0.375`、f3 `0.2916666666666667`、f4 `0`、D `0.5098039215686274`，最大相对差 **0**（门 1e-12）；**子集**（Chr01-03 × clean.3pop.txt）：f2/f3 最大相对差 **1.3e-16**（门 1e-9） |

**偏差**：SE 只在夹具级可比——allel 的 `blen` 按**变异数**分块，reed 按 **bp** 分块（Patterson/admixtools 口径），真实数据上两套块网格不同，故 `compare_fstats.py` 以 `--no-se` 比较估计量；SE 由 reed 的 bp 分块自洽给出（`n_blocks` 记录在输出里）。

### X4.4 LD：`--ld-decay` / `--indep-pairwise`（05 §3.4）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 流式逐对 r² | `LdRequest { decay_max_bp, decay_bin_bp, prune, pairs }`；worker 用 chunk 内 carry 算同 chunk 对；sink 用「本 chunk 的 **head** 窗 × 携带的 **tail** 窗」补全跨 chunk 对；carry 每 chunk 重置，任一对只计一次 | `f0d6d3f` | 区域 `Chr01:1-2,000,000` 在 1/4/56 线程下 md5 相同（427,153 对） |
| 修剪 | `LdPruner` 扩为 `LdWindow::{Variants,Bp}` + plink2 的 step 语义；sink 按位点顺序喂入全量位点（`all_sites`），顺序依赖的决策不进 worker | `f0d6d3f` | 夹具 prune 集合与 plink2 完全一致 |
| 输出 | `PREFIX.ld.decay`（`bin_start bin_end n_pairs mean_r2`）、`PREFIX.ld`（plink2 `.vcor` 列）、`PREFIX.prune.in/.out`（`CHROM_POS`） | `f0d6d3f` | 夹具与区域输出 |
| CLI | `--ld-decay [MAX_BP]`（默认 100k、1 kb 箱）、`--ld-pairs`、`--indep-pairwise W S R2`（W 为位点数或 `W kb`，kb 时要求 step 1） | `f0d6d3f` | `e2e_ld.rs`（2 项）；`cli_contract` 新旗标 |
| 对拍 | `gold/popgen/run_plink2_ld.sh`：plink2 `--r2-unphased` 与 `--indep-pairwise` | `f0d6d3f` | **夹具**：15 对 r² 最大相对差 1.02e-6（plink2 `%g` 6 位，容差 2e-6）、prune 集合相同；**区域**（Chr01:1-2,000,000，427,152 对）：无缺失对、最大相对差 5.0e-6（同为打印精度） |

**偏差**：r² 门取 2e-6 而非 05 §3.4 写的 1e-9——plink2 的 `.vcor` 以 `%g`（6 位有效数字）打印，参考文件分辨率即 5e-7；prune 集合门（集合相等）达标。`--indep-pairwise` 的 kb 窗口与位点窗口都已实现，但只在位点窗口上与 plink2 对拍（夹具 6 位点，kb 窗口无区分度）。

### X4 完成度

X4 的四项（3.1 SFS、3.2 `--fst-global`、3.3 f 统计、3.4 LD）全部落地，每项都有 fixture 逐值门与全量级（子集/区域）对拍，参照物分别为 easySFS、vcftools、scikit-allel、plink2。C-5（拆 `individual.rs`）未做，登记在下方偏差表。

### X4.4-rev `--indep-pairwise` 迁入 `reed filter`（05 §3.4.1，裁决 2026-10-10）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 归属与请求 | `FilterSpec.ld_prune: Option<LdPruneSpec>`（window/kb/step/r2），由共享 `SiteFilterArgs` 的 `--indep-pairwise W[kb] [S] R2` 解析；`reed stats` 对该旗标报错并指向 filter；`reed stats` 侧只留 `--ld-decay/--ld-pairs` | `33c8948` | 手工：`reed filter --indep-pairwise …` 输出剪枝后的 VCF；`reed stats --indep-pairwise …` 报错 |
| pass0 | `reed/src/filter/ld_prune.rs`：与 pass1 同一 `RecordPipeline` 存活判定（head→mask→tail），把每个 contig 的存活位点按序交给 sink 线程上的一个 `LdPruner`；结果是"每 contig 保留位点集合"，pass1 发射前二分查找 | `33c8948` | `e2e_ld::filter_indep_pairwise_keeps_the_plink2_set` |
| 算法 | `LdPruner` 严格照 `plink2_ld.cc`：到达空间的窗口（W 位点或 W bp）、step 推进、到达即剪单态（全 hom-ref / 全 hom-alt / 全杂合 / 无调用）、只对**新进入**的槽位反向最近优先扫描、被剪的 second 仍参与比较、主等位频率高者被剪（1e-6 内平局剪靠后）、contig 末尾的残余窗口也要扫 | `33c8948` | `reed-stats/tests/ld.rs`（7 项） |
| 漏斗 | `after_mask.ld_prune`；`FilterReport.ld_prune = (survivors, pruned)` | `33c8948` | `e2e_ld` 断言 `Some((6, 3))` |
| 对拍 | `gold/popgen/run_plink2_ld.sh`：r² 用 stats 的 `--ld-decay`，prune 集合用 filter 的输出 | `33c8948` | **夹具** `10 1 0.5`：保留位点 `{1,5,6}` == plink2 `prune.in`；**区域**（Chr01:1-2,000,000，8,364 位点）`50 5 0.2`：277 == 277，对称差 **0**；1/4/56 线程正文逐字节相同 |

**偏差**：bp 窗口（`--indep-pairwise 50kb 0.2`）尚未与 plink2 逐 ID 一致（区域内 reed 保留 131 vs plink2 158），窗口推进的 bp 语义仍需对齐；位点窗口（plink2 的主用法与 05 §3.4.1 的门）已逐 ID 相同。

## 未完成 / 偏差

| 项 | 状态 | 说明 |
|---|---|---|
| C-8（部分） | 偏差 | `reed-io` 被 reed-cli 源码直接使用（`Region`、`tag_table`、`bgzf`），未删；C-9 的 7 项里仅 `reed-io` 保留。已由 08 §6 的门面 re-export 覆盖（CLI 现经 `reed::` 使用） |
| 08 §6.3 前端含后端逻辑 | 未做 | `translate/stats.rs` 仍调 `reed::plan::compile` 并自行 `assign_groups`；`translate/filter.rs` 的 `site_spec` 编译同理。移入 `run_filter`/`run_stats` 会改 `run_*` 签名（`StatsRequest` 需携带 `FilterSpec`/pop 文件），留待后续 PR |
| 08 §6.2 单一 `reed::Error` + `kind()` | 未做 | 现有 `FilterError`/`EngineError`/`IoError`/`CoreError` 未合并为 `reed::Error`；退出码仍由前端各自映射。留待后续 PR |
| 08 §6.4.4 cli_contract 往返 | 未做 | 现有往返测试只覆盖 filter；stats/index 的往返测试未补 |
| L-2 完全精选 | 部分 | 门面以 `pub use` 重导出模块（`spec/expr/plan/funnel`、`budget/metrics`、`bgzf/chunk/header`）与类型；`unreachable_pub` 为 crate 级 allow。未逐项把非精选项改 `pub(crate)`（量大，留待后续） |
| C-5 拆 `individual.rs` | 未做 | X4 的四个统计量已落地，但 `individual.rs`（individual F + joint SFS 两个统计量）未拆为 `het.rs` + `sfs/{projected,joint}.rs`；属纯归属调整，留待后续 |
| 联合 SFS（两群体二维谱） | 未做 | 05 §3.1 的两群体二维谱（`JointSfsAcc` 已有实现）未接 CLI；当前 `--sfs` 写一维谱（每群体一个） |
| X5 多等位重算 | 待做 | 05 §4（`--fill-tags` 的逐等位 AC/AF、MAF 次等位口径、HWE 多等位写 `.`） |
| X6 CLI 缺口 | 待做 | `--bed-file`、`--sites-file`、`-q/-Q` 的 `:type` 后缀 |
| LD bp 窗口 | 偏差 | `--indep-pairwise Wkb` 的窗口推进未与 plink2 对齐（位点窗口已逐 ID 相同），见上表 |
| `reed stats --indep-pairwise` | 未做 | 规范 §3.4.1 说 filter/stats 共用；stats 侧需要第二遍，现为显式报错并指向 filter |
| X7 | 已落地 | `.github/workflows/ci.yml`（X2 时创建，数据门条件化） |

> 本轮收口于 green 主干：X0、X1、X2、X3、X7 已落地，未完成项逐条在上表登记。X4 按计划"每个统计量单独 PR"，从 05 §3.1 的 SFS genome-wide 路径开始。

## 基准（2026-10-09 v5 复审基线，56 线程，热缓存；负载与中位数见 `04-performance.md` §1）

（复审时的数字由 `04-performance.md` 给出；实现者在每步完成后按同一命令重测并补登 `bench/suites/e2e/results.tsv`。X0/X1 未触及计算路径，未重跑全量基准。）
