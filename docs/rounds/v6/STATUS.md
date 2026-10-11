# doc_v6 实施状态

> 对照 `06-roadmap.md` 的 Y0–Y7。每项注明 commit 与验证方法（测试名或可复现命令）；偏差单独登记。格式沿用 `docs/rounds/v5/STATUS.md`。

## Y0：lint 门变绿 + 漂移/死代码归零 + 交付性（2026-10-10）

### Y0-A lint 与文档门（E6-1、E6-7、E6-2 的 fmt/doc 两项）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| CI clippy 门 | 修 `reed-filter/src/eval/tail.rs` 的 `doc_lazy_continuation`（合并两段 doc，改列表为 `*` 项）；测试目标的约 60 条 pedantic 逐条修（`assert!`→`assert_eq!`、`()` 模式、cast 用 `try_from`、`write!`/`writeln!` 取代 `push_str(&format!())`、`memchr` 取代逐字节 `filter().count()`、数字字面量加分隔符）；三处必须保留的例外带 `#[allow]` + 理由（表达式/谓词的精确浮点比较、`--version` 构建脚本的 `map_or_else` 形态、worker-panic 测试的显式 `panic!` 消息契约） | `99cb913`、`dc8b709` | `cargo clippy --release --workspace --all-targets -- -D warnings` 退出 0 |
| rustdoc 门 | 修 10 处坏链/冗余链（`bitset.rs`、`count/site_counts.rs`、`count/mod.rs`、`gt_code.rs`、`individual.rs`、`ld.rs`、`error.rs`、`edit_plan.rs`、`plan/mod.rs`、`input.rs`、`job.rs`、`filter_flags.rs`） | `c3cc4c7`、`99cb913` | `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps` 退出 0（唯一 warning 是 bin/lib 同名产物的已知 cargo 问题） |
| fmt 门 | 全树 `cargo fmt`；**偏差**：`docs/rounds/v6/05 §7.3` 写的 "0 diff" 不成立——v6 复审用的 rustfmt 与仓库固定的 1.99 工具链不是同一个，pristine 树在 1.99 的 `cargo fmt --check` 下 158 个文件有 diff。本轮一并格式化 | `99cb913` 等 | `cargo fmt --all --check` 退出 0 |

### Y0-B 漂移清单（05 §4 八条，E6-7）

| # | 位置 | 修正 | commit |
|---|---|---|---|
| 1 | `reed/src/filter/ld_prune.rs:7-9` "contig 各自线程、并行剪枝" | 改为事实：一个 `LdPruner` 在有序 reducer 里串行，66 s 的代价见 04 §6，contig 并行是 Y1 | `c3cc4c7` |
| 2 | `docs/rounds/v5/STATUS.md` X7 行 "本地五个可移植门等价命令全绿" | 划掉并加更正说明（clippy 门当时是红的，用的不是 CI 命令） | `2262ee1` |
| 3 | `reed-io/src/vcf/line_edit.rs:8-10` "无编辑的记录不该走 `edit_vcf_line`" | 改为事实：两处调用点无条件调用（F6-11），逐字复制快速路径排在 Y3 | `c3cc4c7` |
| 4 | `reed-io/src/error.rs:2` 引用 `crate::ReadError` | 改为纯文本 `ReadError`（类型已不存在） | `c3cc4c7` |
| 5 | `reed-core/src/tag_table.rs:91-105` "按声明顺序" | 改为 "sorted by name"（实现是 `names.sort()`） | `c3cc4c7` |
| 6 | `reed-stats/src/individual.rs:44-47` "无逐样本循环" | 改为事实：逐样本循环，词并行改写排在 Y3（F6-11） | `c3cc4c7` |
| 7 | `reed-core/src/count/kernel.rs` SAFETY "dispatcher 断言长度" | 文件已删（A6-22） | `c3cc4c7` |
| 8 | `docs/spec/architecture.md:804`、`reed-engine/src/budget.rs:46-50` "N ≤ 2 走顺序模式" | 两处都注明现状：`run()` 无条件并行，`is_sequential` 无调用者，Y5 裁决 | `c3cc4c7`、`2262ee1` |

### Y0-C 死代码与 unsafe（E6-9、A6-22/23/24/25）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 删除清单 | `reed-core/src/count/kernel.rs`（AVX-512，release 下长度检查不成立 = 可达 UB；348 样本 W=6<8 永不向量化；无生产调用者）、`CountTable`、`math/hypergeom.rs`（与 `reed-stats::sfs::projection_weights` 一问两答且无调用者）、`reed-io/src/field_plan.rs`（无消费者）、`ChunkBytesExport`、`MemBudget::from_half_memory_for_tests` | `c3cc4c7` | `cargo build --release --workspace --all-targets` 0 警告；测试数 382→373（删 3+5+2=10 项、`add_code` 新增 1 项，逐项核对） |
| `CountTable` 删除后的覆盖 | `SiteGroupCounts::add_code` 是生产路径（`stats/job.rs`、`writers/reports.rs` 调用），补 4 项单测（手算、半调用、`code_has_alt`、全 REF） | `c3cc4c7` | `reed-core/src/count/mod.rs` 内联测试 |
| unsafe 守卫 | `allow(unsafe_code)` 从 2 处降到 1 处（只剩 mmap）；`ci/check-doc-drift.sh` 的断言改 `-eq 1`；`docs/spec/architecture.md` §14.4 同步；`bgzf.rs` 的 SAFETY 注释改为 "single audited site" | `c3cc4c7` | 守卫 `ok exactly one allow(unsafe_code) site` |
| `ContigId` 统一（A6-23） | 删 `reed-core/src/contig.rs` 的裸 `u32` 别名；全工作区用 `id.rs` 的 newtype；`reed/src/index.rs` 的 `tid`、`record.key.contig` 的消费点改用 `.get()` | `c3cc4c7` | `cargo check --workspace --all-targets`；gold 对拍与 md5 不变（见 Y0-F） |
| `is_sequential`（A6-25） | 保留（Y5 实现顺序模式或改规范），文档按事实改写 | `c3cc4c7` | — |

### Y0-D 命名与结构（E6-13、E6-14、E6-16、E6-17）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| pixy 词汇出生产代码 | `writers/pixy.rs`→`writers/windows.rs`、`PixyOptions`→`WindowLayoutOptions`、`pixy_kinds`→`window_kinds`、`StatsRequest.pixy`→`layout`、`reed::pixy::write`→`reed::windows::write`；生产 `src/` 的注释也不再出现该词（引用写在文档/测试/gold 脚本里）。测试目录与冻结 oracle 保留参照物名（naming §2.5） | `99cb913`、`dc8b709` | 守卫 `ok no pixy vocabulary in production src`；报告 md5 不变（Y0-F） |
| 测试文件改名 | `m7_fst.rs`→`fst_multi.rs`、`qc_merge_laws.rs`→`summary_merge_laws.rs`、`individual_joint.rs`→`het.rs`+`sfs_joint.rs` | `dc8b709` | 同名测试全绿（6/3/2/1 项） |
| 模块风格统一（T-1） | `reed-cli/src/cli.rs`→`cli/mod.rs`；全树无 `foo.rs`+`foo/` 并存；`docs/spec/naming.md` §2.9 写明约定 | `99cb913`、`2262ee1` | 守卫 `ok no sibling foo.rs + foo/ module pairs` |
| crate 根统一（T-2） | `FilterError` 从 `reed-filter/src/lib.rs` 搬入 `reed-filter/src/error.rs`（四个产品 crate 都有 `error.rs`）；`reed-engine`/`reed-filter` 根加 2–4 个中心类型 `pub use` | `dc8b709` | 守卫 `ok the four product crates have error.rs` |
| `#[doc(hidden)]` 过期（T-3） | `reed-stats/src/lib.rs` 去掉已接 CLI 项的隐藏（`WcMultiAcc`、`ProjectedSfsAcc`、`LdPruner`、`F2Acc…`、`BpBlockJackknife` 等）；只留 `JointSfsAcc`（无 CLI 出口，Y6） | `dc8b709` | `cargo doc -p reed-stats` 索引含 `WcMultiAcc`/`ProjectedSfsAcc`/`LdPruner`/`F2Acc` |
| `ponytail` 标记 | 3 处改为 `TODO` 文本 | `dc8b709` | 守卫 `ok no ponytail markers remain` |
| 公开参数去缩写 | `site_header(rec)`/`decode_record_gt(rec)`→`record`、`deflate_block(src)`→`input`、`bin_range(idx)`/`count_at(idx)`→`index` | `99cb913` | `cargo doc -D warnings` 通过 |

### Y0-E 交付性（E6-10、E6-11、E6-12）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 许可证文本 | 新增 `LICENSE-MIT`（MIT）与 `LICENSE-APACHE`（Apache-2.0 全文）；README 两语言指向两份文本 | `2262ee1` | 守卫 `ok LICENSE-MIT and LICENSE-APACHE exist`；`ls` |
| 构建信息 | `reed-cli/build.rs` 注入 `git describe --always --dirty`、提交日期、`rustc --version`；`--version` 输出 `reed 0.3.0 (64e5dce-dirty 2026-10-10, rustc 1.99.0 (b940084d7 2026-09-28))` | `99cb913` | 运行二进制；守卫 `ok reed-cli has a build.rs for --version provenance` |
| CHANGELOG | `CHANGELOG.md`（Keep a Changelog；Unreleased + 0.3.0 + 0.2.0 三节，从各轮 STATUS 回填） | `2262ee1` | 文件存在 |
| README 四节 | Installation（`cargo install --path crates/reed-cli`）、Input requirements（VCF/BCF/stdin、`##contig` 长度、`-r/-R` 需索引、contig 排序）、Exit codes（0 / 1）、Citation（Nei & Li 1979、Nei 1987、Watterson 1975、Tajima 1989、Hudson 1992、Bhatia 2013、Weir & Cockerham 1984、Patterson 2012、Korunes & Samuk 2021）；中英同步 | `2262ee1` | 人工；README 门面守卫仍绿 |
| `v0.3.0` tag | 打在 Y0 全部提交之后（见本文件末"tag"行） | — | `git tag -l`；`git describe` |

### Y0-F 门与不变性（本轮实测）

| 门 | 命令 | 结果 |
|---|---|---|
| 测试 | `cargo test --release --workspace` | **373 通过 / 0 失败**（52 个目标；382 − 10 删除 + 1 新增） |
| 构建 | `cargo build --release --workspace --all-targets` | 0 警告 / 0 错误 |
| clippy（CI 命令） | `cargo clippy --release --workspace --all-targets -- -D warnings` | 退出 0 |
| rustdoc | `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps` | 退出 0 |
| fmt | `cargo fmt --all --check` | 退出 0 |
| 守卫 | `bash ci/check-doc-drift.sh` | **29/29 ok**（原 23 + Y0 新增 6） |
| gold 正文 md5 | `/tmp/reed-target/release/reed filter clean.vcf.gz <gold 参数>`（`gold/scripts/run_reed.sh`）→ `bcftools view -H \| md5sum` | **`7aedc59ac798729a2674e6e7ef19adfc`**（与 v3–v6 记录逐字符相同） |
| gold 对拍 | `python3 gold/scripts/check_parity.py /tmp/gold_new/reed.vcf.gz --gold target/gold/out/full/gold.vcf.gz` | **PASS**：1,768,759 / 1,768,759 位点，Jaccard 1.000000，348 样本 MATCH，基因型差 mean=0 max=0 |
| 子集数字门 | `reed filter subset.Chr01-03.vcf.gz --min-dp 5 --min-gq 20 --geno 0.1`（与 `--geno 0.2`） | **263,420** / **381,672**（== 三家基线） |
| 报告 md5 | `reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc` | 8 个文件逐一等于 v5 记录：afreq `b2884092…`、hardy `fd8fbf27…`、het `790293d9…`、smiss `5f4fe4f8…`、vmiss `7c325d64…`、idepth `a4b6dd3a…`、ldepth.mean `475ee083…`、stats `189e025b…` |

### Y0-G oracle 引用词

冻结 oracle（`crates/oracle`，仅测试用）的 6 处 `pixy` 引用改为论文引用（Korunes & Samuk 2021/2023），
使 roadmap 的**字面**门 `! grep -rniE "\bpixy" crates/*/src --include=*.rs` 通过（无目录豁免）；守卫同步为字面形式。

新增守卫（`ci/check-doc-drift.sh`）：unsafe 计数 2→1、无 `pixy` 词汇（生产 src）、无 `ponytail`、无 sibling 模块对、四个产品 crate 有 `error.rs`、`build.rs` 构建信息、两份 LICENSE 文本。

**未纳入 Y0-F 的门（留 Y7 或需数据 runner）**：跨线程 sha256、跨格式 md5、BCF 输入全量、`cargo deny`、MSRV 构建、基准回归。

### Y0 偏差登记

| 项 | 说明 |
|---|---|
| fmt 门的既有偏差 | v6 05 §7.3 的 "`cargo fmt --all --check` 0 diff" 是在另一个 rustfmt 版本下得出的；固定工具链 1.99 下 pristine 树有 158 文件 diff。Y0 已全树格式化（`99cb913` 等），此后 fmt 门可信 |
| `reed stats --indep-pairwise` | v5 遗留：stats 侧显式报错并指向 filter；规范 §3.4.1 的"filter/stats 共用"仍待 Y6 |
| `v0.3.0` tag | 见下方"tag"行；tag 指向 Y0 完成的提交，其后的 Y1+ 提交在 tag 之外 |

## Y1：全局分析矩阵化（进行中，2026-10-11）

### Y1 已落地项

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| F6-9（Y1a） | `--ld-decay`/`--ld-pairs` 不再要求 `--pop`（无 `--pop` 时全体样本计为单群体 `all`；`--pop` 仍是可选项） | `d2f579e` | 夹具 + 区域命令退出 0 |
| 剪枝 contig 并行（Y1c） | reducer 只按 contig 收集存活位点，`rayon` 每 contig 一个 `LdPruner`（plink2 的 subcontig 并行；contig 内严格 plink2 顺序语义） | `9bfffa1` | 全量 `--indep-pairwise 50 5 0.2 -Oz1` 235,403 保留 == plink2；1/8/56 线程正文逐字节相同 |
| F6-10 kb 窗口（Y1g） | `LdWindow::Bp` 按 `plink2_ld.cc` `LdPruneNextWindow` 对齐：窗口在触发位点加入前闭合，start 推进到 `bp ≥ bp[触发] − W`；tie 判据用 plink2 的 `kSmallEpsilon`（2⁻⁴⁴）；并修 CLI 单位 bug（`50kb` 此前被读成 50 bp） | `a5cd8a2` | 区域 `50kb 0.2`：158 == plink2 对称差 0；`50 5 0.2`：277 == 277 对称差 0；全量 235,403（4.0 s，pass0_prune 147 ms）；`reed-stats` 9 项单测（2 项 plink2 实测 kb 用例） |
| LD 矩阵化 + 位置区块并行（Y1b/Y1d） | stats 的 LD 不再走 head/tail 携带：worker 按 contig 填站点主序位平面矩阵（样本按群体主序、字对齐——群体对读连续字切片，cohort 和 = 群体和的整数精确相加）；每 contig 矩阵在独立线程做位置区块（4096）并行衰减，`finish` 按 contig 序合并 | `35444fd` | 全量 `--ld-decay 100k`：**8.4 s**（`--pop` 3 群体）/ **6.4 s**（无 `--pop`），门 ≤ 10 s（原 169 s）；cohort 衰减与串行实现逐位相同（区域 3,556,592 对、100/100 bin 全同；夹具 r² PASS）；1/8/56 线程 `.ld.decay`/`.ld` md5 相同（区域 + 全量）；BCF 输入与 VCF 逐对相同（顺带修了旧实现 BCF 路径 contig 不复位的缺陷） |
| 每群体衰减曲线 | `.ld.decay` 增加首列 `pop`：cohort 记 `all`，`--pop` 命名多于一个群体时逐群体各一块；群体曲线只计在该群体样本上有定义的对（对按样本子集计，不按位点） | `35444fd` | `e2e_ld::ld_decay_carries_per_population_curves`（群体计数 ≤ cohort 计数、cohort bin = 15 对）；`.ld.decay` 400 行 = 100 bin ×（all + 3 群体） |

### Y1e/Y1f/Y1h：统一矩阵 pass0（2026-10-11）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 统一 MatrixJob/MatrixSink（Y1e，03 §7） | `reed/src/filter/pass0/{mod,matrix,prune}.rs`：一个 chunk job 对每个 head 幸存记录产 `(4+R)·w` 字的位点行（原始 4 平面 + 每掩码规则一张触发位图），sink 按 contig 拼 `SiteMatrix`（块式，步长校验）并记逐 chunk 的记录位图与全局位点号；`MissingJob`/`SketchJob`/`LdPruneJob` 与 `reed_filter::sketch` 全删，`mid_stage` 改从矩阵按全局位点号读回（计数函数按原始平面 + 规则位图移植，`masked_gt` 等公式逐位保留）；`scan_head_failures`/`RecordView::head_passed` 死面随 sketch 一起删 | `2863b61` | `cargo test --release --workspace` 全绿（移植的 4 项平面计数单测 + 3 项新 e2e）；三策略全基因组正文 md5 一致；子集 1/8/56 线程一致 |
| `--mind` 从 miss 平面按列 popcount | worker 填完位点平面后对 `(miss ∨ mask)` 逐位提取累计 `m0/m1`，单元收口 `called = 2·sites − m0 − m1`（doc_v4 02 回退版 §3.2：分母是 head 幸存者，掩码计缺失）；不再有第二条计数路径 | `2863b61` | 子集三策略 + t1/t56 正文 md5 全等（`fcc6494c…`）；全基因组三策略正文 md5 全等；`mind_counts_post_mask_over_head_survivors` 三策略分支全绿 |
| head 幸存者 + tail 标志语义 | 矩阵收 **head** 幸存者；tail 不再在 pass0 内联（有 `--mind` 时 keep 未定），剪枝用的幸存过滤来自 mid stage 的 per-site tail 判定；无 `--mind` 时 tail 在 pass0 内联（keep 已定），矩阵只收 head∧tail 位点——剪枝序列 = 每条位点规则之后（02 §6 step 6b，plink2 载入后剪枝的顺序） | `2863b61` | `pruner_sees_only_post_tail_sites`：`--geno 0.5` 下 pruner n_seen 3→2、输出 {1} vs 无 tail {1,2}；`mind_counts_ignore_tail_rules`：加 tail 规则后 `--mind` 落样数不变（post-tail 分母会把边界样本翻转） |
| `--mind` + `--indep-pairwise` 组合 | 同一 pass0：矩阵带规则位图，mind 解析 → mid stage（tail over 最终 keep）→ 幸存位图过滤剪枝 → EmitJob 携 `site_keep` 跳过被剪位点；漏斗 `n_output` 扣除被剪数、`n_ld_pruned` 记被剪数（旧代码该组合会因二次 `first.take()` panic） | `2863b61` | `mind_and_indep_pairwise_combine`：9 样本（mind 落 1）+ 每组剪至首位点，输出 [1,6] |
| 流式 contig 剪枝 + 预算（Y1f） | 无 `--mind` 时每 contig 关闭即 `std::thread` 交付剪枝（任一时刻常驻一个矩阵，剪枝与解码重叠）；单 contig 矩阵超 `--memory` 字节预算报 `MemBudgetExceeded`；sketch 形态（整基因组常驻）的预算判断沿用 `choose_for_mind`（位点 × (4+R)·w·8，与 sketch 同式） | `2863b61` | 区域剪枝与 plink2 对称差 0（见下）；全量 235,403 == plink2（与 pre-Y1e 二进制正文 md5 一致） |
| gold 门复跑（Y1h） | gold 正文 md5、子集数字门、报告 md5、plink2 夹具 LD 对拍、区域剪枝对拍全部重测 | `2863b61` | 见下"Y1e 门"表 |

### Y1e 门（2026-10-11 实测，除注明外 56 线程）

| 门 | 结果 |
|---|---|
| 全量 `filter --indep-pairwise 50 5 0.2 -Oz1` | 235,403 保留，正文 md5 == pre-Y1e 二进制（其 plink2 等价在 Y1c 验证）；**3.6 s**（Y1c 基线 4.0 s）；t1 正文 md5 相同 |
| 区域 `50kb 0.2` / `50 5 0.2`（plink2 `p2kb`/`p2` 参考） | 158 == 158、277 == 277，对称差 0 |
| gold 正文 md5（`run_reed.sh`，全量） | **`7aedc59ac798729a2674e6e7ef19adfc`**（与 v3–v6 逐字符相同） |
| 子集数字门 `--min-dp 5 --min-gq 20 --geno 0.1 / 0.2` | 263,420 / 381,672（== 三家基线） |
| 报告 md5（8 文件，全量） | afreq/hardy/het/smiss/vmiss/idepth/ldepth.mean/stats 与 Y0 记录逐字符相同 |
| `gold/popgen/run_plink2_ld.sh` | PASS r²（15 对，1.02e-06）+ PASS prune（3 保留/3 剪除） |
| `--mind` 子集三策略 + t1/t56 | 五份正文 md5 全等（`fcc6494c…`） |
| `--mind` 全基因组三策略 | 三份正文 md5 全等；**6.5/6.5/6.6 s**（single_pass/two_pass/sketch 形态） |
| 测试 | `cargo test --release --workspace` 379 通过 / 0 失败 |
| fmt / clippy（CI 命令）/ rustdoc / 守卫 | 全部退出 0；守卫 29/29 ok |

### Y1 剩余

（无——Y1e/f/h 完成。）

## Y2：I/O 热路径（2026-10-11）

### Y2 已落地项

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| A6-5 deflate 常驻缓冲 | `Deflater` 带常驻 scratch（`len` 恒等于 capacity，只随 bound 增长一次），块输出改为 header/BSIZE/负载/CRC 四段 `extend_from_slice`——每块不再 `resize(bound+8, 0)`（gold 输入 ≈350K 块的 65 KiB memset + realloc 归零） | `0493575` | roundtrip/noodles crosscheck 单测不变；gold 正文 md5 不变 |
| A6-6 inflate 按元数据走块 | `Inflater::inflate_window(block, 精确 ISIZE 窗口)` 免清零（窗口可为上一次用过的旧字节，libdeflate 恰好覆写它报告的量）；`inflate_chunk` 改读 chunk 的 `blocks` 元数据——免二次块头解析、免每块 64 KiB memset，worker 的 scratch 只随最大 chunk 增长一次；校验 Σcomp_len == payload 长度 | `0493575` | 引擎 9 项 pipeline 测试（夹具补上真实 blocks 表）；全部 inflate 单测 |
| A6-7 流式源偏移漂移 | `fill_to_target` 不再 `drain` 空成员：保留在 chunk payload、`blocks` 记 `(size, 0)`、chunk 目标按非空字节计（与 mmap 源同规则）；`carry_offset = compressed_start + cut` 从此等于真实文件偏移（旧实现里被删空成员的字节使后续 chunk 的 `compressed_start` 偏小，stdin/管道输入的虚拟偏移错） | `0493575` | 新增 `streaming_and_mmap_sources_agree_on_chunk_metadata`（空成员夹具：两源逐 chunk 的 offset/blocks 表/first_in_block_offset/payload 字节全等，夹具含 ≥2 个空成员条目）；01-review D2 的两个拼接夹具全绿 |
| A6-8 memrchr | `owned_body_range`/`tail_range` 的标量 `rposition` → `memchr::memrchr`（每 chunk 全文扫描） | `0493575` | seam/记录完整性既有测试 |
| A6-10 GtTable 进程级共享 | filter 侧 `FilterJob`/`MatrixJob` 持 `Arc<GtTable>`：每 run 一份（`native_bcf` 时构造一次），worker 的 `NativeScratch` 只持 Arc 克隆；stats 侧本就每 run 一份 | `a47e43f` | 380 测试；BCF 输入 gold 正文 md5 == VCF |
| A6-11 每记录分配 | `render_info_value` 直用 `format_g6_into`（去掉 String 中转）；`edit_info` 直写 `out`（seen 上栈，>8 条回退堆，结果空才补 `.`）；`parse_pos1`/`parse_qual` 收借用的 contig 字节/`(chrom, pos1)` 元组、`Locus` 只在错误分支构造（site.rs/head.rs/expr/summary.rs/index.rs 的成功路径不再每记录 `from_utf8_lossy().into_owned()`）；`decode_gt_text`/`decode_pair`/`decode_slow`/render.rs `bad()` 泛化到 `LocusSource`，`decode_record_gt` 传 `(contig_name, pos1)`；`render_site_part`/`render_sample_part` 直接用字典里的 contig 字节 | `a47e43f` | g6 gold 渲染测试、INFO/样本编辑测试、BCF 渲染测试全绿；跨格式 md5 门 |
| A6-12 索引写出 | `Indexer::add(name: &[u8])`：CSI 路径（默认）从不构造名字，TBI 才 lossy 转换；`OrderedOutput` 与 `reed index` 的 reducer 侧不再每记录造 String | `a47e43f` | gold 运行含 `-W`，其 `.csi` 照常写出、gold md5 与 parity 不变 |

### Y2 门（2026-10-11 实测，除注明外 56 线程）

| 门 | 结果 |
|---|---|
| gold 正文 md5（`run_reed.sh`，8 线程） | **`7aedc59ac798729a2674e6e7ef19adfc`**（与 v3–v6 逐字符相同）；`check_parity` PASS：348 样本 MATCH，基因型差 mean=0 max=0 |
| 跨格式 | BCF 输入 gold 正文 md5 同 `7aedc59a…`（two_pass_reparse，8.9 s） |
| 子集数字门 | 263,420 / 381,672（== 三家基线） |
| `--mind` 三策略 + t1/t56 | 五份正文 md5 全等（`fcc6494c…`，== Y1e 记录） |
| filter 跨线程 | 子集 t1/t56 正文 md5 相同（`fef4cb07…`） |
| 报告 md5（8 文件，全量） | afreq `b2884092…`、hardy `fd8fbf27…`、het `790293d9…`、smiss `5f4fe4f8…`、vmiss `7c325d64…`、idepth `a4b6dd3a…`、ldepth.mean `475ee083…`、stats `189e025b…` == Y0 记录 |
| stdin 等价 | 子集 8/8、全量 8/8（`cat clean.vcf.gz \| reed stats …` 与文件路径的 8 份输出 md5 相同） |
| plink2 LD 夹具（`gold/popgen/run_plink2_ld.sh`） | PASS r²（15 对，1.02e-06）+ PASS prune（3 保留/3 剪除） |
| filter gold 计时 ×3 | **user 400.2 / 399.5 / 400.0 s（门 ≤ 420，baseline 460 → −60 CPU-s）**；wall 7.70 / 7.64 / 7.67 s（baseline 8.85） |
| stats 五统计 ×3 | user 70.8 / 70.6 / 70.9 s（门 ≤ 68 **未达**，见偏差登记）；wall 1.73 / 1.77 / 1.74 s |
| reports 极差（门 ≤ 3%） | /tmp 3.10 / 3.13 / 3.13 s；tmpfs（/dev/shm）3.09 / 3.10 / 3.11 s——极差 ≤ 1.3% |
| 四件套 | 测试 380 通过 / 0 失败；fmt、clippy（CI 命令）、rustdoc 退出 0；守卫 29/29 ok |

### Y2 偏差登记

| 项 | 说明 |
|---|---|
| stats ≤ 68 s 未达 | 实测 70.6–70.9（baseline 72，−1.4）。04 §3 对 A6-6 的 3–5 CPU-s 估计偏高：旧实现的每块 memset（64 KiB）与实际正文同量（bgzip 满块 ≈ 64 KiB inflated），常驻缓冲真正省掉的是"每 worker 一次性"之外的 ≈15 GB memset ≈ 1.5 CPU-s，加 memrchr 与分配清理后与实测 −1.4 一致；stats 路径没有 deflate（A6-5 的大头在 filter），HWE 缓存（A6-16）在 Y3 |
| reports 双峰（04 §7） | **不复现**：同机 /tmp 与 tmpfs 各 3 次，wall 极差 ≤ 1.3%；`--verbosity 2` 的 reader_ms ≈ 360、reduce_ms ≈ 2.0 s 三次稳定。输出设备无关 → v6 复审时的写回嫌疑排除，双峰归因于当时批次环境（复审在 1.28 负载下测得 3.0–5.4 s）。门按现状达成 |

## Y3：过滤热路径（2026-10-11）

### Y3 已落地项

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| A6-20 | `emit_survivor`：文本路径与原生 BCF 路径的两套"编辑 + 写出"闭包合并为一个方法（site_keep 检查 → fill_info_set → EditPlan → edit_vcf_line → IndexKey → 写出；BCF 路径传入渲染后的整行 `edit_line`） | `731ccf4` | 子集正文 md5 三门不变（见 Y3 门表）；`cargo clippy` 通过 |
| A6-15 | `stale: Vec<&[u8]>` 每 chunk 构建一次（`stale_refs`，run_chunk/run_seam/run_bcf_records/emit 的 run_chunk），不再每幸存记录一次；`filtered.rs` 本就用空表未动。**偏差**：05 写的"预先算在 FilterPlan 里"不可行——`EditPlan.info_remove: &[&[u8]]` 的借用形状从 plan 的 `Vec<Vec<u8>>` 拿不到 `'static`，且 `-x` 可传任意 tag（`translate/filter.rs:301`），每 chunk 构建已消除每记录分配 | `731ccf4` | 同上 |
| A6-16 | HWE 走 `HweCache`：`SiteCounts::hwe_test_input()` 暴露 `(n_called, mac, het)`（degenerate 返回 None，哨兵 p=1 语义逐位保留）；`fill_info_set` 加 cache 参数（worker 级，`--fill-tags HWE`）、`eval_tail` 的 MinHweP（pipeline 字段）、`--hardy` 报告（`ReportScratch` 字段）、mid_stage 每 chunk 一个（决策并行，无锁） | `731ccf4` | gold 正文 md5 与 reports 8 文件 md5 全部不变（cache 命中与未命中返回同一位型）；site_counts/tail 既有单测不变 |
| A6-19 | `Bitset::clear()`：整字 `fill(0)`（padding 保持零不变式）；run_stages 不再逐位 `set(false)` | `731ccf4` | `--mind`/单 pass 子集 md5 不变 |
| A6-17 | 列式掩码：每条 field 规则对整列（连续 u32/f32）一次遍历，首中归属经 `attribute_hit`（bit + mask_rule + mask_hits 语义逐位保留；表达式规则仍先于 field 规则，per-sample 顺序不变）；`eval_mask`/`mask_values` 删除；`--set-gts .` 复用管道的 `(tag, value)` 表原地更新（`FmtValue` 类型与闲置的 `fmt_values` scratch 删除，`SiteValues.fmt` 直接借用该表）——表达式路径无每样本分配；Integer 阈值走 `Op::apply_i64` 整数域比较（u32→i64 精确，序同构） | `a55c560`、`91a59df` | `--mind`（矩阵路径共用 run_stages）子集 md5 `fcc6494c…` == pre-Y3 二进制；set-gts 冒烟 `975e01ca…` == pre-Y3 |
| A6-18 | `CmpOp`/`LogicOp` 解析期一次成型（eval/compare/check 不再热路径匹配算符串 `"&&"`/`"==`）；`InfoIndex`（reed-io::vcf::info，片段字节区间 + round-trip 单测）每记录建一次：head 有任一 INFO/Expr 规则时在 `eval_head` 构建，类型化 INFO 规则与表达式（head 表达式、`--set-gts .` 掩码、tail 表达式）共用；`HeadScratch` 字段 pub 供 reed 侧复用 | `d3ece77` | 子集 gold-args 门（4 tag 表达式 + 8 tag stale）逐字节不变 |

### Y3 门（2026-10-11 实测）

| 门 | 结果 |
|---|---|
| gold 正文 md5（`run_reed.sh` 参数，直接调 reed，56 线程 ×3） | **`7aedc59ac798729a2674e6e7ef19adfc`**（与 v3–v6 逐字符相同，3/3） |
| `check_parity` | **PASS**：1,768,759 / 1,768,759 位点，Jaccard 1.000000，348 样本 MATCH，基因型差 mean=0 max=0 |
| 子集数字门 `--min-dp 5 --min-gq 20 --geno 0.1 / 0.2` | **263,420 / 381,672**（正文 md5 `fef4cb07…` / `371d3a0d…`，== Y2 记录） |
| `stats_site_filters_match_the_filter_then_stats_composition` | 通过（381 项内） |
| 报告 md5（8 文件，全量） | afreq `b2884092…`、hardy `fd8fbf27…`、het `790293d9…`、smiss `5f4fe4f8…`、vmiss `7c325d64…`、idepth `a4b6dd3a…`、ldepth.mean `475ee083…`、stats `189e025b…` == Y0 记录 |
| `--mind 0.1` 子集 | 正文 md5 `fcc6494c…`（== pre-Y3 二进制 == Y1e 记录） |
| `--set-gts .` 基准 | 无每样本分配（`fmt_storage` 原地更新，代码事实）；子集 user 94.7 vs 94.7（pre/post 一致，无回退） |
| filter gold 计时 ×3（同协议 A/B） | **user 395.4 / 395.5 / 395.7 s**（pre-Y3 同协议 398.2 → **−2.6 CPU-s**）；wall 7.62–7.68 s（pre 7.65） |
| 四件套 | 测试 **381** 通过 / 0 失败（+1：`index_get_matches_the_scan`）；fmt、clippy（CI 命令）、rustdoc 退出 0；守卫 29/29 ok |

### Y3 偏差登记

| 项 | 说明 |
|---|---|
| **filter user −30 s 门未达** | 实测 −2.6（398.2 → 395.6，同协议同机 A/B）。逐项消融（子集，261,788 幸存者）：`--fill-tags HWE` 全路径仅 **0.9 CPU-s**（pre-Y3 with-HWE 57.6 vs no-HWE 56.6），全量 ≈5–6，远小于门假设的份额；且 `(n_called, mac)` 工作集超过 8192 memo，clear-on-full 下命中 ≈0（位型不变，数据更均匀时才见效）。掩码的 u32→f64 转换不是瓶颈：整数域比较 395.7 vs 395.4（噪声内），循环受 `masked.get`/`evaluated` 的逐元素分支链与列扫描限制（无 perf 可用的机器上不再猜）。分配类（stale/emit）≈1 CPU-s（no-edit 消融 −1.0）。03/04 对 A6-15…18 的算术在本数据上系统性高估——与 Y2 stats 门（承诺 −3 实测 −1.4）同一模式。下一个大头是 Y4 记录模型（04：`-Ob` 文本往返 = +61 CPU-s），不在 Y3 微项里 |
| A6-15 落点 | 见上表行内说明：per-chunk 构建替代"算在 FilterPlan"，理由是借用形状 + `-x` 任意 tag |
| A6-18 落点 | "预编译为闭包/字节码"实现为解析期 op 枚举（热路径不再匹配算符串）；未做闭包编译（YAGNI） |
| 计时协议注 | `run_reed.sh` 的 `THREADS` 默认 8（config）；与 Y2 基线可比的计时 = 直接调 reed（默认 56 线程），8 线程跑出 wall 66 s / 6× 并行属协议差异，不是回退。正文 md5 与线程数无关（3/3 一致） |

## Y4：`-Ob`/BCF 记录模型（2026-10-11）

### Y4 已落地项

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| A6-9 文本→BCF 编辑内联（03 §8 的记录模型落地形态） | `NativeEncoder::encode_text(line, site, n_samples_in, edit, out)`：INFO 编辑与 GT 掩码/落样在**编码的同一次扫描**里应用——`EditPlan` 的 `info_remove`/`info_set`/`keep_samples`/`mask_gt` 语义逐条对齐 `edit_vcf_line`（stale 丢弃、命中 tag 原位替换、其余按 `info_set` 顺序追加、被掩样本 `./.`、落样跳过）；被替换值走类型化编码（`encode_info_typed`：Integer→`bcf_enc_vint` 宽度、Float 直写原 f32 位、String 走 `render_info_value`），不再"渲染→重解析"。`single_pass.rs` 的 `emit_survivor` 拆成 `EmitSource::{Text, Bcf}` 两路：native `-Ob` 直接 `encode_text`/`encode_bcf`，文本输出与 noodles oracle 仍走 `edit_vcf_line` | 本工作树（待提交） | `encode_text_applies_the_edit_like_the_text_path`（单测：与"edit+noodles"逐字节相等）；`native_and_noodles_encoders_agree`（默认 edit 下 B2/B4 逐字节相等，oracle 门保留） |
| A6-9 BCF→BCF 记录模型 | `NativeEncoder::encode_bcf(rec, dict, edit, out)`：字典下标经输入 `BcfDict` 翻译到输出头（contig/FILTER/INFO/FORMAT key）；未修改的 FORMAT 列**原字节拷贝**（落样按 stride 跳片），只重写 INFO 与 GT 列（GT 从源码直拷，掩码样本写两个 missing、落样跳过，宽度按最大等位重选）；ID/等位原字节拷贝。BCF→任意输出（`-Ob` 与 `-Oz`）都走它：`-Oz` 时不再 `render_record` 整行文本 | 本工作树 | `encode_bcf_matches_the_text_reference`（单测：与"edit+noodles"逐字节相等）；子集 BCF 输入 `-Oz` 正文 md5 == VCF 输入（`fef4cb07…`） |
| A6-14 相位与 `-0.0` | ① `fmt_gt` 增 `vcf44` 形参：头声明 ≥ 4.4 时按 htslib `bcf_format_gt_v2`（`vcf.c:6495-6517`）渲染首等位的前缀 `|`/`/`（仅在推断会错时写），`BcfDict` 解析 `##fileformat` 得到该标志；② `decode_gt` 接受 `|`/`/` 前缀（首等位显式相位）；③ `format_g6_into` 保留 `-0`（htslib `kputd` 语义） | 本工作树 | `fmt_gt_prefixed_phasing_matrix`（11 例码→文本矩阵）；`decode_gt_accepts_the_phasing_prefix`；`g6_matches_htslib_gold_rendings` 增 `(-0.0, "-0")`；`vcf44_phasing_and_negative_zero_survive_the_round_trip`（e2e 夹具，`-Oz` 与 `-Ob` 读回均 == 输入正文）；`gold/fixtures/vcf44_phasing.vcf` 四路对拍 bcftools 逐字节相同（见下） |

### Y4 门（2026-10-11 实测，除注明外 56 线程）

| 门 | 结果 |
|---|---|
| gold 正文 md5（`run_reed.sh`，全量 `-Oz`） | **`7aedc59ac798729a2674e6e7ef19adfc`**（与 v3–v6 逐字符相同） |
| `check_parity` | **PASS**：1,768,759 / 1,768,759 位点，Jaccard 1.000000，348 样本 MATCH，基因型差 mean=0 max=0 |
| 子集数字门 `--min-dp 5 --min-gq 20 --geno 0.1 / 0.2` | **263,420 / 381,672**（正文 md5 `fef4cb07…` / `371d3a0d…`，== Y3 记录） |
| 报告 md5（8 文件，全量） | afreq/hardy/het/smiss/vmiss/idepth/ldepth.mean/stats 与 Y0–Y3 记录逐字符相同 |
| A6-14 夹具四路对拍（`gold/fixtures/vcf44_phasing.vcf`，8 样本，含 `|0/1`、`/0|1`、`0|.`、`|.`、`-0`） | VCF 入 `-Oz`、VCF 入 `-Ob` 读回、BCF 入 `-Oz`、BCF 入 `-Ob` 读回——四者与 `bcftools view -H` **逐字节相同**；pre-Y4 二进制在 BCF 入路径丢前缀（`|0/1`→`0/1`）、`-Ob` 直接报错（`GT allele index is empty or too long`），留档对照 |
| `-Ob` 全量计时（干净段 A/B） | **VCF 入：post 8.02 s / 417.7–418.2 CPU-s vs pre 468.4 CPU-s（−50.7）**；**BCF 入：post 5.34–5.57 s / 273.2 CPU-s vs pre 519.2–521.8 CPU-s（−246）**；`-Oz` gold post 379.0 CPU-s（Y3 记录 395.5，本项顺带消掉 `-Ob` 外的重编码但 `-Oz` 不走编码器，差异属批次） |
| `-Ob` 读回 ≡ `-Oz` 正文 | 子集（VCF 入与 BCF 入）与全量均逐字节相同（bcftools 读回 md5 == `7aedc59a…`） |
| noodles oracle | `native_and_noodles_encoders_agree` 通过（B2/B4 逐字节相等）；`REED_BCF_ENCODER=noodles` 分支保留 |
| 四件套 | 测试 **388** 通过 / 0 失败（+7：4 编码器/GT 单测、2 渲染单测、1 e2e 夹具）；fmt、clippy（CI 命令）、rustdoc 退出 0；守卫 29/29 ok |

### Y4 偏差登记

| 项 | 说明 |
|---|---|
| **全量 `-Ob` ≤ 470 s 门** | 按 **user** 判定达成：VCF 入 417.7–418.2（门 ≤ 470 ✓，−50.7 vs pre 468.4）；BCF 入 273.2（−246）。wall 在批次负载 1.5–16 下呈双峰（同 CPU 下 5.3 s 与 23 s 并存），pre 二进制在同一环境同样双峰（9.9/30.8/27.7）——属环境，不是回退；取安静段最好值 post 5.34 vs pre 9.89 |
| `-Oz` gold 379.0 CPU-s | 与 Y3 记录的 395.5 相差 16.5：`-Oz` 路径不经编码器，属批次环境差异（本次批次 CPU 普遍偏低，见 BCF 入 519→273 的同时 pre 同批次 515–522），不作为 Y4 收益登记 |
| `-Ob` VCF 入的 50.7 CPU-s | 04 §4.2 估 +61，实测 −50.7，落在估计区间内（差值来自 pre 的 `-Ob` 路径 = edit + reparse + encode 三段，模型版把后两段并为一次扫描） |
| `encode_bcf` 的 GT 相位位 | BCF 源码的相位位原样拷贝（不重算）：掩码样本写 `0,0`（两个 missing，htslib `./.` 语义），非掩码样本逐码拷贝。文本路径的 `gt_codes_into` 复刻 htslib 解析（含 4.4 前缀与缺失等位裸相位位），两路在夹具上字节一致 |



| 项 | 状态 | 说明 |
|---|---|---|
| `--mind` + `--indep-pairwise` 的内存预算 | 已知限制 | 组合强制 sketch 形态（矩阵整基因组常驻，剪枝幸存集依赖全基因组 keep）；`choose_for_mind` 的预算判断对此组合不生效，超预算不报错。缓解见 Y6（或后续按 contig 两段式） |
| LD bp 窗口（F6-10） | **已修**（`a5cd8a2`） | 区域 `50kb 0.2` 158 == plink2 对称差 0；02 §3.2 的 M6-4 行随之关闭 |
| C-5、联合 SFS、`--bed-file`/`--sites-file`、前后端收尾三项 | 待做 | 并入 Y6 |
| `is_sequential` / 顺序模式（F6-6） | **已实现**（Y5） | `run()` 按 `budget.is_sequential()` 分派 `run_sequential`（单线程读/解压/计算/写出，复用 `NewlineStitcher` 与 sink）；输出与并行逐字节相同（`sequential_and_parallel_outputs_agree`） |
| `JointSfsAcc` 仍 `#[doc(hidden)]` | 有意 | 无 CLI 出口；Y6 接入联合 SFS 时去掉 |
| 记录模型的范围 | 已知限制 | 落地为 `NativeEncoder` 的 `encode_text`/`encode_bcf` 两个入口（03 §8 的 `Record` 结构体形态未采用：两条源路径各自的扫描即模型，省掉中间表示与一次遍历；`RecordSource`/`RecordSink` trait 待有第三种格式时再抽）。`--set-gts .` 表达式路径经 `plan.mask` 全为 field 规则时仍走 `encode_text`；实测 `--set-gts . -e 'FMT/DP < 10'`（含 `-s` 子集与 `-I`）的 `-Ob` 读回 == `-Oz` 正文（夹具上 4 记录逐字节） |

## Y5：引擎健壮性（2026-10-11）

规格：`docs/spec/v6_y5/README.md`（含与 roadmap 附录的差异表）与 `CODE.md`（难点参考实现）。
实现按规格执行，落地形态与 roadmap 附录不同的四处登记在 §差异。

### Y5 已落地项

| 项 | 内容 | 验证 |
|---|---|---|
| S1 中止（A6-1a） | `AbortSignal`/`AbortHandle` 共享 `Arc<AtomicBool>`，`is_aborted` 读 flag（token 被消费不再回假）；`abort_channel_for(threads)` 容量 `threads+4`、`abort_for(&ThreadBudget)`；六个调用点全换（filter 2、index、stats、reports、matrix）；`AbortHandle::register(&Arc<BytePermits>)` + `trigger` 的 `notify_all`，`BytePermits::acquire` 改无超时 `cv.wait`（§7.5 不再轮询） | `consuming_the_token_does_not_clear_the_flag`、`trigger_wakes_registered_permit_waiters`、`abort_at_eight_times_the_cores_returns_without_leaking_permits`（8×核数 + 64 KiB 预算，< 2 s 返回 `Aborted`） |
| S2 许可 drop guard（A6-1b） | `PermitGuard{Arc<BytePermits>, Arc<Shared>, bytes}`；`InFlightChunk`/`WorkerMsg` 持 guard；reader 三处手写 release 删除；worker 一次解构、失败即 drop；reducer 交付后 drop；`StageMetrics::in_flight_final_bytes` 断言归零 | `permit_guard_returns_its_bytes_on_drop`、`worker_panic_returns_every_permit`（panic 后 2 s 内以 `WorkerPanicked` 收口）、`stage_metrics_split_inflate_and_job_time` 的 `in_flight_final_bytes == 0` |
| S3 接缝上限（A6-2） | `pipeline::MAX_RECORD_BYTES = 1 GiB`（`bcf_input.rs` 改引用）；`NewlineStitcher` 带 `grow` 上限；并行 reducer 与顺序路径共用；删 `stitch_seam`/`finish_fragment` | `stitcher_reassembles_and_caps_fragments`、`stitcher_rejects_a_record_past_the_ceiling`、`stitcher_reassembles_records_across_ten_chunks`（10 B chunk × 100 B 行） |
| S4 顺序模式（A6-3） | `run_sequential`（`Newline` 与 `Raw` 两分支）；`run()` 按 `is_sequential()` 分派；单 chunk 超预算同报 `MemBudgetExceeded`；`budget.rs` 注释与 `architecture.md:804-806` 的 F6-6 现状段改写 | `sequential_and_parallel_outputs_agree`（threads 1/2 vs 8，fold 相同） |
| S5 计量（A6-4） | `StageMetrics` 增 `reader_wait`/`reducer_wait`/`inflate_time`/`job_time`/`deflate_time`/`in_flight_final_bytes`；reader 三处阻塞、reducer `recv` 等待、worker 两段墙钟入原子；deflate 经 `ChunkOutput.deflate_nanos` → `OrderedOutput` → `OrderedSink::deflate_time()`（引擎在 `finish` 前读一次）；`TimingReport` 五字段 + `timing_text` 五行 | `stage_metrics_split_inflate_and_job_time`；子集 `--timing` 实测 `reader_wait_ms`/`reduce_wait_ms`/`inflate_ms`/`job_ms`/`deflate_ms` 全非零（deflate 38.7 s / job 51.3 s） |
| 文档 | `docs/spec/v6_y5/{README,TODO,CODE}.md`；`docs/README.md` 索引补一行 | 人工 |

### Y5 门（2026-10-11 实测）

| 门 | 结果 |
|---|---|
| 测试 | `cargo test --release --workspace`：**398 通过 / 0 失败**（388 → +10：4 引擎集成、2 引擎内联、4 既有调整） |
| gold 正文 md5（`run_reed.sh`，全量 `-Oz`） | **`7aedc59ac798729a2674e6e7ef19adfc`**（与 v3–v6 逐字符相同） |
| `check_parity` | **PASS**：1,768,759 / 1,768,759 位点，Jaccard 1.000000，348 样本 MATCH，基因型差 mean=0 max=0 |
| 子集数字门（`subset.Chr01-03.vcf.gz`） | **263,420 / 381,672**，正文 md5 `fef4cb07943678d45893c5cfdd0c91f9` / `371d3a0d8f1cab0ebb57966cc9a507df`（== Y3/Y4 记录） |
| 子集跨格式 | `-Ob` 读回 == `-Oz` 正文（`fef4cb07…`）；子集 BCF 输入 `-Oz` == 同 md5（263,420） |
| 子集 `--mind 0.1` t1/t56 | 两份正文 md5 `fcc6494c4db3068923599acbc3651819`（== Y1e/Y3 记录） |
| 报告 md5（8 文件，全量） | afreq `b2884092…`、hardy `fd8fbf27…`、het `790293d9…`、smiss `5f4fe4f8…`、vmiss `7c325d64…`、idepth `a4b6dd3a…`、ldepth.mean `475ee083…`、stats `189e025b…` == Y0–Y4 记录 |
| plink2 LD 夹具（`gold/popgen/run_plink2_ld.sh`） | PASS r²（15 对，1.02e-06）+ PASS prune（3 保留/3 剪除） |
| 子集 filter 计时 ×3 | wall 1.121 / 1.771 / 1.118 s（run2 在负载 3.67 上；无回退，Y4 记录同量级） |
| 四件套 | fmt、clippy（CI 命令）、rustdoc 退出 0；守卫 **29/29 ok** |

### Y5 偏差登记

| 项 | 说明 |
|---|---|
| 与 roadmap 附录 S2 的形态 | 附录 `PermitGuard<'a>` 借引用；落地改 `Arc` 持有，避免生命周期参数扩散进 `WorkerMsg`/`ReorderBuffer`/channel。每 chunk 两次原子增减，可忽略 |
| 与 roadmap 附录 S3 的形态 | 附录抽 `SeamStitcher` trait 并把 BCF `FrameChain` 也实现；未做。BCF 走 `RecordFraming::Raw`，接缝在 job 内的链上，reducer 不碰；trait 只有一个实现类型（单实现抽象）。只修上限（A6-2 的真实缺陷） |
| 与 roadmap 附录 S4 的形态 | 附录 `run_sequential(&Abort) -> RunReport` 且只写换行路径；落地为 `&AbortSignal` + 返回 tuple 并含 `Raw` 分支：① 顺序模式 BCF 输入仍是 Raw；② 直接返回 `RunReport` 会绕过 `run` 的 FirstError 收口 |
| 与 roadmap 附录 S5 的形态 | 附录 `ChunkJob::stage_nanos() -> [u64;3]` 回报 inflate/job/deflate；未做。`inflate_chunk` 在引擎 worker（job 看不到），deflate 在 per-worker `Scratch`（job 是 `&J` 共享）。落地：inflate/job 引擎直接测，deflate 随 `ChunkOutput` 到 sink |
| `FrameChain::wait_inner` 仍 10 ms 轮询 | 它同时等 abort/poison/finish 三源，本轮不重排通知路径；中止延迟 ≤10 ms，满足 2 s 门。`BytePermits` 的轮询已按 §7.5 消除 |
| 全量计时未重测 | Y5 不改计算路径（引擎骨架、计量、许可、接缝），子集 ×3 无回退 + gold md5/parity 逐字节不变即达门；全量 filter `-Oz` user 记录沿用 Y4（417.7–418.2） |

## tag

`v0.3.0` 打在 Y0 的收口提交上（`git describe --tags` 从此可解析）。

## 基准（2026-10-10 v6 复审基线，56 线程，热缓存，批次开始负载 1.28；见 `04-performance.md` §1）

| 命令 | wall | user |
|---|---|---|
| filter gold `-Oz6`，VCF 输入 | 8.85 s | 460 s |
| filter gold `-Oz6`，BCF 输入 | 9.00 s | 464 s |
| filter gold `-Ob` | 9.98 s | 521 s |
| filter gold `-Oz1` | 5.64 s | 289 s |
| stats 五统计 VCF / BCF | 1.78 / 1.65 s | 72 / 56 s |
| reports 两族 + summary | 3.0–5.4 s（双峰） | 147 s |
| **Y2 实测（2026-10-11）**：filter gold `-Oz6` VCF | **7.64–7.70 s** | **400 s**（−60） |
| **Y2 实测**：stats 五统计 VCF | 1.72–1.77 s | **70.6–70.9 s**（−1.4，见 Y2 偏差） |
| **Y2 实测**：reports 两族 + summary | **3.09–3.13 s（双峰消失，/tmp 与 tmpfs 一致）** | — |
| **Y3 实测（2026-10-11）**：filter gold `-Oz` VCF，56 线程 ×3 | **7.62–7.68 s** | **395.4–395.7 s**（pre-Y3 同协议 398.2 → −2.6，门 −30 未达，见 Y3 偏差） |
| **Y3 实测**：`--fill-tags HWE`（子集，pre-Y3 消融） | — | 0.9 CPU-s / 261,788 幸存者 |
| stats 边过滤 `--pi` | 3.09 s | 153.6 s |
| index VCF / BCF | 1.63 / 1.12 s | 75 / 44.5 s |
| `--sfs` / `--fst-global` / f 统计 | 1.60 / 1.52 / 1.56 s | 73 / 71 / 71.5 s |
| `--ld-decay 100k` | **8.4 s**（`--pop`）/ **6.4 s** | 413 / 337 s |
| `filter --indep-pairwise 50 5 0.2 -Oz1` | **4.0 s**（Y1c）→ **3.6 s**（Y1e 流式 contig 剪枝） | 181 s |
| `filter --mind 0.1 -Oz1`（三形态） | **6.5 / 6.5 / 6.6 s** | — |
| plink2：pgen 导入 / 从 pgen 剪枝 / 从 pgen 报告 | 27.8 / 0.66 / 1.48 s | 105 / 1.5 / 2.75 s |

> Y0 未触及计算路径（删除的死代码无生产调用者；`memchr` 与 `f64::midpoint` 的替换经 md5/对拍验证逐字节不变），故不重跑全量基准。Y1 起按 04 §9 的纪律重测。
