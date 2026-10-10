# 01 第六轮复审（v5 X0–X7 落地后，2026-10-10）

> 复审对象：`d735802`（doc_v5）… `8709f1f`（SFS 子集对拍）共 23 个提交；工作树干净。
> 方法：构建/测试/clippy/守卫实测；全量基准 3 次中位数 + 线程扫描（04）；gold 全量对拍重跑；五路分 crate 深读审计（engine+core、io 两路交付并经我复核；filter、stats、facade+cli 三路因推理网关超时两次未能交付，这三层以我一手阅读的热路径为准，§3 列出覆盖与局限）；参照实现的源码（plink2、htslib）只引用实际读过的行。

---

## 0. 总评

v5 的路线几乎全部落地：X0–X3、X4.1–X4.4（含剪枝迁入 filter）、X6 的 `-q/-Q :type`、X7 CI 都有提交、测试与对拍数字；gold 正文 md5 在本轮重跑中不变；四个新统计量各只增加 0–2 CPU-s。流式一遍的形态已经在 Amdahl 极限上（filter 56 线程效率 96%，user 在 8/16/32/56 线程上恒定）。

本轮的问题集中在三处，而且都是**结构**问题，不是常数：

1. **全局分析走错了执行模型**：`--indep-pairwise` 66 s、`--ld-decay` 169 s，全部时间在 reducer 线程里串行（`--timing` 的 `reduce_ms 63,967 / 64,436`）；plink2 从内存矩阵剪枝只要 0.66 s 且保留集合与 reed 完全相同（235,403）。模块注释说的"每 contig 自己的线程"没有实现。
2. **CI 的 lint 门是红的**：`cargo clippy --workspace --all-targets -- -D warnings` 在 reed-filter 上 2 个 error（doc lint）+ 测试目标约 60 条 pedantic；STATUS 写 "clippy 0 error" 用的不是 CI 的命令。
3. **一批已证实的浪费与潜伏缺陷**：deflate/inflate 每块 65 KiB 清零、`-Ob` 文本往返 +61 CPU-s、流式 chunk 源在空成员后虚拟偏移漂移、abort 通道按核数而非 `--threads` 定容、计数内核在 release 下不健全（且永不向量化、无调用者）、顺序模式规范有代码无。

---

## 1. 状态核实

| 项 | 结果 |
|---|---|
| 构建 | `cargo build --release --workspace --all-targets` 0 警告 |
| 测试 | **382 通过 / 0 失败**（27 个目标） |
| clippy（CI 命令） | **失败**：`reed-filter/src/eval/tail.rs:75-76` 两处 `doc_lazy_continuation`（`-D clippy::all` 下为 error）+ 测试目标 pedantic ≈ 60 条 |
| 守卫 | `bash ci/check-doc-drift.sh` 23/23 ok |
| gold 全量 | 正文 md5 `7aedc59ac798729a2674e6e7ef19adfc`（两种算法）；`check_parity` PASS，1,768,759 位点，Jaccard 1.0，基因型差 0 |
| 跨格式 | stats 五统计 VCF ≡ BCF（`_pi.txt` `45c2f2e1…`、`_fst.txt` `fa144045…`） |
| 跨线程 | filter 8/16/32/56 线程 user 457–459 s 恒定；stats 70.3–71.7 s 恒定 |
| 远程 | `origin` = github.com/WWz33/reed；`ci/publish.sh` 推送剔除 docs/bench 的子树 |

## 2. STATUS 自述 vs 事实

| STATUS 条目 | 核实 | 结论 |
|---|---|---|
| X0 stash、主干 364 通过 | 现 382 通过（X4/X6 新增 18 项） | ✓ |
| X1a/X1b：C-1/2/3/8/9、前后端三规则、L-1…L-7 | 守卫 23 条全 ok；`reed-cli/Cargo.toml` 只依赖 `reed` | ✓ |
| X2 目录搬迁、X7 CI | `docs/{spec,rounds}`、`ci/`、`bench/suites/`、`target/gold/`；`.github/workflows/ci.yml` 存在 | ✓；但 CI 的 clippy 步骤今天必红（§1） |
| X3 reports ≤ 3.2 s / ≤ 155 CPU-s | CPU 147 s ✓；wall **双峰** 3.0–5.4 s（两组各 3 次），未达"稳定 ≤ 3.2" | 部分 |
| X4.1–X4.3 SFS / fst-global / f 统计 | 全量各 1.5–1.6 s（五统计 1.78 s 之内）；夹具与子集对拍数字在 STATUS | ✓ |
| X4.4 LD："chunk 内 carry + sink 补跨 chunk 对" | 衰减 169 s / 343 CPU-s；chunk 跨度 < 100 kb 窗口 → 绝大多数对在 sink 串行补算 | 正确但**执行模型错** |
| X4.4-rev `--indep-pairwise` 进 filter；`ld_prune.rs:7-9` "contig 并行" | `LdPruneSink` 一个 `Option<LdPruner>` 串行（`ld_prune.rs:318-369`）；pass0 64.4 s 中 reduce 64.0 s | 结果正确（全量 235,403 == plink2），注释不实 |
| X6 `-q/-Q :type` 40/40 | 代码与测试在 | ✓ |
| "clippy 0 error" | CI 命令下 2 error | ✗ |
| `--ld-decay` 可用 | 需要 `--pop`（`translate/stats.rs:439-444` 把它当窗口统计），LD 与群体无关 | 设计瑕疵 |
| gold md5 未重跑（X2 偏差） | 本轮已重跑，不变 | 已补 |

## 3. 审计覆盖与局限

| 层 | 覆盖方式 | 一手核实的关键文件 | 未覆盖 |
|---|---|---|---|
| 引擎 + 核心 | 子审计交付 + 我复核每条 P1/P2 | `pipeline.rs` 全文、`kernel.rs:40-135`、`abort.rs:60-70`、`gt_code.rs`、`planes.rs` | — |
| I/O | 子审计交付 + 我复核偏移漂移、deflate/inflate、line_edit、scanner | `genotype_scan.rs` 全文、`deflate.rs`、`block_writer.rs`、`line_edit.rs:77-216`、`chunk/bgzf.rs:120-140,300-315`、`native_encode.rs:1-45` | 索引写出器与 CSI 内容只按子审计 |
| 过滤 | **子审计未交付**；一手阅读 | `record_pipeline.rs` 全文、`single_pass.rs` 全文、`plan/strategy.rs`、`expr/mod.rs` 的 `eval/value/info_value`、`eval/head.rs` 轮廓、`edit_plan.rs` 的 `fill_info_set`、`sketch.rs:40-75` | 表达式语义与 bcftools 的逐条对照（现有对拍：gold、`-q :type` 40 组、FMT 表达式测试）；词法/语法器 |
| 统计 | **子审计未交付**；一手阅读 + 对拍数字 | `request.rs` 的 push/merge、`window.rs:53-135`、`job.rs` 的 `push_row`/`StatsSink`、`fst.rs` 多群体、`ld.rs` API、`sfs.rs` 在途改动 | 各估计量公式的逐行复核（以 v5 05 §1 的对拍结果为准：pixy/allel/vcftools/plink2/easySFS 逐位一致） |
| 门面 + CLI | **子审计未交付**；一手阅读 | `filter/mod.rs:180-300`、`ld_prune.rs:300-405`、`filtered.rs` 轮廓、`stats/job.rs:915-1000,1296-1500`、`translate/stats.rs` 的泄漏行、`cli.rs` 的统计参数块 | `index.rs`、`reduce.rs`、`mid_stage.rs`、`emit.rs` 的逐行审计 |

所有"最优方案"（02 §2.2、03 §7–9）都只依赖一手核实的事实。

## 4. 发现

| # | 级 | 位置 | 问题 | 证据 | 去向 |
|---|---|---|---|---|---|
| F6-1 | **P1（性能/模型）** | `reed/src/filter/ld_prune.rs:318-369`；`stats/job.rs` 的 `ld_tail` | LD 剪枝与衰减在有序归约点串行：66 s / 169 s；plink2 0.66 s | 04 §6 `--timing`；plink2 实测 | 02 §2.2、03 §7、06 Y1 |
| F6-2 | **P1（工程）** | `reed-filter/src/eval/tail.rs:75-76`；测试目标 | CI 的 clippy 门失败；STATUS 陈述不实 | 05 §1 | 06 Y0 |
| F6-3 | P2（正确性，潜伏） | `reed-io/src/chunk/bgzf.rs:131-135,310` | 流式源丢弃空成员后 `compressed_start` 漂移，`VirtualMap` 的虚拟偏移错；仅 stdin/管道路径可达，无测试 | 03 A6-7 | 06 Y2 |
| F6-4 | P2（健壮性） | `reed-engine/src/abort.rs:63-67`、`pipeline.rs:355-368` | 中止通道按核数定容；失败路径不释放许可 → 可能挂起而非 `Aborted` | 03 A6-1 | 06 Y5 |
| F6-5 | P2（soundness） | `reed-core/src/count/kernel.rs:46,121-133` | release 无长度检查的裸 AVX-512 加载；348 样本下永不向量化；无生产调用者 | 03 A6-22 | 06 Y0（删）或 Y5（修） |
| F6-6 | P2（规范偏差） | `reed-engine/src/budget.rs:46-50`、`pipeline.rs:157`、`docs/spec/architecture.md:804` | 顺序模式承诺未实现 | 03 A6-3 | 06 Y5 |
| F6-7 | P2（性能） | `single_pass.rs:338-345`、`native_encode.rs:131-139` | `-Ob` 文本往返 +61 CPU-s；BCF 输入存活记录渲染→解析→编辑 | 04 §1 | 03 §8、06 Y4 |
| F6-8 | P2（性能） | `deflate.rs:134-135`、`block_writer.rs:183`、`inflate.rs:36,49` | 每块 65 KiB 清零 + 分配，≈ 38 GB memset | 03 A6-5/6 | 06 Y2 |
| F6-9 | P2（设计） | `reed-cli/src/translate/stats.rs:439-444` | `--ld-decay` 被当作窗口统计强制 `--pop` | 实测报错 | 06 Y1 |
| F6-10 | P2（一致性） | `reed-stats/src/ld.rs` 窗口推进 | `--indep-pairwise Wkb` 与 plink2 不一致（131 vs 158） | v5 STATUS 偏差 | 06 Y1 |
| F6-11 | P3（性能） | `single_pass.rs:290,332`；`record_pipeline.rs:345-348,357-372,395-399`；`site_counts.rs:71-81`；`line_edit.rs:123-124,250-302`；`record_pipeline.rs:330`；`chunk/mod.rs:231-246` | 每记录/每样本分配、f64 转换、逐位清零、未缓存 HWE、每 worker GT 表、标量 `rposition` | 03 §3–4 | 06 Y2/Y3 |
| F6-12 | P3（架构） | `single_pass.rs:278-346` vs `:353-412`；`reduce.rs`/`sketch.rs`/`ld_prune.rs` | 两套几乎相同的发射闭包；三套 pass0 | 03 A6-20/21 | 06 Y1/Y3 |
| F6-13 | P3（前后端） | `translate/stats.rs:330,438,507,660`；`main.rs:19`；`translate/filter.rs` 仅 3 个往返测试 | 后端逻辑在前端；无统一 `reed::Error`/`kind()`；stats/index 无往返测试 | v5 STATUS 遗留 | 06 Y6 |
| F6-14 | P3（模型/死面） | `contig.rs:9` vs `id.rs:35`；`math/hypergeom.rs`；`field_plan.rs`；`pipeline.rs:23`；`budget.rs:244-246` | 两个 `ContigId`；重复的投影权重；无消费者的模块与项 | 03 §5 | 06 Y0 |
| F6-15 | P3（文档） | 05 §4 的 8 处 | 注释/规范与代码不符（含 `ld_prune.rs:7-9`） | 05 §4 | 06 Y0 |
| F6-16 | P3（测试） | 05 §3 | 引擎边界、两种源等价、索引内容、fuzz、往返 | 05 §3 | 06 Y7 |
| F6-17 | P3（观测） | reports 模式 | wall 双峰 3.0 / 5.4 s，CPU 恒定 | 04 §7 | 06 Y2 |
| F6-18 | P3（渲染） | `bcf/render.rs:466-489`、`line_edit.rs:260-262` | VCF ≥ 4.4 前缀相位未实现；`-0.0` 打成 `0` | 03 A6-14 | 06 Y4 |
| F6-19 | P3（产品） | `stats/job.rs` 1,936 行；`individual.rs`；联合 SFS；`--bed-file/--sites-file` | v5 的 C-7、C-5、联合谱、窗口层 Regions 路径都未做 | v5 STATUS | 06 Y6 |

## 5. 性能摘要

见 04：filter 8.85 s（不回退）、BCF 输入 9.00 s、stats 1.78 s、index 1.63 / 1.12 s、新统计量 ≈ 五统计成本、剪枝 66.6 s、衰减 169 s、reports 双峰。与 plink2：流式任务领先一到两个数量级，全局分析落后两个数量级（执行模型）。

## 6. 值得保留的做法

- 票据 + 重排环 + 长度前缀链的引擎：形态正确、不变量可证、效率在极限上。
- `RecordPipeline` 的单一判定实现（文本与 BCF 原生两入口共用 `run_stages`），使 pass0/pass1 的存活序列天然一致。
- 可加 partial + `NeumaierSum` + 按 chunk 顺序合并的确定性模型，与跨线程 sha256 门。
- 差分 oracle 纪律（noodles、v0.2 oracle、bcftools 夹具、proptest）与每条 STATUS 行引用测试或命令——这次能逐条核实，仍然靠它。
- 四个新统计量的接入没有破坏任何基线。
