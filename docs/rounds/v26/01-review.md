# 01 复核清单：IO / filter / 并行 / noodles

> 基线核对日期：2026-11-05。文件行号以实施前重新定位为准。`confirmed` 表示静态调用链已确认；回归测试仍需按路线补齐。`verify` 表示旧 review 主张尚缺可复现证据。

## A. 确认的问题

| ID | 优先级 / 状态 | 位置 | 触发与影响 | 最小处理方向 |
|---|---|---|---|---|
| V26-00 | **P0 / confirmed unsoundness** | `reed-io/src/chunk/bgzf.rs::MmapChunkSource::open` unsafe block | safety comment 认为“reed 只读打开、自己不写”即可保证 backing file 不变，但另一进程仍可原地改写或 truncate 同一个文件。`memmap2::Mmap::map` 的安全合同要求映射期间底层文件不被任何进程修改；违反合同属于未定义行为，而不只是读到坏数据/SIGBUS。共享工作目录/NFS/HPC 中文件被替换或截断时可触发。 | **建议默认关闭 mmap 并走 streaming**，除非产品愿意要求输入文件在运行期间不可变，且该前提能由实现/部署保证；仅 advisory lock 不足以约束不合作的写入者。请在 roadmap 拍板采用“安全默认”还是“受约束的 mmap opt-in”，再改 unsafe 入口。 |
| V26-01 | **P0 / confirmed** | `reed-io/src/bgzf/block.rs::is_empty_member`、`chunk/{bgzf,region}.rs`、`chunk/mod.rs::inflate_chunk`、`bgzf/inflate.rs` | block 被 footer 的 ISIZE=0 标记为空时，读取路径会跳过 inflation/CRC 验证；损坏的非空 block 可丢记录。若 ISIZE 被改成很大的数，`inflate_chunk` 会先以该值扩容，之后才由解压器判断内容，可能 OOM。 | 在 reed-io 的共享 BGZF block 层验证完整成员、CRC 与 `ISIZE <= 65536`；只允许验证成功的空成员跳过。stream、mmap、region 复用同一入口。回归覆盖 CRC 错、ISIZE=0 但 payload 非空、ISIZE 越界、合法空成员。 |
| V26-02 | **P0 / confirmed** | `reed-filter/src/semantics.rs::info_af_outcome` | `alts.iter().sum::<u64>()` 没有 checked 加法。输入 `AC=18446744073709551615,1;AN=100` 可在 debug panic，release 回绕后误算 AF。整数列表在开始累加前也可用超长字段消耗内存。 | checked 累加；超出 VCF/实现数值域时返回明确 malformed INFO 错误，不得当作 usable counts。限制/验证 AC 元素数量与整数域。 |
| V26-03 | **P0 / confirmed** | `reed-engine/src/budget.rs::Threads::from_str/workers/tickets`、`abort.rs::abort_channel_for` | 任意非零 `usize` 都能作为线程数。极大值会在 `threads + 4`、`TICKETS_PER_WORKER * workers` 或队列分配处溢出/请求不现实资源。 | CLI 边界给出有依据的上限；budget 派生容量使用 checked arithmetic，所有 CLI 命令复用同一域。 |
| V26-04 | **P1 / confirmed** | `reed-io/src/vcf/info.rs::InfoIndex`、`reed-filter/src/semantics.rs::info_af_outcome`、`expr/mod.rs::info_value` | 重复 INFO 标签被接受。索引/表达式读首个匹配项，而 AC/AN reader 遍历后保留最后一个；例如 `AC=1;AC=3;AN=4` 会使同一记录按不同 AC 值过滤。 | 在共享 INFO 校验处拒绝重复 key；不要让每个 consumer 自选 first/last。检查其他 INFO consumer 一并走该约束。 |
| V26-05 | **P1 / confirmed** | `reed-io/src/bcf/gt.rs::site_header`、`bcf/render.rs::render_site_part` | `n_allele=0` 未在 BCF site header 入口拒绝；renderer 会把它渲染成 REF/ALT 为 `.` 的记录，native GT 路径也将其继续当成 0 allele。上游本地 noodles 0.119 已在 BCF decoder 拒绝无 allele 记录。 | site header 统一要求至少一个 allele（REF）；测试 native 输入、render、GT 三条消费路径均返回 MalformedRecord。 |
| V26-06 | **P1 / confirmed** | `reed-io/src/vcf/line.rs::LineRef::parse/classify_line`、`site.rs::SiteView::parse` | 行分类只验证列数、样本列数和 POS；未拒绝空 CHROM/REF/ALT/FILTER/INFO/FORMAT 等必需字段。像 `chr1\t1\t.\t\tG\t.\tPASS\t.` 可被分类为 record，后续 filter/stats 可能接受或原样输出非法 VCF。 | 在 `reed-io` 增加共享的结构字段验证（区分规范允许的 `.` 与空字段）；所有 VCF record 消费者继续只用 `classify_line`。对照 noodles 0.119 新增的 empty/truncated field 拒绝行为。 |
| V26-07 | **P1 / confirmed** | `reed-engine/src/pipeline.rs::NewlineStitcher`、`reed/src/bcf_input.rs::frame_chunk`、`ThreadBudget::mem` | `--memory` 许可只估算压缩 chunk（`compressed_len * 8`），而跨 chunk 的 VCF fragment 和 BCF `Carried.head` 可各自增长至 1 GiB；carry 不计入许可。输入含超长单记录时，即使指定小预算也可占用远超预算的内存。 | 把 carry 纳入预算或让单记录上限受预算约束；文档同时明确 `--memory` 的真实保证范围。单独 chunk 估算与跨 chunk carry 分别测试。 |
| V26-08 | **P1 / confirmed** | 输出端：`reed/src/filter/mod.rs`、`stats/sink.rs`、`stats/writers/*`、`convert/writers.rs`、`pca/mod.rs`、`dapc/mod.rs` | 多处直接 `File::create` / `fs::write` 覆盖目标。运行中写失败、进程中断或后续输出失败时，旧的完整结果已被截断，新结果留下半文件；多文件报告还可能新旧混杂。 | 先在同目录写临时文件，完成并 flush/sync 后 rename；多文件输出明确是逐文件原子还是整组事务。优先修会覆盖用户既有结果的主输出。 |
| V26-09 | **P1 / confirmed** | `reed-engine/src/budget.rs::physical_memory_bytes` | 默认预算只读 Linux `/proc/meminfo` 的宿主机 `MemTotal/2`，不识别 cgroup/container memory limit；受限容器可能按宿主机大内存放行，随后被 OOM killer 终止。 | 探测 cgroup v2/v1 限额并与宿主机值取小；解析失败保守回退。增加模拟 cgroup 文件的测试。 |
| V26-10 | **P1 / confirmed** | `reed-io/src/header.rs::HeaderFacts`、`reed/Cargo.toml`、workspace `Cargo.toml`、`Cargo.lock` | reed 当前锁 `noodles 0.117.0`，约束 `^0.117` 不会更新到 `.118/.119`；`reed` 直接依赖并暴露 `noodles::vcf::Header`，破坏 reed-io 唯一边界。`reed-io` 还在 dev-dependencies 重复声明 noodles。 | 先隐藏 Header 类型，只通过 reed-io 提供 header metadata/操作；删重复 dev-dependency；单独将 noodles 升到本地已有 `noodles-0.119.0` tag（`9986db41b`）并跑格式兼容测试。**不建议整体重写 CSI/TBI/VCF header/BCF API。** |
| V26-11 | **P2 / confirmed** | `.github/workflows/ci.yml`、`rust-toolchain.toml`、workspace `rust-version` | workspace 声明 MSRV 1.91，仓库工具链固定 1.99，常规 CI 只安装/构建 1.99；依赖更新或语法变化可能无意间抬高 MSRV而 CI 仍绿。 | 增加固定 1.91 的 `cargo check --workspace --all-targets`（或明确下调 MSRV 承诺）；工具链通道与最低版本分开说明。 |
| V26-12 | **P2 / confirmed, flaky** | `crates/reed/src/spill.rs` tests、`convert/job.rs::spill_tests` | `exhausted_candidates_error_without_clobbering` 预测全局计数器的后续 16 个名字；多个测试会增加同一个 `SPILL_COUNTER`，但不是每个递增测试都持 `SPILL_TEST_LOCK`。并行时计数器可在植入文件和调用间推进，测试偶发创建成功。单独/串行通过不能证明并行安全；v25 前后的并行全跑曾出现此失败。 | 测试不依赖全局 counter 的具体值（注入命名器/独立目录状态），或所有消费 counter 的测试统一锁。连续跑默认并行 harness 验收。 |
| V26-13 | **P2 / confirmed** | `reed-io/src/chunk/plain.rs::open_inner`、`chunk/bgzf.rs::scan_header` | 找到 `#CHROM` 前无 header 累积上限。畸形输入缺少 `#CHROM` 且很大时，open 阶段会把 header 全部累积到 Vec，耗尽内存后才报 Header。 | 设明示的最大 header 字节数，超限返回 Header error；测试超限与跨 block 的合法大 header。阈值需考虑大样本/大 header 文件。 |

| V26-18 | **P2 / confirmed scalability optimization** | `reed-engine/src/frame_chain.rs::State::{cells,waiting}` 与 `publish` | 两个 Vec 以全局 chunk id 作索引；每个 chunk publish 后 cell 被消费为 `None`，但 Vec 不回收/复用旧槽位。BCF 长流会保留 O(total chunks) 的 cell 与 waiter 位图；1 MiB chunk 的超长输入会持续增加少量元数据。有限向量不是无界内存 bug，但随文件长度线性增长。 | 若目标负载包含超长 BCF 流，按 ticket window 改成有界 ring；否则保持当前实现。先测 metadata RSS，避免为常规文件引入复杂回收同步。 |
| V26-21 | **P1 / confirmed resource-contract bug** | `reed-reduce/src/backend.rs::faer_gram/faer_sym_eig`、`reed/src/pca/mod.rs::solve` | Dense PCA calls faer with `Parallelism::Rayon(0)`. Faer defines 0 as `rayon::current_num_threads()`; PCA solve runs directly after matrix collection, outside `ThreadBudget::run_auxiliary`, so `--threads 1` or `--threads 2` can still use Rayon’s ambient global pool, and `--threads N` is not the bound used by the eigensolver. | Thread the requested parallelism into reed-reduce and pass `Parallelism::None` or a budgeted Rayon pool matching `ThreadBudget`; do not silently set faer’s process-global setting. Add CLI tests that inspect max active threads for `--threads 1/2/N` and preserve dense result parity. |
| V26-20 | **P1 / confirmed budget undercount** | `reed/src/pca/mod.rs::MatrixUse::peak_bytes`、`reed-reduce/src/backend.rs::faer_sym_eig` | dense estimate 只给 GRM 计一个 `N×N` 阵列；实际 eigensolver 会同时构造 faer 输入矩阵 `am`、特征向量矩阵 `u`、工作 stack，再把 `u` 复制成返回的 `Mat`，原 GRM 与 dosage matrix 此时也仍存活。大样本时多份额外 `N²` 会远超 `--memory` gate。 | 以各 solver 的真实同时存活对象核算峰值；尽可能消费/复用阵列，至少把矩阵和 faer stack 纳入估计。用 RSS 峰值与估算值对照多个 N，明确估算误差边界。 |

## B. v25 候选逐项复核：不把旧 review 当事实

| 候选 | 复核结论 | v26 处理 |
|---|---|---|
| CI-1 gold parity 没有构建/生成 reed 产物 | **确认。** `.github/workflows/ci.yml` 的 data-gates 运行 `reed filter` 做子集计数；随后 `gold/scripts/run_all.sh` 只跑 bcftools/vcftools/plink2，`check_parity.py` 却读取 `out/reed.vcf.gz`。该文件可缺失或来自 self-hosted runner 的旧运行。 | 纳入 P1：data-gates 内 cargo build 当前 checkout，并在 parity 步骤显式运行 `gold/scripts/run_reed.sh`，验证产物来自该 commit。 |
| `N_MISSING` 与 `-g miss` 口径不同 | **不是当前已证实的 bug。** `GtCounts.n_missing` 是 fully-missing genotype 数；`-g miss` 有意把 half-call 也看成 missing，且核心计数另有 `n_half_called`。旧 review 把“统计口径需明确”写成“口径分裂”，证据不足。 | 不作为 bug。补充 docs/spec/CLI 对 fully-missing 与 half-call 的定义和一条对照测试即可。 |
| `--hwe midp` 永远拒绝 | **已明确记录的能力限制，不是缺陷。** CLI spec 明说目前只实现 exact；compiler 给出具体拒绝信息。 | 不列 bug；只有决定实现 mid-p 时才重新开功能项。 |
| reed-reduce 已知答案测试为空 | **旧主张不成立。** `eig.rs` 有小矩阵已知特征值测试，PCA/GRM/Lanczos 也有手算/不变量与 dense/stream 对比测试。 | 不列“空测试”问题；v26 可补独立 oracle/known-answer，但需先指出具体缺失数学性质，不能沿用该标签。 |
| MSRV 1.91 与固定 1.99 | **CI 验证缺口成立**，见 V26-11；不是已证明 1.91 编译失败。 | 增加 MSRV job。 |
| R8 popg bare route / `-o` 契约 | **未确认行为 bug。** 代码目前对会写文件的统计要求 `-o`，对 stdout/table 类目标另有路由；尚未构造唯一复现证明与 CLI 文档冲突。 | 标成 verify，不改变行为；先列出所有 flag 组合并与 cli.md 合约对照。 |
| R9 LD-pruning 同 POS 账目 / r² 边界 | **本轮未做 plink2 复现，不能结论属实。** 代码有测试，但未证明覆盖重复位置与浮点边界。 | 标成 verify；构造最小重复 POS fixture，对照 plink2 的 keep set 和报表计数。 |
| R10 mask attribution 顺序 | v25 已有策略等价 e2e；不再列遗留。 | 关闭。 |
| R13 输出原子性 | 代码仍直接 truncate/write，见 V26-08。 | 确认后纳入。 |
| `stats site_terms` / distance `decode_dose` / kNN helper / 大文件拆分 | 主要是维护与重复成本，本轮未发现足以证明输出错的输入。 | 作为 P2/优化候选，不与格式/资源正确性问题混级。 |

## C. v25 实施后的回归核对

2026-11-05 当前 working tree 中，v25 相关文本 GT、行分类、filter 路由改动已有对应测试；已运行的 `reed-io` / `reed-filter` / `reed-engine`、`reed` lib、`e2e_filter`、`e2e_convert`、`e2e_info`、`e2e_index`、CLI contract 均通过。完整 workspace/release/data-gate 未跑。工作区有用户未提交修改；本轮只写文档，不改这些代码。
