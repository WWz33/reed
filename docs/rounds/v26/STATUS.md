# v26 状态（2026-11-05 立项；2026-11-06 实施完成）

> 本页跟踪 v26 文档项。状态只在代码改动与验收证据齐备后更新；审查发现本身不计"完成"。

| ID | 初审状态 | 实施状态 | 证据 / 验收 |
|---|---|---|---|
| V26-00 mmap safety | confirmed; 方案待拍板 | **已修** | 默认输入路径改为流式（`reed/src/input.rs::open_bgzf`）；`MmapChunkSource::open` 的 `# Safety` 改为调用方契约（要求不可变快照），类型文档写明默认路径不再使用它。实测 7 GB `clean.vcf.gz`：流式 100.8 s / mmap 101.1 s，去掉 `##reed_command` 后**逐字节相同**。采用 roadmap §4 推荐方案 A。 |
| V26-01 BGZF footer/member 验证 | confirmed | **已修** | 新增 `bgzf::block::block_entry`（唯一入口）：校验 BSIZE、`ISIZE <= BGZF_INFLATE_MAX`，空成员用零长窗口让 libdeflate 证明"确实为空"（同时验 CRC）。stream/mmap/region 三源全部走它；`inflate_chunk` 分配前复查表、空成员再证明一次。测试：ISIZE=0 但 payload 非空、ISIZE 越界、坏 CRC、合法空成员、三源一致。 |
| V26-02 AC overflow | confirmed | **已修** | `info_af_outcome` checked 累加，溢出返回 `AcOverflow`；AC 元素数在逐元素解析前先按 `n_alt` 计数拒绝。测试含 `AC=18446744073709551615,1` 与 `u64::MAX` 边界。 |
| V26-03 thread upper bound | confirmed | **已修** | `MAX_THREADS = 1024` 解析期拒绝；`tickets()` 饱和乘；`abort_channel_for` 饱和加。CLI contract 测试覆盖 7 个子命令。采用 roadmap §4 推荐（绝对上限，机读可测）。 |
| V26-04 duplicate INFO | confirmed | **已修** | `validate_info` 拒绝重复标签，`classify_line` 统一入口；`InfoIndex::build` 返回类型化错误；AC/AN reader 改走同一 index（htslib `bcf_get_info_id` 也是首项）。测试覆盖重复 key / 重复 Flag。 |
| V26-05 BCF allele count zero | confirmed | **已修** | `site_header` 与 `render_site_part` 都拒绝 `n_allele = 0`，返回 `MalformedRecord`；renderer 的 `.` 特判删除。GT/render 两条入口测试。 |
| V26-06 empty record fields | confirmed | **已修** | `MalformedLine::Empty(EmptyField)`；空 CHROM/REF/ALT/QUAL/FILTER/INFO/FORMAT 拒绝，`.` 保持合法；判空在样本列一致性之前（否则真缺陷被列数错误掩盖）。 |
| V26-07 carry memory accounting | confirmed | **已修** | `record_ceiling(budget)` = min(1 GiB, `--memory`)；`NewlineStitcher::with_ceiling`；BCF carry 的**两处"先 append 后检查"改为 append 前检查**（原实现可越过 1 GiB，且完全无视预算）。测试：预算 8 B 的 stitcher 拒绝、`record_ceiling` 边界。 |
| V26-08 output atomicity | confirmed | **已修** | 新增 `reed-io::write::atomic::{AtomicFile, write_atomic}`：同目录临时文件 + `create_new` + `sync_all` + `rename`；未提交则 Drop 删除；60 s 以上的陈旧临时文件下次运行清理。替换面：filter/rename/convert `--vcf`（`OrderedOutput::create`）、CSI/TBI、convert 文本格式、报告 sink 的 7 个流式文件、`--ld-pairs`（仍惰性）、PCA 三件、DAPC 四件、SFS、distance 矩阵/表/`.id`、CLI 的 `--report`/`--timing`/`.log`/`.nwk`。实测 `kill -9` 后旧文件逐字节不变；全部命令端到端无临时残留。缺目录仍报错（与旧行为一致）。**多文件组只承诺逐文件原子**，已写入 `docs/spec/cli.md` §12。 |
| V26-09 cgroup-aware default memory | confirmed | **已修** | 探测取 min(host `MemTotal`, cgroup v2 `memory.max`, cgroup v1 `memory.limit_in_bytes`)；`max`/v1 哨兵/0 忽略；路径可注入，夹具树单测覆盖 6 种组合。 |
| V26-10 noodles boundary and 0.119 | confirmed | **已修**（两个独立提交） | ① 升 `noodles 0.119.0`（bcf 0.94 / vcf 0.96），1046 项全绿、无源码改动；删 `reed-io` 重复 dev-dependency。② `HeaderFacts.header` 改为 opaque `VcfHeader`：只暴露 sample/contig/INFO/FORMAT 访问器、`with_samples_pruned`/`with_info_definitions`/`with_pass_filter`/`with_meta_line`/`to_text`/`to_bcf_file_header`/encoder 构造。`reed` 不再直接依赖 noodles（仅 dev-dependency 供集成测试读回索引）。`api_surface` 增加编译门。 |
| V26-11 MSRV CI | confirmed gap | **已修** | CI 新增 `msrv` job（1.91）；本地 `cargo +1.91 check --workspace --all-targets` **实测通过**（1.91.1）。 |
| V26-12 spill test concurrency | confirmed / reproduced | **已修** | 根因是测试预测进程级 `SPILL_COUNTER`。新增 `NameSource` trait + `create_spill_file_from`；三个 SEC-1 测试各持私有 `TestCounter`，测试锁删除。`cargo test -p reed --lib` 默认并行连续 5 次通过。 |
| V26-13 unbounded header accumulation | confirmed | **已修** | `MAX_HEADER_BYTES = 256 MiB`；plain 与 BGZF 头扫描都设上限并报明确错误。测试：无 `#CHROM` 的超长头、4 MiB 合法大头、边界算术。 |
| V26-14 gold parity producer | confirmed | **已修** | `run_all.sh` 在 `REED` 存在时并入 reed 分支；CI data-gates 构建当前 checkout 的二进制、导出 `REED`、显式跑 `run_reed.sh`、断言产物存在后再 `check_parity.py`。区域烟测：reed 分支与参照工具同轮运行并产出 `out/reed.vcf.gz`。 |
| V26-15 R8 popg route | verify | **已核对：一致，非缺陷** | flag 矩阵实测：所有写文件的路由（`--pi/--dxy/--theta/--tajima-d/--fst/--pbs/--f2…/--sfs/--ld-decay/--ld-pairs` + 窗口模式）缺 `-o` 都是 usage error 2；`-o -` 被拒绝并给出 PREFIX 提示；`--sfs`/`--ld-*` 的 bare 路由是唯一的 stdout 路径。 |
| V26-16 R9 LD pruning parity | verify | **已核对：一致，非缺陷** | 构造含 8 个重复 POS 的 40 位点 fixture：reed 与 plink2 的 keep set 都是同一位点（`rs2`）。两等位全反相 fixture：`UNPHASED_R2` 两侧都是 `1`。 |
| V26-17 R13 atomic writes | confirmed | **并入 V26-08 已修** | 见 V26-08。 |
| V26-18 FrameChain metadata growth | confirmed scalability | **已测量，暂不优化** | `Option<Carried>` = 40 B/chunk（`Vec` 的空指针作 niche，无额外字节）。1 TiB BCF（1 MiB chunk）≈ 1.02M chunks ≈ 39 MiB cells + 0.1 MiB waiter 位图。测试钉住该成本。ring 化只对远超当前支持规模的输入有收益，故保留。 |
| V26-20 PCA dense peak undercount | confirmed | **已修** | 估算改为 `4·n·m·8 + 6·n²·8`（pairwise +1 份 n²），覆盖 `x`(2)、`st.w`、faer `wm`、GRM、`am`、`u`、EVD workspace、返回副本。用 `ru_maxrss` 实测拟合：1000/2000/3000 样本 × 2000 位点为 100.0/301.9/627.9 MiB，估计 106.8/305.2/624.0 MiB，误差 ±7%。`auto` 改用完整峰值判断。实测 gate：`--memory 550` 拒绝、`--memory 700` 通过。 |
| V26-21 faer bypasses thread budget | confirmed | **已修** | `reed-reduce::Parallel`（`Serial`/`Threads(n)`）贯穿 `gram_divide`/`sym_eig`/`sym_eig_top_k`、PCA 三入口与 `StreamConfig`；改为 faer 显式 `matmul`/`evd`（原 `*` 运算符读 faer 全局设置）。实测 2500 样本：修复前 `--threads 1` 耗 53.4 s CPU / 21.0 s wall；修复后 23.4 s CPU / 23.5 s wall。满秩输入在 1/2/8 线程下逐字节相同。 |
| v25 false-positive review notes | 已复核 | 关闭/更正 | 不变（HWE midp 是文档化限制；N_MISSING/half-call 是不同定义；reed-reduce 有测试）。 |

## 2026-11-06 复核发现（6 项）及处理

复核从 v25 基线到当前 HEAD 的 IO/VCF/BCF/filter/预算/输出/CI 改动后确认的问题。
逐条修复，每条都先补能失败的测试。

| # | 问题 | 处理 | 证据 |
|---|---|---|---|
| R-1 | `MmapChunkSource::open` 是安全 `pub fn`，但前提要求映射期间文件不被任何进程修改 —— 调用者仍能从安全代码触发 UB | 改为 `pub unsafe fn`，`# Safety` 写明"只读打开不足以建立前提"；内部 `unsafe` 块的 SAFETY 引用该契约；测试调用点显式 `#[allow(unsafe_code)]` + SAFETY 注释 | 仍满足"唯一生产 unsafe 站点"守卫（97/97） |
| R-2 | BCF 首段 carry：`carried.is_empty()` 分支直接 `text[body_end..].to_vec()` 并发布，**完全没有上限检查**；只有"已存在 carry"的分支检查 | 该分支发布前调用 `check_carried(0, crossing, declared_len(...))`，用声明的总长作更紧的界 | 新增 2 个测试；**把修复回退后 2 个测试都失败**，恢复后通过 |
| R-3 | 临时文件清理按"60 s 年龄 + 前缀"删除：会删掉共享前缀的用户文件，也会删掉另一进程仍在写的活动临时文件 | 只清理 `<目标名>.reed-tmp-<本进程 pid>-<纯数字>` 且超龄的文件；其他进程的临时文件一律不动（无法从文件名可靠证明 pid 已死，宁可留陈旧字节也不删） | 测试覆盖 6 类文件：本进程陈旧（删）、本进程新鲜（留）、他进程（留）、共享前缀用户文件（留）、非数字 pid 字段（留）、其他目标的临时文件（留） |
| R-4 | cgroup 探测只读根路径。若进程在子 cgroup 而根为 `max`，真实限额漏掉 | 解析 `/proc/self/cgroup`（v2 `0::/p`、v1 `N:memory:/p`），沿**自己的 cgroup 及其全部祖先**取最小；**并修正 v1 布局** —— 真实 v1 是 `<root>/memory/<rel>/...`（controller 在前），原实现猜的是 `<root>/<rel>/memory/...` | 新增 2 个测试（根无限/子有限、祖先更低、v1、sentinel、不可读回退）；并在本机实证 v1 路径 `/sys/fs/cgroup/memory/system.slice/sshd.service/...` 确实存在 |
| R-5 | `validate_info` 每条记录 `Vec` 分配 + `contains` 的 O(k²)；注释却写"allocation-free and linear" | 改为无分配的"与前缀重比"，注释如实写 O(k²) 及其理由（k 是个位数）；新增独立计数分配器测试 | 新测试：旧实现 **4000 次调用 6000 次分配**，新实现 **0** |
| R-6 | CI 的 data-gate 把 reed 分支跑两次（`run_all.sh` 已并入 + 又显式跑 `run_reed.sh`） | 删除显式的第二次调用，保留 `REED` 导出 + `run_all.sh` 一条路径 | YAML 中 `run_reed.sh` 只剩注释里的一处说明 |

## 2026-11-06 第二轮复核（3 项）

| # | 问题 | 处理 | 证据 |
|---|---|---|---|
| S-1 | 原子清理的**同进程**竞态：同一目标持有一个超过 60 s 未更新的临时文件，再创建该目标的另一个 `AtomicFile` 时仍可能被清理，导致提交失败 | 新增进程内活跃临时文件登记表 `LIVE_TEMPORARIES`（`Mutex<Vec<PathBuf>>`）。创建时在**同一把锁内**完成"清扫 + 登记"，提交/`Drop` 时注销；清扫跳过所有已登记路径。因此"崩溃遗留"和"活跃但安静"可以区分，不再靠年龄猜 | 新增 2 个测试：活跃兄弟临时文件（老化 1 h）在第二次 create 后仍存在且两个都能提交；注销后的陈旧同名文件仍会被清扫。**移除活跃检查后测试失败** |
| S-2 | INFO 查重仍是 O(k²)，标签数无上限 | 按"边界问题"处理，但把上界做成确定性的：新增 `MAX_INFO_TAGS = 4096`，`validate_info` 与 `InfoIndex::build` 在同一条规则上共享（超过即 malformed，且在平方扫描**之前**检查）。真实数据实测最大 7 个标签，余量三个数量级 | 测试：4096 个标签通过、4097 被拒、索引同样拒绝且失败后不遗留表 |
| S-3 | `--threads` 的 1024 上限未进 CLI 规格 | 已写进 `docs/spec/cli.md` §2 的 `--threads` 行（含 0/超限都是 exit 2、错误信息含上限、上限存在理由、远高于受支持部署）；新增 doc-drift 守卫把规格文本与 `MAX_THREADS` 常量绑在一起。顺带修正同一节里 `--solver auto` 仍写旧 `N²·8` 规则的过期描述 | 守卫 98/98；实测 `--threads 1025` 报 `exceeds the supported maximum of 1024`，1024 正常 |

## 本轮验收

- `ci/test-parallel.sh`：**1062 passed / 0 failed**（56 个二进制，17 s）。
- `cargo +stable clippy --workspace --all-targets -- -D warnings`：**0 告警**。
- `ci/check-doc-drift.sh`：**98 项 ok**；`ci/check-structure.sh`：**16 项 ok**。
- `cargo +1.91 check --workspace --all-targets`：通过（MSRV 实测）。
- 真实数据：7 GB `clean.vcf.gz`（3,315,604 记录）filter 结果记录数不变；8 GB `clean.bcf` 走 native BCF 路径通过；BCF→stats 往返通过。
- 未跑：全量 gold/data-gates（本机 plink2 在该子集上另有自身失败，与 reed 无关，故只做了区域烟测）；长时资源压力；完整性能基准回归。

## 复核后追加的偏差说明

R-3 的取舍值得单独记一笔：复核建议"不要跨进程按年龄删除，或至少严格验证文件名并确认所属进程已退出"。采用前半句的**收窄**版本 —— 只删本进程自己的。理由是"确认所属进程已退出"没有可移植且无竞态的实现（`kill(pid,0)` 的 PID 复用、/proc 可见性都依赖平台），而误删已提交数据的代价远大于留下陈旧临时文件。代价：另一进程崩溃留下的临时文件不会被自动清理，需要人工或后续维护。

## 实施中记录的两处方案偏差

1. **V26-13 上限取 256 MiB**（路线只说"设明示上限"）。依据：348 样本头约几 KiB，百万样本约 20 MiB，256 MiB 留足余量且仍能限制单次 open 的最坏分配。
2. **V26-18 不改 ring**。测量显示 1 TiB 输入才 ~39 MiB 元数据，低于任何当前支持负载的噪声水平；改 ring 会引入回收同步复杂度而无实测收益。

## 用户的待拍板项如何落地

roadmap §4 的五项在实施时采用了其推荐默认，并在对应条目与 `docs/spec/cli.md` §12 写明：

1. V26-00 → 默认流式，mmap 保留为显式 API（方案 A）。
2. 输出原子性 → 逐文件 temp+rename，多文件只承诺逐文件（方案 A）。
3. `--threads` 上限 → 绝对 1024（推荐项），未按核数派生。
4. noodles → 先升 0.119（独立提交），再隐藏 `Header` 边界（独立提交）。
5. `--memory` → carry 纳入预算（`record_ceiling`），同时保留 1 GiB 硬上限。

若这些默认需要改变，逐项改动面分别是：`input.rs::open_bgzf`、`write/atomic.rs` 的调用方、`budget.rs::MAX_THREADS`、`header/vcf_header.rs` 的公开面、`pipeline.rs::record_ceiling`。