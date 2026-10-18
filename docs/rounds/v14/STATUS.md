# doc_v14 实施状态

> 本页记录 `06-roadmap.md` C0–C5 的实施状态与证据。本轮是**命令面
> 重组 + 新功能轮**（stats/popg 拆分、分组 help、`--hist`、`rename`），
> 无缺陷修复清单；复审对象是"拆分等价"（0.3.x 输出固化对拍）与
> "实现 vs 参照工具行为"（bcftools stats/reheader/annotate、vcftools
> --singletons）。每步完成即更新。

## 实施状态（C0–C5）

| 项 | 内容 | 状态 | 证据 |
|---|---|---|---|
| C0-1 | 0.3.x stats 输出真值固化（run_split_parity.sh 骨架 + 期望树） | **完成** | 参照二进制 `target/gold/v14/reed-0.3.0-ref`（sha256 `6fe7615a…`，`reed 0.3.0 v0.3.0-93-g4dcbd24`）；22 case 真值树 `target/gold/v14/ref-0.3.0/`（126 文件）；`gold/popgen/run_split_parity.sh record|check` 两模式；check 对 0.3.x：QC 侧全绿、popg 侧如预期红（门有效） |
| C0-2 | 基线墙钟（QC/popg/rename 负载） | **部分**（QC/popg 已测，rename 待 C4） | clean.vcf.gz 3,315,604 记录/348 样本/56 线程：QC 全报告 **3.27 s**（user 2m22s）；popg 全统计（π/dxy/θ/Tajima/Fst/PBS @100k 窗）**2.12 s**（user 1m22s） |
| C0-3 | bcftools 参照真值（直方图节 + reheader/annotate 输出 + 差异登记） | 未开始 | — |
| C1-a | PopgCli + Command::Popg + removed-flag 重定向 | **完成** | `cli/mod.rs` 新增 `PopgCli`（窗口族全 flag）；`StatsCli` 只留 QC + 共享组 + 隐藏重定向字段；`translate/stats.rs::moved_flag` 报错含 `reed popg` 改写 |
| C1-b | translate 拆分 + 库门面改名（run_popgen/run_stats） | **完成** | 新 `translate/{common,popg}.rs` + `stats.rs` 重写；库改名为 `run_popgen`/`run_stats`、`PopgenRequest`/`PopgenReport`/`StatsReportRequest`/`StatsOutcome`；`cis` 守卫同步 |
| C1-c | run_split_parity.sh 转绿（拆分等价） | **完成** | `split parity: PASS (22 cases, byte-equal)`；参照二进制与真值树在 `target/gold/v14/` |
| C1-d | 既有 gold 脚本命令词替换后全绿 | **完成** | `run_convert_parity.sh`、`run_easysfs.sh` 全绿（7 脚本窗口调用改 `popg`） |
| C2-a | GROUPS 常量 + help 模板 + 渲染函数 | **完成** | `GROUPS` + `render_grouped_commands()` + `main.rs` 用 `after_help` 注入；五组渲染实测 |
| C2-b | 分组守卫（变体恰好一组 + help 快照） | **完成** | `groups_cover_every_command_exactly_once` + `grouped_help_lists_the_five_groups` 单测；drift 守卫同步 |
| C3-a | HistPlan + 内建 metric 渲染（qual/dp/dp-site/idl/af） | **完成** | 新 `crates/reed/src/stats/writers/histogram.rs`；内建 metric 复用一个 pass 的 QC 累加器 |
| C3-b | GQ + 通用 INFO/FMT 标签 | **完成** | `--hist gq`（FMT/GQ 行级解码）+ `INFO/<TAG>`/`FMT/<TAG>`（整数稠密 / 浮点 0.1 分箱），header 校验在 translate 前置 |
| C3-c | `--hist` 对拍 + threads {1,4} + bench | **完成** | `gold/popgen/run_hist.sh` 对拍 bcftools 的 QUAL/DP（逐基因型+逐位点）/IDD 四节，全 `SAME`（`hist parity: PASS`）；`e2e_hist.rs` 5 测试含确定性门 |
| C3-d | stats --singletons（C14-3）+ stats `--pop` 分层 QC（2026-10-17 追加进轮） | **完成** | `--singletons` 逐字节对拍 vcftools（fixture + `--multi` 50k）；`--pop`：`PREFIX.afreq.pop`/`.vmiss.pop` 每群体一行（首现序）、`.singletons` 加 POP 列（`.`=不在群体文件）、群体文件/`--force-samples` 规则与 popg 同源（`resolve_groups`）；不给 `--pop` 输出与 0.3.x 逐字节相同；`gold/popgen/run_popqc.sh` fixture PASS（分层值 == 按群体切样本后的 plink2 `--freq/--missing` 值；多等位 ALT 频率为折叠单值）；subset（348 样本/3 群体/1 GB）门与 e2e（`e2e_hist.rs` 4 个新测试：分层正确性、force 规则、cohort 不变性、threads 1v4）全绿 |
| C4-a | RenameCli + 门规则 + 分发 + 分组接入 | **完成** | `cli/mod.rs::RenameCli`（`-h` 让位给 `--header`，`disable_help_flag`）+ `translate/rename.rs` + `Command::Rename`（alias `reheader`） |
| C4-b | VCF 字节编辑四平面 | **完成** | `crates/reed/src/rename/mod.rs`（header/染色体/ID/样本四平面 + 求值序 + 校验表） |
| C4-b' | BCF 头块重写 + 字典重映射 + 快路径 | **未实现**（BCF 输入明确报用法错；登记后续） | — |
| C4-c | 对拍 + e2e + 确定性 | **完成** | `crates/reed/tests/e2e_rename.rs` 8 测试；手工对拍 `bcftools reheader -s`、`annotate --rename-chrs`、`annotate -I`（记录段逐字节 SAME） |
| C4-d | 索引失效 warning + -W | **完成**（bench 行待补） | `rename -W[=csi\|tbi]` 写出按改名后 contig 名的索引（bcftools `view -r` 可查）；输出旁陈旧索引打印 warning 指向 `reed index -f` |
| C5-a | 库 API 改名复核 + README Library use 编译门 | **完成** | 门面根导出 `run_popgen`/`run_stats`/`run_rename`；README 双语 Library use 段已改（编译门随 CI） |
| C5-b | 文档收口 | **完成**（全量 bench 行见 C5-c） | cli.md §1/§7/§12（含 --hist/--singletons/rename -W）；CHANGELOG 0.4.0；README 双语命令面重写；architecture §13.14 + statistics-map 已同步；守卫 87→95 全绿 |
| C5-c | 全量确定性 + 性能终门 | **部分**（功能门全绿；全量 bench 行待补） | `cargo test --workspace` 57 组全绿；`run_split_parity.sh` 22/22 |

## 待拍板（C14，默认已定可推翻）

| 项 | 默认 |
|---|---|
| C14-1 `popg` 主名 + `popgen` visible alias | **已拍板**（2026-10-17） |
| C14-2 裸 `reed stats`（无统计 flag）：报错（默认）or 默认 --summary | 已定默认：报错 |
| C14-3 `--singletons` 进首发轮（默认进）；`--geno-counts` backlog | 已定默认 |
| C14-4 `rename --order-samples` 不进首发轮 | 已定默认 |
| C14-5 `--hist` 浮点 bin 固定边界（bcftools 同款）；参数化 backlog | 已定默认 |
| C14-6 `rename` 后旧索引：warning + 提示 `reed index`（默认重建否） | 已定默认 |

## 差异登记落地状态（B14 → cli.md 差异表，随 C5-b）

| 项 | 状态 |
|---|---|
| B14-1 `--hist` 通用标签（bcftools 无） | 未落地 |
| B14-2 精确 mean/median/p95（bcftools 固定 bin 无分位） | 未落地 |
| B14-3 独立 `.hist` TSV（不嵌 .stats） | 未落地 |
| B14-4 `--chr-style` 有限三态枚举 | 未落地 |
| B14-5 `--set-id` bcftools `%` 模板拼写 + 多等位 ALT 逗号连接 | 未落地 |
| B14-6 `--rm-dup-ids` 加后缀保全部（plink2 丢弃） | 未落地 |
| B14-7 BCF 染色体/样本改名的字典重写 + 纯样本快路径 | 未落地 |

## 测试与守卫（落地后实测）

| 命令 | 结果 |
|---|---|
| `cargo test --workspace` | **全绿**（57 个测试组；含新 `e2e_hist` 8 + `e2e_rename` 14） |
| `cargo clippy --workspace --all-targets -- -D warnings` | **exit 0**（无告警） |
| `bash ci/check-doc-drift.sh` | **全绿**（95 → 96 项；修 4 条改名守卫 + 新增 9 条 v14 守卫） |
| `bash gold/popgen/run_split_parity.sh check` | **PASS**（22 case 逐字节，参照 0.3.x 二进制） |
| `bash gold/popgen/run_convert_parity.sh` | **PASS**（含 `convert --sfs` ≡ `popg --sfs`） |
| `bash gold/popgen/run_easysfs.sh` | **PASS**（rtol 1e-6） |
| `bash gold/popgen/run_hist.sh` | **PASS**（`--hist qual,dp,dp-site,idl` 对拍 bcftools 四节全 SAME） |
| `bash gold/popgen/run_singletons.sh` | **PASS**（`--singletons` 与 vcftools 逐字节） |
| `bash gold/popgen/run_singletons.sh --multi` | **PASS**（50k 记录 / ~13 chunk 的跨 chunk 聚合门） |
| `bash gold/popgen/run_popqc.sh` | **PASS**（`--pop` 分层值 vs 按群体切样本的 plink2 `--freq/--missing`；cohort 文件不变；singletons POP 列） |
| **subset 门**（`subset.Chr01-03.vcf.gz`：1 GB / 348 样本 / ~480k 记录 / Chr01-03，THREADS=8） | `run_split_parity.sh record+check` **PASS（22 case 逐字节）**；`run_hist.sh` **PASS**（QUAL 精确区逐 bin + 加宽尾按总量、DP 双维、IDD）；singletons 0 行（bcftools `+fill-tags AC==1` 独立核验确为 0，reed 与 vcftools 一致）；`run_popqc.sh` **PASS**（pop1/pop2/pop3 分层 afreq/vmiss 逐值等于按群体切样本的 plink2 跑法，多等位折叠差已归一）；分层三文件 threads 1 vs 8 `cmp` 全 SAME |
| subset 线程不变性抽查 | popg 五表 + QC 十文件 threads 1 vs 8 `cmp` 全 SAME |
| rename 对拍（手工） | `bcftools reheader -s` / `annotate --rename-chrs` / `annotate -I` 记录段逐字节 SAME |

## 基线（2026-10-17 实测，`clean.vcf.gz` = 8 GB / 348 样本 / 3,315,604 记录，56 核）

| 负载 | 0.3.x stats | 0.4.0 | 门 |
|---|---|---|---|
| QC 全报告（六 flag） | **3.27 s**（user 2m22s） | 待测 | 退化 ≤ 3% |
| QC + `--hist qual,dp,dp-site,gq,idl,af` | — | 待测 | 增幅 ≤ 10% |
| popg 全统计（π/dxy/θ/Tajima/Fst/PBS） | **2.12 s**（user 1m22s） | 待测 | 退化 ≤ 3% |
| rename 全平面 | —（bcftools 三进程基线待测） | 待测 | ≤ 1/3 bcftools 合计 |

## 对拍覆盖面注记（fixture vs subset）

fixture 只有个位数记录 = **单 chunk**，因此到达不了 `HistAcc::merge`、singletons 的
sink 聚合与 rename 的跨 chunk 去重；上述 subset 门与 `--multi` 门就是为覆盖这些路径加的。
subset 上首跑出的三处 `DIFF` 全部是标签/换算层面的差异（bcftools 溢出桶写作 `>N`、
QUAL 加宽区打印浮点无法还原 bin 边界），计数与总计逐项相同；脚本已把溢出桶归一化、
QUAL 改为"精确区逐 bin + 加宽尾按总量"比较后转绿。`--singletons` 在 subset 上为
0 行是数据事实（独立核验 AC==1 计数为 0），6 行/50,003 行的非空对照分别在
fixture 与 `--multi` 门上。

## 本轮未做（有意）

见 `06-roadmap.md` §4：pca/relatedness、norm/concat/isec/merge/query、
genautocomplete、filter 的 --extract/--write-snplist、stats --sfs 的
显式删除（重定向自然消亡）、rename --order-samples、--hist-bins、
stats/popg 语义调整、嵌套子命令。另登记后续：`rename` 的 BCF 输入/
输出与 `--order-samples`、全量 bench 行（功能门已全绿）。`stats --pop`
分层 QC 已于 2026-10-17 追加完成（C3-d）；`--geno-counts` 仍在 backlog。
