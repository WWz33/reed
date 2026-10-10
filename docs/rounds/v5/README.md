# reed doc_v5：第五轮复审 + 合并版 CLI 规范 + W5/W6 出口设计（2026-10-09 晚）

> 延续 doc_v4（参数语义、W0–W6 路线）。doc_v4 的 W0–W4 已落地（记录在 `docs/rounds/v4/STATUS.md`），W5 在途、W6 未做。
> 本目录给出：落地后的**全项目复审**（01）；把 docs/rounds/v2/02 与 docs/rounds/v4/02 回退版**合成一份当前有效的 CLI 规范**（02）；架构与 `crates/` 残留（03）；新基线（04）；尚未接入 CLI 的统计量的**出口形状**（05）；目录树执行版（06）；路线 X0–X7（07）。
> 本目录不修改代码；实现者按 07 执行并把结果写进 `docs/rounds/v5/STATUS.md`。
>
> **后继修订**：2026-10-10 的第六轮全面审核在 `docs/rounds/v6/`。本目录的 04（性能基线）、07（路线）由 v6 的 04 / 06 取代；03 的残留项与 05 的 LD 执行模型由 v6 的 02 / 03 修订；X0–X7 的实施记录仍以本目录的 STATUS 为准。

## 本轮依据的用户裁决（2026-10-09）

1. `--geno`/`--mind` **保留原样**（plink2 短名，值 = 允许的缺失比例）；撤回 `--max-missing` 家族。→ 已在 W0b 落地。
2. 多群体 Fst：`--fst [wc84|hudson]` 不变，新增**布尔开关 `--fst-global`**（三个群体以上时算合并而非成对；默认成对）；`--fst-pops` 作废。→ 待实现（05 §3.2）。
3. 目录树：只针对 `reed-rs/`；按 06 的树改。→ 待执行（X2）。
4. 不为语义统一而改短名字；语义缺口用条文补。→ 写进 02 §0。
5. 不做 `reed-py`，`reed-cli` 是唯一前端；项目定调为前后端形态。→ 08 §6。
6. **LD 剪枝不照搬 plink 的两段式**（2026-10-10）：`--indep-pairwise` 是 `reed filter` 的位点规则（stats 共用），一条命令直接输出剪枝后的 VCF；`.prune.in/.prune.out` 只是 `--prune-list PREFIX` 的可选副产品；位点选择算法照 plink2 并逐 ID 对拍。→ 02 §3.4、05 §3.4、X4.4。
7. **本项目目前不做多等位**（2026-10-10）：X5 取消；现状即终态（默认 `-M 2` 剔除，放开后要求重算时报错不写错）；05 §4 的设计存档。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 复审：工作树与 HEAD 分别核实（工作树不编译、HEAD 364 测试全绿 0 警告）、STATUS 逐条对照、跨格式 md5、发现 F5-1…F5-11 |
| [02-cli.md](02-cli.md) | **合并版 CLI 规范**：三条命令全部参数（按实测 `--help` 核对）、缺失阈值的精确定义与整数判定、表达式子集、求值顺序、stats 的统计参数/窗口/输出列、`--fst-global`、`--sfs`、退出码、与 `--help` 的差异登记 |
| [03-architecture.md](03-architecture.md) | 建成形态（mmap、BCF 原生过滤输入、长度前缀链、票据）、`unsafe` 策略的两处例外、残余差距 A5-1…A5-5、`crates/` 内残留 C-1…C-9、reports 门的结构修法、守卫扩展 |
| [04-performance.md](04-performance.md) | v5 基线（3 次中位数）与 v4 对照、CPU 账、下一步 P5-1…P5-5、测量纪律修正（负载自抬升） |
| [05-algorithms.md](05-algorithms.md) | 对齐台账；**W5 出口设计**：`--sfs --proj`（genome-wide 路径、dadi/fsc2、easySFS 对拍）、`--fst-global`（vcftools weighted Fst 对拍）、f 统计（Patterson、block jackknife、allel 对拍）、LD（剪枝进 filter、尾窗携带、plink2 对拍，附代码示意）；W6 多等位重算已裁决不做，设计存档 |
| [06-layout.md](06-layout.md) | 目录树执行版（承 docs/rounds/v4/07）：目标树、映射与 `git mv`、要同步改的引用、注意事项 |
| [07-roadmap.md](07-roadmap.md) | X0 解除阻塞 → X1 卫生 → X2 搬迁 → X3 reports 门；X4 统计接入；X5 多等位；X6 CLI 缺口；X7 CI |
| [08-api-surface.md](08-api-surface.md) | **项目定位**：用代码事实回答"像不像 API 库"（结论：CLI 产品，库是内部结构，README 停在 v0.2）；两种终态 A（工具优先，推荐）/ B（API 库）；A 的清单 L-1…L-7 并入 X1 |
| [STATUS.md](STATUS.md) | 实施状态（实现者填写） |

## 结论速览

| 维度 | 结论 |
|---|---|
| 完成度 | doc_v4 W0–W4 全部落地；W5 只有前置调查 + 在途未提交代码；W6 未做 |
| 主干状态 | **工作树不能编译**（在途 W5 改动：`reed-cli` bin E0382、`reed-stats` 测试 E0004）；HEAD `c5b6bfe` 在独立 worktree 验证：364 测试全绿、0 警告、守卫 8/8 |
| 正确性 | stats VCF/BCF 输入输出 md5 相同（与 v4 记录一致）；gold 全量 md5 本轮未重跑，以 W4 记录为准 |
| 性能（56 线程，中位数） | filter VCF `-Oz6` 8.79 s（不回退）；**BCF 输入 8.87 s**（v4 12.1）；stats 五统计 **1.74 / 1.51 s**（v4 2.7 / 3.1）；index **1.60 / 1.05 s**；reports 3.76 s（门 3.2 未达）；边过滤 `--pi` 3.04 s（门边缘） |
| 文档与代码 | STATUS 一处不实（`unsafe` 两处而非一处）；三处注释/帮助文本过期；doc_v2 的三个参数未兑现 |
| 结构 | crate 分层干净；`crates/` 内有 9 条历史残留；目录树方案已给未执行 |

## 需要你确认

| 项 | 推荐 | 备选 |
|---|---|---|
| 在途 W5 代码的处置（X0） | 按 05 §3.1 改成 genome-wide 路径后再提交，期间 `git stash` 让主干可测 | 直接提交现有窗口式 `StatKind::ProjectedSfs` 并补 match 臂（会留下"分窗谱"这个无意义形态） |
| `--sfs` 的默认 | unfolded（REF 为祖先态，与 easySFS `--unfolded` 同）；写浮点，`--sfs-int` 可取整 | 默认 folded（easySFS 的默认） |
| f 统计参数语法 | `--f3 C;A,B`、`--f4 A,B;C,D`（分号分角色组） | 全逗号 `--f3 C,A,B`（靠位置记角色） |
| LD 的默认 | **已裁决**：`--indep-pairwise W[kb] [S] R2` 是过滤规则，一条命令出剪枝后 VCF，`S` 缺省 1（plink2 2.0 语义）；`--ld-decay` 默认 `100k`、1 kb 分箱留在 stats | —（05 §3.4） |
| 边过滤 `--pi` 的门 | 改为 CPU 口径 ≤ 130 s（wall 作参考） | 维持 wall ≤ 3.0 s 并继续优化 |
| 目录搬迁时机（X2） | 与 X1 合成一个"卫生"PR，在 X4 之前 | X4 之后 |
| 项目定位（08） | **已定调：前后端形态**（08 §6）——后端 = `reed` 门面 + 下层 crate，前端 = `reed-cli`（唯一前端），契约 = 02 的参数语法与输出格式；不发布 crate；前端只依赖门面、后端不打 stderr、门面只 `pub use` 三条规则进守卫 | —（`reed-py` 已裁决不做，08 §6.5） |
