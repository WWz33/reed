# 文档索引（当前规范 / 历史）

> docs/rounds/v4/03-architecture.md §2 A4-7 定的规则：每轮文档并存且互相引用，
> 这一页说明**哪一份是当前规范、哪些是历史记录**，以及每轮被谁取代。
> 规则的默认方向：后一轮只修订它明确点到的章节，其余条款继续有效。
>
> 目录搬迁（docs/rounds/v5/06-layout.md）后：`docs/spec/` 是当前规范，
> `docs/rounds/vN/` 是历史轮次（只增不改）。

## 当前规范（按主题）

| 主题 | 当前文档 | 说明 |
|---|---|---|
| CLI 参数语义 | `docs/spec/cli.md` | **合并版 CLI 规范**（= docs/rounds/v5/02），当前唯一有效；v4/v2 的 02 转为历史 |
| 统计出口设计 | `docs/rounds/v5/05-algorithms.md` | SFS / `--fst-global` / f 统计 / LD 的出口形状；对齐台账在 v2/04 + v4/05 |
| 统计与过滤算法 | `docs/rounds/v2/04-algorithms.md` + `docs/rounds/v4/05-algorithms.md` | 方法与文献定义在 v2；对齐状态与缺口台账在 v4/05 |
| 架构 | `docs/spec/architecture.md`（v1 基础）+ `docs/rounds/v2/03-architecture.md` + `docs/rounds/v3/02-bcf.md` + `docs/rounds/v4/03-architecture.md` + `docs/rounds/v5/03-architecture.md` + **`docs/rounds/v6/03-architecture.md`** | v1 长期条款仍有效；BCF 在 v3/02；`unsafe` 策略与前后端边界在 v5/03；**缺陷清单 A6-1…A6-25 与三个最优方案（矩阵化 pass0、`-Ob` 记录模型、列式掩码）在 v6/03** |
| 模型（数据 / 执行 / 统计） | **`docs/rounds/v6/02-models.md`** | 基因型三种表示与 2 位编码的取舍、位平面矩阵作为全局分析的规范输入（含代码示意）、统计模型对拍与偏差台账、M6-1…M6-7 |
| 工程（CI / 测试 / 依赖 / 文档漂移） | **`docs/rounds/v6/05-engineering.md`** | 健康状态、CI 缺项、测试缺口、漂移清单、E6-1…E6-9 |
| 性能 | **`docs/rounds/v6/04-performance.md`**（当前基线）+ `docs/rounds/v5/04-performance.md` + `docs/rounds/v4/04-performance.md` | v6/04：3 次中位数基线、线程效率、deflate 级别研究、与 plink2/bcftools 对照、LD 执行模型问题、P6-1…P6-10；测量纪律在 v5/04 §4 + v6/04 §9 |
| 命名 | `docs/spec/naming.md` + `docs/rounds/v2/05-naming.md` | 版本后缀禁令（v2 05 §4）由 `ci/check-doc-drift.sh` 守卫 |
| 项目定位 | `docs/rounds/v5/08-api-surface.md` | 前后端形态（`reed-cli` 是唯一前端）与门面边界（§6） |
| 路线与验收门 | **`docs/rounds/v6/06-roadmap.md`** | Y0–Y7；v5/07 的 X0–X7 已落地（未完成项并入 Y6），v4/06 的 W0–W6 已完成 |
| 缺陷修复（外部审查核实） | **`docs/rounds/v7/01-review.md`** + `docs/rounds/v7/06-roadmap.md` | 对一份外部审查报告的独立核实与修复规范（D7-1…D7-29）；路线 Z0–Z4。v7 是缺陷修复轮，不重开 v6 的 Y 系列结构项 |
| 线程模型与并行预算 | **`docs/rounds/v8/01-review.md`** + `docs/rounds/v8/06-roadmap.md`** | 全部线程创建点台账与 spec §7.3 的逐条对照；rayon 三阶段与 stats LD sink 脱离 `ThreadBudget` 的修复规范（D8-1…D8-5）；路线 P0–P2。v8 是缺陷修复轮，不重开 v6/v7 的结构项 |
| 引擎错误路径与终止语义 | **`docs/rounds/v9/01-review.md`** + `docs/rounds/v9/02-architecture.md` + `docs/rounds/v9/06-roadmap.md` | 阻塞点—唤醒来源台账；reducer 侧错误/panic **永久挂死**的实测现场（`eu-stack` 5 线程全 `futex_wait`）与规约 W1/W2（唯一错误出口）；通道发送结果语义 C1；输入截断契约（文本丢弃 / index 收录 / BCF 报错）；seam 与 chunk 体的错误一致性 E1；指标口径（peak = 压缩 × `MEM_FACTOR`、reducer 上的 seam 时间漏计）；缺陷 D9-1…D9-8，路线 Q0–Q3。v9 是缺陷修复轮，不重开 v6/v7/v8 的结构项 |
| BGZF/BCF/VCF 格式兼容性 | **`docs/rounds/v10/01-review.md`** + `docs/rounds/v10/02-architecture.md` + `docs/rounds/v10/06-roadmap.md` | 以工作区 htslib/bcftools 源码为基准的格式处理对照：GT 选宽按编码值（N1，修大等位编号被 `as i16` 静默截断）；NULL 描述符非法（N2）；IDX 限界（N3）；BGZF 最小块长（N4）；文本 index 残片承载（N5，对齐 `BcfIndexJob` 的 FrameChain 模式）；多成员 gzip（N6）；PASS 注入（N7）；计数乘法 checked_mul（N8）。缺陷 D10-1…D10-8，路线 R0–R3。v10 是缺陷修复轮，不重开 v6–v9 的结构项 |
| Y5 引擎落地规格 | **`docs/spec/v6_y5/README.md`** + `docs/spec/v6_y5/TODO.md` | A6-1…A6-4 的落地形态（中止原子位 / 许可 drop guard / 接缝上限 / 顺序模式 / 计量）、与 roadmap 附录的差异表、门与执行清单 |
| Y6 前后端落地规格 | **`docs/spec/v6_y6/README.md`** + `docs/spec/v6_y6/CODE.md` + `docs/spec/v6_y6/TODO.md`；批 D 细化 `docs/spec/v6_y6/batch_d/`（README/TODO/CODE） | F6-13/F6-19 的落地形态（reed::Error / 请求内化 / job.rs 拆分 / C-5 / 联合 SFS / `--bed-file`）、与 roadmap 与 08 §6.2 的差异表、门与执行清单；批 D 的真值表/集成点/退出码契约 |
| 实施状态 | **`docs/rounds/v6/STATUS.md`**；缺陷修复轮各自的状态在 **`docs/rounds/v7/STATUS.md`**（Z0–Z4）、**`docs/rounds/v8/STATUS.md`**（P0–P2，含 D8-6）、**`docs/rounds/v9/STATUS.md`**（Q0–Q3，起始全部未修 + 复审基线与实测记录）、**`docs/rounds/v10/STATUS.md`**（R0–R3，起始全部未修 + 复审基线代码证据摘录） | 历史：`docs/rounds/v5/STATUS.md`（X0–X7）、`v4/STATUS.md`（W0–W4）、`v3/STATUS.md`（V0–V6）、`v2/STATUS.md`（M0–M5b/R0–R8） |

## 历史文档（只作记录，不再修订）

| 文档 | 状态 |
|---|---|
| `docs/rounds/v1/design.md` | v1 设计（原 `DESIGN.md`）；除后续轮次明确推翻的章节外仍然有效 |
| `docs/spec/architecture.md` | v1 规范；除 docs/rounds/v2/03 明确推翻的章节外仍然有效 |
| `docs/spec/filter.md` | reed-filter 目标架构（v2，2026-10-05）；执行模型已按 docs/rounds/v2/03 与 v3 落地 |
| `docs/spec/naming.md` | 命名规范；稳定性状态已按 docs/rounds/v5/08 §6 更新 |
| `docs/spec/statistics-map.md` | 统计量到参考实现的映射；对拍结果见各轮 STATUS |
| `docs/rounds/v2/` | 第二轮复审与修订规范；被 v3（BCF/性能/清理）与 v4（参数方向）增量修订 |
| `docs/rounds/v3/` | 第三轮复审与 BCF/性能/清理；V0–V6 已完成，状态在 `v3/STATUS.md` |
| `docs/rounds/v4/` | 第四轮：复审（01）、参数语义规范回退版（02）、架构残留（03）、性能（04）、算法台账（05）、路线 W0–W6（06）、目录树（07）、状态；W0–W4 已落地 |
| `docs/rounds/v5/` | 第五轮：复审（01）、合并版 CLI 规范（02）、架构与前后端边界（03）、性能基线（04）、W5/W6 出口设计（05）、目录树执行版（06）、路线 X0–X7（07）、项目定位（08）、状态 |

## 每轮明确取代的章节（增量，不是整份作废）

| 被取代的条款 | 取代者 |
|---|---|
| v2 02 §3.1 `--mind` / §3.4 `--geno`（plink2 丢弃方向） | v4 02 §3：**名字保留**，语义精确定义（等位计、掩码后、精确边界） |
| v4 02 首版 §5.1–§5.4、§6、§7 的改名（`--max-missing*`） | v4 02（回退版）同名小节：撤回改名 |
| v2 02 §5.2 `--fst-pops A,B,C`（多群体合并 Fst） | v4 05 §2.1 / v5 05 §3.2：作废，改为 `--fst` + 布尔开关 `--fst-global`（默认成对） |
| v2 02（整体）与 v4 02（回退版） | **docs/spec/cli.md**（v5 02 合并版）是当前有效的 CLI 规范 |
| v4 05 §2–§2.1、06 的 W5/W6、07-layout | v5 05（出口设计）、07（X0–X7）、06（目录树执行版） |
| `docs/spec/architecture.md` D12（Python 绑定 M8 之后再做） | v5 08 §6.5：**不做**；项目定调为前后端形态，`reed-cli` 是唯一前端 |
| v2 05 §2.3"`--max-missing` 方向与 `--geno` 相反，不提供别名" | v4 02 §8：结论保留（不提供别名）；`--max-missing` 已不是 reed 的参数 |
| v3 04 §C3"`TagTable` 未做" | `docs/rounds/v3/STATUS.md` 更正为**已做**（v4 01 F4-2） |
| `docs/spec/architecture.md` §13.17 的漏斗 id 示例 | v4 02 §5.3 定名（`after_mask.missing`、`sample.mind`）；§13.17 已就地回改 |
| `docs/spec/architecture.md` §7.7「超出时报 `EngineError::MemBudgetExceeded`」 | **doc_v9 D9-7（待落地）**：改为背压表述——超出时 reader 阻塞，`MemBudgetExceeded` 只用于"单个 chunk 的足迹超过全部预算"。依据是 `docs/rounds/v2/03-architecture.md` §4.2，实现早已如此，spec 未回改 |
| `docs/spec/architecture.md` §7.5「所有阻塞点都用 `select!` 等待数据和 abort」 | **doc_v9 02 §2（待落地）**：追加规约 W2——任何线程的错误与 panic 必须先经 `FirstError` 记录（记录即触发 abort），不得绕过它直接返回。§7.5 的机制本身正确，缺的是"没有旁路"这一条 |
| `v8/STATUS.md:22-23`「`stitcher_reassembles_and_caps_fragments` 失败是预存在的，与本轮无关」 | **doc_v9 D9-8 更正**：预存在属实，但"无关"是误判——它是 D9-3（`NewlineStitcher` 缺 finish 语义）的症状，且暴露了 release-only 验收门对 `debug_assert` 的盲区 |

## 一致性守卫

`bash ci/check-doc-drift.sh` 校验"STATUS 声明 vs 代码事实"的断言：C3 无 noodles、无版本后缀类型、
`--geno`/`--mind` 是真实参数、`--max-missing` 只出现在报错文本里、`reed index` 拒绝 stdin、
`Verbosity` 单一实现、恰两处 `unsafe`、无 `sample.max_missing`/`--fst-pops`/`reed qc`、
前端只依赖门面、后端不打 stderr、门面无 `pub mod`、README 无 v0.2 示例、无 noodles re-export。
文档声明"未做"的项应在该脚本里有一条守卫。
