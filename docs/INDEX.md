# 文档索引（当前规范 / 历史）

> doc_v4 03-architecture.md §2 A4-7：四轮文档（`docs/` 的 v1、doc_v2、doc_v3、doc_v4）并存且互相引用，
> 这一页说明**哪一份是当前规范、哪些是历史记录**，以及每轮被谁取代。
> 规则的默认方向：后一轮只修订它明确点到的章节，其余条款继续有效（doc_v4 02-cli-semantics.md §8 的差异登记是范例）。

## 当前规范（按主题）

| 主题 | 当前文档 | 说明 |
|---|---|---|
| CLI 参数语义 | `doc_v4/02-cli-semantics.md` + `doc_v2/02-cli.md` | v4 §1 的"方向规则"（一切阈值都是保留条件）优先；v4 未提到的章节仍按 doc_v2/02-cli.md |
| 统计与过滤算法 | `doc_v2/04-algorithms.md` + `doc_v4/05-algorithms.md` | 方法与文献定义在 doc_v2；对齐状态与缺口台账在 v4/05 |
| 架构 | `docs/architecture-spec.md`（v1 基础）+ `doc_v2/03-architecture.md` + `doc_v3/02-bcf.md` + `doc_v4/03-architecture.md` | v1 §13 的命名/输出列/提交纪律等长期条款仍有效；BCF 编解码在 doc_v3/02；当前形态与残留在 v4/03 |
| 性能 | `doc_v3/03-performance.md` + `doc_v4/04-performance.md` | 优化清单与测量纪律在 v4/04（含共享机器污染样本的处理） |
| 命名 | `docs/naming-conventions.md` + `doc_v2/05-naming.md` | 版本后缀禁令（doc_v2 05 §4）已由 `scripts/check-doc-drift.sh` 守卫 |
| 路线与验收门 | `doc_v4/06-roadmap.md` | W0–W6；旧的 doc_v2/06、doc_v3/05 均已完成 |
| 实施状态 | `doc_v4/STATUS.md` | 历史：`doc_v3/STATUS.md`（V0–V6）、`doc_v2/STATUS.md`（M0–M5b/R0–R8） |

## 历史文档（只作记录，不再修订）

| 文档 | 状态 |
|---|---|
| `docs/architecture-spec.md` | v1 规范；除 doc_v2/03 明确推翻的章节外仍然有效 |
| `docs/filter-architecture.md` | reed-filter 目标架构（v2，2026-10-05）；执行模型已按 doc_v2/03 与 doc_v3 落地 |
| `docs/statistics-reference-map.md` | 统计量到参考实现的映射；对拍结果见各轮 STATUS |
| `doc_v2/` | 第二轮复审与修订规范；被 doc_v3（BCF/性能/清理）与 doc_v4（参数方向）增量修订 |
| `doc_v3/` | 第三轮复审与 BCF/性能/清理；V0–V6 已完成，状态在 `doc_v3/STATUS.md`（C3 行按 doc_v4 01 F4-2 更正） |
| `doc_v4/` | **当前一轮**：复审（01）、参数语义规范（02）、架构残留（03）、性能（04）、算法台账（05）、路线（06）、状态（STATUS） |

## 每轮明确取代的章节（增量，不是整份作废）

| 被取代的条款 | 取代者 |
|---|---|
| doc_v2 02 §3.1 `--mind`（plink2 丢弃方向） | doc_v4 02 §3.2 `--max-missing-sample`（保留方向、等位计） |
| doc_v2 02 §3.4 `--geno`（plink2 丢弃方向） | doc_v4 02 §3.1 `--max-missing` + `--max-missing-count`（vcftools 方向） |
| doc_v2 05 §2.3"`--max-missing` 方向与 `--geno` 相反，不提供别名" | doc_v4 02 §8：反转，`--max-missing` 成为主名 |
| doc_v3 04 §C3"`TagTable` 未做" | doc_v3/STATUS.md 更正为**已做**（doc_v4 01 F4-2） |
| docs/architecture-spec.md §13.17 的漏斗 id 示例（`after_mask.missing`） | doc_v4 02 §5.3（`after_mask.max_missing` 等）；architecture-spec.md §13.17 已就地更新 |

## 一致性守卫

`bash scripts/check-doc-drift.sh` 校验若干"STATUS 声明 vs 代码事实"的断言（C3 无 noodles、无版本后缀类型、
`--max-missing` 家族存在、`reed index` 拒绝 stdin、`Verbosity` 单一实现）。文档声明"未做"的项应在该脚本里有一条守卫。
