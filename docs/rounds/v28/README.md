# doc_v28：reed-io 内部正确性与输出发布

> 立项 2026-11-10。承接 v27 的运行时审查，但只审 `reed-io` 内部及其必要调用方；
> 重点回答两个问题：还有哪些真实 Bug，哪些通用能力应下沉。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 逐条证据、修正后的定级、撤回项 |
| [02-roadmap.md](02-roadmap.md) | 批次、验收门、明确不做的事 |
| [STATUS.md](STATUS.md) | 实施状态与验收记录 |

## 总结

共 **28 条**：7 条 P0、5 条 P1、9 条 P2、7 条下沉设计（V28-28 是实施期间新发现的 P0）。

**2026-11-10 状态：V28-01…20 与 V28-21…27 全部实施；V28-28 为新发现的第六条 P0，已修复并交叉验证。**
逐条证据、改动和验收见 [STATUS.md](STATUS.md)。

- P0：BCF contig `IDX` 交换、BCF index tid 错位、region 短块静默 EOF、负 GT 变 REF、
  integer vector-end 泄漏为负数、完整性校验前发布旧输出，
  以及 V28-28：`filter -Ob -W` 的 CSI 参考表用行序而记录用字典序——htslib 对两个 contig 都返回 0 条。
- P1：VCF 4.4 首等位 phasing、`NULL,count=0` GT、half-call 策略快路径、零等位 encoder、CRLF 保真。
- 下沉：BCF shared/FORMAT 原语、contig 字典序、`AtomicWriter`、`LineRef` 区间访问器、
  BCF 跨 chunk framing、统一 input open、`VirtualOffset`/`BlockTable` 均已落地。
- 上轮 LD 嵌套线程池主张已撤回：当前代码只有单个 `LdExecutor` pool。
- 上轮 `GT=/0|1` 复现例子已更正：真差异是 `|0/1` 和 `|.`。

## 边界

- 本轮先修正确性和发布时序；下沉只做被修复直接依赖的部分。
- 逐文件原子，不新增多文件事务。
- 不重复 V27-16、V27-17 已落地的 LD 执行器修复。
