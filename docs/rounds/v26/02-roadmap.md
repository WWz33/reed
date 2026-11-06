# 02 v26 路线与验收门

## 1. 顺序

| 步 | 条目 | 前置 | 验收 |
|---|---|---|---|
| 0 | V26-00 mmap safety（需用户拍板） | 无 | 在第 4 节选择方案后实现；unsafe 条件有可审计实现，不只靠注释 |
| 1 | V26-01 BGZF block 验证 | V26-00 | streaming/mmap/region 对畸形 ISIZE、坏 CRC、假空块给同类错误；有效拼接空成员不丢；坏输入不先分配异常大小 buffer |
| 2 | V26-02/04/05/06 INFO/VCF/BCF 记录契约 | 无 | checked AC 溢出；重复 INFO 拒绝；必需空字段拒绝；零 allele BCF 三入口拒绝；text/BCF 与 noodles `.119` 的差分测试 |
| 3 | V26-03/21 thread domain 与 solver parallelism | 无 | CLI 对经用户确认的上限返回 usage；faer dense 路径遵守 `ThreadBudget`，1/2/N 线程实测上限且数值结果不变 |
| 4 | V26-07/09/20 resource accounting | 0, 3 | 长 seam/BCF carry 在设定预算下被限制；cgroup v1/v2 与 host RAM 取小的单测；预算含义和过程峰值界限写清 |
| 5 | V26-08 atomic outputs | 无 | 写入中途注入错误：已有文件保持原样；新输出不留下假成功文件；多文件报告有明确 partial-commit 策略 |
| 6 | V26-10 noodles boundary/update（拍板后） | 2 | 先隐藏 `Header` 跨 crate 暴露并删重复 dev dependency；若确认升级，再锁定 `.119.0`，workspace tests、BCF 编码、CSI/TBI、header fixtures 全过 |
| 7 | V26-11/12/14 CI 与测试稳定 | 无 | data-gate 的 reed 结果由当前 checkout 重建；MSRV job 实际运行；spill 全部测试在默认并行模式重复运行稳定 |
| 8 | Verify-only：R8/R9 与统计语义文档 | 前置实现之外 | R8 flag 矩阵对照 CLI spec；R9 重复 POS / r² 边界对照 plink2；只对复现差异建 bug |
| 9 | P2 清理：`budget::admit`、stats `site_terms`、distance/kNN、V26-18 长流优化与剩余拆分 | 0–8 | V26-18 仅在目标负载确认需要时实施；其余先列调用面，再做纯行为等价重构；不以“文件大”单独拆分 |

## 2. 实施纪律

1. 先加能失败的输入测试，再修共享底层契约；不要在每个调用点继续打 guard。
2. 处理 BGZF 时保持空 BGZF member 与 EOF marker 的兼容；不要把合法的拼接边界误当文件结束。
3. 原子输出按单文件/多文件分别定义承诺；多文件不能声称事务性，除非实现了回滚或目录级替换。
4. noodles 升级与边界收敛分开提交：先做 `.119` 兼容测试，再隐藏类型，便于定位 API 与行为变化。升级版本按 §4 待拍板项执行。
5. v26 当前只立项，不代表这些问题已修。

## 3. 本轮不作为缺陷处理

- `--hwe midp` 当前拒绝：已有 spec 与明确错误信息，是尚未实现的功能。
- `N_MISSING` 统计 fully-missing genotype、`-g miss` 将 half-call 视作 missing：属不同语义；文档需清楚，不要强行合成一个计数。
- 单纯大文件拆分、公式归一和 noodles 整体重写：除非先证明错误或边界收益，否则排在正确性修复之后。

## 4. 需要用户拍板的选择（含建议，不预先选定）

1. **V26-00 mmap 默认策略**：A. 建议默认 streaming，之后可给真正不可变快照提供 opt-in mmap；B. 保留 mmap 默认，但把“输入在整个映射生命周期内不可被任何进程修改”列为硬运行前提。memmap2 的安全合同使“只读打开”不足以单独证明 B 安全。
2. **输出原子性承诺**：A. 建议每个文件先写同目录临时文件、完成后 rename；多文件报告只承诺逐文件原子，并说明失败可出现部分新结果；B. 若需要整组 all-or-nothing，再单独设计 staging directory/回滚。
3. **`--threads` 上限**：建议上限按正式支持的最大 worker 数/集群部署约束定，并在 CLI 文档写明；请确认产品希望支持的最大值，不由本轮任意指定常数。
4. **noodles 升级**：建议在保留现有 BCF/CSI/TBI 实测门的前提下 pin 到本地 `noodles-0.119.0`，并在同一独立提交隐藏 `Header` 类型边界；如需先观察兼容性，可先单做升级验证，不改实现 API。
5. **超长记录与 `--memory`**：建议将 carry 纳入预算，同时保留硬记录上限；请确认 `--memory` 要继续定义为“主要中间数据预算”还是提升为包括 carry 在内的硬上限。当前文档只说它不是进程峰值 RSS 上限。
