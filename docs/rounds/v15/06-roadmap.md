# doc_v15 / 06 路线：F1–F5

> 顺序原则：**先落一行级、零语义争议的修复**（F1），再动统计口径（F2，
> 每处都要 oracle/allel 对拍兜底），再动格式与窗口边界（F3，字节级
> 对拍），规模项（F4）最后（有性能门），收口（F5）。每批完成即更新
> STATUS.md 并跑本批门。

| 批次 | 内容 | 门 |
|---|---|---|
| **F1** 机械修复 | D15-1、D15-5、D15-21、D15-22、D15-27、D15-30、D15-33、D15-13、D15-14、D15-39 | 每条的回归测试修复前红、修复后绿；`cargo test -p reed -p reed-cli`；cli 帮助/契约不变 |
| **F2** 统计口径 | D15-2、D15-11、D15-12、D15-31、D15-18、D15-9、D15-10 | reed-stats/oracle 对拍门 + e2e；`-q/-Q` 修复后 `precac.vcf` 形必须丢弃；Hudson/PBS 输出 0.0；oracle_parity 新增 F2/F3/F4 与零调用组用例 |
| **F3** 格式与窗口 | D15-3、D15-4、D15-6、D15-8、D15-15、D15-19、D15-20、D15-24、D15-25、D15-26、D15-29 | 手工 BGZF 夹具（boundary/idx/seam）与 bcftools 逐字节/逐记录对照；窗口和 = contig 行；convert 组合运行 = 单目标并集 |
| **F4** 规模 | D15-17、D15-28、D15-16（预防性）、D15-7（架构级，预算不足则降级登记） | 性能门（D15-17 单线程 ≥5×；D15-28 ≤2× 对照——实测收在 **2.5×**（3M `-Oz`，对照 0.091s），未达 2×，原因与残留见 STATUS F4 行/02 §D15-28）；RSS 门（D15-7 输出早现） |
| **F5** 收口 | 死代码（D15-36）、merge-law 指向（D15-37）、missing 字段（D15-38）、abort 加锁（D15-32）、两处分配（D15-35 已标点）、cli.md/help/CHANGELOG/docs 索引 | `cargo test --workspace` 全绿；clippy 零告警；gold 脚本不回归（run_split_parity/run_convert_parity/run_hist/run_allel_fstats 可跑部分） |

## 本轮不做

- convert 的完整流式重写若超预算 → D15-7 降级为单目标 spill（STATUS 记账）；
- D15-35 未标点的分配优化（render/native_encode 全量缓冲池化）留给后续轮；其中 `line_edit.rs:format_g6_into` 的零分配已做（2026-10-19 收尾，242.9→182.9 ns/call）；
- D15-36 中涉及 API 公开面的删除（`WindowGrid::bin_of_pos0` 在 reed-core
  公开面上）——只删无调用方且无文档承诺者，其余登记；
- 新功能、命令面调整、v14 结构项重开一律不做。

## 风险

| 风险 | 缓解 |
|---|---|
| D15-19/D15-20 改窗口水位线逻辑，可能引入新的边界回归 | 先固化窗口和 = contig 行的算术门（含重叠 step<size 用例），再改 |
| D15-3 改 header 扫描的跨块读，可能影响 BCF 首记录寻址 | 两侧夹具（H=65279/65280/65281/130560）+ 既有 BCF e2e 全绿 |
| D15-2 改动使既有一处 e2e 期望变化（若有用例依赖陈旧 INFO） | 先跑 `cargo test -p reed --test e2e_filter`，任何变化逐条判定 |
| D15-39 若 writer 不支持 stdout，改"显式拒绝"会破坏 `-o -` 的既有用户 | 查 README/cli.md 是否承诺 stdout；未承诺则拒绝并在 help 写明 |
| D15-8 输出形态改动影响 gold compare_pbs.py | 先读解析器再钉格式，改后跑 gold 脚本 |