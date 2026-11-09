# doc_v27：popg 统计运行时与输出契约审查（2026-11-08 立项）

> 本轮承接 doc_v26 与仓库外 `architecture-review.md`（基线 `7e37a05d77b803e875c4aeca535cf9c392ec14da`）：
> 逐条核对上一轮报告中的 5 条 popg 缺陷，并在三个独立方向（CLI/输出、统计算法、运行时与资源）重扫
> `reed popg` 全链路。只写**已实测复现或可由确定性代码路径证明**的条目；推测一律标 `verify`。
>
> **2026-11-09 更新：第一至十批全部落地**（V27-01…V27-39）。逐条状态、改动与反向验证见
> [STATUS.md](STATUS.md)。01/02 保留立项时的缺陷证据与最小处理方向；实施状态、验收证据和路线进展随批次更新。

范围：`reed-cli` 的 `popg` 前端、`reed/src/stats/`（driver/run/sink/writers）、`reed-stats` 的全部估计器与
`jackknife`/`matrix`/`window`、`reed-engine` 的 `pipeline`/`budget`、`reed-core` 的 `grid`/`count`。
不覆盖 filter/convert/pca/dapc/rename 主路径（除与 popg 共用 writer 或 `prelude` 的接缝）。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 逐条复核：ID、优先级、位置、触发与实测数字、最小处理方向 |
| [02-roadmap.md](02-roadmap.md) | 顺序、验收门、实施纪律、不作为缺陷的旧候选、最终定案与实施记录 |
| [STATUS.md](STATUS.md) | 实施状态（39 条全部落地） |

## 总结

共 **39 条**，其中 **3 条 P0**（可静默毁数据或给出错误数值）、**17 条 P1**、**19 条 P2**。

### P0（数据安全与数值错误）

1. **V27-01** 裸路由 `-o` 可等于输入路径，`exit 0` 静默覆盖输入文件（`refuse_output_eq_input` 存在但 popg 从未调用）。
2. **V27-02** 失败运行会**发布截断的 `.ld`**，覆盖上一轮完好文件——`sink.finish()` 的 commit 早于 `job.finish()` 的完整性校验，违反 doc_v26 V26-08「失败时旧文件逐字节不变」。
3. **V27-03** Hudson Fst / PBS 的分母在 `n < 2` 的位点仍累加，而分子不累加；窗口值被稀释（实测 0.2222 vs 文档/pixy 的 1/3，低覆盖夹具 −39.9%）。文档 04-algorithms.md:56 与 pixy 源码都要求「分子分母同时剔除」，**但现有测试把相反行为钉住了**。

### 上一轮 5 条：全部属实

原报告的 5 条（裸命令忽略 variant-only、LD 单元素尾块、零配对保留旧 `.ld`、`step > size` 误报乱序、全基因组统计强制 contig 长度）逐条实测复现，**数字全部对上**。
本轮对其中三条扩大了范围：

- 零配对问题不限于 `.ld`：`.ld.decay`、`_sfs/`、`.fstats` 在空结果时都保留旧文件（窗口类输出反而正确，因为 writer 无条件调用）。
- LD 单元素尾块是 **4096k+1 整族**（k=1/2/3 实测 4097/8193/12289 全部漏 1 对），不是只有 4097。
- 裸命令忽略 variant-only 的「100 倍」是夹具特定的；本质是 per-pair 口径 vs 完整分母。且裸路由同时产 `pi.bp`/`theta_w.bp`，修复时必须一并处理。

### 方向性结论（正面，避免后续重复排查）

- **分块等价性总体成立**：plain / bgzip `-l 1/2/4/8/9` / `--binary` × `--threads 1/2/3/5/8` × `--memory 8/32/1024`，7 组 × 60 配置下 π/dxy/θ/Tajima/Fst/fstats/contig/genome/LD 全部**逐字节相同**。唯二例外是 V27-06（`step > size`，第三批已修）与 V27-28（多等位绕过升序，第八批已修）。
- **有序 sink 的隐藏输出槽是安全的**：`out` 槽只在全部检查通过后写入，失败路径不残留部分 `SinkOutput`。但 `.ld` 的 commit 早于完整性校验（即 V27-02）。
- **SFS 投影与 easySFS 逐值相同**（1D/2D、folded/unfolded，最大相对差 8.6e-16）。
- **LD 的 `.ld.decay` 跨线程数 md5 一致**（`--threads 2/4/8/16`）。
- **`NeumaierSum` merge 在 30 个种子上逐位相同**。
- **π/dxy 的 `PairAcc` 三列自洽**：`missing = total_pairs·records − comparable` 在所有构造的边界夹具上成立。

### 已拍板（2026-11-08）：pixy 是主旋律

用户裁决：`docs/spec/cli.md` 与 pixy 不一致的地方以 pixy 为准（这与 `cli.md §7.1`
「缺失数据统一采用 Korunes & Samuk 2021（pixy）的逐位点计数法」的既有原则一致）。
据此定案四项，详见 [02-roadmap.md](02-roadmap.md) §4bis：

| 项 | 定案 | 连带 |
|---|---|---|
| V27-03 Hudson/PBS 分母 | 跟 pixy：`n<2` 的位点 num/den 同时剔除 | 改 D15-11 的测试；`compare_pbs.py` 需确保夹具不含 `n<2`（现夹具 `n=8`，已满足） |
| V27-27 `no_snps` | 跟 pixy：只数变异位点 | 与 V27-03 同批 |
| V27-12/V27-13 `\|S\|` 与 `no_sites` | 跟 pixy：只数「该群体有调用」的位点；**dxy 用成对交集** | 改 `04-algorithms.md:51` 的定义；改 `trailing_skipped_site_advances_the_windows` 的分母断言 |
| V27-11 genome `bp` | 新增只供分母用的 `denom_bp`，输出列不变 | 无下游破坏 |
| V27-39 多等位 | **项目不做多等位**：popg 接入 `MultiallelicPolicy::Reject`，多等位在过滤期剔除；`-M > 2` 报 usage 错 | 删除整套「skipped present」机制；D7-11 作废；`trailing_skipped_site_advances_the_windows` 需改构造方式 |

**用户第二条裁决（2026-11-08）**：本项目**不考虑多等位**。这使 V27-12 从「两套 present 定义并存」简化成单一计数器（见 §4bis.6）。

**关键事实**：`n=1` 位点上 scikit-allel 与 pixy 本身互相矛盾
（同一夹具：allel 0.22222222222222224、pixy 0.33333333333333337、reed 现状 = allel）。
所以 reed 不是「算错」，而是跟错了参照；按 pixy 改会同时改变 Hudson Fst 与 PBS 的数值。

V27-05 已于第三批按 V26-08 的「逐文件原子」落地：请求决定输出集合，空结果写表头/零谱/0/NA 行，不做整组事务。

### 已落地（第一至十批，2026-11-09）

| ID | 问题 | 验收 |
|---|---|---|
| V27-01 | 裸路由 `-o` 覆盖输入 | 回归测试；`-o <输入>` → exit 2、输入 md5 不变 |
| V27-02 | 失败运行发布截断 `.ld` | 回归测试（修复前失败）；旧 `.ld` md5/行数不变、无临时残留 |
| V27-08 | SFS 组名路径穿越 | 回归测试（反向验证过）；popg 与 `convert --sfs` 均不再逃出 |
| V27-09 | SFS 一维/联合谱撞名 | 回归测试（修复前失败）；改为明确报错 |
| V27-14 | `--sites-file` 与 `-t` 成并集 | 回归测试（popg 窗口/裸路由 + `stats --freq`）；单族行为不变 |
| V27-15 | `-r` + `-t` 成并集 | 回归测试（反向验证过）；真实数据 reed popg == reed filter == `bcftools view`（1416/492/0） |
| V27-05 | 空结果保留旧输出 | 按请求写表头/零谱/0/NA，旧文件被替换；失败运行仍不发布 |
| V27-06 | `step > size` 误报乱序 | 真实跨块夹具 + seam 单测，乱序仍拒绝 |
| V27-07 | 全基因组统计强制 contig 长度 | SFS/f2/LD 无长度均可运行 |
| V27-20 | `.ld` 不在成功文件清单 | sink 回报实际发布路径，CLI 列出 |
| V27-04 | LD 单元素尾块漏配对 | 4095/4096/4097/4098/8193 边界；4097 位点输出 8,390,656 对，`.ld` == decay 总数 |
| V27-03 | Hudson/PBS 分母混入 n=1 位点 | 跟 pixy：两组 n>=2 才进 num/den；1/3 实测与 pixy 一致 |
| V27-11 | genome `bp` 随 step 变 | `denom_bp` 独立，contig/genome 声明长度、BED 并集；step 不变 |
| V27-12 | variant-only `\|S\|` 定义 | 有调用才占 `|S|`；全缺失审计保留，与 pixy 一致 |
| V27-13 | π/θ 的 `no_sites` 矛盾 | π 逐群体、dxy 成对交集；pixy 逐列一致 |
| V27-27 | `no_snps` 是记录数 | 只数变异位点，跟 pixy |
| V27-39 | 多等位默认门与可见性 | 过滤期剔除；`-M > 2` exit 2；stderr/log 报计数；默认 == `-M 2` |
| V27-16 | LD 线程平方增长 | 单一执行池 + 总预算拆分，峰值 ≤ `--threads`（含 56） |
| V27-17 | LD 出错后不取消不 join | drop 触发 abort、pool join，库级取消测试通过 |
| V27-37 | `--threads 1` 多一线程 | 无辅助池串行，实测 1 线程 |
| V27-18 | 窗口累加器内存 | 流式输出会话 + 标量快照 + raw rows；100k 1bp 窗口约 25.0 MiB，输出不变 |
| V27-19 | BED 内存随 chunk 放大 | 按 contig 共享累加器；50k 区间 64 样本 49.0 MiB，跨线程逐字节一致 |
| V27-10 | variant-only 空窗口写 NA | 空窗为 0，Require/Bypass 仍 NA |
| V27-21…34 | CLI/网格/日志/性能项（除 V27-35） | levels、群体、BED、长度、退出码、日志、升序、临时告警、稀疏窗口与索引均通过回归 |
| V27-35 | 固定坐标 jackknife 网格与等权 SE | ADMIXTOOLS 观察位点动态网格 + 有效位点加权；共享网格与所有切分等价；真实 D 自然估计/SE 全精度一致 |
| V27-38 | 裸路由忽略 variant-only、pair/bp/θ 列量纲错误 | pair 与窗口 avg 对齐；bp 为 `Σsite_mean/L`；θ total=raw、bp=raw/L；缺失夹具与库级回归通过 |
| V27-36 | θ 直方图的零尾存储形状被当成不同统计请求 | 同 group/n_chrom 校验后补零尾；正反合并与整体累加一致；异请求/compact/非法内部形状仍拒绝 |

验收门（第十批实施后）：`cargo test --workspace` 68 个测试二进制 **1126 项全绿**、
`clippy --workspace --all-targets` 0 告警、`check-structure.sh` 16/16、`check-doc-drift.sh` 98 ok / 0 DRIFT（根 README 的线程上限已补齐）、
gold `run_easysfs.sh` PASS（`convert == popg` 逐字节一致）、`run_plink2_ld.sh` PASS、`run_vcftools_fst.sh` PASS、`run_pixy_bed.sh` PASS、`run_allel_fstats.sh` PASS（point estimate）；`run_admixtools_jackknife.sh` 以全精度参考诊断验证 D 自然估计与加权 SE PASS。另有 pixy 低覆盖/全缺失夹具逐列对拍。正常组名输出与修复前逐字节相同；
区域谓词另有真实数据三方对拍 reed popg == reed filter == `bcftools view`（1416 / 492 / 0）。

当前状态：**39 条全部实施完成**。V27-36 是库 API 防御性修复，未声称正常 CLI 可触发 panic；V27-35 改变 SE/Z 与可能的块数，point estimate 保留 reed 自然和/比值和；V27-38 修正裸宽表的量纲，θ total 从旧的平均量改为 raw。

## 测量环境

机器 sri-test-02（Intel Xeon Gold 5520+，56 可用逻辑 CPU），Linux x86_64，rustc 1.99.1。
构建 `cargo build --offline --release -p reed-cli`（`target-dir=/tmp/reed-target`，仓库 `target/` 产物早于源码）。
耗时/内存类数字均为**诊断值**，重复 3–5 次取中位数，不是优化前后对照；线程峰值每 2 ms 采样（可能漏掉更短峰值）。

> 说明：仓库 `target/release/reed` 的 mtime 为 2026-11-03，早于 `crates/` 中 2026-11-06 的源码改动。
> 所有复现均使用按当前源码重新构建的 `/tmp/reed-target/release/reed`。