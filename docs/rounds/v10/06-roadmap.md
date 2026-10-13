# 06 修复路线（R0–R3）

> 原则：**先复现（红测试先行，把"静默损坏/panic/必失败"钉成可跑的门），
> 再修单行级缺陷（一次 PR 一个缺陷族），最后做结构对齐（D10-5）
> 与兼容性收尾**。修复对象是 01-review.md §1 的缺陷清单 D10-1…D10-8；
> 不重开 v6 的 Y 系列、v7 的 Z 系列、v8 的 P 系列、v9 的 Q 系列结构项。
> 规范依据见 02-architecture.md 的对应小节（N1…N8）。每一步都必须
> 满足 §1 通用门。

---

## 1. 通用验收门（每一步都必须满足）【规】

1. **双 profile 全绿**：`cargo test --workspace`（debug）与
   `cargo test --release --workspace` 都通过；
   `cargo clippy --release --workspace --all-targets` 不新增警告。
2. **红 → 绿**：每个 P0/P1 缺陷先提交复现它的失败测试
   （D10-1：放回原阈值 32766 编码断言失败；D10-2/4：`catch_unwind`
   捕获 panic 即红；D10-5：多 chunk 夹具 `run_index` 报错即红），
   再提交修复转绿。测试名描述缺陷本身。
3. **htslib 真值门**：凡输出 BCF 的步骤（R1 的 D10-1、R3 的
   D10-7），用工作区 `bcftools view` 读回写出物，退出码 0 且
   记录数一致。
4. **确定性**：threads ∈ {1, 4} 输出 sha256 相同。
5. **性能**：不回退超过 10%；`bench/suites/e2e/results.tsv`
   追加一行。
6. **同一 PR 内删除被替代的旧实现**，不留并行代码
   （v2 06 §1 纪律 7 沿用）。
7. 本目录相应条目标记完成并注明 commit；
   `docs/rounds/v10/STATUS.md` 同步更新；
   落地一步即在 `ci/check-doc-drift.sh` 登记 01-review §3 的
   对应守卫。

---

## 2. 步骤

### R0 红测试先行（D10-1/2/4/5 的复现夹具）— 门的基础

> R0 把四个"会坏"的情形钉成测试，不修任何代码。四个测试在 R0
> 提交里全红（或 catch 到 panic），后续每步只让自己那一条转绿，
> 防止"修一个顺手破坏另一个"无处可见。

| 项 | 内容 | 门 |
|---|---|---|
| R0-a | D10-1 边界编码测试：等位 62/63、16382/16383、32766/32767 六组，断言 GT 描述符类型与写出字节。**修复前：16383…32766 组断言 INT16（错）** | 测试在 `native_encode.rs` 的 `tests` 模块；放回原阈值 32766 时 `cargo test -p reed-io` 红 |
| R0-b | D10-2 畸形描述符测试：手工拼 `(NULL, count=1)` 的 INFO 字段与 FORMAT 字段（`render.rs` 路径），`catch_unwind` 断言"返回 Err 且无 panic"。**修复前：catch 到 index-out-of-bounds panic** | `cargo test -p reed-io`（debug 才看得见 panic；release 下是 UB 级越界，因此必须在 debug 跑） |
| R0-c | D10-4 畸形块测试：`BSIZE ∈ {0, 3, 17}` 的字节流 × 三条入口（流式 scan、流式 chunk、mmap），断言 `Err` 且无 panic | 构造方式仿 `chunk/bgzf.rs` 的 `compress_lines`，只改 BSIZE 字节 |
| R0-d | D10-5 多 chunk 索引夹具：64 字节小块 + 40 字节 chunk 目标（仿 `chunks_cover_all_records_with_many_seams`），`run_index` 断言成功且与单 chunk 索引字节一致。**修复前：`record line has fewer than 9 columns`** | `cargo test -p reed --test e2e_index` 新增 |

### R1 单行级修复（D10-1 / D10-2 / D10-3 / D10-4）— 缺陷修复的主量

> 四个缺陷的修法都是"一处判断/一个常量"（02 §1–§4），合并为一个
> PR 或四个小 PR 均可；每个缺陷的红绿对照必须在自己的 commit 里
> 可验证（先 revert 修复即红）。

| 项 | 内容 | 规范 | 门 |
|---|---|---|---|
| R1-a | GT 选宽对齐 N1：两处阈值 62/16382（或合并为按编码值扫 min/max 的 helper），`push_gt_code` 上界收紧；R0-a 转绿 + `bcftools view` 读回边界 BCF | 02 §1 | 通用门全部；守卫 #1 |
| R1-b | `typed()` 拒绝 `(NULL, >0)`（N2）；R0-b 转绿 | 02 §2 | 守卫 #2 |
| R1-c | `BcfDict::from_header` IDX 限界（N3）：超界丢弃该 IDX 并 warn、回退自动编号；`IDX=10^9` 测试 + `honors_explicit_idx` 回归 | 02 §3 | 稀疏化属【决】，本轮只落限界 |
| R1-d | `parse_block_size` 最小块长（N4）；R0-c 转绿 | 02 §4 | 守卫 #3 |

### R2 结构对齐：文本 index 残片承载（D10-5）

> 本步是本轮唯一的结构改动，动 `IndexPartial` 的形状与归约器裁决
> 逻辑，但**有现成模板**（同文件 `BcfIndexJob`）。改动面：
> `index.rs` 的 `IndexJob::run_chunk`/`run_seam`、`IndexPartial`、
> `IndexSink::accept`。引擎（`pipeline.rs`）与 v9 W2 出口不动。

| 项 | 内容 | 规范 | 门 |
|---|---|---|---|
| R2-a | `IndexPartial.head/tail` 改承载体（voff + fragment），worker 端用 `VirtualMap` 预映射 voff | 02 §5.1 | R0-d 仍红（形状未接） |
| R2-b | 归约器裁决：`!prev_continues` 解析 head；`prev_continues` 丢弃；tail 暂存为 pending fragment | 02 §5.2 | — |
| R2-c | `run_seam` 用暂存 fragment 的 voff 解析引擎拼接的完整行（替代忽略 `_line`）；畸形 seam 行经 v9 出口报错 | 02 §5.2-3 | R0-d 转绿；seam 畸形夹具报错不挂死 |
| R2-d | 删除 worker 端残片解析旧路径（急切 `parse_record` on 残片） | 02 §5、通用门 6 | 守卫 #6（`run_seam` 使用 `_line`） |

### R3 兼容性收尾（D10-6 / D10-7 / D10-8）+ 文档收口

| 项 | 内容 | 规范 | 门 |
|---|---|---|---|
| R3-a | `plain.rs` 换 `MultiGzDecoder`（N6）；串联 gzip e2e（filter 记录数 = a + b） | 02 §6 | 守卫 #4 |
| R3-b | 输出头注入 PASS（N7）：`filter/mod.rs` 输出头组装处，PASS 位于 FILTER 定义之首；无声明输入 `-Ob` + `FILTER=PASS` 记录 e2e；`bcftools view` 读回 | 02 §7 | 守卫 #5；真值门 |
| R3-c | 四处计数乘法 `checked_mul` 化（N8，【荐】） | 02 §8 | 现有测试全绿 |
| R3-d | spec 回改（02 §9 清单）+ drift 守卫登记 + STATUS 收口 | — | `bash ci/check-doc-drift.sh` 全绿 |

---

## 3. 顺序与依赖

```
R0（红测试）→ R1（单行级四连）→ R2（残片承载）→ R3（收尾）
```

- R1 内部四缺陷互不依赖，可并行分支；
- R2 依赖 v9 的 W2 出口（已落地），不依赖 R1；
- R3 各项互不依赖。**若排期紧张，R3-a/b（P2）可单独成轮**——
  但 D10-6 的静默截断与 D10-7 的报错都面向真实用户输入，不建议
  晚于下一轮。

## 4. 本轮不做什么

| 项 | 原因 |
|---|---|
| `BcfIndexJob` 的任何改动 | 它已是 N5 的模板（01-review §1 D10-5） |
| 读侧 PASS 注入（`BcfDict`） | 【荐】项（02 §7.4），触发面是理论输入，列入后续 |
| 字典稀疏化（HashMap 存储） | 【决】项（02 §3.2），N3 限界已消除 OOM，稀疏化是优化不是修复 |
| 变异类型分类（`bcf_set_variant_type` 对照） | 用户复审已在此前轮次给过结论，本轮不复开 |
| `bcf_record_check` 全量移植 | htslib 的完整记录校验比 N2 宽得多；本轮只补"会 panic/会损坏"的最小集，全量校验若需要应单独立项 |
