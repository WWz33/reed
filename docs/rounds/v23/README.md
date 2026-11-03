# reed doc_v23：全库审查复核（2026-11-03）

> 触发：仓库外审查报告 `../../../../reed-rs-review-2026-11-03.md`（9 个子代理、
> 约 8.9 万行、只读审查）。本轮的职责不是直接修，而是**逐条复核**：
> 哪些是真 bug、哪些是误报、有没有报告没覆盖到的同类位置、修的话代价多大。
>
> 复核方式：当前 HEAD `8b54588`（doc_v22 R2 已提交）；本机 debug
> 构建（`/tmp/reed-target/debug/reed`）、`bcftools 1.24-2`/htslib 1.24-1。
> 崩溃类用最小 fixture 实跑；数值类用独立探针 crate（`/tmp/v23probe`，
> path 依赖 reed-core/reed-reduce，不改仓库）实算；纯内部不变式用代码复核。
>
> 复核阶段**不改源码**。拍板后的实施在同轮完成：A 组 6 条 + V23-08/09 +
> 拍板 5 项已落地，明细见 [02-implementation.md](02-implementation.md)。
>
> 本目录：[01-verification.md](01-verification.md)（逐条证据与判定）、
> [02-implementation.md](02-implementation.md)（实施、测试与验收）、
> [STATUS.md](STATUS.md)（全部条目的判定/建议/状态）。

## 结论速览

审查报告的 **8 条 Must-fix：7 条属实，1 条不成立**。

| # | 审查条目 | 复核判定 | 证据强度 |
|---|---|---|---|
| 1 | BCF `typed`/`typed_int` 无限递归栈溢出 | **属实，P0** | 758 字节损坏 BCF 实测 SIGABRT（exit 134） |
| 2 | BCF CHAR 类型 GT 让渲染/列解码 panic | **属实，P0**（范围需修正） | `convert`/`stats` 实测 panic；`filter` 有守卫、报 malformed |
| 3 | VCF 样本区尾部空列 + 编辑 panic | **属实，P0**（触发条件需修正） | 两个尾随 tab + `-s`/`--min-gq` 实测 panic |
| 4 | AF/MAF 在 AN=0 时 NaN，`!=` 判 true | **不成立** | bcftools 同输入同结果；改 Missing 会反而偏离 bcftools |
| 5 | 表达式嵌套无深度上限，栈溢出 | **属实，P0/P1** | 10,000 层括号实测 SIGABRT（exit 134） |
| 6 | reed-reduce 三处 NaN 排序 panic | **属实，P1**（仅库 API） | 独立探针实测 `isomap` panic；CLI 无调用方 |
| 7 | PCA 特征向量无符号约定 | **属实，P1；已拍板 A**（`top_k` 统一符号规范化，输出符号会变） | 同一矩阵 fallback/faer 实测 PC1/PC2 整体变号 |
| 8 | `round_half_even` 容差把近平局当平局 | **属实，P0/P1** | 探针实测 2/2/32，R 为 1/3/31 |

### 复核中修正/加强的结论

- **①1 触发条件**：单个尾随 tab 不触发；需要多一个**空样本列**（行尾连续
  两个 tab）且计划编辑样本区。带掩码（`--min-gq`）即使没有 `-s` 也会 panic。
- **①2 严重度上调**：负 INFO Float 在 BCF→VCF 渲染时丢负号不是"潜伏"——
  任意带负浮点 INFO 的 BCF（如 `InbreedingCoeff` 一类）经 `reed filter -Ov`
  就会被静默改写（`-1.5e-07` → `1.5e-07`）。`float_into` 的 FORMAT 浮点
  渲染同源。
- **②2 范围修正**：`filter` 的 GT 解码路径**已有**整数类型守卫（干净报
  `BCF GT field is not integer-typed`）；panic 出现在 `convert`/`stats` 的
  渲染路径（`int_element` 从 1 字节元素读 4 字节）。float 列解码路径
  （`float_element_column`）结构上同样危险，但本次 CLI 探针未触发（该输入
  未走原生列解码），需要函数级测试确认。
- **④1 驳回**：`AF!=0` 在 AN=0 位点上 bcftools 也保留（NaN != 0 为真）。
  reed 与 bcftools 逐条一致（`AF>0`/`AF<1`/`AF==0`/`AF!=0` 全对拍）。
  把 NaN 改为 `Val::Missing` 会让 `AF!=0` 只留真 AF>0 位点，反而破坏
  bcftools 兼容；建议只补文档，不改语义。
- **⑥2 范围修正**：`pos_log` 在 `finish_contig` 里 `clear()`，不是"全 contig
  只增不减"；增长上限是**单条 contig 的位点数**。仍是内存问题（大染色体
  × 每线程），但量级按单 contig 算。
- **⑤3 可达性修正**：CLI 拒绝 `--ld-pairs` 单独使用（要求同时给
  `--ld-decay`），"pairs-only 产出一行 decay"只在库 API
  （`LdRequest { decay_max_bp: None, pairs: true }`）可达；CLI 同时给两个
  标志时 `--ld-pairs` 的无上限表**仍然可达**。
- **⑧1 可达性修正**：`isomap`/`grm_isomap`/`knn_from_distances` 在 reed CLI
  没有调用方，只影响库 API；因此从 P0 降为 P1。
- **⑧2 实测**：不只是"可能不同"——同一对称矩阵，fallback backend 的 PC1
  首分量 `-0.99845624390341792`，faer backend 为 `+0.99845624390341725`；
  PC2 同步变号。Cargo 默认 `faer-backend`，切后端或升级 faer 会逐字节改
  输出。

## 实施结果（2026-11-03）

已落地：V23-01…06（P0）、V23-08/09、V23-14/21/31/47/50（拍板项）。
`ci/test-parallel.sh` **893 passed / 0 failed**、clippy 0 告警、doc-drift
97 ok；真实数据 DP/PSC 与负浮点 BCF 均与 bcftools 逐行/逐字节一致。
两处实施偏离拍板文本（V23-05 上限按实测 128 而非 500；V23-03 不改公开
签名），理由见 [02-implementation.md](02-implementation.md)。

## 拍板结果（2026-11-03）

| 编号 | 决定 |
|---|---|
| V23-09 | 按结果正确性改：`SymEigResult::top_k`（及 `eig::SymEig::top_k`）统一做"最大 \|分量\| 取正"；接受输出符号变化；需 V23-48 双 backend 对拍 |
| V23-14 | 选项 C：sink 按 contig 顺序直接追加最终 `PREFIX.ld`；首个 pair 懒创建并写表头（全无 pair 仍不建文件），无临时文件、输出不变、内存 O(单 contig) |
| V23-07 | **丢弃（不修）**：不接受多等位；该条从 v23 修复清单丢弃。默认门已丢弃，显式 `-M >2` 时 AC 按首 ALT（不支持路径） |
| V23-50 | 按推荐 A：保持 AF/MAF 的 NaN 语义（与 bcftools 一致），补文档/测试 |
| V23-21 | 按推荐 A：前缀命令统一拒绝 `-o -` 与空前缀；`filter` 的 stdout 语义不变 |
| V23-47 | 按推荐 A：oracle 保持 complete-case，补缺失 golden，文档记录与 allel 的差异 |

> 以上 5 项已实施（1 项丢弃）；其余条目仍待拍板、未改代码。

## 建议修复顺序

| 优先级 | 项 | 理由 |
|---|---|---|
| P0 | BCF `typed` 递归、CHAR 渲染 panic、尾部空列+编辑 panic、深度嵌套表达式、负 Float 符号、`round_half_even` | 崩溃或静默写错数据；均为小改动（几行到十几行） |
| P1 | 三处 NaN 排序、特征向量符号（已拍板 A）、pass0 多 contig 预算、`read_reference` 整读、`matrices::grm` 全量拷贝、`pos_log` 截断、tql2 不收敛、frame_chain 检查顺序、`standardize` 的 Inf | 内存炸弹/复现性/错误分解 |
| P2 | `--ld-pairs`（已拍板 C）、`-r/-R/-t` 组合、`-o -`/空前缀（已拍板 A）、缓存上限、热路径分配、死代码/精简、oracle 缺失数据 golden（已拍板 A） | 体验/性能/测试护栏，不影响当前正确结果 |
| 驳回/丢弃 | ④1 AF/MAF NaN → Missing（保持）、④2 AC 多等位向量语义（丢弃） | 与 bcftools 一致/不支持多等位 |

## 本轮不做

- 不改任何源码、不跑全量修复后的回归；不把条目写成"已修"。
- 不把"截断风险/性能/精简"条目当作 bug 处理；它们在
  [STATUS.md](STATUS.md) 里单独标注。
- `reed-rs-review-2026-11-03.md` 的原始行号以审查时 HEAD 为准；本轮复核的
  行号已按当前工作区更新（v22 修复已改过部分文件）。
