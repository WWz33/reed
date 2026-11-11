# doc_v29：convert 模块复审与拍板结论

> 立项 2026-11-11。承接 v28 的 `reed-io` 复审，本轮只审 `convert` 及其直接依赖的
> `reed-io` / `reed-stats` 调用面；重点是“哪些是已证实的 Bug、哪些是契约空白、
> 哪些上轮判断需要撤回”。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 逐条证据、状态标记、撤回/降级项、下沉边界 |
| [02-decisions.md](02-decisions.md) | D29 拍板结论、选项与拍板后的实施顺序 |
| [STATUS.md](STATUS.md) | 当前状态：8 项已拍板、Batch 0/1/2 已实施；验证与新增测试锚点 |

## 总结

共 **17 条**：10 条发现（V29-01…10）、4 条下沉设计（V29-11…14）、3 条上轮主张的撤回或降级（W29-01…03）。

**2026-11-11 状态：8 项全部拍板；Batch 0/1/2 已实施并有回归。** 逐条见
[01-review.md](01-review.md)，验证与测试锚点见 [STATUS.md](STATUS.md)。
拍板：D29-01 选 B（默认参考坐标对齐，重开 D15-26）、D29-02 保留 HapMap half-call、
D29-03 SFS 只算参与样本、D29-04 输出冲突预检 exit 2、D29-05 请求即写零谱、
D29-06 多倍体显式警告+计数、D29-07 malformed 记录级一次、D29-08 库内 validate。

### 已实施（10 条发现 + D29-01/D29-06）

- **V29-01**：`--nexus-binary` 在 `FORMAT=DP:GT` 时从样本列起点切 GT，合法 `0/0`、`1/1`
  输出成 `? ?`；现改用 GT 子字段真实 span。
- **V29-02**：`--sfs` 的 `-s` 子集被已排除样本污染；现只解码/校验参与样本，
  `-s` 等价于物理预裁剪（text scanner 与 BCF native 同语义）。
- **V29-03**：输出路径规划不完整——`foo.012` 可被 `--012 -o foo` 覆盖；现运行前
  计算完整输出集合，拒绝重复路径与 input/reference/populations 别名。
- **V29-04**：HapMap 双字母半缺失 `AN`/`NG` 被写成 `0/0`/`1/1`；现保留缺失 slot。
- **V29-05**：header-only VCF + `--sfs` 不写任何 SFS 文件；现请求即写零谱。
- **V29-06**：SFS 的 `n_records` 在多 chunk 合并时只保留第一个 partial；现已合并。
- **V29-07**：MSA contig 名不安全字节（`/` 等）现做百分号编码，避免子目录/撞名。
- **V29-08**：MSA 现校验 REF 与 FASTA（大小写不敏感）、越界、重复 FASTA contig、
  重叠 REF span。
- **V29-09**：同一条 malformed 记录不再按 target 数重复计数。
- **V29-10**：库内 `ConvertRequest::validate` 覆盖 header 无关不变量；`to_args`
  明确为 target-only。
- **D29-01 B**：MSA 改为参考坐标对齐；300 bp 删除夹具由 1299 变 1000。
- **D29-06**：多倍体进入 SFS/`--012` 时显式警告、计数；`--sfs` 跳过该记录，
  `--012` 保持 vcftools 的 `-1`，不再静默。

### 撤回/降级（3 条）

- **W29-01**：MSA 长 REF 输出 1299 是 v15 D15-26 的既定验收契约，**不是 Bug**；
  上轮“G 错位到第 7 位”的说法也不成立。D29-01 已拍板 B：默认改为参考坐标对齐，
  D15-26 被显式取代，属于契约升级而非 Bug 修复。
- **W29-02**：folded 一维 SFS 的 `k=0/n` 端点合并是既定设计与测试锚点，不是 Bug。
- **W29-03**：多倍体在 `--sfs` 被静默跳过、`--012` 写 `-1` 是策略空白；
  D29-06 已拍板 A（显式警告 + 计数，暂不支持）。

### 下沉设计（4 条）

- **V29-11**（未实施）：共享无损 GT 词法与真实 span；V29-01 已直接用现有
  `subfield_span` 修复，词法统一留待后续。
- **V29-12**（部分落地）：`ConvertRequest::validate` 与输出路径计划原语已落地；
  完整 `ConvertPlan` 结构、CLI/target 解析去重未抽出。
- **V29-13**（未实施）：`StatsJob::run_lines` + 统一空 SFS 报告出口；去掉 convert 的
  sanitized 文本重建。
- **V29-14**（未实施）：HapMap 复用 reed-io 的 plain chunk / gzip sniff；spill 注册
  清单收敛为一份。v28 已落地的 `AtomicWriter`、`LineRef` 访问器不重复立项。

## 边界

- 8 项均已拍板；Batch 0/1/2 已实施。V29-11/13/14 与 `ConvertPlan` 的完整抽出是
  后续下沉工作；Batch 4 性能项需先建基准。
- D29-01 B 是显式契约升级：必须同步改 D15-26、测试、CLI 规范与 docsite，
  不能伪装成“修 Bug”。
- 逐文件原子，不新增多文件事务；不把输出计划扩成全局事务层。
- 不统一 convert / PCA / stats 的 GT **策略**，只考虑共享无损词法。
- 不把各格式的筛选门槛、allele 规范化、缺失策略上层抽象为统一框架。
- 不在本轮新增 PLINK、EIGENSTRAT 等新格式；`to_args` 不做完整序列化重设计。
