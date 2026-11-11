# 02 拍板清单与拍板后的实施顺序

> 状态：**8 项全部拍板；Batch 0/1/2 已实施（见 STATUS.md）。**
> 拍板记录（2026-11-11）：D29-01 选 B（默认参考坐标对齐，重开 D15-26）、
> D29-02 选 A（保留 HapMap half-call）、D29-03 选 A（SFS 只算参与样本）、
> D29-04 选 A（输出冲突预检 exit 2）、D29-05 选 A（请求即写零谱）、
> D29-06 选 A（多倍体显式警告 + 计数）、D29-07 选 A（malformed 记录级一次）、
> D29-08 选 A（库内 validate；`to_args` target-only）。

## 拍板清单

### D29-01 MSA 长 REF 的坐标契约（对应 W29-01）——已拍板 B

现状：REF 300 bp 删除时输出宽度 1299（REF block 300 列 + 其后每参考位点 1 列），
由 D15-26 与 `e2e_convert.rs:1225-1227` 锁定。

**拍板：B —— 默认改为参考坐标对齐，重开 D15-26。**
这不是“修一个 Bug”，而是显式契约升级：`D15-26` 的 1299 契约被 v29 取代，现有测试与文档必须同步改写。

| 选项 | 内容 | 结论 |
|---|---|---|
| A | 维持 1299，只文档说明非参考坐标 | 不采用 |
| **B（已拍板）** | 默认参考坐标对齐；REF span 消耗对应参考列 | 采用 |
| C | 默认 1299，另加 `--msa-ref-anchored` | 不采用；默认必须是最科学语义，不把正确性藏在可选开关后 |

**B 的精确规则（实现前必须按此写测试）：**

1. 一条记录 `POS0`、`REF.len() = L` 占用参考区间 `[POS0, POS0 + L)`；下一次参考游标为
   `POS0 + L`，不是 `POS0 + 1`。
2. 位点块宽 `W = max(len(REF), len(ALT_1), ..., len(ALT_n))`，且 `W >= L`（因为 REF 本身
   是候选之一）。每个样本的等位值右侧补 `-` 到 `W`；全缺失仍写 REF（既有 vcf2msa 默认）。
3. 插入（ALT 长于 REF）：仍只消耗 `L` 个参考列，`W - L` 列是其他样本的 gap；
   不得把插入长度加到参考游标上。
4. 重叠记录：下一条 `POS0 < cursor` 时，只有**相同 POS 且相同 REF span**才允许
   last-wins；其余重叠（前一条 REF span 覆盖到后一条 POS）直接报 unsupported overlap，
   不得静默拼接。
5. REF/FASTA：按 `[POS0, POS0 + L)` 做大小写不敏感比较；越界或碱基不一致报错。
   重复 FASTA contig 名在读取参考时拒绝。
6. 输出宽度不再等于“参考长度 + indel 额外列”的猜测值：无插入时等于区域长度；
   有插入时等于区域长度 + 各插入带来的额外列。文档写清列坐标与参考坐标的映射。

**验收示例：**

- 10 bp 参考、POS 2 `REF=AA ALT=A`、POS 6 `REF=A ALT=G`：s0 期望
  `AA-AAGAAAA`（10 列，POS 6 的 G 在第 6 列），不再接受 `AAAAAGAAAA`。
- 300 bp REF 删除 + POS 601 SNP：两样本宽度均为 **1000**；s0 为
  `A` + 299 个 `-` + 其后参考/变异序列；s1 为 300 bp REF + 其后序列。
- 插入：`REF=A ALT=AT` 只消耗 1 个参考列，输出多 1 个插入列；后续位点参考游标不提前。
- 重叠：POS 2 `REF=AAA` 后接 POS 3 任意记录 → 明确报 unsupported overlap；
  POS 2 完全重复且 REF 相同 → last-wins。
- REF/FASTA 不一致、REF 越过 contig 末端、重复 FASTA contig → 明确错误，不写输出。

**文档与对拍影响：**

- `docs/spec/cli.md` 的 `--msa` 语义、`docsite/src/content/docs/convert.md` 与中文版同步改写；
  D15-26 在历史轮次中保留，但标注“v29 已取代”。
- `vcf2msa` 对拍仅保留在参考坐标语义可对齐的简单用例；长 REF/插入的 golden 必须按
  v29 新契约独立手算，不能继续用旧的 1299 结果当参照。
- 不再声称“MSA 无外部 aligner 仍与任意 vcf2msa 输出逐字节一致”；只承诺在
  SNP/简单 indel 域内一致。

### D29-02 HapMap 双字母半缺失（对应 V29-04）——已拍板 A

**拍板：A —— 保留 half-call。** 1 缺失 + 1 已调用写 `i/.` / `./i`，两缺失写 `./.`；
禁止 `i/i` 静默纯合化。

依据“看 bcftools”：仓库 bcftools 的 `view -Ov` 读入 `0/.`、`./1` 时原样保留文本，
不纯合化、不改成 `./.`；`query -f '[%GT]\n'` 同样逐字返回。下游 filter/stats 按既有
half-call 策略解释，不由 HapMap 转换层替它决定。

**拍板后动作**：补 `AA/AG/AN/NA/NG/GN/NN` 与 `-N`/`N+` 表；TASSEL 双向对照进 gold。
两缺失写 `./.`；单缺失保留缺失 slot，未相位次序固定为“已调用在前、缺失在后”。

### D29-03 SFS 的 `-s` 子集语义（对应 V29-02）——已拍板 A

**拍板：A —— 只有被保留且属于该 SFS 群体的样本参与 GT 校验与解码；被排除样本按 missing。**
`--samples` 必须等价于物理预裁剪文件；不能因为被排除样本的 `0/9`、多倍体或超长 token
让保留样本的位点消失。

**拍板后动作**：同时改 `sfs_accepts_gt` 与 stats scanner 入口；验收是
“SFS `-s` 输出 == 物理裁剪后的 SFS 输出”，覆盖被排除样本的越界、多倍体、长前导零 GT。

### D29-04 输出路径冲突与覆盖策略（对应 V29-03、V29-07）——已拍板 A

**拍板：A —— 运行前计算完整输出集合，任何重复路径或与 input/reference/population/log
的别名冲突都 exit 2，不写任何文件。**

包括：`foo.012` 被 `--012 -o foo` 覆盖；MSA 重复 region / flank 收敛到同一路径；
库入口 `ConvertTarget::Vcf{output: Some(input)}`；MSA contig 名含 `/` 或与参考/输入撞名。
保持逐文件原子，不做多文件事务。

**拍板后动作**：`ConvertRequest::validate` + `planned_outputs()`；MSA contig 做文件名
组件校验；CLI 与库入口共用同一计划。

### D29-05 空 SFS 与 header-only 输出（对应 V29-05）——已拍板 A

**拍板：A —— 只要请求了 `--sfs` 就写零谱。**
与 v27-05“请求决定输出集合”一致；header-only、全记录被跳过、空文件三种输入形态
输出集合一致。

**拍板后动作**：请求存在即初始化零 spectrum / 零 report；plain/BGZF/BCF 三输入形态
统一测试，复用同一个空报告构造函数。

### D29-06 多倍体在 SFS / `--012` 的策略（对应 W29-03）——已拍板 A

**拍板：A —— 暂不支持多倍体；进入 SFS/012 时显式警告一次 + 计数，不静默写零谱/`-1`。**
文档写明支持域：diploid + haploid-as-half-call。完整 allele-level 倍性支持另开新能力轮。

**拍板后动作**：SFS 跳过与 `--012` 的 `n != 2` 分支加诊断/计数；测试断言不再出现
“全零谱且无提示”。

### D29-07 malformed 计数语义（对应 V29-09）——已拍板 A

**拍板：A —— 记录级 `malformed` 每行最多 1 次；per-target 明细如需保留则另加字段。**
`malformed` 不得再出现大于 `total_records` 的值，也不随启用 target 数变化。

**拍板后动作**：先加“同一坏记录单/多 target 计数不变量”测试，再改实现；日志文案与
`ConvertReport` 文档同步。

### D29-08 库校验边界与 `to_args` 定位（对应 V29-10）——已拍板 A

**拍板：A —— 库内 `ConvertRequest::validate` 覆盖 header 无关不变量，`ConvertPlan`
补 header 相关部分；`to_args` 明确为 target-only，不承诺完整重放。**

`validate` 至少覆盖：空/重复 target、target 互斥、VCF target 与样本子集的组合、
压缩 level 1–9、outgroup 属于保留样本、MSA reference/region/collapse、SFS projection。
`to_args` 不再被当作执行/校验的序列化格式；日志另记实际样本与群体参数。

**拍板后动作**：先补库直调测试（空 targets、重复 target、`level=0/10/255`、
bin-only outgroup、VCF target + keep mask），再让 CLI 复用同一 validate 结果。

## 拍板后的实施顺序

### Batch 0 — D29-01 契约设计与失败测试

1. 写 v29 MSA 规则小节（§D29-01 的 6 条），把 D15-26 标为被取代。
2. 先写失败测试：10 bp `AA-AAGAAAA`、300 bp 删除 1000、插入、重叠拒绝、
   REF/FASTA 不一致、重复 contig。
3. 更新 `e2e_convert` 的 `msa_long_deletion_is_not_truncated` 与相关注释。

### Batch 1 — 已证实正确性

4. **V29-01**：GT 子字段真实 span；`DP:GT` 三目标 parity。
5. **V29-02 + D29-03 A**：SFS 只校验/解码参与样本；`-s` == 预裁剪。
6. **V29-03 + D29-04 A**：统一输出计划与别名检查；重复 region 提前报错。
7. **V29-04 + D29-02 A**：HapMap 半缺失写 `i/.` / `./i`，表驱动测试。

### Batch 2 — 契约补齐

8. **D29-05 A**：header-only 与全跳过统一 SFS 输出集合。
9. **V29-06**：`merge_sfs` 计数；多 chunk 夹具。
10. **V29-08 + D29-01 B**：MSA REF/FASTA 一致性、重复 contig、重叠 span 规则。
11. **V29-09 + D29-07 A**：malformed 记录级计数。
12. **V29-10 + D29-08 A**：library validate/plan；`to_args` 定位与文档。

### Batch 3 — 下沉（只在 Batch 1/2 需要共用能力时落地）

13. **V29-11**：GT 词法与 span 原语先服务 V29-01，再迁移 PCA/scanner。
14. **V29-12**：`ConvertPlan` / `planned_outputs` 先服务 V29-03，再删 CLI 重复策略。
15. **V29-13**：`StatsJob::run_lines` + 统一空报告；删除 sanitized body 复制。
16. **V29-14**：HapMap chunk/gzip 复用；spill 注册清单合一。

### Batch 4 — 性能/资源（先基准，后改）

- spill 文件身份（fd/`O_NOFOLLOW`）、MSA region 索引 seek、组合 SFS 的 native record view、
  `--memory` 的真实 RSS 口径与临时磁盘预算。
- 这些不因“看起来慢”直接改；先建立 convert 基准/资源回归入口。

## 验收门

- `cargo test -p reed -p reed-cli -p reed-io -p reed-stats`
- 新增回归均需**反向验证**（未修复时失败）：`DP:GT` binary、SFS `-s` 等价、
  MSA 重复 region、HapMap 半缺失、header-only SFS、多 chunk `n_records`、
  MSA 参考坐标对齐（10 bp、300 bp、插入、重叠、REF 不一致）。
- `cargo clippy -p reed -p reed-cli -p reed-io -p reed-stats --all-targets -- -D warnings`
- `cargo fmt --all -- --check`、`git diff --check`、`ci/check-doc-drift.sh`。
- gold：`gold/popgen/run_convert_parity.sh` 记录其二进制来源；长 REF/插入的 MSA
  用手算 golden，不再用 1299 或 vcf2msa 的旧结果当参照。
- 普通双等位、无 `DP:GT`、无重复 region、无半缺失输入的非 MSA 输出逐字节不变；
  MSA 输出按 v29 新契约整体变化，需在发布说明标注。

## 不做

- 不以 Bug 名义静默改 D15-26；D29-01 B 是显式契约升级，必须同步改文档与测试。
- 不把 convert / PCA / stats 的 GT 策略统一成一个 `GtCode` 路径。
- 不新增多文件事务层；不重设计 `to_args` 为全量序列化（D29-08 已选 A）。
- 不在本轮新增 PLINK/EIGENSTRAT 格式或 MSA 外部 aligner。
- 不把“测试全绿”当作 convert 边界已完整的证明。
