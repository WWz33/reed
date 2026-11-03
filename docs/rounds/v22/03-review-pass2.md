# 第二轮评审：字段解析、索引与 stats 消费（doc_v22，2026-11-03）

> 触发：第一轮评审（[02-review.md](02-review.md)）后又提出的
> `GQ` 缺失、`u8` 子字段索引、整数饱和哨兵、`inf/nan` 四项；以及对
> `DP`/`AD`/`GT`/浮点字段的同类问题排查。
>
> 复核方式：当前 HEAD `717de0e`、本机构建 `reed 0.4.0`（debug）、
> `bcftools 1.24-2`/htslib 1.24-1；最小 fixture 在 `/tmp/v22verify`。
>
> **修复已落地（2026-11-03）**：A–D 批全部完成，R2-8/E 批留待后续轮。
> 逐项修复、测试名与验收记录见 [STATUS.md](STATUS.md) 的 "V22-R2"
> 一节；本文保留原始复现与计划，作为修复前基线。

## 结论速览

| ID | 发现 | 判定 | 处理 |
|---|---|---|---|
| V22-R2-1 | `u8` 子字段索引：FORMAT 第 256 个及以后的 GT/DP/GQ/AD 静默变 None | **属实，病态** | 已修（索引改 `usize`） |
| V22-R2-2 | 整数饱和到 `u32::MAX` 与 `INT_MISSING` 重合；AD 求和同样会撞哨兵 | **属实，病态** | 已修（i32 边界 + u64 求和钳制） |
| V22-R2-3 | `parse_float` 接受 `inf`/`nan` | **部分不成立**：GQ 是 Integer，不走此解析器；Float `inf`/`nan` 与 htslib `strtod` 一致 | 不改，文档/测试固化 |
| V22-R2-4 | `--hist dp` 在 AD-only header 上被 CLI 拒绝 | **属实，非病态** | 已修（gate 接受 DP 或 AD） |
| V22-R2-5 | `parse_int` 与 `parse_int_token` 两套整数语义（`+`/`-`） | **属实** | 已修（单一解析器） |
| V22-R2-6 | stats 侧不校验 DP/AD/GQ 类型，Float 字段静默按缺失 | **属实** | 已修（`validate_report_tags`） |
| V22-R2-7 | 真实 Float `nan` 与缺失哨兵重合 | **属实，病态** | 只记录设计选项（§8），后续轮 |
| V22-R2-8 | 负数/有符号 FORMAT Integer 无法比较 | **属实**（表示层限制） | 后续轮（§9），当前负数按 missing |
| V22-R2-9 | 普通“GQ 列缺失”是否会跳过位点 | **不是 bug**（见 §1） | 不改，补测试 |

## 1. 先回答：GQ 列缺失会不会静默跳过全部位点

**普通缺失不会。** 三种形态实测：

| 场景 | `reed filter --min-gq 20` | `reed stats --hist gq` |
|---|---|---|
| header 未声明 `FORMAT/GQ` | exit 2：`FORMAT field GQ is not declared in the header` | exit 2：`--hist "gq": the header declares no FORMAT/GQ` |
| header 声明 GQ，部分记录 FORMAT 无 GQ | 3 条记录全部输出；无 GQ 的记录原样保留 | `#missing` 按缺失样本计数，位点照常扫描 |
| header 声明 GQ，所有记录都无 GQ | 记录全部保留 | 只有 `#missing = n_samples × n_records`，无 bin |

**真正会让 GQ 掩码“对所有位点失效”的是 V22-R2-1**：GQ 在 FORMAT 的第
256 个位置（0-based index 256）时，`index_of` 的 `u8::try_from().ok()`
把它变成 `None`；扫描器把每个样本都当作 GQ 缺失。`--min-gq` 默认
`on_missing = keep`，于是所有位点都不被掩码——**不是位点被丢弃，而是
掩码静默 no-op**。同理 GT 在 index ≥ 256 时，`gt_subfield_index` 也返回
`None`，即使掩码触发也不会写 `./.`。

## 2. V22-R2-1：`u8` 子字段索引

### 复现

构造 257 个 FORMAT key（0-based index 256 处的 key 溢出）：

```text
# fixture A: GT 在 0，GQ 在 256，GQ=10
FORMAT = GT:K0:...:K254:GQ
reed filter --min-gq 20 -Ov   # 输出 GT 仍是 0/1（应掩码为 ./.)
```

```text
# fixture B: GQ 在 0，GT 在 256，GQ=10
FORMAT = GQ:K0:...:K254:GT
reed filter --min-gq 20 -Ov   # 输出 GT 仍是 0/1
```

### htslib 边界

```text
255 个 FORMAT key: bcftools view 通过
256 个 FORMAT key: [E::vcf_parse_format_dict2] FORMAT column ... lists more
                   identifiers than htslib can handle
257 个: 同上
```

即：**bcftools 最多接受 255 个 key**；reed 接受更多但静默错读。属于
"病态输入"，但静默错读比明确报错更危险。

### 受影响位置

| 位置 | 内容 |
|---|---|
| `crates/reed-io/src/vcf/format_layout.rs:157` | `index_of` 把 GT/整数/浮点 key 的位置转 `u8` |
| `crates/reed/src/filter/record_pipeline.rs:736` | `gt_subfield_index` 转 `u8`（滤器 GT 编辑） |
| `crates/reed/src/filtered.rs:180` | 同一 `gt_subfield_index` 逻辑（filtered/view 路径） |

### 计划

1. 首选把 `Layout.gt_index`、`IntLayout.index/fallback_index`、
   `Layout.float_indices`、`gt_subfield_index` 统一改为 `usize`（或 `u16`），
   去掉人为上限；扫描器和行编辑都直接使用位置。
2. 若坚持对齐 htslib：在扫描/记录入口对 `format_text()` 的 key 数 > 255
   报 `MalformedRecord`，所有消费者（filter、stats、rename、convert）
   统一拒绝，而不是各自静默 `None`。
3. 测试：255/256/257 key 三档；GT/GQ/DP/AD 分别放在 index 0 和 256；
   VCF 与 BCF 两条路径（BCF 无该 u8 限制，需确认行为一致）。

## 3. V22-R2-2：整数饱和与缺失哨兵

### 复现

```text
DP=5000000000   --max-dp 100  -> 不掩码（应掩码）
DP=4294967295   --max-dp 100  -> 不掩码（恰好等于 INT_MISSING）
DP=50                         -> 正常掩码
```

`--hist dp` 只统计了 `50`；`5000000000` 和 `4294967295` 都不在直方图里。

AD 求和：

```text
AD=2147483647,2147483647,2147483647 -> 饱和到 u32::MAX = missing -> 不掩码
AD=2147483647,2147483647             -> 和 = u32::MAX-1 -> 正常掩码
```

### 根因

- `genotype_scan.rs:453-465` `parse_int` 用
  `saturating_mul(10).saturating_add(...)`，溢出到 `u32::MAX`，恰好等于
  `INT_MISSING`。
- `genotype_scan.rs:472-483` `parse_int_sum` 跳过等于 `INT_MISSING` 的分量，
  再 `saturating_add`；超范围分量被静默丢弃，求和也可能撞哨兵。
- `bcf/columns.rs:120-138` `int_sum_column` 同样 `saturating_add` 到
  `u32::MAX`。

### htslib 对照

```text
DP=3000000000  -> bcftools 警告 "Extreme FORMAT/DP value ... set to missing"
DP=5000000000  -> 同上，missing
DP=2147483647  -> 保留
```

htslib 以 i32 为界：超出 i32 的文本值直接置 missing，不会饱和成某个
有效值。因此**不能只把饱和点改成 `MAX-1`**：那会让 `5e9` 变成一个巨大的
有效深度，继续偏离 htslib。

### 计划

1. `parse_int` 改为按 `i64` 解析（支持可选 `+`/`-`）：超出 `i32` 范围
   （或负数用于无符号深度字段）→ `INT_MISSING`，与 htslib 一致；
   不再出现 `u32::MAX` 正常值。
2. `parse_int_sum`/`int_sum_column` 用 `u64` 累加，最后 `min(INT_MISSING-1)`，
   保证求和结果永远不是缺失哨兵；单个 AD 分量超出 i32 → 按 htslib
   视为该分量 missing。
3. 文本/BCF 对拍增加边界用例：`i32::MAX`、`i32::MAX+1`、`u32::MAX`、
   三元素 AD 最大值。
4. 若保留 `u32` 列，需在文档中明确"超出 i32 的文本整数按 missing"的
   语义；长期可考虑 `i64`/`Option` 列。

## 4. V22-R2-3：`parse_float` 的 `inf`/`nan`

**判定：不是本轮 bug，不改。** 理由：

- GQ 在 header 中是 Integer，走 `parse_int`，`GQ=inf` 的结果是
  "非数字 → missing"，不会经过 `parse_float`。评审的 GQ 示例用错了
  解析器。
- Float FORMAT 字段确实接受 `inf`/`nan`；但 htslib 的 `strtod` 也接受：

```text
$ bcftools view -H qual_nan.vcf      # QUAL=nan / inf 原样输出
$ bcftools view -H -i 'QUAL<20'      # nan/inf 都不入选
```

因此这是与 htslib 一致的行为，不是兼容性 bug。计划只在 `parse_float`
的文档/测试中固定"`inf` 是值、`nan` 走 missing"的语义；若未来要求严格
VCF Float 语法，再另开决策（会与 htslib 分歧）。

## 5. V22-R2-4：`--hist dp` 的 AD-only 拒绝

### 复现

```text
header: FORMAT/AD Integer，无 FORMAT/DP
$ reed stats ad_only.vcf --hist dp
reed: --hist "dp": the header declares no FORMAT/DP
exit 2

$ reed stats ad_only.vcf --summary    # 正常
$ reed stats ad_only.vcf --depth      # 正常（idepth/ldepth.mean）
```

engine 层的 `IntSpec::depth()` 已经支持 AD 回退；只有 CLI 的
`translate/stats.rs:66-80` 仍要求 header 声明 `FORMAT/DP`。这与 v22
README 的"`--hist dp` 使用同一深度列"直接矛盾，BSA `GT:PL:AD` callset
无法使用该直方图。

### 计划

1. `--hist dp` 的 header gate 改为：`FORMAT/DP` Integer，或
   `FORMAT/AD` Integer（与 filter `--min-dp` 的编译期校验一致）。
2. `--hist dp` 增加真实 CLI 测试（现有 `e2e_hist` 走 engine API，漏了
   CLI gate）。
3. 同一处补 DP/AD 类型校验，见 V22-R2-6。

## 6. V22-R2-5：两套整数解析语义

### 复现

```text
DP=+50:  --hist dp  -> 空直方图（按 missing）
         --hist gq（GQ=+30）-> 30
         --min-dp 10 -> +50 不掩码（值 50 >= 10，结果碰巧对）
         --max-dp 10 -> +50 不掩码（应掩码）
GQ=+30:  --hist gq -> 30
DP=-5:   --min-dp 10 -> 不掩码（应掩码）
```

原因：filter 的 `parse_int` 只接受 ASCII 数字，`+`/`-` 一律 missing；
stats 的 `parse_int_token`（`reports.rs:1122-1128`）用 Rust `u32` 解析，
接受 `+50`、拒绝 `-5`。同一个字段在两条消费路径上语义不同。

htslib 对 `-5`/`+50` 都按整数接受：

```text
$ bcftools view -H -i 'FMT/DP<10' dp_neg.vcf
chr1 10 ... 0/1:-5      # -5 被选中
$ bcftools view -Ob ...  # +50 保留为 50
```

### 计划

1. 收敛为一个整数文本解析器（filter/stats/hist/BCF 编码共用）：
   支持可选 `+`/`-`，按 VCF Integer 语法；不再各写一套。
2. 深度字段（DP/GQ/AD）的负数语义：htslib 保留负值；reed 的 `u32`
   列无法表示。可选方案：
   - A. 深度字段负数 → missing（维持无符号域，文档声明与 htslib 的
     分歧；只解决 `+` 和饱和）；
   - B. 整数列改 `i32`/`i64`（真正支持负数比较，影响 mask/expr/hist，
     改动较大）。
   建议 A 先落 v22，B 作为后续轮（见 V22-R2-8）。
3. 测试：`+50`、`-5`、`+0`、前导零；`--min-dp/--max-dp/--min-gq`、
   `--hist dp/gq/FMT_DP` 的同一值必须给出同一结果。

## 7. V22-R2-6：stats 侧的类型校验

### 复现

```text
header DP=Float, AD=Integer:
  --hist dp / --depth / --summary 都 exit 0；Float DP 静默按缺失
header AD=Float, 记录无 DP:
  --summary/--depth/--hist dp 的 AD 回退静默缺失
header GQ=Float:
  --hist gq 只输出 #missing
```

filter 侧 V22-R7 已经拒绝 Float DP/AD 深度；stats 侧没有对应校验。

### 计划

1. 深度消费者（`--summary` DP/PSC、`--depth`、`--hist dp`）共用一条头部
   校验：DP/AD 若参与深度必须是 Integer；非 Integer → usage error。
2. `--hist gq`：GQ 必须是 Integer；否则 usage error（或明确采用 Float
   直方图，需单独决策）。
3. 测试：DP Float、AD Float、GQ Float 三个 header；filter 与 stats 的
   错误信息/退出码一致。

## 8. V22-R2-7：Float `nan` 与缺失哨兵

`parse_float` 把 `nan` 解析为 `f32::NAN`，列缺失也用 `f32::NAN`；
下游无法区分"合法的 nan 值"与"缺失"。病态输入；当前 `--hist` 的
`parse_float_token` 已把非有限值归为 missing，filter 表达式则保留
NaN 参与比较（比较恒 false）。

计划：低优先级。若要修，需把浮点列改为 `Option<f32>` 或使用 BCF 的
特定位型作缺失哨兵；影响 columns/expr/hist/BCF 渲染，建议单开一轮。

**v22-R2-D 结论（2026-11-03）**：本轮不改，只在 `parse_float` 文档与
`genotype_scan::tests::parse_float_accepts_htslib_inf_and_nan` 固化
"`inf`/`-inf` 是值、`nan` 是 NaN 且与缺失不可分"的语义（文本与 BCF
两条路径一致）。后续轮的设计选项：

| 选项 | 做法 | 代价 |
|---|---|---|
| A | 浮点列改 `Option<f32>` | 文本/BCF 语义最清晰；每值 +4 字节，`ScanOutput`、表达式求值、hist、BCF 渲染全要改 |
| B | 列内用 BCF 的 missing/vector-end 位型（`0x7F800001`/`0x7F800002`）作哨兵 | 保持 4 字节；文本 `nan` 与缺失仍不可分，渲染需特判 |
| C | 维持 `NaN` 哨兵并声明限制（现状） | 零改动；合法 `nan` 被当作缺失 |

## 9. V22-R2-8：负数/有符号 FORMAT Integer（后续轮）

`FMT/<signed-tag>` 的负值在 `parse_int` 里是 missing，在 `int_element_column`
里也被 `value < 0` 归为 missing。若要让 `-10 < -5` 这类表达式和
`--hist FMT/XYZ` 正确，需要把整数列从 `u32` 改为 `i32/i64`。这超出
v22 的"深度回退"范围，记录为后续轮，不在 v22-R2 首批实现。

## 修复批次（完成情况）

| 批次 | 内容 | 优先级 | 状态 / 验收 |
|---|---|---|---|
| V22-R2-A | R2-4 `--hist dp` AD-only gate + CLI 测试；R2-6 stats 类型校验 | P0 | **完成**：AD-only 出直方图；Float DP/AD/GQ 报 usage error（exit 2） |
| V22-R2-B | R2-2 整数解析 i32 边界、`+`、求和防哨兵；R2-5 单一解析器 | P0 | **完成**：边界值文本/BCF 一致；`+50` 两条路径同值；AD 三最大值不撞哨兵 |
| V22-R2-C | R2-1 子字段索引（`usize`/`u16` 或 >255 明确报错），含 GT 编辑两处 | P1 | **完成**：索引全链路 `usize`；GQ/AD/GT 在 index 256+ 正确解码，GT 掩码生效 |
| V22-R2-D | R2-3 文档/测试固化 inf/nan；R2-7 NaN 哨兵（低优先级） | P2 | **R2-3 完成**（测试固化）；R2-7 只记录设计选项（§8） |
| V22-R2-E | R2-8 负数/有符号整数列（后续轮） | 后续 | 未开始；设计约束见 §9 与 06-roadmap |

## 未做

- R2-7（Float NaN 哨兵）与 R2-8/E 批（有符号整数列）留给后续轮，
  设计选项/约束已记录在 §8、§9 与 [06-roadmap.md](06-roadmap.md)。
- `inf`/`nan` 的严格 VCF 语法问题不在 v22-R2 修复范围，保持 htslib
  兼容。
