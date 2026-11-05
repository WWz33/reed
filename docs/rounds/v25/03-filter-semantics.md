# 03 filter 求值语义（V25-04）

> bcftools 的行为（missing_logic 真值表、AC/AN 的 INFO-vs-GT 偏好、cohort 门控）
> 是需要被**调用**的函数，不是需要逐处转录的散文知识。本轮建 `semantics` 模块，
> 让文本表达式、typed INFO 规则、列式 mask、head/tail AF 全部编译成对它的调用。

---

## 1. 表达式路由与 info_trusted

### 现状【据】

- `plan/mod.rs:562/577/881` 用 `expr.references_genotypes()` 决定 head/tail 路由；
- `parser.rs:240-253` 的 match 只覆盖 `Field / Bare / Cmp / Logic / Fmt / Not /
  Subscript`，**`Expr::Info` 落入 `_ => false`**（:251）；
- `head.rs:135`：`SiteValues::from_line(line, pos1, None, true)` 硬编码信任 INFO，
  注释自称 "true is exact here"（引 doc_v20 S58）——对 `Expr::Info` 不成立；
- plan 自己的 `info_trusted` 在 cohort 变化时已被清（`plan/mod.rs:1045`【述】）、
  tail 侧确实执行 S58 门控——**只有 head 表达式漏了**。

### 触发 / 后果

100 样本 VCF，记录带 `INFO/AC=80`（全队列），运行 `-s s1 -e 'INFO/AC>1'`：

- bcftools 在 `bcf_subset` **之前**求值表达式（`vcfview.c:322` 先于 `:341`），
  用全队列的 AC=80 判定 → 保留该记录；
- reed 的 head 规则也在 mask/子集之前，但当时把 `INFO/AC` 路由到了 tail，
  于是读子集后重算的 AC=1 → 记录被错误丢弃；
- 裸写 `AC>1` 解析为 `Expr::Field(Ac)` 路由到 tail，读子集值——这条是
  pre-v25 就与 bcftools 不同的旧分歧，本轮不改口径；
- 与自家单测 `s58_untrusted_info_recomputes_from_genotypes` 宣称的行为自相矛盾。

**2026-11-05 复核结论**：`INFO/` 前缀是纯 INFO 读取，留在 head（与 bcftools
同队列）；裸名仍需基因型解码，留 tail。`references_genotypes` 恢复该区分
（V25-04 曾把两者合一，是回归）。`reads_gt_fields` 保持两种拼写都拒——`--set-gts`
的掩码决策不能依赖基因型或过期 INFO 值，与该路由问题无关。

### 规范【规】

1. `references_genotypes` 删除或委托 `expr_reads_gt_fields`（后者已正确处理
   `Expr::Info`【述】）；plan 三处路由点只保留一棵语义树，路由判定提为 `Expr`
   的关联方法，parser 与 plan 共用；
2. head 求值函数签名带上 `info_trusted: bool`，一律取 plan 单字段；
   `from_line(..., true)` 硬编码与 "true is exact" 注释一并消失；
3. 【荐】`-i/-e 'INFO/AC…'` 在 `cohort_changed` 时同样送 tail（或编译期报错
   引导改用裸名）——两棵语义树二选一，不留"路由看拼写"的暗坑。

### 测试

- S58 用例扩展：`-s` + `INFO/AC / INFO/AN / INFO/AF / INFO/MAF / INFO/MAC` 五个
  基因型衍生 built-in 的 head 路由用例（cohort 变化前后各一）；
- 对拍 bcftools：子集 + `-i/-e 'INFO/AC>1'` 的存活位点集一致。

---

## 2. compare_missing：缺失值比较真值表

### 现状

bcftools 的 missing_logic（`!=` 在缺失时通过、`==` 需两侧缺失、有序比较失败）
三处独立实现（三处存在已【据】确认，行号为【述】）：

| # | 位置 | 语境 |
|---|---|---|
| 1 | `expr/mod.rs:525` `any_pair` | 文本表达式（`-i/-e`） |
| 2 | `eval/head.rs:301` `missing_value_passes` | typed INFO 规则（`-q/-Q` 短路） |
| 3 | `reed/src/filter/record_pipeline.rs:524-560`（INT_MISSING/NaN + `on_missing`） | 列式 mask（`--min-dp` 等） |

主张的分叉用例【述，先 failing 复现】：`--set-gts . -e 'FMT/DP!=3'`（样本 DP 缺失 →
`any_pair` 判 true → mask 生效）vs 等价的
`mask_if [{field:"DP", when:"!= 3"}]`（缺失被跳过 → 不 mask）——同一谓词两种拼写
对缺失样本结果相反。V23-59 只统一了 ①②，③ 的 mask 语境没有 missing 真值表，
只有 `on_missing` 覆盖。

### 规范【规】

reed-filter 新建 `semantics` 模块，缺失值比较只有一份：

```rust
/// bcftools missing_logic + mask 语境的 on_missing 覆盖，单一定义。
pub enum MissingCtx { Expr, InfoRule, ColumnarMask(OnMissing) }
pub fn compare_missing(op: CmpOp, lhs: Option<Val>, rhs: Val, ctx: MissingCtx) -> bool;
```

- 三个调用点全部改走；
- mask 语境（③）的缺省值保持现行为（缺失不触发 mask，除非 `on_missing=mask`），
  与表达式语境（①）的差异**写进 doc**——差异本身是历史契约，不是 bug；
- 真值表全组合（op × lhs 缺失/有值 × ctx）用单测钉死，附 bcftools `filter.c`
  的行号引用。

### 测试

- 先写 failing 用例复现上面的分叉主张——复现则修，不复现则把"等价"钉进测试
  并关闭该子项（R3 处置完毕）；
- 真值表单测 + V23-59 既有用例保持绿。

---

## 3. af_source：AC/AN/AF 的 INFO 优先、GT 回退

### 现状【述】

| 位置 | 实现 | 畸形防护 |
|---|---|---|
| `eval/head.rs:161-229` `allele_counts_info` / `non_reference_af_info` | head 短路（`-q/-Q`） | **无**：`AC=200;AN=100` 算出 freq=2.0 直接判 |
| `eval/tail.rs:108-158` `af_of` | tail（`-m/-M/-v/-V` 的 MAF/MAC） | 有：`alt_sum > an → None` |

两份实现 + 不对称防护：同一畸形记录 head 判 freq=2.0、tail 判"无可用 INFO"，
结论分道。D12-6 / D12-13 / S58 每次补丁都要两处同改。

### 规范【规】

```rust
/// INFO 的 AF/AC 读取 + 畸形防护，单一定义；head/tail 只带策略包装。
pub fn info_af(info: &InfoFacts, n_alt: usize) -> Option<(f64, AltCounts)>;
```

- 读取与防护只有一份；head=None 放行、tail=None 拒绝的策略差异保留（这是
  两阶段设计的本意），但差异点显式写在两个包装里；
- 【荐】head 的 `-Q` 上界对 freq>1 的处理与 tail 对齐（同一个畸形记录同判）。

### 测试

- `AC=200;AN=100` fixture：head `-q` 与 tail `-m` 同判（都当"无可用 INFO"）；
- 既有 D12-6 / S58 用例保持绿。

---

## 4. V25-04b：`when` 的 `<=`/`>=` 回放（R2）

### 现状【据】

`spec/render.rs:126-143` 从 `rule.when` 的字符串前缀（`<`/`>`）决定回放成
`--min-dp` / `--max-dp` / `--min-gq` / `--max-gq`。V23-24 只修了"`<=` 的 `=`
留在值里"（`strip_prefix("<=")`），**闭区间语义仍丢失**：

- 规则 `DP <= 3`（mask DP≤3，保留 DP>3）回放成 `--min-dp 3`（保留 DP≥3）；
- DP=3 的样本在直接求值与 spec→CLI 回放两条路径上结果相反。

同一谓词两条路径、两种结果——P1，不是格式问题。

### 规范【规】

1. `when` 不再用裸 `String` 承载比较语义，存 `(CmpOp, Value)`（或 `Cmp`）；
2. 渲染时由 `CmpOp` 决定可回放的 flag；整数字段把闭区间换算成等价的严格边界
   （`DP <= 3` → `--min-dp 4`），无法精确表达时（非整数域/超出范围）显式报
   "not replayable"（usage），而不是静默降级成开区间；
3. 解析/渲染/求值共用一个 `CmpOp` 映射；JSON 输入兼容旧的 `"when": "<= 3"`
   （解析仍接受），但回放不再假设前缀就是全部语义。

### 测试 / 验收

- proptest：`(CmpOp, Value)` 的 render→parse round-trip 恒等；
- e2e：`DP <= 3` 的 spec 回放与直接 `mask_if` 求值在同一 fixture 上 mask 集
  相同（必须含 DP=3 的边界样本）；
- `spec/render.rs` 不再用 `starts_with('<')/starts_with('>')` 推断比较语义；
- V23-24 的既有 round-trip 用例保持绿。

---

## 5. 拍板：mask 归因顺序（R10）

single-pass（表达式先）与 sketch（声明序）的归因顺序不一致，同一 mask 在两条
策略上可能记到不同 funnel 项，破坏 "byte-identical funnels" 契约。

- **A（荐）**：compile 时固定为声明序，两条策略共用；若改动限于
  `record_pipeline.rs` / `mid_stage.rs` 的循环顺序，并入 V25-04；
- **B**：保持现状，在文档登记"归因顺序不是契约"，并加一条测试钉住差异。

若选 A 且需要改 chunk 结构或输出数值，转 v26，不在本轮硬塞。

---

## 验收（V25-04 整体）

- `INFO/AC` 与裸 `AC` 在 cohort 子集后同队列（对拍 bcftools，存活位点集一致）；
- **三策略等价**：同一 fixture 跑 single-pass / sketch / tail，`-i/-e` 对
  `INFO/AC|AN|AF|MAF|MAC` 与缺失比较（`=`/`!=`/`<`/`<=`/`>`/`>=`）断言存活
  位点集逐字节一致；只测 head 路由不够（F4 的根因是"路由看拼写"）；
- 缺失值比较实现点 3→1（`compare_missing`）；AC/AF 读取 2→1（`info_af`）；
- 路由语义树 2→1（`references_genotypes` 消失或退化为委托）；
- R2 的 round-trip 与边界 e2e 全绿；R10 的 A/B 已在 STATUS 记录；
- `head.rs` 的 "true is exact" 注释与硬编码一起消失；
- `grep 'missing_logic\|bcftools' crates/reed-filter/src/` 的语义引用都指向
  `semantics` 模块一处。

---

## 实施记录（2026-11-04，提交 `00e3a23`）

### 1. 路由与 info_trusted

- `Expr::references_genotypes` 增加 `Expr::Info` 分支（与 `Bare` 同一判定）；
  `expr_reads_gt_fields` 从 `plan/mod.rs` 移到 `Expr::reads_gt_fields`。路由判定
  只剩 `parser.rs` 一棵树，两方法的唯一差异是 `Fmt`（掩码阶段 vs 站点字段），
  已写在文档注释里。
- `eval_head` 增加 `info_trusted: bool` 参数，调用点传 `plan.info_trusted`；
  `from_line(..., true)` 与 "true is exact here" 注释删除。
- **路由选择**：`INFO/AC` 等五个 built-in 在 `[site]` 下与裸名一样被编译期拒绝
  （"belongs in [site_after_mask]"），因为它们的值依赖队列；`INFO/DP` 等普通
  INFO tag 仍留在 head。`[site_after_mask]` 接受两种拼写。测试覆盖五个名字的
  两种拼写（AC/AN/AF/MAF/MAC）与 cohort 子集前后的判定。
- **对拍**：子集 `{s0,s1}` + `AC>2` / `INFO/AC>2` 的记录集相同（都为空，
  recomputed AC=1），全队列 `INFO/AC>2` 保留（INFO AC=4）。

### 2. compare_missing

- `reed_filter::semantics::{compare_missing, MissingSides, MissingCtx}`；三个调用点
  （`any_pair`、`missing_value_passes`、`record_pipeline` 列式 mask）全部改走。
- 真值表单测覆盖 6 个 op × lhs/rhs/both × 3 个 ctx；`Op -> CmpOp` 的映射也在
  semantics 里。
- **R3 处置**：`-e 'FMT/DP!=3'`（表达式语境，缺失满足 `!=` → mask）与
  `mask_if DP != 3`（默认 `Keep` → 缺失不 mask）的差异**复现且是契约**；e2e
  `missing_dp_masks_in_the_expression_context_but_not_by_default_in_mask_if`
  逐样本钉住（DP=3/缺失/4 三种）。同语境内部没有第二套实现。

### 3. info_af

- `semantics::info_af` 读 AC/AN，填 REF=AN−ΣALT，带 `AC` 长度检查（一次性
  warning）与 `ΣAC > AN` 防护；head `non_reference_af_info` 与 tail `af_of`
  都调它。
- **行为变化（有意）**：畸形 `AC=200;AN=100` 在 head 不再算出 AF=2.0；tail 不再
  直接拒绝，而是把它当"无可用 INFO"回退到 genotype counts（与缺失 tag 同径）。
  单测 `malformed_ac_an_is_no_usable_info_for_the_head` /
  `malformed_info_falls_back_to_the_counts_on_the_tail` 钉住。

### 4. R2 回放

- `spec/render.rs` 用 `rule.parsed()` 的 `(Op, Typed)`，不再嗅探首字符；闭区间
  换算成等价严格边界（`DP <= 3` → `--min-dp 4`，`DP >= 20` → `--max-dp 19`）。
- 无法精确表达的形态（`==`/`!=`、小数界）沿用 V22-R6 的"省略"策略，而不是
  静默降级；`to_args` 保持不可失败（provenance 用它），因此没有引入 usage
  错误（与 03 §4 规范的"显式报错"有偏差，记在此）。
- CLI 测试 `closed_mask_bounds_replay_without_changing_the_mask` 验证 DP=0..6
  的 mask 集在 render→parse 前后逐值相同。

### 5. 三策略等价与 R10

- e2e `three_strategies_agree_on_info_builtins_and_missing_comparisons`：
  - INFO built-ins（AC/AN/AF/MAF/MAC）在 SinglePass 与 TwoPassReparse 上逐字节
    相同；
  - 六个缺失比较（`=`/`!=`/`<`/`<=`/`>`/`>=`）在 SinglePass、
    TwoPassReparse、SketchThenEmit 三策略上逐字节相同；
  - **sketch 对"读记录"的 tail 表达式仍按既有契约拒绝**（mid stage 无记录
    字节）。这不是本轮引入的分叉，而是架构限制；测试断言它给出
    "no record bytes" 的 usage 错误，而不是换一个 INFO 源算错。若 v26 要让
    sketch 也支持，需要在 pass1 重算 tail（chunk 结构变化），不属于本轮。
- **R10 拍板：分叉不成立，按 B 落地**。两条策略共用同一个 `RecordPipeline`
  的 `mask_rule`（首个执行序命中），mid stage 的 `rule_hits_into` 再按声明序
  折叠的是同一份 bitmap；e2e
  `mask_attribution_is_strategy_invariant_with_overlapping_rules` 用
  "一个样本同时命中 field 与 expr 规则"的 fixture 断言两策略的
  `funnel`/`n_masked_gt`/`masked_by_rule` 全部相同。归因顺序因此**是**契约
  （首个执行序命中：表达式先、列式 field 后），不是"不确定"。
