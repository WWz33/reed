# 02 修订规范（N1…N14）

> 每条 N 对应 01-review.md §1 的一条或多条缺陷（D12-x），给出**可执行的规范**
> 与**验收标准**。规范只约束"INFO 字段的语义与编码该长成什么样"，量尺是
> 工作区 htslib/bcftools（1.24）与 `docs/spec/cli.md` 的既有承诺；不重开
> v6 的 Y 系列、v10 的结构项。
> 标记：【规】必须；【荐】推荐；【决】需拍板。

---

## §1 策略裁决单点化：尾部表达式否决 sketch（N1 ← D12-1）【规】

策略的**最终**裁决者是驱动（`reed/src/filter/mod.rs:232-265`），因此守卫必须
长在驱动里，不能只留在 `plan::compile`：

1. `FilterPlan` 暴露一个谓词（如 `pub fn has_tail_expr(&self) -> bool`，
   即现有 `plan/mod.rs:740` 的表达式），作为唯一事实源；`compile` 内部改用
   同一谓词。
2. 驱动的 `--mind` 分支：`choose_for_mind` 返回 `SketchThenEmit` 且
   `plan.has_tail_expr()` 时，改选 `TwoPassReparse`（pass1 在记录字节上求值，
   语义与 SinglePass 一致，实测正确）。
3. `--mind + --indep-pairwise`（sketch 强制路径）且 `plan.has_tail_expr()`：
   编译期 `Error::usage`——"尾部表达式与 --mind --indep-pairwise 组合暂不支持
   （sketch 中段没有记录字节）；去掉表达式或去掉 --indep-pairwise"。
   **不静默选错**。
4. 尾部 `-q/-Q`（`MinAf/MaxAf`）读 `values.info`，同属"行字段依赖"：
   第 2 条的门用 `has_tail_expr() || tail 含 MinAf/MaxAf 且 header 声明
   AC+AN` 作判据（sketch 下 INFO 不可见会静默改变 `-q` 口径，见 D12-1 实测）。
   N5 落地后（cohort 门下 `-q` 不再读 INFO），此判据可收窄回 `has_tail_expr`。

**验收**：
- e2e：`--mind 0.9` + `-i 'F_MISSING<=0.5 && QD>10'` 与 SinglePass 输出
  逐字节一致（修复前 sketch 为 0 条即红）。
- e2e：`--mind 0.9 --indep-pairwise … -i '<尾部表达式>'` 退出码 2、报错文本
  含两个参数名。
- 既有 `all_three_strategies_produce_identical_bytes` 的夹具追加一条含 INFO
  读取的尾部表达式规则后仍三策略一致（含 N5 的 `-q` 场景）。

---

## §2 重算门 `has_gt`：无 GT 的记录不动 AC/AN（N2 ← D12-2）【规】

bcftools 的两处同款门（`vcfview.c:343`、`fill-tags.c:763`）移植为 reed 的
逐记录规则：**`counts.has_gt == false` 的记录，AC/AN 既不重算也不从 stale
清单里删除——输入值原样保留**。

1. `fill_info_set`（`eval/edit_plan.rs:21-46`）：`recompute.ac/an` 分支加
   `counts.has_gt` 门（NS/MAF/HWE 同门——它们同样是 GT 派生量，无 GT 时
   bcftools fill-tags 整条记录跳过）。
2. emit 侧（`single_pass.rs:306`/`emit.rs:115` 的 `stale_refs` 消费点）：
   当 `!view.counts.has_gt` 且 stale 中的 AC/AN 来自**自动更新**（区分于用户
   `-x` 显式指名）时，该记录不删这两个标签。实现上给 `EditPlan` 的调用点
   传"本记录生效的 stale 子集"即可，`EditPlan` 形状不变。
   用户显式 `-x INFO/AC` 的删除**不受** `has_gt` 门影响（`annotate -x` 语义
   无条件）。
3. sketch/mid 路径：`matrix.has_gt(local)` 已存在（`mid_stage.rs:324`），
   `SiteEmit.infos` 的填充走同一门。

**验收**：
- 红→绿：`nogt.vcf`（`AC=4;AN=6`、FORMAT=DP）`-s s0,s1` → 输出 `AC=4;AN=6`
  （修复前 `AC=0;AN=0`）。
- 有 GT 记录的自动更新与 `-I` 行为不变（现有 e2e 全绿）；
  `-x INFO/AC` 显式删除在无 GT 记录上仍生效。

---

## §3 编码器整数值域：超界置 missing（N3 ← D12-3）【规】

htslib 解析期规则（`vcf.c:4004-4013`）移植到 reed 的**写入**期：

1. `parse_int_vector`（`native_encode.rs:1115-1123`）：解析成功但
   `v < INT32_MIN_VALUE || v > INT32_MAX_VALUE`（`record.rs:59-61` 常量，
   即 htslib `BCF_MIN/MAX_BT_INT32`）→ 推 missing 哨兵（`i64::MIN`，
   `encode_int_vector` 已按哨兵写各宽度的 missing）。
2. 一次性 `tracing::warn!`（"Integer INFO value out of BCF range set to
   missing"，含 tag 名与首个越界值；每 worker 一次，对齐 htslib 的
   `extreme_int_warned`）。
3. `encode_info_typed` 的 `(Integer, Int(v))` 分支（`:997`）同规则——重算值
   来自 u64 计数，`as i64` 后超界的理论路径一并封死。

**验收**：红测试 `out_of_range_integers_become_missing_in_bcf`
（`XD=3000000000`/`-3000000000`/`-2147483641` 三值 → bcftools 读回均为
`XD=.`）；真值门：与 bcftools 自产 BCF 的读回文本一致；界内值
（`2147483647`/`-2147483640`）不变。

---

## §4 值 × 声明类型的忠实编码（N4 ← D12-4/D12-5）【规】

"一个 INFO 值（源文本字节 或 `InfoValueOut`）× header 声明类型 → BCF 字节"
收成**单一函数**，`encode_info_text` 与 `encode_info_typed` 共用（消除簇 γ
的两份实现）。规则表（htslib `vcf_parse_info` 同款）：

| 声明类型 | 源形态 | 编码 | 依据 |
|---|---|---|---|
| Flag | 无值（裸标签） | `(NULL, 0)` | vcf.c:3947-3948 |
| Flag | **带值** | 按字符串存字面值（`bcf_enc_vchar`） | vcf.c:3949-3951 "if Flag has a value, treat it as a string" |
| Integer/Float | 有值 | 解析（N3 值域 + N13 解析规则） | vcf.c:3990-4056 |
| Integer/Float | 无值（裸标签，畸形但存在的输入） | 单元素 missing | 现状保留（`native_encode.rs:365-372`） |
| String/Character | 有值 | **字面存**，`.` 不特判 | vcf.c:3949-3951 |
| String/Character | 无值 | 零长度字符串 | 现状（`:373`） |

1. 删除 `(InfoType::Flag, _) => NULL` 的无条件分支（`:346`/`:996`）——
   只有"无值"才走 NULL。
2. 删除 `if v == b"." { &[] }` 特判（`:362-364`）——INFO 值的 `.` 是字面
   字符串（整个 INFO 列为 `.` 的缺省记号在更外层已处理，
   `native_encode.rs:308` 的 `info != b"."` 判断）。ID 列的 `.` → 空串
   特判（`:263`）**保留**（htslib 对 ID 确实如此）。
3. `InfoValueOut`（重算值）无 Flag 形态，`encode_info_typed` 的 Flag 分支
   仅在"用户把重算标签声明成 Flag"的畸形组合下触发——维持现状（写 NULL），
   由 N11 的文档说明兜底。

**验收**：红测试三例——`DB=1`、`DB=0`、`SOMESTR=.` 经 `-Ob` 后 bcftools
读回与输入逐字节相同（修复前分别为 `DB`、`DB`、`SOMESTR`）；`-Ov` 输出不变；
真值门：与 bcftools 自产 BCF 读回一致。

---

## §5 `-q/-Q` 的 cohort 门与"陈旧 INFO"警告（N5 ← D12-6）【规】

bcftools 的门（`vcfview.c:341-350`、`filter.c:1471`）：**样本集或基因型被
改动后，INFO/AC、AN 不再可信，`-q` 判定改用（子集/掩码后的）GT 计数**。

1. `FilterSpec` 增加派生谓词（编译期可得）：
   `cohort_changed() = sample_select.is_some() || mind.is_some() ||
   !genotype.mask_if.is_empty() || genotype.set_gts.is_some()`。
2. `compile`（`plan/mod.rs:382-413`）：`cohort_changed()` 为真时**不推**
   `HeadRule::MinAfInfo/MaxAfInfo`（head 的 INFO 短路整层失效，`-q` 全部
   落到尾部 GT 路径——正是 bcftools 的行为序）。
3. `af_of`（`tail.rs:108-134`）：`values.info` 的 INFO 优先分支加同一门
   （门经 `SiteValues` 或 `eval_tail` 参数传入；sketch 的 `from_counts`
   上下文 `info=b"."` 天然走 GT，行为自此三策略一致，与 N1 §4 呼应）。
4. **警告兑现**（`cli.md:78`、v4/02:103）：`!cohort_changed()` 且 header
   声明 AC+AN 时，编译期 `tracing::warn!` 一次——"`-q/-Q` reads INFO/AC,AN,
   which may be stale; use --maf/--mac for genotype-recomputed values"
   （文案与 `warn_stale_info` 同族）。
5. 输出记录的 AC/AN 自动更新（N2 机制）不变——判定与写出自此口径一致，
   "保留依据 AF=0.5、输出写 AF=0"的自相矛盾消失。

**验收**：红测试 `min_af_after_subset_uses_the_genotypes`
（`sub.vcf` `-s s2,s3 -q 0.1` → 0 条，修复前 2 条）；`e2e_af_type.rs`
（无子集，INFO 优先）全绿；`--mind` 场景三策略一致；`-q` + AC/AN header
→ stderr 恰一次警告。

---

## §6 表达式层的类型感知（N6 ← D12-7）【规】

`check` 与求值上下文接 `&TagTable`（已存在，`reed-core/src/tag_table.rs:57`；
`compile` 已持有），取代裸名字列表：

1. **Flag 比较**：`INFO/DB=1`、`DB=1`、`DB!=0`、`DB=0` 按存在性语义求值
   （存在→1、不存在→0，bcftools `filters_set_info_flag`，`filter.c:964-971`；
   文档依据 `doc/bcftools.txt:4067`）。`info_value` 对 Flag 类型标签返回
   `Num(1)/Num(0)` 而非 Missing；`Val::Missing` 仅留给"值缺失"
   （`DP=.`）。裸 `DB`（无比较）维持 false（见第 3 条）。
2. **向量比较**：Number ∈ {A, R, G, .} 或 `Count(>1)` 的标签，比较按
   **any-element**（bcftools `CMP_VECTORS`，`filter.c:2738-2790`：任一对
   满足即真）。`info_value` 返回元素集合或求值处按逗号展开；`AC/AN/AF/MAC/MAF`
   内建字段的 INFO 优先路径（`info_or_gt`，`expr/mod.rs:365-376`）取**首元素**
   维持现状（bcftools 的内建 AC 同样是向量语义，但 reed 的内建字段契约是
   标量——差异登记 cli.md §5，不在本轮扩）。
   **2b. Number=1 声明却给了多值**（畸形输入，如 `XI=5,`）：bcftools 的标量
   setter 读到 **0**——htslib `bcf_unpack` 先把 `v1.i` 初始化为 0，且只在
   `len == 1` 时填充（`vcf.c:4286-4287`），于是多值记录的标量视图是 0 而
   **不是 missing**（实测可区分：`-i 'XI=0'` 命中 `XI=5,` 与 `XI=1,2,3`，
   `-i "XI='.'"` 不命中）。这是 htslib 的实现细节而非设计意图。
   【决】reed 的处置二选一：
   - **(a) 对齐 bcftools**：读作 0。好处是与 bcftools 逐条一致（对拍门可断言）；
     坏处是把一个畸形输入的偶然值固化成语义，且 0 会参与 `XI<5` 之类的比较
     产生反直觉结果。
   - **(b) 当作 missing**（推荐）：`len != 1` 而声明 Number=1 → `Val::Missing`，
     比较恒假（与"缺失值比较为假"的既有规则一致，`cli.md §5`）。好处是语义
     自洽、不固化 htslib 的偶然行为；代价是与 bcftools 在这类畸形输入上分歧
     （需在 cli.md 差异表登记，同 Σac>an 的处理方式）。
   默认选 **(b)**；若用户要求逐条对拍通过，改 (a) 并在差异表删除该行。
3. **裸标签真值**：对齐 bcftools——单个裸 token（无比较符）恒假
   （`filter.c:4352-4380`：`pass_site` 无比较节点不置位）。删除
   `expr/mod.rs:286-292` 的"非零为真"分支，替换为 false + 文档注明；
   或（更友好、仍合规）编译期拒绝裸标签："a bare tag is always false in
   bcftools; write `TAG=1`（flag）or `TAG>0`"——二选一，【决】默认选
   **编译期拒绝**（错误优于静默恒假）。
4. `check` 对 Flag 标签的 `~`/`!~`、向量的下标等未实现形态，照旧
   "not supported yet"（`cli.md §5` 纪律）；新增：Flag 标签参与 `<`/`>`
   时按 bcftools（int 0/1 可比较）放行。
5. `parse_and_check`/`parse_and_check_fmt` 签名相应带 `&TagTable`
   （调用点 `plan/mod.rs:420/:435/:659` 已持有 `tags`）。

**验收**：对拍测试三组（flag3/vec/zero 夹具）——`DB=1`、`DB=0`、`DB!=0`、
`DP4>10`、`DP4<5`、裸 `DP`（含 `-i`/`-e` 两向）与 bcftools 1.24 结果一致；
`cli.md §5` 支持清单回改（Flag 比较、向量 any、裸标签处置）；
既有表达式测试全绿。

---

## §7 量词语义（N7 ← D12-8）【规】

`eval_info`（`head.rs:209-225`）：

```rust
Some(v) => {
    let ok = rule.cmp.apply_f64(v);
    any |= ok;
    all &= ok;
}
```

`any` = ∃ 元素满足比较；`all` = ∀ **可解析**元素满足比较（不可解析元素使
`all` 为假——现状保留，`None => all = false`）。文档（`predicate.rs:168-176`、
`architecture.md:1210`）已如此定义，实现向文档收敛。

**验收**：单测四组合（`>=`/`<=` × any/all，值集 `{1,2}`、`{1,20}`、
含 `.` 元素各一）；修复前 any 两例红。

---

## §8 打包计数上限（N8 ← D12-9）【规】

BCF 记录头的四个计数字段有硬上限（hts-specs §6.3 的位宽；htslib 在文本
解析期即拒绝，`vcf.c:3912`/`:4174`、`vcf_parse_format_dict2`）。reed 在
**写回前**检查：

| 字段 | 位宽 | 上限 | 超限处置 |
|---|---|---|---|
| n_info | 16 | 0xFFFF（65 535） | `MalformedRecord`（"too many INFO fields for BCF"） |
| n_allele | 16 | 0xFFFF | 同上 |
| n_fmt | 8 | 0xFF（255） | 同上 |
| n_sample | 24 | 0xFFFFFF | 同上 |

一个 helper（如 `pack_counts(n_allele, n_info) -> Result<u32, IoError>` 与
`pack_fmt_sample(n_fmt, n_sample)`），四处写回点（`native_encode.rs:336`、
`:479`、`:793`、`:851-853`）全部改走它。BCF 源路径（`encode_bcf`）的输入
计数来自源记录，编辑只会减少 INFO 或增加 `info_set` 条数——检查同样覆盖。

**验收**：红测试三夹具（65 536 INFO 片段、300 FORMAT 键、66 001 等位）——
修复前静默写出坏文件，修复后 `reed filter -Ob` 报 `MalformedRecord` 且退出
码非 0；65 535/255 边界值正常写出且 bcftools 可读回。

---

## §9 ALT 计数单点化（N9 ← D12-10）【规】

1. "ALT 列 → ALT 个数"只有一个实现：`SiteView::n_alt` 的规则
   （`site.rs:129-135`：`.` → 0，否则逗号数 +1）。提取为
   `reed_io::vcf` 的自由函数（如 `count_alts(alts: &[u8]) -> u32`），
   `SiteView::n_alt` 改调它。
2. head 的 `MinAltAlleles/MaxAltAlleles`（`head.rs:82-91`）与表达式的
   `n_alt`（`expr/mod.rs:109`）改用 `count_alts`（head 已有 `alts` 字节；
   `SiteValues::from_line` 同）。
3. `-M 1` 的"0 = 关闭"歧义：`SiteRules::max_alt_alleles` 改
   `Option<u32>`（serde 默认 `None`），translate 直接
   `Some(m.saturating_sub(1))`；plan 按 `Some(max)` 推规则（`max == 0`
   合法：只留不变位点）。`min_alt_alleles` 的 0 语义（关闭）不变。
   `to_args` 渲染同步。

**验收**：红测试 `invariant_sites_count_zero_alts`——`inv.vcf` 上 `-m 2`、
`-M 1`、`-i 'N_ALT=0'`、`-i 'N_ALT=1'` 四例与 bcftools 对拍一致；
spec 往返测试（`to_args`）覆盖 `max_alt_alleles = Some(0)`。

---

## §10 `--fill-tags` 校验（N10 ← D12-11）【规】

1. 支持集常量单点定义（如 `reed_filter::plan::SUPPORTED_FILL_TAGS =
   ["AC","AN","NS","MAF","HWE"]`）；`RecomputeFlags` 的构造改经校验函数：
   大小写不敏感匹配（bcftools `strcasecmp` 同款）并归一为大写；
   未知标签 → `FilterError::Spec`（"unsupported --fill-tags tag {t};
   supported: AC, AN, NS, MAF, HWE"）。校验落在 `compile`（库调用方同样
   受保护），translate 不另做。
2. `cli.md:117` 的承诺收敛：已实现集照列；AF/F_MISSING 移入"【待实现】"
   注记（与 `:type` 后缀当年的处理同款）。
3. 【决】AF/F_MISSING 是否随本轮实现：AF 在双等位门下 = `counts.ac /
   counts.an`（Number=A 单元素），量小；F_MISSING = `n_missing / n_samples`。
   默认**不实现**（本轮只校验 + 文档），实现与否由用户裁决后追加为 N10b。

**验收**：`--fill-tags BOGUS`/`AF` 退出码 2、报错含支持集；`--fill-tags ac`
≡ `AC`；gold 命令输出逐字节不变；`cli.md` 回改后 drift 守卫绿。

---

## §11 默认 `-x` 清单可见化（N11 ← D12-12）【规】

行为不改（gold 对拍的既定默认，v2 02 §3.5），只补可见性：

1. `cli/mod.rs` 的 `-x` help：追加"defaults to the stale list
   AC,AN,AF,MLEAC,MLEAF,ExcessHet,InbreedingCoeff,DP (the gold pipeline's
   `annotate -x`); `-x ''` keeps everything"。
2. README「INFO handling」段与 `cli.md:116` 同款补一句（含 `-x ''` 转义口）。
3. `-I/--no-update` 的 help 注明：与默认 `-x` 叠加时 AC/AN 被删且不回填
   （bcftools `view -I` 无此效果，因为 bcftools 没有默认 `-x`）。
4. 【决】默认改 opt-in（`--drop-stale-default` 或类似）不在本轮；若用户
   裁决改默认，gold 脚本已显式传 `-x`，不受影响，但既有用户脚本行为会变——
   需要 CHANGELOG 显著条目。

**验收**：`reed filter --help` 输出含默认清单；README/cli.md 回改；
drift 守卫（help 文本 grep）。

---

## §12 AC/AN INFO 读取单点化 + 长度校验（N12 ← D12-13）【规】

1. `allele_counts_info`（`head.rs:128-155`）加 `n_alt: u32` 参数：
   `ac.len() != n_alt as usize` → `None`（htslib `bcf_calc_ac` 的
   `ac_len != line->n_allele - 1` → 警告 + 回退 GT，`vcfutils.c:56-64`）。
   首次触发 `tracing::warn!` 一次（"Incorrect number of AC fields"，
   htslib 同款文案族）。
2. `non_reference_af_info`（`head.rs:161-188`）改为基于同一扫描的求和视图
   （或删除，由 `allele_counts_info` 的结果求和——单一实现，簇 δ 收口）。
3. `summary.rs:288-320` 的第三份实现换用同一函数（其现有
   `vals.len() == n_alt` 与 `an >= Σac` 校验并入共享实现）。
4. Σac > an 的处置差异（bcftools `exit(1)` vs reed 回退/丢弃）：维持 reed
   现状，差异写入 cli.md §5 差异表（01-review §2 已核减）。

**验收**：红测试 `ac_length_mismatch_falls_back_to_genotypes`（aclen3 夹具
与 bcftools 对拍）；stats 的 AF/SiS 现有测试全绿；生产代码中
`for (tag, value) in info_fragments` 手工匹配 `b"AC"` 的循环仅剩共享实现一处。

---

## §13 INFO 数值解析单点化（N13 ← D12-14）【规】

> 初审列为【荐】；31 形态 × 3 路径的全量差分（01-review D12-14）显示后果是
> **reed 内部自相矛盾**——同一条记录的过滤判定值、`-Ob` 写出值、`-Ov` 写出值
> 三者不同。内部一致性是正确性问题，不是风格问题，故升为【规】。
> "是否完全对齐 htslib 前缀语义"仍是【决】（第 4 条）。

**1. 单一实现 `parse_info_number`**（放 `reed-io::vcf::info`，与
`info_fragments` 同文件——它是 INFO 值的唯一解释点）：

```rust
/// 按 header 声明的类型解析一个 INFO 值元素。
/// `Ok(None)` = 该元素缺失（`.`、空、不可解析）。
pub fn parse_info_number(text: &[u8], ty: TagType) -> Option<f64>;
```

规则（三种模式，全部**不 trim**——trim 与否是当前矛盾的成因之一；数值
文法不含空白，VCFv4.5 §1.3 的 Float 正则首尾均无 `\s`，Integer 是
"32-bit, signed"）：

| 声明类型 | 规则 |
|---|---|
| Integer | 十进制，可选前置 `-` **或 `+`**（对称于规范 Float 正则的 `[-+]?`，tex:97——`+5` 是合法输入，不是畸形）；越界处置见 N3；其余不可解析 → `None` |
| Float | Rust `str::parse::<f64>` 语法（与规范正则等价：接受 `+5`/`5.0`/`1e3`/`inf`/`nan`，拒绝 `0x10` 与空白）；不可解析 → `None` |
| Flag | 无值 → `Some(1)`（存在性语义，N6）；带值 → 按 String 处理（规范外输入，htslib 恢复策略，N4） |
| String/Character | 不解析为数值（表达式层按字符串比较）；**值可含空格**（VCFv4.5 §1.6.1 tex:431，4.5 新增——String 路径不得 trim/拒绝空格） |

2. **全部消费者改走它**（删除各自的本地解析）：
   - 编码器 `parse_int_vector`/`parse_float_vector`（`native_encode.rs:1115`/`:1126`）；
   - 表达式 `info_value`（`expr/mod.rs:380-407`）——**并且不再退回字节序
     比较**：数值型标签解析失败即 `Val::Missing`（比较恒假），只有 String
     型标签才走字符串比较。这条消除 `abc > 4` 为真的错判（01-review D12-14
     路径 C 第一行）。
   - `Val::as_num`（`expr/mod.rs:153-159`）删掉 `trim()`；
   - head 的 AC/AN 扫描（经 N12 的共享实现）；
   - `summary.rs:296-311`（经 N12 的共享实现）。
   FORMAT 的 `parse_int`（`genotype_scan.rs:421`，u32 饱和 + 缺失哨兵）语义
   不同，**不并**，但两处文档互相引用。

3. **`cli.md §5` 增加"INFO 数值解析：与 htslib 的差异"表**（N13 落地后 reed
   三条路径同值，表中 reed 列即 `parse_info_number` 的行为）。
   **处置**列是这张表的承重部分：**「对齐」= 本轮必须改成与 htslib 相同；
   「不对齐」= reed 的行为更好，永久保持分歧，只登记不改代码；「决」= 待拍板**。
   判定依据见 01-review §2.5。

   | 输入 | htslib（Integer） | htslib（Float） | reed（N13 后） | 差异原因 | 处置 |
   |---|---|---|---|---|---|
   | `5` | `5` | `5` | `5` | — | 一致 |
   | `+5` | `5` | `5` | `5`（整数也接受前置 `+`） | 规范 Float 正则显式允许 `[-+]?`（tex:97）；htslib 同 | **对齐**（N13.1 整数模式接受 `+`，一行改动） |
   | `-5` | `-5` | `-5` | `-5` | — | 一致 |
   | `0005` | `5` | `5` | `5` | — | 一致 |
   | `010` | `10`（十进制） | `10` | `10` | — | 一致 |
   | `5.0` | `5`（前缀解析） | `5` | Integer missing / Float `5` | htslib 停在 `.` | **不对齐** |
   | `1e3` | `1`（前缀解析） | `1000` | Integer missing / Float `1000` | htslib 停在 `e` | **不对齐** |
   | `0x10` | `0`（前缀解析） | `16`（`hts_str2dbl` 认十六进制） | missing | htslib 前缀解析 | **不对齐** |
   | `5abc` | `5`（前缀解析） | `5` | missing | 同上 | **不对齐** |
   | ` 5` / `5 ` | missing / `5` | `5` / `5` | missing（数值文法不含空白，两侧都不 trim） | htslib 浮点容忍前空格 | **不对齐**（reed 三条路径统一为不 trim 是必须；规范数值文法无空白，` 5` 属规范外——注意 §1.6.1 允许的是**值**含空格，主要指 String） |
   | `.` / 空 / `..` | missing | missing | missing | — | 一致 |
   | `inf` / `nan` | missing | `inf` / `nan` | Integer missing / Float `inf`、`nan` | — | 一致 |
   | `abc` | missing | missing | missing（**不再字节序比较**） | reed 表达式曾退回字符串比较 | **对齐**（N13.2：删字节序退回） |
   | `1,2,3`（Number=1） | `0`（`len != 1` 不填 `v1`） | `0` | 见 §6.2b | htslib 的实现细节，非设计 | **决**（默认当 missing，即不对齐） |
   | `1,2,3`（Number=A） | 向量 | 向量 | 向量，比较按 any-element（N6） | reed 曾只比首元素 | **对齐**（N6.2） |
   | `2147483648` | missing（超 `BCF_MAX_BT_INT32`） | `2.14748e+09` | missing（N3） | — | **对齐**（N3） |
   | `-2147483641` | missing（低于 `BCF_MIN_BT_INT32`） | `-2.14748e+09` | missing（N3） | — | **对齐**（N3） |
   | `3000000000` | missing | `3e+09` | missing（N3） | reed 曾回绕成 `-1294967296` | **对齐**（N3） |
   | `Σac > an`（AC/AN） | `exit(1)` 报错 | 同 | 该位点丢弃（不退出） | 批处理稳健性取舍 | **不对齐**（01-review §2.5-B2） |

   表尾注明：reed 的 `-Ov`/`-Oz` 文本输出对 INFO 列**逐字节保真**（§6.7），
   不经解析，故本表只约束"判定"与"BCF 写出"两条路径；htslib 的归一化列
   （`+5`→`5`）**永远不适用于 reed 的文本输出**（01-review §2.5-B1）。

   表中 **5 行「不对齐」是永久分歧**（`5.0`/`1e3`/`0x10`/`5abc`、空格、
   Σac>an；`+5` 一行经 VCFv4.5 对照从"不对齐"转为 N13 顺手对齐——规范
   Float 正则显式允许 `[-+]?`，tex:97），落地时必须在 `cli.md §5` 原样登记
   （含"reed 更严格/更稳健"的理由），否则下一轮复审会把它们重新判成缺陷。

4. 【决】**完全对齐 htslib 前缀解析**（`hts_str2int`/`hts_str2dbl` 语义）：
   默认**不做**——会改变 gold 之外输入的过滤结果、热路径变慢（逐字符扫描
   vs SIMD 友好的严格解析）、gold 数据无此形态。若用户裁决对齐，实现为
   `parse_info_number` 的第二种模式（`ParseMode::Htslib`），规格另立，
   并把上表 reed 列改为与 htslib 同。

**验收**：
- 收口守卫：生产代码中 INFO 值的数值解析只经 `parse_info_number`
  （`grep` 断言 `native_encode.rs`/`expr/mod.rs` 无本地 `parse::<f64>`/
  `parse::<f32>`/`parse_i64` 用于 INFO 值）。
- **内部一致性红测试** `info_value_is_the_same_for_filtering_and_writing`：
  夹具含 `XI= 5`、`XI=1e3`、`XI=5abc`（声明 Integer）——断言"表达式判定为真"
  的记录，其 `-Ob` 写出值经 bcftools 读回后满足同一比较。修复前这三条
  判定为真而写出为 `.`，即红。
- 对拍测试：01-review D12-14 的 31 形态夹具，`-Ob` 与 `-i 'XI>4'` 两条路径
  逐行断言等于上表的 reed 列。
- gold 输出逐字节不变；性能门（06 §1.5）。

---

## §14 死代码：`edit_vcf_line` 的重解析（N14 ← D12-15）【规】

删除 `line_edit.rs:86-94` 的重解析与不可达错误分支。签名收为不可失败
（`-> ()`），三个调用点（`single_pass.rs:378`、`emit.rs:165`、`filtered.rs:177`）
同步去掉 `?`/`map_err`；`IoError::MalformedRecord` 的该产生点从文档注释
（`:77-79`）删除。若保守保留 `Result`，则至少删分支并加
`debug_assert!(LineRef::parse(line.line).is_some())` 注记不变式。

**验收**：`grep -c 'LineRef::parse(line.line)' crates/reed-io/src/vcf/line_edit.rs`
= 0；line_edit 与三个调用方测试全绿；clippy 无新警告。

---

## §15 文档回改清单（随各 N 落地）

| 文档 | 回改 | 关联 |
|---|---|---|
| `docs/spec/cli.md` §3.2（`-q/-Q` 行） | 补 cohort 门语义（子集/掩码后按 GT）+ 警告已兑现 | N5 |
| `docs/spec/cli.md` §3.5（`--fill-tags` 行） | 支持集收敛为 AC/AN/NS/MAF/HWE；AF/F_MISSING 标【待实现】 | N10 |
| `docs/spec/cli.md` §3.5（`-x` 行） | 默认清单 + `-x ''` 转义口 | N11 |
| `docs/spec/cli.md` §5 | Flag 比较、向量 any、裸标签处置、**解析差异表（02 §13.3，含处置列）**、Σac>an 差异 | N6/N13/N12 |
| `docs/spec/cli.md` §5（新增小节「与 htslib 的有意分歧」） | **登记 01-review §2.5-B 的 B1–B5 五条"不改"项**：文本逐字节保真、Σac>an 丢位点不退出、Integer 严格解析（不前缀）、不 trim、Number=1 多值当 missing。每条一句话理由 + 指向 §2.5-B。这是"不改"的**唯一交付物**（守卫 15），不登记则下一轮复审会重判为缺陷 | §2.5-B |
| `docs/spec/cli.md` §3.2（`-m/-M` 行） | `-M 1` 语义（只留不变位点）；`ALT='.'` 计 0 | N9 |
| `README.md`（INFO handling 段） | 默认 `-x` 清单一句话 | N11 |
| `README.zh-CN.md` | 同步 | N11 |
| `CHANGELOG.md` | Unreleased/Changed：`-q` cohort 门、默认 `-x` 可见化、`--fill-tags` 校验（行为变化均显著条目） | N5/N10/N11 |
| `docs/spec/architecture.md` §11.4 | sketch 策略触发条件补"无尾部表达式/行字段依赖"前置 | N1 |
