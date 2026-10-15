# 01 INFO 字段兼容性复审（bcftools/htslib 差分）

> 复审对象：commit `065175c`（doc_v10 全部落地；doc_v11 为门面 API 轮，与本轮
> 主题正交），代码工作区干净。
> 方法：INFO 链路全量通读（reed-io 解析/编辑/BCF 编码 → reed-filter 规则编译与
> 求值 → reed 驱动与策略 → stats 消费点，codegraph 索引辅助定位调用面）+
> 逐条 CLI 夹具实跑 + 下游探针 crate 实跑；真值门为工作区 `bcftools 1.24-2-g270751fc`
> （htslib `1.24-1-gb499e4c7`）与 `../../htslib`（vcf.c、vcfutils.c）、
> `../../bcftools`（vcfview.c、filter.c、plugins/fill-tags.c、doc/bcftools.txt）源码对照。
> 本轮不修改任何代码；修法以规范形式给出（02），实施顺序见 06。
> 证据通道：**【实】**CLI/探针实跑复现；**【据】**reed 或 htslib/bcftools 源码行号。
> 复现夹具存于 `/tmp/inforev/`（复审期），正文内联可重建。

---

## 0. 总评

v10 修的是 BCF/BGZF 的**结构**兼容（描述符、位宽、块长、残片）；本轮把同一把
尺子量到 INFO 字段的**语义**层：值编码、编辑/重算、`-q`/表达式读取。结果呈
清晰的两半：

- **读侧（BCF → 文本）已与 bcftools 逐字节一致**。差分夹具（int8/16/32 边界、
  float `-0`/`1e-05`、String、Flag、`.` 缺失、向量截断共 5 记录）：reed 渲染
  bcftools 写出的 BCF，输出与 bcftools 自渲染 `diff` 为空【实】。v3–v10 的
  `render.rs` 主路径经受住了这轮对照。
- **写侧（文本 → BCF）与过滤语义成族偏差**。共同模式与 v10 相同——**reed 在
  "自己生成的规整数据"上全对，在"别人生成的边界数据"上出错**：超界整数、
  Flag 带值、String `.`、超多 INFO/FORMAT 键；外加三类**语义**分歧：`-q` 的
  INFO 优先级不看 cohort 变化、表达式层拿不到 header 类型、无 GT 记录被重算
  清零。

本轮 15 条缺陷里 12 条能归到 10 个**共性代码根因簇**（§0.5）——同一规则的
多份实现各自漂移（AC/AN 读取 3 份、数值解析 3 份、ALT 计数 3 份、编辑语义
文本/BCF 2 份、打包计数 4 处），或"元数据存在但消费者没查"（`has_tail_expr`、
`has_gt`、`number_a`、cohort 变化）。这是 v10 D10-1（GT 位宽两份拷贝一处错）
的放大版，修法因此以"收口到单一实现"为主，而非逐点打补丁。

分档：

| 档 | 缺陷 | 后果形态 |
|---|---|---|
| P0 | D12-1 sketch 策略下尾部表达式丢行字段 | **静默错误结果**（`--mind` + 含 INFO/QUAL 的 `-i/-e` 丢光记录；`-q` 算错） |
| P0 | D12-2 无 GT 记录的 AC/AN 被子集更新清零 | **静默数据损坏**（输入 `AC=4;AN=6` 输出 `AC=0;AN=0`） |
| P1 | D12-3 Integer INFO 超 int32 回绕 | **静默数据损坏**（`-Ob` 写出 `3000000000 → -1294967296`） |
| P1 | D12-4 Flag 带值丢值、`DB=0` 语义反转 | 静默损坏（`-Ob` 后"清除"变"设置"） |
| P1 | D12-5 String 值 `.` 变空串 | 静默丢值（`-Ov`/`-Ob` 不等价） |
| P1 | D12-6 `-q/-Q` 在 cohort 变化后仍用陈旧 INFO | 与 bcftools 结果相反 + 输出自相矛盾；承诺的警告缺失 |
| P1 | D12-7 表达式无类型信息：Flag 比较恒假、向量只比首元素、裸标签真值相反 | 与 bcftools 结果相反（三类） |
| P1 | D12-8 `quantifier = "any"` 不看比较结果 | 库 API 规则恒真 |
| P2 | D12-9 打包计数（n_info/n_allele/n_fmt）溢出无检查 | 静默写出坏 BCF（双方都读不回） |
| P2 | D12-10 `ALT='.'` 计为 1 个 ALT | `-m/-M/N_ALT` 与 bcftools 相反 |
| P2 | D12-11 `--fill-tags` 不校验标签 | 静默无操作 + 文档承诺落空（AF 反被默认 `-x` 删除） |
| P2 | D12-12 默认 `drop_stale` 不可见 | 裸跑静默删 8 个标签；`-I` 单用丢 AC/AN |
| P2 | D12-13 `allele_counts_info` 不校验 AC 长度 | `-q` 与 bcftools 相反（畸形 AC 输入） |
| P2 | D12-14 INFO 数值解析三份实现、宽严不一 | **reed 内部自相矛盾**（同一记录：过滤判定值 ≠ `-Ob` 写出值 ≠ `-Ov` 写出值）+ 与 htslib 23 处差异 |
| P3 | D12-15 `edit_vcf_line` 死代码重解析 | lean |

已核实**一致**的项（不计入缺陷，作为回归基线）：

| 面 | 证据 |
|---|---|
| BCF→文本 INFO 渲染 | 差分夹具逐字节一致【实】 |
| INFO 浮点 `%g` 渲染 | gold 80 228 个 MAF/HWE 值零失配（`line_edit.rs:258-260` 既有锁定） |
| `-q/-Q` 五种 `:type` | nref/alt1/minor/major/nonmajor 与 bcftools 逐一相同【实】 |
| 子集后 AC/AN 自动更新（有 GT） | `-s s2,s3` → 双方 `AC=0;AN=4`；`-I` → 双方保留输入值【实】 |
| `--fill-tags AC,AN,NS,MAF,HWE` 值口径 | 有 GT 时与 bcftools `+fill-tags` 定义一致（gold 门已锁） |
| v10 全部修复项 | PASS 注入、多成员 gzip、NULL 描述符、块长、GT 位宽——本轮未回归 |

---

## 0.5 共性代码根因簇

> 用户复审问题"这些问题有些是因为共性代码的原因么"的正式回答。
> 15 条缺陷中 12 条落在 10 个簇上；每簇的修法是**收口**（单一实现/把元数据
> 传给消费者），见 02 的对应 N 项。

| 簇 | 共性代码 | 症状 | 规范 |
|---|---|---|---|
| α | 策略在两处决定：`plan/mod.rs:740` 有 `has_tail_expr` 守卫，驱动 `filter/mod.rs:232-265` 的覆盖不查它 | D12-1 | N1 |
| β | `has_gt` 一路传到 `SiteCounts`（`record_pipeline.rs:434`），`-g` 规则查它（`tail.rs:68`），重算不查（`edit_plan.rs:21-46`、`filter/mod.rs:367`） | D12-2 | N2 |
| γ | "值 × 声明类型 → BCF 字节"在 `encode_info_text`（`native_encode.rs:342-375`）与 `encode_info_typed`（`:994-1016`）各写一份，且都没有值域/忠实性规则 | D12-3/4/5/9 | N3/N4/N8 |
| δ | `-q` 的"INFO 优先"在编译期按 header 声明定死一次（`plan/mod.rs:382-413`），运行期不知道 cohort 变过；AC/AN 读取另有三份实现 | D12-6/13 | N5/N12 |
| ε | 表达式层只拿到标签**名字**（`check(expr, info_tags: &[String], …)`，`expr/mod.rs:215`），拿不到 Type/Number | D12-7 | N6 |
| ζ | `eval_info` 的 any/all 求值（`head.rs:209-225`） | D12-8 | N7 |
| η | ALT 计数三份：`site.rs:129-141`（处理 `.`）、`head.rs:82-91`、`expr/mod.rs:109`（都不处理） | D12-10 | N9 |
| θ | 规范/文档承诺与实现之间无校验门（`cli.md:78/:117` vs `recompute.rs:64-66`、`translate/filter.rs:317-323`） | D12-11/12 | N10/N11 |
| ι | "INFO 文本 → 数值"三份解析：`native_encode.rs:1239`（严格 i64，不 trim）、`expr/mod.rs:401`（`str::parse::<f64>`，失败退回字节序比较，且 `:156` 会 trim）、`genotype_scan.rs:421`（饱和 u32），另有 6 处裸 `str::parse`。三者互不相同 → **reed 对同一条记录给出三个值** | D12-14 | N13 |
| κ | `edit_vcf_line` 对已解析的 `LineRef` 再解析（`line_edit.rs:86-94`） | D12-15 | N14 |

---

## 1. 缺陷清单

### D12-1 sketch 策略把尾部表达式放在"无记录字节"的上下文里求值【实】【据】

- **严重度：P0**。后果形态：静默错误结果（正常参数组合丢光记录），且违反
  §11.4 的三策略等价契约。
- **reed 现状**：
  - 计划编译期有守卫：有尾部表达式就选 SinglePass（`plan/mod.rs:740-747`，
    `has_tail_expr`），注释明言"sketch path evaluates tail rules in the mid
    stage, where the record bytes are gone"。
  - 但驱动在 `filter/mod.rs:232-253` 重新选策略：`--mind` 时按内存预算走
    `choose_for_mind`（可能得 SketchThenEmit），**不查 `plan.tail`**；
    `:258-263` 对 `--mind + --indep-pairwise` 更是无条件强制 sketch。
  - mid 阶段的求值上下文是 `SiteValues::from_counts`（`mid_stage.rs:326-331`），
    其行字段全是哑值：`info: b"."`、`qual: None`、`chrom: b""`、`pos: 0`
    （`expr/mod.rs:126-142`）。尾部表达式里任何 INFO/QUAL/行字段读取都拿到
    Missing → 比较恒假。
  - `-q/-Q` 的尾部规则同伤：`af_of` 读 `values.info`（`tail.rs:116`），sketch
    下 `info = b"."` → INFO 路径永远不可用，退回 GT 计数——与 single/two_pass
    的行为不同。
- **实测**（夹具 `mind.vcf`：4 记录，QD=200/200/5/5）：
  ```
  reed filter -i 'F_MISSING<=0.5 && QD>10' -o …                    → 2 条（single_pass）
  reed filter -i 'F_MISSING<=0.5 && QD>10' --mind 0.9 -o …         → 0 条（sketch_then_emit，--report 证实）
  reed filter -i 'F_MISSING<=0.5 && QD>10' --mind 0.9 --strategy single_pass → 2 条
  reed filter -i 'F_MISSING<=0.5 && QUAL>10' --mind 0.9            → 0 条（QUAL 同样丢）
  reed filter --mind 0.9 --indep-pairwise 50 5 0.2 -i '…QD>10'     → 0 条（对照组无表达式 → 1 条）
  reed filter -q 0.1 --mind 0.9（stale.vcf，INFO AC=4/AN=8）        → 0 条（无 --mind → 2 条）
  ```
  既有测试没抓到：`all_three_strategies_produce_identical_bytes`
  （`e2e_filter.rs:403-441`）的 gold_spec 尾部只有 `max_missing`/`min_maf`，
  没有 `TailRule::Expr`。
- **htslib/bcftools 对照**：不适用（reed 特有的执行策略问题）；契约是本仓
  §11.4"三种策略对同一输入必须产出逐字节相同的输出"。
- **修法**：驱动的策略覆盖必须尊重尾部表达式——`--mind` 且 `plan.tail` 含
  `Expr` 时退回 TwoPassReparse（实测正确）；`--mind + --indep-pairwise +
  尾部表达式`无正确执行路径，编译期报错。见 02 §1。
- **验收**：红测试 `mind_with_tail_expression_keeps_records`（修复前 0 条即红）；
  三策略等价测试的夹具扩一条含 INFO 读取的尾部表达式；`--mind --indep-pairwise
  -i 'F_MISSING<=0.5 && QD>10'` 报 usage 错而非静默 0 条。

### D12-2 无 GT 的记录被子集自动更新写成 `AC=0;AN=0`【实】【据】

- **严重度：P0**。后果形态：静默数据损坏（把输入里正确的 INFO 值覆盖成 0）。
- **reed 现状**：
  - `-s` 子集触发自动更新：`filter/mod.rs:363-374`（`subset && !no_update` →
    `recompute.ac/an = true` + AC/AN 进 stale）。
  - 重算值来自 `fill_info_set`（`eval/edit_plan.rs:21-46`），直接写
    `counts.ac/an`——**不看 `counts.has_gt`**。记录无 GT 列时计数全 0。
  - `has_gt` 字段存在且已传播（`site_counts.rs:47`、`record_pipeline.rs:434`、
    `mid_stage.rs:324`），`-g` 规则正确地查了它（`tail.rs:67-85`，注释引用
    bcftools 的 `bcf_get_fmt` null check）——重算是唯一不查的消费者。
- **实测**（夹具 `nogt.vcf`：`AC=4;AN=6`，FORMAT 只有 DP，无 GT，3 样本）：
  ```
  bcftools view -s s0,s1        → AC=4;AN=6（保留）
  reed filter -s s0,s1          → AC=0;AN=0（清零）
  reed filter -s s0,s1 --fill-tags AC,AN → AC=0;AN=0
  reed filter -s s0,s1 -I -x '' → AC=4;AN=6（-I 时正确，证明缺口就在自动更新）
  ```
- **htslib/bcftools 对照**：`vcfview.c:343`
  `if ( args->calc_ac && !bcf_get_fmt(args->hdr,line,"GT") ) update_ac = 0;`
  ——无 GT 不更新；`plugins/fill-tags.c:763` `if ( !fmt_gt ) return;`——无 GT
  不算。reed 自己的 `-g` 实现（`tail.rs:62-66` 注释）已引用同一条规则。
- **修法**：`has_gt == false` 的记录跳过 AC/AN 的"删旧 + 写新"整组操作，
  保留输入值（逐记录门，bcftools 同款）。`fill_info_set` 与 emit 侧的 stale
  过滤都要接 `has_gt`。见 02 §2。
- **验收**：红测试 `gt_less_records_keep_their_ac_an_under_subset`（修复前
  `AC=0;AN=0` 即红）；有 GT 记录的更新行为不变（现有 e2e 全绿）。

### D12-3 `-Ob` 把超出 int32 的 Integer INFO 值静默回绕【实】【据】

- **严重度：P1**。后果形态：静默数据损坏（写出的 BCF 数值错误）。
- **reed 现状**：
  - `parse_int_vector`（`native_encode.rs:1115-1123`）接受任意 i64；
  - `encode_int_vector`（`:1082-1102`）按值域选宽后 `v as i32`（`:1098`）——
    超界值截断不报错；`int_width`（`:1053-1061`）对超 INT32 值域的值也只能
    选 INT32。
  - `[BCF_MIN_BT_INT32, BCF_MAX_BT_INT32]`（`record.rs:59-61`，htslib 保留带
    之外的合法值域）没有在任何写入点被检查——`-2147483641` 这类保留带内的值
    也原样写出。
- **实测**（夹具 `big.vcf`：`XD=3000000000` / `-3000000000` / `2147483647`）：
  ```
  bcftools view -Ob → 读回 XD=.  XD=.  XD=2147483647（伴随一次性警告
                      "Extreme INFO/XD value encountered and set to missing"）
  reed filter -Ob   → 读回 XD=-1294967296  XD=1294967296  XD=2147483647
  reed filter -Ov/-Oz → 逐字节保留（文本路径不解析，无此问题）
  ```
- **htslib 对照**：`vcf.c:4004-4013`——解析期
  `overflow || tmp_val<BCF_MIN_BT_INT32 || tmp_val>BCF_MAX_BT_INT32` →
  `bcf_int32_missing` + 一次性警告。
- **修法**：`parse_int_vector` 对 `[INT32_MIN_VALUE, INT32_MAX_VALUE]`
  （`record.rs` 已有常量）之外的值写 missing 哨兵（`i64::MIN`，编码端已有
  处理），一次性 `tracing::warn!`。见 02 §3。
- **验收**：红测试 `out_of_range_integers_become_missing_in_bcf`（修复前回绕值
  即红）；**htslib 真值门**：`XD=3000000000` 经 reed `-Ob` 后 bcftools 读回
  `XD=.`，与 bcftools 自产 BCF 一致。

### D12-4 声明为 Flag 的标签带值时丢值；`DB=0` 语义反转【实】【据】

- **严重度：P1**。后果形态：静默数据损坏 + 语义反转（"清除"变"设置"）。
- **reed 现状**：`encode_info_text` 的 `(InfoType::Flag, _) => encode_size(s, 0,
  value_type::NULL)`（`native_encode.rs:346`）无条件丢掉值；
  `encode_info_typed`（`:996`）同样。文本路径（`edit_info`）逐字节保留，所以
  `-Ov` 正确、`-Ob` 错——两种输出不等价。
- **实测**（夹具 `flag2.vcf`：`DB=1` 与 `DB` 各一条；`flag0.vcf`：`DB=0`）：
  ```
  输入 DB=1  → bcftools -Ob 读回 DB=1（按字符串保留）；reed -Ob 读回 DB（值丢失）
  输入 DB=0  → bcftools -Ob 读回 DB=0；reed -Ob 读回 DB（= flag 被设置，语义反转）
  输入 DB=1  → reed -Ov 读回 DB=1（文本路径正确）
  ```
- **htslib 对照**：`vcf.c:3949-3951`——
  `} else if ((y>>4&0xf) == BCF_HT_FLAG || (y>>4&0xf) == BCF_HT_STR) { // if
  Flag has a value, treat it as a string` → `bcf_enc_vchar`。即 htslib 对
  "Flag 带值"按字符串编码，读回原样。
- **修法**：Flag 带值时按字符串编码（htslib 语义）；无值 Flag 维持
  `(NULL,0)`。与 D12-5 共用同一个"值 × 声明类型"编码函数（N4）。见 02 §4。
- **验收**：红测试 `flag_with_a_value_round_trips_through_bcf`（`DB=1`/`DB=0`
  两例，修复前即红）；真值门：bcftools 读回与 bcftools 自产一致。

### D12-5 String 型 INFO 值 `.` 经 `-Ob` 变成空串（读回成裸标签）【实】【据】

- **严重度：P1**。后果形态：静默丢值；`-Ov` 与 `-Ob` 输出不等价。
- **reed 现状**：`encode_info_text` 的字符串分支
  `encode_typed_string(s, if v == b"." { &[] } else { v })`
  （`native_encode.rs:362-364`）——把值 `.` 特判成零长度。渲染端
  （`render.rs:155-162`）对 `count==0` 只打键名，于是读回是裸 `SOMESTR`。
- **实测**（夹具 `t.vcf` 第 4 条 `SOMESTR=.`）：
  ```
  bcftools -Ob 读回 → SOMESTR=.（字面保留）
  reed -Ob 读回     → SOMESTR（bcftools 与 reed 读 reed 的文件都得到裸标签）
  reed -Ov          → SOMESTR=.（文本路径正确）
  ```
- **htslib 对照**：`vcf.c:3949-3951` 字符串分支 `bcf_enc_vchar(str, end-val,
  val)` 原样存字节；`.` 只有作为**整个 INFO 列**的缺省记号才有特殊含义，
  作为标签的**值**就是字面 `"."`。ID 列的 `.` → 空串特判（`native_encode.rs:263`）
  是正确的（htslib 对 ID 确实这么做），疑似被顺手抄进了 INFO 字符串分支。
- **修法**：删掉 `if v == b"."` 特判，原样存 `v`。一行。见 02 §4。
- **验收**：红测试 `string_dot_value_round_trips_through_bcf`；真值门同上。

### D12-6 `-q/-Q` 在样本集变化后仍优先用输入的陈旧 INFO/AC、AN【实】【据】

- **严重度：P1**。后果形态：与 bcftools 结果相反 + 输出自相矛盾；承诺的警告缺失。
- **reed 现状**：
  - 编译期只要 header 声明了 AC+AN 就推 head 规则 `MinAfInfo/MaxAfInfo`
    （`plan/mod.rs:382-413`），求值直接读输入 INFO（`head.rs:94-99` →
    `non_reference_af_info`，`head.rs:161-188`）；尾部 `af_of` 同样 INFO 优先
    （`tail.rs:114-127`）。**没有任何 cohort 门**——`-s` 子集、`--mind`、
    基因型掩码之后 INFO 已陈旧，判定却仍用它。
  - 更矛盾的是：输出记录里的 AC/AN 已被自动更新为子集后的值（D12-2 的同一
    机制），于是出现"保留依据 AF=0.5、输出记录写着 AF=0"的自相矛盾。
  - `cli.md:78` 承诺"（INFO/AC 陈旧）取到时打印'陈旧 INFO'警告"
    （v4/02-cli-semantics.md:103 同）；实现只有表达式路径的 `warn_stale_info`
    （`plan/mod.rs:422/:437/:661`），`-q/-Q` 路径没接。
- **实测**（夹具 `sub.vcf`：`AC=4;AN=8`，GT s0/s1=1/1、s2/s3=0/0）：
  ```
  bcftools view -s s2,s3 -q 0.1 → 两条全丢（子集 GT 重算 AF=0）
  reed filter -s s2,s3 -q 0.1   → 两条全留，且输出记录 INFO 为 AC=0;AN=4
  reed filter -q 0.1（无子集）   → 与 bcftools 一致（INFO 优先是双方语义）
  reed filter -q 0.1 …（任何调用）→ 无"陈旧 INFO"警告（-i 'AC>1' 有，-q 没有）
  ```
- **htslib/bcftools 对照**：`vcfview.c:330` 先 `bcf_calc_ac(…INFO|FMT)`，
  但 `:341-350` 一旦有样本子集就 `bcf_subset` 后**用 GT 重算**再做 `-q` 判定；
  `filter.c:1471` 同款门：`if ( !flt->sample_mask_set && bcf_calc_ac(…,
  BCF_UN_INFO)>0 )`——有样本掩码就不信 INFO。
- **修法**：spec 编译期可得 cohort 变化信号（`sample_select.is_some() ||
  mind.is_some() || 掩码规则非空 || set_gts`）；为真时不推 `MinAfInfo/
  MaxAfInfo`、`af_of` 跳过 INFO 优先；INFO 优先真正生效时打一次"陈旧 INFO"
  警告（兑现 cli.md:78）。见 02 §5。
- **验收**：红测试 `min_af_after_subset_uses_the_genotypes`（修复前保留即红）；
  无子集时 INFO 优先行为不变（`e2e_af_type.rs` 全绿）；`-q` + 声明 AC/AN 的
  header → stderr 出现一次警告。

### D12-7 表达式层拿不到 header 类型：Flag 比较恒假、向量只比首元素、裸标签真值相反【实】【据】

- **严重度：P1**。后果形态：与 bcftools 结果相反（三类，均静默）。
- **reed 现状**：
  - `check(expr, info_tags: &[String], format_tags: &[String])`
    （`expr/mod.rs:215`）只拿到名字列表，拿不到 Type/Number——无法区分 Flag，
    也无法按 `cli.md §5` 的纪律"超出子集报 not supported yet"拒掉不支持的形态。
  - **Flag**：`info_value` 对无 `=` 的片段返回 `Val::Missing`
    （`expr/mod.rs:391-393`），任何比较恒假。
  - **向量**：值只取第一个逗号前元素（`expr/mod.rs:397`
    `let first = value.split(',').next()`）。这有两种错法，方向相反：
    声明 Number≠1 时应按 **any-element**（bcftools `CMP_VECTORS`），reed 只看
    首元素 → 漏判；声明 Number=1 而记录给了多值时，bcftools 的标量 setter
    读到 **0**（htslib `bcf_unpack` 在 `len != 1` 时不填 `v1` 联合，而它已被
    初始化为 0，`vcf.c:4286-4287`），reed 仍取首元素 → 多判。实测
    `XI=5,`（声明 Number=1）配 `-i 'XI>4'`：bcftools 丢（读到 0）、reed 留
    （读到 5）；配 `-i 'XI=0'`：bcftools 留、reed 丢（见 D12-14 路径 C）。
    注意 htslib 这个"读作 0"是实现细节而非设计（`XI=.` 的 missing 与
    `XI=5,` 的 0 在 bcftools 里可区分，见 02 §6.2 的【决】）。
    **两种情形都需要 Number 才能做对**，这是 N6 的直接论据。
  - **裸标签真值**：`expr/mod.rs:286-292` 注释声称"matching bcftools'
    truthiness"，实现是"非零数为真"——bcftools 对单个裸 token 恒假。
- **实测**（夹具 `flag3.vcf`：pos10 有 `DB`、pos20 无；`vec.vcf`：
  `DP4=0,0,30,40` 与 `DP4=1,2,3,4`；`zero.vcf`：`DP=0/5/-3`）：
  ```
  -i 'DB=1'    bcftools 留 pos10；reed 留 0 条（'DB!=0'、'INFO/DB=1' 同）
  -i 'DB=0'    bcftools 留 pos20；reed 留 0 条
  -i 'DP4>10'  bcftools 留 pos10（任一元素）；reed 留 0 条（只看首元素 0）
  -i 'DP'      bcftools 留 0 条；reed 留 pos20/30（DP=5/-3 非零）
  -e 'DP'      bcftools 留 3 条；reed 留 1 条
  -i 'DP>0'    双方一致（有比较符的路径是对的）
  ```
- **htslib/bcftools 对照**：
  - Flag：`doc/bcftools.txt:4067`——"1 (or 0) to test the presence (or absence)
    of a flag: `FlagA=1 && FlagB=0`"；实现 `filter.c:964-971`
    `filters_set_info_flag`（存在→1，不存在→0，类型 INT）。
  - 向量：`filter.c:2738-2790` `CMP_VECTORS`——两侧值集合**任一对**满足即真。
  - 裸 token：`filter.c:4352-4380` `filter_test`——`TOK_VAL` 只填值，
    `pass_site` 保持初始 0；只有比较/函数节点会置 `pass_site`。
- **正确先例（本仓库）**：`InfoRule` 路径（`[[site.info]]`）在编译期就拿到
  `TagDef`（`plan/mod.rs:452-492`：类型校验、`number_a`），表达式路径应对齐
  同一信息源（`TagTable` 已存在，`tag_table.rs:57`）。
- **修法**：`check`/求值上下文接 `&TagTable`：Flag 比较按存在性 1/0；
  Number≠1 的标签比较按 any-element；裸标签真值对齐 bcftools（恒假）或编译期
  拒绝；其余未实现形态（下标、函数）维持"not supported yet"。见 02 §6。
- **验收**：三组红测试（flag 比较、向量 any、裸标签）逐一与 bcftools 对拍；
  `cli.md §5` 的支持清单回改。

### D12-8 `[[site.info]]` 的 `quantifier = "any"` 只看可解析性、不看比较结果【实】【据】

- **严重度：P1**。后果形态：库 API 的 Number=A 规则恒真（值可解析即通过）。
- **reed 现状**：`eval_info`（`head.rs:209-225`）：
  ```rust
  Some(v) => { any = true; all &= rule.cmp.apply_f64(v); }
  ```
  `any` 在元素**能解析**时置真，比较结果只喂 `all`。枚举文档
  （`predicate.rs:168-176`）与 `architecture.md:1210` 的语义都是"至少一个
  值**满足比较**"。
- **实测**（探针 crate，`AC=1,2`）：
  ```
  keep_if ">= 10", quantifier any → true（应为 false）
  keep_if "<= 0",  quantifier any → true（应为 false）
  keep_if ">= 10", quantifier all → false（正确）
  ```
- **触发面**：CLI 无 `[[site.info]]` 拼写（INFO 规则经 `-i` 表达式进入
  `site.include`；`spec/render.rs:8` 注明 typed INFO predicates 无 flag 形式），
  但 `InfoPredicate`/`Quantifier` 是公开类型，库调用方直接可踩。
- **修法**：`any |= rule.cmp.apply_f64(v)`（删 `any = true`）。一行。见 02 §7。
- **验收**：红测试 `quantifier_any_requires_a_satisfying_element`（四组合
  ≥/≤ × any/all；修复前 any 两例即红）。

### D12-9 记录打包计数（n_info/n_allele/n_fmt/n_sample）溢出无检查，写出坏 BCF【实】【据】

- **严重度：P2**。后果形态：静默写出双方都读不回的 BCF（极端输入）。
  **夹具有效性说明**（对照 VCFv4.5 后补）：首组夹具用了 65 536 个**重复键**
  `DP=1;DP=1;…`，这违反 §1.6.1（tex:434）"Duplicate keys are not allowed"，
  属规范外输入（bcftools 实测也会原样往返重复键——两边都以宽容恢复处理，
  该夹具测的是计数溢出而非重复键语义）；**缺陷的成立不依赖它**——同轮
  的 65 536 个**互异键**夹具（`XD=0;XD=1;…`，规范内：arbitrary keys
  permitted）复现同一损坏，红测试以互异键夹具为准。
- **reed 现状**：四处直接打包写回：`native_encode.rs:336`/`:793`
  （`n_allele << 16 | n_info`）、`:479`/`:851-853`（`n_fmt << 24 | n_sample`）。
  计数超限（n_info/n_allele > 0xFFFF、n_fmt > 0xFF、n_sample > 0xFFFFFF）时
  高位溢进相邻字段，无报错。
- **实测**：
  ```
  65 536 个 INFO 片段（重复键 DP）→ reed -Ob 成功写出；
    bcftools 读："Bad BCF record … Invalid REF/ALT type 1"；reed 自读：
    "BCF string field is not CHAR-typed"（65 535 个则双方正常）
  300 个 FORMAT 键 → reed -Ob 写出 n_fmt=44（300 & 0xFF），静默截断；
    bcftools 对同输入在解析期报错 "FORMAT column … lists more identifiers
    than htslib can handle"
  66 001 个等位（--max-alleles 70000）→ reed -Ob 写出；bcftools 读：
    "Invalid FILTER type"；reed 自读："BCF integer field has a non-integer type"
  ```
- **htslib 对照**：文本解析期即拒绝——`vcf.c:3912`（`n_info == UINT16_MAX`
  → "Too many INFO entries"）、`vcf.c:4174`（等位同款）、FORMAT 键上限在
  `vcf_parse_format_dict2`。v10 D10-8 已给字段**字节数**加了 checked 乘法，
  本条是同一族缺口在**计数字段**上的残留。
- **修法**：写回前统一检查四个计数上限，超限返回 `MalformedRecord`（带上限
  说明）；文本与 BCF 两条编码路径共用一个 helper。见 02 §8。
- **验收**：红测试 `info_count_overflow_errors_instead_of_corrupting`
  （65 536 片段夹具；修复前写出坏文件即红）；300 FORMAT 键、66 001 等位同款。

### D12-10 `ALT='.'` 被计为 1 个 ALT：`-m/-M/N_ALT` 与 bcftools 相反【实】【据】

- **严重度：P2**。后果形态：不变位点（monomorphic）的过滤结果相反。
- **reed 现状**：ALT 计数有三份实现，只有一份处理 `.`：
  - `SiteView::n_alt`（`site.rs:129-135`）：`alts == b"." → 0`，**正确**；
  - head 规则 `MinAltAlleles/MaxAltAlleles`（`head.rs:82-91`）：逗号数 +1，
    `.` 得 1；
  - 表达式 `n_alt`（`expr/mod.rs:109`）：逗号数 +1，`.` 得 1。
  - 另一独立缺口：`-M 1` 经 translate 变 `max_alt_alleles = 0`
    （`translate/filter.rs:112-114`），而 plan 把 0 当"关闭"
    （`plan/mod.rs:352`）——`-M 1` 整个规则消失。
- **实测**（夹具 `inv.vcf`：pos10 `ALT=.`、pos20 `ALT=G`）：
  ```
  bcftools -m 2 → pos20      reed -m 2 → pos10 pos20
  bcftools -M 1 → pos10      reed -M 1 → pos10 pos20（规则未推入）
  bcftools -i 'N_ALT=0' → pos10   reed -i 'N_ALT=0' → 0 条
  bcftools -i 'N_ALT=1' → pos20   reed -i 'N_ALT=1' → pos10 pos20
  ```
- **htslib/bcftools 对照**：`vcfview.c:418-434`（`-m/-M` 按 `line->n_allele`，
  `ALT=.` 时 `n_allele=1`）；`filter.c:1453-1456` `filters_set_nalt`
  （`n_allele - 1` → 0）。
- **修法**：ALT 计数收口到 `SiteView::n_alt` 的同一规则（`.` → 0），head 与
  表达式共用；`max_alt_alleles` 改 `Option<u32>`（或等价）区分"关闭"与"0"。
  见 02 §9。
- **验收**：红测试 `invariant_sites_count_zero_alts`（`-m 2`/`-M 1`/`N_ALT=0`
  三例与 bcftools 对拍；修复前即红）。

### D12-11 `--fill-tags` 不校验标签：未实现与不存在的都静默通过【实】【据】

- **严重度：P2**。后果形态：静默无操作 + 文档承诺落空 + 与默认 `-x` 叠加成
  "反向操作"。
- **reed 现状**：
  - 实现集只有 AC/AN/NS/MAF/HWE（`RecomputeFlags`，`plan/mod.rs:263-275`；
    `wants()` 精确大小写匹配，`recompute.rs:64-66`）；translate 直接把逗号
    列表塞进 `spec.recompute.tags`（`translate/filter.rs:317-323`），无校验。
  - `cli.md:117` 承诺"重算并写入 AC、AN、**AF**、MAF、NS、HWE、**F_MISSING**…"。
  - 默认 `drop_stale` 含 AF（D12-12）：`--fill-tags AF` 的净效果是**删除**
    输入里已有的 AF。
- **实测**：
  ```
  reed filter --fill-tags AF     → 退出码 0，无警告，输出无 AF（输入 AF=0.5 被默认 -x 删除）
  reed filter --fill-tags BOGUS  → 退出码 0，无警告，无操作
  reed filter --fill-tags ac,an  → 退出码 0，静默无操作（大小写敏感）
  bcftools +fill-tags -t BOGUS   → "the tag \"BOGUS\" is not supported"，exit 1
  ```
- **bcftools 对照**：`plugins/fill-tags.c:500-538`——白名单逐一
  `strcasecmp`，未知标签报错退出；`:541-560` `list_tags` 列出全部支持项。
- **修法**：translate/编译期校验 `tags` ⊆ 支持集（大小写不敏感归一，对齐
  bcftools），未知项 usage 错误并列出支持集；`cli.md:117` 的承诺收敛为已实现
  集，AF/F_MISSING 标注【待实现】或直接实现（AF 在双等位门下 = `ac/an`，量小；
  是否做列为【决】）。见 02 §10。
- **验收**：`--fill-tags AF`、`--fill-tags BOGUS` 退出码 2 且报错列出支持集；
  `--fill-tags ac` 与 `AC` 等价；gold 命令（大写五标签）输出不变。

### D12-12 默认 `drop_stale`（8 标签）对用户不可见；`-I` 单用把 AC/AN 删光【实】【据】

- **严重度：P2**。后果形态：裸跑静默删除用户数据；文档/help 无提示。
- **reed 现状**：
  - `RecomputeSpec::default()` 的 `drop_stale` 默认即 gold 清单
    `["AC","AN","AF","MLEAC","MLEAF","ExcessHet","InbreedingCoeff","DP"]`
    （`recompute.rs:13-19`、`:32-46`）；CLI 不给 `-x` 时**不覆盖**该默认
    （`translate/filter.rs:309-316` 仅在 `Some(list)` 时赋值）。
  - `--help` 对 `-x` 只写 "Remove these INFO tags (bcftools `-x`)"
    （`cli/mod.rs:101-102`）；README:103-105 同样未提默认清单；
    `cli.md:116` 也没有。唯一写明处是 `architecture.md:1201` 的 TOML 示例与
    `gold/scripts/run_reed.sh:15`（gold 管线**显式**传 `-x`，不依赖默认）。
  - `-I` 交互：`reed filter -I`（无子集）→ 默认 `-x` 删 AC/AN、`-I` 又禁止
    重算 → 输出 `INFO=.`；bcftools `view -I` 保留输入值。
- **实测**：
  ```
  输入 DP=30;AF=0.5;MLEAC=1;ExcessHet=0.1;InbreedingCoeff=-0.2;FS=3
  reed filter（无任何 INFO 参数）→ FS=3.0（其余 5 个被删；header 定义仍在）
  reed filter -I                → INFO=.（AC/AN 被删且不回填）
  bcftools view / view -I       → 原样保留
  reed filter -x ''             → 原样保留（转义口存在但未文档化）
  ```
- **定性**：默认清单是 gold 对拍的既定设计（v2 02 §3.5"gold `annotate -x`
  清单为默认"），**行为不改**；缺陷在可见性——`--help`、README、cli.md 都
  必须写明"默认删除以下标签，`-x ''` 关闭"。是否把默认改为 opt-in 列为【决】
  （改默认会破坏 gold 命令的现有拼写习惯与既有用户脚本预期，两个方向都有理由）。
- **修法**：文档三处回改 + `-I` 的帮助文本注明与默认 `-x` 的叠加效果。
  见 02 §11。
- **验收**：`--help` 含默认清单字样；README/cli.md 回改；drift 守卫。

### D12-13 `allele_counts_info` 不校验 AC 长度；AC/AN 读取三份实现口径不一【实】【据】

- **严重度：P2**。后果形态：畸形 INFO 输入下 `-q` 与 bcftools 相反。
- **reed 现状**：AC/AN 的 INFO 读取有三份：
  - `head.rs:128-155` `allele_counts_info`：不比对 AC 元素数与 ALT 数；
  - `head.rs:161-188` `non_reference_af_info`：只求和（求和语义下长度不敏感，
    本身无错，但与上一份重复了扫描逻辑）;
  - `summary.rs:292-320`（stats AF/SiS）：**有**长度校验
    （`:302` `vals.len() == n_alt`）与 `an >= Σac` 校验。
- **实测**（夹具 `aclen3.vcf`：`AC=3`（1 个值）+ 2 个 ALT + `AN=4`，GT 全 0/0）：
  ```
  bcftools view -M 3 -q 0.1 → 警告 "Incorrect number of AC fields" + 丢弃
                              （回退 GT：AF=0）
  reed --max-alleles 3 -q 0.1 → 保留（用 INFO 3/4=0.75）
  ```
  另：`AC=5;AN=4`（Σac > an）时 bcftools `exit(1)` 报错（`vcfutils.c:85-90`），
  reed 静默丢弃该位点——方向一致（都不保留），烈度不同，登记不升级。
- **htslib 对照**：`vcfutils.c:56-64`——`ac_len != n_allele - 1` → 一次性警告
  + `return 0`（调用方回退 GT 计数）。
- **修法**：三份收口成一个实现（放 `reed-filter` 或 `reed-core`，带 `n_alt`
  参数），长度不符 → `None`（回退 GT）；`summary.rs` 换用同一函数。见 02 §12。
- **验收**：红测试 `ac_length_mismatch_falls_back_to_genotypes`（aclen3 夹具，
  修复前保留即红）；stats 的 AF/SiS 现有测试全绿。

### D12-14 "INFO 文本 → 数值"有三份解析实现，宽严互不相同【实】【据】

- **严重度：P2**（初审记 P3，31 形态 × 3 路径的全量差分后上调）。后果形态：
  **同一条记录在 reed 内部得到三种不同的值**——过滤判定、`-Ob` 写出、`-Ov`
  写出互不一致；且其中两种都与 htslib 不同。
- **reed 现状**（三份互不相同的实现）：
  - 编码器（整数）：`parse_i64`（`native_encode.rs:1239-1258`）严格十进制，
    可选前置 `-`，不 trim，拒绝 `+5`/`5.0`/`1e3`/`0x10`/`5abc`
    （其中 **`+5` 是规范内输入**——Float 正则允许 `[-+]?`（tex:97），见
    §2.5-B3 注；拒绝它属 N13 要顺手修的，其余四项是规范外/跨类型形态）。
  - 编码器（浮点）：`str::parse::<f32>`（`:1134`）——Rust 语法，接受
    `inf`/`nan`，不 trim。
  - 表达式：`str::parse::<f64>`（`expr/mod.rs:401`），失败则退回**字节序
    比较**（`compare`，`:452-460`）；且 `Val::as_num` 对字符串值会 `trim()`
    （`:156`）——编码器不 trim，表达式 trim。
  - FORMAT 扫描：`parse_int`（`genotype_scan.rs:421-435`）u32 饱和（语义不同，
    本轮不并）。
  - 另有 `head.rs:141/:146/:174/:179`、`summary.rs:298/:310` 等 6 处裸
    `str::parse`。
- **实测**：31 种书写形态 × 3 条路径的全量差分（夹具 `/tmp/parsecmp/`，
  `XI` Integer / `XF` Float 各一，真值 bcftools 1.24）。

  **路径 A：Integer 标签经 `-Ob`**——10 处不一致：

  | 输入 | bcftools | reed | 成因 |
  |---|---|---|---|
  | `+5` | `5` | `.` | reed 不接受前置 `+` |
  | `5.0` | `5` | `.` | `hts_str2int` 前缀解析到 `.` 停止 |
  | `1e3` | `1` | `.` | 同上，停在 `e` |
  | `0x10` | `0` | `.` | 同上，停在 `x` |
  | `5abc` | `5` | `.` | 同上 |
  | `5 `（尾空格） | `5` | `.` | reed 不 trim |
  | `3000000000` | `.` | `-1294967296` | 值域缺失（D12-3） |
  | `9223372036854775807` | `.` | `-1` | 同上 |
  | `-2147483649` | `.` | `2147483647` | 同上 |
  | `-2147483641` | `.` | `-2147483641` | htslib 保留带（`BCF_MIN_BT_INT32`）内亦置 missing，reed 原样写 |

  **路径 B：Float 标签经 `-Ob`**——4 处不一致（其余 27 形态全同，含
  `inf`/`nan`/`1e400`/`99999999999999999999`，这条路径基本正确）：

  | 输入 | bcftools | reed |
  |---|---|---|
  | `0x10` | `16` | `.` |
  | ` 5`（前空格） | `5` | `.` |
  | `5 `（尾空格） | `5` | `.` |
  | `5abc` | `5` | `.` |

  **路径 C：表达式 `-i 'XI>4'`**——9 处不一致（reed 多保留 9 条）：

  | 输入 | bcftools | reed | 成因 |
  |---|---|---|---|
  | `abc` | 丢 | **留** | 解析失败 → 当字符串 → 字节序 `"abc" > "4"`（`'a'`=0x61 > `'4'`=0x34）为真 |
  | `inf` | 丢 | **留** | Rust `parse::<f64>` 接受 `inf` |
  | `1e3` | 丢 | **留** | bcftools 按 Integer 得 1；reed 按 f64 得 1000 |
  | `3000000000` | 丢 | **留** | bcftools 置 missing；reed f64 解析成功 |
  | `2147483648` | 丢 | **留** | 同上 |
  | `9223372036854775807` | 丢 | **留** | 同上 |
  | `99999999999999999999` | 丢 | **留** | 同上（f64 得 1e20） |
  | ` 5`（前空格） | 丢 | **留** | 表达式 `as_num` 会 trim，编码器不 trim |
  | `5,` | 丢 | **留** | Number=1 声明却给多值：bcftools 标量读作 **0**（htslib `len != 1` 不填 `v1`，`vcf.c:4286`），reed 取首元素 5 |

- **reed 内部自相矛盾（比"与 bcftools 不同"更严重）**：
  ```
  XI= 5   → 表达式判定 5（保留记录）；-Ob 写出 XI=.（missing）；-Ov 逐字节留 " 5"
  XI=1e3  → 表达式判定 1000（保留）；-Ob 写出 XI=.；-Ov 写出 XI=1e3
  XI=5abc → 表达式判定 5（保留）；-Ob 写出 XI=.
  ```
  即：过滤决定所依据的值，与写进输出的值不是同一个。用户无法从任何单一输出
  推断 reed 认为该值是多少。
- **htslib 对照**：`hts_str2int`/`hts_str2dbl`（前缀解析，停在第一个非法字符）
  + 值域检查（超界置 missing，`vcf.c:4004-4013`）。**关键是 bcftools 的表达式
  层始终按 header 声明的类型取值**（Integer 标签走 `filters_set_info_int` →
  `bcf_get_info_int32`，`filter.c:845-861`），所以 `abc` 只会是 missing，永远
  不会退化成字符串比较。reed 的表达式层没有类型信息（D12-7 同根因），只能
  "先试 f64、失败当字符串"，于是多出字节序比较这条路。
- **定性与修法**：完全对齐 htslib 前缀解析会改变 gold 之外输入的过滤结果、
  热路径变慢，收益存疑；**但"三条路径互相一致"是必须的**，与是否对齐 htslib
  无关。最小修法：一个 `parse_info_number`（按声明类型分整数/浮点两种模式）
  供编码器、表达式、head 规则共用，消除内部矛盾；与 htslib 的残余差异登记成
  表进 `cli.md §5`。表达式层的字符串退回比较需要类型信息，随 N6 一并落地。
  是否全对齐 htslib 前缀语义列为【决】。见 02 §13。
- **验收**：三份实现收口为一处（`grep` 守卫）；内部一致性测试——同一夹具上
  "表达式判定为真"的记录，其 `-Ob` 写出值必须满足同一比较（修复前 ` 5`/`1e3`/
  `5abc` 三例即红）；差异表进 cli.md；gold 输出逐字节不变（gold 数据是规整
  十进制，各规则同结果）；性能门（06 §1.5）。

### D12-15 `edit_vcf_line` 对已解析的 `LineRef` 再解析，错误分支不可达【据】

- **严重度：P3（lean）**。
- **reed 现状**：`line_edit.rs:86-94` 重新 `LineRef::parse(line.line)` 并在
  `None` 时报 `MalformedRecord`。但 `LineRef` 字段私有（`line.rs:13`）、唯一
  构造入口是 `parse`（`line.rs:20-32`，成功即 ≥9 个 tab），所以入参
  `line: &LineRef` 恒满足条件——重解析恒 `Some`，错误分支不可达，且每条存活
  记录多付一次 9-tab memchr 扫描（热路径）。
- **修法**：删除重解析块；`edit_vcf_line` 随之不可失败，签名改 `-> ()`
  （三个调用点 `single_pass.rs`/`emit.rs`/`filtered.rs` 同步），或保守保留
  `Result` 仅删死分支。见 02 §14。
- **验收**：`grep -c "LineRef::parse(line.line)" crates/reed-io/src/vcf/line_edit.rs`
  为 0；调用点编译通过；line_edit 既有测试全绿。

---

## 2. 已核减 / 需说明的说法

> 与"是否对齐 htslib/bcftools"有关的核减，**判定与理由统一在 §2.5**
> （那里有判定规则、A/B/C 三类清单和统计）。本表只保留非对齐性质的核减，
> 以及指向 §2.5 对应条目的交叉引用，避免两处判定漂移。

| 说法 | 处置 | 依据 |
|---|---|---|
| "`-Ov` 输出 `XD=+5` 与 bcftools 的 `XD=5` 不一致是缺陷" | **核减 → §2.5-B1（标注不改）**：文本路径逐字节保真是 reed 的既定设计；bcftools 归一化是 record 模型往返的副作用，不是功能 | architecture.md §6.7；本轮实测 |
| "sites-only VCF（8 列）被拒是 INFO 缺陷" | **移出本轮 → §2.5-B8（标注不改）**：`LineRef` 的 9-tab 契约（`line.rs:18-31`）是格式级决定，牵连扫描器/编辑器/编码器/index 全链。独立立项候选（06 §4） | line.rs 注释；06 §4 |
| "表达式应实现 bcftools 的向量下标/函数/算术" | **不做 → §2.5-B7（标注不改）**：维持 `cli.md §5` 的子集纪律；D12-7 只收 flag 比较、向量 any、裸标签三类**已接受语法**的语义，其余照旧报 "not supported yet" | cli.md §5 |
| "`--fill-tags` 应实现 AF/F_MISSING" | **【决】→ §2.5-C3**：N10 的最小修法是校验 + 文档收敛；AF（双等位门下 = ac/an）实现量小但属新功能 | 02 §10 |
| "默认 `drop_stale` 应改 opt-in" | **【决】→ §2.5-C4**：默认清单是 gold 对拍设计（v2 02 §3.5）；N11 最小修法是可见性 | 02 §11 |
| "AC=5;AN=4（Σac>an）应像 bcftools 一样 exit(1)" | **登记不升级 → §2.5-B2（标注不改）**：reed 丢位点继续跑、bcftools 进程退出；两者都不保留该位点，reed 对批处理更稳 | vcfutils.c:85-90；本轮实测 |
| "前缀解析（`+5`/`1e3`/`0x10`/`5abc`）应对齐 htslib" | **核减 → §2.5-B3/B4（标注不改）**：htslib 的前缀解析自己不自洽（`0x10` 在 Integer 得 0、Float 得 16；`5abc` 静默成 5）；reed 三条路径**必须统一**（§2.5-A17），但统一目标是"严格"而非 htslib 的"前缀"。**例外**：`+5` 是规范内输入（Float 正则 `[-+]?`，tex:97），N13 整数模式顺手接受，不属"不对齐"的坚守 | 本轮实测；hts_str2int/hts_str2dbl；VCFv4.5 §1.3 |
| "`-q/-Q` 的 f32 量化是精度缺陷，应改 f64" | **核减 → §2.5-B6（护栏项）**：这是**刻意对齐**——`af_threshold`（`plan/mod.rs:796-798`）复现 bcftools 的 `sscanf("%f")`，改 f64 会在 `AC=34,AN=680` vs `0.05` 边界翻转结果、破坏 gold 对拍 | plan/mod.rs:796-798 注释 |
| "读侧（BCF→文本）也有问题" | **核减**：差分夹具逐字节一致；`render.rs` 的 INFO 路径（`fmt_array`/`fmt_array1`/flag `count==0`）与 htslib `bcf_fmt_array`（vcf.c:3111-3160）、INFO 单值分支（vcf.c:4524-4552）逐分支对照无偏差 | 本轮差分实测 |

---

## 2.5 对齐判定：哪些必须遵从 htslib/bcftools，哪些标注为**不改**

> 本节回答"这些差异里哪些是 reed 更准"。判定规则（按优先级）：
>
> 0. **VCF 规范文本有明文**（本轮参照用户提供的 VCFv4.5 草案
>    `01-VCFv4.5.tex`；4.5 新增条款注明，4.2 已有条款直接适用）→
>    **规范怎么说就怎么办**，优先于任何工具的实现。
> 1. **htslib 有明确设计依据**（源码注释、显式值域检查、两处以上同款门）→
>    **必须对齐**，reed 现状是缺陷。
> 2. **bcftools 文档明文承诺**，而 `cli.md` 又声称"语义同 bcftools 手册"→
>    **必须对齐**（这是对用户的承诺，不是偏好）。
> 3. **reed 违反自己的文档/契约** → 与 bcftools 无关，**必须修**。
> 4. **htslib 的行为是实现细节而非设计**（偶然值、为 C 语言便利）→
>    reed 可自选更自洽的语义，**标注不改**，只登记差异。
> 5. **reed 的行为严格更强**（保真、稳健、拒绝畸形）且不影响互操作 →
>    **标注不改**。
>
> "标注不改"的含义：**永久分歧**，必须写进 `cli.md §5` 差异表，
> **不写代码、不写红测试断言与 bcftools 相同**；下一轮复审不得重新判为缺陷
> （06 §1 通用门据此免除这些项的对拍断言）。
>
> **规范文本核验结果**（4.5 草案，本轮新引入的量尺）：§1.3 Data types
> 明文 Integer 为"32-bit, signed"且 **`[−2³¹, −2³¹+7]` 在 VCF 和 BCF 里都
> 被禁止**（tex:101）——D12-3 因此从"对齐 htslib"升格为"规范明文"；
> Float 有显式正则 `^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$`（允许前置
> `+`、**不含空白**，tex:97）——`+5` 在 Float 是**规范内**输入（reed 已与
> bcftools 一致），` 5` 在数值标签**不匹配任何规范文法**；§1.4.2 明文
> "Flag 类型的 INFO 字段不含 Value 项"（tex:172）——`DB=1` 是**规范外**
> 输入，D12-4 属"畸形输入的恢复策略"对齐；§1.6.1 明文 "Space characters
> are allowed in values"（tex:431，**4.5 新增**）与 "Duplicate keys are
> not allowed"（tex:434）。

### A. 必须对齐（reed 现状是缺陷）

| # | 分歧 | 判定档 | 依据 | 规范 |
|---|---|---|---|---|
| A1 | 超 int32 整数 `-Ob` 回绕 | **0** + 1 | **VCFv4.5 §1.3（tex:101）明文**：Integer 为 32-bit signed，`[−2³¹, −2³¹+7]`（即 `−2147483648…−2147483641`）"cannot be stored in the binary version and therefore are **disallowed in both VCF and BCF**"——reed 写出 `-2147483641` 原样与回绕 `3000000000→−1294967296` 都是写规范禁止的值；htslib 的 missing 处置（`vcf.c:4004-4013`）与规范一致 | N3 |
| A2 | Flag 带值被丢、`DB=0` 语义反转 | 1（恢复策略） | 输入本身**规范外**：§1.4.2（tex:172）"The 'Flag' type indicates that the INFO field **does not contain a Value entry**"。htslib 的处置（`vcf.c:3949-3951` "if Flag has a value, treat it as a string"）是对规范外输入的恢复策略；reed 丢值/反转连恢复都算不上。对齐 htslib 的恢复方式 | N4 |
| A3 | String 值 `.` 变空串 | 1 | `bcf_enc_vchar(str, end-val, val)` 字面存；`.` 作为**单个值**是字面字符串，只有作为**整个 INFO 列**的缺省记号才有特殊含义（§1.6.1 tex:428）。reed 的 `.` 特判是从 ID 列错误外推（ID 的 `.`→空串 htslib 确实做，`native_encode.rs:263` 保留） | N4 |
| A4 | `-q/-Q` 在 cohort 变化后仍用陈旧 INFO | 1（两处同款门） | `vcfview.c:341-350`（`bcf_subset` 后重算）+ `filter.c:1471`（`!flt->sample_mask_set` 才信 INFO） | N5 |
| A5 | 无 GT 记录的 AC/AN 被清零 | 1（两处同款门） | `vcfview.c:343` `!bcf_get_fmt(…,"GT") → update_ac = 0`；`fill-tags.c:763` `if (!fmt_gt) return;`。**reed 自己的 `-g` 已引用同一规则**（`tail.rs:62-66` 注释） | N2 |
| A6 | `ALT='.'` 计为 1 个 ALT | 1 | `filter.c:1453-1456` `filters_set_nalt` 返回 `n_allele - 1`；`vcfview.c:418-434` 的 `-m/-M` 按 `n_allele` | N9 |
| A7 | 打包计数溢出写出坏 BCF | 1 | `vcf.c:3912`（n_info）、`:4174`（n_allele）解析期即拒绝；reed 是写侧，等价位置在写回前 | N8 |
| A8 | Flag 比较 `DB=1`/`DB=0` 恒假 | 2 | `doc/bcftools.txt:4067` "1 (or 0) to test the presence (or absence) of a flag"；实现 `filter.c:964-971` | N6.1 |
| A9 | 向量 INFO 只比首元素 | 2 | `filter.c:2738-2790` `CMP_VECTORS` 双层循环：任一对满足即真 | N6.2 |
| A10 | 表达式对非数值退回字节序比较（`abc > 4` 为真） | 2 + 3 | bcftools 的表达式层**始终按 header 类型取值**（`filter.c:845-861` Integer 标签走 `bcf_get_info_int32`），不存在字符串退回；且 reed 自己的 `cli.md §5` 承诺"与缺失值比较为假" | N13.2 |
| A11 | `quantifier = "any"` 不看比较结果 | 3 | `predicate.rs:168-176` 与 `architecture.md:1210` 都写"至少一个值**满足比较**"——实现向自己的文档收敛，与 bcftools 无关 | N7 |
| A12 | sketch 策略丢尾部表达式记录 | 3 | `architecture.md:1762` 三策略 sha256 相同；`plan/mod.rs:740` 的 `has_tail_expr` 守卫本来就为此存在，只是驱动覆盖时没查 | N1 |
| A13 | 裸标签真值与注释相反 | 3 | `expr/mod.rs:286-287` 注释自称 "matching bcftools' truthiness"；bcftools 单 token 恒假（`filter.c:4352-4380` `pass_site` 无比较节点不置位） | N6.3（【决】：对齐恒假 或 编译期拒绝，默认后者） |
| A14 | `--fill-tags` 不校验标签 | 3 | `cli.md:117` 承诺"定义同 `bcftools +fill-tags`"；bcftools 对未知标签 `exit 1`（`fill-tags.c:530-535`） | N10 |
| A15 | 默认 `-x` 清单在 help/README/cli.md 无提示 | 3 | `cli.md:116` 只写"删除 INFO 标签"，未提默认值；用户无法从任何文档得知裸跑会删 8 个标签 | N11 |
| A16 | `-q/-Q` 缺"陈旧 INFO"警告 | 3 | `cli.md:78` 与 v4/02:103 明文承诺；实现只接在表达式路径（`plan/mod.rs:422/:437/:661`） | N5.4 |
| A17 | **reed 内部自相矛盾**（判定值 ≠ `-Ob` 值 ≠ `-Ov` 值） | 3 | 不是对齐问题，是正确性问题：`XI=1e3` 过滤按 1000 判、`-Ob` 写 `.`、`-Ov` 留 `1e3`。用户无法从任何单一输出推断 reed 认为该值是多少 | N13.1 |
| A18 | AC 长度不符仍按 INFO 判 `-q` | 1 | `vcfutils.c:56-64` `ac_len != n_allele - 1` → 警告 + `return 0`（回退 GT）。**reed 自己的第三份实现已有此校验**（`summary.rs:302`） | N12 |

### B. 标注**不改**（reed 更优，永久分歧，只登记）

| # | 分歧 | 判定档 | reed 现状 | bcftools/htslib | 不改的理由 |
|---|---|---|---|---|---|
| **B1** | `-Ov`/`-Oz` 文本输出的 INFO 列归一化 | 5 | **逐字节保真**（`XD=+5` 原样留） | 归一化（`XD=5`），因为全量走 record 模型往返 | reed 的 §6.7 设计就是"columns 1–7 copy verbatim"。保真**严格强于**归一化：用户写什么就得到什么。bcftools 的归一化是往返的副作用，不是功能 |
| **B2** | `Σac > an`（AC/AN 自相矛盾）的处置 | 5 | 丢该位点，继续处理后续记录 | **`exit(1)` 整个进程退出**（`vcfutils.c:85-90`） | 两者都不保留该位点，但 reed 不因一条畸形记录毁掉整轮批处理。对 3.3M 位点 / 数小时的运行，进程退出是更差的失败模式 |
| **B3** | Integer 的前缀解析（`5.0`/`1e3`/`0x10`/`5abc`；`+5` 见 B3 注） | 4 + 5 | 严格十进制，畸形 → missing | `hts_str2int` 前缀解析：`5abc`→5、`1e3`→1、`0x10`→0 | htslib 的前缀解析**自己也不自洽**：`0x10` 在 Integer 得 0、在 Float 得 16（`hts_str2dbl` 认十六进制）；`5abc` 静默接受成 5 而不告警。`5.0`/`1e3` 在 Integer 标签**不匹配规范的任何数值文法**（§1.3：Integer 是 32-bit signed，Float 正则是另一类型）。reed 拒绝畸形输入并置 missing 更可预测。**注意**：reed 三条路径必须统一（A17），但统一到的目标是"严格"，不是 htslib 的"前缀"。**B3 注（`+5`）**：规范的 Float 正则**显式允许前置 `+`**（tex:97 `[-+]?`），Float 路径 reed 已与 bcftools 一致；Integer 未给文法，但与 Float 对称读，`+5` 应视为合法——N13 的整数模式应**接受可选前置 `+`**（一行改动，消掉路径 A 的第一行分歧），这不属于"不对齐"的坚守范围 |
| **B4** | 空格容忍（` 5` / `5 `，数值标签） | 0 + 4 | 三条路径统一为**不 trim** | Integer 不容忍、Float 容忍前空格（`hts_str2dbl`）——htslib 自己两种类型不一致 | **更正（2026-10-15，对照 VCFv4.5 后）**：初版理由"VCF 规范本身不允许 INFO 值含空格"是**错的**——§1.6.1（tex:431）明文 "Space characters are allowed in values"（**4.5 新增**，changelog），且本轮差分已实测双方对 String 值含空格（`XS=a b c`）处理一致。正确理由：空间的许可是对**值**（主要是 String）的；**数值文法不含空白**——Float 正则（tex:97）首尾均无 `\s`，Integer 是 "32-bit, signed"。` 5` 在数值标签不匹配任何规范文法，属规范外输入；htslib 两种类型的不一致没有设计意图。reed 统一不 trim、missing 处置与规范文法一致，**结论不变、理由更正** |
| **B5** | Number=1 声明却给多值（`XI=5,`） | 0 + 4 | 建议当 **missing**（比较恒假） | 标量读作 **0** | **规范支持 reed 的建议**：§1.4.2（tex:158-159）"if the INFO field contains a single number, then this value must be $1$"——`XI=5,` 在 Number=1 标签上是**规范外输入**。htslib 那个 0 是纯实现细节（`vcf.c:4286` 先 `v1.i = 0`，`len == 1` 分支才填充），无任何注释表明有意。对规范外输入，missing（比较恒假）是最诚实的恢复；0 会参与 `XI<5` 产生反直觉结果。**若用户要求逐条对拍通过**，改为读 0（02 §6.2b 选项 a），代价是把偶然行为固化成契约 |
| **B6** | `-q/-Q` 阈值的 f32 量化 | 5（reed 更精确） | 刻意复现 bcftools 的 `sscanf("%f")` 精度（`af_threshold`，`plan/mod.rs:796-798`，带注释说明边界案例） | 无意为之（C 的 `float` 存储） | 这是**已经对齐**的项，列出是为了防止未来"改进"成 f64：那会在 `AC=34, AN=680` vs `0.05` 这类边界上翻转结果，破坏 gold 对拍 |
| **B7** | 表达式语言子集（下标 `[0]/[*]`、`MAX/MIN/AVG` 函数、算术） | 5 | 报 "not supported yet" | 完整实现 | `cli.md §5` 的既定纪律：**超出子集要报告，绝不用另一种语义静默执行**。这条纪律本身就是 reed 比"部分实现"更可信的地方，不因本轮改变 |
| **B8** | sites-only VCF（8 列） | 4 | 拒绝（`LineRef` 的 9-tab 契约） | 接受 | **移出本轮**（06 §4）：格式级决定，牵连扫描器/编辑器/编码器/index 全链。不是 INFO 缺陷，独立立项候选。在立项前维持"拒绝 + 明确报错"，比部分支持更安全 |

### C. 待拍板（【决】，裁决前不动代码）

| # | 分歧 | 选项 | 规范默认 |
|---|---|---|---|
| C1 | 裸标签（A13） | (a) 对齐 bcftools 恒假 / (b) 编译期拒绝 | **(b)**——错误优于静默恒假 |
| C2 | Number=1 却多值（B5） | (a) 读 0 对齐 bcftools / (b) 当 missing | **(b)**——不把 htslib 偶然行为固化成契约 |
| C3 | `--fill-tags AF/F_MISSING` | (a) 只校验+文档收敛 / (b) 随本轮实现 | **(a)**——实现属新功能，AF 在双等位门下 = `ac/an` 量小，可追加为 N10b |
| C4 | 默认 `drop_stale` | (a) 仅可见化 / (b) 改 opt-in | **(a)**——默认清单是 gold 对拍设计（v2 02 §3.5）；改默认会影响既有用户脚本 |
| C5 | 前缀解析全对齐（B3/B4） | (a) 保持严格 / (b) 加 `ParseMode::Htslib` | **(a)**——改变 gold 之外输入的过滤结果、热路径变慢 |

### D. 统计

- 18 条必须对齐（A1–A18），其中 **A1 有规范明文**（判定档 0：VCFv4.5 §1.3
  的 Integer 值域在 VCF 与 BCF 都被禁止——reed 现状不只是"与 htslib 不同"，
  是写出规范禁止的值）；A2/A3/B5 涉及的输入（Flag 带值、Number=1 却多值）
  本身规范外，对齐属"恢复策略"；其余 2 条属档 2（bcftools 文档承诺）、
  8 条属档 3（reed 违反自己的文档/契约，与 bcftools 无关）。
- **8 条标注不改（B1–B8）**，其中 B4 的理由经 VCFv4.5 对照**更正**（规范
  4.5 起允许值含空格；不 trim 的结论不变，依据改为"数值文法不含空白"）、
  B3 拆出 `+5` 一行转交 N13 顺手对齐（规范 Float 正则显式允许 `+`）、
  B6 是"已对齐、防止未来改坏"的护栏项、B8 是移出本轮的格式级项。
- 5 条待拍板（C1–C5），全部有规范默认值，裁决前不动代码。
- **INFO 解析归属**（影响 A/B 两类的性质）：INFO 值解析**全部是 reed 自写**，
  不经 noodles——`architecture.md:1912` 决策 D1 明文"VCF 样本列/固定列/INFO
  用自写快速扫描（推翻 DESIGN【定】'不写解析代码'的这一部分）"，理由是
  "noodles 每多解一列 FORMAT +28.3%"。noodles 在 `reed-io` 只剩四处：header
  解析（`header.rs:78`，用于建 `TagTable`）、索引写出（csi/tabix）、输出
  header 渲染、BCF 编码器 **oracle**（仅 `REED_BCF_ENCODER=noodles` 下启用，
  `single_pass.rs:38`；生产走 `NativeEncoder`）。
  **这解释了差异的分布**：v10 的结构层能对齐是因为对着 htslib 源码逐条改；
  语义层的 `hts_str2int`/`hts_str2dbl` 从未被参照，且自写了**三份**互不相同的
  实现（D12-14 簇 ι），所以一次性冒出一批。B3/B4 的"不对齐"判定因此是
  有意识的选择，不是遗漏——但前提是 A17 的内部一致性先修好。
- **规范量尺的状态**：本轮参照的 VCFv4.5 是**草案**（用户提供，tex 源）；
  reed 与 htslib 1.24 的 `##fileformat` 实际声明是 4.2/4.3。4.5 新增条款
  （值含空格、percent-encoding）在两边都未实现，本轮只把它们当作**判定
  依据**（判断哪些输入属规范外），不要求实现（06 §4.2）。

---

## 3. 复核门（本轮新增，落地后由 `ci/check-doc-drift.sh` 登记）

> 守卫 1–14 对应 §1 的缺陷修复。**§2.5-B 的 8 条"不改"项没有守卫，且不得
> 新增**——它们是永久分歧，写进 `cli.md §5` 差异表即为落地形式（守卫 15）。
> 下一轮复审若把它们重新判为缺陷，需先推翻 §2.5 的判定规则本身。

1. **D12-1 守卫**：驱动的策略覆盖尊重尾部表达式——`reed/src/filter/mod.rs`
   含对 `TailRule::Expr`（或 plan 暴露的等价谓词）的检查；
   `--mind + --indep-pairwise + 尾部表达式`存在 usage 报错文本。
2. **D12-2 守卫**：`fill_info_set` 或其调用点含 `has_gt` 门（grep）。
3. **D12-3 守卫**：`parse_int_vector` 引用 `INT32_MIN_VALUE`/`INT32_MAX_VALUE`。
4. **D12-4/5 守卫**：`encode_info_text` 不再含 `(InfoType::Flag, _) =>` 的
   无条件 NULL 分支与 `if v == b"."` 空串特判（或两者收口进单一 helper，
   helper 内含 htslib 语义注释）。
5. **D12-6 守卫**：`MinAfInfo`/`MaxAfInfo` 的推入条件含 cohort 门；`-q/-Q`
   路径含"陈旧 INFO"警告文本。
6. **D12-7 守卫**：`expr::check` 的签名含类型信息（`TagTable` 或等价）。
7. **D12-8 守卫**：`eval_info` 的 `any` 分支含 `cmp.apply`（grep
   `any |= ` 或等价形态）。
8. **D12-9 守卫**：`native_encode.rs` 含计数上限检查 helper（单一实现，
   文本/BCF 两路径共用）。
9. **D12-10 守卫**：head 的 ALT 计数不含裸 `memchr_iter(b',', alts)`
   （走共享 helper）；`max_alt_alleles` 为 `Option`（或 `-M 1` 有测试锁定）。
10. **D12-11 守卫**：translate 对 `--fill-tags` 未知标签返回 usage 错误
    （测试名锁定）。
11. **D12-12 守卫**：`--help` 的 `-x` 条目含"默认"字样；README 含默认清单。
12. **D12-13 守卫**：AC/AN INFO 读取单一实现（`grep -c "fn allele_counts"`
    全仓 ≤ 1 处生产定义）。
13. **D12-14 守卫**：INFO 数值解析单一实现；cli.md §5 含解析差异表条目。
14. **D12-15 守卫**：`line_edit.rs` 不含 `LineRef::parse(line.line)`。
15. **"不改"项登记守卫**（§2.5-B）：`cli.md §5` 的差异表含 B1–B5 全部五条
    （文本保真、Σac>an、前缀解析——`+5` 除外，tex:97 属规范内、N13 顺手
    对齐——、空格容忍、Number=1 多值），每条注明"reed 有意为之 + 理由"。
    **这是"不改"的落地形式**——不登记就会在下一轮复审被重新判为缺陷。
    反向断言：测试里不得出现"B 类永久分歧形态与 bcftools 结果相同"的
    对拍（那会把分歧变成契约）；B6（f32 量化）例外——它是对齐项的护栏。
