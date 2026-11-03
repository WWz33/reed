# 审查条目复核：证据、判定与同类位置（doc_v23，2026-11-03）

> 被复核文档：`reed-rs-review-2026-11-03.md`（仓库外）。本文只记录复核，
> 不含代码改动。行号按当前工作区（含 doc_v22 R2 提交 `8b54588`）。
>
> 判定口径：
> - **属实**：按报告给的条件可复现，或代码路径必然触发。
> - **属实（范围修正）**：问题存在，但触发条件/严重度/可达性与报告不同。
> - **不成立**：按报告的条件不可复现，或与参照实现（bcftools/htslib）一致。
> - **代码复核**：只做了静态路径确认，未构造运行时复现。
> - **性能/精简**：不影响正确性，单独分类。

## 批 1 ①：reed-io VCF 文本层

### ①1 尾部空列 + 样本区编辑 panic —— 属实（触发条件修正）

**代码路径**：`crates/reed-io/src/vcf/line_edit.rs:192-241`。`edit_samples`
遍历整个样本区，`sample_idx` 无上界；`SampleMask::contains` /
`Bitset::get` 对越界索引 `assert`（`reed-core/src/bitset.rs:59`）。

**实测**（debug 构建）：

```text
# 单尾随 tab（合法，bcftools 忽略）：不触发
chr1 10 ... GT 0/1 1/1 \t
$ reed filter in.vcf -s s0 -Ov        # exit 0，输出 1 行

# 两个尾随 tab（多一个空样本列）+ 样本编辑：触发
chr1 10 ... GT 0/1 1/1 \t \t
$ reed filter in.vcf -s s0 -Ov
thread '<unnamed>' panicked at crates/reed-core/src/bitset.rs:59:
invariant: bit index out of range
reed: io error: worker panicked on chunk 0
# exit 1

# 掩码路径同样触发（不需要 -s）
$ reed filter trailing_gq.vcf --min-gq 20 -Ov
# 同上 panic（Bitset::get）
```

**修正**：报告说"行以 `...\t` 结尾"即可；实际需要**多一个空列**（行尾
连续两个 tab）才会让循环进入越界字段。单个尾随 tab 在
`field_start = field_end + 1` 后正好等于 `samples.len()`，循环结束。

**其他考虑**：
- 触发面是"任何 `plan.samples_edited()`"：`-s/-S/--keep/--remove` 走
  `keep.contains`，`--min-dp/--min-gq/--set-gts` 走 `mask.get`；两者都会
  越界。
- `edit_samples` 是文本路径；BCF 输入的原生编辑不经过它（但
  `filter` 的 BCF 文本重渲染仍可能走 `LineRef` 编辑，需确认）。
- 修法：给 `edit_vcf_line`/`edit_samples` 传 header 的 `n_samples`，循环内
  `if sample_idx >= n_samples { break; }`（超出部分按 R22-3 契约只能为空）。
  测试：0/1/2 个尾随空列 × `-s`/`--min-gq`。

### ①2 `format_g6_into` 科学计数分支丢负号 —— 属实（严重度上调）

**代码路径**：`crates/reed-io/src/vcf/line_edit.rs:316-394`。指数分支用
`write!(wide, "{:.5e}", v.abs())`，随后只输出 `mantissa`/指数，没有补
`-`；定点分支用 `{v:.precision$}`，符号保留。

**实测**（VCF → `bcftools view -Ob` → `reed filter -Ov`）：

```text
输入 INFO/XF            bcftools view -Ov        reed filter -Ov
XF=-1.5e-7              XF=-1.5e-07              XF=1.5e-07
XF=-1234567             XF=-1.23457e+06          XF=1.23457e+06
XF=-0.00001             XF=-1e-05                XF=1e-05
```

**修正**：报告称"当前调用方恰好不产负值，属潜伏"。实际上
`format_g6_into` 同时是 **BCF→VCF 的 INFO/FORMAT 浮点渲染器**
（`bcf/render.rs:441` 的 `float_into`，以及
`vcf/line_edit.rs:67` 的 `render_info_value`）。任何带负 Float INFO 的
BCF 经 `reed filter -Ov` 都会静默改写；这是输入数据损坏，不只是内部
producer 的潜在问题。`format_g6`（String 版，`line_edit.rs:275`）同源。

**其他考虑**：
- 触发条件是"值取科学计数法"（6 位有效数字四舍五入后指数在 `[-4, 6)` 之外），
  不是所有负数：`-0.123` 走定点分支，符号正确。
- `-0.0` 单独处理（保留负零），不受影响。
- 修法：指数分支在 `mantissa` 前 `if v < 0.0 { out.push(b'-'); }`；
  金标测试补 3 个负值 + `-0.0` + 指数边界（`-0.000099999`/`-999999.5`）。

### ①3 FORMAT 布局缓存 64 上限 —— 性能（代码复核，属实但非 bug）

`format_layout.rs:109` `CACHE_LIMIT = 64`，达到后 `layouts.clear()` +
`index.clear()`（`151-157`）。>64 种 FORMAT 字符串交替的文件命中率归零，
每 record 4 次分配。1000 个不同 FORMAT 的真实 callset 罕见（通常
`GT`、`GT:AD`、`GT:DP:...` 几种），列为性能改进；改 LRU 或提高上限前先
profile。

### ①4 `pruned_header` 零测试 —— 测试缺口（代码复核）

`header.rs:113-128` 是样本子集输出 header 的唯一入口；当前无测试。
建议补两个：保序、`keep` 短于样本数时默认保留。不是运行时 bug。

### ①5 五个 pub 死代码 —— 精简（代码复核）

`is_pass`、`samples_start`、`alt_at`、`info_value`、`format_g6(String)`
在 `crates/` 内除测试/自身外无调用方（`grep` 复核）。删除或
`#[cfg(test)]`。不影响行为。

### ①6 样本列数错误丢弃 found/expected —— 诊断体验（代码复核）

`genotype_scan.rs` 把 `SampleColumnsError::{TooFew, TooMany}` 统一映射为
字符串 `"sample column count disagrees with the header"`，丢掉了
`found/expected`。`SampleColumnsError` 已有字段，格式化即可。建议修，
纯错误信息。

### ①7 空 ALT 计为 1 个 alt —— 属实（畸形输入，低优先级）

`vcf/site.rs:27-33`：`count_alts(b"")` 返回 `commas + 1 = 1`，只有 `.` 特判
为 0。空 ALT 列是畸形记录；若上游放行，`-m/-M`、`N_ALT`、`-q/-Q` 会把
它当 1 个 ALT。建议空与 `.` 都返回 0，并补单测；是否拒绝空 ALT 另议。

## 批 1 ②：reed-io BCF/BGZF/chunk/write

### ②1 `Cursor::typed`/`typed_int` 互递归 —— 属实，P0

**代码路径**：`bcf/record.rs:313-342`：

```rust
pub fn typed(&mut self) -> ... {
    ...
    if count == 15 {
        let ext = self.typed_int()?;      // typed_int -> typed
        ...
    }
}
pub fn typed_int(&mut self) -> ... {
    let (ty, count) = self.typed()?;      // 再次进入 typed
    ...
}
```

连续 `0xF1`（count=15, type=INT8）每字节产生一次 `typed → typed_int → typed`
递归。`format_fields`（`record.rs:209-226`）与 INFO/FILTER 解析均走此路径。

**实测**：从正常 BCF 解压，把第一条记录的 `l_indiv` 改成 200000、indiv 区
填满 `0xF1`，重新 bgzip（压缩后 758 字节）：

```text
$ reed filter corrupt.bcf -Ov
thread '<unknown>' has overflowed its stack
fatal runtime error: stack overflow, aborting
# exit 134 (SIGABRT)
```

**其他考虑**：
- htslib 只允许**一层**扩展：`bcf_dec_typed_int1` 解析扩展描述符后不再
  递归；修复可以照此：`typed()` 在 count==15 时直接读下一字节，并拒绝其
  count 半字节再次为 15（或任何非法 type）。约 5 行。
- 同类入口：`format_fields`、INFO/FILTER 描述符、`decode_record_gt` 的
  key 解析都经过 `typed_int`；修在 `typed` 一处即可全覆盖。
- 不是合法输入；属于"损坏/恶意文件不得让进程 abort"的健壮性门。

### ②2 CHAR 类型 FORMAT 渲染/列解码 panic —— 属实（范围修正）

**代码路径**：
- `bcf/render.rs:390-425` `fmt_gt`：`width = type_width(ty)`；CHAR 宽度 1，
  `int_element(ty, &data[i*width..(i+1)*width])` 走 `_` 分支读 4 字节
  （`record.rs:51-57`）→ 越界 panic。
- `bcf/columns.rs:158-176` `float_element_column`：非 FLOAT 时同样调
  `int_element`（1 字节元素读 4 字节）。
- `bcf/gt.rs:163-171` 的 `decode_slow` **已有**守卫：
  `"BCF GT field is not integer-typed"`。

**实测**（把 GT 的 FORMAT 描述符从 `0x21` 改成 `0x27`，count=2、type=CHAR）：

```text
$ reed filter charfmt.bcf -Ov
reed: io error: malformed record at chr1:10: BCF GT field is not integer-typed
# exit 1，干净错误（filter 的 GT 解码先命中守卫）

$ reed convert charfmt.bcf --012 -o cv
thread '<unnamed>' panicked at record.rs:55:53:
index out of bounds: the len is 1 but the index is 1
reed: worker panicked on chunk 0
# exit 1

$ reed stats charfmt.bcf --freq -o q
# 同一 panic
```

**修正**：报告说 `render_sample_part` 对 GT CHAR panic——在 `filter` 路径
不成立（守卫先触发）；`convert`/`stats` 的渲染路径确实 panic（`fmt_gt`）。
float 列路径（`float_element_column`）在本次 CLI 探针里未触发：带
`FMT/XF` 掩码的 BCF 走的是另一条路径（输出把 CHAR 字节渲染成 NUL，
掩码按 missing 处理），没有进入原生列解码；**函数级**仍是同一越界模式，
建议直接加类型判断（非整数/非浮点 → malformed/NaN），并用函数级测试
钉住，而不是只依赖某条 CLI 路径。

**其他考虑**：
- `fmt_array`/`fmt_array1` 对 CHAR 有专门分支，非 GT 的 CHAR 字段渲染不
  panic；风险集中在 `fmt_gt`（GT）和 `float_element_column`（请求为
  float 的键）。
- 合法 BCF 的 GT 永远是整数类型；这是 malformed 输入的错误处理问题，
  不是合法数据兼容问题。
- 修法：`fmt_gt` 对非 INT8/16/32 返回 `MalformedRecord`（与 gt.rs 同款）；
  `float_element_column` 对非 FLOAT/INT 返回 `NaN`；补最小用例。

### ②3 热路径每记录/每字段堆分配 —— 性能（代码复核）

`render.rs:222/236` 的 `format_fields` Vec、`columns.rs:66` 的字段查找、
`encode_native.rs:382/445` 的 keys Vec。与模块注释 "no per-record
allocation" 不符，属性能债；先 profile 再动。

### ②4 `RegionChunkSource::pending_block` 死机制 —— 精简（代码复核）

`chunk/region.rs:42-43,75-76,128-137,193-196`：全仓库只有
`pending_block: None` 初始化与 `take()`，没有赋值 `Some` 的路径。

### ②5 畸形类型/扩展计数无失败测试 —— 测试缺口（代码复核）

正是 ②1/②2 的触发面；随修复补最小 case 即可。

### ②6 死掉的公开 API —— 精简（代码复核）

`bgzf/deflate.rs:91` `eof_block()`、`write/block_writer.rs:70-77` 的
`OutputFormat::BgzfLevel0`、`encode_native.rs:139-143`、
`encode_noodles.rs:123-128`、`chunk/plain.rs:26,75,103`。是否删由 API 稳定
承诺决定（pub 项删除是 breaking change）。

### ②7 `leak()` 永久泄漏错误消息 —— 低优先级（代码复核）

`encode_noodles.rs:160-164` 每条编码错误 `Box::leak`。当前编码错误是
fatal（每 run 至多一条），不构成无界增长；若未来改成逐记录容错才是问题。

### ②8 `field()` 在引号内匹配 `ID=`/`IDX=` —— 低优先级（代码复核）

`bcf/dict.rs:162-192` 的字符串扫描不识别引号；畸形 header 可能注册错名。
真实文件罕见。

## 批 1 ③：reed-core

### ③1 `round_half_even` 容差判平局 —— 属实，P0/P1

**代码路径**：`reed-core/src/math/r_compat.rs:21-40`：

```rust
if (frac - 0.5).abs() < f64::EPSILON * x.abs().max(1.0) { /* tie */ }
```

**实测**（探针 crate 调真实函数）：

```text
round_half_even(1.4999999999999998) = 2    R: 1
round_half_even(2.5000000000000004) = 2    R: 3
round_half_even(0.35 * 90)          = 32   R: 31
```

R 值的依据：R 的 `round()` 只在二进制值正好等于 `n+0.5` 时取偶；上面三个
都不是精确平局。`0.35*90 = 31.499999999999996`，`xvalDapc` 的训练集大小
因此 off-by-one。

**其他考虑**：
- 修法是精确比较 `frac == 0.5`（`x - floor(x)` 对有限 f64 是精确的），
  不是缩小容差；容差在任何量级上都会把相邻可表示值误判。
- 负数行为当前正确：`-2.5 → -2`、`-0.4 → -0.0`；修复不要改
  `floor`/奇偶逻辑。
- `round_half_even_usize` 的 `r <= usize::MAX as f64` 断言放行 `2^64`
  附近的饱和（报告 ③11）；调用方都是小计数，纯加固。

### ③2 `dot`/`dist2` 长度不匹配静默截断 —— 属实（内部不变式，低优先级）

`reed-core/src/linalg/mod.rs:151-166` 用 `zip`，短的一方截断；同文件
`matvec`/`tmatvec` 有 `assert_eq!`。建议补断言（debug/release 一致），
避免库调用方拿到"合法但错误"的点积。

### ③3 NeumaierSum 文档过度声明 —— 文档问题（代码复核）

`math/sum.rs:1-6` 写"regardless of chunk boundaries 逐位一致"；实际是
"固定分块逐位一致，改分块 ≤4 ulp"。与 spec §8.4 的承诺对齐即可，不改
数值实现。

### ③4 `argsort_desc` 非稳定降序 —— 属实（低优先级）

`linalg/mod.rs:190-197`：升序稳定排序后整体 `reverse()`，并列键的相对
顺序被反转。改为直接 `sort_by(|a,b| v[b].partial_cmp(&v[a])...)` 的稳定
降序。当前调用方是否依赖并列顺序需逐一确认（先 grep 使用点）。

### ③5 `bp_windows` 无溢出保护 —— 加固（代码复核）

`grid.rs:136,144` 的 `start + size`/`start += step`。真实基因组不可达；
CLI 的巨大 `--window` 实测未触发（窗口锚点/循环先结束）。可用
`checked_add`/`saturating_add` 加固，不列为 bug。

### ③6 组文件解析线性扫描 —— 性能（代码复核）

`group.rs:128,158` 每行扫组名；照 `SampleTable` 的 HashMap 模式即可。
组数通常很小，优先级低。

### ③7 `RegionSet::stab` 命中路径每次分配 —— 性能（代码复核）

`grid.rs:270-291` 每次探测 `Vec`；百万位点即百万次 malloc。改
`&mut Vec` 复用。需要 benchmark 确认热点。

### ③8 `HweCache::distribution` 每次克隆 —— 性能（代码复核）

`math/hwe.rs:153-163` 返回 `HweDistribution` 克隆。可返回引用/Arc；缓存
容量 64，克隆成本与 n 相关。先 profile。

### ③9 IUPAC 全表对拍测试 —— 测试缺口（报告已手工比对 62+10 条）

当前实现与 vcf2phylip 2.9 一致；补一条全表断言防回归。

### ③10 `pretty_sequence` 非 10 回退 —— 未实现路径（代码复核）

`r_compat.rs:150-155` 对 `n != 10` 回退为全量整数列；pub API 但
`xvalDapc` 只用 `n=10`。加 `debug_assert`/显式文档，或实现 pretty 步长
梯子并补对拍。

### ③11 `round_half_even_usize` 范围断言 — 见 ③1 其他考虑。

### ③12 测试覆盖总体良好 — 与报告一致；缺口即 ③1/③9/③10。

## 批 2 ④：reed-filter

### ④1 AF/MAF AN=0 的 NaN `!=` —— **不成立（驳回）**

**报告主张**：`AN=0` 时 AF 返回 NaN，`NaN != 0` 为 true，违反"缺失值比较
为假"的契约，建议改为 `Val::Missing`。

**实测**（`AF` 由 GT 计算，两个位点：全缺失 / 0/1+0/0）：

```text
表达式     bcftools view -i     reed filter -i
AF>0       20                   20
AF<1       20                   20
AF==0      （空）               （空）
AF!=0      10 20                10 20
MAF!=0     10 20                10 20
```

bcftools 在同一位点也把全缺失的 AF 当 NaN，`!=` 为 true。**reed 与
bcftools 逐条一致**；改成 Missing 会让 `AF!=0` 只留真 AF>0 的位点，反而
偏离参照实现。spec 的"与缺失值比较为假"指的是 `"."` 字面缺失；NaN 是
计算值。建议只在文档中说明 `AF/MAF` 在 AN=0 时为 NaN，不改语义；若要
改，需要单开决策并同步 bcftools 兼容性说明。

**拍板（2026-11-03）**：按推荐 A——**不改语义**；文档写明 AN=0 时
AF/MAF 为 NaN、比较按 IEEE 语义，并补测试钉住与 bcftools 一致的结果。

### ④2 AC 多等位只取第一个元素 —— 属实

**拍板（2026-11-03）**：**丢弃（不修）**——不接受多等位，该条从 v23 修复
清单丢弃。现状保留：默认 `BiallelicGate` 已丢弃多等位；显式 `-M >2` 时
AC 仍按首 ALT（不支持路径），不实现 any-pair。

**实测**（`AC=0,3` / `AC=3,0`，`-M 3` 放开多等位）：

```text
表达式   bcftools -i   reed -i
AC>0     10 20         20          ← 位点 10 的 AC 向量 [0,3] 被 reed 当第一个 0
XX>0     10 20         10 20       ← 自定义 Number=A 标签走 any-pair，正确
```

**代码**：`expr/mod.rs:390` 的 `Field::Ac` 走 `info_or_gt`（`info_value`
只取第一个逗号元素）的标量视图；
自定义 INFO 标签走 `info_operand`（Number=A → any-pair）。修法：当
`info_trusted` 时 `Field::Ac` 改走 `info_operand`，否则保留 GT 重算路径。
需同时确认 AN/AF/MAC/MAF 是否有同类；AF/AN/MAC/MAF 是 Number=1，无此问题。

### ④3 表达式嵌套无深度上限 —— 属实，P0/P1

**代码**：`expr/parser.rs:295-380` 的 `logic_or → logic_and → unary →
comparison → primary → LParen → logic_or` 纯递归，无计数器。

**实测**：10,000 层括号（20 KB 表达式）：

```text
$ reed filter in.vcf -i "$(python3 -c 'print("("*10000+"1"+")"*10000)')"
thread 'main' has overflowed its stack
fatal runtime error: stack overflow, aborting
# exit 134
```

**其他考虑**：
- 同一深度也会让 AST 的 `Drop` 递归；解析器加 depth 上限（如 500）即可
  同时约束求值和析构。
- 输入是操作者自己给的命令行表达式，不是远端数据；严重度低于文件解析，
  但仍是"用户输入导致 abort"。
- 修法：`Parser` 加 `depth: usize`，`primary` 的 `LParen` 分支递增/递减，
  超限返回 `Err`。测试：499/500/501 层。

### ④4 类型化 INFO 规则把 `TAG=.` 绕过 on_missing —— 属实（语义不一致）

`eval/head.rs:242-274`：标签存在但值为 `.` 时 `parse_number` 返回 None，
比较为 false；只有标签缺席才走 `on_missing`。因此
`DP=.` + `on_missing=pass` 与"DP 缺席"行为不同。建议把 `.`/空值统一走
`on_missing`；补三种形态的用例（缺席 / `.` / 正常值）。

### ④5 `to_args` 把 `<=`/`>=` 渲染成 `= N` —— 属实（可复现性）

`spec/render.rs:114-128`：`trim_start_matches(['<','>'])` 对 `<= 10`
得到 `= 10`，生成 `--min-dp "= 10"`，重放时 CLI 报错/语义错。修法：先匹配
两字符算符，再匹配单字符。测试：`< <= > >=` 四种。

### ④6-④8、④10、④13、④14 —— 性能/清理（代码复核）

- ④6 每条 Expr 头规则每 record 重建 `SiteValues`（`eval/head.rs:129-140`）。
- ④7 `info_or_gt` 扫 INFO 两遍（`expr/mod.rs:435-439`）。
- ④8 字符串/标量比较每 record 分配（`expr/mod.rs:186-193,508-516,693-694`）。
- ④10 `af_of` 非饱和减法（`eval/tail.rs:141`）——改 `saturating_sub`，
  一行；当前输入下未触发下溢。
- ④13 注释重复/矛盾；④14 `stale_for_record` 每 record Vec
  （`eval/edit_plan.rs:65-78`）。

### ④9 reed-filter 无 tests/ —— 测试缺口（属实）

发现 ④1/④2/④4/⑤ 均为 parser/evaluator 风险分支；建议先补 4 个用例
（AC 向量、`TAG=.`、`to_args` 算符、深度上限），再决定是否建 tests/。

### ④11 空 `SampleSelect` 当成队列改变 —— 属实（JSON 边界）

`plan/mod.rs:508-511` 用 `spec.sample_select.is_some()`；空对象反序列化后
`is_some()` 为真，导致 `-q/-Q` 从 INFO 短路切到 GT 重算。改
`is_some_and(|s| !s.is_empty())`。仅 JSON spec 可达。

### ④12 含 `:` 的 contig 无法指定区间 —— 属实（限制）

`spec/region.rs:83-90` 直接 `split_once(':')`。修法：先对 header contig 表
整串匹配（最长匹配）再切分。HLA contig 名是现实场景；优先级中。

## 批 2 ⑤：reed 主 crate

### ⑤1 sketch 形状 pass0 预算只按单 contig —— 属实（代码复核）

`filter/pass0/matrix.rs:917-946`：`close_pending` 只比较
`matrix.size_bytes() > budget`，随后 `self.matrices.push(matrix)`；
`track_chunks`（sketch 形状）把所有 contig 的矩阵常驻，N 条各略低于预算
即总量 N×预算。报告的例子合理。修法：累计 `self.matrices` 的
`size_bytes()` 与预算比较，超限报 `MemBudgetExceeded`；补多 contig fixture
（小预算、3 条 contig）。非 sketch 路径有 `max_prune` 上限，需单独核对
每个 prune 线程的矩阵驻留。

### ⑤2 `--ld-pairs` 表无上限 —— 属实（CLI 可达；已拍板 C）

`sink.rs:248-267` 把每个 contig 的 pair 集合收进 `LdResult.pairs`，按距离
窗口平方增长；`--ld-decay 100k --ld-pairs` 在全基因组上可达 GB 级。
`--ld-pairs` 单独使用已被 CLI 拒绝（见 ⑤3），但"同时给两个标志"仍无上限。

**拍板（2026-11-03）**：**选项 C**——sink 按 contig 顺序直接追加到最终
`PREFIX.ld`：**首个 pair 到达时懒创建并写表头**（全无 pair 仍不创建文件，
与现状一致），每个 contig 的 pair 行按现有顺序写出，不再在内存里聚合成
`Vec`/`String`。不落临时文件，输出与现状逐字节一致，内存 O(单 contig)。
`--ld-pairs` 单独使用仍拒绝；落地补逐字节对比测试。

### ⑤3 LD decay 默认值不一致 —— 属实（仅库 API）

`request.rs:402` 的 `decay_edges()` 在 `decay_max_bp=None` 时用 `bin`；
`sink.rs:247` 用 `unwrap_or(100_000)` 收集 pair。CLI 的
`build_ld_request` 对 `--ld-decay` 总是 `Some`，且拒绝 pairs-only
（"pass --ld-decay as well"），所以 CLI 不可达；库调用方传
`LdRequest { decay_max_bp: None, pairs: true }` 会得到"收集到 100kb、
但 decay 只有 0..bin 一箱"的结果。修法：`pairs=false && None` 时直接
不建 LD plan，或统一 `decay_max_bp` 的缺省。

### ⑤4 `read_reference` 整读 FASTA —— 属实（内存，代码复核）

`convert/job.rs:1514-1517` `fs::read` 整个参考，随后按 contig 再复制进
`Arc<[u8]>`：人类基因组峰值约 2×3 GB。修法：`BufRead` 逐行解析，按 contig
累积；或 mmap。影响 `convert --msa --reference`。

### ⑤5-⑤10 热路径分配/死代码 —— 性能/精简（代码复核）

包括 mid stage 缓冲、每 record `auto` Vec、prune 位图、`PcaScratch` 死字段、
`append_sites` 整块复制、`--rm-dup-ids` HashSet 驻留。均不影响正确性；
`append_sites` 的 `mem::take` 是零风险小改。

## 批 2 ⑥：reed-stats

### ⑥1 LD decay 距离在首 edge 前计为 bin 0 —— 属实（库 API，低）

`ld.rs:164-170`：`Err(0).saturating_sub(1) = 0`。仓库内调用方 edge 从 0 起，
不可触发；pub API 传非 0 起始 edges 时会错。修法：返回
`Option<usize>`/哨兵，`push` 跳过。

### ⑥2 `pos_log` 无界增长 —— 属实（范围修正）

`ld.rs:281,379,438,542`：每个 arrival `push`，`evaluate_window_bp` 只从
`window_start+1` 向后找，旧条目永不访问。但 `finish_contig`
（`ld.rs:542`）会 `clear()`——所以是**单 contig 内**增长，不是全文。
量级：一条 1 亿位点的 contig ≈800 MB/线程；人类最大染色体约 1–2 千万
位点 ≈80–160 MB/线程。修法：窗口推进时把 `[0, window_start]` 截断（
`drain(..new_start)` 或保留偏移），或改用 ring buffer。补大 pos span 单测。

### ⑥3 `WcMultiAcc::push` 每位点 4 次分配 —— 性能（代码复核）

`fst.rs:304-348`：r 通常 ≤10，栈/ scratch 复用可省；报告称约 2× 减速，
需 benchmark 确认。

### ⑥4 merge 大小不匹配 release 静默 —— 属实（防线，低）

`request.rs:289-306`、`window.rs:95-97`、`jackknife.rs:281` 只在
`debug_assert` 检查；release 下不匹配的 `zip` 截断、`match` 落 `_ => {}`。
这些都是内部不变式（同一请求的 partial），CLI 不应触发；建议改
`assert!`/`Result`，防止未来回归。

### ⑥5 两个分支缺数值对拍 —— 测试缺口（代码复核）

WcMulti r≥3 半调用组、F2/F3/F4 merge 律。补 proptest/固定 fixture。

### ⑥6-⑥8 冗余守卫/尺度/测试死代码 —— 清理（低）

pair 与 multi WC84 分量尺度差 2×（比率一致，已有对拍）；PBS 三重守卫；
`projected_sfs` 死代码。不影响正确结果。

## 批 3 ⑦：oracle

### ⑦1 SFS/θW/Tajima complete-case vs allel —— 已拍板 A（测试口径）

`oracle/src/estimators/spectrum.rs:3-7,108-120,268-282` 是整站跳过缺失；
allel 的逐位点删除语义可能不同。报告推断未实测 allel。该口径原本需要
决策：oracle 应复刻 allel（对拍用）还是定义 reed 语义（生产用）；补含
缺失位点的 golden 之前不能确定哪边是"真值"。

**拍板（2026-11-03）**：按推荐 A——oracle 保持 **complete-case**（与生产
实现和 `docs/spec/architecture.md` 的 SFS 口径一致）；补含缺失位点的
golden；在 oracle 文档记录与 allel 逐位点删除语义的差异，不把 allel
当作该族的唯一真值。

### ⑦2 WC84 单组全缺失分支 —— 需实测（oracle）

`oracle/src/estimators/fst.rs:160-171` 保留 c；allel 可能整站 NaN。
建议构造 B 组全缺失 fixture 实跑 allel 后钉死。

### ⑦3-⑦5 golden 覆盖/`is_variant`/hobs-hexp 样本集 —— 测试缺口

均只影响对拍置信度，不是生产数值 bug。`hobs`/`hexp` 用不同样本集是否
有意需作者确认（oracle 文档）。

## 批 3 ⑧：reed-reduce / reed-distance / reed-engine

### ⑧1 三处 NaN 排序 panic —— 属实（库 API，P1）

**代码**：`graph.rs:326`（`isomap`）、`grm_reduce.rs:164`
（`knn_from_distances`）、`grm_reduce.rs:207`（`isomap_from_distances`）
均为 `partial_cmp(...).unwrap()`。`pairwise_distances` 对 NaN 行输出 NaN。
同文件 `knn_graph`（`graph.rs:189-192`）、`umap.rs:190`、`extra.rs:395`
已改 `total_cmp`；`auto_gamma`（`graph.rs:172`）有 `d > 0.0` 过滤，NaN 不
入列，安全。

**实测**（探针调 `reed_reduce::isomap`，矩阵含一 NaN 行）：

```text
thread 'main' panicked at graph.rs:326:52:
called `Option::unwrap()` on a `None` value
# exit 101
```

**可达性**：reed CLI 目前没有 isomap 入口（`grep` 全仓，仅 lib re-export），
所以不是用户命令崩溃，而是库 API。修法：三处 `total_cmp`，各一行；补
与 `s31_knn_graph_sort_tolerates_nan_genotypes` 同款的测试。

### ⑧2 特征向量无符号约定 —— 属实（P1；已拍板 A）

**代码**：`backend.rs:56-69` `top_k` 原样拷贝列；`pca.rs:80-90` `pack`
不规范化。LDA 已做"最大分量取正"（`lda.rs:16-22,489-503`），PCA 没有。

**实测**（探针，同一 6×6 矩阵，两个 backend）：

```text
Fallback: vals=[10.00358088059565, 9.000771282019596]
          pc1 first=-0.99845624390341792  sum=-1.08657588832281937
          pc2 first= 0.05225101460837055 sum=-1.02606657870634099
Faer:     vals=[10.003580880595647, 9.0007712820196]
          pc1 first=+0.99845624390341725  sum=+1.08657588832281804
          pc2 first=-0.05225101460836980 sum=+1.02606657870634255
```

PC1/PC2 整体变号。`reed-reduce/src/lib.rs:7-9` 承诺"Both produce identical
results"；切后端/升级 faer 会改 `.eigenvec` 字节。修法：照 LDA 对每列做
"最大绝对值分量取正"；但这会改变现有输出的符号（包括测试 fixture 的
期望值），需要拍板接受这个 breaking change。

**拍板（2026-11-03）**：按结果正确性改——在 `SymEigResult::top_k` 统一
规范化（最大 \|分量\| 取正，并列取行序最先，全零列不动）；`eig::SymEig::top_k`
同规则。PCA、PCoA、kernel PCA、Isomap、Laplacian、diffusion 共用。接受
`.eigenvec`/embedding 坐标的符号变化，以换取跨 backend 逐位一致；落地
需补 V23-48 的双 backend 对拍。

### ⑧3 frame_chain 看门狗检查顺序 —— 代码顺序属实，误报可达性待定

`frame_chain.rs:160-171` 先比较 `silence_ns`，再检查 cell 是否就绪。
若请求的 cell 已发布但 `last_publish_ns` 已超过 120 s（例如 waiter 线程
被极端调度延迟），会报 "chain stalled" 尽管数据就绪。正常调度下每次
publish 都刷新时间戳，误报需要"已就绪 cell 比上限还老"的组合；建议把
cell 检查移到 watchdog 之前（零风险），或消费成功时刷新时间戳。

### ⑧4 tql2 不收敛静默返回 —— 属实（低概率，防线）

`eig.rs:172-177`：`iter > 100 { break; }` 后照常返回。建议
`debug_assert`/`Result` 或提高迭代上限并在文档说明。构造触发矩阵困难，
但静默错值比报错危险。

### ⑧5 GRM 单线程三重循环 + `gram_divide` 死路径 —— 事实属实（性能）

- `backend.rs:73-94` 的 `gram_divide` 只有自身测试和 lib re-export 调用；
- `pca_dense`（`pca.rs:100`）走 `Grm::mean_impute`
  （`reed-distance/src/grm.rs:36-52`），是标量三重循环；
- 4000×20 万 SNP 的 GRM 慢一到两个数量级。
修法：`Grm::mean_impute` 在 `faer-backend` 下委托 `gram_divide`（注意
`gram_divide` 是 `W W'/m`，与 GRM 的标准化后矩阵一致）；需补
backend 对拍。

### ⑧6 `matrices::grm` 全量拷贝 + 二次扫描 —— 属实（内存炸弹）

`reed-distance/src/matrices.rs:303-315` 分配 `n*m` 的 `x` 拷贝
（2000×100 万 × 8 B = 16 GB），再对每对做 O(m)；`st.w` 本身还在。
修法：直接读 `st.w`，按需要局部乘 `sd`（或预计算每行缩放），
`called_both` 在同一个内积循环里累加。需要保证 pairwise 分母不变；
补与现实现逐值对拍。

### ⑧7 矩阵写出未缓冲 —— 性能（代码复核）

`reed-distance/src/output.rs:60,79-107` 每值 `write!` 到裸 `File`；
`BufWriter` 一行。`write_bin` 已用 Vec 缓冲。

### ⑧8 缺 faer/fallback 与跨线程对拍 —— 测试缺口（属实）

`backend.rs` 只有 residual/trace/gram 近似测试，没有两 backend 的逐值
对比；`e2e_pca` 也不切后端。若修 ⑧2，这个测试是前置。

### ⑧9-⑧13 其它 —— 加固/语义（代码复核）

- ⑨ diffusion_map 平稳向量阈值绝对 `1e-9`（`graph.rs:291`）→ 相对
  `1e-9*λ_max`。
- ⑩ kernel PCA `abs().sqrt()`（`graph.rs:117`、`grm_reduce.rs:71,282`）→
  `max(0.0).sqrt()`；负特征值不应变成实坐标。
- ⑪ `standardize` 只认 NaN（`reed-distance/src/standard.rs:110-152`）：
  Inf 会被当有效值，`mu`/`sd` 变 Inf/NaN，整列静默变 NaN。改为
  `is_finite()`。
- ⑫ LLE 的 M 矩阵 O(n³) 稠密构造（`extra.rs:431-444`）——性能。
- ⑬ `pca::backend()` 文档说"served the most recent call"，实现返回编译期
  常量（`pca.rs:104-107`）——文档修正。

## 批 3 ⑨：reed-cli / reed-graph / bench

### ⑨1 `filter -r` 存在时静默忽略 `-R` —— 属实（P1，静默缺数据）

**实测**（`two.vcf.gz`：chr1/chr2 各一条；`r12.txt` 含两行）：

```text
$ reed filter two.vcf.gz -r chr1 -Ov            → chr1 10
$ reed filter two.vcf.gz -R r12.txt -Ov         → chr1 10, chr2 10
$ reed filter two.vcf.gz -r chr1 -R r12.txt -Ov → chr1 10     ← 丢 chr2
$ reed stats  two.vcf.gz -r chr1 -R r12.txt --freq → afreq 含 chr1+chr2
```

`translate/filter.rs:463-473` 是 `if -r { ... } else if -R { ... }`；stats
两者都生效。修法：合并两个来源（或 clap `conflicts_with` 直接拒绝组合）。
补 `-r/-R/-t/-T` 组合矩阵测试。

### ⑨2 `-r/-R` 与 `-t/-T` 组合截断 `-t` —— 属实

**实测**：`filter -r chr2 -t chr1:10` → 输出 chr2:10（`-t` 应使其为空）。
`translate/filter.rs:340-354` 的目标集合只由一处构造，`-r` 时另一处被
覆盖/截断。最小修法：同时给出时 usage 错；完整修法：交集。

### ⑨3 `-o -` 在前缀命令创建 `-.xxx` 文件 —— 属实（实测）

```text
$ reed stats    ... --freq -o -        → exit 0, 生成 -.afreq
$ reed convert  ... --012  -o -        → exit 0, 生成 -.012/.012.pos/.012.indv/-.log
$ reed distance ... --type ibs -o -    → exit 0, 生成 -.mibs/-.mibs.id
```

`popg` 已拒绝；`stats/convert/pca/distance` 未拒绝。修法：统一
`prefix == "-"` → usage（或真正支持 stdout，但多文件输出没有单一 stdout
语义）。注意 `-o ''` 是另一个问题（⑨5）。

**拍板（2026-11-03）**：按推荐 A——前缀型命令统一拒绝 `-o -` 与空前缀
（`trim().is_empty()`），在打开输入前报 usage；`filter -o -` 的 stdout
语义保持不变。`popg` 的 `-o -` 检查同时上移到参数校验阶段。

### ⑨4 `popg -o -` 检查晚于计算 —— 代码复核（顺序属实）

`translate/popg.rs:144-155` 的检查在 run 之后；应上移到参数校验阶段
（约 61 行）。否则大输入先算完才报 usage。

### ⑨5 `distance -o ''` 未拒绝 —— 属实（实测）

```text
$ reed distance ... -o '' --type ibs
# exit 0，生成 .mibs/.mibs.id（隐藏文件）
```

`distance.rs:50-54` 缺 S19 的空前缀检查。修法同 stats：`trim().is_empty()`
→ usage。

### ⑨6 vs-noodles 基准 BCF 分支开销不可比 —— 基准方法（代码复核）

`bench/suites/vs-noodles/src/main.rs:50-53` 每样本 `series.name(&header)`
字符串查名（348×3.3M ≈ 11 亿次），VCF 分支一次 `select("GT")`。修基准后
重跑基线；不影响生产代码。

### ⑨7 `--format` 对表类 `--type` 静默忽略 —— 属实（低）

`distance.rs:48,283-292`：表类输出固定布局，`--format` 不生效也无提示。
加 usage 错或文档。

### ⑨8 `--contig-lengths` 错误缺文件名 —— 诊断（低）

`common.rs:357` 的错误只带 I/O 错误。补路径。

### ⑨9 tree.rs 内部节点 id 上限文档错 —— 文档（低）

`reed-graph/src/tree.rs:10` 写错；UPGMA 实际最大 id = 2n−1。

### ⑨10 bench 脚本顺序/吞 stderr —— 基准方法（低）

`bench/run_matrix.sh:6-10` 固定先 noo 后 hts、吞 stderr；首轮冷缓存偏差，
panic 留 `n_rec=0`。加 warmup/交替顺序/保留 stderr。

### ⑨ Lean：样本来源冲突逻辑两份实现 —— 精简（代码复核）

`filter.rs:264-317` 与 `common.rs:186-229` 约 40 行 ×2；doc_v21 S32 正是
这类漂移产生过 bug。建议合并到一处。

## 跨区域复核补充（报告未列/需修正）

1. **负 Float 符号**是 BCF→VCF 的数据损坏，不是 latent；应进 Must-fix
   前五。
2. **①1 的触发**需要额外空列；补测试时按"两个尾随 tab"构造。
3. **④1 驳回**，改语义会破坏 bcftools 兼容。
4. **⑧1 仅库 API**，CLI 不可达，优先级低于报告。
5. **⑧2 已实测** backend 间整体变号，不是"可能"。
6. **⑥2** 按单 contig 计，report 的"全 contig"表述需修正。
7. **⑤3** CLI 拒绝 pairs-only；库 API 才有 decay=None 的不一致。
8. **⑤2** 在 `--ld-decay + --ld-pairs` 组合下仍 CLI 可达。
9. **④2** 的同类检查：自定义 Number=A 标签已走 any-pair；`AC` 是唯一
   标量化 built-in；`N_ALT`/`AN`/`AF`/`MAC`/`MAF` 不是向量。
10. **NaN 排序**同类位置已排查：`umap.rs:95`、`auto_gamma`
    （`graph.rs:172`）、`extra.rs` 均有 finite/positive 过滤或
    `total_cmp`；仅 ⑧1 的三处未修。
11. **`-o -`** 实测覆盖 stats/convert/distance；`pca` 按同一前缀辅助函数
    推断同样受影响（未单独跑）。
12. **`format_g6_into`** 也被 `float_into` 用于 FORMAT 浮点渲染；测试应
    同时覆盖 INFO 与 FORMAT 两个入口。
