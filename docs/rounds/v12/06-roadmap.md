# 06 修复路线（R0–R3）

> 原则：**先把每条缺陷钉成可跑的红测试（夹具全部来自 01-review 的实测，
> 可直接重建），再修 P0+P1（一次 PR 一个缺陷族），然后 P2 收口，最后
> lean/文档/守卫收尾**。修复对象是 01-review.md §1 的 D12-1…D12-15；
> 规范依据见 02-architecture.md 的 N1…N14。不重开 v6 的 Y 系列、v7–v10 的
> 结构项、v11 的门面项。每步满足 §1 通用门。

---

## 1. 通用验收门（每一步都必须满足）【规】

1. **双 profile 全绿**：`cargo test --workspace`（debug）与
   `cargo test --release --workspace` 都通过；
   `cargo clippy --release --workspace --all-targets -- -D warnings`
   不新增警告。
2. **红 → 绿**：每个 P0/P1 缺陷先提交复现它的失败测试（测试名描述缺陷
   本身），再提交修复转绿；红测试在修复 commit 前必须实测红过。
3. **htslib 真值门（限 §2.5-A 的项）**：凡改变 BCF 写出字节的步骤（R1 的
   D12-3/4/5、R2 的 D12-9），用工作区 `bcftools view`（1.24）读回写出物：
   退出码 0、INFO 列与 bcftools 自产 BCF 的读回文本一致。
   **例外**：01-review §2.5-B 的"不改"项（B1 文本保真、B2 Σac>an、
   B3 前缀解析、B4 空格容忍、B5 Number=1 多值）**不适用本门**——它们是
   永久分歧，对拍断言必须排除这些形态，否则测试会把分歧钉成契约。
   夹具里出现这些形态时，断言 reed 自己的规范值（02 §13.3 表的 reed 列），
   不断言与 bcftools 相同。
4. **确定性**：threads ∈ {1, 4} 输出 sha256 相同。
5. **性能**：不回退超过 10%；INFO 解析/编码是每记录热路径（N3/N4/N13 都
   触碰），`bench/suites/e2e/results.tsv` 追加一行；348 样本 gold 输入的
   `-Ob` 墙钟时间与本轮基线（v10 STATUS：字节一致、时间不变）对比。
6. **三策略等价**：触碰 tail/策略的步骤（R1-a、R1-f）跑
   `all_three_strategies_produce_identical_bytes` 的扩展夹具（含尾部表达式
   与 `-q`）；SinglePass ≡ TwoPassReparse ≡ SketchThenEmit（或按 N1 报错）。
7. **同一 PR 内删除被替代的旧实现**，不留并行代码（v2 06 §1 纪律 7 沿用）：
   N12/N13 收口时旧的重复实现同 PR 删除。
8. **"不改"项的落地形式是文档**：§2.5-B 的 B1–B5 必须写进 `cli.md §5` 差异表
   （每条注明"reed 有意为之 + 理由"），这是它们的唯一交付物——**不写代码、
   不写对拍测试**（守卫 15）。B6（f32 量化）是护栏项：已对齐，加一条
   "不得改成 f64"的测试注释即可。
9. 本目录相应条目标记完成并注明 commit；`docs/rounds/v12/STATUS.md` 同步
   更新；落地一步即在 `ci/check-doc-drift.sh` 登记 01-review §3 的对应守卫。

---

## 2. 步骤

### R0 红测试先行（D12-1/2/3/4/5/6/7/8/9/10/13/14 的复现夹具）— 门的基础

> R0 不修任何代码。夹具全部是 01-review 实测用过的最小 VCF（内联字符串
> 即可，无需外部文件）；每个红测试的断言是"与 bcftools 1.24 的行为一致"
> （**仅 §2.5-A 的项**）、"三策略一致"、或"reed 三条路径自洽"
> （D12-14 的内部一致性，与 bcftools 无关）。R0 提交里这些测试全红
> （或对 D12-8 用探针断言红）。
>
> **§2.5-B 的"不改"项（B1–B8）在 R0 没有红测试**——它们的交付物是
> `cli.md §5` 差异表（通用门 8、守卫 15），落在 R3-b。若在 R0 给它们写
> 对拍断言，等于把永久分歧钉成契约，下一轮无法再改。

| 项 | 内容 | 门 |
|---|---|---|
| R0-a | D12-1：`mind_with_tail_expression_keeps_records`（4 记录夹具 + `--mind 0.9` + `-i 'F_MISSING<=0.5 && QD>10'`，断言 2 条且与 single_pass 逐字节一致）。**修复前：0 条** | `cargo test -p reed --test e2e_filter` |
| R0-b | D12-2：`gt_less_records_keep_their_ac_an_under_subset`（nogt 夹具，断言输出 `AC=4;AN=6`）。**修复前：`AC=0;AN=0`** | 同上 |
| R0-c | D12-3/4/5：`bcf_value_fidelity`三例（`XD=3000000000` → `XD=.`；`DB=1`/`DB=0` → 字面往返；`SOMESTR=.` → `SOMESTR=.`），bcftools 读回断言。**修复前：回绕/丢值/裸标签** | `cargo test -p reed --test e2e_bcf` + 真值门脚本 |
| R0-d | D12-6：`min_af_after_subset_uses_the_genotypes`（sub 夹具 `-s s2,s3 -q 0.1`，断言 0 条）。**修复前：2 条** | e2e_filter |
| R0-e | D12-7：`expression_flag_vector_bare_parity`（flag3/vec/zero 三夹具 × bcftools 期望值表）。**修复前：flag/向量/裸标签三组全红** | e2e_filter |
| R0-f | D12-8：`quantifier_any_requires_a_satisfying_element`（reed-filter 单测，`AC=1,2` × `>= 10`/`<= 0` × any/all 四组合）。**修复前：any 两例红** | `cargo test -p reed-filter` |
| R0-g | D12-9：`info_count_overflow_errors_instead_of_corrupting`（65 536 片段 / 300 FORMAT 键 / 66 001 等位三夹具，断言 Err）。**修复前：Ok 且文件坏** | e2e_bcf |
| R0-h | D12-10：`invariant_sites_count_zero_alts`（inv 夹具 × `-m 2`/`-M 1`/`N_ALT=0`/`N_ALT=1` 期望表）。**修复前：四例中三例红** | e2e_filter |
| R0-i | D12-13：`ac_length_mismatch_falls_back_to_genotypes`（aclen3 夹具 `-q 0.1`，断言 0 条）。**修复前：1 条** | e2e_filter |
| R0-j | D12-14：`info_value_is_the_same_for_filtering_and_writing`（31 形态夹具；断言"表达式判定为真"的记录其 `-Ob` 写出值满足同一比较 + 逐行等于 02 §13.3 表的 **reed 列**）。**修复前：`XI= 5`/`1e3`/`5abc` 三例判定为真而写出 `.`** | e2e_bcf；**断言 reed 自洽，不对拍 bcftools**（B3/B4 形态属"不改"，通用门 3 例外） |

### R1 P0+P1 修复（D12-1…D12-8）— 主量

> 族内互不依赖，可并行分支；每族一个 PR，红绿对照在自己的 commit 里可验证
> （revert 修复即红）。D12-1（策略）与 D12-6（`-q` cohort 门）有交叠
> （N1 §4 的判据引用 N5），建议 D12-6 先行或同 PR。

| 项 | 内容 | 规范 | 门 |
|---|---|---|---|
| R1-a | D12-1：驱动策略覆盖查 `has_tail_expr`（+ `-q` INFO 判据）；`--mind --indep-pairwise` + 尾部表达式报 usage；R0-a 转绿 | 02 §1 | 通用门 6；守卫 #1 |
| R1-b | D12-2：`fill_info_set` + emit stale 过滤接 `has_gt`；R0-b 转绿 | 02 §2 | 守卫 #2 |
| R1-c | D12-3：`parse_int_vector`/`encode_info_typed` 值域 → missing + 一次性警告；R0-c（整数组）转绿 | 02 §3 | 真值门；守卫 #3 |
| R1-d | D12-4/5：值 × 类型编码收口单一函数（Flag 带值 → 字符串、String `.` 字面存）；R0-c（flag/string 组）转绿 | 02 §4 | 真值门；守卫 #4 |
| R1-e | D12-6：`cohort_changed()` 门（head 不推 INFO 规则 + `af_of` 跳 INFO）+ `-q/-Q` 陈旧警告；R0-d 转绿 | 02 §5 | 守卫 #5；`e2e_af_type` 全绿 |
| R1-f | D12-7：`check`/求值接 `TagTable`；Flag=1/0、向量 any（Number=1 而记录多值 → 整个 missing）、裸标签编译期拒绝（【决】默认）；R0-e 转绿。**与 R2-f 同改 `expr/mod.rs` 的 `info_value`，建议同 PR 或紧邻**（R1-f 加类型分派，R2-f 换解析函数） | 02 §6 | 守卫 #6；bcftools 对拍表全绿 |
| R1-g | D12-8：`any |= cmp.apply_f64(v)`；R0-f 转绿 | 02 §7 | 守卫 #7 |

### R2 P2 收口（D12-9…D12-14）

| 项 | 内容 | 规范 | 门 |
|---|---|---|---|
| R2-a | D12-9：打包计数上限 helper，四写回点接入；R0-g 转绿 | 02 §8 | 守卫 #8 |
| R2-b | D12-10：`count_alts` 单点化（head/expr/SiteView 共用）+ `max_alt_alleles: Option<u32>`；R0-h 转绿 | 02 §9 | 守卫 #9；spec 往返测试 |
| R2-c | D12-11：`--fill-tags` 校验（支持集常量 + 大小写归一 + usage 错）；cli.md:117 回改 | 02 §10 | 守卫 #10；gold 输出不变 |
| R2-d | D12-12：`-x` help/README/cli.md 可见化（含 `-x ''` 与 `-I` 叠加说明） | 02 §11 | 守卫 #11 |
| R2-e | D12-13：AC/AN 读取单点化 + 长度校验 + 一次性警告；summary 第三份实现删除；R0-i 转绿 | 02 §12 | 守卫 #12；stats 测试全绿 |
| R2-f | D12-14：`parse_info_number` 单点化（编码器/表达式/head/summary 共用；删 `as_num` 的 `trim()`；表达式不再退回字节序比较；整数模式**接受可选前置 `+`**——规范 Float 正则允许 `[-+]?`，tex:97，属合法输入非畸形）；cli.md §5 解析差异表（02 §13.3）；R0-j 转绿 | 02 §13 | 守卫 #13；性能门（§1.5）；内部一致性红测试 |

### R3 P3 + 文档与守卫收尾（D12-15）

| 项 | 内容 | 规范 | 门 |
|---|---|---|---|
| R3-a | D12-15：删 `edit_vcf_line` 死分支，签名收不可失败（`single_pass.rs:378`、`emit.rs:165`、`filtered.rs:177` 三调用点同步） | 02 §14 | 守卫 #14；clippy |
| R3-b | 文档收口：02 §15 清单逐项回改 + **`cli.md §5` 差异表登记 §2.5-B 的 B1–B5**（每条注明"reed 有意为之 + 理由"；这是"不改"项的唯一交付物）+ CHANGELOG（`-q` cohort 门、`--fill-tags` 校验、默认 `-x` 可见化、INFO 解析差异表为显著条目）+ drift 守卫全量登记 + STATUS 收口 | 02 §15；01 §2.5-B | `bash ci/check-doc-drift.sh` 全绿；守卫 #15（B1–B5 五条均在差异表内） |

---

## 3. 顺序与依赖

```
R0（红测试） → R1（P0+P1 七连） → R2（P2 六连） → R3（lean/文档）
```

- R1 内部：R1-e（`-q` cohort 门）先于或与 R1-a 同 PR（N1 §4 的判据引用
  N5 的门）；其余互不依赖。
- R1-f 与 R2-f **同改 `expr/mod.rs` 的 `info_value`**（前者加类型分派、后者
  换解析函数并删字节序退回）——建议同一 PR 或紧邻落地，否则后落的一方要
  重做前一方刚改完的函数。
- R2 内部：R2-e（AC/AN 单点化）与 R1-e 触碰同一函数族（`af_of`/
  `allele_counts_info`），建议相邻落地；R2-b 的 `count_alts` 是 R2-e 共享
  实现的入参来源（`n_alt`），先行。
- R2-f（解析单点化）依赖 R1-c/R1-d（编码值处理已收口）与 R2-e（AC/AN
  扫描已收口），排在 R2 末尾。
- **若排期紧张**：R1-a/R1-b（两个 P0，静默错结果/数据损坏）必须最先；
  R1-c/d（BCF 值忠实性）次之（真值门成本低）；R2-d/R3-b（文档）可并入
  任一相邻 PR。

## 4. 本轮不做什么

> 分两类：**「标注不改」= 01-review §2.5-B 的永久分歧**（reed 更优，
> 交付物是 `cli.md §5` 差异表，永不写代码/对拍测试）；**「本轮不做」=
> 范围外或【决】项**（可能后续立项）。

### 4.1 标注不改（§2.5-B，永久分歧）

| 项 | reed 行为 | bcftools/htslib | 不改的理由 |
|---|---|---|---|
| B1 文本输出归一化 | `-Ov`/`-Oz` 对 INFO 列**逐字节保真** | 归一化（`+5`→`5`） | §6.7 设计；保真严格强于归一化，bcftools 的归一化是 record 往返的副作用 |
| B2 `Σac > an` | 丢该位点，继续跑 | **`exit(1)`** 进程退出（`vcfutils.c:85-90`） | 不因一条畸形记录毁掉数小时的批处理 |
| B3 Integer 前缀解析 | 严格十进制，畸形 → missing（**例外**：`+5` 属合法输入——规范 Float 正则允许 `[-+]?`，N13 整数模式顺手接受） | `hts_str2int` 前缀解析（`5abc`→5、`1e3`→1、`0x10`→0） | htslib 自己不自洽（`0x10` 在 Integer 得 0、Float 得 16）；`5.0`/`1e3` 在 Integer 标签不匹配规范任何数值文法（§1.3），拒绝畸形更可预测 |
| B4 数值中的空格 | 三条路径统一**不 trim** | Integer 不容忍、Float 容忍前空格 | **理由更正（对照 VCFv4.5）**：规范 4.5 起允许**值**含空格（§1.6.1 tex:431）——初版"规范不允许空格"是错的；但数值文法（Float 正则 tex:97、Integer"32-bit signed"）首尾无空白，` 5` 属规范外。htslib 两种类型不一致无设计意图。结论不变：不 trim、missing |
| B5 Number=1 却给多值 | 当 **missing**（【决】C2 默认） | 标量读作 **0** | **规范支持 reed**：§1.4.2（tex:158）Number=1 就是"single number"，多值属规范外输入；htslib 那个 0 是实现细节（`vcf.c:4286` 先置 0、`len != 1` 不填），无注释表明有意。missing 是最诚实的恢复 |
| B6 `-q/-Q` 的 f32 量化 | 刻意复现 `sscanf("%f")`（`plan/mod.rs:796-798`） | C 的 `float` 存储（无意为之） | **护栏项**：已对齐，改成 f64 会在 `AC=34,AN=680` vs `0.05` 边界翻转结果、破坏 gold 对拍 |
| B7 表达式语言子集 | 下标/函数/算术报 "not supported yet" | 完整实现 | `cli.md §5` 纪律：超出子集要报告，绝不用另一种语义静默执行 |
| B8 sites-only VCF（8 列） | 拒绝 + 明确报错 | 接受 | 格式级决定（`LineRef` 9-tab 契约牵连扫描器/编辑器/编码器/index 全链）；独立立项候选，见 §4.2 |

### 4.2 范围外 / 待裁决

| 项 | 原因 |
|---|---|
| sites-only VCF（8 列）支持（= B8 的立项） | 独立立项候选，不塞进缺陷修复轮 |
| 表达式语言扩容（下标 `[0]/[*]`、`MAX/MIN/…` 函数、算术）（= B7 的立项） | `cli.md §5` 子集纪律；N6 只收三类**已接受语法**的语义 |
| `--fill-tags AF/F_MISSING` 的实现 | 【决】C3（02 §10.3）：本轮最小修法是校验 + 文档收敛 |
| 默认 `drop_stale` 改 opt-in | 【决】C4（02 §11.4）：默认清单是 gold 对拍设计；本轮只补可见性 |
| 数值解析完全对齐 htslib 前缀语义（= 推翻 B3/B4） | 【决】C5（02 §13.4）：改变现有输入结果、热路径变慢、gold 无此形态 |
| `bcf_record_check` 全量移植 | v10 §4 裁决沿用：只补"会 panic/会损坏/会静默错"的最小集（本轮 N8 即其一） |
| 读侧 PASS 注入（`BcfDict`） | v10 §4 遗留【荐】项，与 INFO 无关，不并轮 |
| `reed stats` 的 INFO 直方图（architecture.md §10 `info.<TAG>`） | 未实现的功能面，非缺陷；M 系列里程碑事项 |
