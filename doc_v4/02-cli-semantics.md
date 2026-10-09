# 02 参数语义统一规范

> 用户裁决（2026-10-09）：**参数语义要对齐一致；不必全部对齐主流软件。缺失阈值用 vcftools `--max-missing 0.8` 的方向；`--maf 0.05` 保持。**
>
> 本文是 doc_v2/02-cli.md 的**修订版**：doc_v2 的"名字相同则语义相同"原则保留，但在"方向"上加一条更高优先级的规则（§1）。
> 其余未提到的参数（`-r/-t/-s/-S/-v/-V/-m/-M/-f/-i/-e/-x/--fill-tags/-O/-o/-W/--pi/--dxy/…`）语义不变。

---

## 1. 方向规则【规】

1. **每个阈值参数都是一条"保留条件"**：下限参数 `X` 的含义是"保留满足 `量 ≥ X` 的对象"，上限参数 `max-X` 的含义是"保留满足 `量 ≤ X` 的对象"。
   不允许任何参数用"丢弃条件"（"量 > X 则丢弃"）表述。
2. **缺失率阈值的值是"必须被调用的比例"**（call rate）。`--max-missing 0.8` = 保留调用率 ≥ 0.8 的位点，即允许至多 20% 缺失。
   这是 vcftools 的约定；名字里的 "max-missing" 读作"对缺失的最大容忍，用保留比例表示"。
3. **同一概念在位点与样本两级使用同一方向、同一分子分母定义**（§3）。
4. 与规则 1–3 冲突的主流名字**不借用**（哪怕语义上有对应物）。违反者只有两个：plink2 的 `--geno`（丢弃方向）与 `--mind`（丢弃方向），本文删除。
5. 边界必须精确：阈值比较用**精确整数判定**（§4），不得通过 `1.0 − x` 之类的浮点换算改变边界。
6. 同一参数在 `reed filter` 与 `reed stats` 中语义、分母、求值阶段、漏斗 id 完全相同（已通过共用的 `SiteFilterArgs` 保证；本次改名同时改两处）。

---

## 2. 现状中的不一致（证据）

| 参数 | 现语义 | 方向 | 来源 | 位置 |
|---|---|---|---|---|
| `--geno FLOAT` | 位点缺失率 > FLOAT 则丢弃（`n_missing/n_kept ≤ m` 保留） | **丢弃** | plink2 | `reed-filter/src/eval/tail.rs:33-36`、`reed-cli/src/translate/filter.rs:116-117` |
| `--mind FLOAT` | 样本缺失率 > FLOAT 则丢弃（`miss/n_records > threshold` 剔除） | **丢弃** | plink2 | `reed/src/filter/mod.rs:256-268` |
| `--maf/--mac/--hwe/-q` | 量 ≥ X 保留 | 保留 | plink2/vcftools/bcftools | `tail.rs:38-47` |
| `--max-maf/--max-mac/-Q` | 量 ≤ X 保留 | 保留 | 同上 | 同上 |
| `--min-dp/--max-dp/--min-gq` | DP ∈ [min,max] 且 GQ ≥ min 的基因型保留（不在区间内的置缺失） | 保留 | plink2 | `translate/filter.rs:88-106` |

gold 管线自己就在两种方向之间换算：`run_vcftools.sh:18` 用 `1−m` 把 `SITE_MISS_MAX=0.10` 翻译成 `--max-missing 0.9`，`run_plink2.sh:29` 与 `run_reed.sh:25` 直接用 `--geno 0.1`。
这正是用户指出的问题：同一条规则两种说法。

**实测基线（2026-10-09，subset.Chr01-03，480,952 位点，348 样本，掩码后）**：

| 规则 | reed `--geno` | vcftools `--max-missing` | plink2 `--geno` |
|---|---|---|---|
| 调用率 ≥ 0.9（`--min-dp 5 --min-gq 20` 之后） | 263,420 | 263,420 | 263,420 |
| 调用率 ≥ 0.8 | 381,672 | 381,672 | — |

三者位点集合相同，所以改方向不改变任何现有结果，只改变参数的说法。改动后必须用同一组数字做回归门（§7）。

---

## 3. 每个阈值的精确定义【规】

记：`N` = 保留的样本数（经 `-s/-S/--keep/--remove` 与样本缺失规则之后）；`L` = 参与样本缺失计算的位点数；
`an` = 位点上保留样本的已调用等位基因数（半调用计 1，全调用计 2，`SiteCounts::an`）；`ac` = ALT 等位数；`n_missing_gt` = 完全缺失的基因型数。

### 3.1 缺失（位点级）

| 参数 | 保留条件 | 分子 / 分母 | 阶段 |
|---|---|---|---|
| `--max-missing X` | `an ≥ X · 2N`（精确整数判定见 §4） | 已调用等位 / 保留样本的等位总数 | 掩码之后、`--maf/--mac/--hwe` 之前（与 vcftools 的 `filter_sites_by_frequency_and_call_rate` 同序，`entry_filters.cpp:822`） |
| `--max-missing-count C` | `2N − an ≤ C` | 缺失等位数 | 同上（vcftools `entry_filters.cpp:918`：`(N_chr − N_non_missing_chr) > max` 则丢弃） |

- **按等位基因计**，不按基因型计：这是 vcftools 的定义（`N_non_missing_chr / get_N_chr()`），也是 `--max-missing-count` 唯一自洽的定义。
  半调用 `0/.` 贡献 1 个已调用等位。plink2 的 `--geno` 按基因型计；两者只在有半调用时不同。本数据实测 **0 个半调用**（Chr01:1–300 kb：351,980 已调用、2,604 缺失、0 半调用，VCF 与 BCF 一致），因此对现有结果没有影响。
- 现实现 `tail.rs:33-36` 用 `n_missing_gt / N`（基因型计）。改为等位计后，无半调用时 `2N − an = 2·n_missing_gt`，比例相同。

### 3.2 缺失（样本级）

| 参数 | 保留条件 | 分子 / 分母 | 阶段 |
|---|---|---|---|
| `--max-missing-sample X` | 样本的已调用等位数 `≥ X · 2L` | 该样本在 L 个位点上的已调用等位 / `2L` | 静态子集之后；在掩码之后统计；`L` = 通过掩码前位点规则（`-v/-V/-m/-M/-f/--snps-only/` 与不依赖基因型的 `-i/-e`）的记录数；在 `--max-missing/--maf/--mac/--hwe` 之前生效 |

- 阶段与分母沿用 V0 已落地的 plink2 口径（掩码后、头部存活记录；`doc_v3/STATUS.md` C0），只改方向与分子（等位计）。
- 样本级规则需要一遍全局归约（现有 pass0/sketch 机制不变）。
- `reed stats` 暂不提供此参数（需要多遍；现状已明确拒绝 `--mind`，改名后同样拒绝并提示用 `reed filter` 先过滤）。

### 3.3 频率与计数（不变）

| 参数 | 保留条件 | 分子 / 分母 |
|---|---|---|
| `--maf X` | `min(ac, an−ac) / an ≥ X` | 次等位 / 已调用等位（掩码后，保留样本） |
| `--max-maf X` | `min(ac, an−ac) / an ≤ X` | 同上 |
| `--mac C` / `--max-mac C` | `min(ac, an−ac) ≥ C` / `≤ C` | 次等位计数 |
| `-q X` / `-Q X` | 非参考等位频率 `≥ X` / `≤ X` | **INFO/AC 与 INFO/AN 同时存在时取 INFO 值**（bcftools 语义），否则 `ac/an`；阈值按 bcftools 量化为 f32 |
| `--hwe P [midp]` | 精确检验 `p ≥ P` | 掩码后、保留样本、双等位位点 |

`-q/-Q` 与 `--maf` 的区别必须在 `--help` 里写明：前者是"非参考等位频率、优先用 INFO"，后者是"次等位频率、总是按基因型重算"。本数据的 INFO/AC 是陈旧的；`-q/-Q` 取到 INFO 值时打印与表达式相同的"陈旧 INFO"警告（现只有表达式有这条警告，§6 F-C3）。

### 3.4 基因型掩码（不变）

| 参数 | 保留条件（否则该基因型置 `./.`） |
|---|---|
| `--min-dp A` / `--max-dp B` | `DP ∈ [A, B]`；DP 缺失不掩码 |
| `--min-gq G` | `GQ ≥ G`；GQ 缺失不掩码 |
| `--set-gts .` + `-i/-e 'FMT/…'` | bcftools `filter -S .` 语义 |

---

## 4. 精确整数判定【规】

阈值是十进制小数，`1 − X` 在二进制浮点中不精确：`1 − 0.8 = 0.19999999999999996`，若实现写成 `missing/N ≤ 1 − X`，`N = 10, missing = 2` 时 `0.2 ≤ 0.19999999999999996` 为假，边界位点被错误丢弃；而 vcftools 保留它（`8/10 ≥ 0.8` 为真）。
为避免这一类错误：

1. 阈值 `X` 解析后量化为有理数 `q / D`，`D = 10^9`，`q = round(X · D)`（覆盖用户会写的任何小数位数；plink2 也做类似的 `(1 + kSmallEpsilon)` 容差，`plink2_filter.cc:4220`）。
2. 位点规则：保留 iff `an · D ≥ q · 2N`（全部 u128 整数运算）。
3. 样本规则：保留 iff `called_alleles · D ≥ q · 2L`。
4. 计数规则：保留 iff `2N − an ≤ C`。

**与两家工具的等价性**：vcftools 保留 iff `an/(2N) ≥ X`（浮点比较，`entry_filters.cpp:820-823`），plink2 丢弃 iff `missing > floor(X' · N · (1+ε))`（`plink2_filter.cc:4220-4246`，`X' = 1 − X`）。
对整数 `N`、实数 `a`：`N − ⌈a⌉ = ⌊N − a⌋`，所以 "`called ≥ ⌈X·N⌉`" 与 "`missing ≤ ⌊(1−X)·N⌋`" 相同；整数判定 2 在 `X·2N` 恰为整数时取等号，与两家一致。
§7 的边界 fixture 必须覆盖 `X·2N` 为整数与非整数两种情况。

---

## 5. CLI 全表（修订后）

### 5.1 `reed filter` 与 `reed stats` 共用的过滤参数

| 参数 | 保留条件 | 阶段 | 变化 |
|---|---|---|---|
| `-r/-R`、`-t/-T` | 区域内 | 载入 | 不变 |
| `-s/-S/--keep/--remove`、`--force-samples` | 样本在集合内 | 载入 | 不变 |
| **`--max-missing-sample X`** | 样本调用率 ≥ X | 掩码后、位点阈值前 | **新；取代 `--mind`**（仅 filter） |
| `-v/-V/--snps-only/-m/-M/-f` | 类型/等位数/FILTER | 掩码前 | 不变 |
| `-i/-e EXPR` | bcftools 表达式 | 按是否依赖基因型分到掩码前/后 | 不变 |
| `--min-dp/--max-dp/--min-gq`、`--set-gts .` | 基因型在区间内（否则置缺失） | 掩码 | 不变 |
| `-q/-Q` | 非参考 AF 在区间内 | 掩码前（INFO 有值时）/后 | 不变；补"陈旧 INFO"警告 |
| **`--max-missing X`** | 位点调用率 ≥ X | 掩码后 | **新；取代 `--geno`** |
| **`--max-missing-count C`** | 位点缺失等位数 ≤ C | 掩码后 | **新**（vcftools） |
| `--maf/--max-maf/--mac/--max-mac` | 次等位频率/计数在区间内 | 掩码后 | 不变 |
| `--hwe P [midp]` | p ≥ P | 掩码后（最后） | 不变 |

### 5.2 删除的参数

| 参数 | 处理 |
|---|---|
| `--geno X` | 报错退出码 2：`--geno was removed: it stated the fraction to drop; use --max-missing {1−X} (the fraction that must be called)`，数值按 `1−X` 以十进制打印（例 `--geno 0.1` → `--max-missing 0.9`） |
| `--mind X` | 同上，提示 `--max-missing-sample {1−X}` |
| `--vcf-min-dp/--vcf-max-dp/--vcf-min-gq` | 保留（与主名同方向、同语义的别名，不违反 §1） |

不提供自动换算的隐藏别名：静默换算正是用户要消除的歧义来源。

### 5.3 漏斗 id

`after_mask.max_missing`、`after_mask.max_missing_count`、`sample.max_missing`（取代 `after_mask.missing`、`sample.mind`）。

### 5.4 `to_args()` 渲染

`FilterSpec::to_args()` 渲染 `--max-missing X`（spec 中存储的就是保留比例 `X`，不再存储 `max_missing: 1−X`）；`##reed_filterCommand` 随之改变。`spec.site_after_mask.max_missing` 字段**改名为 `min_call_rate`**，字段语义与参数一致，杜绝"字段叫 max_missing、值却是 drop rate"的旧状态（现 `TailRules.max_missing` 存的是 `--geno` 的值，`spec/mod.rs:142`）。

### 5.5 `reed stats` 独有参数（不变，列出以备核对方向）

| 参数 | 含义 | 方向性 |
|---|---|---|
| `--pop/--populations FILE` | 群体文件 | — |
| `--pi/--dxy/--theta/--tajima-d/--fst [wc84\|hudson]/--pbs A,B,C [--pbs-normed]` | 统计量 | — |
| `--freq/--hardy/--het/--missing/--depth/--summary` | 报告 | — |
| `--window/--step`、`--levels`、`--skip-empty` | 窗口 | — |
| `--variant-only`、`--bypass-invariant-check` | 输入声明（pixy 语义） | 声明，非阈值 |
| `--contig-lengths FILE`、`--force-samples` | 输入补充 | — |

`--missing` 报告的 `F_MISS` 列仍是"缺失比例"（plink2 `.smiss/.vmiss` 的列定义），它是输出数据而不是阈值参数，不受 §1 约束；帮助文本注明"报告列为缺失比例，阈值参数为调用比例"。

---

## 6. gold 与脚本

- `gold/scripts/run_reed.sh:25`：`--geno "$SITE_MISS_MAX"` → `--max-missing "$(awk -v m="$SITE_MISS_MAX" 'BEGIN{printf "%.10g", 1-m}')"`，与 `run_vcftools.sh:18` 相同的换算（`filters.conf` 继续用 `SITE_MISS_MAX`，它是三家工具的共同来源，不是 reed 的输入）。
- `SAMPLE_MISS_MAX` 在 gold 中用于生成 `samples.keep.txt`，reed 继续用 `-S`，不受影响。
- `gold/popgen` 不涉及这些参数。

---

## 7. 测试门【规】

1. **边界 fixture**（新增 `reed/tests/e2e_missing_boundary.rs`）：
   `N = 10` 与 `N = 7`（`X·2N` 分别为整数与非整数）、`X ∈ {0.8, 0.9, 0.5, 1.0, 0.0}`、含半调用与不含半调用；
   reed `--max-missing X` 的位点集合 == vcftools `--max-missing X`；无半调用的用例同时 == plink2 `--geno 1−X`（含 plink2 的 ε 规则）。
2. **样本级**：fixture 上 reed `--min-gq 20 --max-missing-sample 0.9` 保留的样本集合 == plink2 `--vcf-min-gq 20 --mind 0.1`（V0 的 fixture 改参数后复用）。
3. **回归**：subset 上 `--min-dp 5 --min-gq 20 --max-missing 0.9` → 263,420 位点；`0.8` → 381,672（§2 的实测基线）。
4. **gold**：`run_reed.sh` 改参数后，全量正文 md5 `7aedc59ac798…` 不变，`check_parity.py` PASS。
5. **往返**：`to_args()` 渲染 `--max-missing`，`spec_from_flags(to_args(spec)) == spec`；旧参数 `--geno/--mind` 的报错文本有测试。
6. **精确判定**：单元测试构造 `X = 0.8, N = 10, missing = 2`（等号成立）与 `X = 0.7, N = 10, missing = 3`（`1−0.7 = 0.30000000000000004`）两例，断言保留。

---

## 8. 与 doc_v2/02-cli.md 的差异登记

| doc_v2 条款 | 本文 |
|---|---|
| §3.1 `--mind`（plink2，丢弃方向） | 删除，改 `--max-missing-sample`（保留方向，等位计） |
| §3.4 `--geno`（plink2） | 删除，改 `--max-missing`（vcftools）+ `--max-missing-count` |
| §0 原则 2"名字相同则语义相同" | 保留，但从属于本文 §1 的方向规则 |
| 05-naming.md §2.3"`--max-missing` 方向与 `--geno` 相反，不提供别名" | 反转：`--max-missing` 成为主名，`--geno` 删除 |
