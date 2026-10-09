# 02 参数语义统一规范（2026-10-09 晚：命名回退版）

> 用户裁决（2026-10-09，两条按时间先后）：
> **第一次**：参数语义要对齐一致；不必全部对齐主流软件；`--maf 0.05` 保持。
> **第二次（本文 §5/§8 据它重写）**：缺失阈值只有**位点**与**样本**两级，因此 `--geno`/`--mind` **保留原样**——短名字、plink2 拼写、值 = 允许的缺失比例。撤回 W0 引入的 `--max-missing`、`--max-missing-sample`、`--max-missing-count`。
>
> 本文是 docs/rounds/v2/02-cli.md 的**修订版**；doc_v4 首版（已被 W0 实现）的"改名"部分由本文取代，**但 W0 在实现层的修正全部保留**（精确有理数判定、等位计分子、掩码后求值、头部存活记录作分母）——那些改动用户不可见，且已与 vcftools/plink2 逐位对拍。
> 其余未提到的参数（`-r/-t/-s/-S/-v/-V/-m/-M/-f/-i/-e/-x/--fill-tags/-O/-o/-W/--pi/--dxy/…`）语义不变。

---

## 1. 方向规则【规】

1. **每个阈值参数都是一条"保留条件"**，作用的量分两类：
   - 要**保证**的量（MAF、MAC、HWE p、调用率）：下限参数写作 `≥`；
   - 要**抑制**的量（缺失率、最大 MAF、最大 MAC）：上限参数写作 `≤`。
   不允许"量 > X 则丢弃"这类以动作为中心的表述，也不允许同一个量在位点与样本两级换方向。
2. **缺失率的阈值作用在"缺失率"（要抑制的量）上，保留 plink2 的短名字**：
   `--geno X` = 保留缺失率 ≤ X 的位点（等价地：调用率 ≥ 1−X）；`--mind X` = 保留缺失率 ≤ X 的样本。
   `--geno 0.2` 与 vcftools 的 `--max-missing 0.8` 是同一条规则。
3. **同一概念在位点与样本两级使用同一方向、同一分子分母定义**（§3）。
4. **CLI 路径上不存在 `1−x` 换算**：用户写的数就是存储与比较的数（§4）。唯一出现减法的地方是 §5.2 的撤回参数报错文本，且它用精确十进制渲染。
5. 边界必须精确：比较用**精确整数判定**（§4）。
6. 同一参数在 `reed filter` 与 `reed stats` 中语义、分母、求值阶段、漏斗 id 完全相同（由共用的 `SiteFilterArgs` 保证）。

### 1.1 为什么这里不追求"名字统一"

第一次裁决的理由（"方向要一致"）仍然成立，且**已被规则 1 满足**：`--geno`/`--mind` 是"对缺失率的上限"，与 `--max-maf`/`--max-mac`/`-Q` 同属 `≤` 家族，方向并不冲突——首版把它误判成"方向例外"，才引出改名。真正的取舍只剩**名字**：

- 改名（首版方案）能得到方向统一，代价是 `--max-missing`（13 字符）与 `--max-missing-sample`（20 字符）两个长名字，且与 plink2 的肌肉记忆相反；
- 不改名保住了短名字与 plink2/plink 1.9 的既有习惯，位点/样本两级仍然自洽（同为缺失率、同为 `≤`、同为等位计、同取掩码后）。

用户第二次裁决选后者。本文据此记录：**缺失阈值这一族只有 `--geno`、`--mind` 两个参数，不多造名字。**

---

## 2. 现状与证据

### 2.1 首版（W0）做过的、本次保留的实现层修正

| 修正 | 内容 | 为何保留 |
|---|---|---|
| 精确有理数判定 | 阈值量化为 `q/10^9`，全整数比较 | 与首版无关的技术改进；杜绝边界浮点误判 |
| 等位计分子 | 位点级用"已调用等位数 / 2N"；样本级用"样本已调用等位数 / 2L"（半调用计 1，掩码计 0） | vcftools 的定义；`doc_v4` 首版实测本数据 0 个半调用，故与基因型计结果相同 |
| 掩码后求值 | 位点规则在基因型掩码之后、`--maf/--mac/--hwe` 之前 | 与 vcftools 的 `filter_sites_by_frequency_and_call_rate` 同序（`entry_filters.cpp:822`） |
| 分母 = 头部存活记录 | 样本级规则的 `L` 只数通过 `-v/-V/-m/-M/-f/--snps-only/` 与不依赖基因型的 `-i/-e` 的记录（V0 已定，`docs/rounds/v3/STATUS.md` C0） | plink2 口径，已对拍 |
| 报告形态 | `sample_rule: {rule, threshold, n_dropped}` | W0 新增，与 funnel id 并存（见 §5.3 的登记） |

### 2.2 本次撤回的只有三件事

| 撤回 | 原（W0） | 现 |
|---|---|---|
| 位点级参数名 | `--max-missing X`（保留方向） | `--geno X`（缺失方向，`== --max-missing 1−X`） |
| 样本级参数名 | `--max-missing-sample X` | `--mind X` |
| 计数形式 | `--max-missing-count C` | **删除**（见 §5.2；本族定为两个参数） |

### 2.3 实测基线（2026-10-09，subset.Chr01-03，480,952 位点，348 样本，掩码后）

| 规则 | reed | vcftools | plink2 |
|---|---|---|---|
| 缺失率 ≤ 0.1（`--min-dp 5 --min-gq 20` 之后） | `--geno 0.1` → 263,420 | `--max-missing 0.9` → 263,420 | `--geno 0.1` → 263,420 |
| 缺失率 ≤ 0.2 | `--geno 0.2` → 381,672 | `--max-missing 0.8` → 381,672 | — |

三家集合逐个相同；本表同时是 W0b 的回归门（§7.3）。**回退改名不改变任何数字**，只把说法换回短名字——gold 的 `SITE_MISS_MAX=0.10` 现在可以原样传给 `--geno`，不必再换算。

---

## 3. 每个阈值的精确定义【规】

记：`N` = 保留的样本数（经 `-s/-S/--keep/--remove` 与样本缺失规则之后）；`L` = 参与样本缺失计算的位点数；
`an` = 位点上保留样本的已调用等位基因数（半调用计 1，全调用计 2，`SiteCounts::an`）；`ac` = ALT 等位数。
**缺失等位数** = `2N − an`（位点级）、`2L − called`（样本级）。

### 3.1 缺失（位点级）

| 参数 | 保留条件 | 分子 / 分母 | 阶段 |
|---|---|---|---|
| `--geno X` | `2N − an ≤ X · 2N`（精确整数判定见 §4） | 缺失等位 / 保留样本的等位总数 | 掩码之后、`--maf/--mac/--hwe` 之前（vcftools `entry_filters.cpp:822` 同序） |

- 等价写法：`an ≥ (1−X) · 2N`（§4 给出边界等价性证明）。
- **按等位基因计**，不按基因型计：vcftools `--max-missing` 的 `N_non_missing_chr / get_N_chr()` 即此；半调用 `0/.` 贡献 1 个已调用等位。plink2 的 `--geno` 按基因型计，两者只在有半调用时不同；本数据实测 **0 个半调用**（Chr01:1–300 kb：351,980 已调用、2,604 缺失），VCF 与 BCF 一致，故与两家都逐位相同。
- 计数形式（首版的 `--max-missing-count C`）**不提供**：`--geno X` 配上精确判定已能表达任意阈值（把 `C/2N` 写成小数即可，2N 由头部给出，比较是整数运算）。将来若确有 `C` 的需求，加 `--geno-count INT`，不要复活 `--max-missing-count`。

### 3.2 缺失（样本级）

| 参数 | 保留条件 | 分子 / 分母 | 阶段 |
|---|---|---|---|
| `--mind X` | `2L − called_s ≤ X · 2L` | 该样本在 L 个位点上的缺失等位 / `2L` | 静态子集之后；掩码之后统计；`L` = 头部存活记录数；在位点级阈值之前生效 |

- 只属 `reed filter`。样本级规则需要一遍全局归约（现有 pass0/sketch 机制不变）。
- `reed stats` 不提供：需要多遍（现状即拒绝 `--mind`，报错文本改为"先 `reed filter --mind X`"）。

### 3.3 频率与计数（不变）

| 参数 | 保留条件 | 分子 / 分母 |
|---|---|---|
| `--maf X` | `min(ac, an−ac) / an ≥ X` | 次等位 / 已调用等位（掩码后，保留样本） |
| `--max-maf X` | `min(ac, an−ac) / an ≤ X` | 同上 |
| `--mac C` / `--max-mac C` | `min(ac, an−ac) ≥ C` / `≤ C` | 次等位计数 |
| `-q X` / `-Q X` | 非参考等位频率 `≥ X` / `≤ X` | **INFO/AC 与 INFO/AN 同时存在时取 INFO 值**（bcftools 语义），否则 `ac/an`；阈值按 bcftools 量化为 f32 |
| `--hwe P [midp]` | 精确检验 `p ≥ P` | 掩码后、保留样本、双等位位点 |

`-q/-Q` 与 `--maf` 的区别必须在 `--help` 里写明：前者是"非参考等位频率、优先用 INFO"，后者是"次等位频率、总是按基因型重算"。本数据的 INFO/AC 是陈旧的；`-q/-Q` 取到 INFO 值时打印与表达式相同的"陈旧 INFO"警告。

### 3.4 基因型掩码（不变）

| 参数 | 保留条件（否则该基因型置 `./.`） |
|---|---|
| `--min-dp A` / `--max-dp B` | `DP ∈ [A, B]`；DP 缺失不掩码 |
| `--min-gq G` | `GQ ≥ G`；GQ 缺失不掩码 |
| `--set-gts .` + `-i/-e 'FMT/…'` | bcftools `filter -S .` 语义 |

---

## 4. 精确整数判定【规】

1. 阈值 `X` 解析后量化为有理数 `q / D`，`D = 10^9`，`q = round(X · D)`（覆盖用户会写的任何小数位数；plink2 用 `(1 + kSmallEpsilon)` 达到类似目的，`plink2_filter.cc:4220`）。
2. 位点规则：保留 iff `(2N − an) · D ≤ q · 2N`（全 u128 整数运算）。
3. 样本规则：保留 iff `(2L − called) · D ≤ q · 2L`。

**为什么这里不出浮点陷阱**：首版的陷阱来自 `1 − X`（`1 − 0.8 = 0.19999999999999996`）。回退后用户直接写缺失比例，CLI 路径上没有任何减法，陷阱不可达；精确判定仍然保留，作为"用户写 0.07 就是 0.07"的保证。

**与两家工具的等价性**：vcftools 保留 iff `an/(2N) ≥ 1−X`（浮点比较，`entry_filters.cpp:820-823`），plink2 丢弃 iff `missing > floor(X · N · (1+ε))`（`plink2_filter.cc:4220-4246`）。
对整数 `N`、实数 `a` 有 `N − ⌈a⌉ = ⌊N − a⌋`，故 "`missing ≤ ⌊X·N⌋`"、"`called ≥ ⌈(1−X)·N⌉`" 与整数判定 2 在 `X·2N` 恰为整数时**同时取等号**，三家一致。§7.1 的边界 fixture 必须覆盖 `X·2N` 为整数与非整数两种情况。

---

## 5. CLI 全表（回退后）

### 5.1 `reed filter` 与 `reed stats` 共用的过滤参数

| 参数 | 保留条件 | 阶段 | 变化 |
|---|---|---|---|
| `-r/-R`、`-t/-T` | 区域内 | 载入 | 不变 |
| `-s/-S/--keep/--remove`、`--force-samples` | 样本在集合内 | 载入 | 不变 |
| **`--mind X`** | 样本缺失率 ≤ X | 掩码后、位点阈值前 | **恢复原名**（仅 filter） |
| `-v/-V/--snps-only/-m/-M/-f` | 类型/等位数/FILTER | 掩码前 | 不变 |
| `-i/-e EXPR` | bcftools 表达式 | 按是否依赖基因型分到掩码前/后 | 不变 |
| `--min-dp/--max-dp/--min-gq`、`--set-gts .` | 基因型在区间内（否则置缺失） | 掩码 | 不变 |
| `-q/-Q` | 非参考 AF 在区间内 | 掩码前（INFO 有值时）/后 | 不变；补"陈旧 INFO"警告 |
| **`--geno X`** | 位点缺失率 ≤ X | 掩码后 | **恢复原名** |
| `--maf/--max-maf/--mac/--max-mac` | 次等位频率/计数在区间内 | 掩码后 | 不变 |
| `--hwe P [midp]` | p ≥ P | 掩码后（最后） | 不变 |

### 5.2 撤回的名字

| 撤回的参数 | 处理 |
|---|---|
| `--max-missing X` | 解析期报错，退出码 2：`--max-missing was removed: it stated the fraction that must be called; use --geno {1−X} (the fraction that may be missing)`，`1−X` 以**精确十进制**打印（`--max-missing 0.9` → `--geno 0.1`，不是 `0.09999999999999998`） |
| `--max-missing-sample X` | 同上，提示 `--mind {1−X}` |
| `--max-missing-count C` | 报错：无同义参数；提示用 `--geno` 表达等价的缺失比例（`C/2N`），理由见 §3.1 |
| `--vcf-min-dp/--vcf-max-dp/--vcf-min-gq` | 保留（与主名同语义的别名） |

不提供自动换算的隐藏别名：静默换算正是要消除的歧义来源。这三条报错文本同时是"没有 keep 方向拼写"的守卫（§7.5）。

### 5.3 漏斗 id 与报告

`after_mask.missing`（`--geno`）、`sample.mind`（`--mind`）。
W0 引入的 `after_mask.max_missing`、`after_mask.max_missing_count`、`sample.max_missing` **作废**。

样本级规则在 `FilterReport` 中同时有 `sample_rule: {rule, threshold, n_dropped}`（W0 新增的形态，保留；`rule` 文本改为 `mind`，`threshold` 仍是用户写的 X）。这是 W0 已登记的偏差（`docs/rounds/v3/STATUS.md` 风格的偏差登记见 `docs/rounds/v4/STATUS.md` "未完成/偏差"）：报告只有 head/tail/mask 的记录计数，样本级因此不落在 `funnel` 里；`sample.mind` 作为条目名在 `--verbosity 1` 的漏斗表里指向这一行。

### 5.4 `to_args()` 与 spec 字段

`FilterSpec::to_args()` 渲染 `--geno X` / `--mind X`；spec 里存的就是用户写的缺失比例 `X`，与参数值一一对应（首版"存保留比例、再由渲染取补"的绕行取消，渲染不再做算术）。`##reed_filterCommand` 随之回到 `--geno 0.1 … --mind 0.2`。

| 处 | W0 | 现在 |
|---|---|---|
| CLI | `--max-missing` / `--max-missing-sample` | `--geno` / `--mind` |
| `TailRules.max_missing` | `min_call_rate`（保留比例） | **`max_missing`（缺失比例，即 `--geno` 的值）** |
| 计数规则字段 | 新增 `max_missing_alleles` | 删除 |
| `FilterSpec.min_sample_call_rate` | 由 `mind` 改名而来 | **`mind`** |
| 阈值类型 | `reed_filter::spec::CallRate`（`keeps(called, total)`） | 同名类型或改名 `MissingRate`，方法改 `tolerates(missing, total)`；`complement()`/`render()`（精确十进制）保留，供 §5.2 的报错文本用 |
| 漏斗 id | `after_mask.max_missing` / `sample.max_missing` | `after_mask.missing` / `sample.mind` |

四处（CLI 名、字段名、漏斗 id、`to_args()`）必须同向，往返测试锁定（§7.5）。

### 5.5 `reed stats` 独有参数（不变，列出以备核对方向）

| 参数 | 含义 | 方向性 |
|---|---|---|
| `--pop/--populations FILE` | 群体文件 | — |
| `--pi/--dxy/--theta/--tajima-d/--fst [wc84\|hudson] [--fst-global]/--pbs A,B,C [--pbs-normed]` | 统计量 | `--fst-global`：三个群体以上时算合并（r 群体 WC84）而非成对，默认成对；见 05 §2.1 |
| `--freq/--hardy/--het/--missing/--depth/--summary` | 报告 | — |
| `--window/--step`、`--levels`、`--skip-empty` | 窗口 | — |
| `--variant-only`、`--bypass-invariant-check` | 输入声明（pixy 语义） | 声明，非阈值 |
| `--contig-lengths FILE`、`--force-samples` | 输入补充 | — |

`--missing` 报告的 `F_MISS` 列是"缺失比例"（plink2 `.smiss/.vmiss` 的列定义），与阈值参数 `--geno/--mind` **同为缺失比例、同向**；帮助文本不必再解释两套方向。

---

## 6. gold 与脚本

- `gold/scripts/run_reed.sh:25` 恢复 `--geno "$SITE_MISS_MAX"`，与 `run_plink2.sh:29` 完全相同；**删掉**首版加的 `awk 'BEGIN{printf "%.10g", 1-m}'` 换算。
- `run_vcftools.sh:18` 保持 `--max-missing "$(awk … 1-m)"`（vcftools 侧仍需换算，`filters.conf` 是三家共同来源，不是 reed 的输入）。
- `SAMPLE_MISS_MAX` 在 gold 中用于生成 `samples.keep.txt`，reed 继续用 `-S`，不受影响。
- `gold/popgen` 不涉及这些参数。

---

## 7. 测试门【规】

1. **边界 fixture**（`reed/tests/e2e_missing_boundary.rs`，保留 W0 的用例表）：
   `N = 10` 与 `N = 7`（`X·2N` 分别为整数与非整数）、`X ∈ {0.0, 0.1, 0.2, 0.3, 0.5, 0.714, 0.72, 0.75, 0.8, 0.857, 1.0}`、含半调用与不含半调用；
   reed `--geno X` 的位点集合 == plink2 `--geno X`（含 ε 规则）== vcftools `--max-missing (1−X)`（换算由测试脚本用精确十进制完成）。
2. **样本级**：fixture 上 reed `--min-gq 20 --mind X` 的样本集合 == plink2 `--vcf-min-gq 20 --mind X`（保留集合 {s0,s2} 等 W0 断言，参数名换回）。
3. **回归**：subset 上 `--min-dp 5 --min-gq 20 --geno 0.1` → 263,420；`--geno 0.2` → 381,672（§2.3）。
4. **gold**：`run_reed.sh` 改回后，全量正文 md5 `7aedc59ac798729a2674e6e7ef19adfc` 不变，`check_parity.py` PASS；`##reed_filterCommand` 含 `--geno 0.1`。
5. **往返与守卫**：`to_args()` 渲染 `--geno/--mind`，`spec_from_flags(to_args(spec)) == spec`；三个撤回参数的报错文本各有单测（含 `--max-missing 0.9` → 提示 `--geno 0.1` 的精确十进制）；`ci/check-doc-drift.sh` 的守卫按 §5.2 翻转（`--geno/--mind` 存在、`--max-missing` 只出现在报错文本里）。
6. **精确判定**：单元测试构造 `X = 0.2, N = 10, missing = 2`（等号成立，保留）与 `missing = 3`（剔除），以及 `X = 0.3, N = 10, missing = 3`（`0.3` 无精确二进制表示，仍须保留）。

---

## 8. 与各轮文档的差异登记

| 文档条款 | 本文 |
|---|---|
| doc_v2 02 §3.1 `--mind`（plink2 丢弃方向） | **保留名字**，语义按 §3.2 精确定义（等位计、掩码后、头部存活记录） |
| doc_v2 02 §3.4 `--geno`（plink2） | **保留名字**，语义按 §3.1 精确定义 |
| doc_v2 §0 原则 2"名字相同则语义相同" | 保留 |
| doc_v2 05 §2.3"`--max-missing` 方向与 `--geno` 相反，不提供别名" | 保留其结论（不提供别名）；`--max-missing` 不再是 reed 的参数 |
| doc_v4 首版 02 §1.4"`--geno/--mind` 删除" | **作废**：见 §1.1 |
| doc_v4 首版 02 §5、§5.2、§5.3、§5.4、§6、§7 | **作废**，由本文同名小节取代 |

---

## 9. W0b 实施清单（供实现者）

| 文件 | 改动 |
|---|---|
| `crates/reed-cli/src/cli/filter_flags.rs` | `max_missing*` 三个参数改为隐藏的**报错**参数（方向与 W0 相反），报错函数命名为 `removed_max_missing`、`removed_max_missing_sample`、`removed_max_missing_count`（守卫脚本按这些名字 grep）；恢复可见的 `--geno X` / `--mind X`；`--max-missing-sample` 在 stats 中的拒绝文本改为提示 `reed filter --mind X` |
| `crates/reed-filter/src/spec/{mod,render}.rs` | 字段名回退（`min_call_rate`→`max_missing`、`min_sample_call_rate`→`mind`、删 `max_missing_alleles`）；`TailRule::MaxMissingCount` 删除；`to_args()` 渲染 `--geno/--mind`，不做算术 |
| `crates/reed-filter/src/eval/tail.rs` | 判定式改为 §4.2（缺失方向）；边界单测改名保留 |
| `crates/reed-filter/src/funnel.rs` | 漏斗 id 回退 `after_mask.missing` / `sample.mind`；`sample_rule.rule` 文本 `mind` |
| `crates/reed/src/filter/mod.rs`（pass0/sketch 样本通道） | 只改参数名与前向文本；判定式按 §4.3 |
| `crates/reed-cli/src/translate/filter.rs`、`plan/provenance.rs` | `--geno/--mind` 的往返与 `##reed_filterCommand` |
| `gold/scripts/run_reed.sh` | 删 awk，`--geno "$SITE_MISS_MAX"` |
| `ci/check-doc-drift.sh` | 守卫翻转（§7.5） |
| `docs/spec/architecture.md` §13.17、`docs/README.md` | 漏斗 id 与取代关系按本文回改 |
| `docs/rounds/v4/STATUS.md` | W0 表中涉及 CLI 拼写的行由 W0b 行**取代**（不删原行，注明 superseded） |
