# 02 CLI 规范（v14：当前有效的全部参数）

> **v14（2026-10-17）**：命令面重组——`reed stats` 拆为 QC 命令（本文件的
> §7，`--freq/--hardy/--het/--missing/--depth/--summary/--hist`）与窗口统计
> 命令 **`reed popg`**（别名 `popgen`）；新增 **`reed rename`**（§12）；
> `reed --help` 按组列出命令。窗口族 flag 从 `reed stats` 移走后在旧命令上
> 是**解析期报错 + 改写提示**（§12.3，`--max-missing` 家族同款）。0.3.x 的
> 输出逐字节不变（`gold/popgen/run_split_parity.sh` 22 case）。完整规范见
> `docs/rounds/v14/`。

> 本文把 `docs/rounds/v2/02-cli.md`（命令结构与原则）、`docs/rounds/v4/02-cli-semantics.md` 回退版（缺失阈值的精确语义与整数判定）和 2026-10-09 的三条裁决（`--geno/--mind` 保留原样；`--fst-global` 布尔开关；`--fst-pops` 作废）合成**一份当前有效的规范**，并逐条与 2026-10-09 实测的 `reed --help` 核对。
> 标 **【待实现】** 的条目是规范已定、代码未到；其余全部已实现且有测试。目录搬迁后本文即 `docs/spec/cli.md`（06-layout.md）。

---

## 0. 总原则【规】

1. **一条命令就是全部规范**。没有 `--spec`、`--conf`、TOML、JSON 输入。
2. **参数名 = 参考工具的参数名，语义也相同**。来源优先级：bcftools（I/O、区域、样本、表达式）→ plink2（群体 QC 阈值、基因型阈值、报告）→ vcftools / pixy（窗口统计）。名字相同而语义不同的参数不借用（vcftools `--minDP` 会掩掉 DP 缺失的基因型，不提供）。
3. **统计参数按科学方法命名**，帮助与日志给出发明者与文献，软件名只出现在测试里；缺失数据统一采用 Korunes & Samuk 2021（pixy）的逐位点计数法。
4. **短参数只用 bcftools 中含义相同的字母**。
5. **每个阈值参数都是一条"保留条件"**：要保证的量写下限 `≥`（`--maf/--mac/--hwe`），要抑制的量写上限 `≤`（`--max-maf/--max-mac/-Q/--geno/--mind`）。不存在"量 > X 则丢弃"的表述；CLI 路径上不做 `1−x` 换算（§4）。
6. **不为语义统一而改短名字**（2026-10-09 裁决）：`--geno/--mind` 是 plink2 的拼写，值 = 允许的缺失比例，保留；语义缺口用规范条文补（§4），不用改名补。
7. **参数在处理第一条记录之前全部校验**：标签是否声明、类型是否为数值、统计量所需群体是否存在；出错退出码 2。
8. **库 API 是类型化结构**（`FilterSpec`、`StatsRequest`），CLI 是它的文本表示，`to_args()` 可反向渲染，往返测试锁定。
9. **稳定性承诺只在 CLI 与输出格式**（doc_v5 08 §6）：参数语法与输出文件格式是契约；`reed` 与各 `reed-*` crate 都是内部结构，`publish = false`，Rust 公开项不承诺兼容。前端（`reed-cli`）只依赖门面 `reed`，后端不打 stderr、不认识 clap、不定退出码。

---

## 1. 命令

| 命令 | 对标 | 作用 |
|---|---|---|
| `reed filter` | bcftools `view` + `filter -S .` + `+fill-tags` + `annotate -x`；plink2 的 QC 阈值 | 过滤、基因型掩码、INFO 重算，输出保真的 VCF/BCF |
| `reed stats` | plink2 `--freq/--hardy/--het/--missing`、`bcftools stats`、vcftools 报告 | **QC 报告**：逐位点/逐样本行 + 汇总 + 直方图（§7） |
| `reed popg`（别名 `popgen`） | vcftools 窗口统计、pixy、easySFS、scikit-allel | **窗口群体统计**：π/dxy/Fst/θ/Tajima/PBS/SFS/f 统计/LD 衰减（原 `reed stats` 的统计族） |
| `reed rename`（别名 `reheader`） | bcftools `reheader`、`annotate --rename-chrs/--set-id`；plink2 `--set-*-var-ids/--update-name` | **标识符编辑**：样本名/染色体名/变异 ID/header 行（§12） |
| `reed index` | `bcftools index` | CSI（默认）/ TBI |
| `reed convert` | vcf2phylip 2.9（Ortiz 2019）、vcftools `--012`、TASSEL hapmap（双向）、vcf2msa、easySFS | 转换：系统发育矩阵 / 0-1-2 矩阵 / hapmap / 参考 MSA / 投影 SFS；hapmap → VCF 反向 |

**分组 help（v14 §1）**：`reed --help` 按 Filtering / Statistics / Transform /
File Operations / Miscellaneous 五组列出命令；命令路径本身不变（扁平）。
分组清单是 `cli/mod.rs` 的 `GROUPS` 常量，唯一事实源，守卫测试断言每个
`Command` 变体恰好入一组。

`reed qc` 是 doc_v3 C4 留下的**隐藏且已弃用**的子命令（`cli.rs:49`），`--help` 不列出；v5 路线 X1 删除（03 §4 C-2）。

---

## 2. 通用参数（实测 `--help` 核对）

| 参数 | 含义 | 命令 |
|---|---|---|
| `<INPUT>` / `-` | 输入 VCF / VCF.gz / BCF；`-` = stdin（`index` 只接受文件，索引需要可定位输入） | 三者 |
| `-o, --output FILE\|PREFIX` | filter：输出文件（缺省 stdout），格式来自 `-O` 或扩展名；stats：输出前缀；index：索引路径。别名 `--out` | 三者 |
| `-O, --output-type v\|z\|b\|u[0-9]` | 输出类型与压缩级别（`-Oz6`、`-Ob1`） | filter |
| `-W, --write-index[=csi\|tbi]` | 边写边建索引；**不给则不写**（bcftools 语义；W1 修正了"默认写"的缺陷） | filter |
| `--threads INT` | worker 线程数，默认全部核 | 三者 |
| `--memory MB` | 在途内存预算（MiB），默认物理内存一半 | filter / stats |
| `--tmp-dir DIR` | 临时目录（多遍运行的 stdin 溢写） | filter / stats |
| `--verbosity 0..3` | 0 静默；1 摘要（默认）；2 阶段计时；3 调试（`RUST_LOG`，默认 debug）；>3 报错 | 三者 |
| `--no-version` | 不写任何 `##reed_*` provenance 行 | filter |
| `--no-command` | 只省略 `##reed_command`（argv 含机器路径）；正文不受影响 | filter |
| `--report FILE` | JSON 报告（漏斗、样本规则、计数器） | filter |
| `--timing FILE` | 分阶段计时（open/pass0/mid/pass1、reader/reduce/chain_wait/worker_busy、在途峰值） | filter |
| `-h/--help`、`-V/--version` | — | — |

长度参数接受 `100000`、`100k`、`100kb`、`1M`、`1Mb`。

**与 `--help` 的差异登记**：stats 的 `--tmp-dir` 帮助仍写 "report runs that read twice"，W2 之后报告只读一遍，文本过期（01 F5-6）。

---

## 3. 过滤参数（`reed filter` 全部；`reed stats` 除 INFO 更新与 `--mind` 外全部，用于边过滤边统计）

### 3.1 样本级

| 参数 | 语义 | 来源 |
|---|---|---|
| `-s/--samples [^]LIST`、`-S/--samples-file FILE`、`--keep FILE`、`--remove FILE` | 静态子集；`^` 排除；`--force-samples` 时名单里不存在的样本只警告 | bcftools / plink2 |
| `--mind X`（仅 filter） | 保留样本缺失率 ≤ X 的样本；缺失率 = 缺失等位 / 2L，L = 通过掩码前位点规则的记录数，在掩码之后统计（§4.2） | plink2 |

静态子集与 `--mind` 的关系按 plink2：先静态子集，再动态规则。`reed stats` 不提供 `--mind`（需要多遍），给了就报错并提示先 `reed filter --mind X`。

### 3.2 位点级（掩码前）

| 参数 | 语义 | 来源 |
|---|---|---|
| `-v/--types LIST`、`-V/--exclude-types LIST` | `snps,mnps,indels,other,ref` | bcftools |
| `--snps-only[=just-acgt]` | 只保留 SNP；`just-acgt` 限制 ALT ∈ {A,C,G,T} | plink2 |
| `-m/--min-alleles INT`、`-M/--max-alleles INT` | REF+ALT 等位数上下限；`ALT='.'` 计 0 个 ALT（`-M 1` = 只留不变位点，与 bcftools 一致；doc_v12 N9） | bcftools |
| `-f/--apply-filters LIST` | FILTER 属于列表 | bcftools |
| `-i/--include EXPR`、`-e/--exclude EXPR` | bcftools 表达式子集（§5），可重复，按命令行顺序 | bcftools |
| `-q/--min-af FLOAT`、`-Q/--max-af FLOAT` | 非参考等位频率下限/上限；**INFO 有 AC 与 AN 且样本集未变时取 INFO 值**（bcftools 语义，取到时打印"陈旧 INFO"警告），一旦静态子集/掩码/`--mind`/`--set-gts .` 改变了 cohort，改按（子集/掩码后的）基因型重算——bcftools 同款门（`vcfview.c:341-350`、`filter.c:1471`；doc_v12 N5）；阈值按 bcftools 量化为 f32 | bcftools |

`-q/-Q` 的 `:nref/:alt1/:minor/:major/:nonmajor` 后缀 **【待实现】**（现只有 nref 一种，07 路线 X6）。

### 3.3 基因型掩码

| 参数 | 语义 | 来源 |
|---|---|---|
| `--min-dp A` / `--max-dp B` | `FMT/DP ∉ [A, B]` 的基因型置 `./.`；**DP 缺失不掩码**。别名 `--vcf-min-dp/--vcf-max-dp` | plink2 |
| `--min-gq G` | `FMT/GQ < G` 置缺失；GQ 缺失不掩码。别名 `--vcf-min-gq` | plink2 |
| `--set-gts .` + `-i/-e 'FMT/…'` | 与 `bcftools filter -S . -e EXPR` 相同；只接受 `.` | bcftools |

### 3.4 位点级（掩码后，按基因型重算）

| 参数 | 保留条件 | 分子 / 分母 |
|---|---|---|
| `--geno X` | 位点缺失率 ≤ X（§4.1 精确式） | 缺失等位 / 2N |
| `--maf X` / `--max-maf X` | `min(ac, an−ac)/an ≥ X` / `≤ X` | 次等位 / 已调用等位 |
| `--mac C` / `--max-mac C` | `min(ac, an−ac) ≥ C` / `≤ C` | 次等位计数 |
| `--hwe P [midp]` | 精确检验 `p ≥ P`；`midp` 用 mid-p | 双等位位点 |
| **`--indep-pairwise W[kb] [S] R2`**（裁决 2026-10-10） | LD 剪枝：按 plink2 的算法选出保留位点，**本命令直接输出剪枝后的 VCF/BCF**；位点清单只是 `--prune-list PREFIX` 的可选副产品 | 掩码后的基因型；r² 只取两位点都已调用的样本（05 §3.4） |

`-q/-Q` 与 `--maf` 的区别写在 `--help`：前者"非参考等位频率、优先用 INFO"，后者"次等位频率、总是按基因型重算"。

**LD 剪枝不照搬 plink 的两段式**（用户裁决 2026-10-10："不喜欢 prune.in / prune.out 两段式命令，但第一步选位点我需要"）：

- `reed filter in.vcf.gz --indep-pairwise 50 5 0.2 -Oz -o pruned.vcf.gz` 一条命令完成"选位点 + 提取"；plink2 要 `--indep-pairwise` 再 `--extract plink.prune.in` 两次运行。
- `--prune-list PREFIX`（可选）写 `PREFIX.prune.in` 与 `PREFIX.prune.out`：plink2 格式，变体 ID 一行一个，ID 为 `.` 时用 `chrom_pos`（与报告模式相同的规则），可直接喂给 plink2 / bcftools。不给就不写。
- `reed stats … --indep-pairwise W S R2` 同样可用：统计只在剪枝后的位点上算（边过滤边统计的一种）。
- 反向（拿别处产出的 ID 清单筛位点）是 `--extract FILE`（plink2 拼写），与 `--sites-file`（位置清单）一起归 X6。
- 它是**有顺序依赖的全局位点规则**，走两遍 / sketch 策略（pass0 按 contig 算保留位图，pass1 发射）；漏斗 id `after_mask.ld_prune`。算法、窗口语义、缺失处理、对拍门见 05 §3.4。

### 3.5 INFO 更新（仅 `reed filter`）

| 参数 | 语义 |
|---|---|
| 默认 | 样本被静态子集或 `--mind` 剔除时按 bcftools view 规则更新 `INFO/AC`、`INFO/AN`；只做掩码时不更新（与 `bcftools filter -S .` 一致） |
| `-I, --no-update` | 不更新 |
| `-x, --remove-info LIST` | 删除 INFO 标签；**默认清单**为 gold 管线的陈旧标签 `AC,AN,AF,MLEAC,MLEAF,ExcessHet,InbreedingCoeff,DP`（v2 02 §3.5），`-x ''` 保留全部 |
| `--fill-tags LIST` | 重算并写入 **AC、AN、NS、MAF、HWE**（大小写不敏感，定义同 `bcftools +fill-tags`）；未知标签报 usage 错误并列出支持集（bcftools 同款，`plugins/fill-tags.c:530-535`；doc_v12 N10）。AF/F_MISSING **【待实现】**。**多等位不做**（裁决 2026-10-10）：默认门 `-M 2` 已剔除；显式放开多等位又要求重算时报错，不静默写错（05 §4 存档） |

---

## 4. 缺失阈值的精确定义与整数判定【规】

记 `N` = 保留的样本数；`L` = 参与样本缺失计算的位点数（头部存活记录）；`an` = 位点上已调用等位数（半调用计 1，全调用计 2）；`called_s` = 样本 s 在 L 个位点上的已调用等位数。阈值 `X` 解析后量化为有理数 `q/10^9`（`reed_filter::spec::MissingRate`），比较全程 u128 整数，**没有浮点，没有 `1−x`**。

### 4.1 位点级 `--geno X`
保留 iff `(2N − an) · 10^9 ≤ q · 2N`。等价写法 `an ≥ (1−X)·2N`。**按等位计**（vcftools `--max-missing` 的定义），半调用贡献 1；plink2 按基因型计，两者只在有半调用时不同（本数据 0 个半调用）。阶段：掩码之后、`--maf/--mac/--hwe` 之前。

### 4.2 样本级 `--mind X`
保留 iff `(2L − called_s) · 10^9 ≤ q · 2L`。静态子集之后、掩码之后统计；`L` = 通过 `-v/-V/-m/-M/-f/--snps-only` 与不依赖基因型的 `-i/-e` 的记录数；在位点级阈值之前生效。需要一遍全局归约（pass0 / sketch）。

### 4.3 与参照物的等价性
vcftools 保留 iff `an/(2N) ≥ 1−X`（浮点）；plink2 丢弃 iff `missing > ⌊X·N·(1+ε)⌋`。对整数 `N` 与实数 `a` 有 `N − ⌈a⌉ = ⌊N − a⌋`，所以三者在 `X·2N` 恰为整数时同时取等号。已锁定的门：`e2e_missing_boundary`（N=10/7、含/不含半调用、重复小数边界 0.1428572 等）与 plink2 `--geno X`、vcftools `--max-missing (1−X)` 逐位点一致；子集 `--min-dp 5 --min-gq 20 --geno 0.1` → 263,420、`--geno 0.2` → 381,672（三家相同）。

### 4.4 撤回的拼写
`--max-missing X`、`--max-missing-sample X` 解析期报错并给**精确十进制**换算（`0.9` → `--geno 0.1`）；`--max-missing-count C` 报错并说明本族只有两个参数（计数上限写作 `C/2N`）。退出码 2，无静默别名。

### 4.5 漏斗 id 与报告
`after_mask.missing`（`--geno`）、`sample.mind`（`--mind`）。样本级规则不是 `funnel` 元素，而是报告里的 `sample_rule: {rule: "sample.mind", threshold: "X", n_dropped}`（已登记的形态）。

---

## 5. 表达式（`-i/-e`）

语法语义同 bcftools 手册 EXPRESSIONS，reed 实现其子集，超出子集报 "not supported yet"，不用另一种语义执行：

- 支持：`QUAL FILTER ID CHROM POS TYPE`、`INFO/TAG` 与裸名 `TAG`（裸名先查 INFO）、`FMT/TAG`、`N_ALT N_SAMPLES AC AN AF MAC MAF N_MISSING F_MISSING`；比较 `== = != < <= > >=`、正则 `~ !~`；逻辑 `&& || & |`（对 FORMAT 的含义照 bcftools 区分）、`!`、括号；缺失值 `"."`。
- 未支持：`FMT/` 下标 `[0]/[*]`、通用正则、函数（`MAX/MIN/AVG/…`）。
- `AC/AN/AF/MAC/MAF`：INFO 中存在时取 INFO 值，否则按基因型。检测到表达式使用了 INFO 中已存在的 AC/AN 时打印"陈旧 INFO"警告，提示改用 `--maf/--mac` 或先 `-x INFO/AC,INFO/AN`。
- 与缺失值比较为假（bcftools 语义）。
- 引用基因型量的表达式（`N_MISSING`、`F_MISSING`，以及 INFO 不存在时的 `AC/AN/…`）在掩码之后求值。
- **Flag 按存在性读写**：`DB=1`/`DB=0` 测 flag 的存在/缺失（`doc/bcftools.txt:4067`），`DB!=0` 等价 `DB=1`；`DB` 单独出现不是合法表达式（bcftools 对单个裸 token 恒假，reed 在编译期报错而不是静默丢弃所有记录；doc_v12 N6）。
- **多值 INFO 按任一元素比较**：`Number=A/R/G/.` 或 `Count(n>1)` 的标签逐对比较，任一对满足即真（bcftools `CMP_VECTORS`）；`Number=1` 却给了多个值时按 missing 处理（规范外输入，01-review §2.5-B5）。
- **数值解析**：按 header 声明类型选择文法（Integer 严格十进制含可选 `+`/`-`，Float 为规范正则文法）；解析失败即 missing，不做字节序比较（doc_v12 N13）。

**与 htslib 的有意分歧**（doc_v12 §2.5-B，永久保留，只登记不改）：

| 形态 | htslib | reed | 理由 |
|---|---|---|---|
| `-Ov`/`-Oz` 的 INFO 文本 | 归一化（`+5`→`5`） | 逐字节保真 | §6.7 的编辑原语设计；保真严格强于归一化 |
| `Σac > an`（AC/AN 自相矛盾） | `exit(1)` 退出 | 丢该位点继续 | 一条畸形记录不应毁掉整轮批处理 |
| `5.0`/`1e3`/`0x10`/`5abc`（Integer 标签） | 前缀解析（5 / 1 / 0 / 5） | missing | htslib 的前缀解析自身不自洽（`0x10` 在 Float 得 16）；严格拒绝更可预测 |
| ` 5`/`5 `（数值） | Integer 不容忍、Float 容忍 | 不 trim，missing | 数值文法（§1.3）无空白；htslib 两类型不一致 |
| `Number=1` 却给多值 | 标量读作 0（`vcf.c:4286` 的实现细节） | missing | 不把 htslib 的偶然行为固化成契约（02 §6.2b） |
| `-q/-Q` 阈值的 f32 量化 | C 的 `float` 存储（无意） | 刻意复现 `sscanf("%f")` | 已对齐，防止未来改成 f64 破坏边界案例 |

---

## 6. `reed filter` 的求值顺序与 header

| 步 | 内容 |
|---|---|
| 1 | `-s/-S/--keep/--remove`、`-r/-R/-t/-T` |
| 2 | 载入期位点过滤：`-v/-V/--snps-only/-m/-M/-f` |
| 3 | 不依赖基因型的 `-i/-e` |
| 4 | 掩码：`--min-dp/--max-dp/--min-gq`、`--set-gts .` + `FMT/` 表达式 |
| 5 | `--mind` |
| 6 | `--geno`、`--hwe`、`--maf/--max-maf/--mac/--max-mac`、依赖基因型的 `-i/-e` |
| 6b | `--indep-pairwise`：有顺序依赖的全局位点规则，只看通过第 1–6 步的位点与掩码后的基因型（plink2 的顺序：`--mind/--geno/--maf` 之后）；pass0 算保留位图，pass1 发射 |
| 7 | INFO 更新：`-x` → AC/AN 自动更新 → `--fill-tags` |
| 8 | 输出 |

顺序只影响漏斗计数，不影响结果；唯一影响结果的是"`--mind` 在掩码之后"。

header：`##reed_version`、`##reed_filterCommand=<to_args() 规范化命令>; Date=…`、`##reed_command=<原样 argv>`（`--no-command` 省略；`--no-version` 省略全部）。正文确定，header 含环境信息。

三种策略（single_pass / two_pass_reparse / sketch_then_emit）输出逐字节相同；VCF 与 BCF 输入、`-Oz` 与 `-Ob` 输出正文相同（测试锁定）。

---

## 7. `reed stats`（QC 报告）

> v14：本节的窗口参数、群体族与统计族已搬到 `reed popg`（下一节）；`stats`
> 只保留 QC 报告族、共享选择参数与 `--hist`。

```
reed stats [通用] [选择参数] --freq/--hardy/--het/--missing/--depth/--summary/--hist -o PREFIX <INPUT|->
```

### 7.1 原则
统计参数只写方法；缺失数据按 pixy 逐位点计数；不变位点默认要求在文件里（pixy 的检查），只有变异位点时用 `--variant-only` 声明；窗口从 1 开始。无缺失数据时与 vcftools / scikit-allel / pixy 逐位一致。

### 7.2 统计参数

| 参数 | 统计量 | 方法 | 状态 |
|---|---|---|---|
| `--pi` | π | Nei & Li 1979；缺失：Korunes & Samuk 2021 | 已对拍 pixy（4.8e-14） |
| `--dxy` | dxy | Nei 1987；同上 | 已对拍 |
| `--theta` | Watterson θ | Watterson 1975；逐位点 n | 已对拍 |
| `--tajima-d` | Tajima's D | Tajima 1989；方差按窗口平均调用数算一次（pixy） | 已对拍 |
| `--fst [wc84\|hudson]` | 成对 Fst，所有群体对；默认 `wc84` | Weir & Cockerham 1984；Hudson 1992 / Bhatia 2013 | 已对拍 |
| **`--fst-global`** | 布尔开关：三个群体以上时算**合并**的 r 群体 WC84（一个值），而非成对；默认成对 | Weir & Cockerham 1984（r 群体） | **【待实现】**（05 §3.2；裁决 2026-10-09） |
| `--pbs A,B,C [--pbs-normed]` | PBS，A 为焦点群体；可重复 | Yi et al. 2010；PBSn1：Malaspinas 2016 | 已对拍 allel |
| `--sfs [folded] --proj N[,N…]` | 投影 SFS：每群体一维；两群体起另写每对联合（二维）谱 | Marth 2004；Gutenkunst 2009 | 已对拍 easySFS（`run_easysfs.sh`，含二维联合谱） |
| `--freq` / `--hardy [midp]` / `--het` / `--missing` | 逐位点频率 / HWE 与杂合度 / 个体 F / 缺失率 | Wigginton 2005；Graffelman 2013 | 已对拍 plink2 |
| `--depth` | 逐样本、逐位点平均深度 | — | 已对拍 vcftools |
| `--summary` | bcftools stats 格式 QC 汇总 | — | 已对拍 |
| **`--f2/--f3/--f4/--dstat/--f4-ratio`**、**`--ld-decay [MAX_BP]` / `--ld-pairs`** | f 统计（block jackknife）、LD 衰减 | Patterson 2012；Purcell 2007 | f 统计 X4.3 已落地；LD 衰减见 05 §3.4.2。`--indep-pairwise` **不是** stats 独有的统计量，它是 §3.4 的过滤规则，stats 里的含义是"只在剪枝后位点上统计" |

`--fst-global` 的规则：没有 `--fst` 时报错；与 `hudson` 组合报错（Hudson 没有 r 群体形式）；两个群体时值等于成对 wc84，只提示不报错；输出单独的 `PREFIX_fst_global.txt`。`--fst-pops` **作废**（doc_v2 02 §5.2 那一行由本表取代）。

### 7.3 群体、窗口与声明

| 参数 | 语义 | 状态 |
|---|---|---|
| `--pop FILE`（别名 `--populations`） | `sample<TAB>population`，无表头；缺失样本默认报错（pixy），`--force-samples` 改为跳过 | 已实现 |
| `--window SIZE`、`--step SIZE` | 窗口与步长（默认等于窗口）；网格 `[1,w]`、`[1+s,w+s]`… | 已实现 |
| `--levels window,contig,genome` | 输出层级，默认 window | 已实现 |
| `--skip-empty` | 不输出没有位点落入的窗口（默认输出并写 `NA`） | 已实现 |
| `--contig-lengths FILE` | header 缺 `##contig` 长度时补充（`.fai` 两列）；缺失即报错，不编造 | 已实现 |
| `--variant-only` | 声明只有变异位点：窗口内未出现的位置视为单态且全部调用，计入 π/dxy/θ 分母 | 已实现（== vcftools `--window-pi`） |
| `--bypass-invariant-check` | 不含不变位点也不声明时强行只统计出现的位点（高估 π/dxy，警告） | 已实现 |
| `--bed-file BED`、`--sites-file FILE` | BED3（0-based 半开）给定窗口，区间可重叠、各自成窗；`CHROM<TAB>POS`（1-based）只统计列出的位点；`--bed-file` 与 `--window/--step` 互斥 | 已实现（BED 与 pixy `--bed_file` 逐列对拍，`gold/popgen/run_pixy_bed.sh`） |

### 7.4 输出（`-o PREFIX`；坐标 1-based 闭区间；缺失写 `NA`）

| 统计 | 文件 | 列（实测 writer） |
|---|---|---|
| π | `PREFIX_pi.txt` | `pop chromosome window_pos_1 window_pos_2 avg_pi no_sites count_diffs count_comparisons count_missing` |
| dxy | `PREFIX_dxy.txt` | `pop1 pop2 chromosome window_pos_1 window_pos_2 avg_dxy no_sites count_diffs count_comparisons count_missing` |
| Fst | `PREFIX_fst.txt` | `pop1 pop2 chromosome window_pos_1 window_pos_2 avg_wc_fst\|avg_hudson_fst no_snps`（列名取 pixy 的拼写 `avg_wc_fst`，使下游 pixy 脚本可直接用） |
| Fst global | `PREFIX_fst_global.txt` | `pops chromosome window_pos_1 window_pos_2 avg_wc_fst no_snps`【待实现】 |
| θ | `PREFIX_watterson_theta.txt` | `pop chromosome window_pos_1 window_pos_2 avg_watterson_theta no_sites raw_watterson_theta no_var_sites weighted_no_sites` |
| Tajima's D | `PREFIX_tajima_d.txt` | `pop chromosome window_pos_1 window_pos_2 tajima_d no_sites raw_pi raw_watterson_theta tajima_d_stdev` |
| PBS | `PREFIX_pbs.txt` | `pop1 pop2 pop3 chromosome window_pos_1 window_pos_2 pbs no_snps fst_12 fst_13 fst_23` |
| SFS | `PREFIX_sfs/dadi/<pop>-<n>.sfs`、`PREFIX_sfs/fastsimcoal2/<pop>_MAFpop0.obs`；联合（二维，每对一份）：`PREFIX_sfs/dadi/<A>-<B>.sfs`、`PREFIX_sfs/fastsimcoal2/<prefix>_jointMAFpop<idxB>_<idxA>.obs` | dadi / fastsimcoal2 标准格式，含 mask 行（easySFS 实测布局） |
| freq / hardy / het / missing | `PREFIX.afreq .hardy .het .smiss .vmiss`；`--pop` 时另写 `.afreq.pop`/`.vmiss.pop`（每群体一行，§12.2） | plink2 列布局 |
| depth | `PREFIX.idepth`、`PREFIX.ldepth.mean` | vcftools 列布局 |
| summary | `PREFIX.stats` | bcftools stats 文本 |
| f 统计 / LD | `PREFIX.fstats`、`PREFIX.ld.decay`、`PREFIX.ld`（`--ld-pairs`） | 05 §3.3–3.4。剪枝清单 `.prune.in/.prune.out` 不是 stats 的输出，只由 `--prune-list PREFIX` 产生（filter 与 stats 都可给） |
| 日志 | `PREFIX.log` | 命令、版本、每个统计量的方法与文献、过滤漏斗、计数器 |

### 7.5 例子

```bash
reed stats allsites.vcf.gz --pop pops.txt --pi --dxy --fst --tajima-d --window 10k -o out --threads 56
reed stats clean.vcf.gz --pop pops.txt --variant-only --pi --fst --fst-global --window 100k --step 50k -o out
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc
reed stats clean.vcf.gz --pop pops.txt --freq --missing --singletons -o qc   # 另写 .afreq.pop/.vmiss.pop，singletons 加 POP 列
reed stats clean.vcf.gz --pop pops.txt --variant-only --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 --pi --window 100k -o out
reed stats clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out
```

---

## 8. `reed index`

`reed index [-c|-t] [-o FILE] [-f] [--threads N] <FILE>`；默认 CSI，`-t` TBI；stdin 明确拒绝（"needs a seekable file"）。输出与 bcftools 的索引在 `bcftools view -r` 下逐字节等价（测试锁定）。

---

## 8b. `reed convert`（doc_v13）

转换命令：一次读入，多目标输出（同一 worker 池，输出与线程数无关）。至少给一个转换目标，否则 usage 错。

| 参数 | 语义 | 来源 |
|---|---|---|
| `--sfs [folded] --proj N[,N…] --sfs-int --pop FILE` | 投影 SFS：与 `reed stats --sfs` **逐字节等价**（同前缀同输入下 `diff -r` 为空）；`--pop` 缺省全体一组 | easySFS |
| `--phylip` / `--fasta` / `--nexus` | 核苷酸矩阵（relaxed PHYLIP 布局 / FASTA / NEXUS）；三格式共享一份列模型，**`--phylip` 是 opt-in**（与 vcf2phylip 默认写相反，B13-1） | vcf2phylip |
| `--nexus-binary` | 二元 NEXUS（SNAPP；仅双等位、二倍体 GT 参与，其余 `?`） | vcf2phylip |
| `--outgroup NAME` | 该样本写在每个矩阵第一行；不在（子集后的）样本集 → usage 错（B13-4） | vcf2phylip |
| `--min-samples-locus N` | 位点门：有 GT 的**样本**数下限（phylo 默认 4，超样本数自动截断；`--012`/`--hapmap` 只在显式给出时应用） | vcf2phylip |
| `--resolve-iupac` | 杂合随机拆解为单字母；**确定性**（seed 0，threads {1,4} 输出相同），与 vcf2phylip 的每次不同有意分歧（B13-3） | vcf2phylip |
| `--write-used-sites` | 写 `PREFIX.used_sites.tsv`（仅核苷酸流接受的位点） | vcf2phylip |
| `--012` | `PREFIX.012`（每样本一行 `<原header下标>\t0/1/2/-1`）、`PREFIX.012.pos`、`PREFIX.012.indv`；双等位门（`>2` 等位跳过），**半调用与 haploid 写 `-1`**（vcftools 写 1，B13-8） | vcftools |
| `--hapmap` | `PREFIX.hmp.txt`：11 固定列 + 样本列；`alleles` 是 REF/ALT **声明集**（TASSEL 写观察集，B13-5）；格取 TASSEL 读侧字母表（A/C/G/T/R/Y/S/W/K/M/N/-，其余 → `N`） | TASSEL |
| `--msa --reference FASTA [--region CHR:B-E] [--flank N] [--collapse-suffix STR]` | 参考驱动的全位点 MSA：每 contig/区域一个 `PREFIX_<contig>[_<s>-<e>].fasta`，无记录处写 REF；indel 位点确定性填充（不调外部对齐器，B13-6）；样本名默认不折叠（B13-7） | vcf2msa |
| `--vcf [-O v\|z[level]] [--no-version] [--no-command] [-o FILE]` | **反向转换**：输入是 hapmap（`.hmp.txt`），输出 VCF。与每个 VCF 输入目标互斥（门规则 6），不接受 `-s/-S/--keep/--remove`。`chrom`/`pos`/`rs#`/`alleles` 直抄（`ID` = `rs#`，REF/ALT 按声明顺序，B13-11）；格按 IUPAC 反查展开为 GT（`R`→`0/1` 等），仅 `N`/`?` 是缺失（`./.`）；`-`/`+`/`0` 按 `alleles` 列查下标；`alleles=NA` → 全 `./.` 行。`-O` 缺省 `v`（纯文本），`-o` 以 `.gz`/`.bgz` 结尾时缺省 `z6`；`b`/`u`（BCF）报 usage 错。头无 `##contig`（TASSEL 同款），`##reed_version`/`##reed_command` 由 `--no-version`/`--no-command` 分别关闭（B13-12） | TASSEL |
| `-o PREFIX` | 输出前缀（缺省：输入文件名去扩展）；`PREFIX.log` 写五计数与参数。`--vcf` 下 `-o` 是**输出文件路径**（缺省 stdout） | reed |

**与参照工具的有意分歧（B13 登记，不写对拍断言）**：B13-1 `--phylip` opt-in；B13-2 畸形行跳过 + 计数 + 一次性警告（vcf2phylip 打印整行）；B13-3 resolve 固定 seed；B13-4 outgroup 报错；B13-5 hapmap 声明集；B13-6 msa 不调外部对齐器；B13-7 msa 样本名不折叠；B13-8 `--012` 半调用/haploid = `-1`；**B13-9** AMBIG 表未命中（REF/ALT 含表外字母）vcf2phylip **直接崩溃**，reed 跳过 + 计数；**B13-11** hmp→VCF 的 `ID`/REF/ALT 直抄声明值（TASSEL 的 `ExportUtils` 会重排/重命名）；**B13-12** hmp→VCF 头不写 `##contig`（TASSEL 也不写）且 `##reed_*` 可关。

**对拍门**：`gold/popgen/run_convert_parity.sh`（vcf2phylip 五文件逐字节 / vcftools `--012` 干净夹具逐字节 + B13-8 断言 / TASSEL 读回 / vcf2msa 序列相等 / convert ≡ stats 的 SFS 树 / hapmap↔VCF 往返 + TASSEL 读 reed 的 VCF）；可选 `DATA.vcf.gz` 参数追加全量档（同一负载跑参照与 reed 各一遍并 diff）。**实测**（clean.vcf.gz，8 GB / 348 样本 / 331 万记录，28 核）：五格式一次跑 vcf2phylip 2.9 = 27 min 27 s，`reed convert` = 10.4 s（≈158×，五文件逐字节相同；v13 STATUS §基线）。

---

## 9. 退出码与消息

0 成功；1 运行期失败（带 `contig:pos`）；2 用法错误（参数、表达式、标签、群体缺失、不变位点检查未过、撤回的参数）。错误 `reed: …`（Display，不带引号）；警告 `reed: warning: …`。

---

## 10. 库 API 与往返

`FilterSpec`（`spec/{cohort,mask,missing_rate,predicate,recompute,region,render}.rs`）与 `StatsRequest` 是类型化结构；`to_args()` 渲染为 §3/§7 的参数且 `parse(render(x)) == x`（`translate::filter::tests::spec_to_args_round_trips_through_the_flags` 等）。clap 结构体、单位解析、错误渲染只在 `reed-cli`。

---

## 12. v14 新增面（doc_v14）

### 12.1 `reed popg`（窗口群体统计）

原 `reed stats` 的窗口统计族，flag 拼写、默认值、输出文件与列布局全部不变
（§7.2–§7.4 的统计族与 §7.3 的窗口/清单参数即为 `popg` 的参数面）。
`--pop FILE` 在窗口模式必需（现状不变）。`--sfs` 的等价契约更新为
`reed popg --sfs` ≡ `reed convert --sfs` 逐字节（原 `stats --sfs`）。

### 12.2 `reed stats` 的 QC 面

`--freq/--hardy/--het/--missing/--depth/--summary` 的拼写、输出文件与列布局
零变化。新增 `--hist LIST`：逗号分隔的内建指标（`qual`、`dp`、`dp-site`、
`gq`、`idl`、`af`）与通用 `INFO/<TAG>`/`FMT/<TAG>`；每个 metric 写
`PREFIX.<m>.hist` TSV（列 `bin_lower bin_upper count percent`；整数指标附
精确 `#mean`/`#median`/`#p95` 行，缺失观测附 `#missing` 行）。`idl` 用有
符号长度（缺失为负、插入为正，bcftools IDD 约定）。标签必须在 header 声明
且为 Integer/Float，否则用法错（§0.7）。`--hist` 可单独使用（裸
`reed stats` 仍是用法错，但给了 `--hist` 即为合法请求）。

另有 `--singletons`：逐等位（含 REF）输出 singleton（计数 1）与 private
doubleton（计数 2 且有纯合携带者），列 `CHROM POS SINGLETON/DOUBLETON
ALLELE INDV`（vcftools 布局，`PREFIX.singletons`）。

`--pop FILE`（`--force-samples` 同义沿用）分层 QC：`--freq`/`--missing` 各
另写一份按群体的行——`PREFIX.afreq.pop`（列
`#CHROM ID POP REF ALT ALT_FREQS OBS_CT`）与 `PREFIX.vmiss.pop`（列
`#CHROM ID POP MISSING_CT OBS_CT F_MISS`），每群体一行、按文件首现序；
`--singletons` 加第六列 `POP`（携带者的群体，不在群体文件中的样本写
`.`，计数仍按全体样本）。分层值 == 把该群体样本单独取出后的 plink2
`--freq/--missing` 值（`gold/popgen/run_popqc.sh`，逐群体对拍；多等位
ALT 频率为折叠单值，与 `.afreq` 的 `CollapseAlt` 约定一致）。不给
`--pop` 时全部输出与 0.3.x 逐字节相同；给了 `--pop` 但没有
`--freq/--missing/--singletons` 之一时为用法错（不静默忽略）。多等位位点
不拆分（`CollapseAlt` 语义：ALT_FREQS 是全部 ALT 等位的合并频率）——
`stats` 全家族一致，非分层独有。

### 12.3 搬迁 flag 的重定向

§7 窗口族与统计族 flag 出现在 `reed stats` 上时，是**解析期用法错（exit 2）
+ 完整改写提示**，例如：

```
reed: --pi moved to `reed popg`: the windowed population-genetics statistics
      live there now; run `reed popg <input> --pop FILE --pi ...`
```

无静默别名、无转发执行（同 `--max-missing` 家族，§4.4）。

### 12.4 `reed rename`

四个平面，求值序 header 替换 → 染色体 → ID → 样本 → header meta 编辑；
只编辑 header 行、CHROM 列、ID 列与样本名，**从不解码基因型**（结构性守卫：
`crates/reed/src/rename/` 不引用 GenotypeScanner/GtCode），记录其余字节直通。

| 平面 | flag |
|---|---|
| header | `-h/--header FILE`、`-f/--fai FILE`、`--add-meta LINE`、`--remove-meta KEY` |
| 染色体 | `--rename-chrs FILE`（`##contig` 行与 CHROM 列同改）、`--chr-style plain\|ucsc\|mt-is-m` |
| ID | `--set-id PATTERN`（`%CHROM %POS %REF %ALT`；ID 用改名后的 CHROM）、`--set-missing-var-ids PATTERN`、`--update-name FILE`、`--idspace-to SEP`、`--rm-dup-ids`、`--clear-ids` |
| 样本 | `-s "old new,…"`、`-N/--samples-file FILE`、`--replace-sample-name PATTERN` |

输出 `-o FILE`（默认 stdout）、`-O v\|z[level]`（BCF 输入/输出暂不支持）、
`-W[=csi\|tbi]`（filter 语义：默认不写索引；索引按改名后的 contig 名建）。
输出旁若存在陈旧索引会打印 warning 提示 `reed index -f`。
校验（§4.3）：改后样本名重复、名字不存在、`--set-id` 与
`--set-missing-var-ids` 同给、`--clear-ids` 与其他 ID flag 同给、
`--idspace-to` 与 `--set-id` 同给、无任何编辑 flag —— 全部用法错。

## 11. 本轮核对出的差异（登记到 01-review.md）

| # | 差异 | 处理 |
|---|---|---|
| ~~D-1~~ | `--bed-file`、`--sites-file`（doc_v2 §5.3）已实现（Y6） | 完成 |
| D-2 | `-q/-Q` 的 `:type` 后缀未实现 | 路线 X6 |
| D-3 | stats `--tmp-dir` 帮助文本过期（"read twice"） | 路线 X1 |
| D-4 | `reed qc` 隐藏子命令仍在 | 路线 X1（03 C-2） |
| D-5 | `--fst-pops` 仍出现在 `reed-stats/src/lib.rs:47` 注释与 doc_v2 02 §5.2 | 路线 X4 随 `--fst-global` 清理 |
| D-6 | Fst 列名文档写 `avg_wc84_fst`、代码写 `avg_wc_fst` | 本文以代码（pixy 拼写）为准 |
| D-7 | 在途 X4.4 把 `--indep-pairwise` 放在 `reed stats` 并以 `.prune.in/.out` 为产物（照 05 §3.4 旧文） | 2026-10-10 裁决改为 §3.4 的过滤规则 + 一条命令出剪枝后 VCF；`--ld-decay/--ld-pairs` 留在 stats 不变（05 §3.4.3） |
