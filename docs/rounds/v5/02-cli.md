# 02 CLI 规范（v5 合并版：当前有效的全部参数）

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
| `reed stats` | vcftools 统计、pixy、plink2 `--freq/--hardy/--het/--missing`、easySFS、`bcftools stats` | 所有统计一个命令，一个统计量一个参数 |
| `reed index` | `bcftools index` | CSI（默认）/ TBI |

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
| `-m/--min-alleles INT`、`-M/--max-alleles INT` | REF+ALT 等位数上下限 | bcftools |
| `-f/--apply-filters LIST` | FILTER 属于列表 | bcftools |
| `-i/--include EXPR`、`-e/--exclude EXPR` | bcftools 表达式子集（§5），可重复，按命令行顺序 | bcftools |
| `-q/--min-af FLOAT`、`-Q/--max-af FLOAT` | 非参考等位频率下限/上限；**INFO 有 AC 与 AN 时取 INFO 值**（bcftools 语义，本数据的 INFO/AC 是陈旧的，取到时打印"陈旧 INFO"警告），否则按基因型；阈值按 bcftools 量化为 f32 | bcftools |

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
| `-x, --remove-info LIST` | 删除 INFO 标签 |
| `--fill-tags LIST` | 重算并写入 AC、AN、AF、MAF、NS、HWE、F_MISSING…，定义同 `bcftools +fill-tags`。**多等位不做**（裁决 2026-10-10）：默认门 `-M 2` 已剔除；显式放开多等位又要求重算时报错，不静默写错（05 §4 存档） |

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

## 7. `reed stats`

```
reed stats [通用] [过滤参数] [--pop FILE] [窗口参数] <统计参数…> -o PREFIX <INPUT|->
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
| **`--sfs [folded] --proj N[,N…]`** | 投影 SFS（1 群体一维，2 群体二维） | Marth 2004；Gutenkunst 2009 | **【待实现】**（05 §3.1，W5 进行中） |
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
| `--bed-file BED`、`--sites-file FILE` | BED 给定窗口；只统计列出的位点 | **【待实现】**（doc_v2 承诺，X6） |

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
| SFS | `PREFIX_sfs/dadi/<pop>-<n>.sfs`、`PREFIX_sfs/fastsimcoal2/<pop>_MAFpop0.obs`（二维：`<A>-<B>-<nA>-<nB>.sfs`、`jointMAFpop1_0.obs`） | dadi / fastsimcoal2 标准格式，含 mask 行【待实现】 |
| freq / hardy / het / missing | `PREFIX.afreq .hardy .het .smiss .vmiss` | plink2 列布局 |
| depth | `PREFIX.idepth`、`PREFIX.ldepth.mean` | vcftools 列布局 |
| summary | `PREFIX.stats` | bcftools stats 文本 |
| f 统计 / LD | `PREFIX.fstats`、`PREFIX.ld.decay`、`PREFIX.ld`（`--ld-pairs`） | 05 §3.3–3.4。剪枝清单 `.prune.in/.prune.out` 不是 stats 的输出，只由 `--prune-list PREFIX` 产生（filter 与 stats 都可给） |
| 日志 | `PREFIX.log` | 命令、版本、每个统计量的方法与文献、过滤漏斗、计数器 |

### 7.5 例子

```bash
reed stats allsites.vcf.gz --pop pops.txt --pi --dxy --fst --tajima-d --window 10k -o out --threads 56
reed stats clean.vcf.gz --pop pops.txt --variant-only --pi --fst --fst-global --window 100k --step 50k -o out
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc
reed stats clean.vcf.gz --pop pops.txt --variant-only --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 --pi --window 100k -o out
reed stats clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out        # 【待实现】
```

---

## 8. `reed index`

`reed index [-c|-t] [-o FILE] [-f] [--threads N] <FILE>`；默认 CSI，`-t` TBI；stdin 明确拒绝（"needs a seekable file"）。输出与 bcftools 的索引在 `bcftools view -r` 下逐字节等价（测试锁定）。

---

## 9. 退出码与消息

0 成功；1 运行期失败（带 `contig:pos`）；2 用法错误（参数、表达式、标签、群体缺失、不变位点检查未过、撤回的参数）。错误 `reed: …`（Display，不带引号）；警告 `reed: warning: …`。

---

## 10. 库 API 与往返

`FilterSpec`（`spec/{cohort,mask,missing_rate,predicate,recompute,region,render}.rs`）与 `StatsRequest` 是类型化结构；`to_args()` 渲染为 §3/§7 的参数且 `parse(render(x)) == x`（`translate::filter::tests::spec_to_args_round_trips_through_the_flags` 等）。clap 结构体、单位解析、错误渲染只在 `reed-cli`。

---

## 11. 本轮核对出的差异（登记到 01-review.md）

| # | 差异 | 处理 |
|---|---|---|
| D-1 | `--bed-file`、`--sites-file`（doc_v2 §5.3）未实现 | 路线 X6 |
| D-2 | `-q/-Q` 的 `:type` 后缀未实现 | 路线 X6 |
| D-3 | stats `--tmp-dir` 帮助文本过期（"read twice"） | 路线 X1 |
| D-4 | `reed qc` 隐藏子命令仍在 | 路线 X1（03 C-2） |
| D-5 | `--fst-pops` 仍出现在 `reed-stats/src/lib.rs:47` 注释与 doc_v2 02 §5.2 | 路线 X4 随 `--fst-global` 清理 |
| D-6 | Fst 列名文档写 `avg_wc84_fst`、代码写 `avg_wc_fst` | 本文以代码（pixy 拼写）为准 |
| D-7 | 在途 X4.4 把 `--indep-pairwise` 放在 `reed stats` 并以 `.prune.in/.out` 为产物（照 05 §3.4 旧文） | 2026-10-10 裁决改为 §3.4 的过滤规则 + 一条命令出剪枝后 VCF；`--ld-decay/--ld-pairs` 留在 stats 不变（05 §3.4.3） |
