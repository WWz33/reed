# 02 CLI 规范（取代 v1 §11.1 的 TOML Spec 与 §18 CLI）

> 用户裁决：
> 1. 2026-10-06 第一次：**不要 conf、不要 TOML**，采用 bcftools / vcftools / plink2 风格的命令行参数；
> 2. 2026-10-06 第二次：**统计参数要最直白最简单；其他所有行为与主流程序、主流算法对齐。reed 不做算法创新，只做架构创新。**
>
> 3. 2026-10-06 第三次：**统计参数只报告科学方法（谁发明的方法），不报告软件名**；缺失数据的处理以 pixy 为默认。
>
> 由此得到本规范的总原则：**reed 不发明统计方法**。每个统计量都是有文献的经典方法，缺失数据统一用 pixy 的方法；
> reed 提供的是速度（多核、单遍、保真），不是新的统计定义。

---

## 0. 总原则【规】

1. **一条命令就是全部规范**。不存在 `--spec`、`--conf`、TOML、JSON 输入。
2. **参数名 = 参考工具的参数名，语义也必须相同**。来源优先级：
   bcftools（通用 I/O、区域、样本、表达式）→ plink2（群体 QC 阈值、基因型阈值、`--freq/--hardy/--het/--missing`）
   → vcftools / pixy（窗口统计）。**名字相同而语义不同的参数一律不借用**（例如 vcftools 的 `--minDP` 会把 DP 缺失的基因型
   也掩掉，与 plink2/bcftools 相反，因此不提供这个别名，§3.3）。
3. **统计参数以科学方法命名**（§5.1）：每个统计量对应一个有文献的方法；缺失数据统一采用 Korunes & Samuk 2021（pixy）的方法。
   软件名不出现在参数中，只出现在测试的对拍里（无缺失数据时与 vcftools/scikit-allel/pixy 逐位一致）。
4. **短参数只用于 bcftools 中含义相同的字母**。
5. **参数在处理第一条记录之前全部校验**：标签是否声明、类型是否为数值、统计量所需的群体是否存在。出错退出码 2。
6. **库 API 仍是类型化结构**（`FilterSpec`、`StatsRequest`，builder 构造）；CLI 是它的文本表示，可反向渲染（§7）。

---

## 1. 命令集合（解决【决 CLI-1】）

| 命令 | 对标 | 作用 |
|---|---|---|
| `reed filter` | bcftools `view` + `filter -S .` + `+fill-tags` + `annotate -x`；plink2 的 QC 阈值 | 过滤、基因型掩码、INFO 重算，输出保真的 VCF/BCF |
| `reed stats` | vcftools 的统计输出、pixy、plink2 `--freq/--hardy/--het/--missing`、easySFS、`bcftools stats` | **所有统计合并为一个命令**：一个统计量一个参数，参数按科学方法命名 |
| `reed index` | `bcftools index` | 建 CSI（默认）/ TBI 索引 |

裁决理由：vcftools 和 plink2 都是"一个程序、每个统计量一个参数"，用户最熟悉，也最直白。
`reed stats --summary` 产出 `bcftools stats` 格式的 QC 报告，所以不需要单独的 `qc` 命令，也就不存在与 bcftools 的
`stats` 撞名的问题。上一版 doc_v2 中的 `reed popgen` **取消**。legacy 二进制 `reed-filter` 删除。

---

## 2. 通用参数（所有命令）

| 参数 | 含义 | 来源 |
|---|---|---|
| `<FILE>` / `-` | 输入；`-` = stdin；格式按魔数探测（BGZF VCF、VCF、BCF、普通 gzip） | bcftools |
| `-o, --output FILE` | 输出文件；`reed stats` 中为输出前缀（与 vcftools `--out`、plink2 `--out` 同义，`--out` 作为别名） | bcftools / vcftools / plink2 |
| `-O, --output-type {v,z,b,u}[0-9]` | 输出格式与压缩级别；未给出时按扩展名推断，冲突时报错 | bcftools |
| `-W, --write-index[=csi\|tbi]` | 边写边建索引，默认 CSI | bcftools |
| `--threads INT` | 线程总数，默认全部核；无短参数 | bcftools / plink2 |
| `-r/-R` | 区域（需要索引） | bcftools |
| `-t/-T` | 区域（流式，`^` 排除） | bcftools |
| `-s/-S` | 样本子集（`^` 排除） | bcftools |
| `--keep FILE` / `--remove FILE` | 样本子集（每行一个样本名；plink2 的 FID/IID 两列格式也接受） | plink2 / vcftools |
| `--force-samples` | 不存在的样本只警告 | bcftools |
| `--memory MB` | 内存预算（MiB） | plink2 |
| `--tmp-dir DIR` | 溢写目录 | — |
| `--no-version` | 不写 `##reed_*` header 行 | bcftools |
| `--verbosity INT` | 0 静默 … 3 调试，默认 1 | bcftools |
| `-h/--help`、`--version` | — | — |

长度参数接受 `100000`、`100k`、`100kb`、`1M`、`1Mb`。

---

## 3. 过滤参数（`reed filter` 全部支持；`reed stats` 支持除 INFO 更新外的全部，用于边过滤边统计）

### 3.1 样本级

| 参数 | 语义（与来源完全一致） | 来源 |
|---|---|---|
| `-s/-S`、`--keep/--remove` | 静态子集 | bcftools / plink2 |
| `--mind FLOAT` | 样本缺失率 > FLOAT 的样本剔除。**plink2 语义**：在基因型掩码（§3.3）之后计算，分母是载入阶段（区域、`-v/-m/-M/--snps-only` 等载入期过滤）之后的全部位点 | plink2 |

注意：gold 管线的 keep 名单是在**原始**基因型上、全量统一算出来的，与 plink2 的 `--mind` 语义不同。
gold 对拍时，reed 与另外三个工具一样直接使用 `-S samples.keep.txt`（§8），而不是 `--mind`。

### 3.2 位点级（掩码前）

| 参数 | 语义 | 来源 |
|---|---|---|
| `-v/-V LIST` | `snps,mnps,indels,other,ref` | bcftools |
| `-m/-M INT` | 等位基因总数上下限 | bcftools（plink2 `--min-alleles/--max-alleles` 作为别名，语义相同） |
| `--snps-only[=just-acgt]` | 只保留 SNP | plink2 |
| `-f LIST` | FILTER 属于列表 | bcftools |
| `-i/-e EXPR` | bcftools 表达式（§4） | bcftools |
| `-q/-Q FREQ[:type]` | 等位频率下限/上限（`:minor`、`:alt`、`:nref`、`:major`） | bcftools view。**注意**：与 bcftools 一样，INFO 中有 AC/AN 时使用 INFO 值（本数据的 INFO/AC 是陈旧的）；需要按基因型重算请用 `--maf` |

### 3.3 基因型掩码

| 参数 | 语义 | 来源 |
|---|---|---|
| `--min-dp INT` / `--max-dp INT` | `FMT/DP < INT` / `> INT` 的基因型置缺失；**DP 缺失时不掩码** | plink2 `--vcf-min-dp/--vcf-max-dp`（`plink2_import.cc:901-923`，缺失值 `ScanInt32` 失败即通过），plink2 名字作为别名 |
| `--min-gq INT` | `FMT/GQ < INT` 置缺失；GQ 缺失时不掩码 | plink2 `--vcf-min-gq`（别名） |
| `--set-gts .`（仅 `reed filter`）配合 `-i/-e` 的 `FMT/` 表达式 | 与 `bcftools filter -S . -e EXPR` 完全相同 | bcftools `filter --set-GTs`。bcftools 的短参数 `-S` 在 view 中表示样本文件，reed 合并两个命令后把 `-S` 留给样本文件，这里只提供长参数（§6 注） |

不提供 vcftools 的 `--minDP/--maxDP/--minGQ`：vcftools 在 DP/GQ 缺失时把值当作 −1，于是**会掩码**
（`vcf_entry.cpp:540-577`），与上表相反。

### 3.4 位点级（掩码后，由基因型重算）

| 参数 | 语义 | 来源 |
|---|---|---|
| `--geno FLOAT` | 位点缺失率上限 | plink2 |
| `--maf FLOAT` / `--max-maf FLOAT` | 次等位频率 | plink2 |
| `--mac INT` / `--max-mac INT` | 次等位计数 | plink2 |
| `--hwe P [midp]` | 精确检验；`midp` 修正 | plink2 |

### 3.5 INFO 更新（仅 `reed filter`）

| 参数 | 语义 | 来源 |
|---|---|---|
| （默认） | 样本被 `-s/-S/--keep/--remove/--mind` 剔除时，按 bcftools view 的规则更新 `INFO/AC`、`INFO/AN`；只做基因型掩码时**不**更新（与 bcftools filter -S . 一致） | bcftools |
| `-I, --no-update` | 不更新 | bcftools |
| `--fill-tags LIST` | 重算并写入（AC、AN、AF、MAF、NS、HWE、F_MISSING…），定义与 `bcftools +fill-tags` 相同 | bcftools |
| `-x, --remove-info LIST` | 删除 INFO 标签 | bcftools annotate |

---

## 4. 表达式（`-i/-e`）

**语法与语义与 bcftools 相同**（bcftools 手册 "EXPRESSIONS" 一节），reed 实现其中的一个子集，超出子集的写法报
"not supported yet"，而不是用另一种语义执行：

- 支持：`QUAL`、`FILTER`、`ID`、`CHROM`、`POS`、`TYPE`、`INFO/TAG` 与裸写的 `TAG`（bcftools 规则：裸名先查 INFO）、
  `FMT/TAG`、`N_ALT`、`N_SAMPLES`、`AC`、`AN`、`AF`、`MAC`、`MAF`、`N_MISSING`、`F_MISSING`；
  比较 `== = != < <= > >=`、正则 `~ !~`；逻辑 `&& || & |`（bcftools 中两者对 FORMAT 的含义不同，reed 照搬）、`!`、括号；
  下标 `[0]`、`[*]`；缺失值 `"."`。
- 内建量 `AC/AN/AF/MAC/MAF` 的取值规则照搬 bcftools：**INFO 中存在时取 INFO 值，否则按基因型计算**。
  这正是本数据集的陷阱（INFO/AC 是陈旧的）。reed 不改语义，而是在检测到"表达式使用了 INFO 中已存在的 AC/AN"时打印警告，
  提示改用 `--maf/--mac`（plink2 语义，总是按基因型计算）或先 `-x INFO/AC,INFO/AN`。
- 与缺失值比较为假（bcftools 语义）：`-i 'QD>=4'` 会去掉没有 QD 的位点，`-e 'QD<4'` 会保留它们。
- 引用按基因型计算的量（`N_MISSING`、`F_MISSING`，以及 INFO 中不存在时的 `AC/AN/AF/MAC/MAF`）的表达式在掩码之后求值，
  等价于 bcftools 管线中"先 `filter -S .` 再 `view -i`"。

---

## 5. `reed stats`

```
reed stats [通用参数] [过滤参数] [--pop FILE] [窗口参数] <统计参数…> -o PREFIX <FILE|->
```

### 5.1 原则：参数只写科学方法，不写软件名【规】（用户裁决 2026-10-06 第三次）

- 统计参数按**方法本身**命名，帮助文本与日志给出**方法的发明者与文献**，不出现软件名。
- 经典统计量在各软件之间公式相同；差异只来自三件事：缺失数据处理、不变位点的分母约定、窗口网格（04-algorithms.md §2）。
  reed 对这三件事各只采用一种做法：
  1. **缺失数据**：Korunes & Samuk 2021（pixy）的逐位点计数法——π、dxy 用"差异数之和 / 可比较对数之和"，
     Watterson θ 与 Tajima's D 按逐位点实际调用数计算（Tajima's D 的方差按窗口平均调用数计算一次，即 pixy 2.x 采用的做法）。
     这是 reed 的唯一默认，不提供其他缺失数据模型。
  2. **不变位点**：默认要求输入含不变位点（all-sites VCF），并像 pixy 一样检查；只有变异位点时用 `--variant-only` 声明，
     此时窗口内未出现的位点视为"单态且所有样本都已调用"，计入分母。
  3. **窗口网格**：从位置 1 开始，`[1, w]`、`[1+s, w+s]`…（pixy/vcftools 的网格）。
- 无缺失数据、且窗口相同时，reed 的结果与 vcftools、scikit-allel、pixy 逐位一致（这条作为对拍门写进 06-roadmap.md）。

### 5.2 统计参数

| 参数 | 统计量 | 方法（日志中写出的引用） |
|---|---|---|
| `--pi` | 群体内核苷酸多样性 π | Nei & Li 1979；缺失数据：Korunes & Samuk 2021 |
| `--dxy` | 群体间绝对分化 dxy | Nei 1987；缺失数据：Korunes & Samuk 2021 |
| `--fst [wc84\|hudson]` | Fst（所有群体对），默认 `wc84` | `wc84`：Weir & Cockerham 1984；`hudson`：Hudson, Slatkin & Maddison 1992，按 Bhatia et al. 2013 的比值之和形式 |
| `--fst-pops A,B,C` | 多群体合并的一个 Fst | Weir & Cockerham 1984（r 个群体） |
| `--theta` | Watterson θ | Watterson 1975 |
| `--tajima-d` | Tajima's D | Tajima 1989 |
| `--pbs A,B,C` | PBS（A 为焦点群体），基于 Hudson Fst | Yi et al. 2010；`--pbs-normed`：Malaspinas et al. 2016（PBSn1） |
| `--sfs` + `--proj N[,N…]` | 位点频率谱（dadi / fastsimcoal2 文件），投影到每个群体 N 个单倍型 | 超几何投影：Marth et al. 2004；Gutenkunst et al. 2009 |
| `--freq` | 逐位点等位频率 | — |
| `--hardy [midp]` | 逐位点 HWE 精确检验、观测/期望杂合度 | Wigginton, Cutler & Abecasis 2005；`midp`：Graffelman & Moreno 2013 |
| `--het` | 逐样本纯合度与近交系数 F（矩估计） | — |
| `--missing` | 逐样本与逐位点缺失率 | — |
| `--depth` | 逐样本与逐位点平均深度 | — |
| `--summary` | QC 汇总报告（bcftools stats 文本格式，便于 MultiQC） | — |

`--fst` 只有两种取值，因为 WC84 与 Hudson 是两个不同的估计量（不同的科学方法），而不是同一方法的不同实现。

### 5.3 群体、窗口与输入声明

| 参数 | 语义 |
|---|---|
| `--pop FILE` | 两列 `sample<TAB>population`，无表头（`--populations` 作为别名） |
| `--window SIZE` | 窗口大小 |
| `--step SIZE` | 步长，默认等于窗口 |
| `--bed-file BED` | 用 BED 给定窗口 |
| `--sites-file FILE` | 只统计列出的位点 |
| `--variant-only` | 声明输入只有变异位点：π、dxy、θ 的分母补上窗口内未出现的位点（视为单态、全部调用） |
| `--bypass-invariant-check` | 输入不含不变位点、也不声明 `--variant-only` 时，强行按"只统计出现的位点"计算（与 pixy 同名参数的行为相同；会高估 π 和 dxy，日志中给出警告） |

Fst、PBS、Tajima's D 是比值统计量，单态位点对分子分母的贡献都是 0，所以不受 `--variant-only` 影响。

### 5.4 输出

`-o PREFIX`（也接受 `--out`）。

| 统计 | 文件 | 列 |
|---|---|---|
| π | `PREFIX_pi.txt` | `pop chromosome window_pos_1 window_pos_2 avg_pi no_sites count_diffs count_comparisons count_missing` |
| dxy | `PREFIX_dxy.txt` | `pop1 pop2 chromosome window_pos_1 window_pos_2 avg_dxy no_sites count_diffs count_comparisons count_missing` |
| Fst | `PREFIX_fst.txt` | `pop1 pop2 chromosome window_pos_1 window_pos_2 avg_wc84_fst\|avg_hudson_fst no_snps` + 分量列 |
| θ | `PREFIX_watterson_theta.txt` | `pop chromosome window_pos_1 window_pos_2 avg_watterson_theta no_sites raw_watterson_theta no_var_sites weighted_no_sites` |
| Tajima's D | `PREFIX_tajima_d.txt` | `pop chromosome window_pos_1 window_pos_2 tajima_d no_sites raw_pi raw_watterson_theta tajima_d_stdev` |
| PBS | `PREFIX_pbs.txt` | `pop1 pop2 pop3 chromosome window_pos_1 window_pos_2 pbs no_snps fst_12 fst_13 fst_23` |
| SFS | `PREFIX_sfs/dadi/*.sfs`、`PREFIX_sfs/fastsimcoal2/*.obs` | dadi / fastsimcoal2 标准格式 |
| freq / hardy / het / missing | `PREFIX.afreq`、`PREFIX.hardy`、`PREFIX.het`、`PREFIX.smiss`、`PREFIX.vmiss` | plink2 的列布局（社区最通用的表格格式） |
| depth | `PREFIX.idepth`、`PREFIX.ldepth.mean` | vcftools 的列布局（plink2 没有对应输出） |
| summary | `PREFIX.stats` | bcftools stats 文本格式 |
| 所有运行 | `PREFIX.log` | 命令行、版本、每个统计量的方法与文献、过滤漏斗、计数器 |

- 窗口统计沿用 pixy 的列布局：既然默认方法来自 pixy，下游已有的 pixy 绘图/汇总脚本可以直接使用。
- 坐标为 1-based 闭区间；缺失值写 `NA`。
- 表格布局只是文件格式约定，不涉及统计方法。

### 5.5 例子

```bash
# 默认：含不变位点的 VCF
reed stats allsites.vcf.gz --pop pops.txt --pi --dxy --fst --tajima-d --window 10k -o out --threads 56

# 只有变异位点的 VCF
reed stats clean.vcf.gz --pop pops.txt --variant-only --pi --fst --tajima-d --window 100k --step 50k -o out

# 逐位点/逐样本报告与 QC 汇总
reed stats clean.vcf.gz --freq --hardy --het --missing --summary -o qc

# 边过滤边统计
reed stats clean.vcf.gz --pop pops.txt --variant-only --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 \
  --pi --tajima-d --window 100k -o out
```

---

## 6. `reed filter`

```
reed filter [通用参数] [过滤参数] [--report FILE] <FILE|->
```

注：`bcftools filter -S/--set-GTs` 与 `bcftools view -S/--samples-file` 的短参数冲突。reed 把两个命令合并成一个，
因此 `-S` 保留给样本文件（view 语义），设置基因型用长参数 `--set-gts .`；配合 `-i/-e` 的 `FMT/` 表达式使用时，语义与
`bcftools filter -S . -e EXPR` 完全相同。

**求值顺序**：与 bcftools 管线 `view → filter -S . → view → +fill-tags → annotate -x` 一致，与 plink2 的阈值顺序一致。
对于互相独立的位点过滤，顺序不影响最终结果，只影响漏斗中每条规则的计数；会影响结果的只有"`--mind` 在掩码之后计算"这一点（§3.1）。

| 步 | 内容 |
|---|---|
| 1 | `-s/-S/--keep/--remove`、`-r/-R/-t/-T` |
| 2 | 载入期位点过滤：`-v/-V/--snps-only/-m/-M/-f` |
| 3 | 不依赖基因型的 `-i/-e`（按命令行顺序） |
| 4 | 基因型掩码：`--min-dp/--max-dp/--min-gq`、`--set-gts .` + `FMT/` 表达式 |
| 5 | `--mind` |
| 6 | `--geno`、`--hwe`、`--maf/--max-maf/--mac/--max-mac`（plink2 顺序），以及依赖基因型的 `-i/-e` |
| 7 | INFO 更新：`-x` → AC/AN 自动更新 → `--fill-tags` |
| 8 | 输出 |

**输出 header**：`##reed_version`、`##reed_filterCommand=<规范化命令>; Date=…`（仿 `##bcftools_viewCommand`）。

**gold 等价命令**（`gold/config/filters.conf` 默认值）：

```bash
reed filter clean.vcf.gz -S samples.keep.txt \
  --snps-only=just-acgt -M 2 -f PASS -i 'QUAL>=30' \
  -e 'QD<4 || FS>50 || MQ<40 || SOR>3' \
  --min-dp 5 --max-dp 100 --min-gq 20 \
  --geno 0.1 --maf 0.05 \
  -x INFO/AC,INFO/AN,INFO/AF,INFO/MLEAC,INFO/MLEAF,INFO/ExcessHet,INFO/InbreedingCoeff,INFO/DP \
  --fill-tags AC,AN,NS,MAF,HWE \
  -Oz -W -o clean.filtered.reed.vcf.gz --threads 56
```

---

## 7. 库 API 与 CLI【规】

- 删除：`FilterSpecV2` 的 TOML 反序列化、`read_spec_toml`、`compat::from_gold_conf`、`toml` 依赖、`--spec/--conf`、`REED_THREADS`。
- 保留：类型化 `FilterSpec`、`StatsRequest`（builder）。
- 新增：`reed_filter::expr`（bcftools 表达式子集的解析器、校验器、求值器）；`FilterSpec::to_args()`/`StatsRequest::to_args()`
  规范化渲染；往返性质测试 `parse(render(x)) == x`。
- clap 结构体、单位解析、错误渲染在 `reed-cli` crate。

---

## 8. gold 与对拍

- `gold/scripts/run_reed.sh`：与 `run_bcftools.sh` 相同的方式，从 `filters.conf` 生成 §6 的 reed 命令，并且和另外三个工具一样使用
  `run_all.sh` 统一算出的 `samples.keep.txt`。
- 统计对拍（新增 `gold/popgen/`）：软件只作为测试的参照物——默认方法对 pixy 2.2.4；`--variant-only` 的 π 对 vcftools `--window-pi`
  （两者定义相同，包括缺失数据）；无缺失数据的 fixture 上同时对 vcftools、scikit-allel、pixy 三家；PBS 对 scikit-allel；
  SFS 对 easySFS；hardy/freq/het/missing 对 plink2；summary 对 bcftools stats。见 04-algorithms.md §5、06-roadmap.md R7。

---

## 9. 退出码与消息

0 成功；1 运行期失败（带 `contig:pos` 定位）；2 用法错误（参数、表达式、标签、群体缺失、不变位点检查未通过）。
错误：`reed <cmd>: error: …`；警告：`reed <cmd>: warning: …`。

---

## 10. 从现有 CLI 迁移

| 现有 | 新 |
|---|---|
| `reed filter --spec f.toml` / `--conf filters.conf` | 参数（§6）；gold 用 `run_reed.sh` |
| `reed filter -o x.vcf.gz --write-index --compress-level 6` | `reed filter -Oz6 -W -o x.vcf.gz` |
| `reed stats -g pops --size 100000 --step 50000 --levels …`（默认计算 51 列） | `reed stats --pop pops --pi --fst … --window 100k --step 50k -o PREFIX`（显式选择统计量） |
| `REED_THREADS` | `--threads` |
| `reed-filter` 二进制 | 删除 |
