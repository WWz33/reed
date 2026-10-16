# 02 修订规范（N13-1…N13-11）

> 本文给出 `reed convert` 子命令（Phylo 多格式 + SFS 出口迁移 + 012/hapmap/msa
> + hapmap→VCF 反向）的**可执行规范**：语义模型（§1–§6，每条对照 vcf2phylip 2.9 /
> easySFS / vcftools / TASSEL / vcf2msa 实测语义）、并行执行模型（§7）、门面 API
> （§8）、验收标准（§9）、文档回改清单（§10）；N13-11 为反向目标的独立编号。
> 量尺：`docs/spec/cli.md` §0 总原则、`docs/spec/architecture.md` §7 融合式
> chunk 流水线、`docs/rounds/v5/05-algorithms.md` §3.1（SFS 出口形状）、
> `docs/rounds/v8`（线程预算纪律）、doc_v11（门面公开面纪律）。
> 标记沿用全仓约定：【规】必须；【荐】推荐；【决】需拍板。
> 缺陷与差异清单见 01-review.md（D13-1…D13-8 / B13-1…B13-12 / C13-1…C13-6）。

---

## §1 命令面：`reed convert`（N13-1）【规】

### 1.1 子命令与命令组

`Command` 枚举新增第四个变体；命令面按**转换目标**分组（cli.md §0.1"一条
命令就是全部规范"，参数名照抄参考工具）：

```
reed convert <INPUT> [--sfs [folded]] [--proj N[,N...]] [--sfs-int]        # SFS 组
             [--phylip] [--fasta] [--nexus] [--nexus-binary]               # Phylo 组
             [--outgroup NAME] [--min-samples-locus N] [--resolve-iupac]
             [--write-used-sites]                                          # Phylo 调节
             [--012]                                                       # 基因型矩阵组
             [--hapmap]                                                    # 基因型矩阵组
             [--msa] [--reference FILE] [--region CHR:B-E] [--flank N]
             [--collapse-suffix STR]                                       # MSA 组
             [--vcf] [-O v|z[level]] [--no-version] [--no-command]         # 反向组
             [通用参数：--threads/--memory/--verbosity/-o]
```

**门规则**【规】：

1. 六组（SFS / Phylo / 调节 / 012-hapmap / MSA / 反向）**至少给一个转换
   目标**；一个都没有 → usage 错误（退出码 2），报错文本列出可选项
   （"nothing to convert: give --sfs, a phylo format (--phylip/--fasta/
   --nexus), --012, --hapmap, --msa or --vcf"）。不落入"读全文件、写零
   文件"的空跑。
2. `--phylip` 是**opt-in 开关**，默认关闭——与 vcf2phylip 的 `-p`（默认写
   PHYLIP）**语义相反但拼写对齐**。理由：cli.md §0.7"参数在处理第一条记录
   之前全部校验"下，`reed convert` 的哲学是"要什么写什么"，暗写一个用户
   没要的文件违反最小意外；本条是 B13-1 的落地，对拍断言按"指定等价集"
   书写（01 §5-B13-1）。
3. `--outgroup/--min-samples-locus/--resolve-iupac/--write-used-sites` 不给
   Phylo 目标时 → usage 错误（"needs a phylo output"）。同 §1.3 的单义性
   纪律：拒绝无消费者的参数。
4. `--sfs` 必须给 `--proj`（沿用 stats 侧既有校验与报错文本）；`--proj`
   不带 `--sfs` → usage 错误。**若用户想跑 SFS 又顺手要 phylip**：允许
   `--sfs --proj 20 --phylip`（两目标并行执行，§7）。
5. `--msa` 必须给 `--reference FILE`（FASTA，须与 VCF 同参考）；`--reference`
   不带 `--msa` → usage 错。`--region`/`--flank`/`--collapse-suffix` 只在
   `--msa` 时接受（同第 3 条纪律）。
6. `--vcf`（hapmap→VCF）**独占运行**：它要求输入是 hapmap（`.hmp.txt`），
   而其余目标要求 VCF 输入；同给 `--vcf` 与任何 VCF 目标（或它们的调节
   参数）→ usage 错（"cannot mix --vcf (hapmap input) with VCF-input
   targets"）。输出走 filter 的 `-O` 语义：`-O v|z[0-9]`（`b`/`u` 报
   usage——BCF 记录模型对合成文本无收益）；不给 `-O` 时从 `-o` 扩展名推
   （`.gz` → z6），再不给 → 文本写 stdout。

### 1.2 SFS 侧：`--sfs` 与 stats 完全等价迁移【规】

`reed convert --sfs` 的一切语义与 `reed stats --sfs` 逐字节相同：参数
（`[folded]`、`--proj N[,N…]` 广播/逐群体、`--sfs-int`）、`--pop` 群体文件
语义、输出布局（`PREFIX_sfs/dadi|fastsimcoal2/…`）、`.log` 行、退出码。
因此：

- **`reed stats --sfs` 不删除、不隐藏**——保持既有脚本（`gold/popgen/
  run_easysfs.sh` 及用户工作流）不破坏；本命令是它的镜像入口。落地时用
  **共享翻译函数**（CLI 结构体 → `SfsRequest`）保证两个入口不漂移，并加
  往返测试锁死（§9.4）。
- 语义唯一事实源是 `docs/rounds/v5/05-algorithms.md` §3.1 与
  `crates/reed/src/stats/request.rs` 的 `SfsRequest`。本规范不重述其算法
  （超几何投影、folded、`m < n` 丢点计数），只迁移出口。

### 1.3 短字母不占用【规】

SFS 组照抄 stats 侧长名（`--sfs/--proj/--sfs-int`），**不占短字母**；
Phylo 组不引入 vcf2phylip 的短字母（`-i/-f/-n/-b/-r/-w/-p/-m/-o`——`-o`、
`-m` 已是 reed 全局参数，其余字母低频且 reed 无兼容包袱，cli.md §0.4 的
"短字母只用 bcftools 同义字母"不强制借 Python 脚本的短名）。`-m` 保留给
将来"最小等位计数"类 bcftools 同义参数。**`--012` 是数字开头的长名**
（照抄 vcftools 的 flag 拼写；clap 的长名允许数字开头，`--012` 直接可用）。

### 1.4 参数表（新 CLI 全量）

| 参数 | 语义 | 来源 |
|---|---|---|
| `--sfs [folded]` | 投影 SFS（等价 `stats --sfs`；`folded` 折叠） | easySFS / Marth 2004 |
| `--proj N[,N…]` | 投影大小（染色体数），一个值广播 | easySFS |
| `--sfs-int` | 谱值取整 | easySFS `--dtype int` |
| `--phylip` | relaxed PHYLIP（opt-in） | vcf2phylip |
| `--fasta` | FASTA 矩阵 | vcf2phylip |
| `--nexus` | NEXUS（DNA）矩阵 | vcf2phylip |
| `--nexus-binary` | 二元 NEXUS（SNAPP；仅二倍体 GT 参与） | vcf2phylip |
| `--outgroup NAME` | 该样本写成第一行 | vcf2phylip |
| `--min-samples-locus N` | 位点保留门（默认 4） | vcf2phylip |
| `--resolve-iupac` | 杂合随机拆解 | vcf2phylip |
| `--write-used-sites` | 写通过位点的坐标表 | vcf2phylip |
| `--012` | 0/1/2/-1 基因型矩阵三件套（§3.8） | vcftools `--012` |
| `--hapmap` | TASSEL hapmap `.hmp.txt`（§3.9） | TASSEL ExportUtils |
| `--msa` | 参考序列全位点 MSA，FASTA 输出（§3.10） | vcf2msa.py |
| `--reference FILE` | 参考基因组 FASTA（`--msa` 必选） | vcf2msa `-f` |
| `--region CHR:B-E` | MSA 区域（可重复；缺省全 contig） | vcf2msa `-r/-R` |
| `--flank N` | 区域两侧加侧翼碱基数（默认 0） | vcf2msa `-F` |
| `--collapse-suffix STR` | MSA 样本名按首个 STR 折叠（默认不折叠） | vcf2msa `.` 折叠（B13-7） |
| `--vcf` | hapmap 输入转回 VCF（record→record，§3.12）；配 `-O v\|z[0-9]` | TASSEL hmp→VCF |
| `-o, --output PREFIX` | 输出前缀（cli.md §2 既有语义） | reed |
| `--threads/--memory/--verbosity` | cli.md §2 既有语义 | reed |

## §2 转换模型：一次请求多目标（N13-2）【规】

### 2.1 `ConvertRequest`（门面类型）

```rust
/// One conversion target.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ConvertTarget {
    /// Projected SFS (SfsRequest + int dtype); semantics == stats --sfs.
    Sfs { request: SfsRequest, int_dtype: bool },
    /// Nucleotide alignments (any nonempty subset of phylip/fasta/nexus).
    Phylo { formats: PhyloFormats, options: PhyloOptions },
    /// Binary NEXUS (SNAPP), independent site stream and own counters.
    NexusBin,
    /// 0/1/2/-1 genotype-count matrix (vcftools --012, §3.8).
    CountMatrix,
    /// TASSEL hapmap .hmp.txt (§3.9).
    Hapmap,
    /// Reference-sequence whole-site MSA (vcf2msa, §3.10).
    Msa { reference: PathBuf, regions: Vec<Region>, flank: u32, collapse_suffix: Option<String> },
    /// A hapmap input converted back to VCF (§3.12): `v` (text) or `z`
    /// (BGZF at `level`). Exclusive with every VCF-input target (§2.4b).
    Vcf { compressed: bool, level: u8 },
}
```

`ConvertRequest { targets: Vec<ConvertTarget>, … }`——`targets` 非空由
驱动校验（§1.1）。**每目标独立**：计数器、`used_sites`、错误各自独立，
一目标的格式错误不牵连其他目标（除首错即停的全局契约，§7.3）。

### 2.2 SFS 目标 = stats 引用，不复制

`ConvertTarget::Sfs` 内嵌 `SfsRequest`；驱动侧转换是**对 stats 现有
`StatsJob` 的复用**（§7.5），不是新实现。收尾由既有
`reed::sfs::write_projected_sfs` 写文件——**零拷贝复用，零新代码**。
`Populations` 装载沿用 stats 侧（`--pop` 缺省 = 全体一组）。

### 2.3 Phylo 组：一个模型、多格式渲染

`--phylip/--fasta/--nexus` 三开关解析进 `PhyloFormats`（bitflags 风格
结构体）：

```rust
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PhyloFormats { pub phylip: bool, pub fasta: bool, pub nexus: bool }
```

**共享同一份逐位点模型**（列模型，§3）：一次列转换、三种**渲染**。
对拍断言：三格式由同一模型渲染 → `fasta` 序列 == `phylip` 行 == `nexus`
行（去填充后逐字符相等），作为内部不变式测试（§9.2）。

### 2.4 NEXUS-binary 与 Phylo 组的关系【规】

二元 NEXUS **不**进 `PhyloFormats`：它的位点门（仅双等位）与状态字母表
（0/1/2/?）都独立于核苷酸列模型（01 §1.4）。两个模型在 chunk 阶段各算
各的列，共享同一次行解析与位点门（`min_samples_locus` 两模型都先过，
01 §3（歧义点 4）已核实 vcf2phylip 的位点门先于两种模型）。`--outgroup`
对它生效（vcf2phylip 同款）。

### 2.4b 012 / hapmap / msa：三种新目标与共享行视图（D13-7）【规】

- `CountMatrix`（`--012`）与 `Hapmap` 是**位点行模型**（每通过位点一行
  或一值），与 Phylo 的**列模型**（每样本一列）转置关系不同——但三者
  共享同一份**行级 GT→等位索引解析**（D13-3 的原始索引视图，双等位判定
  在索引映射之前）。禁止出现第四份独立的行级 GT 解析（01 §6 守卫 11）。
- `Msa` 是**参考驱动**目标：行解析只发生在"该参考位置有 VCF 记录"时，
  无记录位点直接复制 REF（§3.10）；与其他目标共享行视图但**不共享位点
  门**（msa 无 min-samples-locus 概念，无记录=REF 是它的缺失语义）。
- 三者的 partial 形态：`CountMatrix` 为 per-sample `Vec<i8>` 列块（渲染
  时转置成行）；`Hapmap` 为**行缓存**（每通过位点一行 String/字节块，
  顺序 append——行模型天然 chunk 序拼接，无转置）；`Msa` 为区域切片的
  字节块。全部加法合并、chunk 序归约（§7.2 同款确定性论证）。
- **`Vcf` 是第五种形态，且与前四种根本不同**（D13-8）：它不是"折叠成
  矩阵再收尾渲染"，而是 **record→record 流式翻译**（hapmap 每行 → 一条
  VCF 记录），输出经 `-O` 的 `OrderedOutput`/`BlockWriter` 边归约边写
  （§7.7），不进 drain 的矩阵渲染。因此 `Vcf` 目标**独占一次运行**：
  它要求输入是 hapmap（不是 VCF），与 SFS/Phylo/012/hapmap/msa 互斥
  （后者都吃 VCF 输入）；同给两类 → usage 错（§1.1 门规则 6）。

### 2.5 样本子集与 GT 掩码

`reed convert` 携带 `--samples/--samples-file/--keep/--remove`（复用
`filter_flags::SampleArgs`），**不携带** `--mind`（stats 侧同款 usage 错
误 + 指示先 `reed filter --mind`）、不携带 `--min-dp/--max-dp/--min-gq`
（01 §5-C13-1【决】；对齐需要重开推敲，不塞进首发轮）。转换前 GT 不掩
码、子集后照实输出——vcf2phylip 没有 DP/GQ 概念，多给反而偏离参照。
**`--vcf` 除外**：hapmap 输入的样本集就是列头，无法子集化，给
`-s/-S/--keep/--remove` → usage 错（§1.1 门规则 6）。

---

## §3 Phylo 语义规范（N13-3）【规】

本节把 vcf2phylip 2.9 的行为钉成规范。凡"vcf2phylip 如此"处都有 01-review
§1 的行号对照；凡偏离处都在 B13 表里登记。

### 3.1 位点保留门（chunk 层）

按序判定，任一步不过即弃点并计数：

1. **非空、非 `#` 行**，tab 列数 == N+9（`is_anomalous`，vcf2phylip:82）；
   列数错 → 一次性警告后跳过（**行为**对齐 vcf2phylip 的打印跳过，不
   因一条坏行失败——cli.md §9 运行期错误是"带 contig:pos 的失败"，而
   vcf2phylip 对畸形行是宽容跳过；**错误哲学差异属 B13-2**，reed 有意
   对齐参照而非自家默认）。
2. **每样本基因型**：`./.` 开头计缺失；有基因型样本数 ≥ `--min-samples-locus`
   （默认 4；自动 min(N, m) 截断，vcf2phylip:229）。shallow 计数入 `.log`。3. **SNP 检验**（`is_snp`，vcf2phylip:90）：REF 长度 1，且 `ALT` 里逗号
   分隔的等位**个数 == 逗号数+1 且每个等位长度 1**（实现把 `<NON_REF>` 替换
   为 REF 后检验——GVCF 兼容；`-` 替换为 `*` 后进字母表）。非 SNP → mnp
   计数（`.log`）。
4. **列转换**（§3.3/§3.4）遇无法映射的等位（AMBIG 未命中）或等位索引越界
   → 该行按 malformed 跳过 + 计数（B13-9：参照在 AMBIG 未命中处**崩溃**，
   在索引越界处跳过；reed 统一为跳过）。

`.log` 写五个计数（对齐 vcf2phylip 的 stdout 报告，化为 reed 的 `.log`
通道，01 §1.6）：

| 计数 | vcf2phylip 原文 |
|---|---|
| `total_records` | Total of genotypes processed |
| `excess_missing` | Genotypes excluded … missing data |
| `mnps` | excluded for being MNPs |
| `accepted` | SNPs that passed the filters |
| `biallelic_for_bin` | Biallelic SNPs selected for binary NEXUS |

### 3.2 等位字母表

REF/ALT 每等位：**大写化、`-`→`*`、`<NON_REF>`→REF**（替换后进字母表；
`*` 参与 IUPAC 消歧，01 §3 歧义点 1）。字母表本身**不做 ACGTN 白名单预判**
——参照的拒绝机制是"GT 实际用到的等位组合不在 AMBIG 表 → 整行 malformed"
（`nt_dict` 收 REF/ALT 原文，01 §3 歧义点 3：REF 非纯 ACGT 的等位只有
出现在被保留样本的 GT 里才可能触发）。因此 reed 的列转换以"组合不在 AMBIG 表"（B13-9）与"等位索引超出
记录声明"为 malformed 判据，等位字母只做 `-`→`*`、大写化与
`<NON_REF>` 替换。

### 3.3 核苷酸列（phylip/fasta/nexus 共用）

对每个**保留样本**（子集后）：

1. GT 文本的首 `:` 字段，去掉 `/`、`|`，得等位索引串（如 `02`）。
2. 每索引查字母表得核苷酸字母集，**排序去重拼接**（如 `{G, *}` → `*G`）。
3. 查 AMBIG 表（IUPAC + `*` 组合的 62 项全表，规约进 reed-core 常量表
   `iupac.rs`；缺失 `./.` → `N`）得输出字母；`--resolve-iupac` 时改为
   **从等位索引串均匀随机取一**再查表（每样本每位点独立伯努利）。
4. 拼成列，逐位点 append。

**随机源**【规】：`--resolve-iupac` 用 reed 已有的确定性 seed 机制
（`StdRng::seed_from_u64`，seed 来自 `FilterOptions`/请求字段），**默认
seed = 0** → 同输入同输出（sha256 复现，§9.5）。`--seed` 参数**不加**
（01 §5-B13-3：确定性纪律，参数面最小）。文档必须写明"随机拆解是
确定性的、seed 固定 0"（B13-3 登记：与 vcf2phylip 的"每跑不同"有意
分歧，理由 = reed 的可复现是硬约束 architecture.md §5.10）。

### 3.4 二元 NEXUS 列（`--nexus-binary`）

- 只认可**双等位、二倍体、相位未知**的 GT 字面（`0/0,0|0→0`；
  `0/1,0|1,1/0,1|0→1`；`1/1,1|1→2`；`./.,.|.→?`）。相位已知不影响
  （vcf2phylip 的 GEN_BIN 表含 `0|0` 等），其余一切（多等位、非二倍体、
  缺失单倍型）→ `?`。
- 位点门：仅 `len(ALT)==1` 的记录参与（vcf2phylip:333 `len(record[4])==1`）。
- SNAPP 上游软件普遍不接受列内混 `?` 的全缺失列；本规范不做行级剔除
  （对齐 vcf2phylip）。

### 3.5 `--outgroup`

名字必须精确匹配保留样本名；`-o s1,s2` 取第一个（逗号/分号截断，
vcf2phylip:213）；不在样本集 → usage 错（**偏离** vcf2phylip 的静默忽略，
B13-4：cli.md §0.7 校验哲学，错比静默好）。效果：该样本行写在矩阵第一行
（phylip/nexus/fasta/bin.nexus 同步）。

### 3.6 输出布局

对拍固定为 vcf2phylip 实测布局（01 §1.5）：

| 格式 | 文件 | 头 | 行格式 |
|---|---|---|---|
| phylip | `PREFIX.phy` | `N L\n` | `name` + 右填充到 `max_name_len+3` + 序列 |
| fasta | `PREFIX.fasta` | — | `>name\n` 序列 |
| nexus | `PREFIX.nexus` | `#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=n NCHAR=L;\n\tFORMAT DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n` | 同 phylip 填充 |
| bin.nexus | `PREFIX.bin.nexus` | `…DATATYPE=SNP MISSING=? GAP=- ;\nMATRIX\n` | 同上 |
| 尾 | nexus/bin | `;\nEND;\n` | — |
| used_sites | `PREFIX.used_sites.tsv` | `#CHROM\tPOS\tNUM_SAMPLES\n` | 通过位点逐行 |

- `PREFIX` 默认：**输入文件名去扩展**（`in.vcf.gz` → `in`），由 `-o` 覆盖。
- phylip/nexus 的填充宽 = max name len + 3（vcf2phylip:389）。
- fasta 每样本两行，序列单行（整条一行，不折行）。
- 尾部 `;\nEND;\n` 是 vcf2phylip 的原文（`.bin.nexus` 亦然）。
- **临时文件与 `.minN` 后缀不移植**（D13-4）：`min_samples_locus` 直接是
  `--min-samples-locus` 的值，进 `.log`；reed 不在磁盘上做转置。

### 3.7 `--min-samples-locus` 与既有缺失语义的关系【规】

它是 vcf2phylip 语义的**整数样本数门**（vcf2phylip:100：有 GT 的**样本**
数，不是等位数）。与 reed 既有参数的关系在 `--help` 写明：
`--geno 0.8` ≈ `--min-samples-locus ⌈0.8·N⌉`（近似换算，分母语义不同：
`--geno` 的分母是子集后全部样本）。不提供换算、不提供别名——两个语义并
存于不同参照工具，reed 原样呈现（cli.md §0.2）。

### 3.8 `--012`：基因型计数矩阵（vcftools `--012` 语义，01 §1b.1）【规】

三个文件（vcftools 布局）：

| 文件 | 内容 |
|---|---|
| `PREFIX.012` | 矩阵：每保留样本一行；行首无索引列以外的内容（vcftools 的首列是样本下标，**保留**该布局：`<row_idx>\t<vals…>`）；每通过位点一值 |
| `PREFIX.012.pos` | 每通过位点一行 `CHROM\tPOS` |
| `PREFIX.012.indv` | 每保留样本一行样本名（矩阵行序） |

编码【规】：

| GT 形态 | 输出 | 对照 vcftools（vft-conv:445） |
|---|---|---|
| `0/0`（纯合 REF） | `0` | 相同 |
| `0/1`、`1/0`（杂合） | `1` | 相同 |
| `1/1`（纯合 ALT） | `2` | 相同 |
| `./.`（全缺失/无 GT） | `-1` | 相同 |
| 半调用 `0/.`、`./1` | **`-1`** | vcftools 为 `1`（B13-8：reed 按半信息=缺失） |
| haploid `0`、`1` | **`-1`** | vcftools 为 `0`/`1`（同 B13-8 口径：单倍型不构成二倍体计数；与 stats 侧 `HaploidAsHalfCall` 策略一致） |

- **位点门**：双等位（`n_alleles ≤ 2`，含 `ALT='.'` 的单等位记录 → 全 0
  列，对齐 vft-conv:436-440 的 `> 2` 判定）；`--min-samples-locus` 同门
  生效（vcftools 无此门，是 reed 的加叠—— Phylo 组不给时默认 1，即不
  加叠，保持 vcftools 裸行为；给了 `--min-samples-locus` 则四组共用同门，
  见 §3.1）。
- **`-1` 的语义**在 `.log` 写明（"half-calls and haploid calls emit -1;
  vcftools emits 1/0/1"——B13-8 的可见化）。
- **复用面**：计数矩阵可从既有 `SiteCounts` 之外按样本物化——但注意
  `--012` 需要**逐样本**值（矩阵每行一样本），`SiteCounts`（`crates/
  reed-core/src/count/site_counts.rs:20`）是位点聚合，不适用；逐样本编码
  从 D13-7 的多等位原始索引视图取（双等位判定先于索引映射），chunk 内
  append 进 `Partial` 的 per-sample `Vec<i8>` 列块（§7.2 同款）。

### 3.9 `--hapmap`：TASSEL hapmap（01 §1b.2）【规】

单文件 `PREFIX.hmp.txt`（tab 分隔）：

- **头行**：`rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\t
  assayLSID\tpanelLSID\tQCcode\t<sample…>`（11 固定列 + 样本名，
  ExportUtils.java:127-128）。
- **每位点一行**：
  - `rs#`：VCF ID；`.` → `chrom:pos`（与 reed 报告模式同款规则）。
  - `alleles`：**REF+ALT 声明集**，`/` 连接（多等位全列；B13-5——TASSEL
    写观察集，与 chunk 归约冲突）。
  - `chrom`/`pos`：原值；`strand` 恒 `+`；`assembly#`…`QCcode` 7 列恒
    `NA`（TASSEL 写侧同款，ExportUtils.java:176-177）。
  - 每样本一格：**TASSEL 读侧字母表内的单字母**（A/C/G/T/N/R/Y/S/W/K/M
    与 `-`，其余一律 → `N`，B13-10——实测 TASSEL 读不回小写软掩蔽码与
    三碱基码 V/H/D/B）。缺失 `./.` → `N`；半调用/haploid → 组合字母
    （如 `*A` → `N`，因小写码被收敛）。
- **位点门**：SNP 检验（§3.1 第 3 步）+ `--min-samples-locus`（同门）。
- **多等位**：`alleles` 列多值；样本格按 GT 实际等位组合查 AMBIG（组合
  无键 → malformed 跳行，同 §3.1 第 4 步）。
- 对拍：TASSEL **读回**验证（`BuilderFromHapMap` 读 reed 输出：导入 +
  再导出成功、site 数保持；`gold/popgen/run_convert_parity.sh` 的门 3），
  不 diff TASSEL 自产文件（其 alleles 列语义与命名都不同，B13-5）。

### 3.10 `--msa`：参考序列全位点 MSA（vcf2msa 语义，01 §1b.3）【规】

- **输入**：`--reference FILE`（FASTA）+ VCF；`--region CHR:B-E`（1-based
  闭区间，可重复），缺省全 contig；`--flank N` 侧翼。**区域约定**：
  vcf2msa 的区域端点是排他的（`-r chr1:1-10` 输出 9 个碱基，实测），
  reed 用自然闭区间——对拍时映射 `B-(E+1)`（01 §1b.3）。
- **输出**：**每 contig（或每区域）一个 FASTA 文件**：`PREFIX_<contig>
  [_<start>-<end>].fasta`，记录头 `><contig>[:start-end]_<sample>`，序列
  = 区域内**逐参考碱基**的样本序列（无记录处写 REF——与 Phylo 组的
  "只写变异位点"根本不同，帮助文本写明）。
- **逐位点规则**（对齐 vcf2msa:160-247 的次序）：
  1. 有 VCF 记录：按 `genotype_resolve` 语义取样本等位组合 → IUPAC；indel
     （等位长 > 1 或 `*`）走确定性填充（B13-6）：本位点各样本值取
     max len，短的以 `-` 补齐（vcf2msa 在 clustalo **缺失时直接崩溃**，
     其填充回退只在 MUSCLE 返回非零时才走到——reed 的填充路径是唯一
     可移植选择，01 §1b.3）。
  2. 无记录样本（GT 缺失或无该样本列）：写 REF 字母。
  3. 掩蔽（首版**不实现** mpileup 掩蔽——vcf2msa 的 `-m` 需 per-sample
     mpileup 文件，属另一输入面；`--min-dp` 的 VCF 内 DP 掩蔽同样不进
     （C13-1 同族）；后续按需立项）。
- **样本名折叠**：默认不折叠（B13-7）；`--collapse-suffix .` 显式开启
  （首个 `.` 前缀折叠，vcf2msa 同款）。
- **确定性**：区域划分只由 `--region/--flank` 决定；contig 输出文件按
  头文件 contig 序。
- **chunk 并行**：区域内的逐位点循环按 chunk 切（chunk 含 REF 段的
  内存映射切片——参考 FASTA 由 reader 线程一次性载入共享 `Arc<[u8]>`，
  worker 按区域偏移取 REF 段）；无记录位点零解析成本（直接 REF 复制）。

### 3.11 未纳入的 vcftools 转换格式（登记，不做）【规】

`--plink/--plink-tped`（PED/MAP、TPED/TFAM）：**【决】C13-5**，默认不进
首发轮（plink2 直读 VCF；`--012` + `.pos`/`.indv` 已覆盖行/列下标需求）。
`--IMPUTE/--BEAGLE-GL/--BEAGLE-PL/--ldhat/--ldhat-geno/--ldhelmet`：
**【决】C13-6**，本轮不做（相位/似然专门格式；IMPUTE/LDhat 的"整位点
零容忍丢弃"与 chunk 归约正交，需单独设计；06 §4 登记）。语义台账已在
01 §1b.4 钉死，后续立项直接引用。

---

### 3.12 `--vcf`：hapmap → VCF（反向转换，01 §1b.6）【规】

**输入**：hapmap 文本（`.hmp.txt`，11 固定列 + 样本列；单字母 IUPAC 或
双字母，按首行样本区字符数 ÷ taxa 数自动判定，对齐 TASSEL
BuilderFromHapMap.java:112-124）。**输出**：VCF 文本（`-Ov`）或 BGZF
（`-Oz`，级别同 filter 的 `-O`）。

**头**（reed 风格，B13-12）：

```
##fileformat=VCFv4.3
##reed_version=<crate version>   # 除非 --no-version
##reed_command=<argv>            # 除非 --no-version/--no-command
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  <samples...>
```

**不写 `##contig` 行**——TASSEL 自己的 hmp→VCF 也不写（C0 实测输出无任何
`##contig`），头因此在 open 时即完整可知，流式写出无需预扫。下游需要
contig 长度时走 `--contig-lengths`（与消费 TASSEL 导出的 VCF 同款约束）。
不写 TASSEL 那五个未使用的 FORMAT 声明（AD/DP/GQ/PL）与三个 INFO 声明
（NS/DP/AF）——B13-12。

**逐行翻译**（record→record，D13-8）：

| VCF 列 | 来源 | 规则 |
|---|---|---|
| `CHROM` | hapmap `chrom` | 原样转写；**不强制 sorted**（无 `##contig` 头后乱序不失真；下游索引工具自会拒绝，与消费 TASSEL 导出物一致） |
| `POS` | hapmap `pos` | 整数（非整数 → malformed 跳过，对齐读侧 :242-244） |
| `ID` | hapmap `rs#` | 原样（**不**像 TASSEL 重写为内部 siteName，B13-11） |
| `REF` | hapmap `alleles` 第 1 个 | `A/G` → REF=`A`（**声明序**，非 major-allele 猜测，B13-11） |
| `ALT` | `alleles` 其余 | `A/G` → `G`；`A/G/T` → `G,T`；单等位（`A`）→ ALT=`.`；`-`（gap）原样列出，GT 用其下标 |
| `QUAL` | — | 恒 `.` |
| `FILTER` | — | 恒 `PASS` |
| `INFO` | — | 恒 `.` |
| `FORMAT` | — | 恒 `GT` |
| 每样本 GT | hapmap 格 | 见下 |

**GT 编码**（核心，01 §1b.6 的 IUPAC 反查）：

1. 把该位点 `alleles` 列拆成等位字母表 `L`（声明序，如 `A/G` → `[A,G]`）。
2. 单字母格：缺失码只有 **`N` 与 `?`** → `./.`（TASSEL 写侧用 `N` 表缺失，
   C0 实测；`?` 是 GEN_BIN 族的缺失记号）。其余字母查 IUPAC 展开表
   （`ambig_expand`：`R`→{A,G}、`A`→{A}、`-`→{-}、`+`→{+}、`0`→{-,+}；
   小写按大写处理——比 TASSEL 读侧宽容）得到等位字母集合：
   - 集合的每个字母在 `L` 中的下标即等位编号；**纯合**（单字母集）→
     `i/i`，**杂合**（双字母集）→ 下标升序 `i/j`。
   - 展开超过两个字母（V/H/D/B，二倍体不可表示）或字母不在 `L`
     （`alleles` 列与实际格不一致）→ malformed 跳过 + 计数（B13-2 口径）。
3. 双字母格（`AG`/`AA`）：两字符各查单字母 → 两个等位编号 → `a/b`
   （TASSEL 读侧把第二字符当 allele1，:258-266；reed 对齐此序，输出
   非相位 `a/b`，因为 hapmap 双字母不携带相位）。
4. 相位：hapmap 无相位 → 一律 `a/b`（`/`），不写 `|`。

**`alleles` 为 `NA` 或缺列**（TASSEL 的无数据位点，:147-150）：REF=`.`、
ALT=`.`、全样本 `./.`（对齐 TASSEL 的 nAlleles==0 分支，ExportUtils:486-496）。

**往返保真契约**（B13-11 的可验证面）：`reed convert in.vcf --hapmap` →
`reed convert out.hmp.txt --vcf` 在**双等位 SNP、无缺失、REF 在 alleles
首位**的位点上，REF/ALT/GT 与原 VCF **逐记录相同**（`e2e_convert` 的
round-trip 测试锁定；多等位/缺失/indel 位点是 hapmap 的有损域，测试排除）。

**确定性**：逐行翻译无跨行状态（头在 open 时一次写全，无 `##contig` 累积），
与线程数无关；`-Oz` 的 deflate 级别固定（缺省 6，`-Oz<level>` 可改）。

## §4 SFS 迁移与 stats 的关系（N13-4）【规】

### 4.1 行为等价契约

`reed convert --sfs` ≡ `reed stats --sfs` **逐字节**（含 `.log`、退出码、
目录布局）。验收（§9.4）：

```
reed convert in.vcf.gz --pop p.txt --sfs --proj 20 -o A
reed stats  in.vcf.gz --pop p.txt --sfs --proj 20 -o B
diff -r A_sfs B_sfs   # 必须为空
```

### 4.2 stats 侧不删不改

`stats --sfs` 的解析、校验、报错全部保留（原样继续工作）。共享的方式是
**翻译层共享**：`translate/stats.rs` 里 SFS 参数块（`args.sfs/proj/
sfs_int` + pop 解析）提取为 `fn build_sfs_request(...)`，`translate/
convert.rs` 与 `translate/stats.rs` 共同调用；两个 CLI 结构体的字段名、
帮助文本一致（复制即漂移的教训：doc_v12 §0.5 元模式 1，AC/AN 三份实现
各自漂移——SFS 出口不做第三个）。

### 4.3 往返测试【规】

`ConvertCli → ConvertRequest → to_args()` 的往返测试锁死两个入口：
`convert --sfs` 的 args 渲染必须与 `stats --sfs` 的渲染**同名同序**
（`--sfs folded --proj 20,20 --sfs-int` 的子序列）。`SfsRequest` 的
`to_args` 若 stats 侧没有，则补上（对称 v11 的 `Levels`/`StatKind`
既有 `to_args` 面）。

---

## §5 错误与报告通道（N13-5）【规】

### 5.1 错误面

- 请求构造失败（门规则、`--proj` 缺失、outgroup 不在样本集、未知格式名）
  → `ErrorKind::Usage`（退出码 2）。
- 运行期失败（I/O、BcfNative 解码、维度爆内存预算、`--vcf` 的非 hapmap
  输入）→ `ErrorKind::Runtime`（退出码 1，带 contig:pos）。非 hapmap 输入
  在头扫描处报 `header error: not a hapmap file: …`，与 `filter`/`stats`
  吃到非 VCF 输入同口径（"给错格式的输入"是数据错误不是用法错误）。
- 畸形记录：**跳过 + 计数**（§3.1），一次性警告（对齐 v12 的
  "每 worker 一次"模式）。`--verbosity 0` 静默。
- 首错即停遵守 v9 W1/W2（错误先经 `FirstError` 记录、通道旁路禁止）——
  convert 的工作体是纯函数，错误传播走既有 chunk 管道语义。

### 5.2 `.log`（verbosity ≥ 1 时写 `PREFIX.log`）

```
command、version（既有节）
convert targets: sfs(unfolded,proj=20), phylip, nexus-binary
counters: total_records=… excess_missing=… mnps=… accepted=… biallelic_for_bin=…
outgroup: NAME (first taxon) / none
```

SFS 侧 `.log` 沿用 stats 现有行（`skipped m<n` 计数等），两目标合并写
进同一 `.log`（`targets` 行标注各自归属）。

### 5.3 无 `--report`

首版不写 JSON 漏斗报告（convert 不是 filter，无漏斗）。`.log` 是唯一
审计面。后续如 stats 侧要收编 convert 计数，再开报告（06 §4 不做）。

---

## §6 命名与文档（N13-6）【规】

| 对象 | 名字 | 依据 |
|---|---|---|
| 子命令 | `convert` | 用户命名 |
| CLI 结构体 | `ConvertCli` | 既有 `FilterCli/StatsCli/IndexCli` 模式 |
| 门面类型 | `ConvertRequest`、`ConvertTarget`、`ConvertReport`、`PhyloFormats`、`PhyloOptions` | v11 门面纪律 |
| 驱动函数 | `run_convert(path, &request, &budget) -> ConvertReport`、`run_convert_opened` | `run_stats` 对称 |
| 翻译模块 | `translate/convert.rs` | `translate/{filter,stats,index}.rs` 模式 |
| 核心模块 | `reed-core/src/iupac.rs`（AMBIG/GEN_BIN 常量 + 单测） | L3 纯函数 |
| 渲染模块 | `crates/reed/src/convert/writers.rs`（phylip/fasta/nexus/bin 渲染，纯函数 + 单测） | stats `writers/` 模式 |
| 模型模块 | `crates/reed/src/convert/mod.rs`（`ConvertJob` 等，crate 内 `mod convert`） | stats 模式 |

- **版本后缀禁令**（v2 05 §4）：类型不带 `V13`/`Convert2` 等后缀。
- 帮助文本：`--min-samples-locus` 的帮助引用 Ortiz 2019（DOI
  10.5281/zenodo.2540861）；`--sfs` 引用 Marth 2004 / Gutenkunst 2009
  （stats 同款）。**软件名只出现在测试与 gold 脚本里**（cli.md §0.3）。
- **署名与许可**【规】：vcf2phylip 是 GPLv3（工作区
  `../../vcf2phylip/LICENSE`）。reed 是 MIT/Apache 双许可——**不复制其
  代码**，只实现其行为（本规范即行为规格；AMBIG 表是事实数据不受版权
  保护，但表项的呈现以 IUPAC 标准的核苷酸编码为准，注明"per IUPAC
  nucleotide ambiguity codes"）。README 的 Credits 段落引 Ortiz 2019
  原文（对齐既有 pixy/easySFS 的引用方式）。CHANGELOG 的 convert 条目
  同样引 Ortiz 2019。

---

## §7 并行执行模型（N13-7）【规】

> 与 vcf2phylip 的单线程 readlines-chunk 循环（vcf2phylip:281）相比，
> reed convert 的价值主张就是并行。模型完全复用既有融合式 chunk 流水线
> （architecture.md §7），本节只钉"复用哪几块、多目标怎么并"。

### 7.1 总体形态

```
[reader 1 线程]  解 BGZF → 按 chunk 切（固定压缩字节，对齐 BGZF 块）
[worker N 线程]  每 chunk：一次行解析 → 各目标各自折叠列/谱
                 Partial = { phylo: Vec<u8>列块, bin: Vec<u8>列块,
                             sfs: SfsPartial(既有), counts: CounterVec }
[reducer 1 线程] 按 chunk 序合并 → 列块按序拼接；SFS 加法合并；
                 计数器求和；水位线无（无窗口）→ 全部收尾在 drain
[main]           收尾：序列拼接→ 渲染 → 写文件（4 类）
```

- **chunk 边界只由输入决定**（§7.2 原则 5 既有）→ 输出与线程数无关。
- worker 内每行**一次解析、多目标消费**：`LineRef` 只 parse 一次，SFS
  计数器与 phylo 列函数各取所需——**单次解析多消费者**，解析成本摊销。

### 7.2 Phylo 的 chunk 并行折叠

列模型天然可并行：每个 worker 把本 chunk 的通过位点列**按序**折叠成
**列块**（每样本一列的紧凑 `Vec<u8>`，或先写"位点 × 样本"的转置小块）。
Partial 结构：

```rust
pub(super) struct PhyloPartial {
    /// 本 chunk 的通过位点数
    n_sites: u64,
    /// 转置列块：per-sample 列字节（样本数 × n_sites）
    cols: Vec<Vec<u8>>,      // 仅 Phylo 目标
    /// bin 列字节（bin.nexus 目标；独立字母表）
    bin_cols: Vec<Vec<u8>>,  // 仅 NexusBin 目标
    /// used_sites 行（CHROM/POS/n_samples）
    used: Vec<UsedSite>,
    /// 计数器（§3.1 的五个）
    counters: ConvertCounters,
}
```

reducer 按 chunk 序 `concat`：列块 append（`Vec::extend`），`n_sites`、
计数器加法合并——**加法结合 + 顺序归约 ⇒ 与线程数无关的确定性**
（对齐 v6/02 §2.1 的位平面矩阵 partial 合并模式；无 seam：列转换不跨
chunk 边界依赖前文，与 v9 的 seam 分析一致，天然 seam-free）。

### 7.3 首错即停与预算

- worker 内单行错（I/O/解码）→ `EngineError` → 既有 abort/`FirstError`
  通道（v9 W1/W2）——**convert 不新增任何旁路错误出口**。
- `ThreadBudget`：reader 1 + workers N + reducer 1，全部在既有
  `ThreadBudget` 内记名（v8 D8 纪律：无匿名线程、无脱离预算的 rayon
  花园）。**多目标不乘线程**：`N` 个目标共享同一组 worker——chunk 内
  串行消费多目标（单次解析的摊销即在于此），不是每目标一个池。
- 内存：列块在途预算走既有 `MemBudget`（peak = 压缩 × MEM_FACTOR 的
  口径，v9 指标纪律）；输出序列在收尾渲染时一次性物化（≈ N_samples ×
  L 字节；对 3 M 位点 × 700 样本 ≈ 2.1 GB——超出典型预算时收尾渲染
  改流式逐行从列块拉取，**不整块物化**序列；首版即实现流式渲染，不
  先做整块再优化）。

### 7.4 确定性【规】

1. 同输入 + 同参数（含目标集）+ 同 seed（恒 0）⇒ 全部输出文件 sha256
   相同，threads ∈ {1, 4} 验证（§9.5）。
2. chunk 序归约 + 加法合并 ⇒ 列序 = 输入位点序，天然稳定。
3. `--resolve-iupac` 的随机消费按"位点 × 样本"唯一索引派生（列序号
   进 hash），与线程划分无关。

### 7.5 与 `StatsJob` 的复用关系

`ConvertJob` 实现 `ChunkJob`（`reed-engine::pipeline::run` 的入参），
组合 `StatsJob` 的 SFS 能力：

- **方案（定）**：`ConvertJob { sfs: Option<StatsJob>, phylo: Option<…>,
  bin: Option<…>, count: Option<…>, hmp: Option<…>, msa: Option<…> }`——
  SFS 目标在场时内部持有一个 `StatsJob`（仅 SFS
  kinds 的空窗请求），`run_chunk` 先调 `sfs.run_chunk(...)` 再做各目标
  折叠；`StatsPartial` 的 sfs 字段拆出嵌进 `ConvertPartial`。
  **不改 `StatsJob` 公开面**（stats 的调用方零感知）；convert 侧通过
  组合复用，stats 侧零改动。若实践中发现 `StatsJob` 的窗口字段不可
  精简复用（空窗仍造 grid），fallback 方案 B：把 SFS 折叠逻辑（`push_sfs`
  一族，run.rs:494-560）提为 `reed-stats` 的独立可组合单元（`SfsFold`
  trait），stats 与 convert 各自的 job 调它——**实现期二选一，规范不锁
  死组合方式，只锁死"单一折叠实现"**（D13-5 的守卫：`push_sfs` 不得
  出现第二份拷贝）。
- 012/hapmap/msa 三个新目标是 `ConvertJob` 的平级成员（不是 `StatsJob`
  的组合物）：`run_chunk` 内对每行**一次解析**后按目标消费（§2.4b），
  多目标共享 worker 池不乘线程（§7.3）。

### 7.6 收尾渲染与写出

drain 阶段（main 线程或 reducer 完成后）：

1. 合并后的列块 → 按格式渲染（`convert/writers.rs` 纯函数：
   `render_phylip(names, cols, …)` / `render_fasta` / `render_nexus` /
   `render_nexus_bin` / `render_012`（列块转置成行）/ `render_hapmap`
   （行缓存直写）/ `render_msa`（区域字节块直写））；outgroup 重排 =
   行序数组一次置换，不动物理列。
2. 文件写出：`PREFIX.phy/.fasta/.nexus/.bin.nexus` + `PREFIX_sfs/…`（SFS
   目标直接调 `reed::sfs::write_projected_sfs`）+ `PREFIX.used_sites.tsv`
   + `PREFIX.012`/`PREFIX.012.pos`/`PREFIX.012.indv`（CountMatrix）+
   `PREFIX.hmp.txt`（Hapmap）+ `PREFIX_<contig>[_<s>-<e>].fasta`（Msa，
   每 contig/区域一文件）。
3. `.log` 收尾写（verbosity ≥ 1）。

### 7.7 `--vcf` 的 record→record 输出路径（D13-8）

hapmap→VCF **不走** §7.1–§7.6 的"折叠 + drain 渲染"：它是流式翻译，
输出边归约边写，复用 filter 的 `-O` 机制：

```
[reader 1]  hapmap 文本按 1 MiB 切 chunk（PlainChunkSource 同款；头两行
            —— ##SAMPLE 注释 + 列头 —— 在 open 时扫描，得样本名与
            one/two-letter 判定；嗅探用的首行数据行回放给 source，
            不丢记录）
[worker N]  每 chunk：逐行翻译成 VCF 记录字节，worker-local BlockWriter
            （-Oz 时 deflate 成块；-Ov 时裸字节）→ ChunkOutput
[reducer 1] 按 chunk 序 OrderedOutput.accept_chunk（块序 = 输入序）
[main]      OrderedOutput.finish（EOF 块 / flush）
```

- **头先写**：`##fileformat`/`##reed_*`/`##FORMAT=GT`/`#CHROM` 行在
  open 时即完整可知（样本名来自 hapmap 列头；**无 `##contig`**，见 §3.12），
  `OrderedOutput.write_header` 一次写出，随后纯流式；头不随 chunk 累积。
- **复用面**：`reed_io::write::block_writer::BlockWriter`（deflate/裸字节）、
  `ordered_output::OrderedOutput`（块序写出 + EOF）、`ChunkJob`/`OrderedSink`
  （引擎）。**不新造压缩路径**（D13-8）。
- **独立驱动**：`run_convert` 检测到唯一目标是 `Vcf` 时走
  `convert::hapmap::run_hapmap_to_vcf`（自己的 job/sink），不与矩阵折叠的
  `ConvertJob` 混跑（§1.1 门规则 6 的互斥保证一次运行只有一类）。

---

## §8 门面 API（N13-8）【规】

`crates/reed/src/lib.rs` 增补（v11 纪律：根导出、请求按引用、可失败
`to_args`、无 `pub mod`）：

```rust
// convert 驱动与类型（08 §6.2 的 convert 对应物）
pub use crate::convert::{
    ConvertCounters, ConvertReport, ConvertRequest, ConvertTarget, MsaRegion,
    PhyloFormats, PhyloOptions, run_convert, run_convert_opened,
};
// `ConvertTarget::Vcf` 携带 `(compressed, level, output, command_line,
// no_version)`——filter 侧没有可复用的公开输出类型（它从 `-O` + 输出路径
// 在 CLI 翻译层推 format）。hapmap→VCF 首版只支持 `v`（文本）与 `z`
// （BGZF），不支持 `b`/`u`（BCF 记录模型对"从 hapmap 合成"无收益）；
// 给 `b`/`u` → usage 错（§1.1 门规则 6）。
```

- `run_convert(path: Option<&Path>, request: &ConvertRequest, budget:
  &ThreadBudget) -> Result<ConvertReport, Error>`——请求**按引用**（v11
  A11-2 教训直接沿用）。
- `ConvertReport { pub files: Vec<PathBuf>, pub counters:
  Option<ConvertCounters>, pub sfs: Option<SfsReport> }`——`counters`
  在无 Phylo 目标时为 `None`（SFS-only 跑法不造 Phylo 计数器）；
  `sfs` 直接复用既有 `SfsReport`（门面已有类型，零新增）。
- `to_args(&[...])`：`ConvertRequest` 渲染回 CLI 拼写（provenance 记录
  用），可失败（`Result`）——v11 A11-3 同款。
- **`to_args` 不调 `default_kinds` 类 panic 路径**（v11 A11-3 教训）：
  convert 无 kinds 概念，渲染函数从枚举纯分派。
- README「Library use」追加 `run_convert` 三行示例，**进编译门**
  （v11 的 README 示例测试纪律）。

## §9 验收标准（N13-9，全条目汇总）

| 门 | 内容 | 落点 |
|---|---|---|
| §9.1 语义对拍 | `gold/popgen/run_vcf2phylip.sh`：fixture + 全量（clean.vcf.gz）× 格式 × IUPAC/resolve/used_sites，diff 全等（resolve 跑法降级为结构断言，06 §5）；真值 = 工作区 `../../vcf2phylip/vcf2phylip.py` 2.9 | C2 |
| §9.2 内部不变式 | 三格式同模型：fasta 行 == phylip 行 == nexus 行（去填充）；bin 独立字母表 | C2 |
| §9.3 行为分档表 | 01 §5 的 reed↔参照工具差异表（B13-1…8），每条进 `cli.md`（B13 登记）| C2/C3 |
| §9.8 新格式对拍 | `gold/popgen/run_convert_parity.sh`：012 对 vcftools（fixture 半调用/杂合分行）、hapmap 对 TASSEL 读回验证、msa 对 vcf2msa（无 indel 区域 diff 全等 + indel 区域结构断言）| C2b |
| §9.9 反向对拍 | 同脚本门 6：VCF → hapmap → VCF 记录逐行相同（无损域）+ `-Oz` 帧 + TASSEL 读 reed 的 VCF（`-importGuess`，站点数保持）| C4 ✅ |
| §9.4 SFS 等价 | convert ≡ stats 逐字节（diff -r _sfs）；往返测试双入口同名同序 | C3 |
| §9.5 确定性 | threads {1,4} sha256 相同（含 resolve-iupac、多目标混合） | 每步 |
| §9.6 双 profile | `cargo test --workspace` debug+release 全绿；clippy -D warnings 无新告警 | 每步 |
| §9.7 性能 | 01 §8 **实测**：全量 158× vcf2phylip（28 核）且五文件逐字节；012/SFS 全量档低优先未测 | C3 ✅ |

## §10 文档回改清单（N13-10）

| 文件 | 回改 |
|---|---|
| `docs/spec/cli.md` | §1 命令表加 `reed convert` 行；§7.2/§7.4 的 SFS 行加"或 `reed convert`"；新增 §11「`reed convert`」全参数表 + B13 差异表 |
| `README.md` / `README.zh-CN.md` | 命令一览 + Credits 引 Ortiz 2019 + Library use 的 run_convert 示例 |
| `CHANGELOG.md` | `### Added`：`reed convert`（SFS/Phylo/多格式、并行模型、等价声明） |
| `docs/README.md` | 当前规范表加 doc_v13 行 |
| `ci/check-doc-drift.sh` | 追加守卫（01 §6 清单 1–10） |
| `docs/spec/architecture.md` | L2 分层图加 convert 驱动模块行（§7.5 的组合关系） |

## §11 决策默认表

| 决策 | 默认 | 待拍板 |
|---|---|---|
| C13-1 GT 掩码族（`--min-dp/--max-dp/--min-gq`）进 convert | **不进**（最小门；GT 掩码语义与"列字母表"无交集，进则引 §3.1 门外的一层掩码模型） | §2.5 |
| C13-2 `--min-samples-locus` 是否复用 `--geno` 语义 | **不复用**（整数样本数门，vcf2phylip 原语义；帮助注明与 `--geno` 的近似换算） | §3.7 |
| C13-3 `--pop` 是否必选（SFS 目标时） | **沿用 stats**：缺省全体一组，不强制 `--pop` | §1.2 |
| C13-4 msa 样本名是否折叠 `.` 后缀 | **不折叠**（vcf2msa 的折叠是作者管线约定；`--collapse-suffix` 显式开启） | §3.10 |
| C13-5 `--plink/--plink-tped` 是否进首发轮 | **不进**（plink2 直读 VCF；012 已覆盖下标需求） | §3.11 |
| C13-6 `--IMPUTE/--BEAGLE-*/--ldhat*/--ldhelmet` 是否实现 | **本轮不做**（相位/似然专门格式，后续立项） | §3.11、06 §4 |
| B13-1 `--phylip` opt-in（vs vcf2phylip 默认写） | **opt-in**（对拍断言按指定等价格式集书写） | §1.1 |
| B13-2 畸形行跳过不报错退出 | **跳过+计数**（对齐参照；cli.md §9 运行期错误保留给 I/O/解码） | §3.1 |
| B13-3 resolve-iupac 固定 seed | **seed=0 确定性**（sha256 门；文档写明） | §3.3 |
| B13-4 outgroup 不在样本集 → usage 错 | **报错**（cli.md §0.7 校验哲学） | §3.5 |
| B13-5 hapmap `alleles` 列取声明集 | **REF+ALT 声明**（观察集需全体样本频率排序，与 chunk 加法归约冲突；TASSEL 读侧只按 `/` 切分不受影响） | §3.9 |
| B13-6 msa indel 的 gap 摆放 | **确定性填充**（不引 ClustalO 外部进程；与 vcf2msa 的 MUSCLE 失败回退一致） | §3.10 |
| B13-7 msa 样本名折叠默认 | **不折叠**（= C13-4 的行为面） | §3.10 |
| B13-8 `--012` 半调用/haploid 输出 | **`-1`**（vcftools 输出 `1`；reed 按 GtCode 半信息=缺失口径，`.log` 写明） | §3.8 |
| B13-11 hmp→VCF 的 REF/ID 来源 | **照抄 hapmap 列**（`alleles[0]`→REF、`rs#`→ID、GT 按声明序）；不学 TASSEL 的 major-as-REF + 内部 siteName 重写 | §3.12 |
| B13-12 hmp→VCF 的头 | `VCFv4.3` + 仅 `##FORMAT=GT`（不写 TASSEL 的五个空 FORMAT 声明） | §3.12 |
| B13-9 AMBIG 表未命中 | **跳过 + 计数**（参照在表外字母处**崩溃**（KeyError 未捕获）、索引越界处跳过；reed 统一） | §3.1/§3.3 |
| B13-10 hapmap 格字母表 | **收敛到 TASSEL 读侧集合**（A/C/G/T/N/R/Y/S/W/K/M/-，其余 → `N`；参照写侧自己也读不回小写/三碱基码） | §3.9 |
