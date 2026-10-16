# 01 语义对照与差异台账（reed convert vs vcf2phylip 2.9 / easySFS / TASSEL）

> doc_v13 是**新功能轮**（与前十二轮的缺陷修复性质不同）：为 `reed` 增加
> `convert` 子命令——vcf2phylip 的 Phylo 多格式转换（并行化）+ SFS 出口从
> `stats --sfs` 的等价迁移；第二轮增补 `--012`/`--hapmap`/`--msa`，第三轮
> 增补 hapmap→VCF 反向（`--vcf`）。本文是**语义对照台账**：把各参照工具的
> 行为逐条钉死（§1–§4，全部带源码行号【据】），列出 reed 的有意分歧
> （§5 分档表：D13 缺陷类 / B13 有意分歧类 / C13 待拍板类）与守卫清单
> （§6）。修法与 API 规范见 02，实施顺序见 06。
>
> 参照物：`../../vcf2phylip/vcf2phylip.py`（Ortiz 2019，v2.9，2023-07-07，
> GPLv3——**只实现行为，不复制代码**，见 02 §6）；easySFS（已对拍门
> `gold/popgen/run_easysfs.sh` 既有）；vcftools（`--012`）、TASSEL 5
> （hapmap 双向）、vcf2msa.py（第二轮/第三轮增补，§1b）。全部结论来自
> 源码通读、既有对拍脚本与参照工具实测（C0 基线）。

---

## §0 为什么是这两个工具、为什么并成一个子命令

- **vcf2phylip** 是系统发育分析的标准前处理（Ortiz 2019，DOI
  10.5281/zenodo.2540861）：VCF → PHYLIP/FASTA/NEXUS/二元 NEXUS（SNAPP）。
  它是**单线程** Python（`readlines(50000)` 逐块，vcf2phylip:281），README
  自述 20 GB VCF ~27 分钟——reed 的 chunk 并行流水线正是此负载的形状。
- **SFS 出口**已在 `reed stats --sfs` 落地（doc_v5 §3.1，对拍 easySFS 全
  绿），但出口挂在一个"窗口统计"命令下面语义别扭：SFS 是 genome-wide 的
  格式转换（VCF → dadi/fsc2 谱文件），与"转换"同类。用户指令（2026-10-16）
  明确：sfs convert 迁移到新子命令。
- 两者共享同一执行形状：**读一遍 VCF、按位点折叠、收尾写一组格式文件**
  ——不是 filter 的"逐记录变换输出"，也不是 stats 的"窗口聚合"。这正当化
  一个新的命令面，而不是把 phylip 塞进 stats（stats 的参数表已经 60+ 行）。

## §1 vcf2phylip 语义台账（【据】全部为 vcf2phylip.py 行号）

### §1.1 参数与门

| 行号 | 行为 |
|---|---|
| :150-215 | argparse：`-i` 必选；`--output-folder`（默认 `./`，不存在则建）；`--output-prefix`（默认输入名去扩展）；`-m/--min-samples-locus`（默认 4）；`-o/--outgroup`（默认空）；`-p` 关闭默认的 PHYLIP；`-f/-n/-b/-r/-w` 开关 |
| :213 | outgroup 取 `,` 与 `;` 的第一个切分（`"a,b"` → `a`） |
| :229 | `min_samples_locus = min(num_samples, m)`：m 超过样本数自动截断 |
| :241 | 输出前缀追加 `.min<N>` 后缀（N = 生效的 min_samples_locus） |

### §1.2 输入与位点门（主循环 :281-341）

| 行号 | 行为 |
|---|---|
| :67-79 | 样本名从 `#CHROM` 行取第 10 列起，`./` 后缀剥掉（`record[i].replace("./", "")`） |
| :82-87 | `is_anomalous`：列数 != 样本数+9 的行**打印 "Skipping malformed line" 后跳过**（不退出） |
| :90-97 | `is_snp`：REF 长 1 且 `alt.count(",") + 1 == len(alt) - alt.count(",")`（即每个逗号分隔等位长 1）；`<NON_REF>` 先替换为 REF（GVCF） |
| :100-108 | `num_genotypes`：`record[i]` 以 `.` 开头（`./.`、`.|.`、`.`）计缺失；返回有 GT 的**样本数**（不是等位数） |
| :281 | `readlines(50000)`：50 000 行一块读（大文件优化，非语义） |
| :302 | `num_samples_locus < min_samples_locus` → shallow 计数、弃点 |

### §1.3 列转换（核苷酸）

| 行号 | 行为 |
|---|---|
| :111-133 | `get_matrix_column`：`nt_dict`：`0`→REF（`-`→`*`、大写化）、`.`→`N`、ALT 逗号切分逐个入表（`-`→`*`、`<NON_REF>`→REF） |
| :117-119 | GT 取首 `:` 字段，去 `/` 和 `|`，等位索引→字母，`sorted(set(...))` 拼接 |
| :122 | 键不在 AMBIG（未知核苷酸组合）→ 整行返 `"malformed"` → 打印跳过 |
| :126-128 | `resolve_IUPAC=False`：查 AMBIG 表；`True`：`random.choice(geno_num)` 均匀取一等位索引再查表 |
| :29-47 | AMBIG 表 62 项：IUPAC（M/R/W/S/Y/K/V/H/D/B/N）+ `*` 参与的组合（小写 = 有 N 或 `*` 参与）；`"*"` 单独 → `"-"`；`ACGT` 全 → `N` |

### §1.4 列转换（二元 NEXUS）

| 行号 | 行为 |
|---|---|
| :49-61 | GEN_BIN：`0/0,0|0→0`、`0/1,0|1,1/0,1|0→1`、`1/1,1|1→2`、`./.,.|.→?`；**其余一切 GT 字面 → `?`** |
| :333 | 只有 `len(record[4]) == 1`（ALT 恰一字符 = 双等位）的记录进 bin 流（计数 snp_biallelic） |

### §1.5 输出（转置与写出 :346-498）

| 行号 | 行为 |
|---|---|
| :246-256 | 先写"竖排"临时文件（`.tmp` / `.bin.tmp`），转置阶段逐行取第 s 列（**O(N·L) 的重复 IO**，见 D13-2） |
| :370-371 | phylip 头：`"{n_samples} {n_accepted}\n"`（relaxed 版，无交错块） |
| :377-379 | nexus 头：`#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=n NCHAR=L;\n\tFORMAT DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n` |
| :382-384 | bin.nexus 头：同上但 `DATATYPE=SNP MISSING=? GAP=-` |
| :412 | 名字列填充宽 = 最长样本名 + 3（首次定义） |
| :395-407 | outgroup 样本先写（三格式 + bin 各自转置一遍临时文件） |
| :486/:490 | 尾：`;\nEND;\n` |
| :429/:455 | fasta：`>name\n` + 序列单行 |
| :495-497 | 删临时文件（:495 `.tmp`、:497 `.bin.tmp`） |
| :354-357 | used_sites：`#CHROM\tPOS\tNUM_SAMPLES\n` + 通过位点行（NUM_SAMPLES = 该位点有 GT 的样本数） |

### §1.6 stdout 统计（:344-353）

五计数（进入 reed 的 `.log`，02 §3.1）：total / excess-missing / mnp /
accepted / biallelic-for-bin；每 500 000 行进度打印（:294-295）。

## §1b 其他参照工具语义台账（第二轮增补，2026-10-16）

> 用户点名 `--012`（vcftools）、hapmap（TASSEL）、msa（vcf2msa.py），源码
> 均在工作区。以下行号：【v】= `../../vcftools/src/cpp/variant_file_format_convert.cpp`
> （下称 vft-conv）、【vft-set】= `vcf_entry_setters.cpp`、【vft-par】=
> `parameters.cpp`、【t】= `../../tassel/src/main/java/net/maizegenetics/dna/snp/`
> （下写 `ExportUtils.java` / `io/BuilderFromHapMap.java` 相对 t）、【m】=
> `../../vcf2msa.py/vcf2msa.py`。

### §1b.1 vcftools `--012`（0/1/2 矩阵）【据】

| 位置 | 行为 |
|---|---|
| vft-par:183 | flag 名就是 `--012`（`output_012_matrix`） |
| vft-conv:377-384 | 三个文件：`PREFIX.012`（矩阵）、`PREFIX.012.pos`（CHROM\tPOS 每行一位点）、`PREFIX.012.indv`（每行一样本名） |
| vft-conv:436-440 | **双等位门**：`N_alleles > 2` → one_off_warning（"Only outputting biallelic loci"）+ 跳过记录（`> 2` 意味着**单等位 ALT='.' 的记录保留**，输出全 0 列） |
| vft-conv:445 | GT pair `(-1,-1)`（双缺失或无 GT）→ `-1`；`(0,0)`→`0`；`(1,1)`→`2`；**其余一切（含半调用、haploid）→ `1`** |
| vft-conv:396-400 | **每样本一个临时文件**（mkstemp），矩阵行以样本为行主序拼出——与 vcf2phylip 同型的"临时文件转置"，D13-2 同款禁止 |
| vft-set:68-95 | GT 解析：无 `/`\|`|` → haploid（ploidy=1，第二等位记 `.`→-1）；多倍体直接 error |
| vft-conv:441 | `.pos` 行在**双等位门之后**写——被跳过的多等位位点不进 `.pos`（矩阵列与 pos 行严格对齐） |

**关键语义坑**（C0 实测确认）：半调用 `0/.` 的 pair 是 `(0,-1)` → 落
`1`（与杂合同码）；`./.` 是 `(-1,-1)` → `-1`；haploid `1` 是
`(1,-1)` → `1`、haploid `0` 是 `(0,-1)` → `1`（**注意 haploid REF 也落
1**，不是 0）。**reed 规范见 02 §3.8（B13-8：半调用/haploid 一律 `-1`）。**
另：vcftools 的 `.012` 行首是**原始 header 下标**（`ui`），子集后不重编
（reed 对齐此行为，02 §3.8 已写明）。

### §1b.2 TASSEL hapmap（`.hmp.txt`）【据】

| 位置 | 行为 |
|---|---|
| ExportUtils.java:127-128 | 11 个固定列头：`rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode`，然后样本名（tab 分隔） |
| ExportUtils.java:153-174 | 每位点一行：`rs#`=位点名；`alleles`=**按基因型实际出现的等位**排序去重以 `/` 连接（可为 1 或 3+ 个；无数据时写 `NA`） |
| ExportUtils.java:176-177 | `chrom`=contig 名、`pos`=物理位置；`strand` 恒 `+`；`assembly#`…`QCcode` 7 列恒 `NA` |
| ExportUtils.java:180-199 | 每样本一格：默认（diploid=false 分支）**单字母 IUPAC**（`diploidAsString`）；`diploid=true` 时**双字母**（如 `AG`） |
| io/BuilderFromHapMap.java:53-57 | 读侧 `NUM_HAPMAP_NON_TAXA_HEADERS = 11`；列序 `rs#(0) alleles(1) chrom(2) pos(3)` |
| io/BuilderFromHapMap.java:258-266 | 读侧双字母按 `[i+1] [i]` 顺序入核（**忽略相位**）；单字母按 IUPAC 展开为 diploid |

**关键语义坑**（C0 实测确认）：TASSEL 的 `alleles` 列是**从基因型推出的
观察等位**，且排序是**按频率**（实测 chr1:3 的 `G/A` 而非 REF/ALT 序）；
写侧用 `AlleleFreqCache.allelesSortedByFrequencyNucleotide` 全体样本一起
排。reed 的写侧规范（02 §3.9）取**REF+ALT 声明集**（多等位全列）而非观察
集：差异登记 B13-5（观察集依赖全体样本物化，与 chunk 流水线的确定性归约
冲突）。`rs#` 实测是 TASSEL 内部命名（`S1_1` 之类），reed 取 `chrom:pos`
（ID 为 `.` 时）或 VCF ID——写侧命名无对拍真值，读回验证只查 taxa/site
数与可解析性。**读侧字母表更窄**（实测）：`getNucleotideDiploidByte` 的
哈希只有 A/C/G/T/N/R/Y/S/W/K/M/-/+/0/X/Z 等**大写**项——小写软掩蔽码
（`a`/`m`）与三碱基码（V/H/D/B）**读回即报错**（"illegal value"）。reed
的 hapmap 格因此**收敛到该子集**（其余 → `N`，02 §3.9），这是对 vcf2phylip
AMBIG 表的有意收窄。另：`strand` 写 `+`；`assembly#`…`QCcode` 7 列写
`NA`（TASSEL 自己的导出同款）。

### §1b.3 vcf2msa.py（参考序列 + 全位点 MSA）【据】

| 位置 | 行为 |
|---|---|
| m:56-79 | **必须参考 FASTA**（`-f`）+ VCF（`-v`，bgzip+tabix）；**逐参考碱基**循环（`for nuc in range(spos, epos)`），非逐记录 |
| m:155-170 | 样本名按第一个 `.` 前缀折叠（`SampleID.variant` → `SampleID`）；每 contig 输出 `contig_<name>.fasta`，记录头 `>contig:start-end_sample` |
| m:160-175 | 每位点：先查 VCF 有无该位置记录；**无记录的样本写 REF**（从 VCF 的 REF 或参考序列取）——与 vcf2phylip 的"只写变异位点"根本不同 |
| m:424-476 | `genotype_resolve`：多等位取**排序去重拼接** → `reverse_iupac_case`（AG→R 等，保留大小写）；indel（长度>1 或 `*`）单独分支：`--indel` 时优先取 indel、多个 indel 取**最短** |
| m:479-517 | `reverse_iupac_case`：IUPAC 表含大小写两份（小写=软掩蔽保留） |
| m:179-247 | 长度不齐（indel）时调 **ClustalO 外部进程**对齐；失败回退原值；掩蔽样本写 `N` 重复至 maxlen、缺失样本写 `-` 重复 |
| m:100-115 | mpileup 掩蔽（可选）：per-sample mpileup 深度 < `--cov` 的位置写 `N` |
| m:138-141（实测） | **区域端点是排他的**：`-r chr1:1-10` 输出 9 个碱基（`range(spos, epos)` 的 `epos = end−1`）；reed 的 `--region B-E` 是自然闭区间（输出 E−B+1 个），对拍映射为 vcf2msa `B-(E+1)` |
| m:198-207（实测） | indel 位点调 **clustalo**；缺失时 `FileNotFoundError` **不被捕获**（只捕 `ValueError`）→ 直接崩溃。B13-6 的"确定性填充"因此不是可选优化而是**唯一可移植路径** |
| m:649-657 | 参数：`-f/-v/-c/-m/-d/-r/-R/-F/--indel/--force` |
| m LICENSE | **GPLv3**——与 vcf2phylip 同款"只实现行为不复制代码"约束 |

**关键语义坑**：msa 是**全位点（含单态/参考位）**输出，每 contig 一文件、
样本序列长度 = 区域长度（含 indel 对齐后的 gap）。ClustalO 外部依赖是
对齐 gap 摆放的唯一裁决者——reed 首发不移植外部对齐（D13-7：indel 位点
的 gap 摆放用确定性规则——ALT 等位字面长度差以 `N`/`-` 填充对齐 maxlen，
即 m:247-268 的 MUSCLE 失败回退路径），B13-6 登记。vcf2msa 依赖随机
region 遍历（`vfh.fetch`），无确定性隐患；但样本名 `.` 折叠语义特殊
（B13-7 登记，reed 默认**不折叠**，给 `--collapse-suffix .` 显式开关
【决】C13-4）。

### §1b.4 vcftools 其余转换格式盘点（同文件，顺手钉死）

| flag | 文件 | 位点门 | 每格语义 | 备注 |
|---|---|---|---|---|
| `--plink`（vft-conv:10） | `.ped` + `.map` | 双等位（>2 跳过+警告，:108-112） | 等位字母两列/样本（缺失 `0 0`；haploid 第二列按相位补同字母或 0，:151-167）；ped 前 6 列 FID IID 0 0 0 0（:46-47）；chrom 映射 chrN→N/X/Y/XY/MT（:54-66，`--chrom-map` 覆盖） | PLINK 1 经典格式；plink2 已直读 VCF，**优先级低** |
| `--plink-tped`（:219） | `.tped` + `.tfam` | 同上 | 转置：每位点一行，前 4 列 chrom rs-id cM pos 后跟 2N 等位列 | 同上 |
| `--IMPUTE`（:505；vft-par:255） | `.impute.legend` + `.impute.hap` + `.impute.hap.indv` | 双等位 + **全相位 + 零缺失**（任一样本缺失/未相位→整位点丢弃，:546-570） | legend：`ID pos allele0 allele1`（ID `.`→`chrom-pos`）；hap：每样本两列 0/1 | 相位数据专用 |
| `--BEAGLE-GL/--BEAGLE-PL`（:1267；vft-par:184-185） | `.BEAGLE.GL`/`.BEAGLE.PL` | 双等位 + FORMAT 必须有 GL/PL 标签 | 头 `marker alleleA alleleB` + 每样本 3 列（GL 取 10^x）；marker=`chrom:pos` | 似然导出，非 GT 转换 |
| `--ldhat`（:616；vft-par:274，**自动 phased_only**） | `.ldhat.sites` + `.ldhat.locs` | 双等位；**必须 --chr 单染色体**（vft-par:717） | sites 头 `2N n_sites 1`；每 haplotype 两行（`>name-k`）；碱基 0/1/?；locs 头 `n_sites max_pos/1000.0 L` + 每位点 pos/1000.0（4 位小数） | LDhat rhos 需要位置 |
| `--ldhat-geno`（:799） | 同上 | 同上 | sites 头 `N n_sites 2`；每样本一行；**genotype 编码**：00→0、11→1、01/10→2、haploid(相位`\|`)→0/1、缺失→`?`（:905-937） | 注意与 --012 的编码**不同**（这里杂合是 2） |
| `--ldhelmet`（:983；vft-par:275，自动 phased_only+remove_indels） | `.ldhelmet.snps` + `.ldhelmet.pos` | 双等位 + indel 剔除 | snps **无计数头**，直接 `>name-k` + 碱基行（等位字母 ACGT，缺失 `N`）；pos 每行 POS（整数，4 位小数格式） | — |

**结论**（进 06 的排期分档）：`--012` 进首发（用户点名、无外部依赖、
与 GtCode 语义对齐）；hapmap 进首发（用户点名）；msa 进首发但 indel
对齐降级为确定性填充（B13-6）；`--plink/--plink-tped` 低优先级（plink2
直读 VCF）【决】C13-5；`--IMPUTE/--BEAGLE-GL/--ldhat*/--ldhelmet` 属
相位/似然专门格式，**不在本轮**（06 §4）——其中 IMPUTE/LDhat 的"整位点
零容忍丢弃"语义与 reed 的 chunk 归约正交，可后补。

### §1b.5 工作区其他候选（查过源码后排除）

| 工具/目录 | 排除原因 |
|---|---|
| bcftools `+vcf2smpl`/plugins | bcftools 工作区无此插件；且 reed 已对拍 bcftools 全套 I/O |
| scikit-allel / sgkit / pixy | 消费格式是数组/表格，不是独立"转换目标"；`--freq` 等已有 stats 出口 |
| flashpca/GCTA/admixtools/treemix | 消费 PLINK/文本格式（.ped/.tped/EIGENSTRAT），是 `--plink` 的下游；treemix 的 allele freq 输入 = `--freq` 已有出口 + gzip 脚本 |
| SNPGenie | 消费标准 VCF，不是转换目标 |
| angsd | beagle 似然格式消费方；`--BEAGLE-GL` 同族，归 06 §4 |
| noodles/rust-htslib | 库非工具；reed 已用其 BCF 路径 |

### §1b.6 TASSEL hapmap → VCF（导出侧真值，C0 实测）【据】

> 用户追问"hapmap 应该可以相互转换"——写侧（VCF→hmp）已落地（§1b.2），
> 读侧（hmp→VCF）是真值门。实测命令：
> `run_pipeline.pl -importGuess f.hmp.txt -export vcf_back -exportType VCF`
> （夹具 `c012_clean.vcf`，17 位点 × 8 样本）。

| 位置 | 行为 |
|---|---|
| ExportUtils.java:226 | `##fileformat=VCFv4.0`（不是 4.2/4.3） |
| :228-231 | **无参考**时写 `##Tassel=<ID=GenotypeTable,Version=5,Description="Reference allele is not known. The major allele was used as reference allele">` |
| :232-244 | 恒定声明 `FORMAT` 的 GT/AD/DP/GQ/PL 与 `INFO` 的 NS/DP/AF（**即使输出只有 GT**） |
| :253-257 | `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT` + 样本名（`.trim()`） |
| :260-267 | **REF = 主等位**（hapmap 不存 REF；`p.getAllele(WHICH_ALLELE.Reference)` 在无参考时回落 major），`knownVariants` 存在（来自 hapmap 的 `alleles` 列）时改用 `knownVariants[0]` |
| :488-496 | `REF`/`ALT` 来自 `knownVariants`（`/` 切分）；无 knownVariants 时按频率排序等位，`-`（gap）不进 ALT 而记 indelIndex |
| :539-548 | **无等位位点**：`REF ALT QUAL FILTER INFO` = `. . . PASS .`，全样本 `./.` |
| :554-557 | `QUAL` 恒 `.`；`FILTER` 恒 `PASS` |
| :486-492 | `ID` = `gt.siteName(site)`，即 hapmap 的 `rs#` 列 |
| :474-486 | `CHROM`/`POS` = hapmap 的 `chrom`/`pos` 列 |
| BuilderFromHapMap.java:112-124 | **单/双字母自动判定**：首行样本区字符数 ÷ taxa 数 ≈ 2.0 → 单字母（IUPAC），≈ 3.0 → 双字母（`AA`/`AG`） |
| :258-266（读侧） | 双字母按 `getDiploidValue(input.charAt(i+1), input.charAt(i))` 入核——**第二字符是 allele1**（即 `AG` 读作 G/A） |
| :242-244（读侧） | `pos` 必须是整数，否则 `IllegalArgumentException` |
| NucleotideAlignmentConstants:140-158 | 读侧字母表：A/C/G/T/+/-/Z/N/X/R/Y/S/W/K/M/0（大写）；小写与 V/H/D/B → `UNDEFINED_DIPLOID_ALLELE` → **读回即报错**（B13-10 的根因） |

**实测要点（C0）**：TASSEL 导出的 VCF 里 `CHROM` 是 **`1`/`2`**（hapmap 的
`chrom` 列原样），`ID` 是 **`S1_1`**（TASSEL 内部 siteName，**不是** reed
写的 `chr1:1`）——因为 TASSEL 的 `-importGuess` 会重写 rs# 为自己的命名。
`REF` 与 `ALT` 是**主/次等位**而非原 VCF 的 REF/ALT（pos3 原 `A/G`
→ TASSEL `G/A`，pos5 原 `A/G,T` → `T/A,G`），GT 相应重编号（原 `0/1` →
`1/0`）。**这是有损往返的根源**：hapmap 不存 REF，任何 hmp→VCF 都必须
重新选 REF。

**reed 规范（02 §3.12）**：不模仿 TASSEL 的内部 siteName 重写与 major-as-REF
猜测——reed **照抄 hapmap 的列**（`rs#`→ID、`chrom`→CHROM、`pos`→POS、
`alleles` 的**第一个**→REF、其余→ALT），GT 按"等位在该位点 `alleles`
列表中的下标"编码（与 TASSEL 的 knownVariants 分支同构，:488-496），
IUPAC 格按 AMBIG 反查展开（`R`→`A/G` 两等位 → `0/1`，顺序取 `alleles`
列的声明序，**不做字母序重排**）。`QUAL`=`.`、`FILTER`=`PASS`、`INFO`=`.`、
`FORMAT`=`GT`（不写 TASSEL 那五个未使用的 FORMAT 声明）。缺失格只有
**`N` 与 `?`** → `./.`（TASSEL 写侧用 `N` 表缺失，C0 实测）；`-`/`+`/`0`
是等位字母，按 `alleles` 列查下标（`0` 展开为 `-` 与 `+` 两个等位字母，
仍是等位下标）；`alleles=NA` 的无数据行 → REF/ALT/QUAL/FILTER/INFO = `. . . PASS .`
且全样本 `./.`（对齐 TASSEL 的 nAlleles==0 分支，:539-548）；双字母格按
TASSEL 读侧序（第二字符是 allele1，:258-266）解析后输出非相位 `a/b`；
展开超过两字母（V/H/D/B）或字母不在 `alleles` 声明集 → malformed 跳过 +
计数（B13-2 口径，参照读侧直接报错）。

## §2 easySFS / stats --sfs 现状（迁移对象）

| 事实 | 位置 |
|---|---|
| 请求类型 `SfsRequest { projections: Vec<u32>, folded: bool }`；报告 `SfsReport { spectra, projections, folded, n_sites, n_skipped_m_lt_n, joint }` | `crates/reed/src/stats/request.rs:435-470` |
| 折叠与 partial：`SfsPartial::new`、`push_sfs`（每行消费，超几何权重经 `ProjectionCache`） | `stats/run.rs:494-560` |
| 写出：`write_projected_sfs` / `write_joint_sfs`（dadi + fsc2 布局） | `stats/writers/sfs.rs:26-123` |
| CLI 解析与校验（`--proj` 缺失报 usage） | `reed-cli/src/translate/stats.rs:715-763` |
| 门面导出：`SfsRequest`/`SfsReport` 已在根 | `crates/reed/src/lib.rs:89-94` |
| 对拍门：`gold/popgen/run_easysfs.sh`（unfolded+folded × 1D+joint，rtol 1e-6） | gold/popgen |

## §3 语义歧义点核实（实现前必须知道的坑）

> 通读源码发现的、README 没写的语义细节。规范（02 §3）已按"对齐参照"
> 逐条落定，这里留证据。

1. **ALT 含 `<NON_REF>` 或 `-` 的处理顺序**（:97, :113-114）：先 `replace
   ("-", "*")` 再 `replace("<NON_REF>", REF)`、REF 也 `-`→`*` 且**大写化**。
   `-` 是 GATK 删除符号，进字母表成 `*`，参与 IUPAC 消歧（`*A`→`a`）。
2. **`is_snp` 的判定是"REF 长度 1 且 ALT 每个等位长度 1"**（:95）——
   ALT 数不限（多等位 SNP 保留，GT 里的 2/3 号等位照样进字母表）；单核苷
   酸 ALT 与 `*`/`<NON_REF>` 组合也过。**reed 的 `--snps-only`（just-acgt）
   语义不同**——convert 不复用它，见 §5 B13 分歧表。
3. **AMBIG 未命中是崩溃，不是跳过**（实测修正，C0 复现）：`nt_dict` 收
   REF/ALT 原文（:112），GT 用到的等位组合不在 AMBIG 表时，:129
   `AMBIG[geno_nuc]` 的 KeyError **落在 try 块外**（try 只包
   `nt_dict[j]` 的 KeyError）→ vcf2phylip **直接崩溃**。实测：`REF=R` +
   GT `0/0` 即 `python3 vcf2phylip.py` 退出 traceback（`KeyError: 'R'`）。
   等位**索引**越界（GT 写 `2/2` 而只有 1 个 ALT）才走 `"malformed"` 跳过。
   **reed 规范（B13-9）**：两种情况都跳过 + `malformed` 计数 + 一次性警告
   （不崩溃、不刷屏），对拍断言排除这两种记录。
4. **min-samples-locus 数的是样本不是染色体**（:100-108）：`0/.` 半调用
   不以 `.` 开头 → **算有基因型**。与 reed `--geno`（缺失**等位**/2N）分
   母不同——02 §3.7 的帮助文本已写明。
5. **`-p` 是"关闭默认"而 reed 的 `--phylip` 是"开启"**：B13-1，02 §1.1。
6. **bin.nexus 的 `?` 不区分"GT 缺失"与"非双等位样本 GT"**（GEN_BIN 兜底）。
7. **outgroup 静默忽略不存在的名字**（:396 `if outgroup in sample_names`
   无 else）——reed 改报错，B13-4。
8. **used_sites 只记核苷酸流**（:327-329 在 `fasta or nexus or not
   phylipdisable` 分支内）：bin-only 跑法不写 used_sites（02 §3.6 对齐）。

## §4 reed 侧可复用面盘点

| 能力 | 位置 | 复用方式 |
|---|---|---|
| chunk 流水线（reader/workers/reducer、首错即停、预算、指标） | `reed-engine::pipeline` | 直接作为 `ConvertJob` 的执行壳（02 §7） |
| `LineRef` 逐行解析 + GT 扫描（`GtCode`/`GenotypeScanner`） | `reed-core`、`stats/run.rs:567` | **注意**：phylo 列转换需要**等位索引 → 字母**的多等位视图，而 `GtCode` 是 CollapseAlt 后的双等位位平面——**不能直接用 GtCode**，需保留原始等位索引（02 §7.2 的列折叠在行文本上自取 GT 字段，或扩展 scanner 的原始视图）。这是 D13-3 的正面表述 |
| SFS 折叠（`push_sfs` + `ProjectionCache`） | `stats/run.rs:494-560` | 组合复用（02 §7.5），零拷贝不漂移 |
| SFS 写出 | `stats/writers/sfs.rs` | 直接调用 |
| 群体装载（`--pop`、`Populations`） | `stats/request.rs:256-…` | 直接复用 |
| 样本子集参数组 | `reed-cli/cli/filter_flags::SampleArgs` | CLI flatten 复用 |
| 输入打开（VCF/gz/BCF/stdin 统一） | `crates/reed/src/input.rs` | `run_convert_opened` 对称 |
| 确定性/退出码/verbosity 通道 | 全仓既有 | 直接沿用 |

## §5 差异分档表（D13 / B13 / C13）

### D13（不做就是错的——对拍会挂）

| 项 | 内容 | 规范 |
|---|---|---|
| D13-1 | **三格式同模型不变式**：fasta/phylip/nexus 必须同一份列模型渲染，禁止三份独立转换实现（doc_v12 §0.5 元模式 1 的预防性约束） | 02 §2.3、§9.2 |
| D13-2 | **禁止临时文件转置**：vcf2phylip 的 `.tmp` 竖排 + 逐行取列是 O(N·L) 重复 IO（每样本重读全临时文件，:430-446）；reed 在 partial 内存里列分块、收尾一次转置 | 02 §7.2、§7.6 |
| D13-3 | **多等位 GT 视图**：`GtCode`（CollapseAlt）不可用于核苷酸列；列折叠必须取**原始等位索引**（0/1/2/3…）映射字母 | 02 §3.3、§4 |
| D13-4 | **`.minN` 前缀后缀不移植**：reed 用 `-o PREFIX` 显式前缀 + `.log` 计数，输出名不带隐藏参数痕迹 | 02 §3.6 |
| D13-5 | **SFS 折叠单一实现**：convert 侧不得出现 `push_sfs` 的第二份拷贝（stats/convert 漂移预防） | 02 §7.5、§4.2 |
| D13-6 | **`--sfs` 两入口等价**：convert ≡ stats 逐字节（含校验错误文本）；翻译层共享单一 `build_sfs_request` | 02 §4 |
| D13-7 | **012/hapmap/msa 共享位点门与行解析**：`--012`/`--hapmap` 复用同一 chunk 行解析与 `--min-samples-locus` 门；禁止第四份独立的行级 GT 解析（D13-3 同族：多等位原始索引视图一处实现，phylo/012/hapmap/msa 四目标共用） | 02 §3.8/§3.9/§3.10 |
| D13-8 | **hapmap→VCF 是 record→record 映射，不是矩阵折叠**：两格式都站点为主序，逐行翻译（IUPAC 反查展开成 GT），复用 `-O` 的 `OrderedOutput`/`BlockWriter` 写出（不新造压缩路径）；头**无 `##contig`**（TASSEL 同款）因此 CHROM 顺序不参与任何头状态，原样转写不强制 sorted（下游索引工具自会拒绝乱序，与消费 TASSEL 导出物同款约束） | 02 §3.12、§7.7 |

### B13（有意分歧——登记进 cli.md 差异表，不写对拍断言）

| 项 | reed | 参照 | 理由 | 规范 |
|---|---|---|---|---|
| B13-1 | `--phylip` opt-in（默认不写） | PHYLIP 默认写（`-p` 关闭） | reed"要什么写什么"；对拍按指定等价集 | 02 §1.1 |
| B13-2 | 畸形行跳过+计数+一次性警告 | 打印整行跳过 | 3 M 行级别的 stdout 刷屏不可用；`--verbosity 0` 可静默 | 02 §3.1 |
| B13-3 | `--resolve-iupac` 固定 seed=0，确定性可复现 | 每次运行不同（全局 random） | reed 确定性是硬约束（architecture.md §5.10）；文档写明 | 02 §3.3 |
| B13-4 | outgroup 名不在样本集 → usage 错（退出码 2） | 静默不重排 | cli.md §0.7 校验哲学 | 02 §3.5 |
| B13-5 | hapmap 的 `alleles` 列取 **REF+ALT 声明集** | TASSEL 写"观察等位集"（按全体样本频率排序，ExportUtils.java:157-174） | 观察集需全体样本物化后频率排序，与 chunk 加法归约冲突；声明集确定且信息等价（下游 TASSEL 读侧只按 `/` 切分） | 02 §3.9 |
| B13-6 | msa 的 indel 位点 gap 摆放用**确定性填充**（ALT 字面长度差，`N`/`-` 填到 maxlen） | vcf2msa 调 **ClustalO 外部进程**对齐（m:198-207），失败才走同款填充 | 不引外部进程依赖；填充路径与 vcf2msa 的 MUSCLE 失败回退一致；文档写明"非多序列对齐" | 02 §3.10 |
| B13-7 | msa 样本名**不折叠** `.` 后缀（`--collapse-suffix` 显式开关） | vcf2msa 恒按第一个 `.` 折叠（m:56-60） | 折叠是作者个人管线约定（README 明言 "might not be the desired behavior"）；默认保守，显式开启 | 02 §3.10 |
| B13-8 | `--012` 的半调用 `0/.` 与 haploid 输出 `-1`（缺失） | vcftools 输出 `1`（"其余一切→1"分支，vft-conv:445；haploid REF 也落 1，C0 实测） | vcftools 的三态分支把半调用、haploid 与杂合混同；reed 按既有 `GtCode` 语义（半调用=半信息）置缺失更诚实；对拍断言排除这些样本（fixture 分行） | 02 §3.8 |
| B13-9 | AMBIG 表未命中的记录**跳过 + 计数**（REF/ALT 含表外字母，或等位索引越界） | 表外字母 **KeyError 崩溃**（try 块只护索引查找，vcf2phylip:129）；索引越界走 malformed 跳过 | 3M 记录级输入上崩溃丢整轮；reed 统一为"跳过+计数+一次性警告"，与畸形行同口径（B13-2） | 02 §3.1/§3.3 |
| B13-10 | hapmap 格收敛到 TASSEL 读侧字母表（A/C/G/T/N/R/Y/S/W/K/M/-，其余 → `N`） | TASSEL 写侧输出 vcf2phylip 风格的全 IUPAC（含小写软掩蔽码） | TASSEL **自己读不回**小写与三碱基码（实测 "illegal value: a"）；写一个自己生态读不了的文件没有意义 | 02 §3.9 |
| B13-11 | hmp→VCF **照抄 hapmap 列**（`rs#`→ID、`alleles[0]`→REF、其余→ALT，GT 按 `alleles` 声明序编码） | TASSEL 重写 `ID` 为内部 siteName（`S1_1`）、用 **major allele** 当 REF 并重编号 GT | hapmap 不存 REF，任何 hmp→VCF 都有损；reed 选**可预测**的"声明序"而非 TASSEL 的频率猜测——同输入同输出，且 `--hapmap`→`--vcf` 往返在双等位位点上**保真**（REF/ALT/GT 复原），TASSEL 路径不保证 | 02 §3.12 |
| B13-12 | hmp→VCF 写 `##fileformat=VCFv4.3` + 仅 `##FORMAT=<ID=GT>`（reed 输出头风格） | TASSEL 写 `VCFv4.0` + GT/AD/DP/GQ/PL 五个 FORMAT 声明（实际只填 GT） | 未使用的 FORMAT 声明是噪声；reed 的输出头对齐自家 filter（`##reed_*` provenance 可选）；VCF 版本取 reed 现行 4.3 | 02 §3.12 |

### C13（待拍板——默认已定，可推翻）

| 项 | 默认 | 备选 | 规范 |
|---|---|---|---|
| C13-1 | convert 不带 GT 掩码族（`--min-dp/--max-dp/--min-gq`） | 带上（转换前掩码，对齐 plink2 风格） | 02 §2.5 |
| C13-2 | `--min-samples-locus` 独立参数（整数样本门） | 复用 `--geno`（分数等位门） | 02 §3.7 |
| C13-3 | SFS 目标 `--pop` 可缺省（全体一组，沿 stats） | 必选 | 02 §1.2 |
| C13-4 | msa 不折叠样本名 `.` 后缀 | 对齐 vcf2msa 默认折叠 | 02 §3.10 |
| C13-5 | `--plink`/`--plink-tped` 不进首发轮 | 进（vcftools 全套平移） | 02 §3.11、06 §4 |
| C13-6 | `--impute/--beagle-gl/--beagle-pl/--ldhat*/--ldhelmet` 不做（本轮） | 后续轮按需求单独立项 | 06 §4 |

## §6 守卫清单（落地后进 `ci/check-doc-drift.sh`）

1. `Command` 枚举恰四变体（Filter/Stats/Index/Convert）；`reed --help`
   列出 convert。
2. `translate/` 下 SFS 请求构建只有一个共享函数（`build_sfs_request`），
   `translate/stats.rs` 与 `translate/convert.rs` 都调它（D13-6）。
3. `push_sfs` 全仓唯一实现点（D13-5）。
4. `reed-core` 含 `iupac` 模块；AMBIG 表项数 = 62、GEN_BIN = 10（防抄漏）。
5. convert 无 `--mind`（出现即解析错误，对齐 stats 侧做法）。
6. `--phylip` 默认 false（帮助文本不暗示默认写 PHYLIP）。
7. `resolve_iupac` 输出确定性测试存在且 threads ∈ {1,4} 双跑。
8. `run_convert` 请求按引用、根导出、README 示例在编译门内（v11 纪律
   自动延伸）。
9. `gold/popgen/run_vcf2phylip.sh` 存在且 CI 侧可寻址（gold 目录登记）。
10. convert 的渲染层无临时文件 API 调用（`tempfile`/`NamedTempFile` 不进
    convert 模块；D13-2）。
11. 多等位原始索引视图（D13-3/D13-7）全仓单一实现：phylo 列、012 列、
    hapmap 格、msa 解析四消费者共用的行级 GT→等位索引函数只有一个
    定义点（`grep` 计数守卫）。
12. `--012` 半调用置 `-1` 的单测存在（B13-8 的契约面）；fixture 里半调用
    样本与杂合样本分行（对拍可分别断言）。
13. `--hapmap` 写出的 `alleles` 列来自 REF+ALT 声明（有"含未观察 ALT"
    的测试位点；B13-5 的契约面）。
14. AMBIG 未命中/索引越界走 malformed 跳过路径（`e2e_convert` 的 fixture
    含 `REF=R` 与越界索引两条；B13-9 的契约面）。
15. hapmap 格的字母表收敛有测试（B13-10 的契约面，`hapmap_layout_and_alphabet`）。
16. `--vcf` 目标存在且与 VCF 输入目标互斥（同给 → usage 错，退出码 2）；
    非 hapmap 输入 → 运行期 header 错（退出码 1，与 `filter`/`stats` 吃到
    非 VCF 输入同口径）；`ConvertTarget::Vcf` 在门面根导出（D13-8）。
17. hmp→VCF 的 IUPAC 反查表（`R`→`A/G` 等）与 AMBIG 互为逆运算有测试
    （`iupac` 模块的 `ambig_expand` 单测 + `e2e_convert` 的双等位往返保真）。
18. hmp→VCF 的头形态锁定：**无 `##contig`**（TASSEL 也不写）且
    `##reed_version` 存在、`--no-version` 可关（B13-12 的契约面；
    `hapmap.rs` 里不得出现字面 `"##contig`）。

## §7 对拍设计（真值门）

`gold/popgen/run_vcf2phylip.sh`（新）：

- **fixture**：`gold/fixtures/phylo.vcf`（新；8 样本 × ~20 位点，覆盖：
  单/多等位 SNP、MNP、indel、`*` ALT、`<NON_REF>`、`./.` 全缺失列、
  haploid GT、phased GT、REF 非纯 ACGT 一条、半调用 `0/.`、列数错误行）。
  期望输出由 vcf2phylip 2.9 现跑生成并**人工核对后**入库（gold 惯例：
  make_fixtures.py 一并生成）。
- **全量**：`clean.vcf.gz`（348 样本）phylip + fasta + nexus + bin +
  resolve + used_sites 六跑，`diff -r` 全等（B13 分歧点用**等价集**参数
  匹配：reed 显式给 `--phylip`；vcf2phylip 不给 `-p`——两边格式集相同）。
- **SFS 迁移门**：`run_easysfs.sh` 增跑一遍 `reed convert --sfs` 入口
  （同 fixture 同容差），与 stats 入口 diff -r 空。
- 门文本：`Phylo parity: PASS` / `SFS parity via convert: PASS`。

`gold/popgen/run_convert_parity.sh`（新，第二轮增补的三个格式）：

- **fixture**：`gold/fixtures/convert012.vcf`（半调用/杂合/haploid/多等位/
  ALT='.' 分行覆盖 B13-8 场景）与 `gold/fixtures/convert_hmp.vcf`（含"未
  观察 ALT"位点覆盖 B13-5）；msa 另需 `convert_msa.fa` 参考 + 小区域 VCF。
- **真值**：`--012` 对 vcftools 工作区版（`vcftools --012`，注意 fixture
  排除半调用样本或断言按 reed 契约）；`--hapmap` 对 TASSEL `ExportUtils`
  的**读侧**验证（用 TASSEL BuilderFromHapMap 读回 reed 的 .hmp.txt 断言
  样本数/位点数/等位——写侧无独立 CLI 可对拍，读回验证是 TASSEL 生态
  的实际消费方式）；msa 对 vcf2msa.py 无 indel 区域 diff 全等、含 indel
  区域结构断言（B13-6）。
- **SFS 三入口**：convert/stats 双入口等价门沿用。

`gold/popgen/run_convert_parity.sh` 门 6（第三轮增补，hapmap↔VCF 反向）：

- **fixture**：脚本内联 4 位点双等位无缺失 VCF（往返无损域）。
- **真值**：VCF → `--hapmap` → `--vcf`，`grep -v '^#'` 后的记录行与输入
  **逐行相同**（B13-11 的保真面）；缺省 `-O` 必须是纯文本、`-Oz` 必须
  gzip 帧（`1f 8b`）；TASSEL 5 `-importGuess` 读 reed 的 VCF 再导出，
  站点数保持（B13-12 的读侧契约）。
- 门文本：`hapmap ↔ VCF round trip and TASSEL read-back of reed's VCF`。

## §8 性能（实测基线，2026-10-16）

**全量对拍与性能同一负载**：`clean.vcf.gz`（8 GB BGZF、**348 样本、
3 315 604 记录**、`Chr01` 起），28 核，phylip+fasta+nexus+bin+used_sites
五输出一次跑（vcf2phylip 用 `-f -n -b -w`）。

| 负载 | 参照（单线程） | reed convert（28 核） | 实测加速 |
|---|---|---|---|
| 五格式一次跑，全量 | vcf2phylip 2.9：**27 min 27 s**（wall；与之相符 README 的 20 GB/27 min 自述——瓶颈是逐样本重读 1.16 GB 临时文件的转置，vcf2phylip:430-446） | **10.4 s**（wall；user 3 m 33 s） | **≈ 158×**（门 ≥8×；即使按 reed 单线程 CPU 时间 ≈3.5 min 比也是 ≈8×） |
| 输出体量 | `.phy` 1 153 835 076 B、`.bin.nexus` 1 153 835 173 B、`used_sites` 62 242 937 B | **五文件与 vcf2phylip 逐字节相同**（`cmp` 全等） | 门 §9.1 全量档完成 |
| 峰值内存 | — | 8.2 GiB（= mmap 输入页 ~8 GB，可回收 + 列矩阵 ~2.3 GB + worker 往返缓冲 ~0.2 GB；`--memory` 只约束拷贝式源） | 无硬门 |
| 线程确定性 | — | threads 28 vs 4：全量 `.phy`/`.bin.nexus` `cmp` 相同 | 门 §9.5 |

`--012`／SFS 全量的参照基线未测（低优先）；fixture 档的逐字节门已在
`run_convert_parity.sh` 内。数字随本表进入 STATUS。

> 注：`R0 基线步`字样已随轮次编号统一改为 `C0`（v13 用 C0–C3）。
