# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed 对 VCF/BCF 文件计算群体遗传学统计量并执行硬过滤，输出与原文件逐字节一致的
格式。纯 Rust 实现，默认用满所有核心，不需要系统库：BGZF 与 BCF 编解码直接编译
进二进制，工具链能跑的地方它就能跑。它同时也是 Rust 库（`reed` crate），其他
工具可以直接调用同一套统计，不必起子进程。

## 安装

```bash
cargo install --path crates/reed-cli     # 需要 Rust 1.91+
reed --version
```

## 用法

```bash
# 全基因组与窗口 π、dxy、Fst、θ、Tajima's D
reed popg allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out

# 只含变异位点的输入
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --fst --window 100k --step 50k -o out

# 位点级与样本级 QC 报告 + bcftools-stats 汇总，一次扫描
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# 同一遍出直方图：内建指标与任意数值 INFO/FORMAT 标签
reed stats clean.vcf.gz --hist qual,dp-site,gq,INFO/DP -o qc

# 标识符编辑：样本名、染色体名、变异 ID、header 行
reed rename in.vcf.gz -s "s0 sampleA" --rename-chrs chrs.txt \
    --set-missing-var-ids '%CHROM:%POS' -Oz -o renamed.vcf.gz

# 位点与基因型过滤，输出带索引
reed filter in.vcf.gz -o out.vcf.gz -Oz6 --geno 0.1 --maf 0.05 \
    --min-dp 5 --min-gq 20 -W

# LD 剪枝
reed filter in.vcf.gz --indep-pairwise 50kb 5 0.2 -Oz6 -o pruned.vcf.gz

# CSI 索引
reed index in.vcf.gz

# 系统发育矩阵、0/1/2 矩阵、hapmap 与参考 MSA
reed convert in.vcf.gz --phylip --fasta --nexus --nexus-binary --012 --hapmap -o matrices
reed convert in.vcf.gz --msa --reference ref.fa --region chr1:1-100000 -o aln

# 反向：hapmap 转回 VCF（-Oz 写 BGZF；不给 -o 时流到 stdout）
reed convert panels.hmp.txt --vcf -Oz -o panels.vcf.gz

# 投影 SFS，与 `reed popg --sfs` 逐字节相同
reed convert clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out
```

`reed <command> --help` 列出全部选项。

## reed popg

窗口群体遗传统计学。统计开关任意组合，每个开关产出自己的文件。

| 开关 | 输出 |
|---|---|
| `--pi` | 每群体 π，`PREFIX_pi.txt` |
| `--dxy` | 每群体对 dxy，`PREFIX_dxy.txt` |
| `--fst [wc84\|hudson]` | 全部群体对的 Fst，`PREFIX_fst.txt` |
| `--fst-global` | 全群体 WC84 Fst，`PREFIX_fst_global.txt` |
| `--theta` | 每群体 Watterson θ，`PREFIX_watterson_theta.txt` |
| `--tajima-d` | 每群体 Tajima's D，`PREFIX_tajima_d.txt` |
| `--pbs A,B,C` | 指定三重 PBS，`PREFIX_pbs.txt`；`--pbs-normed` 产 PBSn1 |
| `--sfs` | 投影 SFS，落在 `PREFIX_sfs/dadi/` 与 `PREFIX_sfs/fastsimcoal2/`；需要 `--proj N` |
| `--f2 A,B`、`--f3 C;A,B`、`--f4 A,B;C,D`、`--dstat A,B;C,D`、`--f4-ratio A,B;C,D;X` | Patterson f 统计量 + 区块刀切，`PREFIX.fstats` |
| `--ld-decay [MAX_BP]` | 1 kb 分箱的配对 r²，`PREFIX.ld.decay`；`--ld-pairs` 加逐对表 `PREFIX.ld` |

窗口选项：

| 开关 | 含义 |
|---|---|
| `--window <SIZE>` | 窗口大小（bp）；接受 `100k`、`1M`（默认 100000） |
| `--step <STEP>` | 窗口步长（默认等于窗口大小） |
| `--levels <LIST>` | 输出的层级：`window`、`contig`、`genome`（默认 `window`） |
| `--skip-empty` | 丢弃无位点落入的窗口 |
| `--pop <FILE>` | `样本<TAB>群体` 行、无表头；窗口模式必需 |
| `--variant-only` | 按全变异输入处理：文件不含的位置计为纯合且全部样本有调用 |
| `--bypass-invariant-check` | 只在文件含有的位点上计算（π/dxy 分母有偏） |
| `--contig-lengths <FILE>` | `.fai` 式 `名称<TAB>长度`，用于缺 `##contig` 行的头 |

`popgen` 是可见别名。这些窗口开关在 0.3.x 属于 `reed stats`；现在在
`reed stats` 上使用会得到指出改写命令的解析期报错（不是静默别名）。

## reed stats

位点级与样本级 QC 报告。每个开关产出自己的文件。

| 开关 | 输出 |
|---|---|
| `--freq`、`--hardy`、`--het`、`--missing`、`--depth` | `PREFIX.afreq`、`PREFIX.hardy`、`PREFIX.het`、`PREFIX.smiss` + `PREFIX.vmiss`、`PREFIX.idepth` + `PREFIX.ldepth.mean` |
| `--summary` | bcftools-stats 布局的 QC 汇总，`PREFIX.stats` |
| `--hist <LIST>` | 每个指标一份 `PREFIX.<m>.hist` TSV |

`--hist` 用逗号分隔内建指标（`qual`、`dp`、`dp-site`、`gq`、`idl`、`af`）
与任意数值标签（`INFO/<TAG>` 或 `FMT/<TAG>`）。每份文件含 `bin_lower`、
`bin_upper`、`count`、`percent` 列；整数指标另附精确 `#mean`、`#median`、
`#p95` 行。标签未在 header 声明、或不是 Integer/Float 类型，都是用法错误。

`--pop FILE` 在这里可选（分层 QC）；`popg` 的窗口模式必需。与
`--freq/--missing/--singletons` 之一同给时分层：`PREFIX.afreq.pop` 与
`PREFIX.vmiss.pop` 每群体一行（按首现序；数值等于把该群体样本单独取出后
的 plink2 跑法），`PREFIX.singletons` 加 POP 列（携带者不在群体文件时写
`.`）。群体文件之外的输出不变；`--pop` 不与这三个 flag 之一同给则是用法
错误。多等位 ALT 频率保持折叠（与 `.afreq` 一致）。

## reed rename

标识符编辑：样本名、染色体名、变异 ID 与 header 行。全程不解码基因型，
记录逐字节直通。

| 开关 | 含义 |
|---|---|
| `-s "old new,…"`、`-N FILE` | 样本改名（文件里是 `old<TAB>new`） |
| `--replace-sample-name PATTERN` | 字面改名（`pre*` 去前缀、`*suf` 去后缀，否则删除该文本；正则元字符报错） |
| `--rename-chrs FILE` | `old<TAB>new` 染色体名；`##contig` 行与 CHROM 列同改 |
| `--chr-style plain\|ucsc\|mt-is-m` | 染色体命名规范化（无 `chr`、带 `chr`、`chrM`） |
| `--set-id PATTERN` | 用 `%CHROM %POS %REF %ALT` 模板重写全部 ID |
| `--set-missing-var-ids PATTERN` | 只补 `.` 的 ID |
| `--update-name FILE`、`--idspace-to SEP`、`--rm-dup-ids`、`--clear-ids` | ID 映射、换分隔符、重复加后缀、全部置 `.` |
| `-h FILE`、`-f FILE`、`--add-meta LINE`、`--remove-meta KEY` | 换 header、按 `.fai` 补 `##contig` 长度、增删 header 行 |

输出：`-o FILE`（默认 stdout）、`-O v|z[level]`。`reheader` 是可见别名。
改成已存在的名字、名字不存在、或没有任何编辑请求，都是用法错误。

## reed filter

选择，沿用 bcftools：`-r/-R`、`-t/-T`、`-s/-S`、`--keep`、`--remove`
（`--force-samples` 把未知样本名降级为警告）、`-v/-V`、`--snps-only`、
`-m/-M`、`-f`、`-i/-e`。

位点过滤，在掩码之后应用：`-q/-Q/--min-af/--max-af`（可带
`:nref|alt1|minor|major|nonmajor` 后缀）、`--geno`、`--maf`、`--max-maf`、
`--mac`、`--max-mac`、`--hwe <P [midp]>`。

基因型掩码：`--min-dp`、`--max-dp`、`--min-gq` 从 FMT/DP 与 FMT/GQ 掩码；
`-i/-e` 的 FORMAT 表达式配 `--set-gts .` 精确掩掉命中样本，与
`bcftools filter -S .` 一致。被掩基因型对后续位点过滤按缺失计。

样本过滤：`--mind <FLOAT>` 丢弃缺失率超阈值的样本。
`--indep-pairwise <W[kb] [S] R2>` 用 plink2 算法按 LD 剪枝并输出保留位点；
窗口 `W` 可带 `kb` 后缀，步长 `S` 默认 1，`R2` 是 r² 上限。

INFO 处理：样本集变化时 reed 更新 `INFO/AC` 与 `INFO/AN`（无 GT 列的记录
保留输入值）；`-I/--no-update` 关闭。`-x/--remove-info LIST` 删除 tag——默认
删 gold 管线的陈旧清单（`AC,AN,AF,MLEAC,MLEAF,ExcessHet,InbreedingCoeff,DP`），
`-x ''` 保留全部。`--fill-tags LIST` 重算 AC、AN、NS、MAF 或 HWE（大小写不
敏感；不支持的标签报错，与 bcftools 一致）。

输出：`-O v|z|b|u[0-9]` 选择普通 VCF、BGZF VCF、BCF 或未压缩 BCF，可带压缩
级别。`-o` 指定文件（默认 stdout），`-W[=csi|tbi]` 写索引，`--report` 写 JSON
运行报告，`--timing` 记录各阶段耗时，`--no-version`、`--no-command` 抑制
`##reed_*` 来源信息行。

## reed index

`-c/--csi`（默认）或 `-t/--tbi`，`-o` 改写到别处，`-f` 覆盖已有索引。输入
必须是文件：索引需要可寻址的源。

## reed convert

一次并行读入、多目标输出（各目标共享同一 worker 池，输出与 `--threads`
无关）。语义严格对齐参照工具（vcf2phylip 2.9（Ortiz 2019）、vcftools
`--012`、TASSEL hapmap、vcf2msa），`gold/popgen/run_convert_parity.sh`
逐字节校验。

| 参数 | 输出 |
|---|---|
| `--phylip` / `--fasta` / `--nexus` | `PREFIX.phy` / `.fasta` / `.nexus` 核苷酸矩阵（relaxed PHYLIP 布局、IUPAC 简并码；`--phylip` 是 opt-in，与 vcf2phylip 相反） |
| `--nexus-binary` | `PREFIX.bin.nexus`（SNAPP；双等位二倍体位点，其余 `?`） |
| `--012` | `PREFIX.012`（每样本 0/1/2/-1）、`.012.pos`、`.012.indv` |
| `--hapmap` | `PREFIX.hmp.txt`（TASSEL 布局） |
| `--vcf [-O v\|z[level]]` | **反向转换**：hapmap（`.hmp.txt`）输入转回 VCF——`ID` 取 `rs#`、REF/ALT 按声明顺序、IUPAC 格展开为 `GT`（`R` → `0/1`）；`-O z` 写 BGZF，不给 `-o` 时流到 stdout |
| `--msa --reference ref.fa [--region CHR:B-E]` | `PREFIX_<contig>.fasta`：逐参考碱基的全位点比对 |
| `--sfs --proj N[,N…] [--pop FILE]` | `PREFIX_sfs/dadi`、`PREFIX_sfs/fastsimcoal2`，与 `reed popg --sfs` 相同 |
| `--outgroup NAME`、`--min-samples-locus N`、`--resolve-iupac`、`--write-used-sites` | 矩阵调节（见 `--help`） |

`--hapmap` 与 `--vcf` 是同一对格式的两个方向：在无损域（双等位、二倍体、
无缺失的 SNP）上 VCF → hapmap → VCF 往返逐记录还原，且 TASSEL 5 能读回
reed 写的 VCF。

与参照工具的有意分歧（登记在 `docs/spec/cli.md` §8b）：畸形记录跳过并计数
（不崩溃、不刷屏）；`--resolve-iupac` 确定性（固定 seed）；`--012` 的半调用
与 haploid 写 `-1`（vcftools 写 1）；hapmap 的 `alleles` 列写声明集；
MSA 的 indel 用确定性填充（不调外部对齐器）；hapmap→VCF 的头不写
`##contig`（TASSEL 也不写），`##reed_*` 溯源行可关。

8 GB VCF（348 样本、331 万记录）上五个系统发育输出与 vcf2phylip 2.9
**逐字节相同**，耗时 **10 秒** vs 其 27 分钟（28 核）——逐样本临时文件转置
换成内存列块。

## 输入要求

- VCF、VCF.gz 或 BCF，按文件头识别。`-` 从 stdin 读 VCF 文本；`reed index`
  必须给文件。
- 头必须声明 `##contig=<ID=…,length=…>`。窗口统计按声明长度裁剪每个 contig
  的末窗，索引需要 contig 顺序。头缺这些行时用 `--contig-lengths FILE`。
- 记录必须按 contig 与位置排序。索引与窗口统计依赖这一点。
- `-r/-R` 需要输入旁有 `.csi` 或 `.tbi`（用 `reed index` 建）。普通扫描不需要。
- 基因型必须是二倍体。reed 默认拒绝多等位记录（`-M 2`），`0/.` 这类半调用
  按等位基因级处理，`--geno` 与 `--mind` 的计数与 plink2 一致。

## 输出

`filter` 与 `index` 写出 `-o` 指定的文件。`stats` 把 `-o` 当前缀，每个统计
一份文件（见上表）。不给 `-o` 时，窗口统计输出 TSV 到 stdout，报告文件不写。

来源信息行：`filter` 的输出头部带 `##reed_version`、`##reed_command`、
`##reed_filterCommand`、`##reed_filter`；`stats` 的写进 `PREFIX.log`。

窗口表统一形态：每个窗口每群体（或群体对）一行，列为 `chromosome`、
`window_pos_1`、`window_pos_2`、估计值与其背后的计数（π/dxy 的 `no_sites`、
`count_diffs`、`count_comparisons`、`count_missing`；Fst、PBS、f 统计量的
`no_snps`）。分母为零的估计值打印 `NA`。

## 作为库使用

`reed` crate 是公开 API：typed 请求进，结构化报告出，writer 把报告写成与
CLI 相同的文件。

```toml
[dependencies]
reed = { path = "crates/reed" }
```

```rust
use reed::{PopgenRequest, Populations, StatKind};

let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::File("pops.txt".into()),
    ..PopgenRequest::default()
};
let report = reed::run_popgen(Some("in.vcf.gz".as_ref()), &request, &reed::budget::ThreadBudget::default())?;
let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false);
```

转换请求同形；目标类型化，报告列出写出的全部文件：

```rust
use reed::{ConvertRequest, ConvertTarget, PhyloFormats, PhyloOptions};

let request = ConvertRequest {
    targets: vec![ConvertTarget::Phylo {
        formats: PhyloFormats { phylip: true, fasta: true, nexus: false },
        options: PhyloOptions::default(),
    }],
    prefix: "matrices".into(),
    ..ConvertRequest::default()
};
let report = reed::run_convert(Some("in.vcf.gz".as_ref()), &request, &reed::budget::ThreadBudget::default())?;
```

手上已经有 `样本 -> 群体` 映射的调用方可以直接传入，不必落一个群体文件：

```rust
let request = PopgenRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::Assignment(vec![
        ("sample1".into(), "popA".into()),
        ("sample2".into(), "popB".into()),
    ]),
    ..PopgenRequest::default()
};
```

门面 re-export 所有需要的东西：`run_filter`、`run_popgen`、`run_stats`、
`run_rename`、`run_index`、`run_convert`，请求类型（`PopgenRequest`、
`FilterSpec`、`StatsReportRequest`、`RenameSpec`、`ConvertRequest`），
报告（`PopgenReport`、`FilterReport`、`StatsOutcome`、`RenameReport`、
`IndexOutcome`、`ConvertReport`），
writer（`reed::windows`、`reed::reports`、`reed::histogram`、`reed::sfs`、
`reed::fstats`、`reed::ld`、`reed::summary`），以及单一 `reed::Error`
（`kind()` 区分用法错误与运行期失败）。`PopgenRequest::to_args(&group_names)`
与 `FilterSpec::to_args()` 把请求渲染回 CLI 拼写，程序化构造请求的工具也能
记录等价的命令行；CLI 拼写不出来的请求（per-bp kind、`Regions` grid）由
`PopgenRequest::to_args` 返回用法错误。只有这个门面是稳定性承诺；
`crates/` 下的下层 crate 是内部实现，`reed::Error` 通过封口变体把它们挡在
签名之外——用 `kind()` 分类，用 `Display`/`source()` 读消息。

## 已知限制

- 多等位记录被拒绝。`-M` 提到 2 以上时 reed 不会重算逐等位统计，默认保持拒绝。
- LD 剪枝在 `filter` 里，需要第二遍扫描。`reed stats` 与 `reed popg` 都拒绝
  `--indep-pairwise` 并指向 `reed filter`。
- f 统计量的标准误用按碱基对的区块刀切，scikit-allel 按变异数；估计值一致，
  标准误不可比。
- Tajima's D 按窗口平均等位数处的方差做分母（pixy 约定），与 vcftools 的固定
  2N 不同。

## 引用

reed 实现的是已发表的估计量，请引用它使用的方法：

- π：Nei & Li 1979；dxy：Nei 1987；θ：Watterson 1975；D：Tajima 1989。
- Hudson Fst：Hudson 1992 与 Bhatia et al. 2013；WC84 Fst：Weir & Cockerham 1984。
- f 统计量：Patterson et al. 2012；PBS：Yi et al. 2010。
- π、dxy、θ、D 的缺失数据处理：Korunes & Samuk 2021。

## 许可

MIT OR Apache-2.0。见 [LICENSE-MIT](LICENSE-MIT) 与 [LICENSE-APACHE](LICENSE-APACHE)。
