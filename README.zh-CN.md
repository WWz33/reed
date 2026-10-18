# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed 是用 Rust 写的群体遗传学 API 库。它在 VCF/BCF 上计算 π、dxy、Fst、θ、
Tajima's D 等，并做过滤、LD 剪枝、索引与格式转换。

## 安装

```bash
cargo install --path crates/reed-cli     # 需要 Rust 1.91+
reed --version
```

## 用法

```bash
# 窗口化的 π、dxy、Fst、θ、Tajima's D
reed popg allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out

# 只含变异位点的输入
reed popg clean.vcf.gz --pop pops.txt --variant-only \
    --fst --window 100k --step 50k -o out

# QC 报告与汇总
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# 直方图
reed stats clean.vcf.gz --hist qual,dp-site,gq,INFO/DP -o qc

# 改样本名、染色体名，补缺失变异 ID
reed rename in.vcf.gz -s "s0 sampleA" --rename-chrs chrs.txt \
    --set-missing-var-ids '%CHROM:%POS' -Oz -o renamed.vcf.gz

# 位点与基因型过滤，输出带索引
reed filter in.vcf.gz -o out.vcf.gz -Oz6 --geno 0.1 --maf 0.05 \
    --min-dp 5 --min-gq 20 -W

# LD 剪枝
reed filter in.vcf.gz --indep-pairwise 50 5 0.2 -Oz6 -o pruned.vcf.gz

# CSI 索引
reed index in.vcf.gz

# 系统发育矩阵、0/1/2 矩阵、hapmap 与参考驱动的 MSA
reed convert in.vcf.gz --phylip --fasta --nexus --nexus-binary --012 --hapmap -o matrices
reed convert in.vcf.gz --msa --reference ref.fa --region chr1:1-100000 -o aln

# hapmap 转回 VCF（-Oz 写 BGZF；不给 -o 时输出到 stdout）
reed convert panels.hmp.txt --vcf -Oz -o panels.vcf.gz

# 投影 SFS
reed convert clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out
```

## 命令

- `filter` .. 选择、基因型掩码、QC 阈值、LD 剪枝
- `stats` .. 位点级与样本级 QC：freq、hardy、het、missing、depth、summary
- `popg` .. 窗口化统计：π、dxy、Fst、PBS、SFS、f 统计、LD
- `convert` .. 系统发育矩阵、基因型矩阵、hapmap、MSA、SFS
- `rename` .. 标识符编辑：样本名、染色体名、变异 ID、header 行
- `index` .. 建 CSI/TBI 索引

运行 `reed` 列出全部命令；`reed <command> --help` 列出全部选项。

## reed popg [OPTIONS] FILE

按窗口、按 contig 或全基因组的群体遗传统计量。

### 统计选项

| 选项 | 说明 |
|---|---|
| `--pi` | 每群体 π（`PREFIX_pi.txt`） |
| `--dxy` | 每群体对 dxy（`PREFIX_dxy.txt`） |
| `--fst [wc84\|hudson]` | 全部群体对的 Fst，默认 wc84（`PREFIX_fst.txt`） |
| `--fst-global` | 与 `--fst` 同给：全部群体的全局 WC84 Fst（`PREFIX_fst_global.txt`） |
| `--theta` | 每群体 Watterson θ（`PREFIX_watterson_theta.txt`） |
| `--tajima-d` | 每群体 Tajima's D（`PREFIX_tajima_d.txt`） |
| `--pbs A,B,C` | 指定三重的 PBS，可重复；`--pbs-normed` 写 PBSn1（`PREFIX_pbs.txt`） |
| `--sfs [folded]` | 投影 SFS，默认 unfolded；需要 `--proj N`（`PREFIX_sfs/dadi/`、`PREFIX_sfs/fastsimcoal2/`） |
| `--f2 A,B`、`--f3 C;A,B`、`--f4 A,B;C,D`、`--dstat A,B;C,D`、`--f4-ratio A,B;C,D;X` | Patterson f 统计量，5 Mb 区块刀切（`PREFIX.fstats`） |
| `--ld-decay [MAX_BP]` | MAX_BP 内（默认 100k）的配对 r²，1 kb 分箱（`PREFIX.ld.decay`）；`--ld-pairs` 另写逐对表（`PREFIX.ld`） |

### 窗口选项

| 选项 | 说明 |
|---|---|
| `--window SIZE` | 窗口大小（bp），接受 `100k`、`1M`；默认 100000 |
| `--step STEP` | 窗口步长；默认等于窗口大小 |
| `--levels LIST` | 输出行级：`window`、`contig`、`genome`；默认 `window` |
| `--skip-empty` | 丢弃没有位点落入的窗口 |
| `--pop FILE` | `样本<TAB>群体` 行、无表头；窗口模式必需 |
| `--variant-only` | 声明全变异输入：文件不含的窗口位置计为纯合且全部样本有调用 |
| `--bypass-invariant-check` | 只在文件含有的位点上计算（高估 π/dxy） |
| `--contig-lengths FILE` | `.fai` 式 `名称<TAB>长度`，用于缺 `##contig` 行的头 |

## reed stats [OPTIONS] FILE

位点级与样本级的 QC 报告。

| 选项 | 说明 |
|---|---|
| `--freq` | 逐位点等位频率（`PREFIX.afreq`） |
| `--hardy` | HWE 精确检验，观测与期望杂合度（`PREFIX.hardy`） |
| `--het` | 逐样本近交系数 F（`PREFIX.het`） |
| `--missing` | 逐样本与逐位点缺失率（`PREFIX.smiss`、`PREFIX.vmiss`） |
| `--depth` | 逐样本与逐位点平均深度，读 FMT/DP（`PREFIX.idepth`、`PREFIX.ldepth.mean`） |
| `--summary` | QC 汇总（`PREFIX.stats`） |
| `--hist LIST` | 每个指标一份直方图 `PREFIX.<m>.hist`：内建 `qual`、`dp`、`dp-site`、`gq`、`idl`、`af`，或任意数值 `INFO/<TAG>`/`FMT/<TAG>`；整数指标附 `#mean`/`#median`/`#p95` |
| `--singletons` | 逐等位 singleton 与私有 doubleton（`PREFIX.singletons`） |
| `--pop FILE` | `样本<TAB>群体` 行、无表头；与 `--freq`/`--missing` 同给时写分层文件 `PREFIX.afreq.pop`、`PREFIX.vmiss.pop`（每群体一行），与 `--singletons` 同给时 `PREFIX.singletons` 加 POP 列 |

`PREFIX.dp.hist`：

```
bin_lower	bin_upper	count	percent
5	5	12	1.234568
6	6	48	4.938272
#mean	8.104938
#median	7
#p95	21
```

## reed rename [OPTIONS] FILE

编辑样本名、染色体名、变异 ID 与 header 行。

| 选项 | 说明 |
|---|---|
| `-s "old new,…"`, `-N FILE` | 样本改名（文件里是 `old<TAB>new` 行） |
| `--replace-sample-name PATTERN` | 字面改名（`pre*` 去前缀、`*suf` 去后缀，否则删除该文本） |
| `--rename-chrs FILE` | `old<TAB>new` 染色体名；`##contig` 行与 CHROM 列同改 |
| `--chr-style plain\|ucsc\|mt-is-m` | 染色体名规范化（无 `chr`、带 `chr`、`chrM`） |
| `--set-id PATTERN` | 用 `%CHROM %POS %REF %ALT` 模板重写全部 ID |
| `--set-missing-var-ids PATTERN` | 只补 `.` 的 ID |
| `--update-name FILE`、`--idspace-to SEP`、`--rm-dup-ids`、`--clear-ids` | ID 映射、换分隔符、重复加后缀、全部置 `.` |
| `-h FILE`、`-f FILE`、`--add-meta LINE`、`--remove-meta KEY` | 换 header、按 `.fai` 补 `##contig` 长度、增删 header 行 |

输出：`-o FILE`（默认 stdout）、`-O v|z[level]`。

## reed filter [OPTIONS] FILE

对 VCF/BCF 做子集、过滤与掩码；LD 剪枝；写出 VCF 或 BCF。

| 选项 | 说明 |
|---|---|
| `-r LIST` / `-R FILE` | 区域；需输入旁有 `.csi`/`.tbi`，用 `reed index` 建 |
| `-t LIST` / `-T FILE` | 目标位点；流式扫描，不需索引 |
| `-s/-S`、`--keep`、`--remove` | 样本子集；`--force-samples` 把未知样本名降级为警告 |
| `-v/-V`、`--snps-only` | 按变体类型选位点；`--snps-only` 只留 SNP |
| `-m/-M INT` | 记录等位数下限/上限；默认 `-M 2`，拒绝多等位记录 |
| `-f LIST` | 只留 FILTER 列含 LIST 中值的位点 |
| `-i/-e EXPR` | 表达式纳入/排除位点 |
| `-q/-Q/--min-af/--max-af FLOAT[:type]`、`--geno`、`--maf`、`--max-maf`、`--mac`、`--max-mac`、`--hwe <P [midp]>` | 位点过滤，掩码之后应用：等位基因频率上下限（`:type` 取 `nref\|alt1\|minor\|major\|nonmajor`）、缺失率、MAC、HWE |
| `--min-dp`、`--max-dp`、`--min-gq` | 从 FMT/DP、FMT/GQ 掩码基因型；被掩基因型对位点过滤按缺失计 |
| `--set-gts .` | 配 `-i/-e`，把命中的 FORMAT 基因型置 `.` |
| `--mind FLOAT` | 丢弃缺失率超过 FLOAT 的样本 |
| `--indep-pairwise W[kb] [S] R2` | LD 剪枝，写出保留位点：`W` 裸数字是位点个数、带 `kb` 后缀是区间（步长固定 1）；`S` 默认 1；`R2` 是 r² 上限 |
| `-I/--no-update` | 不更新 `INFO/AC`、`INFO/AN`（默认随样本集更新；无 GT 列的记录保留输入值） |
| `-x/--remove-info LIST` | 删除 INFO tag；默认删 `AC,AN,AF,MLEAC,MLEAF,ExcessHet,InbreedingCoeff,DP`，`-x ''` 保留全部 |
| `--fill-tags LIST` | 重算 AC、AN、NS、MAF 或 HWE（大小写不敏感） |
| `-O v\|z\|b\|u[0-9]`、`-o FILE`、`-W[=csi\|tbi]` | 输出类型（可带压缩级别）、输出文件（默认 stdout）、写索引 |
| `--report`、`--timing`、`--no-version`、`--no-command` | JSON 运行报告、各阶段耗时、抑制 `##reed_*` 行 |

## reed index [OPTIONS] FILE

| 选项 | 说明 |
|---|---|
| `-c/--csi` | 写 CSI 索引（默认） |
| `-t/--tbi` | 写 TBI 索引 |
| `-o FILE` | 索引写到别处 |
| `-f` | 覆盖已有索引 |

输入必须是文件。

## reed convert [OPTIONS] FILE

把 VCF/BCF 转成其他格式；hapmap 输入可转回 VCF。

| 参数 | 输出 |
|---|---|
| `--phylip` / `--fasta` / `--nexus` | `PREFIX.phy` / `.fasta` / `.nexus` 核苷酸矩阵 |
| `--nexus-binary` | `PREFIX.bin.nexus` |
| `--012` | `PREFIX.012` + `.012.pos` + `.012.indv` |
| `--hapmap` | `PREFIX.hmp.txt` |
| `--vcf [-O v\|z[level]]` | hapmap 输入转回 VCF |
| `--msa --reference ref.fa [--region CHR:B-E]` | `PREFIX_<contig>.fasta` |
| `--sfs --proj N[,N…] [--pop FILE]` | `PREFIX_sfs/dadi`、`PREFIX_sfs/fastsimcoal2` |
| `--outgroup NAME`、`--min-samples-locus N`、`--resolve-iupac`、`--write-used-sites` | 调节项，见 `--help` |

差异记录在 `docs/convert-parity.md`。

## 表达式（`-i/-e`）

`-i` 保留为真的位点，`-e` 丢弃为真的位点。

| 类别 | 内容 |
|---|---|
| 列 | `CHROM`、`POS`、`ID`、`QUAL`、`FILTER`、`TYPE` |
| 标签 | `INFO/TAG`、`FMT/TAG`；裸 `TAG` 先查 INFO |
| 计算量 | `N_ALT`、`N_SAMPLES`、`N_MISSING`、`F_MISSING`；`AC`、`AN`、`AF`、`MAC`、`MAF`（INFO 有值取 INFO，否则按基因型） |
| 运算符 | `==`、`=`、`!=`、`<`、`<=`、`>`、`>=`；正则 `~`、`!~`；`&&`、`\|\|`、`!`、括号 |
| 缺失 | `"."`；与缺失的比较为假 |

```
reed filter in.vcf.gz -i 'QUAL>30 && MAF>=0.05' -Oz -o out.vcf.gz
reed filter in.vcf.gz -e 'FMT/DP<5' --set-gts . -Oz -o out.vcf.gz

# GATK SNP 硬过滤
reed filter in.vcf.gz -e 'QD<2 || QUAL<30 || SOR>3 || FS>60 || MQ<40 || MQRankSum<-12.5 || ReadPosRankSum<-8' -Oz -o out.vcf.gz
# GATK indel 硬过滤
reed filter in.vcf.gz -e 'QD<2 || QUAL<30 || FS>200 || ReadPosRankSum<-20' -Oz -o out.vcf.gz
```

## 输出

`filter`、`rename`、`convert` 与 `index` 写出 `-o` 指定的文件；`stats` 与
`popg` 把 `-o` 当前缀（`stats` 无 stdout 模式）。`popg` 另写 `PREFIX.log`。

`PREFIX_pi.txt`：

```
pop	chromosome	window_pos_1	window_pos_2	avg_pi	no_sites	count_diffs	count_comparisons	count_missing
popA	chr1	1	10000	0.001234	95	12	812	3
```

分母为零的估计值打印 `NA`。

## 作为库使用

`reed` crate 是公开 API。

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

已有 `样本 -> 群体` 映射的可以直接传入，不必写群体文件：

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

crate 根导出：`run_filter`、`run_popgen`、`run_stats`、
`run_rename`、`run_index`、`run_convert`，请求类型（`PopgenRequest`、
`FilterSpec`、`StatsReportRequest`、`RenameSpec`、`ConvertRequest`），
报告（`PopgenReport`、`FilterReport`、`StatsOutcome`、`RenameReport`、
`IndexOutcome`、`ConvertReport`），writer（`reed::windows`、`reed::reports`、
`reed::histogram`、`reed::sfs`、`reed::fstats`、`reed::ld`、`reed::summary`），
以及单一 `reed::Error`。`PopgenRequest::to_args(&group_names)` 与
`FilterSpec::to_args()` 把请求渲染回 CLI 拼写；CLI 拼写不出来的请求（per-bp
kind、`Regions` grid）返回错误。只有 crate 根有稳定性承诺，`crates/` 下
的下层 crate 是内部实现；错误用 `Error::kind()` 分类，消息走
`Display`/`source()`。

## 引用

reed 实现的是已发表的估计量，请引用它使用的方法：

- π：Nei & Li 1979；dxy：Nei 1987；θ：Watterson 1975；D：Tajima 1989。
- Hudson Fst：Hudson 1992 与 Bhatia et al. 2013；WC84 Fst：Weir & Cockerham 1984。
- f 统计量：Patterson et al. 2012；PBS：Yi et al. 2010。
- π、dxy、θ、D 的缺失数据处理：Korunes & Samuk 2021。

## 许可

MIT OR Apache-2.0。见 [LICENSE-MIT](LICENSE-MIT) 与
[LICENSE-APACHE](LICENSE-APACHE)。
