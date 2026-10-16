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
reed stats allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out

# 只含变异位点的输入
reed stats clean.vcf.gz --pop pops.txt --variant-only \
    --fst --window 100k --step 50k -o out

# 位点级与样本级 QC 报告 + bcftools-stats 汇总，一次扫描
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# 位点与基因型过滤，输出带索引
reed filter in.vcf.gz -o out.vcf.gz -Oz6 --geno 0.1 --maf 0.05 \
    --min-dp 5 --min-gq 20 -W

# LD 剪枝
reed filter in.vcf.gz --indep-pairwise 50kb 5 0.2 -Oz6 -o pruned.vcf.gz

# CSI 索引
reed index in.vcf.gz
```

`reed <command> --help` 列出全部选项。

## reed stats

统计开关任意组合，每个开关产出自己的文件。

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
| `--freq`、`--hardy`、`--het`、`--missing`、`--depth` | `PREFIX.afreq`、`PREFIX.hardy`、`PREFIX.het`、`PREFIX.smiss` + `PREFIX.vmiss`、`PREFIX.idepth` + `PREFIX.ldepth.mean` |
| `--summary` | bcftools-stats 布局的 QC 汇总，`PREFIX.stats` |

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

样本与区域选择在 `stats` 和 `filter` 中一致：`-r/-R` 区域（需要索引）、
`-t/-T` targets（流式）、`-s/-S`、`--keep`、`--remove`。位点过滤（`-v/-V`、
`--snps-only`、`-m/-M`、`-f`、`-i/-e`、`-q/-Q`、`--geno`、`--maf/--max-maf`、
`--mac/--max-mac`、`--hwe`）决定统计能看到什么。

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
use reed::{Populations, StatKind, StatsRequest};

let request = StatsRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::File("pops.txt".into()),
    ..StatsRequest::default()
};
let report = reed::run_stats(Some("in.vcf.gz".as_ref()), &request, &reed::budget::ThreadBudget::default())?;
let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false);
```

手上已经有 `样本 -> 群体` 映射的调用方可以直接传入，不必落一个群体文件：

```rust
let request = StatsRequest {
    kinds: vec![StatKind::Pi { group: 0, per_bp: false }],
    populations: Populations::Assignment(vec![
        ("sample1".into(), "popA".into()),
        ("sample2".into(), "popB".into()),
    ]),
    ..StatsRequest::default()
};
```

门面 re-export 所有需要的东西：`run_filter`、`run_stats`、`run_reports`、
`run_index`，请求类型（`StatsRequest`、`FilterSpec`、`ReportRequest`），
报告（`StatsReport`、`FilterReport`、`ReportsOutcome`、`IndexOutcome`），
writer（`reed::windows`、`reed::reports`、`reed::sfs`、`reed::fstats`、
`reed::ld`、`reed::summary`），以及单一 `reed::Error`（`kind()` 区分用法错误
与运行期失败）。`StatsRequest::to_args(&group_names)` 与
`FilterSpec::to_args()` 把请求渲染回 CLI 拼写，程序化构造请求的工具也能记录
等价的命令行；CLI 拼写不出来的请求（per-bp kind、`Regions` grid）由
`StatsRequest::to_args` 返回用法错误。只有这个门面是稳定性承诺；
`crates/` 下的下层 crate 是内部实现，`reed::Error` 通过封口变体把它们挡在
签名之外——用 `kind()` 分类，用 `Display`/`source()` 读消息。

## 已知限制

- 多等位记录被拒绝。`-M` 提到 2 以上时 reed 不会重算逐等位统计，默认保持拒绝。
- LD 剪枝在 `filter` 里，需要第二遍扫描。`reed stats` 拒绝 `--indep-pairwise`
  并指向 `reed filter`。
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
