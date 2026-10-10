# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed 对 VCF、VCF.gz 与 BCF 文件计算群体遗传学统计量并执行硬过滤，把工作切分到
所有核心上。纯 Rust 实现：传输层用 noodles，BGZF 与 BCF 编解码是 reed 自己的，
不链接 GPL 代码。`reed` 二进制是产品，`crates/` 下的 crate 是它的内部结构；
只有 CLI 参数与输出文件格式是稳定性承诺。

## 安装

需要 Rust 1.91 或更新版本，仓库在 `rust-toolchain.toml` 里固定 1.99。

```bash
cargo install --path crates/reed-cli
reed --version   # reed 0.3.0 (64e5dce 2026-10-10, rustc 1.99.0)
```

也可以就地构建：

```bash
cargo build --release -p reed-cli     # 产物在 target/release/reed
```

不需要安装任何系统库：BGZF 与 BCF 编解码全部由 Rust 源码编译进来。

## 快速开始

```bash
# 全位点文件上的全基因组与窗口 π、dxy、Fst、θ、Tajima's D
reed stats allsites.vcf.gz --pop pops.txt --pi --dxy --fst --theta --tajima-d \
    --window 10k -o out

# 仅变异位点输入，成对与全局 Fst
reed stats clean.vcf.gz --pop pops.txt --variant-only \
    --fst --window 100k --step 50k -o out

# 逐位点与逐样本 QC 报告加 bcftools-stats 汇总，一遍读取
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# 位点与基因型过滤，并给输出建索引
reed filter in.vcf.gz -o out.vcf.gz -Oz6 --geno 0.1 --maf 0.05 --min-dp 5 --min-gq 20 -W

# LD 剪枝：一条命令直接写出剪枝后的 VCF
reed filter in.vcf.gz --indep-pairwise 50kb 5 0.2 -Oz6 -o pruned.vcf.gz

# CSI 索引
reed index in.vcf.gz
```

## 命令

| 命令 | 作用 | 对标 |
|---|---|---|
| `reed filter` | 硬过滤、基因型掩码、INFO 重算；保真的 VCF/BCF 输出 | bcftools `view`、`filter -S .`、`+fill-tags`、`annotate -x`；plink2 的 QC 阈值 |
| `reed stats` | π、dxy、Fst、Watterson θ、Tajima's D、PBS、SFS、f 统计、LD、逐位点与逐样本报告、bcftools-stats 汇总 | vcftools、pixy、plink2、easySFS、scikit-allel、`bcftools stats` |
| `reed index` | CSI（默认）或 TBI | `bcftools index` |

参数以 `reed <命令> --help` 为准，完整规范见 `docs/spec/cli.md`。

### reed stats

统计参数，可任意组合，每项各写一个输出文件：

| 参数 | 输出 |
|---|---|
| `--pi` | 每群体 π，`PREFIX_pi.txt` |
| `--dxy` | 每群体对 dxy，`PREFIX_dxy.txt` |
| `--fst [wc84\|hudson]` | 全部群体对的 Fst，`PREFIX_fst.txt`（方法体现在列头） |
| `--fst-global` | 全部群体的全局 WC84 Fst，`PREFIX_fst_global.txt` |
| `--theta` | 每群体 Watterson θ，`PREFIX_watterson_theta.txt` |
| `--tajima-d` | 每群体 Tajima's D，`PREFIX_tajima_d.txt` |
| `--pbs A,B,C` | 指定三元组的 PBS，`PREFIX_pbs.txt`；加 `--pbs-normed` 输出 PBSn1 |
| `--sfs` | 投影 SFS，写在 `PREFIX_sfs/dadi/` 与 `PREFIX_sfs/fastsimcoal2/`；需要 `--proj N` |
| `--f2 A,B`、`--f3 C;A,B`、`--f4 A,B;C,D`、`--dstat A,B;C,D`、`--f4-ratio A,B;C,D;X` | Patterson f 统计与分块 jackknife，`PREFIX.fstats` |
| `--ld-decay [MAX_BP]` | 1 kb 分箱的成对 r²，`PREFIX.ld.decay`；加 `--ld-pairs` 另写逐对表 `PREFIX.ld` |
| `--freq`、`--hardy`、`--het`、`--missing`、`--depth` | `PREFIX.afreq`、`PREFIX.hardy`、`PREFIX.het`、`PREFIX.smiss` + `PREFIX.vmiss`、`PREFIX.idepth` + `PREFIX.ldepth.mean` |
| `--summary` | bcftools-stats 布局的 QC 汇总，`PREFIX.stats` |

窗口参数：

| 参数 | 含义 |
|---|---|
| `--window <SIZE>` | 窗口大小（bp），接受 `100k`、`1M`（默认 100000） |
| `--step <STEP>` | 窗口步长（默认等于窗口大小） |
| `--levels <LIST>` | 输出层级：`window`、`contig`、`genome`（默认 `window`） |
| `--skip-empty` | 跳过没有位点落入的窗口 |
| `--pop <FILE>` | `样本<TAB>群体` 行，无表头；窗口模式与所有按群体的统计量都需要 |
| `--variant-only` | 声明输入仅含变异位点：文件中没有的位置按单态且全部调用处理 |
| `--bypass-invariant-check` | 只在文件含有的位点上计算，接受 π/dxy 分母偏高 |
| `--contig-lengths <FILE>` | `.fai` 风格的 `name<TAB>length` 文件，用于头里没有 `##contig` 行时 |

π、dxy、θ 与 Tajima's D 需要包含单态位点的分母。默认情况下 reed 会检查输入是否
带有不变位点，没有就报错停下；`--variant-only` 与 `--bypass-invariant-check`
是绕过这一检查的两种声明方式。

样本与区域选择在 `stats` 与 `filter` 中一致：`-r/-R` 区域（需要索引）、
`-t/-T` 目标（顺序扫描）、`-s/-S`、`--keep`、`--remove`。位点过滤
（`-v/-V`、`--snps-only`、`-m/-M`、`-f`、`-i/-e`、`-q/-Q`、`--geno`、
`--maf/--max-maf`、`--mac/--max-mac`、`--hwe`）决定统计量看到哪些位点。

### reed filter

选择，沿用 bcftools：`-r/-R`、`-t/-T`、`-s/-S`、`--keep`、`--remove`
（`--force-samples` 把未知样本名降级为警告）、`-v/-V`、`--snps-only`、
`-m/-M`、`-f`、`-i/-e`。

位点过滤，在掩码之后执行：`-q/-Q/--min-af/--max-af`，可带
`:nref|alt1|minor|major|nonmajor` 后缀；`--geno`、`--maf`、`--max-maf`、
`--mac`、`--max-mac`、`--hwe <P [midp]>`。

基因型掩码：`--min-dp`、`--max-dp`、`--min-gq` 依据 FMT/DP 与 FMT/GQ 掩码；
`-i/-e` 的 FORMAT 表达式配 `--set-gts .` 精确掩掉命中的基因型，等同
`bcftools filter -S .`。被掩码的基因型在后续位点过滤中按缺失计。

样本过滤：`--mind <FLOAT>` 丢弃缺失率超过阈值的样本。
`--indep-pairwise <W[kb] [S] R2>` 按 plink2 算法做 LD 剪枝并输出保留的位点；
窗口 `W` 可带 `kb` 后缀，步长 `S` 默认 1，`R2` 是 r² 上限。

INFO 处理：样本集改变时 reed 更新 `INFO/AC` 与 `INFO/AN`；`-I/--no-update`
关闭该行为，`-x/--remove-info LIST` 删除标签，`--fill-tags LIST` 重算标签。

输出：`-O v|z|b|u[0-9]` 选普通 VCF、BGZF VCF、BCF 或未压缩 BCF，可带压缩级别；
`-o` 指定文件（默认 stdout）；`-W[=csi|tbi]` 写索引；`--report` 写 JSON 运行报告；
`--timing` 记录分阶段耗时；`--no-version` 与 `--no-command` 关闭 `##reed_*`
溯源头行。

### reed index

`-c/--csi`（默认）或 `-t/--tbi`，`-o` 写到别处，`-f` 覆盖已有索引。输入必须是
文件：建索引需要可寻址的输入。

## 输入要求

- VCF、VCF.gz 或 BCF，按文件头几个字节自动识别。`-` 表示从 stdin 读 VCF 文本；
  `reed index` 只接受文件。
- 头里必须有 `##contig=<ID=…,length=…>` 声明。窗口统计要把 contig 的最后一个
  窗口裁剪到声明长度，索引需要 contig 顺序。头里没有这些行时用
  `--contig-lengths FILE` 补上。
- 记录必须按 contig 与位置排序。索引与窗口统计都依赖这一点。
- `-r/-R` 需要输入旁的 `.csi` 或 `.tbi`，用 `reed index` 生成；顺序扫描不需要。
- 基因型按二倍体处理。默认下 reed 拒绝多等位记录（`-M 2`）；对 `0/.` 这类
  半调用按等位计数，因此 `--geno` 与 `--mind` 的口径与 plink2 相同。

## 输出

`filter` 与 `index` 写 `-o` 指定的文件。`stats` 把 `-o` 当作前缀，每个统计量
各写一个文件，名称见上面的表。不给 `-o` 时，窗口统计以 TSV 打到 stdout，
报告类文件不写。

`filter` 的输出会带溯源头行（`##reed_version`、`##reed_command`、
`##reed_filterCommand`、`##reed_filter`）；`stats` 把这些写进 `PREFIX.log`。

窗口表格结构一致：每个群体（或群体对）每个窗口一行，含 `chromosome`、
`window_pos_1`、`window_pos_2`、估计值，以及它背后的计数（π 与 dxy 有
`no_sites`、`count_diffs`、`count_comparisons`、`count_missing`，Fst、PBS 与
f 统计有 `no_snps`）。分母为 0 的估计值打印 `NA`。

## 验证

正确性以参照工具为准，而不是自洽。下表的每一行在 `gold/` 下都有脚本，数字
来自 2026-10-10 的复审（记录在 `docs/rounds/v6/`）。

| 检查 | 参照 | 结果 |
|---|---|---|
| π、dxy（含不变位点） | pixy 2.2.4 | 28,470 个窗口，最大相对误差 4.8e-14 |
| Watterson θ、Tajima's D | pixy | 3.1e-14 / 8.9e-11 |
| Fst wc84 / hudson（成对） | pixy | 2.9e-15 / 2.2e-15 |
| 全局 WC84 Fst | vcftools weighted Fst | 2.2e-16 |
| PBS | scikit-allel | 逐列精确 |
| 投影 SFS | easySFS、dadi 2.4.4 | 夹具 1e-15，子集 5.7e-11 |
| f2/f3/f4/D 估计值 | scikit-allel | 1.3e-16 |
| LD r² 与剪枝集合 | plink2 | 全量保留位点逐个相同（235,403） |
| `--freq/--hardy/--het/--missing` | plink2 | 逐行一致 |
| `--depth` | vcftools | 逐行一致 |
| `--summary` | `bcftools stats` | 九节逐行一致 |
| gold 过滤规则 | bcftools | 1,768,759 位点，`check_parity` PASS，正文 md5 `7aedc59a…` |
| `--geno/--mind` 边界 | plink2、vcftools | N=10/7 与半调用上一致；子集 263,420 / 381,672 |
| BCF 对 VCF、mmap 对流式 | bcftools | 正文逐字节相同，CSI 相同 |

线程数不影响输出：测试在 1、2、3、7、56 线程上对三种执行策略分别哈希，结果相同。

## 性能

全量 `clean.vcf.gz`（7.0 GB，3,315,604 条记录，348 样本），56 线程，热缓存，
三次取中位数：

| 命令 | 墙钟时间 |
|---|---|
| `filter` gold 规则，VCF 输入，`-Oz6` | 8.85 s |
| `filter` gold 规则，BCF 输入，`-Oz6` | 9.00 s |
| `filter` gold 规则，`-Ob` 输出 | 9.98 s |
| `stats` π/dxy/θ/D/Fst，VCF / BCF | 1.78 s / 1.65 s |
| `stats` 报告 + 汇总 | 3.0 s 到 5.4 s |
| `index` VCF.gz / BCF | 1.63 s / 1.12 s |

读取端对普通文件做 mmap，并在 BGZF 块边界上切分压缩块，线程之间不需要重新
对齐边界。deflate 是 `filter` 的主要开销。小内存机器上用 `--memory <MB>`
限制在途缓冲。

## 已知限制

- 不支持多等位记录。把 `-M` 放宽到 2 以上时 reed 不会重算每个等位的统计量，
  默认因此把这类记录挡在外面。
- LD 剪枝在 `filter` 里，需要它的第二遍扫描。`reed stats` 会拒绝
  `--indep-pairwise` 并指向 `reed filter`。
- `--ld-decay` 目前也要求 `--pop FILE`，尽管 LD 与群体标签无关。
- `--window` 当前必须显式给 `--step`，帮助文本里"默认等于窗口大小"的默认值
  目前不被接受。
- f 统计的标准误按碱基分块 jackknife，scikit-allel 按变异数分块。估计值一致，
  标准误不可比。
- Tajima's D 用窗口平均等位数的方差作分母（pixy 口径），与 vcftools 固定 2N
  的做法不同。

## 构建与测试

```bash
cargo test --release --workspace
cargo build --release --workspace --all-targets
cargo clippy --release --workspace --all-targets -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
cargo fmt --all --check
bash ci/check-doc-drift.sh
```

lint 门覆盖库与二进制目标。测试目标还有未清的 pedantic lint，所以对
`--all-targets` 开 `-D warnings` 会报出来。

## 引用

reed 实现的是已发表的方法，并与上面的参照工具对拍。若在基于这些方法的工作中
使用 reed，请引用原始论文：

- Nei & Li 1979（π）；Nei 1987（dxy）；Watterson 1975（θ）；Tajima 1989（D）。
- Hudson 1992 与 Bhatia et al. 2013（Hudson Fst）；Weir & Cockerham 1984（WC84 Fst）。
- Patterson et al. 2012（f 统计）；Yi et al. 2010（PBS）。
- Korunes & Samuk 2021（π、dxy、θ、D 的缺失数据处理）。

## 仓库结构

```
crates/   八个 crate；唯一的二进制由 reed-cli 构建
gold/     验证：夹具与对拍脚本
bench/    度量套件与结果日志
docs/     规范（spec/ 为当前，rounds/ 为历史）
ci/       漂移守卫与 CI 入口
```

## 许可

MIT OR Apache-2.0。见 [LICENSE-MIT](LICENSE-MIT) 与
[LICENSE-APACHE](LICENSE-APACHE)。
