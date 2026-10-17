# 01 语义对照与差异台账（doc_v14：命令面重组 vs bcftools / vcftools / plink2 / seqkit）

> 证据基准：reed 侧 = 2026-10-17 工作树（v13 落地后，commit `4dcbd24`），
> 四命令实读；参照侧 = 本工作区源码副本（bcftools 1.9 dev、vcftools
> 0.1.18、plink-ng 2.0.0-a.7LM、seqkit 2.10.1 实测 `--help`）。

---

## §1 reed 命令面现状（拆分前的完整核实）

### §1.1 四个子命令与 flag 计数【据】

| 命令 | 定义处 | flag 面 | 备注 |
|---|---|---|---|
| `filter` | `cli/mod.rs:76`（FilterCli） | Region/Sample/Site 三共享组 + `--mind` + INFO 族（`-I/-x/--fill-tags`）+ 输出族（`-O/-o/-W`）+ 运行族 | v12 后稳定 |
| `stats` | `cli/mod.rs:194`（StatsCli） | Region/Sample/Site 三共享组 + QC 族（`--freq/--hardy/--het/--missing/--depth/--summary`）+ 窗口族（`--window/--step/--levels/--skip-empty/--bed-file/--sites-file`）+ 群体族（`--pop/--force-samples/--variant-only/--bypass-invariant-check/--contig-lengths`）+ 统计族（`--pi/--dxy/--theta/--tajima-d/--fst/--fst-global/--pbs/--f2/f3/f4/dstat/f4-ratio/--block-size/--ld-decay/--ld-pairs/--sfs/--proj/--sfs-int`） | **拆分对象**：QC 族 vs 其余全部 |
| `index` | `cli/mod.rs:167` | `-c/-t/-o/-f` + 运行族 | 稳定 |
| `convert` | `cli/mod.rs:376`（ConvertCli，v13） | 九目标（sfs/phylip/fasta/nexus/nexus-binary/012/hapmap/msa/vcf）+ Sample 组 + 调节族 + `-o/-O` | 本轮只改 help 归组 |

### §1.2 库层早已存在的分界线【据】

- 门面函数：`run_stats`（窗口统计，`stats/driver.rs:28`）与
  `run_reports`（QC 报告，`stats/writers/reports.rs:733`）是两个独立入口，
  请求类型分别是 `StatsRequest` 与 `ReportRequest`——**CLI 的 stats 是
  一个命令装两个门面**。
- 翻译层分叉：`translate/stats.rs` 里 reports 路线与
  `windowed_route || sfs || ld_decay` 路线的判定（v13 后 SFS 双入口）。
- 累加器基建：`QcReportAcc`（`reed-stats/src/summary/report.rs:26`）每个
  QC run 已积累 QUAL（`QualAcc`，SNV/ts/tv/indel 分列）、FMT/DP、INFO/DP、
  IDD（插入/缺失各一）、AF 直方图；`IntHistogram`
  （`summary/histogram.rs:14`）提供精确整数分位与可合并合并。
  **`--hist` 的主体是渲染层**，新解码只有 FMT/GQ 与通用标签。

### §1.3 共享参数组机制【据】

`filter_flags.rs` 定义 `RegionArgs`（`-r/-R/-t/-T`）、`SampleArgs`
（`-s/-S/--keep/--remove`）、`SiteFilterArgs`（`-v/-V/--snps-only/-m/-M/
-f/-i/-e/--min-dp/--max-dp/--min-gq/--geno/--maf/--mac/--hwe/--indep-pairwise`
等）；`FilterCli`/`StatsCli`/`ConvertCli`（仅 Sample 组）各自 flatten。
语义分叉在物理上不可能（同一结构体字段），`cli_contract.rs` 锁共享 help
文本。**本轮判据成文：选择参数共享，动作参数各归**（02 §5）。

---

## §2 参照工具台账

### §2.1 bcftools（25 个注册子命令，`main.c:79-193`【据】）

`main.c` 的命令表注册 25 个子命令（index、tabix、annotate、concat、
convert、head、isec、merge、norm、plugin、query、reheader、sort、view、
call、consensus、cnv、csq、filter、gtcheck、mpileup、polysomy、roh、
stats、som；其中 tabix 隐藏不广告、plugin/polysomy/som 受编译门控，
help 常显约 21 个）。
bcftools 自身用**注释行分组**（"VCF/BCF manipulation" / "VCF/BCF
analysis"）+ 扁平命令——与本轮方向一致（分组是 help 呈现，不是命令嵌套）。

与本轮相关的命令语义：

| bcftools 命令/flag | 行号 | 语义 | 对 reed 的映射 |
|---|---|---|---|
| `reheader -h FILE` | `reheader.c:708` | 整体替换 header（记录字节直通，不解码） | `rename -h` 同名同义 |
| `reheader -s PAIRS` / `-N FILE` | `reheader.c:710-712` | 样本改名（"old new" 对，逗号分隔或文件） | `rename -s/-N` |
| `reheader -f FAI` | `reheader.c:706` | 按 .fai 更新 `##contig` 行 | `rename -f`（对 reed 尤其值：窗口统计硬性要求 contig 长度） |
| `annotate --rename-chrs FILE` | `vcfannotate.c:3951` | 染色体映射 `old<TAB>new`，同时改头与记录 | `rename --rename-chrs` |
| `annotate -I/--set-id EXPR` | `vcfannotate.c:3934` | 模板写 ID（`%CHROM:%POS` 等） | `rename --set-id` |
| `stats` | `vcfstats.c` | QC 汇总：SN/TSTV/SiS/AF/QUAL/ST/DP/IDD/PSC 节 | reed `stats --summary` 已对拍；`--hist` 把其中直方图节独立成 TSV |
| `stats --html` 直方图数据 | vcfstats.c plot 节 | QUAL/DP/IDD 直方图仅嵌在 .stats 文本里 | B14-2：reed 独立 `.hist` TSV + 精确分位 |
| `query -f FORMAT` | `vcfquery.c` | 自定义列提取 | 不进本轮（backlog，v15+ 候选） |

**bcftools 没有的**（reed 差异化）：任意 `INFO/TAG`/`FMT/TAG` 直方图、
整数精确分位数（bcftools 固定 bin 无分位）、`GQ` 直方图。

### §2.2 vcftools（110 个长 flag 的单命令，`parameters.cpp` 去重实测【据】）

vcftools 是"一个命令装下所有统计"的反面教材——本轮拆分的动机佐证。
QC 族拼写对照（reed 已实现，零改动）：

| vcftools | reed stats（拆分后不变） |
|---|---|
| `--freq` / `--freq2` | `--freq`（plink2 afreq 布局） |
| `--hardy` | `--hardy`（plink2 布局） |
| `--het` | `--het`（plink2 布局） |
| `--missing-indv/--missing-site` | `--missing`（smiss+vmiss） |
| `--depth` / `--site-depth` / `--site-mean-depth` | `--depth`（idepth+ldepth.mean） |
| `--singletons` | C14-3（拟进本轮，vcftools 布局：CHR/POS/A/R/A1/A2/TYPE/POP） |

### §2.3 plink2（378 flag，实测 `--help` 2026-01 版【实】）

与本轮相关：

| plink2 flag | 语义 | 对 reed 的映射 |
|---|---|---|
| `--update-name FILE` | 映射改 variant ID | `rename --update-name` |
| `--set-all-var-ids` / `--set-missing-var-ids` | 模板 ID（`@`/`#`/`$r`/`$a` 占位） | `rename --set-id` / `--set-missing-var-ids`（bcftools `%` 模板拼写优先） |
| `--idspace-to SEP` | ID 分隔符替换 | `rename --idspace-to` |
| `--rm-dup` | 去重复（ID 维度） | `rename --rm-dup-ids`（加后缀而非丢弃，B14-6） |
| `--rename-chrs FILE` | 染色体改名 | 同 bcftools 拼写 `--rename-chrs` |
| `--output-chr MT/0/M/n` | 染色体命名约定 | `rename --chr-style`（枚举子集，B14-4） |
| `--update-ids/--update-map/--update-alleles` | 改样本 ID/坐标/等位 | **不进**（D14-4：改坐标破坏排序不变量；改等位须重算 INFO） |
| `--freq x/--missing/--hardy/--het` | QC 报告 | reed stats 已有，popg 分层 QC 用 `--within` 思路的 `--pop` |
| `--pca` / `--make-king` | PCA / 亲缘 | **不进本轮**（v15 候选，02 §1 网格预留位） |

### §2.4 seqkit 2.10.1（分组 help 的形态参照【实】）

seqkit `--help` 的分组机制：顶层命令全部扁平（`grep`/`rmdup` 是顶层
动词），分组只出现在 help 的标题行（"Commands for Set Operation:" 等
8 个 "Commands for …" 组 + 1 个 "Additional Commands" 组，共 9 组）。
**这是本轮分组 help 的直接范本**：命令路径不变、可发现性靠
分组标题。seqkit 的 `rename`（去重 ID）/`replace`（正则改名）动词先例
支持 `rename` 命名。`genautocomplete`（shell 补全生成）是
Miscellaneous 组的成员先例（reed 用 `clap_complete`，本轮不做）。

---

## §3 拆分判据与搬迁表

**判据一句话**：报告行的单位是**位点或样本** → `stats`；是**群体（对/
三元组）× 窗口或全基因组** → `popg`。

### §3.1 搬迁表（现 StatsCli 的每个 flag 去向）【规】

| 去向 | flag | 备注 |
|---|---|---|
| **stats（留）** | `--freq --hardy --het --missing --depth --summary` | 拼写与输出布局零变化 |
| **stats（新）** | `--hist LIST` `--geno-counts`（backlog，C14-3）`--singletons`（默认进轮，C14-3） | 02 §3 |
| **stats（共享组）** | Region/Sample/Site 三组 | 照旧 flatten |
| **stats（保留降级）** | `--pop FILE` `--force-samples` | **不是新增**：现 stats 已有这两个 flag（窗口模式必需），拆分后 stats 侧语义降级为"分层 QC 可选"；popg 侧为窗口模式必需（现状）。同名 flag 双命令并存，语义按命令定（分层 vs 必需） |
| **stats（hidden 报错族）** | `--mind` `--max-missing-sample` | 现状即 hidden 解析报错（`unsupported_stats_mind`/`removed_max_missing_sample`），拆分后保持 |
| **stats（运行/输出族）** | `-o/--output` `--threads` `--tmp-dir` `--memory` `--verbosity` | 两命令各自保留（`--tmp-dir` 帮助文本过期的 D-3 项顺手修正） |
| **popg（搬）** | `--pi --dxy --theta --tajima-d --fst --fst-global --pbs --pbs-normed --f2 --f3 --f4 --dstat --f4-ratio --block-size --ld-decay --ld-pairs --sfs --proj --sfs-int` | 拼写不变，仅换所属命令 |
| **popg（搬）** | `--variant-only --bypass-invariant-check --contig-lengths` | 清单族（`--pop/--force-samples` 见上：popg 侧照旧） |
| **popg（搬）** | `--window --step --levels --skip-empty --bed-file --sites-file` | 窗口族 |
| **popg（共享组）** | Region/Sample/Site 三组 | 照旧 |
| **popg（运行/输出族）** | `-o/--output` `--threads` `--tmp-dir` `--memory` `--verbosity` | 同 stats |
| **遗留 stats 上的群体 flag** | 上三行的全部 | 解析期报错 + 改写提示（removed-flag 模式，§3.3） |

> 搬迁表与 `StatsCli` 的 47 个字段（2026-10-17 `awk` 全量提取核对）一一
> 对应；`input` 为位置参数两边各有。

边界裁决（已定）：

- `--ld-decay/--ld-pairs` → popg。证据：LD 平面是 group-major
  （`translate/stats.rs:409`，裸给时全体样本为一群体），行 = 群体×距离
  bin，符合判据。
- `--indep-pairwise`（Site 组成员）→ 两命令照旧都有（filter 语义 =
  剪枝后输出；stats/popg 语义 = 在剪枝后位点上统计）。现状 stats 报
  "not supported yet"（`translate/stats.rs:274`），拆分时该报错文本改为
  引导 popg。
- `--summary` → stats。行 = 全文件汇总，无群体维度。
- `--sfs` → popg（搬）+ convert 既有入口不动（D13-6 双入口契约更新为
  popg/convert 双入口，stats 遗留入口报错引导）。

### §3.2 命名（C14-1，已拍板：`popg`）

主名 `popg`（population genetics 的缩写；用户 2026-10-17 拍板），
visible alias `popgen` 兜底自释拼写。help 的 about 写 "windowed
population-genetics statistics"。`qc` 名**不采用**（cli.md §1 已登记
其弃用史）。

### §3.3 遗留重定向（removed-flag 模式）【规】

旧 `reed stats in.vcf.gz --pi --pop pops.txt` 在 0.4.0 起：

```
reed: --pi moved to reed popg: run 'reed popg in.vcf.gz --pop pops.txt --pi …'
```

exit 2（usage 错，与 `--max-missing` 家族同款：解析期报错、给精确改写、
无静默别名）。QC flag 与 popg flag 拼写冲突为零（两集合不相交，搬迁表
可全量核实）。

---

## §4 reed 侧可复用面盘点

| 资产 | 位置 | 本轮用途 |
|---|---|---|
| `IntHistogram`/`QualDist`/`AfAcc`/IDD 累加器 | `reed-stats/src/summary/{histogram,report}.rs` | `--hist` 直接消费；新增 GQ 与通用标签累加进 `QcReportAcc`（或并立 `HistRequest` 计划） |
| `QcReportAcc` 每run 已积累的直方图 | `report.rs:26-53` | qual/dp/dp-site/idl/af 六 metric 零新计算，纯渲染 |
| `run_reports`/`ReportRequest` | `reports.rs:733/38` | 拆分后 stats 驱动器（改名 `run_stats`，02 §6） |
| `run_stats`/`StatsRequest` | `driver.rs:28` | popg 驱动器（改名 `run_popgen`） |
| removed-flag 解析器模式 | `filter_flags.rs:13-40` | popg flag 在旧 stats 上的重定向 |
| `clap` help 模板 + `Command` 元数据 | `cli/mod.rs` | 分组渲染器：`get_subcommands()` 遍历 + `GROUPS` 名单分桶（02 §1） |
| `edit_vcf_line`/`LineRef` 字节编辑 | `reed-io`（v12 后 `edit_vcf_line` 无失败路径） | `rename` 的 CHROM/ID 列替换主体 |
| header 渲染（VCF 文本头 + BCF 头块） | `reed-io/src/header.rs` | `rename -h/-f/--add-meta/--remove-meta` 与 BCF 字典重写 |
| `MissingRate` 精确有理数模式 | `reed-filter spec` | （风格参照）重定向报错给精确改写命令 |
| `open_regions`/facts | `reed/src/input.rs` | rename 的输入开路与样本名校验 |
| 输出/压缩路径（`-O v|z|b|u` + `-W`） | filter 既有 | `rename -o/-O` 直接复用 |
| `cli_contract.rs` 契约测试框架 | `reed-cli/tests` | 分组守卫、重定向报错文本、`--hist`/rename 契约 |

---

## §5 差异分档表（D14 / B14 / C14）

### D14（不做就是错的——对拍会挂或契约破裂）

| # | 约束 | 规范处 |
|---|---|---|
| D14-1 | **QC flag 零改动**：`--freq/--hardy/--het/--missing/--depth/--summary` 的拼写、输出文件与列布局、退出码在拆分后逐字节不变（gold `run_vcftools.sh`/`run_bcftools_af.sh`/`compare_reports.py` 只换命令词即须全绿） | 02 §2.2 |
| D14-2 | **重定向不是别名**：遗留 stats 上的群体 flag 解析期报错 + 给出含完整改写命令的错误文本，无静默转发（`--max-missing` 家族模式，cli.md §4.4 同款纪律） | 02 §2.4 |
| D14-3 | **SFS 双入口契约更新**：`convert --sfs` ≡ `popg --sfs` 逐字节（原 stats --sfs 的等价门改挂 popg）；`build_sfs_request` 仍单一实现（D13-6 延续），`translate` 层不得出现第二份 | 02 §2.3 |
| D14-4 | **`rename` 不改 POS/REF/ALT/基因型**：只编辑 header 行、CHROM 列、ID 列、样本名；改坐标破坏排序不变量、改等位须重算 INFO，两者都是 out of scope（plink2 `--update-map/--update-alleles` 不收） | 02 §4.1 |
| D14-5 | **`rename --rename-chrs` 必须同时改 `##contig` 行与记录 CHROM 列**：只改一头会使索引语义与头的 contig 表断裂（bcftools annotate 同款双改） | 02 §4.3 |
| D14-6 | **分组渲染单源**：`GROUPS` 常量是唯一分组清单；守卫测试断言 `Command` 每变体恰好入一组；README 命令表由同一常量生成或以测试锁定与 help 一致（防 v12 §0.5"多份实现各自漂移"） | 02 §1.4 |
| D14-7 | **`--hist` 累加器定律**：新累加器（GQ、通用标签）满足 merge 结合律/交换律并入既有 chunk 流水线；threads ∈ {1,4} 输出 sha256 相同（§8.4 定律沿用） | 02 §3.4 |
| D14-8 | **help 分组不引入嵌套子命令路径**：`reed filter …` 等调用路径一个字符不变（cli.md §0.9 契约的直接推论；违反 = 破坏全部既有脚本） | 02 §1.1 |

### B14（有意分歧——登记 cli.md 差异表，不写对拍断言）

| # | 分歧 | 与参照的差异 |
|---|---|---|
| B14-1 | `--hist` 通用标签（`INFO/TAG`/`FMT/TAG` 任意数值标签出直方图） | bcftools stats 只对固定字段（QUAL/DP/IDD）出直方图且嵌在 .stats 文本里 |
| B14-2 | 整数 metric 附精确 mean/median/p95（`IntHistogram::quantile` 精确分位） | bcftools 无分位数、固定 bin |
| B14-3 | `--hist` 输出独立 `PREFIX.<m>.hist` TSV（bin_lower,bin_upper,count,percent），不嵌 .stats | bcftools 嵌入式 |
| B14-4 | `--chr-style` 有限枚举（`plain|ucsc|mt-is-m` 三态），不做任意映射 | plink2 `--output-chr` 四态 + 数字模式；有限枚举避免语义蔓延 |
| B14-5 | `--set-id` 多等位 ALT 用 `,` 连接（`%ALT` 展开为 `A,B`）；模板占位符用 bcftools `%` 拼写 | plink2 用 `@`/`#`/`$a` 拼写——不并轨，bcftools 优先（cli.md §0.2 优先级） |
| B14-6 | `--rm-dup-ids` 重复 ID 加 `_1/_2` 后缀保全部；plink2 `--rm-dup` 丢弃。丢弃位点属于"选择"，加后缀属于"编辑"——rename 只做后者 | plink2 |
| B14-7 | `rename` 对 BCF 输入的染色体/样本改名须重写字典与头块（全量重编码）；纯样本改名可走快路径（只重写头块、记录直通） | bcftools reheader 同款快路径，保住卖点 |

### C14（待拍板——默认已定，可推翻）

| # | 决策 | 默认 |
|---|---|---|
| C14-1 | ~~`popg` 主名 + `popgen` alias，还是反向~~ | **已拍板（2026-10-17）：`popg` 主名 + `popgen` visible alias** |
| C14-2 | 裸 `reed stats in.vcf.gz`（无统计 flag）报错 or 默认 `--summary` | 报错（保持显式家风；bcftools stats 无选项全做的肌肉记忆让位） |
| C14-3 | `--singletons`/`--geno-counts` 进不进首发轮 | `--singletons` 进（SiSAcc 原料已有，渲染层）；`--geno-counts` 留 backlog |
| C14-4 | `rename --order-samples`（按清单重排样本列）进不进首发轮 | 不进（样本列重排是 merge/concat 轮的前置，单独做语义半残） |
| C14-5 | `--hist` 浮点标签 bin 边界 | bcftools 同款固定边界（QUAL 0.1 步长先例）；参数化（`--hist-bins`）留 backlog |
| C14-6 | `rename` 后旧 `.csi/.tbi` 失效：警告 or 默认重建 | 警告 + 提示 `reed index`（与 filter 的 `-W` 显式哲学一致，不隐式写文件） |

---

## §6 守卫清单（落地后进 `ci/check-doc-drift.sh`，编号续 v13 的 18）

| # | 守卫 |
|---|---|
| 1 | `Command` 枚举变体 = {Filter, Stats, Poph, Index, Convert, Rename}，`popgen` 是 stats/popg 之外的 alias 不占变体（doc_v14 guard 1） |
| 2 | 每个变体恰好出现在 `GROUPS` 一组（doc_v14 guard 2/D14-6） |
| 3 | 遗留群体 flag 在 stats 上有 removed-flag 解析器（`--pi` 不得作为 stats 的合法 flag 静默存在）（doc_v14 guard 3/D14-2） |
| 4 | `run_popgen`/`run_stats`（新含义）在门面根导出，旧 `run_reports` 名不复存在（doc_v14 guard 4） |
| 5 | `build_sfs_request` 仍单一实现且 popg/convert 共用（D14-3，续 v13 guard 2） |
| 6 | `--hist` 契约测试存在：metric 集合、通用标签、TSV 列头（doc_v14 guard 6） |
| 7 | `--hist` 累加器定律测试存在（merge 结合/交换）（D14-7） |
| 8 | `rename` 不解码基因型：`crates/reed/src/rename/` 无 GenotypeScanner/GtCode 引用（D14-4 的结构性守卫） |
| 9 | `--rename-chrs` 同改头与记录的 e2e 存在（D14-5） |
| 10 | `rename` 拒绝改 POS/REF/ALT 的负向测试存在（或：相关 flag 不在 RenameCli 定义中——grep 守卫） |
| 11 | 重定向报错文本含 `reed popg`（doc_v14 guard 11） |
| 12 | `--singletons`（若 C14-3 进轮）对拍 vcftools 的 gold 脚本存在 |
| 13 | 分组 help 快照测试：`reed --help` 输出含五组标题（doc_v14 guard 13） |
| 14 | stats 的 `--pop` 分层 QC 有 e2e（分层 missingness/freq 列） |

---

## §7 对拍真值门设计（gold 脚本）

| 门 | 脚本 | 真值 | 断言 |
|---|---|---|---|
| 拆分等价 | `gold/popgen/run_split_parity.sh`（新） | 0.3.x 的 `reed stats` 输出 vs 0.4.0 的 `reed popg` | QC 六报告 + 窗口统计全部逐字节（六 flag 组合矩阵） |
| `--hist` 内建 | `gold/popgen/run_hist.sh`（新） | `bcftools stats` 的 QUAL/DP/IDD 直方图节 | reed `.hist` TSV 的 bin/计数与 bcftools 文本节一致（注意：bcftools 的 DP 节把 per-genotype 与 per-site 两维并列在同一节——列 4/5 基因型、列 6/7 位点，实测 `bcftools stats` 输出核实——因此对拍断言按列映射到 reed 的 `dp` 与 `dp-site` 两个 metric；af 除外——bcftools bin 语义按 init_iaf，B14-2 分位不比） |
| `--hist` 通用标签 | e2e（`crates/reed/tests/e2e_hist.rs`） | 构造 INFO/DP 浮点、FMT/GQ | 累加器定律 + 手算 bin 断言 |
| `rename` 样本/染色体/ID | `gold/popgen/run_rename.sh`（新） | `bcftools reheader -s`、`bcftools annotate --rename-chrs`、`--set-id` | 输出 VCF 逐字节（bcftools 重排序头行的差异登记为等价集比较，B14 表） |
| `rename` BCF 路径 | e2e | VCF 路径输出转 BCF 再读回 | 记录与头逐项相等 |
| `--singletons` | `gold/popgen/run_singletons.sh`（新，若 C14-3 进轮） | vcftools `--singletons` | 列布局/值一致（POP 列接 `--pop`） |
| 全量确定性 | 沿用 e2e 纪律 | — | threads {1,4} sha256（hist/rename 全 flag） |

bcftools 侧注意：`reheader -s` 不重排头行、`annotate --rename-chrs` 会
重写 contig 行顺序——等价集比较（v13 B13 同款 diff -r 语义放宽）而不是
裸 diff；差异点进 B14 表补登。

---

## §8 性能基线（C0-2 实测后填入 STATUS）

| 负载 | 0.3.x `stats` | 0.4.0 拆分后 | 门 |
|---|---|---|---|
| QC 全报告（`--freq --hardy --het --missing --depth --summary`）clean.vcf.gz | 待测 | 待测 | 拆分纯搬运，墙钟退化 ≤ 3% |
| QC + `--hist qual,dp,dp-site,gq,idl,af` | —（新） | 待测 | 相对 QC-only 增幅 ≤ 10%（GQ 解码为主要新成本） |
| `popg` 全统计（π/dxy/Fst/θ/Tajima/PBS + LD decay） | 待测（0.3.x stats 同命令） | 待测 | 同上 ≤ 3% |
| `rename`（染色体+ID+样本全平面）clean.vcf.gz | —（新） | 待测 | 对标 `bcftools annotate --rename-chrs --set-id` + `reheader -s` 三进程合计的 ≤ 1/3 墙钟（单遍 vs 三遍 IO） |

每步跑 `bench/suites/e2e/results.tsv` 追加一行（v13 惯例）。
