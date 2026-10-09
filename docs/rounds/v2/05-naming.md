# 05 命名规范增量（对 v1 §13 的修订）

> v1 §13 中没有在这里提到的规则全部继续有效。本文只列修订、废止与新增项。

---

## 1. 废止

| v1 条款 | 状态 | 原因 |
|---|---|---|
| §13.15 配置键（TOML） | **废止** | 不再有配置文件（02-cli.md §0） |
| §13.14 中 `--spec`、`-c/--conf` 相关内容 | **废止** | 同上 |
| §13.1 "层名集合"中缺少 `cli` 的执行情况 | 重申 | 二进制必须放在 `reed-cli` crate，`reed` 库不依赖 clap |

---

## 2. CLI 参数命名（取代 v1 §13.14）【规】

1. **借用名字的前提是语义完全一致**。借用顺序：bcftools → plink2 → vcftools（02-cli.md §0）。
2. reed 自有参数使用 kebab-case 长参数（`--min-dp`、`--fill-tags`、`--sfs-missing`）。
3. **别名只给语义完全相同的参数**：plink2 的 `--vcf-min-dp/--vcf-max-dp/--vcf-min-gq/--min-alleles/--max-alleles`、
   vcftools/plink2 的 `--out`、pixy 的 `--populations`。vcftools 的 `--minDP/--maxDP/--minGQ`（DP/GQ 缺失时会掩码）、
   `--max-missing`（方向与 `--geno` 相反）语义不同，**不提供别名**。
4. **统计参数用统计量的通用名字**，不带窗口、不带工具名：`--pi`、`--dxy`、`--fst`、`--theta`、`--tajima-d`、`--pbs`、`--sfs`、
   `--freq`、`--hardy`、`--het`、`--missing`、`--depth`、`--summary`。窗口用 `--window/--step/--bed-file`；方法取值用作者式名字（`--fst wc84|hudson`、`--pbs-normed`），不出现软件名。
   不采用 vcftools 的"统计量 + 窗口"合一写法（`--window-pi 100000`），因为同时请求多个统计量时会重复写窗口。
5. 短参数与 bcftools 完全一致：

| 短参数 | 长参数 | 用途 |
|---|---|---|
| `-o` | `--output` | 输出 |
| `-O` | `--output-type` | 输出类型 |
| `-W` | `--write-index` | 写索引 |
| `-r` / `-R` | `--regions` / `--regions-file` | 索引区域 |
| `-t` / `-T` | `--targets` / `--targets-file` | 流式区域 |
| `-s` / `-S` | `--samples` / `--samples-file` | 样本 |
| `-i` / `-e` | `--include` / `--exclude` | 表达式 |
| `-f` | `--apply-filters` | FILTER |
| `-v` / `-V` | `--types` / `--exclude-types` | 变异类型 |
| `-m` / `-M` | `--min-alleles` / `--max-alleles` | 等位数 |
| `-x` | `--remove-info` | 删除 INFO |
| `-I` | `--no-update` | 不更新 AC/AN |
| `-q` / `-Q` | `--min-af` / `--max-af` | 频率阈值（bcftools view，含 `:minor` 等类型） |

6. **禁止**：`-v` 表示 verbose、`-t` 表示 threads、`-q` 表示 quiet。`--pop`、`--window` 等统计参数不设短参数（避免与 bcftools 将来的字母冲突）。
7. 环境变量不参与行为配置（删除 `REED_THREADS`）。日志级别使用 `RUST_LOG`（tracing 惯例）或 `--verbosity`。

---

## 3. 命令名【规】

`reed filter`、`reed stats`（全部统计，vcftools/plink2 模式）、`reed index`（02-cli.md §1，CLI-1 已裁决）。
库中驱动的路径与命令名一致：`reed::filter::run`、`reed::stats::run`、`reed::index::run`
（不写 `run_filter`，避免模块名重复）。

---

## 4. 类型与文件名的修订【规】

| 规则 | 说明 |
|---|---|
| 禁止版本与新旧后缀 | 标识符与文件名中不得出现 `V2`、`Legacy`、`Compat`、`New`、`Old`（`FilterSpecV2` → `FilterSpec`，`spec/legacy.rs`、`plan/legacy.rs`、`spec/compat.rs` 删除）。版本只存在于 crate 版本号与 git 历史中 |
| 统计请求类型 | `StatKind` → `StatSpec`（与 v1 §8.1 一致：Spec = 用户声明）；`StatKey` → `StatId` |
| 新增核心类型 | `GtCode`、`GenotypeScanner`、`FormatLayoutCache`、`SiteView`、`HeaderBundle`、`TagTable`/`TagId`、`BlockWriter`、`OrderedOutput`、`IndexKey`、`Tickets` |
| 文件名 | `vcf/sample_scan.rs` → `vcf/genotype_scan.rs`（它扫描的是基因型与 FORMAT 列，不只是"样本"）；`vcf/gt.rs` 删除 |
| 驱动目录 | `reed/src/filter/`、`reed/src/stats/`（每个命令一个目录，按 03 §8 拆文件） |
| 统计量模块 | `reed-stats/src/<统计量>.rs`，按方法命名：`pi.rs`、`dxy.rs`、`fst_wc84.rs`、`fst_hudson.rs`、`watterson_theta.rs`、`tajima_d.rs`、`pbs.rs`、`sfs_projection.rs`；模块文件中不出现软件名（软件名只出现在测试与对拍脚本中） |
| 方法枚举 | 以作者与年份命名：`FstMethod { WeirCockerham1984, Hudson1992 }`（CLI 取值 `wc84`、`hudson`）、`PbsMethod { Yi2010, Malaspinas2016 }`；日志中的引用字符串由枚举的 `citation()` 给出 |
| CLI crate | `reed-cli/src/args/<command>.rs`（clap 结构体）、`reed-cli/src/translate/<command>.rs`（参数 → Spec） |

---

## 5. 漏斗规则标识（修订 v1 §13.17）

点分小写，层级与 02-cli.md §6 的求值顺序表对应：

```
sample.list            sample.keep           sample.remove
region.targets         region.regions
site.types             site.alleles          site.snps_only        site.filter
site.include#1         site.exclude#2        （序号为该参数在命令行中的出现顺序）
genotype.min_dp        genotype.max_dp       genotype.min_gq       genotype.set_gts#1
sample.mind（掩码后）   after_mask.geno       after_mask.hwe        after_mask.maf        after_mask.mac        after_mask.include#3
```

报告 JSON 中每个 `-i/-e` 规则同时带上 `expr` 字段（表达式原文的规范化形式）。

---

## 6. 输出文件与列（补充 v1 §12、§13.16）

| 对象 | 规则 |
|---|---|
| 窗口统计 | `PREFIX_<stat>.txt`（`pi, dxy, fst, watterson_theta, tajima_d, pbs`），列布局沿用 pixy（02-cli.md §5.4）；1-based 闭区间；缺失值 `NA` |
| 逐位点/逐样本报告 | `PREFIX.afreq/.hardy/.het/.smiss/.vmiss`（plink2 的列布局）、`PREFIX.idepth/.ldepth.mean`（vcftools 的列布局） |
| SFS | `PREFIX_sfs/dadi/`、`PREFIX_sfs/fastsimcoal2/` |
| QC 汇总 | `PREFIX.stats`（bcftools stats 文本格式） |
| 运行日志 | `PREFIX.log`（仿 plink2） |
| VCF header 行 | `##reed_version`、`##reed_filterCommand=…; Date=…`（仿 bcftools）。`##reed_ascertainment` 废止 |

## 7. 变量命名补充（补充 v1 §13.7）

| 概念 | 规范名 |
|---|---|
| 一个样本的基因型编码 | `code`；切片 `codes` |
| 每群体的编码直方图 | `code_hist` |
| 字节预算许可 | `permit`；许可数量 `n_permit_bytes` |
| 规范化参数序列 | `canonical_args` |
| 表达式的语法树 / 编译后的谓词 | `expr_ast` / `predicate` |
