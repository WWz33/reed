# 02 架构与规范修订（doc_v14）

> 本文只写**变更面**：分组的 help、拆分后的两个统计命令、`--hist`、
> `rename`、库 API 改名与共享判据。凡未提及的行为（filter/convert 全部
> 语义、三共享参数组、退出码、错误文本家规）一律沿用
> `docs/spec/cli.md` 现行条文。标记：【规】必须；【荐】推荐；【决】
> 待拍板（01 §5 C14 表）。

---

## §1 分组 help（seqkit 式）

### §1.1 形态与不变量【规】

- 命令**保持扁平**：单一 `Command` 枚举、`reed <verb> …` 调用路径一个
  字符不变（D14-8）。分组只存在于 `reed --help` 的呈现层。
- `reed --help` 的 Commands 段由默认单列表改为**按组渲染**：

```
Commands for Filtering:
  filter           hard filtering: selection, masking, QC thresholds, LD pruning

Commands for Statistics:
  stats            per-site and per-sample QC: freq, hardy, het, missing, depth, summary, histograms
  popg             windowed population genetics: pi, dxy, Fst, PBS, SFS, f-stats, LD decay

Commands for Transform:
  convert          export to phylogenetic matrices, 012, hapmap, MSA, SFS
  rename           edit identifiers (samples, contigs, variant IDs) and header lines

Commands for File Operations:
  index            build a CSI/TBI index

Commands for Miscellaneous:
  help             print this message or the help of the given subcommand
```

- `--version`/`-h` 行为不变；子命令 `--help` 完全不变（分组只在顶层）。

### §1.2 组的定义（目标态网格）【规】

| 组 | 本轮成员 | 预留位（后续轮，仅登记不实现） |
|---|---|---|
| Filtering | `filter` | — |
| Statistics | `stats`, `popg` | `pca`, `relatedness`（v15 候选） |
| Transform | `convert`, `rename` | `norm`, `query`（v15+ 候选） |
| File Operations | `index` | `concat`, `isec`, `merge`（v15+ 候选） |
| Miscellaneous | `help` | `genautocomplete`（clap_complete，可随时加） |

网格意义：新动词落组即回答"它是什么形状的工具"；**Statistics 组的准入
判据 = 报告行的单位**（位点/样本 → stats；群体×窗口/全基因组 → popg；
矩阵输出 → pca/relatedness 独立动词），Transform 组的准入判据 =
**逐记录编辑/转换但不清除**（不改存活集，改表示）。

### §1.3 实现机制【规】

- `reed-cli/src/cli/mod.rs` 新增常量：

```rust
/// The help groups (doc_v14 02 §1.2): the single source of the grouping;
/// a test asserts every Command variant appears in exactly one group.
const GROUPS: &[(&str, &[&str])] = &[
    ("Filtering", &["filter"]),
    ("Statistics", &["stats", "popg"]),
    ("Transform", &["convert", "rename"]),
    ("File Operations", &["index"]),
    ("Miscellaneous", &["help"]),
];
```

- 顶层 `Cli` 设 `help_template`：用模板标记位 + 一个小渲染函数遍历
  `Command::command().get_subcommands()`，按 `GROUPS` 名单分桶，每命令
  取 clap 自身的 `get_about()`（**不手抄帮助文本**——防 v12 D-3 类
  帮助文本过期漂移）；桶外变体守卫测试拦截。
- `popg` 的 alias `popgen` 不入组名单（alias 不是变体）；`reheader`
  作为 `rename` 的 visible alias 同理。

### §1.4 守卫【规】

1. 单元测试：`GROUPS` 展开去重后的集合 == `Command` 变体名集合
   **∪ {help}**（clap 的 `help` 子命令是自动生成的、不占 `Command`
   变体，但出现在 `--help` 的组渲染里；多一个、少一个、重复都挂）。
2. `cli_contract.rs` 快照：`reed --help` 含五个组标题、每组行数正确。
3. README（双语）命令表与 `GROUPS` 同源（生成或测试锁定一致）。

---

## §2 命令面拆分：`stats` 与 `popg`

### §2.1 命令定义【规】

| 命令 | 驱动门面 | 行的单位 | `--pop` |
|---|---|---|---|
| `reed stats <INPUT> [选择组] [QC flags] [-o PREFIX]` | `run_stats`（新含义，原 `run_reports`） | 位点 / 样本 / 全文件汇总 | **可选**（分层 QC） |
| `reed popg <INPUT> [选择组] [--pop FILE] [清单/窗口参数] [统计参数…] [-o PREFIX]` | `run_popgen`（原 `run_stats`） | 群体（对/三元组）× 窗口 / 全基因组 | 窗口模式必需（现状不变） |

flag 清单与搬迁表见 01 §3.1；`popg` 主名 + `popgen` visible alias
（C14-1 已拍板，2026-10-17）。

### §2.2 stats（拆分后）的参数面【规】

- **继承不动**：`--freq --hardy --het --missing --depth --summary`
  （拼写、输出文件、列布局逐字节不变，D14-1）；
  Region/Sample/Site 三共享组照旧 flatten。
- **新增**：
  - `--hist LIST`（§3）；
  - `--pop FILE`、`--force-samples`（**保留降级**：现 stats 已有的两个
    flag，拆分后 stats 侧语义为"分层 QC 可选"——`--missing`/`--freq`/
    `--singletons` 增加按群体的行/列；不给 `--pop` 时行为与 0.3.x 完全
    相同）；
  - `--singletons`（C14-3 默认进轮：`PREFIX.singletons`，vcftools 列
    布局 CHR POS A R A1 A2 TYPE POP，POP 列接 `--pop`，无 `--pop` 时
    为 `.`）。`--geno-counts` 留 backlog（C14-3）。
- **无统计 flag 的裸调用**：usage 错，提示可选 flag 列表（C14-2 默认；
  bcftools stats 的"无选项全做"不沿用）。

### §2.3 popg（拆分后）的参数面【规】

- 搬迁 flag 的拼写、默认值、互斥规则、输出文件与列布局全部不变
  （01 §3.1 搬迁表全量列出）；`--indep-pairwise` 在 Site 组照旧存在，
  报错文本从 "not supported yet … feed its output to reed stats" 改为
  指向 popg 自身后续支持或 filter（现状语义不变，只改引导文本）。
- `--sfs` 等价契约：`popg --sfs` ≡ `convert --sfs` 逐字节（D14-3，
  原 `stats --sfs` 的门改挂 popg；`translate::build_sfs_request`
  单一实现不动）。

### §2.4 遗留重定向【规】

旧 `stats` 上的全部搬迁 flag（01 §3.1 后两栏）成为 stats 的 hidden
removed-flags：解析期报错（exit 2）：

```
reed: --pi moved to reed popg: run 'reed popg <input> --pop pops.txt --pi …'
reed: --window moved to reed popg: run 'reed popg <input> --window 100k …'
```

错误文本给出**可执行的完整改写**（把原命令行里的该 flag 连同已知参数
原样搬运，模式同 `--max-missing → --geno` 的精确换算提示）；无静默
别名、无转发执行（D14-2）。

### §2.5 兼容与版本【规】

- 版本 0.4.0；`docs/spec/cli.md` 拆为 "stats（QC）" 与 "popg" 两节，
  重定向表进 §11 差异登记；CHANGELOG 逐 flag 列搬迁；README 双语命令
  面重写。
- 库 API 同轮改名（§6）——库无兼容承诺（cli.md §0.9），但
  `to_args()` 往返测试须同步改造：`StatsRequest::to_args` 渲染的
  前缀变为 `popg`。

---

## §3 `--hist`（stats）

### §3.1 请求语法【规】

```
--hist LIST        # 逗号分隔；元素 = 内建 metric 或 INFO/TAG 或 FMT/TAG
```

- 内建 metric：`qual`（QUAL，0.1 bin，SNV/indel 合并计数，与 `--summary` 的 QUAL 节同源）、`dp`
  （FMT/DP，逐基因型）、`dp-site`（INFO/DP，逐位点）、`gq`（FMT/GQ，
  逐基因型，**新增解码**）、`idl`（indel 长度，ins/del 分列）、`af`
  （等位频率，bcftools `init_iaf` 的 bin 语义）。
- 通用标签：`INFO/<TAG>`、`FMT/<TAG>` 任意数值标签（Integer/Float；
  Flag/String/Character 报 usage 错）；`FMT/` 聚合全体样本观测。
- 重复元素去重；空元素忽略；`--hist` 与 `--summary` 可同跑（独立文件，
  一次扫描双消费——v13 多目标共享行解析的纪律延续）。
- **可路由**（荐）：`--hist` 不要求其他 QC flag，可单独
  `reed stats in.vcf.gz --hist qual,gq -o out`（裸 stats 报错规则的
  例外：`--hist` 即统计请求）。

### §3.2 输出布局【规】

`PREFIX.<m>.hist`（`<m>` = metric 名或 `INFO_TAG`/`FMT_TAG` 形态）：

```
bin_lower	bin_upper	count	percent
0	1	12345	0.456789
```

- 整数 metric（dp/dp-site/gq/idl 及整数标签）：`bin_lower == bin_upper`
  为精确值行；文件尾附三行汇总 `#mean`/`#median`/`#p95`（`#` 前缀注释
  行，`IntHistogram::quantile` 精确分位，B14-2）。
- 浮点 metric（qual/af 及浮点标签）：固定边界（C14-5 默认 bcftools
  同款：qual 0.1 步长到 100 + overflow 桶；af 按 init_iaf bin）。
- 溢出桶命名 `>max`；缺失观测（`.）`不入任何桶但计入 `#missing` 行。

### §3.3 计划与解码【规】

- `HistPlan` 在翻译期完成：解析 LIST、按 header 校验每个标签存在且
  为数值类型（cli.md §0.7"参数在第一条记录前全部校验"沿用），内建
  metric 检查对应字段存在（无 FMT/GQ 声明 → usage 错，同 `--depth`
  的既有行为）。
- worker 内：`--summary` 已解码的 QUAL/INFO-DP/FMT-DP/IDD 直接复用
  （零重复解析）；GQ 与通用标签走 GenotypeScanner 的字段计划扩展
  （`--depth` 独立解码 DP 的先例）。

### §3.4 累加器与确定性【规】

新累加器（GQ、通用 `IntHistogram`/`FloatHist`）满足 merge 结合/交换律
（D14-7）；`IntHistogram::merge` 既有实现沿用；e2e 断言 threads
∈ {1,4} 输出 sha256 相同。

### §3.5 B14 登记

`--hist` 的差异（通用标签、独立 TSV、精确分位、浮点固定边界）全部
登记 `docs/spec/cli.md` 差异表（B14-1/2/3、C14-5），对拍门只比
bcftools 有的部分（§7 of 01）。

---

## §4 `rename`（新子命令）

### §4.1 定位与硬边界【规】

`rename` 是**标识符编辑器**：只碰 header 行、CHROM 列、ID 列、样本名。
**禁止**：解码基因型（结构性守卫：`crates/reed/src/rename/` 不得引用
GenotypeScanner/GtCode）；修改 POS/REF/ALT/FILTER/INFO/FORMAT（D14-4）；
清除断言（"选择"语义归 filter）。Transform 组成员（§1.2 判据）。

```
reed rename <INPUT> [编辑平面 flags…] [-o FILE] [-O v|z|b|u[0-9]] [-W[=csi|tbi]]
```

`reheader` 为 visible alias（bcftools 用户可发现性）。

### §4.2 四个平面与 flag 全集【规】

| 平面 | flag | 语义（来源） |
|---|---|---|
| header | `-h, --header FILE` | 整体替换 header（bcftools reheader `-h`；新头须含全部样本名） |
| header | `-f, --fai FILE` | 按 .fai 补/改 `##contig` 行（bcftools reheader `-f`；reed 窗口统计的前置刚需） |
| header | `--add-meta LINE`（可重复） | 追加 `##key=value` 行 |
| header | `--remove-meta KEY`（可重复） | 删除全部 `##KEY=` 行（剥 `##reed_command`、caller 来源） |
| 染色体 | `--rename-chrs FILE` | `old<TAB>new` 映射；**同时改 `##contig` 行与记录 CHROM 列**（D14-5；bcftools annotate 同名） |
| 染色体 | `--chr-style plain\|ucsc\|mt-is-m` | 约定转换三态：加/去 `chr` 前缀、`MT↔M`（B14-4；plink2 `--output-chr` 的有限子集） |
| ID | `--set-id PATTERN` | 模板写全部 ID：`%CHROM %POS %REF %ALT`（`%ALT` 多等位逗号连接，B14-5；bcftools annotate `-I` 拼写） |
| ID | `--set-missing-var-ids PATTERN` | 只填 `.` 的 ID（plink2 同名拼写） |
| ID | `--update-name FILE` | `old<TAB>new` 映射改 ID（plink2 同名） |
| ID | `--idspace-to SEP` | 把 ID 中现分隔符（`:`) 换成 SEP（plink2 同名；GCTA/treemix 的冒号过敏解药） |
| ID | `--rm-dup-ids` | 重复 ID 加 `_1/_2…` 后缀保全部（B14-6；seqkit rename 语义） |
| ID | `--clear-ids` | 全部置 `.` |
| 样本 | `-s, --samples PAIRS` | `"old new,old2 new2"`（bcftools reheader `-s`） |
| 样本 | `-N, --samples-file FILE` | 改名对文件（bcftools reheader `-N`） |
| 样本 | `--replace-sample-name PATTERN` | 字面改名：`pre*` 去前缀、`*suf` 去后缀、其余整体删除；带正则元字符（`^$[]()|+?{}\`）报用法错（不引入正则依赖） |

输出族 `-o/-O/-W` 与运行族（`--threads/--memory/--verbosity/--no-version/
--no-command`）复用 filter 的拼写与语义。**不收** Region/Sample/Site
选择组（编辑无选择语义；先筛再改变管道 `reed filter … | reed rename …`）。
`--order-samples` 不进首发轮（C14-4）。

### §4.3 求值顺序与校验【规】

```
1. -h 整体替换（含样本表）      5. ID：--update-name → --set-id/
2. 染色体：--rename-chrs 与        --set-missing-var-ids → --idspace-to
   --chr-style（先映射后约定，      → --rm-dup-ids → --clear-ids
   二者可叠加；同改头与记录）     6. 样本：-s/-N → --replace-sample-name
3. -f/--fai 覆盖 contig 长度      7. --add-meta/--remove-meta
4. 全部校验（下表）               8. 写出
```

校验（全部在第一条记录前，exit 2）：

| 规则 | 行为 |
|---|---|
| 改名后样本名重复 / 新名与未改名样本冲突 | 报错 |
| `-s/-N/--update-name/--rename-chrs` 映射里的名字在输入中不存在 | 报错（`--force-samples` 风格降级为警告不引入——rename 是编辑器，拼错名字静默跳过比 filter 危险；0.4.0 先严后松） |
| `--set-id` 与 `--set-missing-var-ids` 同给 | 报错（互斥） |
| `--clear-ids` 与其他 ID 平面 flag 同给 | 报错 |
| `--idspace-to` 与 `--set-id` 同给 | 报错（先模板再换分隔符语义含混；用户应直接在模板里写 SEP） |
| 无任何编辑 flag | usage 错（`--help` 引导） |
| 输出为 BCF 且染色体改名 | 允许（字典重写，B14-7），但头块全量重编码路径须 e2e 锁定 |

### §4.4 写出与索引【规】

- VCF/VCF.gz 输入 → 字节编辑路径（`edit_vcf_line`：CHROM/ID 列替换
  + 头行编辑；记录其余字节直通）。BCF 输入 → 头块重写 + 记录重编码
  （染色体重名时 contig 字典重映射）。纯样本改名（BCF 与 VCF 均是）
  走快路径：只重写头块、记录体直通（B14-7，bcftools reheader 同款
  卖点）。
- 旧 `.csi/.tbi` 在写出时**不自动失效处理**：检测到输入旁存在索引且
  染色体/位点 ID 改动影响其语义时打印 warning 提示 `reed index`
  （C14-6 默认）；`-W` 照 filter 语义为新输出建索引。

### §4.5 provenance【规】

`##reed_version`/`##reed_command`/`##reed_rename`（编辑摘要：各平面
计数）进输出头；`--no-version/--no-command` 复用 filter 语义。

---

## §5 共享参数组判据（成文）【规】

**选择参数共享，动作参数各归**：

| 组 | 内容 | 归属 |
|---|---|---|
| `RegionArgs`/`SampleArgs`/`SiteFilterArgs` | 回答"统计/输出看见哪些位点与样本" | filter、stats、popg（现状照旧）；convert 只收 Sample 组（v13 现状） |
| 各命令动作 flag | 回答"做什么" | 各命令自持 |

`rename` 不收任何选择组（§4.2）；后续 `norm`（Transform）同样不收，
`concat/isec/merge`（File Operations）按各自语义定义（后续轮规范）。

---

## §6 库 API 改名【规】

| 旧 | 新 | 说明 |
|---|---|---|
| `run_stats` / `run_stats_opened` | `run_popgen` / `run_popgen_opened` | 驱动 popg；`StatsRequest` → `PopgenRequest`（`StatKind` 等内核类型名不动——它们是估计器词汇不是命令词汇） |
| `run_reports` / `run_reports_opened` | `run_stats` / `run_stats_opened` | 驱动 stats；`ReportRequest` → `StatsReportRequest`，`ReportsOutcome` → `StatsOutcome` |
| （新增） | `run_rename` / `RenameRequest` / `RenameOutcome` | 门面根导出；`to_args()` 可失败（v11 A11-3 纪律） |
| （新增） | `HistRequest`（`ReportRequest` 的扩展字段或并立类型，实现期定） | 累加器在 `reed-stats::summary` 扩展 |

改名同轮完成（库无兼容承诺）；`to_args` 往返测试全部更新；README 的
Library use 一节示例改写并进编译门（v13 惯例）。

---

## §7 错误、退出码与日志【规】

- 退出码家规不变：usage 2 / 运行期 1；`reed: …` / `reed: warning: …`。
- 新错误文本族：
  - 重定向：`--pi moved to reed popg: run '…'`（D14-2）；
  - `--hist`：未知 metric、非数值标签、标签未声明；
  - `rename`：§4.3 校验表。
- `popg` 的窗口模式写 `PREFIX.log`（方法与文献、漏斗、计数器）；QC 报告
  （`stats`）不写日志——0.3.x 的 QC 路径即如此，拆分未改（实测参照
  二进制只写窗口模式的日志）。

## §8 验收门（实施期逐条）

1. 双 profile 全绿 + clippy `-D warnings` + `check-doc-drift.sh`
   （含本轮守卫 1–14）。
2. 拆分等价门：`run_split_parity.sh` 全绿（01 §7）。
3. 确定性：`--hist` 与 `rename` 全 flag 组合 threads {1,4} sha256 相同。
4. 既有面不回归：v13 的 `run_convert_parity.sh`、`run_easysfs.sh`、
   `run_pixy*.sh`、`run_vcftools*.sh`、`run_plink2_ld.sh` 全绿（部分
   脚本命令词替换 stats→popg，diff 语义不变）。
5. 性能门：01 §8 表的退化上限。
6. README（双语）/cli.md/CHANGELOG/架构 spec §13.14（子命令清单）同轮。

## §9 文档回改清单

| 文档 | 改动 |
|---|---|
| `docs/spec/cli.md` | §1 命令表加 popg/rename；§7 拆为 stats/popg 两节；§11 登记 B14 与重定向 |
| `README.md` / `README.zh-CN.md` | 命令面重写（分组表 + stats/popg/rename 用例 + `--hist` 表） |
| `CHANGELOG.md` | 0.4.0：逐 flag 搬迁清单 + 新 flag + 库 API 改名 |
| `docs/spec/architecture.md` §13.14 | 子命令清单更新（filter/stats/popg/index/convert/rename + 组） |
| `docs/spec/statistics-map.md` | 落库计划补 `--hist`（§10.2 落地）与 `--singletons` |
