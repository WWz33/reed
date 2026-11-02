# reed doc_v22：真实 callset 的 I/O 兼容与深度回退（2026-11-03）

> 触发：真实 BSA-seq callset `bsaseq/raw.bsa.vcf.gz`（4 样本、`GT:PL:AD`、
> 无逐样本 `FMT/DP`）上的失败，以及仓库外对三类"bcftools 能读、reed 报错"
> 输入形态的复现清单。
> 方法：每条先用最小 fixture 在本机构建 `reed` 上复现，再用
> `bcftools 1.24-2`（htslib 1.24-1）确认参考行为；修复后补失败测试，
> 并在真实 `raw.bsa.vcf.gz` 与真实 `subvcf.vcf.gz` 上回归。
>
> 本目录：README（结论与决策）、[01-review.md](01-review.md)（逐条复现与证据）、
> [STATUS.md](STATUS.md)（修复表与测试名）、[06-roadmap.md](06-roadmap.md)
> （架构收口决策 D22-1 与 R22 收口记录）。

## 结论速览

| # | 报告的问题 | 判定 | 说明 |
|---|---|---|---|
| 1 | sites-only VCF（8 列）完全不能读 | **属实（高）** | `LineRef::parse` 要求 9 个 tab；`filter`/`index`/`stats` 全部报 `record line has fewer than 9 columns`。bcftools `view -G` 正常输出该形态 |
| 2 | 单倍体 GT 被默认策略拒绝且无 CLI 开关 | **属实（中）** | `filter` 的 `PloidyPolicy::RequireDiploid` 写死在标志层；spec 已有 `allow_haploid`，但没有标志可设。`stats` 本来就用 `HaploidAsHalfCall` |
| 3 | 尾随空样本列导致整文件报错 | **属实** | `more sample columns than the header declares`；bcftools 静默忽略行尾空列 |
| 4 | （触发项）无 `FMT/DP` 时 DP 掩码与深度报告恒空 | **属实（高）** | `--min-dp` 直接报 `FORMAT field DP is not declared`；`stats --summary` 的 DP/PSC 与 `--depth`/`--hist dp` 无数据。`report.rs` 的 "AD fallback" 注释从未实现且口径写错 |
| 5 | BCF `n_fmt == 0`（文本 FORMAT `.`）渲染列错位 | **属实** | 渲染成"空 FORMAT 列 + 多一个空样本列"，htslib 回读报 `FORMAT '' is not defined in the header`；现按 `vcf_format` 输出 FORMAT `.` + 每样本 `.`。真实 `bcftools view -Ob` 生成的 `FORMAT=.` BCF 上 `reed filter -Ov` 与 `bcftools view` 数据行逐字节相同 |

## 关键决策

| 决策 | 结论 |
|---|---|
| 深度的定义 | 记录内：FORMAT 有 `DP` 键 → 该样本的 `DP`（值为 `.` 不回退）；无 `DP` 键 → 该样本 `AD` 向量之和（bcftools `calc_sample_depth`）。AD 向量全量求和（不按 `Number=R` 截断），`.` 跳过，向量尾截断，全缺失记 missing |
| 回退放在哪一层 | **IO 层**：`IntSpec { key, fallback }`；文本扫描器与 BCF 原生解码共用同一 spec。深度列由消费方显式构造 `IntSpec::depth()`（DP，回退 AD 求和），IO 不再按字段名推断（`spec_for_key` 已删除，doc_v22 R22-1）。filter/stats 不再各自解析文本 |
| `--min-dp/--max-dp` 的兼容性 | 只在"该记录 FORMAT 无 DP 键"时回退；有 DP 的 callset 与 plink2/vcftools/bcftools 行为逐位不变。`docs/spec/cli.md` 已标注 AD 回退为 reed 扩展 |
| 深度 0 的处理 | bcftools 只在 `dp > 0` 时计入 DP 直方图/PSC 平均深度（`vcfstats.c:1119`）；reed 原先也计入 0，本轮一并修正 |
| `--ploidy` 默认值 | 保持 `require-diploid`（不改变既有 callset 的语义），新增 `--ploidy allow-haploid` 显式放宽；`stats` 的既有 `HaploidAsHalfCall` 不动 |
| sites-only | `LineRef` 接受 ≥8 列；`FORMAT`/样本区为空。filter/index/stats/rename/convert 的 9 列下限改为 8。sites-only 输出保持 8 列（`-Ov/-Oz/-Ob` 均验证） |
| 尾随空列 | 仅忽略**超过 header 样本数**且为空的列；非空额外列仍报 malformed（不掩盖真实损坏） |
| BCF `n_fmt == 0` 的渲染 | 按 htslib `vcf_format`：FORMAT 列 `.`，每样本一个 `.`（输出列数 = 8 + 1 + n_sample）；`n_fmt == 0 && n_sample == 0`（sites-only BCF）不写 FORMAT/样本区，保持 8 列（V22-R1）。真实 `bcftools view -Ob` 的 `FORMAT=.` 文件对拍通过 |
| `--depth` 的 vcftools 兼容 | 深度列同 `--min-dp`（DP，记录无 DP 键时 AD 之和）；有 DP 的输入与 vcftools 输出不变，文档标注扩展 |

## 改动文件

| 层 | 文件 | 内容 |
|---|---|---|
| io | `crates/reed-io/src/vcf/format_layout.rs` | `IntSpec` / `IntLayout`、`FormatLayoutCache::new_with_specs`；回退策略由调用方显式声明 |
| io | `crates/reed-io/src/vcf/genotype_scan.rs` | `new_with_specs`、`parse_int_sum`、8 列与尾随空列处理 |
| io | `crates/reed-io/src/vcf/line.rs` | `LineRef` 支持 8 列；`format_text`/`samples_text`/`n_sample_cols`/`samples_start` 的空区语义 |
| io | `crates/reed-io/src/vcf/line_edit.rs` | sites-only 记录输出保持 8 列 |
| io | `crates/reed-io/src/bcf/columns.rs` | `decode_format_columns` 按显式 `IntSpec` 解码：主键标量、回退键求和 |
| io | `crates/reed-io/src/bcf/render.rs` | `n_fmt == 0` 按 htslib `vcf_format` 渲染 FORMAT `.` + 每样本 `.`；sites-only 保持 8 列；FORMAT 描述符遍历与 `columns.rs` 共用 `record::format_fields` |
| filter | `crates/reed-filter/src/plan/mod.rs` | 仅声明 `FORMAT/AD` 时 DP 掩码规则合法 |
| filter | `crates/reed/src/filter/record_pipeline.rs` | 文本扫描器按编译计划的显式 `IntSpec` 建列 |
| filter | `crates/reed-cli/src/cli/{mod.rs,filter_flags.rs}` | `--ploidy`、`--min-dp/--max-dp` 帮助 |
| filter | `crates/reed-cli/src/translate/filter.rs` | `--ploidy` → `policy.ploidy` |
| stats | `crates/reed/src/stats/writers/{summary.rs,reports.rs,histogram.rs}` | 深度列统一、`dp > 0` |
| stats | `crates/reed-stats/src/summary/report.rs` | 文档口径修正 |
| 其余 | `crates/reed/src/{index,genotype,rename,convert,filtered}.rs` 等 | 9 列错误信息与 8 列下限 |
| 文档 | `docs/spec/cli.md`、`CHANGELOG.md` | 深度/ploidy/sites-only 语义 |

## 验收

- `ci/test-parallel.sh`（全量工作区）：**860 passed / 0 failed**（56 个测试二进制）。
- `cargo clippy --workspace --all-targets`：0 告警。
- 真实数据：`raw.bsa.vcf.gz` chr1:1-200000 上 `--min-dp 10` 掩码 1111 个原调用
  + 174 个原缺失（与逐样本 AD 求和判据一致）；`--summary` 的 DP/PSC 与
  `bcftools stats -s -` 逐行相同。`subvcf.vcf.gz` 全量 `--summary` 仍逐行一致。
- BCF `FORMAT=.`：`bcftools view -Ob` 生成的真实文件上，`reed filter -Ov`
  数据行与 `bcftools view` 逐字节相同。