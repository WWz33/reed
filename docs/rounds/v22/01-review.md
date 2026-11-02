# 复现与证据（doc_v22）

> 基线：本机构建 `reed 0.4.0`（HEAD `02109e9` + 本轮改动）。
> 参考实现：`bcftools 1.24-2-g270751fc` / htslib 1.24-1-gb499e4c7。
> 所有 fixture 都在 `/tmp/bugcheck` 下用最小文本构造；真实数据为
> `/sri/home/sri2025201067/bsaseq/raw.bsa.vcf.gz`。

## 0. 先回答：INFO/FORMAT 解析在哪

原始文本的解析**不在一个文件里**，但每个格式只有一条解码路径：

| 关注点 | 位置 |
|---|---|
| 固定列 + INFO/FORMAT/样本区切片 | `crates/reed-io/src/vcf/line.rs`（`LineRef`，只有偏移，不解析内容） |
| 位点视图（CHROM/POS/ALT/allele） | `crates/reed-io/src/vcf/site.rs` |
| INFO 片段与索引 | `crates/reed-io/src/vcf/info.rs` |
| FORMAT 布局缓存（GT/整数/浮点列、本轮加 `IntSpec`） | `crates/reed-io/src/vcf/format_layout.rs` |
| GT 与 FORMAT 列的单遍解码 | `crates/reed-io/src/vcf/genotype_scan.rs` |
| BCF 原生解码/渲染 | `crates/reed-io/src/bcf/{record,render,columns,gt,dict,encode_native,encode_noodles}.rs` |
| 头部 TAG 声明表 | `crates/reed-core/src/tag_table.rs` |
| 头部驱动的规则校验（INFO/掩码） | `crates/reed-filter/src/plan/mod.rs` |
| 表达式解析（`INFO/…`、`FMT/…`） | `crates/reed-filter/src/expr/` |

所以：**VCF 文本的字段解析集中在 `reed-io/src/vcf/`**（`line` → `site`/`info` →
`format_layout` → `genotype_scan`），BCF 有平行的 `reed-io/src/bcf/` 解码器；
`reed-filter` 只做头部校验和规则编译，不重复解析记录文本。

## 1. sites-only VCF（8 列）

### 复现

```bash
cat > sites_only.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="d">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10	.	A	G	50	PASS	DP=12
chr1	20	.	C	T	60	PASS	DP=20
EOF
reed filter sites_only.vcf -Ov -o out.vcf
reed stats  sites_only.vcf --freq -o out
bgzip -c sites_only.vcf > sites_only.vcf.gz && reed index sites_only.vcf.gz
```

修前三条命令均为：`record line has fewer than 9 columns`（`index` 在压缩后
同样失败）；`bcftools view` 正常输出 8 列记录。

### 根因与修法

- `LineRef::parse` 收集 9 个 tab，缺少第 9 个就返回 `None`——把"记录截断"和
  "合法的 sites-only 形态"混为一谈。
- `LineRef` 记录实际 tab 数：7 = sites-only（8 列）、8 = 有 FORMAT 无样本、
  9 = FORMAT + 样本。`format_text()`/`samples_text()` 在无对应列时返回空；
  `n_sample_cols()` 为 0；`samples_start()` 指向行尾。
- `GenotypeScanner` 在 `n_samples == 0` 时跳过样本循环（此前会把空样本区
  当成一个空样本，随后报样本数列不一致）。
- `edit_vcf_line` 在无 FORMAT/样本区时只写 8 列；过滤输出（VCF 与 BCF）不再
  带多余的空列。
- 其余 `record line has fewer than 9 columns` 的报错与测试同步改为 8 列下限。

### 修后

```text
$ reed filter sites_only.vcf -Ov -o out.vcf
$ awk -F'\t' '!/^#/{print NF}' out.vcf
8
8
$ reed index sites_only.vcf.gz
reed index: 2 records -> sites_only.vcf.gz.csi in 0.0s
$ reed stats sites_only.vcf --freq --summary -o out
```

`out.afreq` 为 `ALT_FREQS=NA, OBS_CT=0`（无样本时的合法输出），`--summary`
正常写出。

## 2. 单倍体 GT

### 复现

```bash
cat > gt_shapes.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="g">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1
chr1	100	.	A	G	50	PASS	.	GT	1	0
chr1	200	.	A	G	50	PASS	.	GT	0/1	1/1
EOF
reed filter gt_shapes.vcf -Ov -o out.vcf
reed filter --help | grep -i ploidy   # 修前：无输出
```

修前：`policy rejected genotype at chr1:100: haploid genotype under
PloidyPolicy::RequireDiploid`。`bcftools view` 接受整份文件。

### 根因与修法

- `FilterSpec.policy.ploidy` 早已支持 `allow_haploid`（serde 名），
  `stats` 的 `stats_policies()` 也用 `HaploidAsHalfCall`；只有 `reed filter`
  的标志层没有入口，`FilterCli` 直接落到 `PolicyRules::default()`
  （`RequireDiploid`）。
- 新增 `--ploidy require-diploid|allow-haploid`（默认保持
  `require-diploid`），翻译到 `spec.policy.ploidy`。

### 修后

```text
$ reed filter gt_shapes.vcf --ploidy allow-haploid -Ov -o out.vcf
#CHROM ... FORMAT s0 s1
chr1 100 ... GT 1 0
chr1 200 ... GT 0/1 1/1
```

`reed stats gt_shapes.vcf --freq` 修前修后均通过（策略本就是半调用）。

## 3. 尾随空样本列

### 复现

```bash
printf '...header...\tFORMAT\ts0\ts1\nchr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\t\n' > trail_tab.vcf
reed filter trail_tab.vcf -Ov -o out.vcf
```

修前：`malformed record at chr1:10: more sample columns than the header
declares`。bcftools `view` 静默丢掉行尾空列并继续。

### 根因与修法

`GenotypeScanner::scan` 在 `sample_count >= n_samples` 时无条件报错。改为：
超出 header 样本数的**空**列忽略（继续扫描），非空额外列仍报 malformed——
不把真正的"多了一列基因型"当成可忽略噪声。

### 修后

命令通过；记录按原字节保留行尾 tab（未编辑时逐字节直通）。

## 4. 无 FMT/DP（触发项）：深度回退

### 复现

```bash
reed filter raw.bsa.vcf.gz --min-dp 10 -Ov -o out.vcf
# 修前：invalid filter spec: FORMAT field DP is not declared in the header
reed stats raw.bsa.vcf.gz --summary -o qc
# 修前：DP/PSC 节为空；模块注释却声称 "AD fallback: the maximum non-ref AD"
```

bcftools `stats -s -`：`calc_sample_depth` 在记录 FORMAT 无 DP 键时用该样本
AD 向量之和；`dp > 0` 才计入直方图与平均深度。

### 修法

- `IntSpec { key, fallback }`：主键读一个整数，回退键的逗号分隔向量求和；
  深度列由调用方显式构造 `IntSpec::depth() = DP（有则用），否则 AD 之和`
  （D22-1 收口后 `spec_for_key` 已删除）。
- `reed-io/vcf/genotype_scan.rs` 与 `reed-io/bcf/columns.rs` 的
  `decode_format_columns` 共用同一语义；filter 文本/BCF 两条路径一致。
- `--min-dp/--max-dp` 的 plan 校验：header 无 `FORMAT/DP` 时允许用
  Integer 的 `FORMAT/AD`。
- stats 的 DP/PSC、`--depth`、`--hist dp` 都读这条深度列；`dp == 0` 不再入
  直方图（bcftools `if (dp>0)`）。

### 修后（真实数据 chr1:1-200000）

```text
$ reed filter raw.bsa.vcf.gz -r chr1:1-200000 --min-dp 10 -Ov -o out.vcf
masked=1285 kept=3075        # 1111 个原调用 + 174 个原缺失
# 逐样本 AD 求和判据：called_low=1111, called_ok=3075, missing=174
$ reed stats raw.bsa.vcf.gz -r chr1:1-200000 --summary -o qc
$ diff <(grep -E '^DP|^PSC' qc.stats) <(bcftools stats -s - ... | grep -E '^DP|^PSC')
（无输出：逐行相同）
```

`subvcf.vcf.gz` 全量 `--summary` 的 DP/PSC 仍与 bcftools 逐行一致（有 DP 的
callset 行为不变）。

## 5. BCF `n_fmt == 0`（文本 FORMAT `.`）渲染

### 复现

```bash
printf '##fileformat=VCFv4.2\n##contig=<ID=Chr01,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
Chr01\t20\t.\tA\t.\t.\t.\t.\t.\t.\t.\n' > fmt_dot.vcf
bcftools view -Ob -o fmt_dot.bcf fmt_dot.vcf
reed filter fmt_dot.bcf -Ov -o out.vcf
```

`bcftools` 写入的 BCF 记录 `n_fmt == 0`、`n_sample == 2`（文本 FORMAT 为
`.`）。期望渲染（`bcftools view -Ov`）：

```text
Chr01	20	.	A	.	.	.	.	.	.	.
FORMAT=. + s1=. + s2=.        # 10 列
```

修前 reed 的 `render_sample_part` 在 `n_fmt == 0` 时先 push 了 FORMAT 前的
`\t`，随后循环 `0..=n_sample` 再写 `\t.`：FORMAT 列本身为空、样本列多出一
列（11 列），htslib 回读报 `FORMAT '' is not defined in the header`。

### 修法

按 htslib `vcf_format`：`n_fmt == 0` 时先写 FORMAT 列的 `.`，再按
`n_sample` 写 `\t.`（每列自带分隔符，不是"FORMAT 之后再补一个"）。
文件：`crates/reed-io/src/bcf/render.rs`（`render_sample_part`）。

### 修后

```text
$ reed filter fmt_dot.bcf -Ov -o out.vcf
$ diff <(grep -v '^#' out.vcf) <(bcftools view -Ov fmt_dot.bcf | grep -v '^#')
（无输出：数据行逐字节相同）
```

测试：`bcf::render::tests::fmt_no_format_keys_render_as_missing`（手工构造
`n_fmt == 0` 记录，断言输出 `Chr01\t20\t.\tA\t.\t.\t.\t.\t.\t.\t.\n`）。