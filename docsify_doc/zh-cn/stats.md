# stats

## 用法

```
reed stats [OPTIONS] --freq --hardy ... -o PREFIX <INPUT>
```

位点级与样本级 QC 报告。单遍读取，`-o` 当作输出前缀，报告写文件。不给 `-o`、一个报告都没选、`--pop` 没配能分层的报告，都是用法错。

## 参数

### 报告

```
    --freq                  每位点等位频率 → PREFIX.afreq
    --hardy                 HWE 精确检验与观察/期望杂合度 → PREFIX.hardy
    --het                   每样本近交系数 F → PREFIX.het
    --missing               每样本与每位点缺失率 → PREFIX.smiss、PREFIX.vmiss
    --depth                 基于 FMT/DP 的平均深度 → PREFIX.idepth、PREFIX.ldepth.mean
    --summary               QC 汇总 → PREFIX.stats
    --hist LIST             直方图，一个指标一个文件 → PREFIX.<指标>.hist。内置 qual、
                            dp（FMT/DP）、dp-site（INFO/DP）、gq、idl（有符号 indel
                            长度）、af；也收 INFO/<TAG>、FMT/<TAG>
    --singletons            逐等位 singleton 与私有 doubleton → PREFIX.singletons
    --sites-file FILE       CHROM<TAB>POS（1-based），只报告列出的位点
-o, --output PREFIX         输出前缀；报告写文件，不写 stdout。别名 --out
```

### 群体分层

```
    --pop FILE              样本<TAB>群体 两列，无表头。给 --freq/--missing 加分层文件，
                            给 --singletons 加 POP 列。别名 --populations
    --force-samples         群体文件里的未知样本只警告并跳过；默认报错
```

### 输入与区域

```
<INPUT>                     VCF/VCF.gz/BCF；- 读 stdin
-r, --regions LIST          只要这些区域，逗号分隔；需要输入旁边的索引
-R, --regions-file FILE     区域清单文件，格式同 filter；需要索引
-t, --targets LIST          同样的区域写法，流式，不需要索引
-T, --targets-file FILE     区域清单文件；行首 ^ 表示排除。不需要索引
```

### 全局参数

```
    --threads INT           工作线程数，默认全部可用核
    --memory MB             在途内存预算，MiB，默认物理内存一半
    --tmp-dir DIR           为兼容接受；单遍报告不溢写
    --verbosity INT         0 静默、1 摘要（默认）、2 阶段耗时、3 调试
```

### 会被拒绝的参数

```
-s, --samples / -S, --samples-file / --keep / --remove
                            报告模式不收，报错要求先跑 reed filter
    --mind                  单独的错误：样本规则需要先跑一遍全量数据，
                            提示先 reed filter --mind X
    --indep-pairwise        单独的错误：剪枝需要第二遍扫描，
                            提示先 reed filter --indep-pairwise
-v, --types / -V, --exclude-types / --snps-only
-m, --min-alleles / -M, --max-alleles / -f, --apply-filters
-i, --include / -e, --exclude / -g, --genotype
    --min-dp / --max-dp / --min-gq / --set-gts
                            （掩码三件套的别名 --vcf-min-dp / --vcf-max-dp /
                            --vcf-min-gq）
-q, --min-af / -Q, --max-af / --geno / --maf / --max-maf / --mac /
--max-mac / --hwe
                            位点规则与掩码同样不收，报同一条 report mode 错误
```

report mode 错误信息：

```
reed: report mode (--freq/--hardy/--het/--missing/--depth/--summary) does not apply
-s/-S/--keep/--remove or site/threshold filters yet; run `reed filter` first and feed its output to reed stats
```

`--indep-pairwise` 的错误信息单独给出剪枝建议。要统计过滤后的子集，先 `reed filter` 出文件，再喂给 `stats`。

## 输出格式

`qc/raw.singletons`（文件自带表头，一行一个等位；S 表示只有 1 份拷贝，D 表示 2 份拷贝且携带者纯合）：

```
CHROM	POS	SINGLETON/DOUBLETON	ALLELE	INDV
chr1	501	D	G	s3
```

`qc/raw.vmiss`：

```
#CHROM	ID	MISSING_CT	OBS_CT	F_MISS
chr1	rs2	1	10	0.1
chr1	rs4	1	10	0.1
chr1	rs8	1	10	0.1
```

`qc/raw.dp.hist`（分桶就是 DP 的取值本身，只写非空桶）：

```
bin_lower	bin_upper	count	percent
2	2	126	2.520000
3	3	47	0.940000
4	4	58	1.160000
5	5	70	1.400000
```

多等位位点不拆分：`.afreq` 的 ALT_FREQS 是全部 ALT 等位合并的频率（`CollapseAlt` 语义），全家族一致。

## HWE 精确检验

`--hardy` 用的检验与 `filter --hwe` 同一实现：给定已调用二倍体数 `n` 与次要等位拷贝数 `mac`，杂合子数 `nh` 的条件分布为 `P(nh) ∝ n!/(hm!·nh!·hM!)·2^nh`（`hm = (mac−nh)/2`），在对数域按 `P(nh+2)/P(nh) = 4·hm·hM/((nh+1)(nh+2))` 递推，归一化后双尾 p 值取所有概率不超过观测值的总和。单态位点 p = 1。mid-p 变体（`p_mid = p_exact − 0.5·p_obs`）在库中已实现，尚未暴露为命令行选项。

## 直方图

内置指标：

| 指标 | 取值 |
|---|---|
| `qual` | QUAL |
| `dp` | `FMT/DP`（基因型级） |
| `dp-site` | `INFO/DP`（位点级） |
| `gq` | `FMT/GQ` |
| `idl` | 有符号 indel 长度（缺失为负、插入为正） |
| `af` | 等位频率 |

任何 header 声明为 Integer 或 Float 的 `INFO/<TAG>`、`FMT/<TAG>` 也能用。分桶边界就是该量本身的单位：QUAL 按质量值、`af` 按频率。浮点指标的边界是浮点；整数指标文件末尾附三行精确分位数：

```
#mean	31.012800
#median	30
#p95	58
```

有缺失观测时附 `#missing` 行；超出累加器量程的整数观测写 `>max` 行且不进分位数。标签必须在 header 声明且为数值类型，否则用法错。

## 群体分层

`--pop FILE`（`样本<TAB>群体`，无表头）给三个报告加分层数据：

- `--freq` 另写 `PREFIX.afreq.pop`（列 `#CHROM ID POP REF ALT ALT_FREQS OBS_CT`），每群体一行；
- `--missing` 另写 `PREFIX.vmiss.pop`（列 `#CHROM ID POP MISSING_CT OBS_CT F_MISS`）；
- `--singletons` 加第六列 POP（携带者的群体；不在群体文件里的样本记 `.`，计数仍按全体样本）。

```
qc.afreq.pop
#CHROM	ID	POP	REF	ALT	ALT_FREQS	OBS_CT
chr1	rs2	popA	A	C	0.5	10
chr1	rs2	popB	A	C	0.5	8

qc.vmiss.pop
#CHROM	ID	POP	MISSING_CT	OBS_CT	F_MISS
chr1	rs2	popA	0	5	0
chr1	rs2	popB	1	5	0.2
```

分层值等于把该群体的样本单独取出后的对应报告值。群体文件里的未知样本默认报错，`--force-samples` 改为跳过并警告。

## 例子

```bash
mkdir -p qc

# 全套 QC
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc/raw

# 直方图选阈值
reed stats raw.vcf.gz --hist qual,dp-site,gq,INFO/DP -o qc/raw

# 逐等位 singleton
reed stats extra.vcf.gz --singletons -o qc/raw

# 分层等位频率与缺失率
reed stats clean.vcf.gz --pop pops.txt --freq --missing -o qc/by_pop

# 只看一段区域的深度
reed stats raw.vcf.gz -r chr1:1-1000000 --depth -o qc/chr1

# 只统计清单上的位点
reed stats raw.vcf.gz --sites-file candidate_sites.txt --freq -o qc/cand

# 先过滤再统计
reed filter raw.vcf.gz --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 -o clean.vcf.gz -Oz
reed stats clean.vcf.gz --freq --missing -o qc/clean
```
