---
title: filter
---

## 用法

```
reed filter [OPTIONS] <INPUT>
```

硬过滤：位点与样本选择、基因型掩码、QC 阈值、LD 剪枝。默认只保留双等位位点（`-M 2`），输出写 stdout，`-o` 改写到文件。只输出 VCF，不做统计汇总。

## 参数

### 输入与区域

```
<INPUT>                     VCF/VCF.gz/BCF；- 读 stdin。多遍运行（--mind、LD 剪枝）
                            先溢写到 --tmp-dir
-r, --regions LIST          只要这些区域，逗号分隔：chr1、chr1:100、chr1:100-200、
                            chr1:100-、chr1:-200。需要输入旁有 .csi/.tbi 索引
-R, --regions-file FILE     区域清单文件，一行一个。三列及以上按 BED 读（0-based 半开），
                            两列按 CHROM POS（1-based），一列按区域串。需要索引
-t, --targets LIST          同样的区域写法，流式处理，不需要索引
-T, --targets-file FILE     区域清单文件，格式同 -R；行首 ^ 表示排除。不需要索引
```

`#` 注释行与空行跳过。

### 样本选择

```
-s, --samples [^]LIST       只要这些样本，逗号分隔；^ 前缀取补（^A,B 是除 A、B 外全要）
-S, --samples-file FILE     样本名单放文件，一行一个
    --keep FILE             同 -S，另一种拼写
    --remove FILE           要丢掉的样本名单，一行一个
    --force-samples         名单里出现 header 没有的样本名时只警告，不报错
    --mind FLOAT            缺失率高于 FLOAT 的样本丢掉。按缺失等位 / 2·位点计，
                            统计范围是掩码之后、通过掩码前规则的位点；在
                            --geno/--maf/--hwe 之前执行；边界值保留
```

`-S` 与 `--keep` 同给、`-s` 与 `-S/--keep` 同给、`-s ^…` 与 `--remove` 同给都是用法错。

### 位点选择

```
-v, --types LIST            保留这些变异类型：snps、mnps、indels、other、ref
-V, --exclude-types LIST    丢掉这些变异类型
    --snps-only[=just-acgt] 只留 SNP。裸写是任意单碱基替换；=just-acgt 额外要求每个
                            REF/ALT 都在 A/C/G/T 内
-m, --min-alleles INT       等位数下限，REF+ALT 计
-M, --max-alleles INT       等位数上限，REF+ALT 计。默认 2；-M 1 只留不变位点
-f, --apply-filters LIST    FILTER 列在 LIST 里的记录才留，逗号分隔多个值
-i, --include EXPR          保留满足表达式的记录，可重复（见过滤表达式页）
-e, --exclude EXPR          丢掉满足表达式的记录，可重复
-g, --genotype [^](hom|het|miss)
                            按基因型组成选位点：hom 存在纯合变异、het 存在杂合、
                            存在缺失；可重复，之间是 AND；^ 取反。半调用按缺失算
```

### 基因型掩码

```
    --min-dp INT            FMT/DP 低于 INT 的基因型当作没测到；DP 本身缺失的不受影响。
                            别名 --vcf-min-dp
    --max-dp INT            FMT/DP 高于 INT 的基因型当作没测到。别名 --vcf-max-dp
    --min-gq INT            FMT/GQ 低于 INT 的基因型当作没测到；GQ 缺失的不动。
                            别名 --vcf-min-gq
    --set-gts .             配合 -i/-e：FORMAT 表达式按样本逐个求值，命中的基因型置为
                            缺失，位点本身保留。只接受 .
```

掩码掉的基因型按缺失参与下面所有位点阈值。

### 频率与缺失阈值

```
    --geno FLOAT            缺失等位比例超过 FLOAT 的位点丢掉，按掩码之后计。按等位计：
                            半调用算 1 个已调用等位，全缺失的基因型算 2 个缺失；
                            边界值保留
    --maf FLOAT             次要等位频率低于 FLOAT 的位点丢掉（掩码之后）
    --max-maf FLOAT         次要等位频率上限
    --mac INT               次要等位计数下限
    --max-mac INT           次要等位计数上限。别写 --max-mac 0：0 是内部的"关闭"值，
                            会被拒绝；要 MAC ≤ 1 写 1
-q, --min-af FLOAT[:type]   等位频率下限。优先读 INFO/AC、INFO/AN（样本集没被改过时），
                            否则按基因型现算。:type 选口径：nref 非参考等位（默认）、
                            alt1 第一个 ALT、minor 次要、major 主要、nonmajor 非主要
-Q, --max-af FLOAT[:type]   等位频率上限，:type 同 -q
    --hwe P [midp]          HWE 精确检验 p 值低于 P 的位点丢掉；可选 midp 用 mid-p 值
```

`-q/-Q` 优先用 INFO 里的 AC/AN，`--maf` 总是按基因型现算。样本子集或掩码改动了样本集后，`-q/-Q` 改按基因型重算并提示"陈旧 INFO"。

### LD 剪枝

```
    --indep-pairwise W[kb] [S] R2
                            窗口、步长、r² 上限。窗口按位点数计（写 50kb 则按碱基距离，
                            此时步长必须是 1）。窗口内 r² 超过 R2 的位点对剪掉次要等位
                            频率更低的那个。输出就是剪枝后的 VCF
```

### 输出

```
-O, --output-type v|z|b|u[0-9]
                            v 纯文本 VCF、z BGZF VCF、b BCF、u 未压缩 BCF；
                            可带压缩级别（-Oz6、-Ob1）
-o, --output FILE           输出文件；不给写 stdout。指向输入文件本身会被拒绝。
                            不给 -O 时按扩展名猜（.gz 压缩、.bcf 走 BCF）。
                            别名 --out
-W, --write-index[=csi|tbi] 边写边建索引：-W 或 -W=csi（默认）、-W=tbi；不给不建。
                            只对 BGZF 输出有效，且必须给了 -o
    --report FILE           JSON 运行报告：每条规则的累计通过数、掩码归因、
                            --mind 丢弃数、LD 剪枝数、耗时
    --timing FILE           各阶段毫秒数
    --no-version            header 里不写 ##reed_* 溯源行
    --no-command            只省掉 ##reed_command 行（含机器路径的原始 argv），
                            其余溯源行保留
    --tmp-dir DIR           多遍运行的 stdin 溢写目录
    --threads INT           工作线程数，默认全部可用核
    --memory MB             在途内存预算，MiB，默认物理内存一半
    --verbosity INT         0 静默、1 摘要（默认）、2 阶段耗时、3 调试
```

### INFO 标签

```
-I, --no-update             不更新 INFO/AC、INFO/AN。默认样本集变了就按新样本集重算
-x, --remove-info LIST      删掉这些 INFO 标签。默认删 AC,AN,AF,MLEAC,MLEAF,
                            ExcessHet,InbreedingCoeff,DP；-x '' 表示一个都不删
    --fill-tags LIST        按当前样本集重算并写入：AC、AN、NS、MAF、HWE，
                            大小写不敏感；不支持的标签是用法错
```

## 求值顺序

每条记录按固定顺序过八个阶段。顺序决定漏斗计数，不影响最终结果；唯一影响结果的是 `--mind` 在掩码之后：

1. 样本子集与区域/位置选择
2. 载入期位点过滤：`-v/-V/--snps-only/-m/-M/-f`
3. 不依赖基因型的 `-i/-e`
4. 基因型掩码：`--min-dp/--max-dp/--min-gq`、`--set-gts .` 加 FORMAT 表达式
5. `--mind`
6. `--geno`、`--hwe`、`--maf/--max-maf/--mac/--max-mac`、依赖基因型的 `-i/-e`
7. LD 剪枝（`--indep-pairwise`，需要全局位点顺序）
8. INFO 更新（`-x` → AC/AN 自动更新 → `--fill-tags`）与输出

需要第一遍全局信息的规则（`--mind`、LD 剪枝）先扫一遍输入，第二遍再产出。

`--max-missing` 这类旧拼写没有静默别名：`--max-missing 0.9` 被拒绝，错误信息要求改用 `--geno 0.1`（含义相反，0.9 是必须调用的比例）。

## 例子

```bash
# 区域查询，需要先 reed index 建索引
reed filter raw.vcf.gz -r chr1,chr2:1000-200000 -R targets.bed -o region.vcf.gz -Ov

# 位点类型与表达式组合
reed filter raw.vcf.gz --snps-only=just-acgt -M 2 -f PASS -i 'QUAL>30' -e 'F_MISSING>0.1' -o snp.vcf.gz -Oz

# 样本子集，顺带按缺失率丢样本
reed filter raw.vcf.gz -s sampleA,sampleB --mind 0.2 -o subset.vcf.gz -Oz

# 基因型掩码：低深度、高深度、低 GQ 都当作没测到
reed filter raw.vcf.gz --min-dp 5 --min-gq 20 --max-dp 100 -o masked.vcf.gz -Oz

# 按表达式掩码，用 --set-gts 把命中样本的基因型置缺失
reed filter raw.vcf.gz -e 'FMT/DP<5 || FMT/GQ<20' --set-gts . -o masked.vcf.gz -Oz

# 位点阈值
reed filter clean.vcf.gz --geno 0.1 --maf 0.05 --hwe 1e-6 -o qc.vcf.gz -Oz

# LD 剪枝；输出就是剪枝后的 VCF
reed filter clean.vcf.gz --indep-pairwise 50 5 0.2 -Oz6 -o pruned.vcf.gz

# 重算 INFO 标签，不删任何旧标签
reed filter raw.vcf.gz --fill-tags AC,AN,NS,MAF,HWE -x '' -o annotated.vcf.gz -Oz

# 带 JSON 报告、分阶段计时和索引
reed filter raw.vcf.gz -o clean.vcf.gz -Oz6 -W --report clean.json --timing clean.timing
```

跑完 stderr 一行摘要加过滤漏斗：

```
reed filter: 500 -> 131 records (10 of 10 samples), strategy single_pass, 0.0s
funnel:
  site.biallelic_gate          500
  after_mask.missing           177
  after_mask.maf               131
  after_mask.hwe               131
```

`funnel` 的每行是一条规则过后的剩余位点数。`--mqf`、`--min-mean-dp` 这类未实现的参数在解析期就报错，不会跑一半才失败。
