# convert

## 用法

```
reed convert [OPTIONS] --phylip --012 ... -o PREFIX <INPUT>
reed convert --vcf [OPTIONS] -o FILE <HAPMAP>
```

格式转换。一次读入，多目标输出；`-o` 缺省取输入文件名去掉扩展名。`--vcf` 是唯一的反向转换（hapmap 转 VCF），与所有 VCF 输入目标互斥。

## 参数

### 转换目标

给哪个写哪个，可以同时给多个。

```
    --phylip                PHYLIP 矩阵（relaxed 布局） → PREFIX.phy
    --fasta                 FASTA 矩阵 → PREFIX.fasta
    --nexus                 NEXUS（DNA）矩阵 → PREFIX.nexus
    --nexus-binary          二元 NEXUS（SNP 数据类型），只收双等位、二倍体位点 →
                            PREFIX.bin.nexus
    --012                   0/1/2/-1 基因型剂量矩阵 → PREFIX.012、
                            PREFIX.012.pos、PREFIX.012.indv
    --hapmap                Hapmap 矩阵（11 个固定列加每样本一列） → PREFIX.hmp.txt
    --msa                   参考驱动的全位点 MSA，每条 contig 或每个区域一个文件。
                            需要 --reference → PREFIX_<contig>[_<start>-<end>].fasta
    --sfs                   投影 SFS，与 reed popg --sfs 共用实现。需要 --proj →
                            PREFIX_sfs/dadi、PREFIX_sfs/fastsimcoal2
    --vcf                   反向：hapmap 输入转回 VCF；与其他目标互斥 → -o FILE
```

### 矩阵内容

```
    --outgroup NAME         这个样本写在比对的第一行
    --min-samples-locus N   一个位点至少 N 个样本有基因型才收。系统发育目标默认 4；
                            --012 与 --hapmap 只在给了本参数时才设门
    --resolve-iupac         把杂合位点定为一个碱基；抽签种子固定（0），
                            重复运行结果一致
    --write-used-sites      另写 PREFIX.used_sites.tsv，列出实际进入矩阵的位点及样本数
```

### MSA

```
    --reference FASTA       参考基因组；缺它时 --msa 报错。必须与 VCF 匹配
    --region CHR:B-E        只做这一段（1-based 闭区间），可重复；默认整条 contig
    --flank INT             每个区域两侧各扩 INT bp，默认 0
    --collapse-suffix STR   样本名在 STR 第一次出现处截断；截断后撞名的样本被拒绝
```

### SFS

```
    --sfs[=folded]          展开谱（默认，REF 为祖先）；--sfs folded 折叠为次要等位谱
    --proj N[,N...]         每个群体投影到的染色体数；一个值广播，多个值按群体首现
                            顺序对应
    --sfs-int               投影值取整
    --pop FILE              群体文件，样本<TAB>群体，无表头；不给时全部样本算一组。
                            别名 --populations
    --force-samples         群体文件里的未知样本只警告并跳过
```

### hapmap 转 VCF

```
    --vcf                   输入是 .hmp.txt（gzip/BGZF 自动识别，stdin 也可以）
-O, --output-type v|z[0-9]  纯文本 VCF 或 BGZF；没有 BCF 输出
    --no-version            不写 ##reed_* 溯源行
    --no-command            不写 ##reed_command 行
```

### 输入与样本

```
<INPUT>                     VCF/VCF.gz/BCF（--vcf 时是 hapmap）；- 读 stdin
-s, --samples [^]LIST       样本子集，逗号分隔，^ 前缀取补
-S, --samples-file FILE     样本名单文件，一行一个
    --keep FILE             同 -S，另一种拼写
    --remove FILE           要丢掉的样本名单
```

### 输出与全局参数

```
-o, --output PREFIX         输出前缀；缺省取输入文件名去掉扩展名。
                            --vcf 时是输出文件。别名 --out
    --threads INT           工作线程数，默认全部可用核
    --memory MB             在途内存预算，MiB，默认物理内存一半
    --verbosity INT         0 静默、1 摘要（默认）、2 阶段耗时、3 调试
```

## 通用规则

- 系统发育矩阵和 hapmap 只接受 REF 与全部 ALT 都是单字符的位点；多字符位点跳过并计数。`--012` 只跳过多等位位点。
- 一个目标的门只影响那个目标：系统发育流跳过的位点仍会进入 `--012`、hapmap、MSA。
- 畸形记录（列数不对、POS 坏、等位字母不认识）跳过并计数，只在 stderr 警告一次，总数写进 `PREFIX.log`。

## 系统发育矩阵

纯合写碱基本身，杂合写 IUPAC 简并码（`A/G→R`、`C/T→Y`、`A/C→M`、`G/T→K`、`A/T→W`、`C/G→S`；三碱基 `V/H/D/B`；四碱基 `N`），缺失写 `N`。

`PREFIX.phy`（relaxed 布局，首行是 `样本数 位点数`）：

```
3 5
sampleA     ACGTA
sampleB     ACGTC
sampleC     NGTTA
```

真实运行的样子（10 个样本、122 个位点，序列行截断显示）：

```
10 122
s1    MAGGRCAAMAMTARGNRNCTTCNRWRCGWACANGTTAGACRWGCGCMMCTWTMAAAAT…
s2    MRRAACACAGCWTGRARMAAWNTATGCNNTNGANWATAACRAGCGCCMAAWWCNANGT…
s3    MGAGNCCACGCTAAAMRMATTMTAARMGAWAARGWAAAGAGTRANCANMWNTACTTRT…
…
```

`PREFIX.fasta`：

```
>s1
MAGGRCAAMAMTARGNRNCTTCNRWRCGWACANGTTAGACRWGCGCMMCTWTMAAAATNANNTAATWTMWMAACAMMCANARATMMACTRRAAWWGAACC…
```

`PREFIX.nexus`：

```
#NEXUS

BEGIN DATA;
	DIMENSIONS NTAX=10 NCHAR=122;
	FORMAT DATATYPE=DNA MISSING=N GAP=- ;
MATRIX
s1    MAGGRCAAMAMTARGNRNCTTCNRWRCGWACANGTTAGACRWGCGCMMCTWTMAAAAT…
s2    MRRAACACAGCWTGRARMAAWNTATGCNNTNGANWATAACRAGCGCCMAAWWCNANGT…
;
END;
```

`PREFIX.bin.nexus`（`--nexus-binary`，`0/0→0`，`0/1→1`，`1/1→2`，其余含 `./.` 写 `?`）：

```
#NEXUS

BEGIN DATA;
	DIMENSIONS NTAX=10 NCHAR=122;
	FORMAT DATATYPE=SNP MISSING=? GAP=- ;
MATRIX
s1    102212001012012?1?2222?111221020?2220202112222112212100002?0??200212111002…
s2    111002020221221011001?20222??2?20?1020021022222100112?0?222022211?00?2100…
;
END;
```

## 基因型矩阵（`--012`）

写三个文件。行首下标是样本在原 VCF header 里的序号；缺失、半调用和单倍体一律写 `-1`：

```
PREFIX.012                                    PREFIX.012.pos    PREFIX.012.indv
0  1 0 2 2 1 2 0 0 1 … 2 -1 0 -1 -1 2 0 0 0 1  chr1  400         s1
1  1 1 1 0 0 2 0 2 0 … 0 0 1 0 0 1 0 1 2 2 0   chr1  800         s2
2  1 2 0 2 -1 2 2 0 2 … 0 -1 0 0 0 1 0 1 0 1   …                 …
```

## Hapmap（`--hapmap`）

`PREFIX.hmp.txt`，11 个固定列加每样本一列：

```
rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode	s1	s2	s3	s4	s5	s6	s7	s8	s9	s10
rs2	A/C	chr1	400	+	NA	NA	NA	NA	NA	NA	M	M	M	C	A	M	N	A	C	M
rs4	A/G	chr1	800	+	NA	NA	NA	NA	NA	NA	…	…	…	…	…	…	…	…	…	…
```

读法：第一行里 `s4` 是 `C`，即 A/C 纯合 ALT；`s1` 是 `M`，A/C 杂合；`s7` 是 `N`，缺失。几个要点：

- 格子用 IUPAC 单字母，再做一次字母表收缩：只保留 `A C G T N R Y S W K M -`，其余（小写字母、三碱基码）写 `N`。
- `rs#` 取 VCF 的 ID；ID 为 `.` 时写 `chrom:pos`。
- `alleles` 列写 REF/ALT 的声明集（斜杠连接、去重、`ALT='.'` 跳过、`*` 写 `-`）。
- strand 写 `+`，中间六个元数据列写 `NA`。

`--write-used-sites` 时另写的 `PREFIX.used_sites.tsv` 长这样：

```
#CHROM	POS	NUM_SAMPLES
chr1	400	9
chr1	800	9
```

## 参考 MSA（`--msa`）

每条 contig（或每个 `--region`）写一个 `PREFIX_<contig>[_<start>-<end>].fasta`，记录头是 `>contig:start-end_样本`：

```
>chr1:1-1000_s1
AAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA…
>chr1:1-1000_s2
AAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA…
```

位点上的样本值：单字符等位并入 IUPAC 码；indel 等位保留字面文本，全是 indel 时取最短的（平局按字典序）；完全缺失的基因型写参考碱基。各样本的值右填充 `-` 到该位点的最大长度，不调外部比对器。输入必须按 contig 排序；同一位置出现多条记录时取最后一条。

## 投影 SFS（`--sfs`）

与 `reed popg --sfs` 共用同一实现，同前缀同输入时输出一致。

dadi 布局（`PREFIX_sfs/dadi/<pop>-<n>.sfs`，三行：维度、谱值、mask）：

```
7 unfolded
2.7357142857142858 10.776190476190477 26.226190476190478 35.642857142857146 30.773809523809526 14.073809523809524 1.7714285714285716
1 0 0 0 0 0 1
```

fastsimcoal2 布局（`PREFIX_sfs/fastsimcoal2/<pop>_DAFpop0.obs`，展开谱用 DAF、折叠谱用 MAF）：

```
1 observation
d0_0	d0_1	d0_2	d0_3	d0_4	d0_5	d0_6
2.7357142857142858 10.776190476190477 26.226190476190478 35.642857142857146 30.773809523809526 14.073809523809524 1.7714285714285716
```

两个群体起另写每对的二维联合谱（`dadi/<A>-<B>.sfs`、`fastsimcoal2/<base>_jointMAFpop<idxB>_<idxA>.obs`）。上面输出对应的目录结构：

```
out_sfs/dadi/            popA-6.sfs  popB-6.sfs  popA-popB.sfs
out_sfs/fastsimcoal2/    popA_DAFpop0.obs  popB_DAFpop0.obs  out_jointDAFpop1_0.obs
```

## Hapmap 转 VCF（`--vcf`）

翻译规则：

- `CHROM/POS/ID` 直抄 chrom、pos、rs# 列；REF/ALT 按声明顺序直抄 `alleles` 列，不翻转、不重排。`alleles=NA` 的行整行基因型写 `./.`。
- 格子展开为 GT：`R → 0/1` 这类简并码反查展开；只有 `N` 和 `?` 是缺失（`./.`）；展开出的碱基必须是声明等位（大小写不敏感）；`-`、`+`、`0` 按 `alleles` 列查下标。
- 输出头带 `##fileformat=VCFv4.3` 和 `##FORMAT=<ID=GT>`，不写 `##contig` 行。
- 不接受样本子集和任何 VCF 输入目标的参数。

## 例子

```bash
# 一次写全部核苷酸矩阵与剂量矩阵
reed convert clean.vcf.gz --phylip --fasta --nexus --nexus-binary --012 --hapmap -o matrices

# 系统发育矩阵：指定外群、位点样本数门、杂合定碱基、记录用到的位点
reed convert clean.vcf.gz --phylip --outgroup sampleZ \
    --min-samples-locus 10 --resolve-iupac --write-used-sites -o matrices

# 参考驱动 MSA，只做一段并两侧各扩 500 bp
reed convert clean.vcf.gz --msa --reference ref.fa --region chr1:100000-200000 --flank 500 -o aln

# 投影 SFS
reed convert clean.vcf.gz --pop pops.txt --sfs --proj 20,20 -o out

# hapmap 转回 VCF
reed convert panels.hmp.txt --vcf -Oz6 -o panels.vcf.gz
```
