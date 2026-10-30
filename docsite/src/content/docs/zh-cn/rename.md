---
title: rename
---

## 用法

```
reed rename [OPTIONS] <INPUT>
```

标识符编辑：样本名、contig 名、变异 ID、header 行。别名 `reheader`。默认输出到 stdout，基因型列不参与任何编辑，数据内容不会改变。

## 参数

### header

```
-h, --header FILE           用文件内容替换整个 header；文件里的 #CHROM 行必须与输入的
                            样本完全一致（同名同序），否则报错。-h 在本命令里属于
                            --header，帮助要用 --help
-f, --fai FILE              从 .fai（名字<TAB>长度）设置 ##contig 长度：已有的改
                            length=，header 里没有的 contig 在 #CHROM 行之前插入新行
    --add-meta LINE         在 header 末尾追加 ##key=value 行，可重复
    --remove-meta KEY       删除所有 ##KEY= 行，可重复
```

### contig

```
    --rename-chrs FILE      旧名<TAB>新名 映射文件（空白分隔与 # 注释都接受）；
                            ##contig 行与记录的 CHROM 列一起改。旧名必须是 header
                            里已有的 contig
    --chr-style plain|ucsc|mt-is-m
                            规范化 contig 名：plain 去掉 chr 前缀（1、MT）；ucsc 加
                            chr（chr1、chrMT）；mt-is-m 加 chr 且线粒体写作 chrM。
                            先应用 --rename-chrs 的映射，再套风格
```

改名后若出现重复的 contig 名会被拒绝（会产生重复 `##contig` 行并合并两个 contig 的记录）。输出旁已有陈旧索引时打印警告，提示 `reed index -f` 重建。

### ID

```
    --set-id PATTERN        重写所有 ID。模板支持 %CHROM %POS %REF %ALT，%% 是字面 %，
                            如 %CHROM:%POS:%REF:%ALT
    --set-missing-var-ids PATTERN
                            同上，但只填 ID 为 . 的位点
    --update-name FILE      按 旧ID<TAB>新ID 映射改 ID；运行结束时任何一个没匹配到的
                            键都报错（拼错的键不能静默跳过）
    --idspace-to SEP        把 ID 里的 : 替换为 SEP
    --rm-dup-ids            重复 ID 加后缀 _1、_2…；. 不参与去重
    --clear-ids             所有 ID 置 .
```

### 样本

```
-s, --samples PAIRS         内联改名对，逗号分隔，每对两个空白分隔的名字：
                            -s "old new,old2 new2"
-N, --samples-file FILE     旧名<TAB>新名 每行一对
    --replace-sample-name PATTERN
                            字面量改写（无正则）：pre* 去前缀，*suf 去后缀，其余情况
                            删除文本出现的每一处。正则元字符报用法错
```

改名对的校验：旧名必须存在；新名不能为空；改名后撞名报错。

### 输出

```
-o, --output FILE           输出文件；默认 stdout。不给 -O 时按扩展名推断压缩。
                            别指向输入文件本身。别名 --out
-O, --output-type v|z[0-9]  纯文本 VCF 或 BGZF；可带压缩级别。BCF 输入先转 VCF 再
                            处理，-O b/-O u 被拒绝
-W, --write-index[=csi|tbi] 边写边建索引。索引按改名后的 header contig 表构建
                            （header 被替换而没有任何 ##contig 行时拒绝）；
                            需要 -o FILE 与 BGZF 输出
    --no-version            不写 ##reed_* 溯源行
    --no-command            不写 ##reed_command 行
    --threads INT           工作线程数，默认全部可用核
    --memory MB             在途内存预算，MiB，默认物理内存一半
    --verbosity INT         0 静默、1 摘要（默认）、2 阶段耗时、3 调试
```

## 求值顺序

编辑按固定顺序应用：header 替换 → contig 改名 → ID 操作 → 样本改名 → header 元数据增删。ID 操作里 `%CHROM` 用的是改名之后的 contig 名。

## 冲突校验

以下组合是用法错：

- `--set-id` 与 `--update-name`（模板会重写全部 ID，映射永远匹配不上）；
- `--clear-ids` 与任何其他 ID 参数；
- `--set-id`/`--set-missing-var-ids` 与 `--idspace-to`（把分隔符写进模板）；
- 一个编辑参数都不给。

## 输出细节

header 末尾插入 `##reed_version`、`##reed_command`、`##reed_rename=<触碰的平面>`（`--no-version`/`--no-command` 可关）。结束行摘要：记录数、contig 数、样本数、改动的 ID 数与耗时。

## 例子

```bash
# 样本改名 + contig 改名 + 补缺失 ID，一条命令
reed rename in.vcf.gz -s "s0 sampleA s1 sampleB" \
    --rename-chrs chrs.txt --set-missing-var-ids '%CHROM:%POS' \
    -Oz -o renamed.vcf.gz -W

# 批量样本改名（文件）
reed rename in.vcf.gz -N sample_map.tsv -Oz -o renamed.vcf.gz

# 去掉测序批次前缀
reed rename in.vcf.gz --replace-sample-name 'Batch3_' -Oz -o clean_names.vcf.gz

# ID 去冒号（有的下游工具不收冒号）+ 去重
reed rename in.vcf.gz --idspace-to '_' --rm-dup-ids -Oz -o ids.vcf.gz

# 全部 ID 重写为 chrom:pos:ref:alt
reed rename in.vcf.gz --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ids.vcf.gz

# 换 header、补 contig 长度、删旧的参考注释行
reed rename in.vcf.gz -h new_header.txt --remove-meta INFO -Oz -o newh.vcf.gz
reed rename in.vcf.gz -f ref.fa.fai -Oz -o contigs.vcf.gz

# contig 名规范化为 UCSC 风格
reed rename in.vcf.gz --chr-style ucsc -Oz -o ucsc.vcf.gz
```
