# 过滤表达式

`filter` 的 `-i/--include` 与 `-e/--exclude` 接受表达式：`-i` 保留表达式为真的位点，`-e` 丢弃它们。两者都可重复，按命令行顺序依次应用。

## 值

| 类别 | 写法 | 求值时机 |
|---|---|---|
| 记录列 | `CHROM` `POS` `ID` `QUAL` `FILTER` `TYPE` | 掩码前 |
| INFO 标签 | `INFO/TAG` 或裸名 `TAG`（先查 INFO） | 掩码前 |
| FORMAT 标签 | `FMT/TAG` | 只在 `--set-gts .` 的按样本掩码里有意义 |
| 计算值 | `N_ALT` `N_SAMPLES` `N_MISSING` `F_MISSING` | `N_MISSING`、`F_MISSING` 在掩码后 |
| 等位统计 | `AC` `AN` `AF` `MAC` `MAF` | INFO 里有就取 INFO；样本集被改过就按基因型现算 |
| 字面量 | 数字（含负数）、字符串、`"."`（缺失） | |

`AC/AN/AF/MAC/MAF` 的现算公式：`AF = ac/an`，`MAC = an − max(ac, an−ac)`，`MAF = min(af, 1−af)`。表达式用到 INFO 里已有的 AC/AN 时会提示"陈旧 INFO"，建议改用 `--maf/--mac` 或先 `-x INFO/AC,INFO/AN`。

## 运算符

| 类别 | 运算符 |
|---|---|
| 比较 | `==` `=` `!=` `<` `<=` `>` `>=` |
| 正则 | `~` 匹配、`!~` 不匹配（只支持字面模式加 `^`/`$` 锚；正则元字符 `.*+?[](){}\|` 报错） |
| 逻辑 | `&&` `\|\|` `!`，括号分组；优先级 `!` > 比较 > `&&` > `\|\|` |

语义要点：

- 与缺失值的比较一律为假。`QUAL=="."` 匹配 QUAL 缺失的位点；`QUAL>30` 对缺失 QUAL 的位点为假。
- flag 型 INFO 标签按存在性读写：`DB=1` 测存在，`DB=0` 测缺失。裸 flag（单独一个 `DB`）是编译错。
- 多值 INFO（`Number=A/R/G/.` 或计数 >1）按"任一元素满足即真"比较。
- `Number=1` 的标签给了多个值按缺失处理；数值解析失败按缺失处理，不做字符串比较。

## 未支持

- FORMAT 下标 `FMT/DP[0]`、`FMT/DP[*]`（解析到即报错）；
- 任意算术（没有 `+ - * /`）；
- 通用正则；
- `--set-gts .` 的表达式混用 FORMAT 值与基因型派生的位点值；FORMAT 标签必须 `Number=1` 且数值型。

## 例子

```bash
# QUAL 与次要等位频率双阈值
reed filter in.vcf.gz -i 'QUAL>30 && MAF>=0.05' -Oz -o out.vcf.gz

# 深度不足的基因型置缺失，位点保留
reed filter in.vcf.gz -e 'FMT/DP<5' --set-gts . -Oz -o masked.vcf.gz

# 同上，用 GQ 联合判断
reed filter in.vcf.gz -e 'FMT/DP<5 || FMT/GQ<20' --set-gts . -Oz -o masked.vcf.gz

# 变异类型与 ID
reed filter in.vcf.gz -i 'TYPE="snp" && ID!="."' -Oz -o snps_named.vcf.gz

# FILTER flag 存在性
reed filter in.vcf.gz -e 'LowQual=1' -Oz -o keep.vcf.gz

# QUAL 缺失的位点单独处理
reed filter in.vcf.gz -e 'QUAL=="."' -Oz -o defined_qual.vcf.gz

# 两组硬过滤（indel 与 SNP 阈值不同）
reed filter in.vcf.gz -e 'QD<2 || QUAL<30 || SOR>3 || FS>60 || MQ<40 || MQRankSum<-12.5 || ReadPosRankSum<-8' -Oz -o snp.vcf.gz
reed filter in.vcf.gz -e 'QD<2 || QUAL<30 || FS>200 || ReadPosRankSum<-20' -Oz -o indel.vcf.gz
```
