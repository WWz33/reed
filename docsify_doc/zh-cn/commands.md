# 命令与全局约定

`reed --help` 按组列出命令：Filtering（`filter`）、Statistics（`stats`、`popg`、`pca`）、Relatedness（`distance`）、Transform（`convert`、`rename`）、File Operations（`index`）。`reed <command> --help` 列出某个命令的全部参数。

## 全局参数

每个子命令都接受：

```
<INPUT>                 输入 VCF/VCF.gz/BCF；- 读 stdin（index 拒绝 stdin，
                        索引需要可定位的文件）
    --threads INT       工作线程数，默认全部可用核
    --memory MB         在途内存预算，单位 MiB，默认物理内存一半
    --verbosity INT     0 静默；1 摘要（默认）；2 阶段计时；3 调试日志
                        （走 RUST_LOG，默认 debug）；大于 3 报错
```

仅 `filter` 有：

```
    --tmp-dir DIR       多遍运行的 stdin 溢写目录
    --report FILE       JSON 运行报告：过滤漏斗、样本规则、计数器、策略
    --timing FILE       分阶段计时（key-value 文本）
```

`filter` 和 `rename` 另有 `--no-version`（省掉全部 `##reed_*` 溯源行）和 `--no-command`（只省 `##reed_command`，即含机器路径的原始 argv）。

## 长度与数值写法

长度参数（`--window`、`--step`、`--block-size`、`--flank`、`--ld-decay`、`--indep-pairwise` 的窗口）接受 `100000`、`100k`、`100kb`、`1M`、`1Mb`，后缀大小写不敏感，k 是 10³，M 是 10⁶。缺失率阈值（`--geno`、`--mind`、`--maf`、`--hwe`）是小数，范围 [0,1]。

## 共享参数组

三组参数在多个命令间拼写与语义一致：

- 样本子集：`-s/--samples [^]LIST`（逗号分隔，`^` 前缀取补）、`-S/--samples-file FILE`、`--keep FILE`、`--remove FILE`、`--force-samples`。`filter`、`stats`、`popg`、`pca`、`distance`、`convert` 使用；`-S` 和 `--keep` 同给、`-s` 与 `-S/--keep` 同给、`-s ^…` 与 `--remove` 同给都是用法错。
- 区域选择：`-r/--regions LIST`、`-R/--regions-file FILE`（需要索引）；`-t/--targets LIST`、`-T/--targets-file FILE`（流式，不需要索引）。文件格式自动识别：三列及以上按 BED 读（0-based 半开），两列按 `CHROM POS`（1-based），一列按区域串（`Chr01`、`Chr01:100-200`）。`#` 注释行与空行跳过；`-T` 文件里 `^` 前缀表示排除。
- 位点级过滤与阈值：`-v/-V`、`--snps-only`、`-m/-M`、`-f`、`-i/-e`、`-q/-Q`、`--geno`、`--maf/--max-maf`、`--mac/--max-mac`、`--hwe`。`filter` 全部实现；`popg` 里有同名的位点规则，含义是"边过滤边统计"；`stats` 的报告模式拒绝这些参数，要求先跑 `filter`。

## 退出码与错误格式

| 退出码 | 含义 |
|---|---|
| 0 | 成功 |
| 1 | 运行期失败（错误信息带 `contig:pos` 定位） |
| 2 | 用法错误：参数冲突、标签未声明、群体不存在、表达式编译失败、不变位点检查未通过、已撤回的参数 |

错误写 stderr，格式 `reed: <原因>`，警告格式 `reed: warning: <原因>`。解析期报错的参数给出改写建议，例如 `--max-missing 0.9` 会被拒绝并提示 `use --geno 0.1 (the fraction that may be missing)`；没有静默别名。

## 求值顺序（filter 管线）

`filter` 的每条记录按固定顺序过八个阶段，顺序决定漏斗计数，不影响最终结果；唯一影响结果的是 `--mind` 在掩码之后：

1. 样本子集与区域/位置选择
2. 载入期位点过滤：`-v/-V/--snps-only/-m/-M/-f`
3. 不依赖基因型的 `-i/-e`
4. 基因型掩码：`--min-dp/--max-dp/--min-gq`、`--set-gts .` 加 FORMAT 表达式
5. `--mind`
6. `--geno`、`--hwe`、`--maf/--max-maf/--mac/--max-mac`、依赖基因型的 `-i/-e`
7. LD 剪枝（`--indep-pairwise`，需要全局位点顺序）
8. INFO 更新（`-x` → AC/AN 自动更新 → `--fill-tags`）与输出

需要第一遍全局信息的规则（`--mind`、LD 剪枝）会先扫一遍输入，第二遍再产出输出。
