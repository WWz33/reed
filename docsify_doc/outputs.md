# 输出与日志

输出文件名、数值格式、日志与报告的公共约定都收在这一页；每个命令自己特有的输出文件和列，写在各自命令页里。

## -o 的两种含义

| 命令 | `-o` 含义 | 缺省行为 |
|---|---|---|
| `filter`、`rename` | 输出文件路径 | stdout |
| `index` | 索引文件路径 | 写在输入旁 |
| `convert` | 输出前缀（`--vcf` 时是文件路径） | 取输入文件名去扩展名；stdin 无 `-o` 报错 |
| `stats`、`popg`、`pca`、`distance` | 输出前缀 | 无缺省，空缺报错 |

`stats` 和 `popg` 没有 stdout 模式；`-o` 为空串报错（会在工作目录写出隐藏文件）。

## 数值格式

- 零分母、无覆盖、空窗口的估计量写 `NA`，不写 0。
- 坐标是 1-based 闭区间；BED 输入按 BED 的 0-based 半开约定读取并转成 1-based。
- 浮点输出多位小数（矩阵 9 位），分位数与比例原精度输出。

## 溯源 header

`filter` 和 `rename` 会在 header 里插入：

```
##reed_version=<版本>
##reed_command=<原始 argv>        # --no-command 省略
##reed_filterCommand=<规范化命令>  # 仅 filter
##reed_filter_spec=<规范化 JSON>   # 仅 filter
##reed_filter_sha256=<规范 JSON 哈希>
##reed_rename=<触碰的平面>         # 仅 rename
```

`--no-version` 省掉全部 `##reed_*` 行；`--no-command` 只省原始 argv 那行（它含机器路径）。

## PREFIX.log

`popg` 的窗口模式与 `convert` 都写 `PREFIX.log`，内容包括：reed 版本、完整命令行、每个统计量的方法与文献（如 π 的估计式与缺失数据处理方式）、不变位点模式、窗口与记录计数、`convert` 的目标与计数器（总记录数、缺样本跳过数、多字符位点数、畸形行数）。

## filter 的 JSON 报告

`--report FILE` 输出一轮过滤的完整账目：

```json
{
  "n_input_records": 500,
  "n_output_records": 373,
  "n_samples_in": 10,
  "n_samples_out": 10,
  "funnel": [
    { "rule": "site.biallelic_gate", "n_passed": 500 },
    { "rule": "after_mask.missing", "n_passed": 373 }
  ],
  "n_masked_gt": 610,
  "masked_by_rule": [
    { "rule": "genotype.mask.DP", "n_masked": 639 }
  ],
  "strategy": "single_pass",
  "elapsed_secs": 0.0
}
```

漏斗按规则列出累计通过数，规则 id 与[求值顺序](commands.md)的阶段对应（`site.regions`、`after_mask.missing`、`after_mask.ld_prune` 等）。`--mind` 的样本规则单列在 `sample_rule`；LD 剪枝计数在 `ld_prune`。`--timing FILE` 另给各阶段毫秒数。

## 退出码

| 码 | 含义 |
|---|---|
| 0 | 成功 |
| 1 | 运行期失败（错误带 `contig:pos` 定位） |
| 2 | 用法错误（参数冲突、标签未声明、表达式编译失败、群体缺失等） |
