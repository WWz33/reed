# 并行评审与修复（doc_v22，2026-11-03）

> 触发：R22 收口后的六切片并行评审（文本 IO、BCF、filter+CLI、stats、
> 测试/CI/gold/文档、跨切面），全部只读。本文件记录评审结论、本机复核证据
> 与处理结果；修复完成后逐项回填状态。
>
> 复核方式：本机 `reed 0.4.0`（debug 构建）+ `bcftools 1.24-2` +
> `pixy 2.2.4`，最小 fixture 在 `/tmp/v22verify`；未改评审前的代码。
>
> 第二轮（GQ/u8 索引/整数哨兵/`--hist dp` AD-only）见
> [03-review-pass2.md](03-review-pass2.md)，全部并入 v22，当前**待修**。

## 结论

评审 18 项中 **17 项属实**（含 5 项 must-fix），1 项（#2）按用户要求
跳过——`.gitignore`/提交项不在本轮修复范围。其余全部修复。

## 分级与处理

| # | 级别 | 结论 | 位置 / 证据 | 处理 |
|---|---|---|---|---|
| 1 | must | **属实**：sites-only BCF 渲染多一列 | `bcf/render.rs:199`；`bcftools view -Ob sites.vcf` 后 `reed filter -Ov` 输出 9 列（末尾多 `.`），htslib 为 8 列 | 已修 |
| 2 | must | **属实（提交风险）**：新测试/文档/gold 脚本被 ignore | `crates/reed-io/.gitignore:1`、根 `.gitignore:19-20`；`git status --ignored` 确认 | **按用户要求不修** |
| 3 | must | **属实**：`ci/test-parallel.sh` 可假绿 | `ci/test-parallel.sh:17,30-47`；`JOBS=abc` 实测 exit 0 + 0 passed；假 cargo 模拟集成测试编译失败也 exit 0 | 已修 |
| 4 | must | **属实**：`--hist gq`/`FMT/TAG` 尾随空列虚增 missing | `stats/writers/reports.rs:694`；2 样本 4 个 GQ + 行尾 tab → `#missing 1`，对照无此行 | 已修 |
| 5 | must | **属实**：`run_pixy.sh` 默认输入必红 | `gold/popgen/run_pixy.sh:10-13,25-33`；pixy 2.2.4 对 variant-only 输入报 `no invariant sites`，exit 1 | 已修 |
| 6 | should | **属实**：`to_args` 丢 ploidy/depth | `spec/render.rs:30-193`；`--ploidy allow-haploid` 输出头无 `##reed_filterCommand`；`depth=false` DP 掩码渲染成 `--min-dp` | 已修 |
| 7 | should | **属实**：DP Integer + AD Float 时深度掩码静默失效 | `plan/mod.rs:627-657`；实测无 DP 的记录 `0/1:3.5,4.5` 未被掩码且无诊断 | 已修 |
| 8 | should | **属实**：`rename` 拒绝 sites-only | `reed/src/rename/mod.rs:411-417`；`--rename-chrs` 实测 exit 2，`cli.md:56` 却承诺支持 | 已修 |
| 9 | should | **属实**：`--min-dp` 无 CLI 层回退测试 | e2e 直接构造 `depth: true`；`stdin_streaming.rs:492` 是 DP-present 输入 | 已修 |
| 10 | should | **属实**：`bcf_vs_noodles` 绝对路径 + 文档称“仓库内” | `tests/bcf_vs_noodles.rs:203`；`CHANGELOG.md:67`、`06-roadmap.md:57` | 已修 |
| 11 | should | **属实**：FORMAT 描述符遍历重复，注释理由不成立 | `bcf/columns.rs:64-76` vs `bcf/render.rs:213-224`，两边都先收 `fields` | 已修 |
| 12 | should | **属实**：v22 文档指向已删除/搬移路径 | `v22/README.md:29,42,46,48`、`01-review.md:19,156` 仍写 `spec_for_key`/`native_encode.rs` | 已修 |
| 13 | should | **属实**：`--hist dp` 的 dp>0 与 AD 回退无直接测试 | `e2e_hist.rs:86` 只覆盖 DP-present fixture | 已修 |
| 14 | nice | **属实**：五份 plan→IntSpec 映射拷贝 | `filter/mod.rs:708`、`single_pass.rs:203`、`pass0/matrix.rs:143,1020`、`filtered.rs:130` | 已修 |
| 15 | nice | **属实**：死/投机 API（`IntSpec::scalar` 实际有测试调用，保留） | `SampleColumnsError::found/expected/Display` 0 调用；`FormatLayoutCache::new` 仅测试；`IntMode::Sum` 无生产调用 | 已修（删除 `IntMode`/`sum`/`FormatLayoutCache::new`/error 访问器；保留 `scalar`） |
| 16 | nice | **属实**：尾随空列规则在 `scan()` 内两份 | `genotype_scan.rs:165` 与 `:192-200`；validate 通过后循环内非空分支不可达 | 已修 |
| 17 | nice | **属实**：v22 CLI 两例只断言退出码 | `cli_contract.rs:1105-1110,1129-1133` | 已修 |
| 18 | nice | **属实**：陈旧文字 / 断链 | `e2e_stats.rs:677`；`CHANGELOG.md:74,81`；`v22/README.md:58`；`cargo doc` 报 `super::encode::bcf_header_text` 未解析 | 已修 |

## 复核证据摘录

### #1 sites-only BCF

```text
$ bcftools view -Ob -o sites.bcf sites.vcf
$ reed filter sites.bcf -Ov -o out.vcf
$ grep -v '^#' out.vcf | awk -F'\t' '{print NF}'   # 9 9
$ bcftools view -H sites.bcf | awk -F'\t' '{print NF}'  # 8 8
```

### #3 test-parallel 假绿

```text
$ JOBS=abc ci/test-parallel.sh
xargs: invalid number "abc" for -P option
test binaries: 0 passed, 0 failed, 0s wall, abc jobs
$ echo $?   # 0
```

假 cargo 返回构建失败、但输出一个测试可执行文件（`--doc` 成功）时：
`test binaries: 83 passed, 0 failed`，exit 0。

### #4 hist missing

```text
2 样本 × 2 记录，4 个 GQ，第二条记录行尾多一个 tab：
$ reed stats trail_gq.vcf --hist gq -o out
$ grep missing out.gq.hist
#missing        1
去掉行尾 tab 的对照文件没有 #missing 行。
```

### #5 pixy default

```text
$ python3 -m pixy --stats pi dxy --vcf variant_only.vcf.gz --populations pops.txt \
    --window_size 100000 --output_folder pixy_default --n_cores 1
ValueError: The provided VCF appears to contain no invariant sites (ALT = ".")
exit 1
```

### #7 DP+Float AD

```text
header: DP=Integer, AD=Float
$ reed filter in.vcf --min-dp 10 -Ov
chr1 10 ... GT:AD  0/1:3.5,4.5  1/1:1.5,2.5   # 无 DP 键 → 未掩码
chr1 20 ... GT:DP  ./.:5        ./.:6          # 有 DP → 正常掩码
```

### #8 rename sites-only

```text
$ reed rename sites.vcf --rename-chrs map.txt -o out.vcf
reed: the input has no sample columns (a sites-only VCF); the CHROM/ID planes need 9-column records
exit 2
```

### #18 断链

```text
$ cargo doc -p reed-io --no-deps
warning: unresolved link to `super::encode::bcf_header_text`
  --> crates/reed-io/src/bcf/encode_native.rs:14:7
```

## 不修项

- **#2**：`.gitignore` 吞掉新文件属于提交流程问题，按用户要求跳过；提交时
  对 `crates/reed-io/tests/format_columns_vs_bcf.rs`、`docs/rounds/v22/`、
  `gold/scripts/run_parallel.py` 等新文件使用 `git add -f`。

## 修复记录

逐项结果见 [STATUS.md](STATUS.md) 的“V22-R 评审修复”表。
