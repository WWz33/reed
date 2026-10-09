# reed

[English](README.md) | [简体中文](README.zh-CN.md)

**reed** 是一个用纯 Rust 编写的多线程 VCF/BCF 工具，用于群体遗传学统计与
硬过滤。它是一个**命令行产品**：`reed` 二进制是交付物，各 crate 是其内部结构；
只有 CLI 参数与输出文件格式是稳定性承诺（doc_v5 08 §6）。

不依赖 htslib、不含 GPL：传输层基于 noodles，加上 reed 自己的 BGZF 与 BCF 编解码。

## 命令

| 命令 | 作用 | 对标 |
|---|---|---|
| `reed filter` | 硬过滤、基因型掩码、INFO 重算；保真的 VCF/BCF 输出 | bcftools `view` + `filter -S .` + `+fill-tags` + `annotate -x`；plink2 的 QC 阈值 |
| `reed stats` | 一条命令下的全部统计量：π、dxy、Fst、Watterson θ、Tajima's D、PBS、SFS、逐位点/逐样本报告、bcftools stats 汇总 | vcftools、pixy、plink2、easySFS、`bcftools stats` |
| `reed index` | CSI（默认）或 TBI | `bcftools index` |

参数以 `reed <命令> --help` 为准；完整规范见合并版 CLI 文档
（`docs/rounds/v5/02-cli.md`，搬迁后为 `docs/spec/cli.md`）。

```bash
# 全位点输入的窗口统计，56 线程
reed stats allsites.vcf.gz --pop pops.txt --pi --dxy --fst --tajima-d --window 10k -o out

# 仅变异位点；窗口与全基因组 Fst
reed stats clean.vcf.gz --pop pops.txt --variant-only --pi --fst --window 100k --step 50k -o out

# 逐位点/逐样本报告 + bcftools stats 汇总，一遍读取
reed stats clean.vcf.gz --freq --hardy --het --missing --depth --summary -o qc

# 硬过滤，掩码后位点门 + 索引
reed filter in.vcf.gz -o out.vcf.gz -Oz6 --geno 0.1 --maf 0.05 --min-dp 5 --min-gq 20 -W
```

`-o` 对 `filter`/`index` 是文件，对 `stats` 是输出**前缀**（每个统计量各写一个文件）。
统计参数按科学方法命名（`--pi`、`--tajima-d`、`--fst wc84|hudson`）；软件名只出现在
对拍套件里。

## 验证与性能

正确性以参照工具为准，而非自洽。π/dxy 对 pixy（最大相对误差 4.8e-14）、
Fst 对 pixy（2.9e-15）、报告对 plink2/vcftools 逐行一致、汇总对 `bcftools stats`
九节逐行一致；gold 过滤规则与 bcftools 对拍（1,768,759 位点，`check_parity` PASS）。
输出跨线程数逐字节一致。

全量 `clean.vcf.gz`（7.0 GB，3,315,604 条，348 样本），56 线程、热缓存、3 次中位数：
`filter` `-Oz6` 8.79 s；`stats` 五统计 1.74 s（VCF）/ 1.51 s（BCF）；
reports + summary 3.76 s；`index` 1.60 s（VCF.gz）/ 1.05 s（BCF）。

## 仓库结构

```
crates/   八个 crate；唯一的二进制是 reed-cli
gold/     验证：夹具与对拍脚本
bench/    度量套件与结果日志
docs/     规范（spec/ 为当前，rounds/ 为历史）
ci/       漂移守卫与 CI 入口
```

## 构建与测试

```bash
cargo test --release --workspace
cargo build --release --workspace --all-targets
bash ci/check-doc-drift.sh
```

## 许可

MIT OR Apache-2.0。
