# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed 是用 Rust 写的群体遗传学命令行工具。它读 VCF/BCF，按窗口或全基因组计算 π、dxy、Fst、θ、Tajima's D、PBS、SFS、f 统计量与 LD 衰减；也做硬过滤、LD 剪枝、QC 报告、PCA/DAPC、个体间矩阵、格式转换与索引。

## 安装

```bash
cargo install --path crates/reed-cli
reed --version
```

需要 Rust 1.91 或更高版本。

## 快速开始

```bash
# 建索引
reed index raw.vcf.gz

# 硬过滤：先掩码浅深度基因型，再按位点阈值过滤
reed filter raw.vcf.gz -o clean.vcf.gz -Oz6 --min-dp 5 --min-gq 20 --geno 0.1 --maf 0.05 -W

# 只含变异位点的文件上做窗口统计
reed popg clean.vcf.gz --pop pops.txt --variant-only --pi --fst --window 100k -o out

# QC 报告与直方图
mkdir -p qc
reed stats raw.vcf.gz --freq --missing --depth --hist qual,dp-site,gq -o qc/raw

# LD 剪枝，然后 PCA
reed filter clean.vcf.gz --indep-pairwise 50 5 0.2 -Oz6 -o pruned.vcf.gz
reed pca pruned.vcf.gz --pcs 20 -o pcs

# 个体间矩阵与建树
reed distance clean.vcf.gz --type grm -o rel
reed distance clean.vcf.gz --type 1-ibs --tree nj --bootstrap 100 -o tree
```

## 文档

完整手册在 [`docsify_doc/zh-cn/README.md`](docsify_doc/zh-cn/README.md)：安装、每个命令一页的完整参数、逐步的快速开始、过滤表达式、输出约定与 Rust API。本地预览：

```bash
docsify serve docsify_doc
```

## 引用

使用到的方法请引用对应文献，完整列表见 [`docsify_doc/zh-cn/citation.md`](docsify_doc/zh-cn/citation.md)。

## 许可

MIT OR Apache-2.0。见 [LICENSE-MIT](LICENSE-MIT) 与 [LICENSE-APACHE](LICENSE-APACHE)。
