# reed

[English](README.md) | [简体中文](README.zh-CN.md)

reed 是用 Rust 写的群体遗传学命令行工具。它读 VCF/BCF，按窗口或全基因组计算 π、dxy、Fst、θ、Tajima's D、PBS、SFS、f 统计量与 LD 衰减；也做硬过滤、LD 剪枝、QC 报告、PCA/DAPC、个体间矩阵、格式转换与索引。

## 安装

从 [GitHub Releases](https://github.com/WWz33/reed/releases) 下载对应平台的压缩包，运行时无需安装 Rust、Cargo、htslib
或 bcftools。

| 平台 | 压缩包后缀 | 最低系统要求 |
|---|---|---|
| Linux x86_64 | `x86_64-unknown-linux-gnu-glibc217.tar.gz` | glibc 2.17+，内核 3.2+ |
| Linux ARM64 | `aarch64-unknown-linux-gnu-glibc217.tar.gz` | glibc 2.17+，内核 4.1+ |
| macOS Apple Silicon | `aarch64-apple-darwin-macos11.tar.gz` | macOS 11.0 |
| macOS Intel | `x86_64-apple-darwin-macos11.tar.gz` | macOS 11.0 |

压缩包名称为 `reed-v<版本>-<后缀>`。同时下载 `SHA256SUMS`，Linux 使用
`sha256sum --ignore-missing --check SHA256SUMS`，macOS 使用
`shasum -a 256 --ignore-missing --check SHA256SUMS` 校验，然后解压并将 `reed`
放到 PATH。每个包内附许可证和 `BUILD-INFO.txt`，记录目标平台、源码提交与编译器版本。
macOS 二进制尚未签名或公证。

也可以在源码 checkout 中安装：

```bash
cargo install --locked --path crates/reed-cli
reed --version
```

源码编译支持 Rust 1.91 或更高版本，使用当前环境已启用的 Rust 工具链。编译随依赖提供的
libdeflate 需要 C 编译器，不需要另行安装 libdeflate 或 htslib。开发分支 CI 使用当前
stable 跑完整检查；main 仅检查 stable 和 Rust 1.91.1 MSRV 下的源码构建。
正式发布的二进制固定用 Rust 1.91.1 构建，并在四个平台上做冒烟检查。

`--threads N` 接受 1–1024。发布包使用 CPU 基础指令集，不要求与构建 runner 相同的 CPU 型号。

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

完整手册在 [`docsite/`](docsite/)：安装、每个命令一页的完整参数、逐步的快速开始、过滤表达式、输出约定与 Rust API。用 Astro Starlight 构建，发布在 <https://wwz33.github.io/reed/>。本地预览：

```bash
cd docsite && npm install && npm run dev
```

## 引用

使用到的方法请引用对应文献，完整列表见 [`docsite/src/content/docs/zh-cn/citation.md`](docsite/src/content/docs/zh-cn/citation.md)。

## 许可

MIT OR Apache-2.0。见 [LICENSE-MIT](LICENSE-MIT) 与 [LICENSE-APACHE](LICENSE-APACHE)。
