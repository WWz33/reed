---
title: 安装
---

## 编译要求

[GitHub Releases](https://github.com/WWz33/reed/releases) 提供 Linux x86_64/ARM64（glibc 2.17 或更高）和 macOS Intel/Apple Silicon（macOS 11.0 或更高）的预编译包。运行时不需要 Rust、Cargo、htslib 或 bcftools。macOS 二进制尚未签名或公证。

源码编译需要 Rust 1.91 或更高版本，使用当前环境已启用的 Rust 工具链。编译随依赖提供的 libdeflate 需要 C 编译器，不需要另行安装 libdeflate 或 htslib。CI 使用当前 stable 测试，并检查 Rust 1.91.1 MSRV；正式发布的二进制固定用 Rust 1.91.1 构建。

## 安装发布包

下载对应平台的压缩包，解压后将 `reed` 可执行文件放到 `PATH`：

| 平台 | 压缩包后缀 | 最低系统要求 |
|---|---|---|
| Linux x86_64 | `x86_64-unknown-linux-gnu-glibc217.tar.gz` | glibc 2.17+，内核 3.2+ |
| Linux ARM64 | `aarch64-unknown-linux-gnu-glibc217.tar.gz` | glibc 2.17+，内核 4.1+ |
| macOS Apple Silicon | `aarch64-apple-darwin-macos11.tar.gz` | macOS 11.0 |
| macOS Intel | `x86_64-apple-darwin-macos11.tar.gz` | macOS 11.0 |

解压前用 `SHA256SUMS` 校验：

```bash
sha256sum --ignore-missing --check SHA256SUMS
```

macOS 使用 `shasum -a 256 --ignore-missing --check SHA256SUMS`。每个包都包含许可证和记录源码提交、编译器版本的 `BUILD-INFO.txt`。

## 从源码安装命令行

```bash
cargo install --locked --path crates/reed-cli
reed --version
```

`cargo install` 把二进制装到 `~/.cargo/bin`，确认它在 PATH 里。该命令使用当前已启用的 Rust 工具链。`reed --version` 的长格式带构建信息：

```
reed 0.4.0 (<提交> <构建日期>, rustc 1.91.1)
```

括号里是 git 提交、构建日期和 rustc 版本，从 git tarball 构建时显示 `unknown`。

也可以只构建不安装，二进制在 `target/release/reed`：

```bash
cargo build --release
```

## 运行时资源

- 线程：`--threads N`，默认 `std::thread::available_parallelism()`，即全部可用核。
- 内存：`--memory MB`，单位 MiB，默认物理内存的一半。预算约束的是中间数据，不等于进程峰值的硬上限。
- 临时目录：`filter` 在 stdin 输入需要多遍处理时会把输入溢写到 `--tmp-dir`（默认系统临时目录）。

## 作为库引用

只需要库、不要命令行的话：

```toml
[dependencies]
reed = { path = "crates/reed" }
```

路径依赖适合仓库内引用。稳定承诺只在 crate 根（`crates/reed`），`crates/` 下层的 crate 是内部实现，类型和签名随时可能变。用法见[库调用](../library/)。
