---
title: 安装
---

## 编译要求

Rust 1.91 或更高版本。仓库根目录的 `rust-toolchain.toml` 记录了开发用的工具链版本，装好 rustup 之后进入仓库目录会自动切换。

reed 是纯 Rust 实现，VCF/BCF 解析、BGZF 编解码、索引写入都在 crate 内部完成，编译和运行不需要额外的 C 库依赖。

## 从源码安装命令行

```bash
cargo install --path crates/reed-cli
reed --version
```

`cargo install` 把二进制装到 `~/.cargo/bin`，确认它在 PATH 里。`reed --version` 的长格式带构建信息：

```
reed 0.4.0 (v0.3.0-112-g53729ef 2026-10-28, rustc 1.99.1)
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
