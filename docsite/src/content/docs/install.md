---
title: Installation
---

## Build requirements

Prebuilt archives are available from [GitHub Releases](https://github.com/WWz33/reed/releases) for Linux x86_64/ARM64 (glibc 2.17 or newer) and macOS Intel/Apple Silicon (macOS 11.0 or newer). Runtime use does not require Rust, Cargo, htslib or bcftools. The macOS binaries are not signed or notarized.

Source builds require Rust 1.91 or newer and use the Rust toolchain active in your environment. A C compiler is needed for bundled libdeflate; no separately installed libdeflate or htslib is needed. CI tests current stable and checks the 1.91.1 MSRV; tagged release binaries are built with 1.91.1.

## Install a release archive

Choose the matching release archive, extract it, and put its `reed` executable on `PATH`:

| Platform | Archive suffix | Minimum system |
|---|---|---|
| Linux x86_64 | `x86_64-unknown-linux-gnu-glibc217.tar.gz` | glibc 2.17+, kernel 3.2+ |
| Linux ARM64 | `aarch64-unknown-linux-gnu-glibc217.tar.gz` | glibc 2.17+, kernel 4.1+ |
| macOS Apple Silicon | `aarch64-apple-darwin-macos11.tar.gz` | macOS 11.0 |
| macOS Intel | `x86_64-apple-darwin-macos11.tar.gz` | macOS 11.0 |

Verify the archive against `SHA256SUMS` before extracting it:

```bash
sha256sum --ignore-missing --check SHA256SUMS
```

On macOS, use `shasum -a 256 --ignore-missing --check SHA256SUMS` instead. Each archive contains the licenses and `BUILD-INFO.txt` with its source commit and compiler version.

## Install the command line from source

```bash
cargo install --locked --path crates/reed-cli
reed --version
```

`cargo install` puts the binary in `~/.cargo/bin`; make sure that directory is on PATH. The long form of `reed --version` carries build information:

```
reed 0.4.0 (<commit> <build-date>, rustc 1.91.1)
```

The parentheses hold the git commit, build date and rustc version; a build from a git tarball shows `unknown` instead.

To build without installing, the binary lands in `target/release/reed`:

```bash
cargo build --release
```

## Runtime resources

- Threads: `--threads N`, default `std::thread::available_parallelism()`, i.e. every available core.
- Memory: `--memory MB`, in MiB, default half of physical memory. The budget constrains intermediate data, not the process peak.
- Temporary directory: when `filter` needs multiple passes over stdin, it spills the input to `--tmp-dir` (default: the system temporary directory).

## Using it as a library

To use the library without the command line:

```toml
[dependencies]
reed = { path = "crates/reed" }
```

A path dependency suits use inside the repository. The stability promise covers the crate root (`crates/reed`) only; crates below `crates/` are internal implementation whose types and signatures may change at any time. See [Library use](../library/).
