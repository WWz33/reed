# Installation

## Build requirements

Rust 1.91 or newer. `rust-toolchain.toml` in the repository root records the toolchain used for development; with rustup installed, entering the repository directory switches to it automatically.

reed is a pure Rust implementation: VCF/BCF parsing, BGZF encoding and decoding, and index writing all live inside the crate. Building and running need no extra C library dependencies.

## Install the command line from source

```bash
cargo install --path crates/reed-cli
reed --version
```

`cargo install` puts the binary in `~/.cargo/bin`; make sure that directory is on PATH. The long form of `reed --version` carries build information:

```
reed 0.4.0 (v0.3.0-112-g53729ef 2026-10-28, rustc 1.99.1)
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

A path dependency suits use inside the repository. The stability promise covers the crate root (`crates/reed`) only; crates below `crates/` are internal implementation whose types and signatures may change at any time. See [Library use](library.md).
