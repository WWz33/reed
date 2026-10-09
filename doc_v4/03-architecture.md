# 03 架构现状与残留（V0–V6 之后）

> doc_v2/03-architecture.md 的目标架构已基本建成。本文记录建成形态、与目标的残余差距、以及本轮新增的结构性建议。

---

## 1. 建成形态

### 1.1 crate 图（8 个）

```
reed-core     模型 + 内核 + 数学 + TagTable             依赖：bitflags、serde、thiserror、tracing
reed-io       格式、BGZF、chunk 源、VCF/BCF 解码与编辑、写出、索引    依赖：reed-core、noodles、libdeflater、memchr、flate2、bstr
reed-engine   chunk 流水线、票据、字节许可、重排、中止、长度前缀链     依赖：reed-core、reed-io、crossbeam-channel
reed-stats    统计量累加器 + QC summary（原 reed-qc）                 依赖：reed-core、serde
reed-filter   Spec、表达式、plan、求值、sketch、provenance（无 noodles）依赖：reed-core、reed-io、memchr、serde、sha2
reed          驱动：filter/stats/index、输入统一入口、BCF 包装         依赖：以上 + rayon、noodles、bstr、serde_json
reed-cli      二进制 `reed`：clap、单位、翻译                         依赖：reed + 各 crate
reed-core-v02-oracle  冻结的 v0.2（dev-dependency）
```

- **reed-filter 的纯度目标已达成**（doc_v2 03 §7、doc_v3 C3）：`Cargo.toml` 不含 noodles/clap/toml/crossbeam；`plan::compile` 接收 `reed_core::tag_table::TagTable`（由 `reed_io::header::tag_table` 从 header 构建）。`doc_v3/STATUS.md` 仍记为"未做"（01 F4-2）。
- reed-filter 仍依赖 reed-io（`LineRef`、`edit_vcf_line` 等在 eval 中使用）。L3 "只依赖 reed-core" 的最终目标未达成，但依赖方向正确（无环），不是阻塞项。

### 1.2 关键机制（已验证）

| 机制 | 位置 | 状态 |
|---|---|---|
| BAPP（块对齐并行解析） | `reed-io/src/chunk/bgzf.rs` | reader 直读进 chunk 缓冲（V3 优化） |
| 真正的票据 + 字节许可 | `reed-engine/src/pipeline.rs:180-262`、`budget.rs` | 在途 ≤ tickets，`in_flight_peak` 随线程数线性（56 线程 912 MB） |
| 长度前缀链（BCF 定界） | `reed-engine/src/frame_chain.rs`、`reed/src/bcf_input.rs:129-200` | 先发布后解码、定向唤醒、毒化、120 s 看门狗 |
| 统一记录管线 | `reed/src/filter/record_pipeline.rs` | single/two-pass/sketch 与 `FilteredJob` 共用 |
| 输出层 | `reed-io/src/write/{block_writer,ordered_output,index}.rs` | 记录可跨块、按序 append、边写边建 CSI/TBI |
| BCF 原生编码器 | `reed-io/src/bcf/native_encode.rs` | 与 noodles 编码器逐字节一致（测试锁定） |
| BCF 原生 GT 路径 | `reed-io/src/bcf/gt.rs`（64 K 查表） | 只用于无过滤计划的 stats |

---

## 2. 残余差距

| # | 目标（doc_v2/v3） | 现状 | 建议 |
|---|---|---|---|
| A4-1 | 消费者面向格式无关的 `SiteFields`/`GenotypeSource`（doc_v3 02 §2.3） | 没有抽象成 trait；BCF 通过"渲染成文本"复用文本消费者，只有 stats 有原生入口（`BcfRecordFn`） | 不强求 trait；但把 `RecordPipeline` 的输入从 `&[u8]` 文本行改为"`SiteView` + `ScanOutput`"两段式接口，让 BCF 解码器直接产出这两样（01 F4-5） |
| A4-2 | `FilteredJob` 向内层传语义对象 | 传重写后的文本行 | 内层 job 接收 `(SiteView, &[GtCode], columns)`（01 F4-7） |
| A4-3 | `reed stats` 报告两族一次读取 | 读两遍 | 合并 partial（01 F4-6） |
| A4-4 | 命名：禁止版本后缀 | `FilterSpecV2` 15 处 | 改名（01 F4-12） |
| A4-5 | `--timing` 覆盖所有阶段 | 只有 open/pass0/mid/pass1 + worker_busy | 加 reader 的 `next_chunk` 累计、reducer 的 `accept` 累计、链等待累计（BCF），否则"reader 受限"只能靠推断 |
| A4-6 | 原生 BCF 手术（doc_v3 B3b） | 未做（文档标注可选） | 维持不做；A4-1 落地后 BCF 输入的额外成本只剩输出侧的渲染+编码 |
| A4-7 | 文档树 | `docs/architecture-spec.md`（v1）、doc_v2、doc_v3、doc_v4 并存，互相引用 | 新增 `docs/INDEX.md`：一页说明哪个文档是当前规范、哪些是历史；每个旧目录 README 顶部加"已被 doc_vN 取代的章节"列表 |

---

## 3. 本轮新增的结构性建议

### 3.1 参数到 spec 的字段命名与方向一致（配合 02）

`TailRules.max_missing: Option<f64>` 现存丢弃比例（`--geno` 的值）。改为 `min_call_rate: Option<f64>`（保留比例），并新增 `max_missing_alleles: Option<u64>`（`--max-missing-count`）、`FilterSpec.min_sample_call_rate: Option<f64>`（取代 `mind`）。字段名、CLI 名、漏斗 id、`to_args()` 四处同向，往返测试锁定。

### 3.2 阈值判定抽成一个类型

```rust
/// A "keep if called fraction ≥ X" threshold, stored as a rational q/10^9
/// so boundary cases are exact (02-cli-semantics.md §4).
pub struct CallRate { q: u64 }
impl CallRate {
    pub fn parse(text: &str) -> Result<Self, SpecError>;      // 0 ≤ X ≤ 1
    pub fn keeps(&self, called: u64, total: u64) -> bool;     // called·D ≥ q·total (u128)
    pub fn render(&self) -> String;                            // shortest decimal
}
```

位点与样本规则共用它；plink2 的 ε 规则只在对拍测试里出现。

### 3.3 `--verbosity` 的统一实现

| 级 | filter | stats | index |
|---|---|---|---|
| 0 | 无 stderr 输出 | 同 | 同 |
| 1（默认） | 一行摘要 + 漏斗表 | 一行摘要 | 一行摘要 |
| 2 | + `--timing` 内容到 stderr | 同 | 同 |
| 3 | + `RUST_LOG=debug` 的 tracing | 同 | 同 |

实现位置：`reed-cli/src/translate/*.rs` 统一一个 `Verbosity` 枚举；驱动不打印，返回报告，由 CLI 决定打印什么（现 `reed/src/filter/report.rs` 的 `eprint_report` 应移到 CLI）。

### 3.4 reader 吞吐（BCF stats 的门仍未达）

doc_v3 STATUS 的分析成立：reader 单线程 3.4–3.7 GB/s，64 KiB 粒度 `read` 的 syscall + 拷贝是主项。建议按顺序试：
1. `posix_fadvise(SEQUENTIAL)` + 把 `BufReader` 容量提到 8 MiB（零风险，先测）；
2. reader 用 `read_exact` 一次读入整个 chunk 的压缩字节（当前按块读；块头已知 BSIZE，可以先读 1 MiB 再在内存里切块头，只有最后一块需要补读）；
3. mmap 输入文件，chunk 持 `Arc<Mmap>` + 区间（doc_v3 O11）。
每步以 `--timing` 的 reader 累计（A4-5）为门。

### 3.5 长度前缀链的看门狗（01 F4-14）

把 `STALL_LIMIT` 从"单 cell 等待 120 s"改为"链上任何发布事件的全局静默 ≥ 120 s"，并在错误中带上等待的 chunk id、已等待时长、最后一次发布的 chunk id。NFS 卡顿时不再误判为死锁。

### 3.6 文档漂移的防范

- `doc_vN/STATUS.md` 每一行的"验证"列必须引用一个测试名或一条可复现命令；文档声明"未做"的项，CI 用一条 grep 守卫（例如 `grep -q noodles crates/reed-filter/Cargo.toml && echo "C3 still pending"`），避免 F4-2 这类不一致。
- `--help` 文本来自 clap 结构体的 doc 注释：为 `<INPUT>` 等共用项建一个常量，三处共用（F4-13）。
