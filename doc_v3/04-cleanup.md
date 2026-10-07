# 04 语义修正、结构清理与 CLI 一致性

> 本文列出不涉及 BCF、也不以性能为目的的修正项。C 开头的条目各自独立，可以穿插在 BCF 工作之间完成。

---

## 1. 语义与结构

### C0 `--mind` 的实现与声明的语义不一致【P1，对齐】

- **现状**：CLI 帮助与 doc_v2 02 §3.1 都写 "plink2 语义：掩码之后计算，分母是载入期过滤后的位点"；
  但实现（`reed/src/filter/reduce.rs`、`sketch.rs` 的 `cohort_miss`）统计的是**全部记录上的原始（掩码前）缺失**，
  包括头部规则淘汰的记录——这是 gold 管线的口径，不是 plink2 的。
- **plink2 的实际行为**（源码核实）：VCF 导入时 `--vcf-min-dp/--vcf-min-gq` 先掩码；`MindFilter`（`plink2.cc:2061`）在
  `EnforceGenoThresh`（`plink2.cc:2707`）之前执行，分母为载入期过滤（区域、`--snps-only`、`--max-alleles`、`--var-filter` 等）后的全部位点。
- **修法**：`--mind` 按 plink2 计算——统计**通过头部规则（与基因型无关的位点规则）的记录**上、**掩码之后**的样本缺失；
  在尾部规则（`--geno`、`--hwe`、`--maf`…）之前生效。三种执行策略同时修改并保持输出一致。
  gold 对拍继续用 `-S samples.keep.txt`（gold 的原始口径），不受影响。
- **门**：与 plink2 在同一 fixture 上 `--vcf-min-gq 20 --mind 0.1` 保留的样本集合一致；三策略输出一致。

### C0b 样本规则互斥导致参数被静默忽略【P1】

- **现状**：`reed-cli/src/translate/filter.rs:23-36` 用 `if / else if` 链把 `-s`、`-S/--keep`、`--remove`、`--mind` 变成**一个** `CohortRule`
  （`spec.cohort: Option<CohortRule>`）。同时给出 `-S keep.txt --mind 0.2` 时 `--mind` 被忽略；`-S` 与 `--keep` 同时给出时 `--keep` 被忽略；
  没有任何提示。
- **修法**：`cohort` 改为两部分——静态子集 `sample_select: Option<SampleSelect>`（`-s/-S/--keep/--remove` 的组合，按 bcftools/plink2 规则求交与差）
  与动态规则 `mind: Option<f64>`；两者可以同时存在，顺序是先静态子集、后 `--mind`（plink2 同序）。
  同类参数互相冲突时（例如 `-s` 与 `-S` 同时给出）报错，而不是择一。`-s ^A,B` 的排除语法照 bcftools 实现（现在只有区域支持 `^`）。
- **门**：组合用例的保留样本集合与 bcftools（静态部分）、plink2（`--keep` + `--mind`）一致；`to_args()` 往返测试覆盖组合。

### C1 删除清单 v2 残留（F6）

删除 `reed-filter/src/plan/ascertainment.rs`、`provenance.rs` 中的 `reed_ascertainment` 行、`reed-core/src/manifest.rs`，
修改 `reed/tests/e2e_filter.rs:325` 的断言（doc_v2 04 §1 的裁决）。

### C2 reed-stats 的公开面只保留"有 CLI 出口且有参照物对拍"的方法（F7）

| 现导出 | 处理 |
|---|---|
| `NeiMultiAcc`、`HudsonFiniteAcc` | 删除（doc_v2 04 §1 已撤下） |
| `FstMatrixAcc`、`WcMultiAcc` | 保留但 `pub(crate)`，直到 `--fst-pops` 接入并与 vcftools 多群体 WC 对拍 |
| `LdDecayAcc`、`LdPruner`、`pair_sums` | `pub(crate)`，直到 R8 的 LD 接入并与 plink2 `--r2-unphased`/`--indep-pairwise` 对拍 |
| `F2Acc/F3Acc/F4Acc`、`f4_ratio`、`BlockJackknife` | `pub(crate)`，直到 R8 接入并与 scikit-allel `admixture.py` 对拍 |
| `ProjectedSfsAcc`、`JointSfsAcc` | `pub(crate)`，直到 `--sfs --proj` 接入并与 easySFS 对拍 |
| `IndividualFAcc` | 已被 `--het` 使用，保留 |

原则：库的公开 API 不出现没有参照物的统计方法（"不做算法创新"的 API 投影）。

### C3 reed-filter 去掉 noodles（F8）

在 reed-core 定义 `TagTable`（`TagId`、名字、`Number`、`Type`），由 reed-io 从 header 构建；`plan::compile` 只接收 `TagTable`。
完成后 reed-filter 只依赖 reed-core（doc_v2 03 §7），表达式与计划可以脱离 I/O 单测。

### C4 合并 QC 入口（F9）

- 只保留 `reed stats --summary`；`reed qc` 在一个版本内作为隐藏别名并打印弃用警告，之后删除。
- reed-qc crate 并入 reed-stats 的 `summary/` 模块（doc_v2 03 §8），`reed/src/{qc.rs,qc_report.rs,reports.rs,pixy.rs,stats.rs}`
  按 doc_v2 03 §8 归入 `reed/src/stats/`（`job.rs`、`sink.rs`、`writers/*.rs`）。

### C5 类型化 `EngineError`（F11）

按 doc_v2 03 §4.4 落地（`Job(Box<dyn Error + Send + Sync>)` 保留 source 链），在 BCF 工作之前完成（02-bcf.md B0）。

### C6 测试用临时文件

本机 `/tmp` 中残留的 `reed-stdin-*.vcf` 来自测试运行。O5（03-performance.md）修复后，测试中再加一条断言：进程结束后临时目录中没有 reed 的文件。

### C7 clippy pedantic

约 700 条警告中，`expect()`（404 条）应逐条区分：内部不变量保留（消息以 `invariant:` 开头），输入相关条件改为返回错误
（doc_v2 06 §4.3）；`cast_possible_truncation` 类逐条处理。其余风格类可以在 `clippy.toml` 中按项豁免，避免噪声淹没有用的警告。

---

## 2. CLI 一致性（对照 doc_v2 02-cli.md）

### 2.1 需要修正

| 现状 | 规范 | 处理 |
|---|---|---|
| stats 的群体文件主名是 `-g/--groups`，`--pop/--populations` 只是隐藏别名 | `--pop`（`--populations` 别名） | 主名改为 `--pop`；`-g/--groups` 改为隐藏别名并在一个版本后删除 |
| stats 不接受过滤参数与 `-r/-R/-t/-T/-s/-S/--keep/--remove` | 接受 §3 的全部过滤参数（INFO 更新类除外） | 把 filter 的参数组抽成共用的 clap `Args`（`filter_flags.rs`），stats 复用 |
| `--pbs` 是布尔值，固定取前三个群体 | `--pbs A,B,C`（可重复） | 改为接受群体名 |
| 缺 `--fst-pops`、`--sfs`、`--proj`、`--bed-file`、`--sites-file` | 有 | 随 R7/R8 接入；未接入前不在 `--help` 中出现 |
| `--contig-lengths-path` | `--contig-lengths` | 改名（旧名隐藏别名） |
| `--levels window,contig,genome` | 规范未列出 | 作为输出聚合参数保留，补进 02-cli 文档（它只是同一累加器的不同合并层级，不涉及统计方法） |
| `index -t/--type csi\|tbi` | bcftools：`-c/--csi`（默认）、`-t/--tbi`；`-t` 在 bcftools index 中就是"建 TBI" | 改为 `-c/--csi`、`-t/--tbi`、`-o/--output`、`-f/--force`、`--threads`（O4 之后有意义）；`-n/--nrecords`、`-s/--stats` 可选 |
| filter 的 `--compress-level` 与 `-Oz[0-9]` 重复 | 只用 `-O` | 删除 `--compress-level`（或作为隐藏别名） |
| filter/stats 的 `-o, --out` | `-o, --output`，`--out` 作为别名 | 主名改为 `--output` |
| `-O` 只接受 `v\|z` | `v\|z\|b\|u` | 随 BCF（02-bcf.md B2） |
| 缺 plink2 别名 `--vcf-min-dp/--vcf-max-dp/--vcf-min-gq` | 有 | 加别名 |
| 缺 `-q/-Q`（bcftools `--min-af/--max-af`） | 有 | 实现（语义同 bcftools：INFO 有 AC/AN 时取 INFO 值） |
| 缺 `--memory`、`--tmp-dir`、`--no-version`、`--verbosity`、`--timing` | 有 | 实现；`--tmp-dir` 与 O5 一起做，`--timing` 见 O10 |
| `qc` 子命令 | 并入 `stats --summary` | C4 |

### 2.2 过时的帮助文本

| 位置 | 现文本 | 应为 |
|---|---|---|
| filter/stats/qc/index `<INPUT>` | "Input BGZF VCF (.vcf.gz)" | "Input VCF/VCF.gz/BCF, or `-` for stdin"（BCF 在 B1 之前写明暂不支持） |
| filter `-o` | "Output VCF (.vcf.gz) — omit for plain VCF on the path/stdout" | "Output file (default: stdout); format from `-O` or the extension" |
| stats `-o` | "Write the TSV here instead of stdout" | "Output prefix: PREFIX_pi.txt, PREFIX.afreq, …" |
| `--threads` | "Worker threads (default: all cores)" | "Total threads (default: all cores)"（实现是总线程数） |
| `-O` | "b/u are refused until the BCF encoder lands" | B2 之后改为正式说明 |

### 2.3 规范的一处补充

doc_v2 02-cli §4 规定：表达式中的 `AC/AN/AF/MAC/MAF` 取到 INFO 中已有的值时打印警告。**实现中没有这条警告**
（`reed-filter/src/expr/*`、`reed-cli/src/translate/filter.rs` 中均无）。需要补上，因为它正是本数据集的陷阱（INFO/AC 陈旧），
警告应提示改用 `--maf/--mac`（按基因型重算）或先 `-x INFO/AC,INFO/AN`。
