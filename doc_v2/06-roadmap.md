# 06 修复与重构路线

> 原则：**先修正确性，再拆结构，最后加功能**。
> 本数据集上性能已经超额完成（stats 2.65 s，filter 14.2 s），后续每一步的性能门只要求"不回退超过 10%"。

---

## 1. 通用验收门（每一步都必须满足）【规】

1. `cargo test --release --workspace` 全绿；`cargo clippy --workspace --all-targets` 不新增警告。
2. 每个 P0 缺陷**先写出复现它的失败测试，再修复**（红 → 绿），测试名描述缺陷本身。
3. 过滤对拍：全量 `bcftools query` 的 GT 与 AC/AN/NS md5 与 gold 一致（01-review §1.2 的命令写成 `xtask parity`）。
4. 确定性：threads ∈ {1, 2, 3, 7, 56} 输出 sha256 相同（filter 与 stats 都要测）。
5. 性能：全量 stats ≤ 3.0 s（3 群体 π/dxy/Fst/θ/Tajima's D）、filter ≤ 15.6 s（当前值 +10%），并在 `bench/e2e/results.tsv` 追加一行。
6. **对齐**：凡涉及统计数值的改动，`gold/popgen` 对拍（04-algorithms.md §5）必须通过。
7. 同一 PR 内删除被替代的旧实现，不留下"以后再删"的并行代码。
8. doc_v2 中相应条目标记为已完成，并注明 commit。

---

## 2. 步骤

### R0 立即修复（不改结构，最小改动）

| 项 | 内容 | 回归测试 |
|---|---|---|
| D1 | 真正的票据：reader 取令牌，reducer 按序交付后归还 | 现有失败用例转绿；新增"chunk 0 延迟 500 ms"测试 |
| D2 | BGZF 空块跳过，只在底层 EOF 结束 | 拼接文件 fixture；在每个块偏移处插入 EOF 块的穷举测试 |
| D11 | 块按 65,280 字节切，记录可以跨块；`deflate_block` 超限报错 | 合成 3,000 样本（单行 > 100 KB）输出，`bgzip -t` + `bcftools view` 读回 |
| D12 | `-o x.bcf` 显式报错（BCF 输出在 R6 实现） | CLI 测试 |
| D4 | tid 用 u32，未知 contig 报错；`--write-index` 默认 CSI，TBI 遇到 ≥ 2^29 的位置报错 | 70,000 个 contig 的合成文件；位置 6e8 的合成记录用 CSI 查询 |
| D14 | `parse_allele` checked 算术 | 单元测试 |
| D19 | 补登 M4/M5/M5b 的性能行 | — |

### R1 删除 legacy，统一 header 与输入入口

- 按 03 §2 删除；legacy 测试中有价值的用例迁移到 `reed-stats/tests/`。
- `HeaderBundle` + `reed::input::open`；驱动不再调用 `record_header`/`filter_header`。
- stdin（`-`）对 BGZF VCF 可用。
- 门：§1 全部；工作区代码量预计减少约 6,000 行。

### R2 统一记录解码（03 §3）

- `GtCode`、`GenotypeScanner`、`FormatLayoutCache`、`SiteView`；stats、filter（单遍/sketch/reduce）、qc 全部迁移。
- 策略在扫描器中执行；filter 报告的 `policy_counters` 开始有值（修复 D13）。
- 门：filter 输出与 R1 结束时 sha256 相同；stats TSV 与 R1 结束时相同（文件布局变化前）；差分测试对 noodles 通过；fuzz 60 s 无失败。

### R3 引擎

- 按字节的票据预算、单一代码路径、去掉 watcher 线程、类型化 `EngineError`（03 §4）。
- 门：§1；内存上界测试（在途字节 ≤ 预算）。

### R4 输出层

- `BlockWriter`/`OrderedOutput`/索引写出统一（03 §5）；filter 两个 job 的重复代码删除；`reed/src/filter.rs` 按 03 §8 拆分。
- 门：§1；输出 sha256 与 R3 相同（块边界改变时除外：块切分规则从"按记录"改为"按 65,280 字节"会改变压缩字节，
  此时以解压后正文 md5 相同为门）。

### R5 CLI 重写（02-cli.md）

- 新建 `reed-cli` crate；clap 参数组；单位解析；表达式解析器与校验（`reed-filter::expr`）；`FilterSpec::to_args` 往返测试；
  `##reed_filterCommand`；删除 TOML/conf 路径与 `toml` 依赖。
- 命令：`reed filter`、`reed stats`（按科学方法命名的统计参数、`--variant-only`、写出器，02-cli.md §5）、`reed index`；reed-qc 并入 reed-stats 的 `summary/`，成为 `--summary`。
- `gold/scripts/run_reed.sh` 翻译 `filters.conf`；`compare_tools.py` 增加 reed 执行器。
- 门：§1；gold 等价命令（02-cli.md §6）输出的 GT/标签与 gold 一致；`reed stats --summary` 与 `bcftools stats` 九节逐行一致；
  表达式解析器对 bcftools 手册中的例子逐条对拍（同一输入下 `bcftools view -i/-e` 与 `reed filter -i/-e` 选出的位点集合相同）。

### R6 输入/输出补全（03 §6）

- 未压缩 VCF、普通 gzip、BCF 输入；`-r/-R` 索引区域；BCF 输出（`-Ob/-Ou`）。
- 门：同一数据以四种输入格式运行 filter/stats，结果一致；BCF 输出经 `bcftools view` 读回后与 VCF 输出一致。

### R7 统计对齐（04-algorithms.md）

- 先做输入契约（D7 contig 长度不编造、D8 位置升序、D9 pixy 群体文件读法）——属于正确性，可提前到 R2 之后。
- 搭建 `gold/popgen`：在合成 fixture 与全量数据上运行 pixy 2.2.4、vcftools、scikit-allel、plink2、easySFS、bcftools stats，产出期望输出。
- 按 04 §3 的唯一定义重写各统计量（经典方法 + pixy 的缺失数据处理 + `--variant-only` 分母补全）；撤下 complete-case、清单门、reed 自定义网格、PBS 交集语义（04 §1）；
  `reed-core/src/manifest.rs` 与 `##reed_ascertainment` 删除。
- 实现 pixy 的不变位点检查与 `--bypass-invariant-check`。
- 门：04 §5 对拍表的每一行都有用例并通过；特别是"无缺失数据时 vcftools、scikit-allel、pixy、reed 四者一致"的用例必须存在；
  全量上 `--variant-only --pi` 对 vcftools `--window-pi`、`--freq/--hardy/--het/--missing` 对 plink2 通过。

### R8 更多参考实现

只增加"已有主流实现可对拍"的统计量：LD（plink2 `--r2-unphased`/`--indep-pairwise`、vcftools `--geno-r2`）、
f 统计（Patterson et al. 2012，参照 scikit-allel `admixture.py`）、联合 SFS（Gutenkunst et al. 2009，参照 easySFS）。
HWE 记忆化、位平面 LD 等属于实现层优化，可以做，但数值必须与参照物一致。参数仍按方法命名。

---

## 3. 顺序与依赖

```
R0 ──► R1 ──► R2 ──► R3 ──► R4 ──► R5 ──► R6
              │                     │
              └──► R7（输入契约先行；对齐可与 R3–R4 并行）
                                    └──► R8
```

输入契约（D7、D8、D9）属于正确性 P0，建议在 R2 完成后立即插入；统计对齐需要 R5 的参数与写出器，可以先在库层实现并对拍。

---

## 4. 流程纪律（复审中暴露的问题）【规】

1. **测试红着不提交**。`pipeline.rs` 的确定性测试是高频偶发失败（复审时 5 次运行失败 4 次），却没有阻止后续里程碑的提交。
   对并发测试，CI/本地脚本应该重复运行 20 次，以暴露偶发失败。
2. **每个里程碑都入基准日志**（v1 §17.3）。
3. **"invariant"只能用于程序内部不变量**。由输入数据决定的条件（块大小、contig 名、位置范围、样本数）一律返回错误，
   不能 `expect`/`assert`/`as` 截断。clippy 的 `cast_possible_truncation` 警告必须逐条处理。
4. **不新增统计语义**：任何新统计量、新估计器、新缺失处理方式，都必须先指明它复现的主流工具与源码位置；指不出来就不做。
5. **新实现落地时删除旧实现**，不允许出现第二份 GT 解析器、第二份输出写出器。
6. **规范与实现的偏差必须记录**：实现者选择偏离 spec 时（例如 stats 驱动不用 `GroupTable`），在 doc_v2 中登记并说明理由，
   而不是悄悄绕过。
