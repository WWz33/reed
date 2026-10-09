# 07 路线（X0–X7）

> 通用门沿用各轮：测试全绿、gold 对拍（正文 md5 `7aedc59ac798729a2674e6e7ef19adfc` + `check_parity` PASS）、三策略与跨格式 md5、性能不回退并入 `bench/.../results.tsv`（含负载与中位数）、同 PR 删除旧实现、更新 `docs/rounds/v5/STATUS.md`（每行引用测试名或可复现命令）。

```
X0 解除工作树阻塞 ──► X1 crate 与文档卫生（C-1/2/3/8/9、A5-2/3、unsafe 守卫）──► X2 目录搬迁（06）
                                                                                   ├─► X3 reports 门（03 §5）
X4 统计接入（SFS → --fst-global → f 统计 → LD；含 C-5）可与 X1–X3 并行
X5 多等位重算（05 §4）
X6 CLI 缺口（--bed-file、--sites-file、-q :type）
X7 CI（.github/workflows + ci/ 守卫）——X2 之后立刻
C-4（noodles 编码器退出生产路径）与 C-7（job.rs 拆分）各自单独 PR，任何时间
```

| 步 | 内容 | 门 |
|---|---|---|
| **X0** | 工作树当前 `cargo test --workspace` 编译失败（`accumulator_laws.rs` 未覆盖在途的 `StatAcc::ProjectedSfs`）。要么补齐 match 臂并按 05 §3.1 把 SFS 改成 genome-wide 请求路径后提交，要么 `git stash` 保主干干净 | `cargo test --release --workspace` 全绿；`cargo build --release --workspace --all-targets` 0 警告 |
| **X1** | 03 §4 的 C-1（删 M3 `qc.rs`）、C-2（删 `reed qc`）、C-3（`qc_report.rs` → `summary.rs`、`.qc.json` → `.summary.json`）、C-8（reed-cli 未用依赖）、C-9（`crates/oracle`）；A5-3 注释清理（`sample.max_missing`、`--fst-pops`、stats `--tmp-dir` 帮助）；A5-2 `unsafe` 陈述改为两处 + 守卫；**08 §4 的 L-1…L-7**（README 重写为 v0.3 现实、门面公开面收口、去 noodles re-export、8 个 crate `publish = false`、description 更新、examples 移走、嵌入 API doctest）；**08 §6 的前后端边界**：`reed-cli` 只依赖 `reed`（门面 re-export `spec / expr / budget / report`）、后端 5 处 `eprintln!` 改 `tracing::warn!`、`plan::compile` 与 `assign_groups` 移入 `run_*`、单一 `reed::Error` 带 `kind()` 供退出码 | 测试全绿；全部输出 md5 不变；`bash ci/check-doc-drift.sh` 含新守卫全 ok（含 08 §6.4 四条）；`cargo doc --no-deps -p reed` 的公开项与 08 §6.2 一致；`cargo test --doc -p reed` ≥ 3 |
| **X2** | 06-layout.md §4 的 `git mv` 全部执行：`docs/{spec,rounds}`、`ci/`、`bench/suites/`、`gold/{work,out}` → `target/gold/`（`REED_SCRATCH`）、`docs/rounds/v1/design.md` → `docs/rounds/v1/`；§5 的引用一次 sed；`docs/spec/cli.md` = docs/rounds/v5/02 | `cargo test` 全绿；`gold/scripts/run_all.sh` 在新路径跑通并 md5 不变；`! grep -rq 'doc_v[0-9]/' docs/README.md` |
| **X3** | 03 §5：reports job 一开始就用 `tracking_alleles()` 扫描器，summary 累加器消费 `ScanOutput`，删自建扫描 | 8 个报告文件 md5 不变；全量 ≤ 3.2 s、CPU ≤ 155 s（负载 < 1，3 次中位数） |
| **X4** | 05 §3 按序：3.1 SFS（genome-wide 路径、dadi/fsc2 写出、`run_easysfs.sh`）→ 3.2 `--fst-global`（`StatKind::FstGlobal`、`_fst_global.txt`、`run_vcftools_fst.sh`；作废 `--fst-pops` 的注释与 doc_v2 行）→ 3.3 f 统计（`.fstats`、block jackknife、`compare_fstats.py`）→ 3.4 LD（尾窗携带、`LdPruner` 按 contig 串行、`run_plink2_ld.sh`）；C-5 随 SFS 一起拆 `individual.rs` | 每项：fixture 逐值门 + 全量门（05 §3 各小节）；通过后去掉 `#[doc(hidden)]`、`--help` 列出、`cli_contract` 快照更新 |
| **X5** | 05 §4 多等位逐等位重算 | `multiallelic.vcf` 与 `bcftools +fill-tags` 逐字节一致；gold md5 不变 |
| **X6** | `--bed-file`（BED 窗口，等价于 pixy `--bed_file`）、`--sites-file`、`-q/-Q` 的 `:type` 后缀 | pixy BED 模式逐列一致；bcftools `-q 0.05:minor` 位点集合一致 |
| **X7** | `.github/workflows/ci.yml`：`cargo test --release --workspace`、`clippy`、`bash ci/check-doc-drift.sh`、子集数字门（263,420 / 381,672） | CI 绿；守卫脚本在 CI 失败即阻断 |

X0 先于一切，当天完成；X1 与 X2 合成一个"卫生"PR；X4 的每个统计量单独 PR。
