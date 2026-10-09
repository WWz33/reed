# 01 第五轮复审（doc_v4 W0–W4 落地后，2026-10-09 晚）

> 复审对象：`61d1a14`（doc_v4）之后的 10 个提交，到 `c5b6bfe`（W5 前置调查）；以及工作树的未提交改动。
> 方法：对照 `doc_v4/STATUS.md` 逐条核实（代码事实、grep 守卫、`--help`）；工作树与 HEAD 分别构建测试（HEAD 在独立 worktree，同一 target 配置）；全量基准 3 次取中位数（04）；跨格式 md5。本机无 `perf`、无 `/usr/bin/time`。

---

## 0. 总评

doc_v4 的 W0–W4 全部落地且质量与前几轮一致：每项有 commit、测试名、数字；两条用户裁决（`--geno/--mind` 保留、`--fst-global`）在代码里的部分已经执行到位（W0b），守卫脚本 8 条全 ok。性能上 W3（mmap）与 W4（BCF 原生过滤输入）兑现了 doc_v4 的全部目标：stats 五统计 1.7 s / 1.5 s（BCF ≤ VCF），BCF 输入过滤 8.9 s（与 VCF 持平），index 也顺带快了三分之一。

本轮的问题分三类：

1. **主干卫生**：复审当时工作树**不能编译**（W5 在途改动未提交，§1）；STATUS 关于 `unsafe` 的陈述与代码不符（两处而非一处）。
2. **设计待定的功能**：W5 的统计出口在 doc_v4 只写了累加器没定形状，实现者按窗口统计的模子开了头，与 SFS 是 genome-wide 量的事实冲突（F5-3）；本轮在 05 把 SFS / `--fst-global` / f 统计 / LD 的出口形状全部定下来。
3. **历史残留**：`crates/` 内两代并存与名字漂移（03 §4 C-1…C-9）、目录树（06）、doc_v2 承诺未兑现的三个参数（02 §11）。

---

## 1. 状态核实

| 项 | 结果 |
|---|---|
| 工作树 | `git status`：`crates/reed-stats/src/{request,sfs}.rs` 有未提交的 W5 改动（`StatKind::ProjectedSfs`、`bins` 扩到 `n+1`）；本轮文档改动在 `doc_v4/*`、`docs/INDEX.md`；`doc_v4/07-layout.md`、`doc_v5/` 未跟踪 |
| 工作树构建 | `cargo build --release --workspace --all-targets` **失败**：`reed-cli` bin `error[E0382]`（use of moved value）；`cargo test --workspace` **编译失败**：`reed-stats` 测试 `accumulator_laws.rs:60-118` 的 `match` 未覆盖 `StatAcc::ProjectedSfs`（`E0004`） |
| HEAD `c5b6bfe`（独立 worktree `/tmp/v5_head`，`CARGO_TARGET_DIR=/tmp/reed-target-head`，已清理） | `cargo test --release --workspace` **364 通过 / 0 失败**；`cargo build --release --workspace --all-targets` **0 警告** → STATUS "最终验证"在 HEAD 上成立 |
| clippy（工作树可编译部分） | 0 error；pedantic 提示 ≈ 112 条（`assert_eq` 22、`u64→u32` 截断 15、`format!` 追加 14、naive bytecount 11、长字面量 11…）沿 doc_v3 C7 豁免 |
| `bash scripts/check-doc-drift.sh` | 8 条全 ok（C3 无 noodles、无 reed-qc、无版本后缀类型、`--geno/--mind` 真实、`--max-missing` 只在报错、漏斗 id、index 拒 stdin、`Verbosity` 单一） |
| `--help` 核对 | filter / stats / index 三份与 02-cli.md 逐项一致，差异登记在 02 §11（D-1…D-6） |
| 规模 | 8 crate，src ≈ 32.9k 行（reed 7.5k、io 7.0k、filter 4.2k、stats 4.0k、core 4.0k、oracle 2.9k、cli 1.9k、engine 1.5k） |

---

## 2. STATUS 自述 vs 事实

| STATUS 条目 | 核实 | 结论 |
|---|---|---|
| W0b：`--geno/--mind` 恢复、`--max-missing*` 报错、`MissingRate`、漏斗 id `after_mask.missing`/`sample.mind` | `filter_flags.rs:144-151`、`cli.rs:81-86`、`funnel.rs:11,16`、守卫 | ✓ |
| W1：`--verbosity`、index 拒 stdin、`INPUT_HELP`、去 `-g`、`FilterSpec`、`--no-command`、报告模式 `-t/-T` | `--help` 三份；`cli_contract.rs` | ✓ |
| W2：reports 一次读取 3.58 s；边过滤 `--pi` 2.99 s（门 ≤ 3.0） | 本轮 3.76 s / 3.04 s（04 §1） | 一次读取 ✓；**reports 门未达**（STATUS 已登记）；`--pi` 门在边缘（F5-11） |
| W3：mmap reader；stats 1.71 / 1.47 s；BCF ≤ VCF | 1.74 / 1.51 s；index 也从 2.4/2.2 → 1.60/1.05 s（STATUS 未写） | ✓，且收益比记录的更广 |
| W3 偏差："唯一一处 `#[allow(unsafe_code)]`" | `grep -rn 'allow(unsafe_code)'`：**两处**——`reed-core/src/count/kernel.rs:16`（模块级，AVX-512 内核，6 个 `unsafe`，doc_v1 §14.4 的批准例外）与 `reed-io/src/chunk/bgzf.rs:361`（mmap） | ✗ 陈述不实（F5-2） |
| W4：BCF 原生过滤输入；全量 9.06 s；user −26% | 8.87 s / 465 s；正文 md5 与 VCF 输入相同（STATUS 记录，本轮未重跑 gold 全量写出） | ✓ |
| W5："未做，前置调查完成" | 工作树有在途代码（未提交、未入 STATUS），且测试不编译 | 状态落后于代码（F5-1、F5-3） |
| gold 全量 md5 `7aedc59ac798…`、`check_parity` PASS | 本轮**未重跑**（需写出 ≈ 2 GB 并跑 3 家对拍，以 STATUS 的 W4 重跑为准） | 未独立核实，标注 |

---

## 3. 正确性核实（本轮）

| 检查 | 结果 |
|---|---|
| stats 五统计 VCF vs BCF 输入 | `_pi.txt` md5 `45c2f2e1632076fa601eb5575d6288f3`（与 doc_v4 记录相同）、`_fst.txt` `fa144045…`，两种输入相同 |
| 三策略 / 跨格式 / 跨 reader 门 | 由 HEAD 的 364 项测试覆盖（`e2e_filter`、`e2e_bcf`、`determinism`、`e2e_missing_boundary` 等 9 个 e2e 文件） |
| 子集数字门 263,420 / 381,672 | STATUS W0b 记录；`e2e_missing_boundary` 在 HEAD 通过 |
| 守卫 | 8/8 ok |

---

## 4. 性能

见 04-performance.md。要点：BCF 输入过滤 12.1 → 8.87 s；stats 五统计 2.7/3.1 → 1.74/1.51 s；index 2.4/2.2 → 1.60/1.05 s；reports 5.3 → 3.76 s（门 3.2 未达）；filter VCF `-Oz6` 8.79 s 不回退。测量纪律修正见 04 §4。

---

## 5. 发现

| # | 级 | 位置 | 问题 | 证据 | 修法 |
|---|---|---|---|---|---|
| F5-1 | **P0（主干）** | 工作树 `crates/reed-stats/src/{request,sfs}.rs`（未提交） | `reed-cli` bin 与 `reed-stats` 测试都不编译；主干处于"不能测"的状态 | §1 两条错误 | 07 X0：补齐 `accumulator_laws.rs` 的 match 臂并按 05 §3.1 完成后提交，或 `git stash`；CI 第一条就是 `cargo test` 编译通过（X7） |
| F5-2 | P2（文档/策略） | `doc_v4/STATUS.md` W3 偏差行；`bgzf.rs:356` 注释 | "唯一 unsafe" 不实：还有 `count/kernel.rs` 的模块级允许 | §2 | 文本改为"两处"；守卫断言 `allow(unsafe_code)` 恰好 2 处（03 §2） |
| F5-3 | P2（设计） | 在途 `StatKind::ProjectedSfs { group, n, folded }` | SFS 是 genome-wide 一个谱，被放进窗口统计的 `StatKind`/分箱模型；`--levels window` 会产出无意义的分窗谱，且 dadi/fsc2 写出需要收尾合并而非逐窗写行 | diff | 05 §3.1：`StatsRequest.sfs` 独立 genome-wide 路径，partial 并列上交、sink 收尾写出；保留 `bins = n+1` 的改动 |
| F5-4 | P2（性能） | `reed/src/stats/writers/qc_report.rs:153` | `--summary` 自建 `tracking_alleles()` 扫描器重扫 GT，reports 两族 3.76 s / 186 CPU-s，门 ≤ 3.2 s 未达 | 04 §1 | 03 §5：summary 消费 `ScanOutput`（X3） |
| F5-5 | P3（文本漂移） | `reed-filter/src/funnel.rs:67,72`、`reed-stats/src/lib.rs:47`、stats `--tmp-dir` 帮助 | 注释仍写 `sample.max_missing`、`--fst-pops`、"report runs that read twice" | grep | X1 清理；守卫加禁词 |
| F5-6 | P3（残留） | `reed/src/stats/writers/qc.rs`、`reed-cli` `QcCli`/`translate/qc.rs`/`Command::Qc` | M3 死代码与隐藏的 deprecated `reed qc` | 03 §4 C-1/C-2 | X1 删除 |
| F5-7 | P3（残留） | `reed-io/src/bcf/encode.rs`、`reed/src/filter/single_pass.rs:31` | noodles 编码器仍在生产路径作可选项 | 03 §4 C-4 | 单独 PR 移出生产路径 |
| F5-8 | P3（规范缺口） | `reed stats --help` | doc_v2 承诺的 `--bed-file`、`--sites-file`、`-q/-Q :type` 未实现 | 02 §11 D-1/D-2 | X6 |
| F5-9 | P3（测量） | 本轮基准脚本 | 逐次采样的 `load_1m` 被自己的 56 线程运行抬到 9–42，"负载 > 1 不入日志"会剔除全部样本 | 04 §4 | 批次开始采样 + 3 次极差 ≤ 3% 判据 |
| F5-10 | P3（结构） | `crates/`、顶层目录 | C-3/5/6/7/8/9 与 06 的目录树都未执行 | 03 §4、06 | X1 / X2 |
| F5-11 | P3（门） | W2 边过滤 `--pi` 门 ≤ 3.0 s | 本轮 3.045 / 3.042 / 3.005 s，STATUS 记 2.99 "达标"；实际在门边缘 | 04 §1 | 门改写为 CPU 口径（≤ 155 → 130 s，P5-2），wall 作参考 |

---

## 6. 值得保留的做法

- `--timing` 的 reader/reduce/chain_wait/worker_busy 分解：本轮所有"瓶颈在哪"的判断都来自它。
- 每条 STATUS 行引用测试名或命令：本轮能逐条核实，正是因为这条纪律；唯一不实的陈述（unsafe）恰好是没有引用守卫的那条。
- 用户裁决先落文档再落代码：W0 → W0b 的往返只花了一个提交，因为语义层（精确判定、等位计）与拼写层是分开的。
- HEAD 在独立 worktree 验证：工作树不干净时仍能给出可信的"主干状态"。
