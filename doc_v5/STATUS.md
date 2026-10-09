# doc_v5 实施状态

> 对照 `07-roadmap.md` 的 X0–X7。每项注明 commit 与验证方法（测试名或可复现命令）；偏差单独登记。格式沿用 `doc_v4/STATUS.md`。

## 已完成

### X0 解除工作树阻塞（07 X0；01 F5-1 / 03 A5-4）

复审当时工作树有两处编译错误（`reed-cli` bin `E0382`、`reed-stats` 测试 `E0004`）。处置：在途 W5 代码是**窗口式** `StatKind::ProjectedSfs`，与 05 §3.1 的 genome-wide 路径冲突（01 F5-3），且其 CLI 默认（folded）、输出路径（`PREFIX.sfs`）、文件格式均与 05 §3.1 不符；按 07 X0 认可的"`git stash` 保主干干净"一支，把在途代码整体 stash，留待 X4 按设计重做。

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| 在途 W5 代码 stash | `git stash push -u -m "w5-inflight-window-sfs (superseded by 05 §3.1 genome-wide path; to rework in X4)" -- crates Cargo.lock`（stash@{0}）；备份 `/tmp/v5_w5_inflight_tracked.patch`、`/tmp/v5_w5_writers_sfs.rs.bak` | 本提交（仅文档） | `git status`：`crates/`、`Cargo.lock` 已干净，仅余文档改动 |
| 主干可测 | 主干回到 `c5b6bfe` 的代码 | 本提交 | `cargo test --release --workspace` **364 passed / 0 failed**；`cargo build --release --workspace --all-targets` **0 警告** |
| 守卫基线 | — | 本提交 | `bash scripts/check-doc-drift.sh` **8/8 ok** |

## 未完成 / 偏差

| 项 | 状态 | 说明 |
|---|---|---|
| X1 | 待做 | crate 与文档卫生（C-1/2/3/8/9、A5-2/3、unsafe 守卫） |
| X2 | 待做 | 目录搬迁（06） |
| X3 | 待做 | reports 门（03 §5） |
| X4 | 待做 | 统计接入（SFS → `--fst-global` → f 统计 → LD；含 C-5）；SFS 按 05 §3.1 的 genome-wide 路径重做（stash@{0} 可参考，但其窗口形态不采用） |
| X5 | 待做 | 多等位重算（05 §4） |
| X6 | 待做 | CLI 缺口（`--bed-file`、`--sites-file`、`-q :type`） |
| X7 | 待做 | CI（`.github/workflows` + `ci/` 守卫） |

## 基准（2026-10-09 v5 复审基线，56 线程，热缓存；负载与中位数见 `04-performance.md` §1）

（复审时的数字由 `04-performance.md` 给出；实现者在每步完成后按同一命令重测并补登 `bench/e2e/results.tsv`）
