# doc_v6 实施状态

> 对照 `06-roadmap.md` 的 Y0–Y7。每项注明 commit 与验证方法（测试名或可复现命令）；偏差单独登记。格式沿用 `docs/rounds/v5/STATUS.md`。

## 已完成

（实现者填写）

## 未完成 / 偏差

| 项 | 状态 | 说明 |
|---|---|---|
| Y0–Y7 | 待做 | 见 `06-roadmap.md` |

## 基准（2026-10-10 v6 复审基线，56 线程，热缓存，批次开始负载 1.28；见 `04-performance.md` §1）

| 命令 | wall | user |
|---|---|---|
| filter gold `-Oz6`，VCF 输入 | 8.85 s | 460 s |
| filter gold `-Oz6`，BCF 输入 | 9.00 s | 464 s |
| filter gold `-Ob` | 9.98 s | 521 s |
| filter gold `-Oz1` | 5.64 s | 289 s |
| stats 五统计 VCF / BCF | 1.78 / 1.65 s | 72 / 56 s |
| reports 两族 + summary | 3.0–5.4 s（双峰） | 147 s |
| stats 边过滤 `--pi` | 3.09 s | 153.6 s |
| index VCF / BCF | 1.63 / 1.12 s | 75 / 44.5 s |
| `--sfs` / `--fst-global` / f 统计 | 1.60 / 1.52 / 1.56 s | 73 / 71 / 71.5 s |
| `--ld-decay 100k` | 169 s | 343 s |
| `filter --indep-pairwise 50 5 0.2 -Oz1` | 66.6 s | 289 s |
| plink2：pgen 导入 / 从 pgen 剪枝 / 从 pgen 报告 | 27.8 / 0.66 / 1.48 s | 105 / 1.5 / 2.75 s |
