# doc_v4 实施状态

> 对照 `06-roadmap.md` 的 W0–W6。每项注明 commit 与验证方法（测试名或可复现命令）；偏差单独登记。格式沿用 `doc_v3/STATUS.md`。

## 已完成

（实现者填写）

## 未完成 / 偏差

（实现者填写）

## 基准（2026-10-09 复审基线，56 线程，负载 0.3，热缓存）

| 命令 | wall | user |
|---|---|---|
| filter gold `-Oz6`（VCF 输入） | 8.8 s | 455 s |
| filter gold `-Ob` | 9.9 s | 519 s |
| filter gold，BCF 输入 | 12.1 s | 630 s |
| stats 五统计 VCF / BCF | 2.7 / 3.1 s | 71 / 56 s |
| stats 边过滤 `--pi` | 3.4 s | 170 s |
| reports 两族 | 5.3 s | 232 s |
| index VCF / BCF | 2.4 / 2.2 s | 77 / 46 s |
