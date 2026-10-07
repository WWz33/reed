# doc_v3 实施状态（2026-10-08）

> 对照 `05-roadmap.md` 的 V0–V6。每项注明 commit 与验证方法；偏差单独登记（格式沿用 doc_v2/STATUS.md）。
> 全量门：`cargo test --release --workspace` 全绿；`cargo clippy --release --workspace` 无 deny 级问题。

## 已完成

### V0 语义修正（C0、C0b）

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| C0 | `--mind` 改为 plink2 语义：统计**通过头部规则**的记录上、**掩码之后**的样本缺失；分母 = 头部存活的记录数。`MissingJob`（SinglePass/TwoPass 的 pass0）改走与单遍相同的 `RecordPipeline`（含掩码字段扫描），sketch pass0 同步改为 head 存活 + post-mask | 见下 | 与 plink2 对拍：fixture 上 `--vcf-min-gq 20 --mind 0.1` 两边都保留 {A,C}（B 掩码 2/10、D 缺失 3/10 被剔除，C 恰为 0.1 保留）；另有头部规则收缩分母的回归用例（s0 旧口径 91/100→剔除、新口径 1/10→保留），三策略输出一致 |
| C0b | `cohort: Option<CohortRule>` 拆为 `sample_select: Option<SampleSelect>`（`-s/-S/--keep/--remove`，keep ∩ ¬remove）与 `mind: Option<f64>`；先静态子集、后 `--mind`（plink2 同序）。`-s ^A,B` 排除语法；同极性参数（`-s`+`-S`、`-s`+`--keep`、`-s ^`+`--remove`、重复 `-s`）报错而非静默择一 | 见下 | 组合用例：`-S keep.txt --mind 0.2` 保留集合正确；`to_args()` 往返测试覆盖组合与冲突报错 |
| 防御 | `eval_head`/`eval_tail` 的计数器长度与规则数不一致时断言失败（此前 zip 静默跳过规则，导致 pass0 漏评估头部规则） | 见下 | 该断言在开发中直接定位了 pass0 funnel 未初始化的缺陷 |

开发中发现并修复的缺陷（非文档列出）：pass0 的 `FunnelPartial` 未按规则数初始化，`eval_head` 的 `zip` 空转使头部规则被跳过（SinglePass/TwoPass 的 `--mind` 因此按旧口径统计）；pass0 的扫描器未带计划 FORMAT 列，掩码未参与缺失计数。

### 待办

- V1 O1–O5、V2 B0、V3 B1、V4 B3a、V5 B2、V6 B3b；清理 C1–C4、C7 与 CLI 一致性（04 §2）。
