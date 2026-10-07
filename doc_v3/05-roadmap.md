# 05 路线（V0–V6）

> 延续 doc_v2/06-roadmap.md 的通用验收门（测试全绿、gold 对拍、确定性、性能不回退、入基准日志、同 PR 删除旧实现、更新状态文档），
> 本文只列新增步骤与各自的门。状态记录写在 `doc_v3/STATUS.md`（格式沿用 doc_v2/STATUS.md）。

---

## 顺序

```
V0 语义修正（C0、C0b）
 └─► V1 串行瓶颈与资源（O1–O5）
      └─► V2 B0：格式无关的消费接口 + 类型化错误（C5）
           └─► V3 B1：BCF 输入
                └─► V4 B3a：BCF → 文本渲染（BCF 输入的 filter 可用）
                     └─► V5 B2：-Ob/-Ou 输出
                          └─► V6 B3b：原生 BCF 手术；O6/O7 扫描与编辑提速
清理项 C1–C4、C7 与 CLI 一致性（04 §2）可穿插在任意两步之间。
```

理由：V0 是正确性问题，先修；V1 改动小、收益确定（4 处串行段，单项最多 10 倍）；BCF 之前先做 B0，否则 BCF 会复制一遍所有消费者。

---

## 各步的门

| 步 | 内容 | 门 |
|---|---|---|
| **V0** | C0 `--mind` 按 plink2 实现；C0b 静态子集与 `--mind` 分开建模、冲突参数报错、`-s ^…` | 与 plink2 的样本集合一致（`--vcf-min-gq` + `--mind` 组合）；gold 用 `-S` 的对拍不变；三策略输出一致 |
| **V1** | O1 报告在 worker 中计算与格式化；O2 报告按请求装配；O3 缓冲池 + mimalloc；O4 并行 `reed index`；O5 stdin 流式 | 输出逐字节不变；03-performance.md §4 的目标逐项入日志；`/tmp` 无残留 |
| **V2** | `SiteFields`/`GenotypeSource`/`RecordRef`，所有消费者迁移；类型化 `EngineError` | 全部现有输出逐字节不变；性能不回退 |
| **V3** | BCF 输入（02-bcf.md §2）：header 与字典、长度前缀链、GT 查表、数值列；`reed stats`、`reed index` 支持 BCF | `clean.bcf` 与 `clean.vcf.gz` 的 stats 输出逐字节相同；差分测试对 noodles；链的切分穷举测试；`stats --pi` 在 BCF 上不慢于 VCF.gz |
| **V4** | BCF → 文本渲染（htslib 格式）；`reed filter` 支持 BCF 输入（VCF 输出） | `clean.bcf` 渲染正文 ≡ `bcftools view -H clean.bcf`；BCF 与 VCF 输入的 filter 输出正文相同 |
| **V5** | `-Ob/-Ou`：显式 IDX 的 header、noodles 记录编码、0 级 BGZF、BCF 的 CSI | `bcftools view` 读回与 VCF 输出一致；字典删除/新增用例；`-Ob` 额外 wall ≤ 30% |
| **V6** | 原生 BCF 手术（B3b）；O6 单遍扫描、O7 按段拷贝；O9 sketch 第二遍（先测量） | 各自的微基准与全量对照；输出不变 |

---

## 与 doc_v2 未完成项的关系

| doc_v2 未完成项 | 安排 |
|---|---|
| R5：`FMT/` 下标 `[0]/[*]`、通用正则 | 保持"明确报 not supported yet"，需求出现再做；不阻塞本路线 |
| R6：BCF 输入与输出 | V2–V6 |
| R7：`--sfs --proj` CLI、多等位逐等位计数 | V1 之后任意时间；`--sfs` 接入时同时完成 C2 中对应累加器的公开 |
| R8：LD、Patterson f 统计、联合 SFS 的 CLI 与对拍 | 同上；每接入一项就把 C2 中对应的 `pub(crate)` 改为公开 |
