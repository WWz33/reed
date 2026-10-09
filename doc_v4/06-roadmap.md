# 06 路线（W0–W6）

> 通用门沿用 doc_v2/06 §1 与 doc_v3：测试全绿、gold 对拍、三策略与跨格式 md5、性能不回退并入日志（含负载与中位数，04 §4）、同 PR 删除旧实现、更新 `doc_v4/STATUS.md`。

```
W0 参数语义统一（02）──► W1 CLI/文档一致性 ──► W2 多付的工作（O-1、O-2、O-5）
                                              └─► W3 reader 吞吐（O-3）──► W4 BCF 原生过滤输入（O-4）
W5 统计接入（SFS → --fst-pops → f 统计 → LD）可与 W2–W4 并行
W6 多等位逐等位重算
```

| 步 | 内容 | 门 |
|---|---|---|
| **W0** | `CallRate` 类型（03 §3.2）；`--max-missing`、`--max-missing-count`、`--max-missing-sample`；删除 `--geno/--mind`（报错 + 换算提示）；spec 字段改名 `min_call_rate`/`max_missing_alleles`/`min_sample_call_rate`；漏斗 id；`to_args()`；`run_reed.sh` 换算；分子改为等位计 | 02 §7 全部六条：边界 fixture 与 vcftools/plink2 一致；subset 263,420 / 381,672；gold 正文 md5 `7aedc59ac798…` 不变且 `check_parity.py` PASS；往返测试；精确判定单测 |
| **W1** | `--verbosity` 统一（03 §3.3）；`reed index` 拒绝 stdin 并改帮助；`<INPUT>` 帮助统一；去 `-g` 短名；`FilterSpecV2` → `FilterSpec`；`doc_v3/STATUS.md` 修正 C3；`docs/INDEX.md`；`results.tsv` 补登 V0–V6 | CLI 快照测试；grep 守卫（03 §3.6） |
| **W2** | O-1 报告一次读取；O-2 `FilteredJob` 传 codes；O-5 `--timing` 扩展 | 输出 md5 不变；reports 两族 ≤ 3.2 s；边过滤 `--pi` ≤ 3.0 s |
| **W3** | O-3 reader 吞吐三步（fadvise/大缓冲 → 整 chunk 读 → mmap），每步独立 A/B | `--timing` reader 累计下降；BCF stats ≤ VCF stats |
| **W4** | O-4 `RecordPipeline` 的 BCF 原生输入 | BCF 输入 filter 正文 md5 不变；worker_busy −25%；全量 ≤ 9.5 s |
| **W5** | 统计接入（05 §2 顺序）+ `gold/popgen` 对应脚本 + fixture 生成器 | 每项逐位/逐行对拍通过后才去掉 `#[doc(hidden)]` |
| **W6** | 多等位逐等位重算（05 §3） | 与 `bcftools +fill-tags` 在多等位 fixture 上逐字节一致；默认门行为不变 |

W0 单独成一个 PR，先于一切：它改变用户可见的参数，越早落地越少迁移成本。
