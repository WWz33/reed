# 06 路线（W0–W6）

> 通用门沿用 docs/rounds/v2/06 §1 与 doc_v3：测试全绿、gold 对拍、三策略与跨格式 md5、性能不回退并入日志（含负载与中位数，04 §4）、同 PR 删除旧实现、更新 `docs/rounds/v4/STATUS.md`。

```
W0 语义与精确判定（02 首版）──► W0b CLI 拼写回退（02 回退版）──► W1 CLI/文档一致性 ──► W2 多付的工作（O-1、O-2、O-5）
                                                                    └─► W3 reader 吞吐（O-3）──► W4 BCF 原生过滤输入（O-4）
W5 统计接入（SFS → `--fst` 的合并形式 → f 统计 → LD）可与 W0b–W4 并行
W6 多等位逐等位重算
```

| 步 | 内容 | 门 |
|---|---|---|
| **W0** | 首版：阈值类型的精确有理数判定、等位计分子、掩码后求值、头部存活记录作分母、漏斗机制、`to_args()` 往返。**首版的 CLI 改名（`--max-missing*`）与 `--geno/--mind` 的删除已由 W0b 取代，其余全部保留** | 与实现层有关的条目；数字回归（263,420 / 381,672）在 W0b 复核 |
| **W0b** | CLI 拼写回退（02 回退版 §5/§9）：恢复可见的 `--geno X` / `--mind X`（值 = 缺失比例）；`--max-missing`/`--max-missing-sample`/`--max-missing-count` 改为解析期报错 + 精确十进制换算；删 `--max-missing-count`；spec 字段名（`max_missing`/`mind`）、漏斗 id（`after_mask.missing`/`sample.mind`）、`to_args()`、`##reed_filterCommand` 回退；`run_reed.sh` 删 awk；守卫脚本按 02 §7.5 翻转；`docs/rounds/v4/STATUS.md` 的 W0 CLI 行标注 superseded | 02 §7 全部六条（回退版）：边界 fixture 与 plink2 `--geno X` / vcftools `--max-missing 1−X` 一致；subset `--geno 0.1/0.2` → 263,420 / 381,672；gold 正文 md5 `7aedc59ac798…` 不变且 `check_parity.py` PASS；往返测试；三个撤回参数的报错文本单测 |
| **W1** | `--verbosity` 统一（03 §3.3）；`reed index` 拒绝 stdin 并改帮助；`<INPUT>` 帮助统一；去 `-g` 短名；`FilterSpecV2` → `FilterSpec`；`docs/rounds/v3/STATUS.md` 修正 C3；`docs/README.md`；`results.tsv` 补登 V0–V6 | CLI 快照测试；grep 守卫（03 §3.6） |
| **W2** | O-1 报告一次读取；O-2 `FilteredJob` 传 codes；O-5 `--timing` 扩展 | 输出 md5 不变；reports 两族 ≤ 3.2 s；边过滤 `--pi` ≤ 3.0 s |
| **W3** | O-3 reader 吞吐三步（fadvise/大缓冲 → 整 chunk 读 → mmap），每步独立 A/B | `--timing` reader 累计下降；BCF stats ≤ VCF stats |
| **W4** | O-4 `RecordPipeline` 的 BCF 原生输入 | BCF 输入 filter 正文 md5 不变；worker_busy −25%；全量 ≤ 9.5 s |
| **W5** | 统计接入（05 §2 顺序）+ `gold/popgen` 对应脚本 + fixture 生成器 | 每项逐位/逐行对拍通过后才去掉 `#[doc(hidden)]` |
| **W6** | 多等位逐等位重算（05 §3） | 与 `bcftools +fill-tags` 在多等位 fixture 上逐字节一致；默认门行为不变 |

W0 与 W0b 合成一个 PR，先于一切：它们决定用户可见的参数名，越早定稿越少迁移成本。
