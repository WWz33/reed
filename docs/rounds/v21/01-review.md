# 01 复审：真实数据上核实评审（M1–M11）

> 对象：仓库外报告 `reed-rs-review-2026-10-27.md`（HEAD `fa8131f`）。
> 方法：本机构建 `reed 0.4.0`（`rust-toolchain.toml` 钉的 `1.99` 工具链目录
> 缺 `rustc`，用 `stable` 1.99.1 等价构建），在**真实**文件上重跑每条：
> `subvcf.vcf.gz`（348 样本 / 3000 记录 / 20 contig 头、仅 `Chr05` 有记录）、
> `clean.vcf.gz`（7.0 GB，628,192 BGZF 块）、`clean.bcf`、`clean.3pop.txt`。
> 每条给出复现命令与实测输出【据】。

---

## 1. 评审属实的 must-fix（10/10）

| 报告# | 真实数据复现 | 判定 |
|---|---|---|
| 1 pca/dapc/distance 不能读 BCF | `reed pca sub.bcf` → `malformed record at <二进制>:0: record line has fewer than 9 columns`；同源 VCF 正常（348 samples / 3000 sites）。`distance`、`pca --dapc` 同 | 真 |
| 2 `.tbi` 区域用头顺序取 tid | `.tbi`（tabix，1 ref）→ `invalid reference sequence ID: 4` | 真（"CSI unaffected" 有误，见 §3） |
| 3 `-T/-R` 传 BED 静默空跑 | `-T real.bed` → 0 条、exit 0（bcftools 124 条）；`-R` → `region contig "Chr05\t16000000\t16100000" is not in the header` | 真 |
| 4 `-T CHROM:POS` 被拒 | `-t Chr05:16004641` → `region is missing '-'`（bcftools 1 条）；两列 `-T` 同样静默 0 条 | 真 |
| 5 空 CHROM 首行触发 790 GB 分配 | 首行 CHROM 置空 + `-W` → `memory allocation of 790273982464 bytes failed`，rc 134 | 真 |
| 6 零样本产出非法 VCF | `--remove` 全部 348 样本 → 头 8 列、记录 10 列、exit 0；bcftools 警告并写 8 列 | 真 |
| 7 `rename --clear-ids` 覆盖 `--set-missing-var-ids` | 3000 个 ID 全变 `.` | 真 |
| 8 UPGMA tie 与 hclust 相反 | 报告的 4×4 矩阵：reed 第二次合并 `(s1,s3)`、根高 2.5；R `hclust` 为 `(s0,s2)→+s3→+s1`、根高 2.6667。整数矩阵穷举：`older` 规则 25/25 命中 R，`newer` 19/25 | 真 |
| 9 popg PBS 错标/panic | 真实 pop 文件把 `pop2` 的 119 样本改名 → `--force-samples --pbs pop1,pop2,pop3` **panic** `index out of bounds: the len is 2 but the index is 2`（报告原话）；五组场景下数值错位、标签相同、exit 0 | 真 |
| 10 `pca --solver dense` 不计 n² | `--memory` 以 **MiB** 计；2000×50、`--solver dense --missing meanimpute` 在 16 MiB 预算下通过旧门（2nm=1.6MB）并运行，实际 dense 峰值 ~400 MiB / 2000 样本时栈溢出中止 | 真 |

## 2. 真实数据翻出的新 bug：**M3**

报告第 299 行把 reviewer 报的 `-r` → `input ends inside a record` 判为误报。

- `subvcf.vcf.gz` 上按真实记录区间取 40 个：**36 个失败（90%）**；同批在
  `clean.vcf.gz` 上 0 个。上面 §1 表里 #2 的 `-r Chr05` 也落在其中。
- 机理：`RegionChunkSource::next_chunk` 读到 `end_coffset` 所在块为止（整块），
  该块膨胀文本常常止于半条记录；`pipeline` 收尾 `stitcher.finish()` 见到残留
  尾巴即报 `TRUNCATED_INPUT`。
- 决定性实验（排除索引因素）：把 `subvcf` 的**同一份文本**打成两种 BGZF，
  用 **reed 自己的** `reed index` 建索引：
  | 打块方式 | bcftools `-r` | reed `-r` |
  |---|---|---|
  | 固定 60000 字节切（跨行） | 正常 | `input ends inside a record (truncated file)` |
  | 行边界切 | 正常 | 正常 |
  | 原始 subvcf（2/20 块行对齐） | 正常 | `input ends inside a record` |
- 为什么报告漏判：唯一测过的文件是 `clean.vcf.gz`，它的块**恰好**块块行对齐
  （`bcftools annotate` 重写所致）；`bcftools view -r` 抽子集、`bgzip` 直压
  都会产出跨行块，正是 `subvcf` 的形态。

## 3. 对评审结论的两处更正

1. **#2 的"CSI is unaffected"不成立**：`subvcf` 自带的 `.csi`（20 ref，头顺序）
   与头一致，tid 落到正确 ref，于是失败变成 §M3 的截断；而 `bcftools index -c`
   重建的 `.csi` 只含出现过的 contig（1 ref），与头错位 → `invalid reference
   sequence ID: 4`。即 `.csi` 同 `.tbi` 一样会因 tid 错位失败，取决于索引里
   是否有完整的名字表。
2. **§2 的"Not reproducible"**：见 M3。

## 4. 其余 should-fix 抽查

第 11–22、29、30、32、33、34 条在真实数据上逐条复现（`--hwe 0.05 exact` 静默、
`--hwe midp` 用 "M7" 报错、`--maf 0.4 --max-maf 0.2` 静默空、`-Ob99` 双前缀双
`io error`、`-W` 静默、`-o ''` 写隐藏文件、`theta` 空窗打 0 而 `pi` 打 NA、
D-stat 零分母丢行、`--snps-only=just-acgt` 空操作、`--max-mac 0` 语义、`-S`+
`--keep` 忽略后者、空值 no-op、`pca/distance` 无 `--verbosity`）。修法见
[STATUS.md](STATUS.md) 的 S 系列表。
