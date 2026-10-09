# gold — 硬过滤黄金标准 (filter gold standard)

**目的**：用三个权威工具（bcftools / vcftools / plink2）**并行**跑同一套硬过滤规范，
产出一个**黄金标准 VCF**，供 reed 未来自带的过滤器层做**逐位点对拍**（differential
testing）——reed 过滤输出必须与金标准位点集合一致、同名样本基因型零不一致，才算通过。

这是 reed 的**验证分区**（与 `bench/` 同级）：夹具、脚本与配置**随仓库跟踪**，
生成物一律进 `target/gold/`（`$REED_SCRATCH`），不入库（`docs/rounds/v5/06-layout.md`）。

---

## 1. 目录结构

```
gold/
├── README.md
├── config/filters.conf        # ★ 过滤规范唯一真源 (shell 可 source, python 可解析)
├── fixtures/                  # 小体积边界夹具（随仓库跟踪）
├── scripts/
│   ├── lib.sh                 # 公共：路径/配置/工具检查/样本名单/统计
│   ├── move_scratch.sh        # 一次性：把历史 work/ out/ 搬进 $REED_SCRATCH
│   ├── make_subset.sh         # 切染色体子集（先小集定型参数）
│   ├── profile_input.py       # 画像：密度 + 各字段分布 + 过滤漏斗
│   ├── run_bcftools.sh        # bcftools 分支 (common + INFO → gold)
│   ├── run_vcftools.sh        # vcftools 分支 (common)
│   ├── run_plink2.sh          # plink2 分支 (common)
│   ├── run_reed.sh            # reed 分支（对拍对象）
│   ├── run_all.sh             # 并行编排 + 汇总 + 三方对拍
│   ├── compare_tools.py       # 三方一致性报告
│   └── check_parity.py        # ★ 对拍门：任意过滤 VCF vs gold (原名 parity.py, 2026-10-05 更名, 见 docs/spec/naming.md §6.2)
├── popgen/                    # 群体统计对拍（pixy/vcftools/plink2/…）
└── (生成物不在仓库内)
```

生成物位置：`$REED_SCRATCH/out`（产物）与 `$REED_SCRATCH/work`（中间件），
`REED_SCRATCH` 默认 `target/gold`（仓库根下）。首次迁移跑
`bash gold/scripts/move_scratch.sh`。

---

## 2. 过滤规范

规范分两层，全部参数集中在 `config/filters.conf`：

| 层 | 过滤 | bcftools | vcftools | plink2 |
|---|---|---|---|---|
| common | 双等位 SNP (ACGT) | `-m2 -M2 -v snps` | `--remove-indels --min/max-alleles 2` | `--snps-only just-acgt --min/max-alleles 2` |
| common | FILTER==PASS | `-f PASS` | `--remove-filtered-all` | `--vcf-filter` |
| common | QUAL ≥ q | `-i QUAL>=q` | `--minQ` | `--var-min-qual` |
| common | 基因型 DP∈[lo,hi] → 缺失 | `filter -S .` | `--minDP/--maxDP` | `--vcf-min-dp/--vcf-max-dp` |
| common | 基因型 GQ ≥ g → 缺失 | `filter -S .` | `--minGQ` | `--vcf-min-gq` |
| common | 位点缺失 ≤ m | `F_MISSING<=m` | `--max-missing 1-m` | `--geno m` |
| common | MAF ≥ f | `MAF>=f` (fill-tags) | `--maf` | `--maf` |
| common | HWE p ≥ h | `HWE>=h` (fill-tags) | `--hwe` | `--hwe` |
| common | 样本缺失 ≤ s | `-S keep.txt` | `--keep keep.txt` | `--keep keep.txt` |
| **INFO** | QD/FS/MQ/SOR | `-e 'QD<.. \| FS>.. \| ..'` | ✗ | ✗ |

- **gold** = bcftools 的 common+INFO 全层输出（`out/gold.vcf.gz`）。
- 三方对拍在 **common 层**（`out/bcftools.common.vcf.gz` vs vcftools vs plink2）。
- **样本保留名单**由 `run_all.sh` 在**全量输入**上统一算一次，三工具共享同一文件，
  消除 plink2 `--mind` 在小区间上误删样本的问题。

---

## 3. 本数据集的关键事实与踩坑（务必读）

输入：`../clean.vcf.gz`（348 样本、20 条 Chr、**3,315,604 位点、全为 SNP、0 indel、0 多等位、全 PASS**；
GATK HaplotypeCaller/GenotypeGVCFs 产出）。**已按 MAF 0.05 / max-missing 0.8 过滤过。**

1. **INFO/AC,AN,AF,MLEAC,MLEAF 是陈旧的**：348 样本是大队列的子集（样本名不连续），
   INFO 未重算。例：POS Chr01:358 实际 AN=672，INFO 却写 AN=466。
   → 依赖 INFO/AC 的 MAF/AF 过滤会**用错值**。bcftools 分支必须在 `+fill-tags` 前
   `annotate -x` 剔除这些陈旧标签，让 MAF/HWE 由基因型重算（vcftools/plink2 本就重算）。

2. **bcftools 表达式在 FMT 字段上的 `||` 有 bug**：`-e 'FMT/DP<5 || FMT/GQ<20' -S .`
   会把 99.9% 基因型置缺失；单条件正常。必须用**按位 `|`**：
   `-e 'FMT/DP<5 | FMT/DP>100 | FMT/GQ<20'`。位点级 `&&` 正常。

3. **`+fill-tags` 的 MAF/HWE 是正确的**（与 plink2/vcftools 逐位点一致到 1e-6），
   前提是上一条的陈旧 INFO 已剔除。`hardy` 插件未编译，但 fill-tags 已提供 HWE。

4. **plink2 `--export vcf` 把 `Chr01` 写成数字码 `1`**：分支里用输入 header 的 contig 名
   `bcftools annotate --rename-chrs` 还原，三工具才能用同一 `(CHROM,POS,REF,ALT)` 键对拍。

5. **`bcftools gtcheck` 默认做全对全样本比较**（且默认 `-u PL,GT` 的分数非零）：
   同名样本对拍必须 `-u GT,GT -P <pairs>`（每行 `sample<TAB>sample`）。

6. **本数据是自交物种（大豆）**：全基因组纯合子过量，HWE p<1e-6 的位点占 ~91.5%，
   掩码前后都如此（非掩码假象）。**HWE 过滤对自交物种不适用**，默认关闭。

---

## 4. 输入画像（Chr01-03 子集，480,952 位点）

密度：Chr01 2.56 / Chr02 2.83 / Chr03 4.33 sites/kb（平均间距 231–391 bp）。

| 过滤层 | 保留 | 说明 |
|---|---|---|
| type+PASS+QUAL≥30 | 100.0% | 空操作（QUAL 最小 1201，SNP_FILTER 已在 GATK 侧） |
| 位点缺失≤0.10（掩码后） | 54.2% | 有效 |
| MAF≥0.05 | 97.1% | 近似空操作（输入已 MAF 0.05） |
| HWE≥1e-6 | 8.5% | 自交物种，不适用 |
| **累积** | **3.47%** | HWE 占主导 |

字段分位数：QUAL p50=5.1e4；位点缺失(掩码后) p50=0.089/p75=0.175；
MAF p50=0.24；逐基因型 DP p50=18、GQ<20 占 11.8%。

**据此的默认参数**（`config/filters.conf`）：真正起作用的是 GQ/DP 掩码 + 位点缺失；
QUAL/PASS/双等位/样本缺失 为空操作守卫；HWE 默认 `0`（关闭），若需全参数含 HWE 设 `HWE_MIN=1e-6`。

---

## 5. 用法

```bash
export PATH=$HOME/biosoft:$HOME/biosoft/htslib-1.23:$PATH
export BCFTOOLS_PLUGINS=$HOME/biosoft/bcftools_dir/plugins
export REED_SCRATCH=target/gold        # 生成物去向（默认即此，可覆盖）

cd ~/reed/reed-rs

# (1) 先切子集 + 画像，据分布定参数
bash gold/scripts/make_subset.sh Chr01,Chr02,Chr03 "$REED_SCRATCH/work/subset.Chr01-03.vcf.gz"
python3 gold/scripts/profile_input.py --input "$REED_SCRATCH/work/subset.Chr01-03.vcf.gz" --outdir "$REED_SCRATCH/out"

# (2) 子集上跑金标准（三工具并行 + 三方对拍）
INPUT=$PWD/$REED_SCRATCH/work/subset.Chr01-03.vcf.gz OUTDIR=$PWD/$REED_SCRATCH/out/subset.Chr01-03 bash gold/scripts/run_all.sh

# (3) 全量跑（默认输入 ../clean.vcf.gz）
bash gold/scripts/run_all.sh

# (4) 对拍：把任意「过滤后 VCF」与金标准比
python3 gold/scripts/check_parity.py "$REED_SCRATCH/out/gold.vcf.gz" --gold "$REED_SCRATCH/out/gold.vcf.gz"
python3 gold/scripts/check_parity.py <reed输出的过滤VCF> --json
```

参数覆盖（环境变量优先）：

```bash
MAF_MIN=0.10 HWE_MIN=1e-6 SITE_MISS_MAX=0.2 bash gold/scripts/run_all.sh
```

### 产物（均在 `$REED_SCRATCH/out`）

| 文件 | 含义 |
|---|---|
| `out/gold.vcf.gz` | 黄金标准 (= bcftools full 层) |
| `out/bcftools.common.vcf.gz` | bcftools common 层（用于三方对拍） |
| `out/vcftools.vcf.gz` / `out/plink2.vcf.gz` | 另两工具 common 层 |
| `out/concordance.{txt,tsv,md}` | 三方位点集合/基因型一致性报告 |
| `out/samples.keep.txt` | 共享样本保留名单 |
| `out/*.stats.json` | 各工具位点数/样本数/生效规范 |
| `out/profile.{md,json}` | 输入画像 |

`check_parity.py` 判定：位点集合完全一致 **且** 同名样本基因型零不一致 → `PASS`(退出码 0)，否则 `FAIL`(1)。

---

## 6. 与 reed 的衔接

DESIGN §6「验证策略」要求金标准对拍、计数逐位一致。本目录提供过滤层的对应金标准：
reed 的 L1.5 过滤器层（DESIGN §11/§16 立项）实现后，把其输出交给 `check_parity.py`，
即可验证「过滤即改变 ascertainment」的清单降级规则是否与权威工具一致。

参考实现（只读，不翻译）：vcftools(LGPL)、plink2(GPLv3) 仅行为对照；bcftools(MIT/Expat)
的过滤表达式内核可借鉴（DESIGN §8 许可证政策）。
