# gold/popgen — 统计对拍（软件只作参照物，不进产品代码）

按 `docs/rounds/v2/04-algorithms.md` §5：

| 对拍 | 命令 | 判据 |
|---|---|---|
| π/dxy，pixy 默认（all-sites 输入） | `run_pixy.sh` 的 `pixy_default/` | 逐列一致 |
| π/dxy，`--bypass-invariant-check` | `run_pixy.sh` 的 `pixy_bypass/` | 相对误差 ≤ 1e-9 |
| π，`--variant-only` | `run_vcftools.sh` + `compare.py --vcftools --vcftools-only` | 容差 1e-4：vcftools 逐位点用 n/(n−1) 有限样本校正，pixy/reed 用比值之和，含缺失数据的窗口实测最大相对差 1.05e-5；量级远小于任何真实算法错误 |

```bash
# 1) 参照物（pixy 2.2.4 / vcftools；本机源码装在 ../../.. 的 pixy、scikit-allel）
./run_pixy.sh ../../../clean.vcf.gz ../../../clean.3pop.txt /tmp/popgen_full 100000
./run_vcftools.sh ../../../clean.vcf.gz ../../../clean.3pop.txt /tmp/popgen_full 100000
# 2) reed 输出（两种模式）
reed stats ../../../clean.vcf.gz --pop ../../../clean.3pop.txt --pi --dxy --window 100k \
  --bypass-invariant-check -o /tmp/popgen_full/reed_bypass
reed stats ../../../clean.vcf.gz --pop ../../../clean.3pop.txt --pi --dxy --window 100k \
  --variant-only -o /tmp/popgen_full/reed_vonly
# 3) 比较
python3 compare.py /tmp/popgen_full/reed_bypass /tmp/popgen_full
python3 compare.py /tmp/popgen_full/reed_vonly /tmp/popgen_full --vcftools
```

all-sites 输入（不变位点写作 `ALT="."`）上 pixy 默认模式与 reed 默认模式（不加
`--variant-only`/`--bypass-invariant-check`）逐列一致；reed 在只有变异位点的输入上
默认报错，与 pixy 的 `--bypass_invariant_check` 检查一致。

## fixture 生成器（doc_v4 05 §4）

```bash
python3 make_fixtures.py [--out DIR]     # 默认 gold/popgen/fixtures
```

固定、可手算的小 VCF（无随机数），供各对拍脚本共用：

| fixture | 内容 | 用途 |
|---|---|---|
| `missing_rate.vcf` | 10 样本，位点缺失率 0/0.1/0.2/0.3/0.5/0.9/1.0，另有半调用 | `--geno/--mind` 边界（plink2 `--geno`、vcftools `--max-missing 1−X`） |
| `sfs.vcf` + `sfs.pops` | 4 样本 1 群体，ALT 计数 1,1,2,2,3,3,4,4,5,6,7（共 8 条染色体） | 投影 SFS（easySFS `--proj 6 --unfolded --dtype int -a` → dadi `7 unfolded` / `1 3 3 2 1 1 0`） |
| `multiallelic.vcf` | 双等位对照 + 2/3 个 ALT + 一个缺失基因型；INFO/AC,AN 故意陈旧 | 多等位逐等位重算（`bcftools +fill-tags`，doc_v4 W6） |
| `fstats.vcf` + `fstats.pops` | 3 群体 ×2 样本 | Patterson f2/f3/f4、D（scikit-allel `admixture.py`） |
| `ld.vcf` + `ld.pops` | 6 样本，两个连锁块 | LD r²（plink2 `--r2-unphased`/`--indep-pairwise`） |

## 逐位点/逐样本报告（plink2）

```bash
reed stats in.vcf.gz --freq --hardy --het --missing -o /tmp/rep
plink2 --vcf in.vcf.gz --freq --hardy --het --missing --out /tmp/pl2
python3 compare_reports.py /tmp/rep /tmp/pl2
```

plink2 会把染色体名改写成数字编码（Chr01 → 1），并在 VCF 的 ID 缺失时合成
`chrom_pos`（reed 用同样的规则），因此比较以（染色体编号, 位置）为键。

依赖：pixy 2.2.4、scikit-allel 1.3.13（+ dask）、vcftools、bcftools、plink2。

已知且已记录的差异（04 §4 明确规定的约定，非数值错误）：

1. **末窗坐标**：reed 把每个 contig 的末窗截断到 contig 长度，pixy 打印整窗大小。
   数值不受影响，比较以 `(pop, chrom, window_pos_1)` 为窗口标识。
2. **末窗分母**：`--variant-only` 下 reed 只把 contig 长度内的缺失位置计入分母；
   vcftools/pixy 按整窗（100 kb）计数，把 contig 末端之外的位点也算作单态。
   因此 vcftools 对比显式排除每个 contig 的末窗并打印排除行数（每群体每 contig 一行）。
   （D7 要求 header 必须有 contig 长度，正是为了这条截断规则。）
