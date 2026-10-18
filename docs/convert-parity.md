# 参照工具差异与对拍（convert 及其他）

> README 只列参数与输出文件；参照工具的语义细节、有意分歧与对拍方式集中在本页。
> 规范条文仍是 `docs/spec/cli.md`（§8b convert、§11 差异登记）与
> `docs/rounds/v13/02-architecture.md`；本页是面向用户的汇总。

## convert 与参照工具

对齐目标：vcf2phylip 2.9（Ortiz 2019）、vcftools `--012`、TASSEL hapmap、
vcf2msa。登记在 `docs/spec/cli.md` §8b 的有意分歧：

- 畸形记录跳过并计数（不崩溃、不刷屏）；
- `--resolve-iupac` 确定性（固定 seed，不取系统随机数）；
- `--012` 的半调用与 haploid 写 `-1`（vcftools 写 1）；
- hapmap 的 `alleles` 列写声明集；
- MSA 的 indel 用确定性填充（不调外部对齐器）；
- hapmap→VCF 的头不写 `##contig`（TASSEL 也不写），`##reed_*` 溯源行可关。

布局与字母表：PHYLIP 用 relaxed 布局与 IUPAC 简并码；`--phylip` 是 opt-in
（vcf2phylip 默认输出 phylip）；binary NEXUS 只用双等位二倍体位点，其余写
`?`；hapmap→VCF 时 IUPAC 格展开为 `GT`（`R` → `0/1`）。

往返：在无损域（双等位、二倍体、无缺失的 SNP）上 VCF → hapmap → VCF 逐记录
还原，TASSEL 5 能读回 reed 写的 VCF。

对拍：`gold/popgen/run_convert_parity.sh` 跑六种组合并逐字节比较；8 GB VCF
（348 样本、331 万记录）上五个系统发育输出与 vcf2phylip 2.9 逐字节相同，
耗时 10 秒 vs 其 27 分钟（28 核）——逐样本临时文件转置换成了内存列块。

## 估计量的参照差异

- **f 统计量的标准误**：reed 按碱基对做区块刀切，scikit-allel 按变异数；
  估计值一致，标准误不可比。
- **Tajima's D 的分母**：按窗口平均等位数处的方差（pixy 约定），与 vcftools
  的固定 2N 不同。
- **LD 剪枝**：需要第二遍扫描；`--indep-pairwise` 只存在于 `reed filter`，
  `reed stats` 与 `reed popg` 拒绝它并指向 `reed filter`。

## 对拍脚本

| 脚本 | 对拍对象 |
|---|---|
| `gold/popgen/run_split_parity.sh` | 0.3.x 冻结输出（拆分等价，fixture + subset） |
| `gold/popgen/run_convert_parity.sh` | vcf2phylip / vcftools `--012` / TASSEL / vcf2msa |
| `gold/popgen/run_easysfs.sh` | easySFS（SFS，含二维联合谱） |
| `gold/popgen/run_pixy.sh`、`run_pixy_bed.sh` | pixy（π/dxy，含 BED 窗口） |
| `gold/popgen/run_vcftools.sh`、`run_vcftools_fst.sh` | vcftools（窗口统计、Fst） |
| `gold/popgen/run_plink2_ld.sh` | plink2（LD） |
| `gold/popgen/run_bcftools_af.sh` | bcftools（`-q/-Q` 等位频率） |
| `gold/popgen/run_allel_fstats.sh` | scikit-allel（f 统计量） |
| `gold/popgen/run_hist.sh`、`run_singletons.sh`、`run_popqc.sh` | bcftools stats / vcftools / plink2（v14 新面） |

跑法与判据写在各脚本头部注释；脚本与 fixture 清单见 `gold/popgen/README.md`。