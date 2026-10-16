# reed doc_v13：`convert` 子命令——Phylo 多格式 + SFS 迁移 + 012/hapmap/msa + hapmap↔VCF 双向（2026-10-16）

> 本目录是 `docs/rounds/v12/` 之后的**第十三轮规范**，性质是**新功能轮**
> （非缺陷修复）：为 reed 新增 `reed convert` 子命令——vcf2phylip 的
> PHYLIP/FASTA/NEXUS/二元 NEXUS 转换（并行化）+ SFS 出口从 `stats --sfs`
> 的等价迁移，统一门面 API（`ConvertRequest`/`run_convert`）。
> 初版 doc 只写规范不改代码（用户指令 2026-10-16：先写 doc）；随后按用户
> 点名分三轮增补 `--012`/`--hapmap`/`--msa` 与 `--vcf`（hapmap 反向），
> 并全部落地（状态见 STATUS.md）。
> 方法：参照工具源码通读（`../../vcf2phylip/vcf2phylip.py` v2.9、
> `../../easySFS/easySFS.py`）+ reed 侧可复用面盘点（§4）+ 既有对拍脚本
> （`gold/popgen/run_easysfs.sh`）确认 SFS 行为基线。

---

## 背景

reed 目前有三个子命令（filter/stats/index）。两个转换类需求缺失：

1. **系统发育输入**：vcf2phylip（Ortiz 2019，DOI 10.5281/zenodo.2540861）
   是 VCF → PHYLIP/FASTA/NEXUS/二元 NEXUS（SNAPP）的事实标准工具，但它是
   单线程 Python（README 自述 20 GB VCF ~27 分钟；其转置阶段对每个样本
   重读一遍临时文件，O(N·L) 重复 IO）。reed 的融合式 chunk 流水线正是
   此负载的形状：列折叠是每行纯字节操作、无跨位点依赖（天然 seam-free）。
2. **SFS 出口位置**：`--sfs --proj` 挂在 `reed stats`（窗口统计命令）下
   语义别扭——投影 SFS 是 genome-wide 的格式转换（VCF → dadi/fsc2 谱文件），
   与"转换"同类。迁移到 `convert` 后 stats 入口保留（兼容承诺）。

两者共享同一执行形状（读一遍 VCF、按位点折叠、收尾写一组格式文件），
正当化一个新命令面，而不是继续膨胀 stats 的参数表。

**第二轮增补（2026-10-16，同日）**：用户点名再加三个格式——`--012`
（vcftools 的 0/1/2 基因型矩阵）、`--hapmap`（TASSEL 的 .hmp.txt）、
`--msa`（vcf2msa.py 的参考序列全位点 MSA），源码都在工作区
（`../../vcftools`、`../../tassel`、`../../vcf2msa.py`）。顺手盘点并钉死
了 vcftools 的其余转换格式（plink/plink-tped/IMPUTE/BEAGLE-GL/LDhat/
LDhelmet，01 §1b.4 台账 + C13-5/6 裁决不进首发轮）。命令面从三组扩为
**五组**（SFS / Phylo / 调节 / 012-hapmap / MSA），执行模型不变（新目标
是 ConvertJob 的平级成员，D13-7 保证行级解析单一实现）。

**第三轮增补（2026-10-16，同日）**：用户点名"hapmap 应该可以相互转换"
——补 `--vcf` **反向目标**（hapmap → VCF，路线 C4）。真值来自 TASSEL
`ExportUtils.writeToVCF` 与 `BuilderFromHapMap` 的读侧（01 §1b.6，含 C0
实测的 hmp→VCF 输出）。形态与前五目标根本不同：**record→record 流式
翻译**（D13-8），不折叠成矩阵，复用 filter 的 `-O` 写出机制（无新压缩
路径）；`--vcf` 与所有 VCF 输入目标互斥（门规则 6）。新增 B13-11/12 两条
有意分歧与守卫 16–18。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 语义对照台账：vcf2phylip 2.9 逐行号语义（§1）、**vcftools --012 / TASSEL hapmap（含 hmp→VCF 反向真值 §1b.6）/ vcf2msa / vcftools 其余格式与排除清单（§1b.1–§1b.6，第二/三轮增补）**、easySFS/stats 现状（§2）、语义歧义点核实（§3）、reed 可复用面（§4）、差异分档 D13-1…8 / B13-1…12 / C13-1…6（§5）、守卫清单 1–18（§6）、对拍真值门设计（§7）、性能**实测基线**（§8：全量 158× 且五文件逐字节） |
| [02-architecture.md](02-architecture.md) | 修订规范 N13-1…N13-11：命令面与门规则（§1）、ConvertRequest/多目标模型（§2）、Phylo 语义（位点门/IUPAC/双等位 bin/输出布局，§3）、**hapmap→VCF 反向目标（§3.12）**、SFS 迁移等价契约（§4）、错误与 .log（§5）、命名与许可（§6）、并行执行模型（chunk 折叠/partial 合并/确定性，**§7.7 反向目标的 record→record 流式**）、门面 API（§8）、验收门（§9）、文档回改（§10） |
| [06-roadmap.md](06-roadmap.md) | 实施路线 C0–C4：真值与基线 → phylip 垂直切片 → 多格式渲染 → **012/hapmap/msa 三连（C2b，第二轮增补）** → SFS 迁移+并行验证+文档 → **hapmap→VCF 反向（C4，第三轮增补）**；通用验收门、"本轮不做"、风险与开放问题 |
| [STATUS.md](STATUS.md) | 实施状态（C0–C4 全部完成）+ 差异登记落地状态 + 基线表（已实测） |

## 结论速览

| 维度 | 结论 |
|---|---|
| 命令面 | `reed convert <INPUT> --sfs / --phylip / --fasta / --nexus / --nexus-binary / --012 / --hapmap / --msa / --vcf [调节参数]`；六组至少给一个目标，否则 usage 错；`--phylip` 是 opt-in（B13-1）；`--msa` 须配 `--reference`；`--vcf` 读 hapmap 输入，与 VCF 输入目标互斥（门规则 6） |
| 语义模型 | Phylo 三格式共享一份列模型（一次转换、三种渲染，D13-1 不变式）；bin.nexus 独立字母表（0/1/2/?）但共享行解析与位点门；多等位 GT 用**原始等位索引**（GtCode 的 CollapseAlt 不可用，D13-3）；`--vcf` 是 record→record 翻译（D13-8），不折叠 |
| 并行 | 完全复用融合式 chunk 流水线：worker 内一次行解析多目标消费，partial = 列块 + SfsPartial + 计数器，reducer 按 chunk 序加法合并 → 天然确定性与线程无关；多目标共享同一 worker 池不乘线程；首错即停走既有 FirstError 通道（v8/v9 纪律）；`--vcf` 用同一 run 驱动 + `OrderedOutput` 流式（§7.7） |
| SFS 迁移 | `convert --sfs` ≡ `stats --sfs` **逐字节**（含 .log/退出码/目录布局）；stats 入口不删不改；翻译层共享单一 `build_sfs_request`（D13-6），折叠单一实现（D13-5），往返测试锁死双入口 |
| 双向 | `--hapmap` 与 `--vcf` 互为反向：无损域（双等位、二倍体、无缺失）上往返逐记录还原；TASSEL 5 两个方向都读得回（脚本门 3/6） |
| API | `ConvertRequest { targets: Vec<ConvertTarget> }` + `run_convert(path, &request, &budget)`；请求按引用、`to_args` 可失败、根导出——v11 门面纪律直接沿用 |
| 有意分歧 | 12 条（B13-1…12：phylip opt-in、畸形行跳过+计数、resolve 固定 seed=0、outgroup 报错、hapmap alleles 声明集、msa indel 确定性填充、msa 样本名不折叠、012 半调用/haploid=-1、AMBIG 未命中跳过（参照崩溃）、hapmap 字母表收敛、hmp→VCF 直抄声明值、hmp→VCF 无 `##contig` 且 `##reed_*` 可关）——全部登记 cli.md §8b，分歧点对拍降级为分 fixture/读回验证 |
| 待拍板 | 6 条（C13-1…6：GT 掩码族不进、min-samples-locus 独立、--pop 可缺省、msa 样本名不折叠、plink 不进首发、相位/似然格式不做）——默认已定，可推翻 |
| 许可 | vcf2phylip 与 vcf2msa 均为 GPLv3：**只实现行为，不复制代码**；AMBIG 表以 IUPAC 标准编码为准；Credits 引 Ortiz 2019（vcf2phylip）、vcftools、TASSEL、Chafin（vcf2msa） |

## 与 v11、v12 的关系

- **v11（门面 API 轮）**：convert 的门面面（§8）直接执行 v11 纪律——
  请求按引用（A11-2 教训）、`to_args` 可失败（A11-3 教训）、报告类型根
  导出（N6 先例）、README 示例进编译门。
- **v12（INFO 兼容性轮）**：正交（convert 不触碰 INFO 编辑/表达式/编码
  路径）；但 §0.5 的两个元模式直接塑造本轮设计——"多份实现各自漂移"
  → D13-1（三格式单一模型）、D13-5（SFS 折叠单一实现）、D13-6（翻译层
  共享）是预防性约束；"元数据存在但消费者没查"不适用（新代码无遗留
  元数据）。
- **v8/v9（线程/错误路径轮）**：§7 并行模型全部在 ThreadBudget 记名、
  无匿名线程、无 FirstError 旁路；列折叠无 seam（不触发 v9 的接缝分析
  分支）。
- **v5/08 §6.2（门面分工）**：`reed-cli` 仍是唯一前端；convert 的全部
  新代码在门面 `reed` 与 `reed-core`（IUPAC 常量属 L3 纯函数），无第
  三前端。

## 决策状态

| 决策 | 状态 |
|---|---|
| C13-1 GT 掩码族（--min-dp/--max-dp/--min-gq）是否进 convert | **待拍板**；默认不进（02 §2.5） |
| C13-2 `--min-samples-locus` 独立参数还是复用 `--geno` | **待拍板**；默认独立（02 §3.7） |
| C13-3 SFS 目标 `--pop` 是否必选 | **待拍板**；默认可缺省（沿 stats，02 §1.2） |
| C13-4 msa 样本名是否折叠 `.` 后缀 | **待拍板**；默认不折叠（02 §3.10） |
| C13-5 `--plink/--plink-tped` 是否进首发轮 | **待拍板**；默认不进（02 §3.11） |
| C13-6 `--IMPUTE/--BEAGLE-*/--ldhat*/--ldhelmet` 是否实现 | **待拍板**；默认本轮不做（06 §4） |
| B13-1…8 八条有意分歧 | 已定（等价集对拍 + cli.md 差异表登记） |
| B13-9…12 四条有意分歧 | 已定并落地（B13-9 跳过路径、B13-10 字母表收敛、B13-11/12 反向目标抄写规则） |
| ConvertJob 组合 StatsJob（方案 A）还是提取 SfsFold（方案 B） | **实现期二选一**（02 §7.5：规范只锁"单一折叠实现"，不锁组合方式） |

## 标记约定

【规】必须；【荐】推荐；【决】需要拍板。
证据通道：**【据】**vcf2phylip.py / easySFS.py / reed 源码行号；
**【实】**gold 对拍脚本实测（对拍门随 C0-3 建立）。
