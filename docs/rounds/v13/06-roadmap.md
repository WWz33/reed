# 06 实施路线（C0–C4）

> 原则：**先钉真值（fixture + gold 脚本），再垂直切第一格式（phylip 最简，
> 对拍全等），然后一次请求多目标（fasta/nexus 渲染层 + bin 独立模型），
> 最后 SFS 迁移与并行收尾验证**。每步满足 §1 通用门。本轮是新功能轮：
> 没有红测试先行（无缺陷可红），替代物是**对拍真值门**——每步以
> vcf2phylip 2.9 / easySFS 的实测输出为绿。
>
> 不重开 v6–v12 的任何结构项；不触碰 filter/stats 既有语义（stats --sfs
> 只加共享翻译，不改行为）。

---

## 1. 通用验收门（每步）【规】

1. **双 profile 全绿**：`cargo test --workspace`（debug）与
   `cargo test --release --workspace`；
   `cargo clippy --release --workspace --all-targets -- -D warnings`
   不新增警告。
2. **确定性**：本步新增的每个输出，threads ∈ {1, 4} 下 sha256 相同
   （含 `--resolve-iupac`、多目标混合跑法）。
3. **既有面不回归**：`gold/` 既有对拍（filter/stats/SFS）全绿——
   stats 侧只允许"翻译函数提取"造成的纯移动 diff，行为零变化
   （`run_easysfs.sh` 重跑确认）。
4. **退出码与消息**：usage 错 2 / 运行期错 1，错误文本按 02 §5.1。
5. **性能**：每步跑 `bench/suites/e2e/results.tsv` 追加一行；convert
   全量跑法对照 C0-2 基线，最终门见 02 §9.7。
6. **文档同步**：落地一步即更新本目录 STATUS 与 `ci/check-doc-drift.sh`
   对应守卫（01 §6）；`docs/spec/cli.md` 的 convert 节随 C1 落地（02 §10）。

## 2. 步骤

### C0 真值与基线（不写产品代码）

| 项 | 内容 | 门 |
|---|---|---|
| C0-1 | `gold/fixtures/phylo.vcf`（01 §7 的覆盖清单）+ `make_fixtures.py` 生成逻辑 + vcf2phylip 2.9 的期望输出入库（六种跑法：默认集/phylip-only/fasta+nexus/bin/resolve/used_sites，人工核对） | fixture 文件 + 期望文件齐 |
| C0-1b | 新格式 fixtures：`convert012.vcf`（半调用/haploid/杂合分行，01 §7）+ `convert_hmp.vcf`（含未观察 ALT）+ `convert_msa.fa` + 小区域 VCF；期望输出：vcftools `--012` 现跑、TASSEL 读回断言脚本、vcf2msa 无 indel 区域现跑 | fixture + 期望/脚本齐 |
| C0-2 | 基线实测：vcf2phylip 在 clean.vcf.gz 的墙钟（单线程）与 reed `stats --sfs` 的墙钟，记入 STATUS §基线 | STATUS 表两行 |
| C0-2b | 新格式基线：vcftools `--012` 与 vcf2msa（小区域）在 clean.vcf.gz 的墙钟 | STATUS 表两行 |
| C0-3 | `gold/popgen/run_vcf2phylip.sh` + `compare` 辅助（diff -r 语义：目录级全等，B13 分歧用等价集参数） | 脚本可跑（此时 reed 侧失败即红，作未来门） |
| C0-3b | `gold/popgen/run_convert_parity.sh`（012/hapmap/msa 三门，01 §7 第二段） | 脚本可跑（同上） |

### C1 核心垂直切片：`--phylip` 单格式（最简路径）

| 项 | 内容 | 规范 |
|---|---|---|
| C1-a | `reed-core/src/iupac.rs`：AMBIG（62 项）+ GEN_BIN（10 项）常量 + 纯函数 + 单测（含"组合不在表 → malformed 跳过"，B13-9） | 02 §3.2/§3.3/§3.4 |
| C1-b | `ConvertCli`（门规则、`--outgroup/--min-samples-locus` 校验）+ `translate/convert.rs` 骨架 + `Command::Convert` 分发；`--phylip` 单目标 | 02 §1.1 |
| C1-c | `crates/reed/src/convert/`：`ConvertJob`（chunk 内行解析 → 位点门 → 列折叠 `PhyloPartial`）+ reducer 合并 + `render_phylip` + 写出；`run_convert`/`run_convert_opened` 门面 | 02 §2/§3/§7/§8 |
| C1-d | e2e：`phylo.vcf` → `.phy` 与 C0-1 期望 diff 全等；`.log` 五计数正确 | 02 §9.1 |

### C2 多格式渲染 + 调节参数（同一模型）

| 项 | 内容 | 规范 |
|---|---|---|
| C2-a | `PhyloFormats` 三开关 + `render_fasta`/`render_nexus`（纯函数）；三格式同模型不变式测试（02 §9.2） | 02 §2.3/§3.6 |
| C2-b | `--nexus-binary`：独立 bin 列折叠（`GEN_BIN` + 双等位门）+ `render_nexus_bin`；outgroup 对 bin 生效 | 02 §2.4/§3.4 |
| C2-c | `--resolve-iupac`（确定性 seed=0）+ `--write-used-sites` + outgroup 重排 | 02 §3.3/§3.5/§3.6 |
| C2-d | 对拍门全开：六跑法全等（C0-3 脚本转绿）；threads {1,4} sha256 | 02 §9.1/§9.5 |

### C2b 新格式三连：`--012` / `--hapmap` / `--msa`（第二轮增补）

> 三者共享 C1-c 定型的行视图（D13-7：多等位原始索引单一实现，四目标
> 共用）；每项独立可切，顺序 012 → hapmap → msa（复杂度递增：行模型 →
> 行模型 → 参考驱动）。012/hapmap 复用 `--min-samples-locus` 门与既有
> 校验面；msa 引入 `--reference` 装载（FASTA 读入 + 区域解析）。

| 项 | 内容 | 规范 |
|---|---|---|
| C2b-a | `--012`：per-sample `Vec<i8>` 列块折叠（B13-8 半调用/haploid → `-1`）+ `render_012`（转置）+ `.pos`/`.indv` 写出；fixture `convert012.vcf` 对拍 vcftools（半调用样本分行断言） | 02 §3.8、01 §1b.1 |
| C2b-b | `--hapmap`：行缓存模型 + `render_hapmap`（11 固定列 + 声明集 alleles，B13-5）+ TASSEL `BuilderFromHapMap` 读回验证门 | 02 §3.9、01 §1b.2 |
| C2b-c | `--msa`：`--reference` 装载（`Arc` 共享）+ 区域/侧翼解析 + 逐参考碱基折叠（无记录=REF）+ indel 确定性填充（B13-6）+ 每 contig/区域文件写出；对拍 vcf2msa（无 indel diff 全等） | 02 §3.10、01 §1b.3 |
| C2b-d | 三格式 threads {1,4} sha256 + 混合跑法（`--012 --hapmap --msa` 同跑）确定性；drift 守卫 11–13 登记 | 02 §9.5/§9.8 |

### C3 SFS 迁移 + 并行验证 + 文档

| 项 | 内容 | 规范 |
|---|---|---|
| C3-a | `build_sfs_request` 从 `translate/stats.rs` 提取共享；`ConvertTarget::Sfs` 装载；`ConvertJob` 组合 SFS 折叠（方案 A 组合 `StatsJob`，受阻则方案 B 提 `SfsFold`，02 §7.5） | 02 §4/§7.5 |
| C3-b | 等价验收：convert ≡ stats 逐字节（diff -r）+ 往返测试 + `run_easysfs.sh` 双入口 | 02 §4.1/§4.3、§9.4 |
| C3-c | 混合目标跑法（`--sfs --phylip` 同跑）确定性与计数器独立性；性能门：≥ 8× vcf2phylip（clean.vcf.gz、28 核）、SFS-only 退化 ≤ 10% | 02 §9.5/§9.7 |
| C3-d | 文档收口：cli.md §11 convert 节 + B13 差异表（B13-1…8）；README（双语）+ Credits（Ortiz 2019 / vcftools / TASSEL / Chafin）+ Library use 示例进编译门；CHANGELOG；architecture.md 分层图；drift 守卫（01 §6 清单 1–13） | 02 §10 |

### C4 hapmap → VCF 反向转换（第三轮增补，用户点名"hapmap 应该可以相互转换"）

> 真值已在 01 §1b.6 钉死（TASSEL ExportUtils.writeToVCF + BuilderFromHapMap
> 读侧，C0 实测含 hmp→VCF 实跑输出）。形态与前四目标根本不同：
> record→record 流式翻译（D13-8），复用 `-O` 写出机制（§7.7）。

| 项 | 内容 | 规范 | 状态 |
|---|---|---|---|
| C4-a | hapmap 输入探测与装载：`convert::hapmap` 模块——头两行扫描（`##SAMPLE` 注释 + 11 固定列 + 样本名；one/two-letter 判定按首行样本区字符数 ÷ taxa，TASSEL :112-124 同款）、1 MiB 行对齐 chunk 切分 | 02 §3.12/§7.7 | **完成**（`read_header` + 嗅探行回放 + `HapmapSource`） |
| C4-b | 逐行翻译：`chrom`/`pos`/`rs#`/`alleles` 直抄（B13-11）；IUPAC 格 → GT（`reed-core::iupac` 加 `ambig_expand` 逆表）；缺失码（仅 N/?）→ `./.`；`-`/`+`/`0` 是等位字母按 `alleles` 列查下标；`alleles=NA` → 全 `./.` 行；malformed 跳过+计数 | 02 §3.12、01 §1b.6 | **完成**（`translate`/`cell_gt`；one/two-letter 两态） |
| C4-c | 输出：`--vcf` 目标 + `-O v\|z[level]`（`b`/`u` 报 usage）；头在 open 时完整（无 `##contig`，TASSEL 同款）；`OrderedOutput` 流式写出；CLI 门规则 6（与 VCF 目标互斥） | 02 §1.1/§7.7 | **完成**（`run_hapmap_to_vcf` + `run_vcf_cli`；`-o` 扩展名推断 z6；stdout 流式） |
| C4-d | 验收：双等位无缺失位点的 `--hapmap` → `--vcf` **往返保真** e2e（REF/ALT/GT 逐记录还原）；TASSEL 可读性（reed 产的 VCF 被 `-importGuess` 接受）；threads {1,4} sha256；drift 守卫 16–18 | 02 §3.12 往返契约 | **完成**（e2e 往返 + threads 1==8；脚本门 6 TASSEL 读回 4 站点；守卫 16–18 全绿） |

## 3. 顺序与依赖

```
C0（真值/基线） → C1（phylip 垂直切片） → C2（多格式+调节）
       → C2b（012/hapmap/msa 三连） → C3（SFS+并行+文档）
       → C4（hapmap→VCF 反向，独立可切）
```

- C1 是唯一"新执行面"步骤：job/partial/render 的骨架在此定型，C2/C2b/C3
  只加渲染与目标，不动流水线——风险集中在一处。
- C2b 依赖 C1-c 的行视图与位点门实现（D13-7 单一实现点在 C1-c 建立，
  C2b 三目标作为消费者接入）；C2b-a 先行（行模型最简），C2b-c 的参考
  装载是独立新面（FASTA 读入器，C1/C2 未触碰）。
- C3-a 依赖 C1-c 的 `ConvertJob` 形状；若选方案 B（提 `SfsFold`）会触碰
  `stats/run.rs`，必须同 PR 重跑 stats 全量对拍（通用门 3）。
- C4 只依赖 C2b-b（`--hapmap` 写侧）做往返门，代码面独立（新输入探测 +
  新写出路径），可与任何步骤并行。
- **排期紧张时**：C1+C2-d 是最小可发布面（Phylo 三格式+bin）；C2b 三项
  相互独立可任选（012 最常用先做）；C3 可独立后补——两入口等价契约
  保证用户工作流无感；C4 完全独立。

## 4. 本轮不做什么

| 项 | 原因 |
|---|---|
| `reed stats --sfs` 的删除/隐藏 | 兼容承诺（cli.md §0.9）；convert 是镜像入口不是搬家 |
| GT 掩码族进 convert | C13-1【决】；语义与列字母表无交集 |
| `--min-samples-locus` 换算成 `--geno` | C13-2【决】；两个参照语义并存，原样呈现 |
| `--seed` 参数 | 确定性纪律（B13-3）；固定 seed=0 |
| 多维（>2D）联合 SFS | 沿 v5 §3.1 既定边界（easySFS 同界） |
| sites-only VCF（8 列） | v12 §4.2 既定范围外 |
| PHYLIP 交错格式（interleaved）、NEXUS traits/的距离块 | vcf2phylip 只写 relaxed 单块；对齐参照 |
| convert 的 JSON `--report` | 02 §5.3；`.log` 是首版审计面 |
| `--plink`/`--plink-tped`（PED/MAP、TPED/TFAM） | C13-5【决】；plink2 直读 VCF，`--012`+`.pos`/`.indv` 覆盖下标需求；语义台账已在 01 §1b.4 钉死，立项直接引用 |
| `--IMPUTE`/`--BEAGLE-GL/PL`/`--ldhat`/`--ldhat-geno`/`--ldhelmet` | C13-6【决】；相位/似然专门格式；IMPUTE/LDhat 的"整位点零容忍丢弃"与 chunk 归约正交需单独设计；台账在 01 §1b.4 |
| `--vcf` 的样本子集（`-s/-S/--keep/--remove`） | hapmap 的样本集就是列头，没有"子集后重写头"的参照行为；给这些参数 → usage 错（门规则 6） |
| `--vcf` 的 BCF 输出（`-Ob/-Ou`）与 `##contig` 补全 | BCF 记录模型对"从 hapmap 合成"无收益；TASSEL 的 hmp→VCF 也不写 `##contig`，下游要长度走 `--contig-lengths`（B13-12） |
| msa 的 mpileup 掩蔽（vcf2msa `-m`）与 VCF 内 DP 掩蔽 | 另一输入面/与 C13-1 同族；后续按需立项 |
| `--region` 作为**输入过滤器**（非 msa 区域定义） | msa 的 `--region` 是输出区域定义（02 §3.10）；输入级区域过滤走 `reed filter` 管线，不重复 |

## 5. 风险与开放问题

| 风险 | 缓解 |
|---|---|
| `GtCode` CollapseAlt 与多等位列视图的张力（D13-3/D13-7） | C1-c 在行文本上自取 GT 字段（列折叠独立于 stats 扫描器），不强行复用 `GtCode`；若性能不达标再评估扩展 scanner 原始视图 |
| `StatsJob` 组合可能带进空窗 grid 的无谓开销 | 02 §7.5 方案 B 兜底（提 `SfsFold`）；C3-a 实测后二选一 |
| 大样本 × 高位点数的收尾渲染内存（700×3M ≈ 2.1 GB 序列） | 02 §7.3 已定：首版即流式渲染（逐行从列块拉取），不整块物化 |
| vcf2phylip 期望输出里包含其 stdout 的随机性（resolve 跑法） | C0-1 的 resolve 期望用固定种子无法从 vcf2phylip 获得——resolve 跑法的对拍降级为**结构断言**（长度、字母表 ⊆ {A,C,G,T,N}、杂合位点全部拆解为纯合字母），不 diff 全文；B13-3 已登记 |
| TASSEL hapmap 无独立写出 CLI 可对拍 | 写侧对拍降级为**读回验证**（TASSEL `BuilderFromHapMap` 读 reed 输出断言 taxa/site/等位）+ `alleles` 列的声明集断言（B13-5 守卫 13） |
| msa 全位点输出的 I/O 体量（区域长 × 样本数） | 流式渲染（02 §7.3 同款）；无记录位点零解析（REF 复制）；`--region` 强制用户显式给区域时全 contig 输出给一次性警告 |
| `--012` 对拍受 B13-8（半调用/haploid）分歧影响 | fixture `convert012.vcf` 把半调用/haploid 样本与正常样本**分行**（不同位点），对拍断言分别写：正常位点 diff vcftools 全等，半调用位点断言 reed 契约（`-1`） |
| hmp→VCF 的 REF 选择在任何实现下都有损（hapmap 不存 REF） | B13-11 定为"照抄声明序"——同输入同输出、可预测；**无损域**（双等位、REF 在 `alleles` 首位、无缺失）由往返 e2e 与脚本门 6 双断言，多等位/缺失/indel 明确排除在保真契约外（02 §3.12） |
| 头嗅探（首行数据行判单/双字母）可能吞掉一条记录 | `read_header` 把嗅探行存进 `HapmapHeader.first_row`，`HapmapSource` 在新 source 首块回放——脚本门 6 的 4 站点全等即此路径的回归断言 |
