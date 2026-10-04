# reed 设计决策记录 (DESIGN.md)

> 来源:2026-10-04 grilling 会话 + 14 个参考项目源码调研 + 自有数据基准测试。
> 标记约定:【定】= 用户已明确拍板;【荐】= 建议(可改);【据】= 有当日实测/源码证据。

## 1. 定位

- 纯 Rust **库**(lib 优先,CLI/绑定后置)【定】。
- 目标:群体遗传学统计的**规范实现**——每个 `(统计量, 估计器)` 是完整不可拆的定义:精确公式 + 分母约定 + 缺失数据处理 + 文献引用 + 金标准对比测试。
- 输入格式:**仅 VCF/BCF**(CSI 索引;PGEN 明确砍掉)【定】。gVCF 风格 BCF(含 `<NON_REF>`/单态记录)属于 VCF/BCF 范畴,不砍。
- 倍性:v0.1 二倍体;多等位位点:v0.1 拒绝(显式报错,不许静默丢弃)。
- 正确性优先,速度第二;速度阶段的主路线 = 2-bit SWAR(见 §7)。

## 2. 架构:三层

```
L1 传输   noodles 0.117 单后端(惰性记录流 + header 元数据),不写解析代码【定·实测支持】
L2 物化   窗口缓冲 WindowBuffer:把记录解码进定长窗口,内存 O(窗口)
L3 统计   估计器 = 封死的完整定义;纯函数 over 不可变窗口;计数式结果
```

后端选型证据【据】:自有数据(348 样本 × 3,315,604 位点)

| 负载 | noodles | rust-htslib | 说明 |
|---|---|---|---|
| clean.bcf 8GB(8 个 FORMAT 字段) | **110.4 s** | 129.1 s | noodles 惰性只解 GT 列胜出 |
| clean.vcf.gz 7GB(同字段) | **187.3 s** | 302.9 s | 同上 |
| clean_bg.vcf.gz(GT:PS 两字段) | 90.7 s | **49.2 s** | htslib 块级解包便宜,反超 |

机制:htslib 一碰 FORMAT 就解包整块(`bcf1_t` 注释明说 lazy=按记录块);noodles 按列解。popgen 负载(GT 为主)下 noodles 架构正确且实测不落下风。所有输出计数两库与 bcftools **逐位一致**(金标准校验模板已跑通)。

## 3. 数据模型(窗口缓冲)

- `contig_len` 是**一等公民**,缺失即构造/读取报错(Beagle 输出丢 `##contig`,实测案例)。禁止静默用 pos 差近似。
- 基因型表示(v0.1):`Vec<i8>`,每样本 2 个等位基因,-1 = 缺失;平行 `Bits` 保存缺失与**逐基因型相位位**(实测:GT 与 PGT 物理定相分离,相位必须逐基因型而非数组级布尔)。
- 统计量缺失处理(属估计器定义,不可配置),按族固定:
  - π/dxy/Fst/PBS/f-统计:**逐等位基因删除**(allele-level deletion);
  - SFS/θ_W/Tajima's D:**complete-case**(组内任一等位基因缺失即跳过该位点并计数);
  - HWE/杂合度/LD:**完整基因型个体**(任一等位基因缺失的个体被排除)。
- 单态记录(`is_variant=false`)允许存在于缓冲,分母估计器(逐样本对)依赖它。

## 4. Stats Model(核心设计)

三根正交轴,禁止互相污染:

### 4.1 输入清单 SiteCompleteness(数据事实,一次声明)
```rust
enum Ascertainment { Complete, Filtered, Unknown }
enum SiteCompleteness { WithInvariants, VariantOnly(Ascertainment) }
```
- 决定“文件里有什么”,在窗口构造时生效;每个估计器从清单推导可用性。
- `VariantOnly(Unknown)` 对分母敏感统计量**默认拒绝**(pixy 同款立场);放行必须显式声明。【定】
- 清单与文件事实可以矛盾(用户谎报);库提供 `DetectedFacts` 扫描辅助(has_monomorphic, contigs_with_lengths),不做自动魔法。

### 4.2 估计器:封死的定义(公式+分母+缺失语义+引用不可拆)
- 用户选**估计器**,永远选不到“错的分母配错的公式”。
- 可用性规则(实现于 `availability()`):

| 估计器 | WithInvariants | VariantOnly(Complete) | VariantOnly(Filtered) | VariantOnly(Unknown) |
|---|---|---|---|---|
| π / dxy · PerPairComparable(pixy 语义) | 可用 | 有偏(不变位点缺席,分母只数变异位点) | 拒绝 | 拒绝 |
| π / dxy · AccessibleBp(allel 语义,分母=窗口 bp) | 可用 | 可用 | 有偏(被滤掉的真变异不可辨) | 拒绝 |
| Fst · Hudson / WC84 / Nei86 | 可用 | 有偏提示 | 有偏提示 | 有偏提示 |
| PBS · Yi2010 / PbsN1(Fst 之上) | 可用 | 有偏提示 | 有偏提示 | 有偏提示 |

Fst 特例依据:Hudson 的分子分母在同一位点集合上自洽,缺失不变位点不会结构性破坏它(只引入 ascertainment);π/dxy 的分母是位点普查,缺不变位点即分母错误。故 Fst/PBS 的门槛是提示,π/dxy 的门槛是拒绝。
注:θ_W 的分母同为窗口 bp、分子是分离位点计数,与 AccessibleBp 同理应在 VO(Unknown) 下拒绝——**已知缺口,v0.2 修**(spectrum 族尚未接入 availability 门,详见 §12 审计记录)。SFS/Tajima/HWE/LD/f-统计按 Fst 类处理(有偏提示)或按位点内统计豁免,接入门时逐一判定。

### 4.3 结果:永不裸浮点
所有结果携带精确整数计数(分子/分母/缺失对)+ 清单偏差标记 + 窗口坐标;`value()` 是派生的 f64。
- 作用:数字可审计;空分母自动 `None` 而非垃圾值;v0.2 位并行优化以“计数逐位一致”为回归锚点。

### 4.4 守卫
- 窗口**永不跨 contig**。
- 小样本/小窗口:`comps == 0` → `value() = None`;Tajima's D 的 `min_sites` 守卫(allel 先例)在实现时同样内置。【荐】

### 4.5 命名
对齐 sgkit 词汇(`call_genotype`/`variant_position` 精神)为【荐】,非【定】;Rust 侧类型名以本库自洽为准,文档记录映射。

## 5. API 铁律

- 迭代器/回调式访问;**禁逐样本随机访问**(noodles `series.get(i)` 实测 O(n²),正确姿势 `series.iter(header)` 一趟流过)。
- 窗口网格:size + step(滑动窗口 = step < size),从 contig 首个记录所在窗口到末个记录所在窗口,空窗口照发(网格一致),`value()=None`。
- 读取器:`for_each_window(path, opts, |window| ...)` 回调式流(真正流式);`read_windows()` 便利包装(小数据/测试)。

## 6. 验证策略

- **金标准**:同一 fixture 上,allel / sgkit / pixy 与本市实现四方对比(sgkit 自身即用 `test_divergence__windowed_scikit_allel_comparison` 对拍 allel,先例已验证可行)。
- **计数校验**:整数计数逐位一致是硬门;浮点值容差 1e-12。
- **回归锚点**:v0.1 i8 实现 → v0.2 SWAR,分子/分母/缺失计数不得变化。
- 测试工程: `reed/bench/`(noodles vs rust-htslib vs bcftools 三方计数校验,已跑通)保留作 I/O 层基准与兼容性看门狗。

## 7. 路线图

- **v0.1(骨架 + 统计族全量,2026-10-04)**:WindowBuffer / PopulationView / SiteCompleteness / 读取器(VCF、VCF.gz、BCF);统计族实现状态见 §12。
- **v0.2**:2-bit SWAR 物化(GCTA/plink2 位并行模式,一 u64 = 32 样本);Fst Weir-Cockerham 84 + Nei 86(估计器扩展);Tajima's D(min_sites 守卫,需 Watterson θ);区域并行(contig 级 rayon + CSI);流式 Iterator API。
- **v0.3**:snputils 语义的 f2/f3/f4/D-stat + block jackknife(BSD 白名单,可移植);LD r²;HWE;PCA(flashpca 思想)。
- 数据侧待办:管线重跑 gVCF(`GenotypeGVCFs --include-non-variant-sites`)以获得 `WithInvariants` 输入,解锁逐样本对估计器全量语义。

## 8. 许可证政策(翻译白名单)

- 可翻译/移植(带署名):MIT — noodles、rust-htslib、htslib(本体 MIT/Expat)、scikit-allel、pixy、vcflib、vcfstats;Apache-2.0 — sgkit;BSD-3 — snputils。
- **只读禁区**(只看代码记坑 + 按论文/规范重写 + 行为对比,禁翻译):vcftools(LGPL)、plink-ng(GPLv3)、flashpca/GCTA/SNPGenie(GPL)。PGEN 格式自 spec(`pgen_spec.tex`)clean-room,不受 GPL 约束。

## 9. 参考项目地图(14 家;源码级统计族记录见 [`docs/statistics-reference-map.md`](docs/statistics-reference-map.md))

11 个并行子代理已对 14 家逐一完成**源码级统计族审读**(file:line + 公式摘录 + 分母/缺失语义,6 处代表引用抽查属实),汇总为参考实现地图:逐族跨项目对照、金标准认定(Hudson=allel 正典、WC84=allel+plink2+GCTA 三家一致、f-统计=snputils/allel、HWE mid-p 看 plink2)、与 reed 落库计划(f-统计无偏校正+块刀切、HWE mid-p、LD 缺失语义选项等,见地图 §4)。

模型先例:noodles/rust-htslib/htslib(记录流)、allel/pixy(位点数组)、sgkit(xarray 数据集,无 contig_length 之教训)、plink2/PGEN、GCTA(uint64 块 + 缺失位向量)、flashpca(稠密矩阵) 。

## 10. 工程环境事实(本机)

- rustc:用 rustup stable(1.99);conda 的 1.90 低于 noodles 0.117 MSRV(1.91)。
- crates 源:清华 sparse(aliyun 索引滞后,noodles-vcf 0.94 缺失);noodles 本体用本地 path 依赖。
- 构建产物必须放本地盘(`target-dir` 已在 `.cargo/config.toml` 指到 /tmp)——/sri 并行文件系统不支持构建脚本的并发删除。
- rust-htslib 构建备查:LIBCLANG_PATH + BINDGEN_EXTRA_CLANG_ARGS(见 bench/)。本库 v0.1 不依赖 htslib。
- codegraph 已索引全部 14 个参考仓库。

## 11. 未决/远期

- Beagle 文件缺 `##contig` 的兜底(外部 lengths 文件注入)【荐,未实现】。
- htslib 作为可选第二后端的门接口(Reader trait 抽象)何时引入——先观察 v0.1 使用。
- 多线程基准(htslib set_threads vs noodles multithreaded_reader)与 contig 级 rayon 的实测对比。

## 12. 统计族全表与实现状态(2026-10-04)

金标准工程:`/tmp/allelvenv`(scikit-allel 1.3.13)+ `/tmp/gen_golden.py`;HWE 精确检验另用 Levene-Haldane 分布独立枚举交叉验证。27 个测试全绿。

| 族 | 统计量(估计器) | 状态 | 金标准来源 |
|---|---|---|---|
| 频率 | allele counts(分组聚合);SFS 折叠/展开 | ✅ 实现 | 手算 |
| 多样性 | π(PerPairComparable / AccessibleBp);Watterson θ;Tajima's D(min_sites 守卫) | ✅ 实现 | scikit-allel |
| 分化 | dxy(两估计器);Fst:Hudson / Weir-Cockerham 84 / Nei 86;PBS:Yi2010 / PBSn1 | ✅ 实现 | scikit-allel(WC84 组件 a=0.125,b=0,c=0.75 逐位一致;PBS 两定义与 `allel.pbs` normed 开关精确对拍) |
| HWE-近交 | HWE 精确检验(Wigginton 2005,sgkit 可译实现);Hobs/Hexp/Fis 逐位点 | ✅ 实现 | Levene-Haldane 独立枚举 |
| LD | Rogers–Huff r²(逐对 + 窗口汇总) | ✅ 实现 | scikit-allel |
| f-统计 | f2 / f3 / f4 / D-stat(点估计) | ✅ 实现;block jackknife 未(计划 v0.3,snputils 语义) | 手算 |
| 结构 | PCA / GRM / relatedness / IBD | ⏳ 未实现(需全基因组矩阵积累层) | — |
| 选择-相位 | EHH / iHS / XPEHH / nSL / Garud H | ⏳ 未实现(需相位扫描层;Garud H 需单倍型哈希) | — |
| 密码子 | πN/πS / dN/dS | ⏳ 远期(依赖注释) | — |

封死的定义细节(实现即文档):SFS/θ/Tajima 仅用组内 **complete-case 位点**(缺失位点跳过并计数);π/dxy 的 AccessibleBp 分子为位点均值之和、分母为窗口 bp;Hudson Fst 为**逐位点归一化求和**(allel `hudson_fst` 语义:num_s = mpd_AB − (mpd_A+mpd_B)/2,den_s = mpd_AB,位点上任一项无定义则跳过);WC84 组件按位点、双等位基因求和后取 a/(a+b+c),n = an//2 地板语义、s² 与 n_C 各带 1/(r−1) 因子(allel fst.py 逐行对齐);PBS 的 Fst 为窗口级逐位点归一化 Hudson、clip 到 [0, 0.99999](T 上限 ≈11.5)、T=−ln(1−Fst)(与 scikit-allel 完全一致);f-统计/D 为逐位点频率乘积的均值(ALT 频率,逐等位基因删除)。

### 12.1 三代理审计记录(2026-10-04,数学/边界/契约三视角)

审计发现并经本人独立复核(数值实验对拍 allel)证实的错误,**已全部修复**:

1. **Hudson Fst 用了池化比值而非逐位点归一化**(CRITICAL):等覆盖时两式等价,不等覆盖/缺失数据下偏差可达 35 倍(实测 −0.0307 vs −1.0847)。修复:fst.rs/pbs.rs 改为 allel 逐位点语义。原有金标准 fixture(等分 4 样本)无法照亮此错——**测试盲区教训:金标准必须包含不等覆盖 + 缺失数据的随机 fixture**(已补 `tests/unequal_groups.rs`,3 项对拍 allel 逐位一致)。
2. **WC84 的 s² 与 n_C 各缺 1/(r−1) 因子,且 n 未按 allel 的 `an//2` 地板**(CRITICAL):同数据 a/b/c 三项全偏(b 差 9 倍)。修复后与 allel 逐位一致。
3. **PBS clip 上界用了 1−ε 而非 allel 的 0.99999**(MEDIUM):Fst=1 窗口的 T 上限 36 vs 11.5。已改。
4. **Fst/PBS 的 PairCounts.missing 恒为 0、PBS 丢弃 between 计数**(契约违规):审计线索失真。已补。
5. **pos ≥ contig_len 的记录被静默丢弃**(契约违规):现已显式报错 `PositionBeyondContigLength`。
6. **全缺失窗口的 AccessibleBp 返回 Some(0.0)**:现改为 None。
7. **D-stat 分母守卫 `!= 0` → `> 0`**。已改。

审计指出、**评估后暂不采信/降级**的项:AccessibleBp 分母用窗口跨度而非真实 accessible-bp 计数(设计如此——WithInvariants 语义下位点齐全;位级掩码是 v0.3 的 `is_accessible` 特性);Tajima 方差字段与 d 的重复计算(等价,重构时机 v0.2);结果字段全 pub 可事后篡改(封死靠约定与文档,v0.2 考虑 accessor);`min_sites` 为调用方参数(allel 同款 API,保留)。

**遗留(v0.2 待办,按优先级)**:spectrum/hwe/ld/fstats/allele_counts 接入 availability 门(θ_W 需拒绝 VO(Unknown));六种窗口级结果类型补 contig/start/end;`DetectedFacts` 扫描辅助;reader 的 NotDiploid/Multiallelic/OutOfOrder 错误路径测试;空窗口发射与多 contig 测试。

### 12.2 边界代理审计(2026-10-04,第二轮,对修复后代码复审)

8 项发现经本人验证(1/2 用构造实验,3 用 debug 构建)证实并**已全部修复**:

1. **`.` 样本列使整个流中止**(CRITICAL):VCF 规范合法的 `.` 样本列被 noodles 解码为 `Ok(None)`,decode_gt 误判为致命错误;BCF 的 end-of-vector GT 同理。修复:解码为 `./.` 缺失基因型(`reed-io/src/lib.rs`)。回归:`dot_sample_column_decodes_as_missing`。
2. **heterozygosity 的 p 分母口径不一致**(CRITICAL):分子按逐等位基因删减计 alt,分母却用完整基因型数×2 —— 半缺失基因型(如 `1/.`)存在时 hexp/Fis 静默错误(实测 S0=0/1,S1=1/. 给出 hexp=0 而应为 4/9)。修复:p = alt/(ref+alt)(`hwe.rs`),hobs 保持个体级。回归:新增断言。
3. **ld_r2_summary 空窗口 u64 下溢**(MEDIUM):k=0 时 `k*(k−1)/2` 在 debug 构建 panic(reader 按设计发射空窗口)。修复:k≥2 守卫(`ld.rs`)。
4. **Hudson 位点纳入规则与声明的 allel 语义偏差**(MEDIUM):原实现要求两组 n≥2 才计入 den;allel 语义是 den 只要 between 有定义(n≥1)即计入,num 才需双侧 n≥2。已对齐(注意 pixy 是"两侧都丢"的第三种语义,记录在案供将来对比)。
5. **WindowBuffer::new 未强制 end ≤ contig_len 与 pos 严格升序**(MINOR):已补 `EndBeyondContigLength`/`PositionNotAscending` 校验。
6. **非连续 contig 记录静默复用 ContigId 并重置网格**(MINOR):VCF 规范要求同 contig 记录连续;已改为显式报错 `ContigRevisit`。回归:`contig_revisit_is_rejected`。
7. **空窗口发射成本 O(gap/step)**(MINOR,设计取舍):巨大 gap + 小 step 会产生大量空窗口;v0.2 考虑 step 下限或跳过选项。
8. **PopulationView 与 window 的样本数不绑定**(MINOR,API 契约洞):错配时索引 panic 而非报错;v0.2 加 `validate_against(window)`。

另修:PBS 逐对 `pair_mpd` 三次重复调用合并为一次;`DStatResult` 分母守卫 `!=0`→`>0`;pbs.rs 陈旧 clip 文档。**驳回**的项:LD 下溢在 release 下被 0·MAX 抵消不算崩溃(但 debug 下成立,仍修);`.` 列报错误导性(修复后该路径已不存在)。

测试基线:37/37 全绿(8 套件),clippy 0 警告。
