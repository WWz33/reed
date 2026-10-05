# reed 设计决策记录 (DESIGN.md)

> 来源:2026-10-04 grilling 会话 + 14 个参考项目源码调研 + 自有数据基准测试。
> 标记约定:【定】= 用户已明确拍板;【荐】= 建议(可改);【据】= 有当日实测/源码证据。

## 1. 定位

- 纯 Rust **库**(lib 优先,CLI/绑定后置)【定】。
- **API 库形态约束(2026-10-04 显式升格,全文生效)**【定】:
  - 一切设计决策先问"调用方拿这个 API 怎么组合",再问"命令行用户怎么用"——后者置后甚至不做;
  - 公共 API 的承诺面最小化:配置对象是值类型(`Clone`/`Debug`),错误走 `Result` 且类型化(不许构造期校验缺位、运行期 panic);回调/迭代器风格优先于"库替用户做主"的默认值;
  - **不提供"一键全统计"式入口**——组合权在调用方(§4.2 同一原则的 API 投影);
  - 对外部工具语义的借鉴只取**算法内核与求值器 AST**(如 bcftools `filter_t`),不照搬其 CLI 接口形态;
  - 参考项目的取舍以此过滤:评估某特性时常问"CLI 工具需要它,库也需要吗"(例:输出格式/报告生成不属于库;过滤器/统计量属于)。
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

**拓扑正交原则(2026-10-04 补,详见 §15)**:popgen 统计量按数据拓扑分三类,禁止互相迁就——

1. **位点普查型**(π/dxy/Fst/SFS/…):分母是窗口 bp 或成对计数,定长窗口是其自然载体,走 `WindowBuffer`;
2. **路径延伸型**(EHH/iHS/XP-EHH/nSL/Garud H/ROH):衰减跨度由数据涌现(500kb~3Mb),**窗口装不下**(源码级裁决:allel iHS 全染色体单趟 ssl 流式扫描、无窗口概念,`selection.py:324` + `opt/stats.pyx:726`;vcflib iHS 整染色体 string `substr`,`iHS.cpp:179-253`,且其文档自认窗口法是"重叠+丢弃+去重"的补丁,`iHS.cpp:57`)——走未来的 `HaplotypeStream` 通道【定·源码裁决】;
3. **全基因组累加型**(PCA/GRM/KING/LD matrix):需要 O(N²) 样本空间外积的跨 contig 累加器——走未来的 `MatrixAccumulator` 层。

后端选型证据【据】:自有数据(348 样本 × 3,315,604 位点)

| 负载 | noodles | rust-htslib | 说明 |
|---|---|---|---|
| clean.bcf 8GB(8 个 FORMAT 字段) | **110.4 s** | 129.1 s | noodles 惰性只解 GT 列胜出 |
| clean.vcf.gz 7GB(同字段) | **187.3 s** | 302.9 s | 同上 |
| clean_bg.vcf.gz(GT:PS 两字段) | 90.7 s | **49.2 s** | htslib 块级解包便宜,反超 |

机制:htslib 一碰 FORMAT 就解包整块(`bcf1_t` 注释明说 lazy=按记录块);noodles 按列解。popgen 负载(GT 为主)下 noodles 架构正确且实测不落下风。所有输出计数两库与 bcftools **逐位一致**(金标准校验模板已跑通)。

**⚠ 选型的字段访问模式边界(2026-10-04 补,§16.2 实测)**:上表结论仅在 **GT-only 访问**下成立。rust-htslib 读取路径每记录执行 `bcf_unpack(BCF_UN_ALL)`(`rust-htslib/src/bcf/mod.rs:213-218`),其耗时与调用方访问哪些字段无关;noodles 按需解列,每多解一列多付一列的钱。clean.bcf 实测(同日同缓存态):GT-only 时 noodles 110.5s vs hts 129.1s(noodles 快 14.5%);**GT+GQ 时 noodles 142.5s(+28.3%) vs hts 129.1s(零增量)——htslib 反超 9.2%,选型翻转**。因此任何"多解 FORMAT 字段"的特性(GQ/DP 过滤,§11 L1.5)立项前必须先重测选型;Reader trait 双后端抽象(§11 第一条)因此升级为该特性的前置依赖。

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
- **v0.2(2026-10-04 实验裁决后刷新,优先级据 §12.3;P0/P1/P2 已于同日落地,见 §12.4)**:
  - **P0 正确性修复**:PBS 的 `PairCounts.missing` 补记账(实验证实恒 0,含被 `continue` 跳过位点);`watterson_theta.value()` 空/全缺失窗口守卫(实验证实返回 `Some(0.0)`,违反 §4.4 铁律)。
  - **P1 结构统一**:九种结果类型补 contig/start/end(三种逐位点型补 contig);`PbsResult`/`LdSummary` 补 `value()`;`WattersonResult` 拆 total/per-bp 双口径;spectrum/hwe/ld/fstats 接入 availability 门(θ_W 拒绝 VO(Unknown))。
  - 原有:2-bit SWAR 物化(GCTA/plink2 位并行模式,一 u64 = 32 样本);区域并行(contig 级 rayon + CSI);流式 Iterator API。
- **v0.3(统计量扩展,优先级据 §14.3)**:f-统计无偏校正 + block jackknife(snputils BSD);**joint SFS**(allel `sf.py:203-330` MIT);**f4-ratio**(snputils `fstats.py:1222` BSD);**个体近交 F**(allel `hw.py:103`);**LD D/D'**(Lewontin 1964 论文重写,vcftools/plink2 只读参照);HWE mid-p(Graffelman & Moreno 2013);LD dosage 物化 + 整数算子(实测 1.2x,动机=确定性非精度)。
- **v0.4:选择族 = HaplotypeStream 通道**(§15):`ehh_decay` + `ihs` 先行(allel cython 核 MIT 可译,ssl 增量模式),XP-EHH/nSL/Garud H/单倍型多样性随后;ROH HMM(allel `roh.py` MIT)同层。
- **v1.0:MatrixAccumulator 层**:GRM(Yang 2010 论文重写)、PCA(Halko 2011)、KING(Manichaikul 2010 论文);密码子族(πN/πS,SNPGenie 论文规则重写,依赖注释层)。
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
- **PBS 三群体交集语义 vs fst() 成对语义的分歧**(2026-10-04 实验证实:同一数据 `fst(A,B)` 得 `num/den=-0.0833/0.5`,`pbs().fst_ab=None`):reed 用三群体 complete intersection 保证三支位点集对称,allel 用三对独立 `moving_hudson_fst` 最大化信息利用。取舍待定:文档化现行为,或提供 Intersection/PerPair 开关【荐:先文档化,开关随 v0.3 一并评估】。
- **基因型不确定性族(GPAT 系:bFst/pFst/pVst/hapLrt/gl-XPEHH)是否引入**:需 GL/PL 输入通道,与"v0.1 只支持 GT"的定位冲突;若引入是输入层级别的决策【未决】。
- **`standardize_by_allele_count` 分箱标准化的宿层**(iHS/nSL 后处理必经,allel `selection.py:1165`、vcflib `normalize-iHS`):需整染色体/大区段全量值,应作为独立后处理函数随 v0.4 选择族一并设计【未决】。
- **L1.5 过滤器层立项**(2026-10-04 定稿，架构见 [`docs/filter-architecture.md`](docs/filter-architecture.md)):过滤器层是真实空白(grep 证实全仓库无过滤基础设施)。**库定位下的动机重述**:reed 是 API 库,不是 CLI——若无过滤器层,调用方必须在脚本里手动预过滤再传路径,被迫产生中间落盘/子进程开销与"清单申报是否如实反映过滤"的记忆负担;库内过滤的价值是**消除文件侧效应、保持 ascertainment 声明与数据事实同源**(声明和过滤出自同一处代码)。**立项裁决(同日经历三轮实测与用户裁决闭环)**:(a) **noodles 单后端定稿**:noodles 0.117 自带 `MultithreadedWriter/Reader`(实测 348 样本负载 5.2 GB/s,4.6x 加速),写侧不依赖 htslib;读侧 GT+GQ 场景 noodles 慢 9.2%(§16.2)接受为既定成本,htslib 全线移除以换免构建环境/FFI 负担;(b) **过滤即改变 ascertainment**:升级为 Ascertainment v2 结构体(sites/invariants/genotypes/cohort 四维),计划展开时自动推导降级规则;(c) **LD pruning 为独立 ordered-carry pass**(§16.1-2);(d) **群体维度全线移除**:群体结构属 reed-core(PopulationView),过滤器层只认 cohort(keep-list 组合已覆盖);(e) **Emitter 默认 stdout-VCF**:Unix 过滤器语义,落盘格式为 CLI 便利层;(f) **记录流透传而非列式批次**:未触碰字段零重序列化,保真 by construction。金标准与对拍门见 `gold/`。

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

### 12.3 第三轮复审:实验裁决审计(2026-10-04,对评审意见本身的复核)

第三轮评审提出的 10 项意见经**可运行实验**复核(验证 crate 建于 `/tmp/reed-verify`,依赖 reed-core 公共 API,不触碰本仓库;HWE 基准用 Python `Fraction` 精确有理数独立实现 Levene-Haldane 递推作零误差对照)。裁决:**8 项证实、3 项推翻/修正**——评审意见本身也要被审计,性能直觉尤其不可信,一切以实测为准【据】。

**实验证实的项(待修复,v0.2 P0/P1)**:

1. **PBS 的 `PairCounts.missing` 恒为 0**(HIGH,§12.1 第 4 项的修复漏了 PBS 侧):构造 C 群体全缺失位点,`within_a/b/c.missing` 与三个 `between_*.missing` 全部为 0;部分缺失(`[0,-1]`)同样不累计。且 `pbs.rs:166-169` 的 `continue` 在计数前跳过,缺失信息被彻底丢弃——与 `fst.rs` 每窗口完整记录 `total−comps` 的行为不一致。§12.1-4 只修了 Fst。
2. **空/全缺失窗口的 `watterson_theta.value()` 返回 `Some(0.0)`**(HIGH):`n_sites=0, window_bp>0` 时返回 `Some(0.0)`(同窗口 pi(AccessibleBp) 正确返回 None,因 pi 检查了 `comparable>0` 而 watterson 只检查 `window_bp>0 && a_n>0`);全缺失窗口同病。违反 §4.4 "空分母→None" 铁律,下游会把"未测序区"误读为"严格单态区"。同族 tajima_d 反而无此问题(variance=0 → None),说明这是孤立漏洞。修法:`n_sites>0` 守卫。
3. **PBS 与 fst() 的位点集分歧**(HIGH,语义决策):C 全缺失时 `fst(A,B)` 单独算得 `num/den=−0.0833/0.5`,同数据 `pbs().fst_ab=None`(三群体交集为空)。原因:`pbs.rs:166` 用 complete intersection,allel 的 `pbs()`(`selection.py:1314-1319`)是三对独立 `moving_hudson_fst`。取舍已记入 §11 未决。
4. **元数据断层**(grep 证实):`SfsResult`/`WattersonResult`/`TajimaDResult`/`LdSummary`/`F2F3F4DStat`/`AlleleCountsResult` 九种窗口级类型无 contig/start/end;`HweSiteResult`/`HetSiteResult`/`LdPair` 只有 `pos:u32` 无 contig。
5. **`value()` 缺失**:`PbsResult`(`pub pbs`)与 `LdSummary`(`pub mean_r2`)是唯二无 `value() -> Option<f64>` 的窗口级类型。
6. **availability 门禁断层**:仅 dxy/pi/fst/pbs 四模块有 `availability()`;spectrum/hwe/ld/fstats/allele_counts 五模块无(与 §4.2 已知缺口一致,θ_W 应拒绝 VO(Unknown))。
7. **fstats 朴素无 Patterson 校正**(源码级):allel `admixture.py:38,84,141` 的 `h_hat=ac0·ac1/(an·(an−1))`、`f2=(a−b)²−ha/sa−hb/sb`、`f3=(c−a)(c−b)−hc/sc` 在 reed `fstats.rs` 中均无减项。
8. **WC84/Nei86 硬编码 r=2**:`fst.rs:168,178,202` 的 `(2.0−1.0)` 字面量与固定两项平均,K≥3 全局 Fst 无法直接泛化。

**实验推翻/修正的项(评审错误,以实测为准)**:

9. **LD "20~50 倍提升"严重夸大→实测 1.2x**(1000 位点×100 样本×10% 缺失,release×3 轮:214ms vs 172ms,物化本身仅 0.23ms):`dosage()` 实为两次数组索引+比较,并无直觉中的"解包开销";真热点是 O(k²·S) 浮点乘加本身,物化省不掉。正确预期 **1.2~2x**,动机应表述为**确定性+速度**(对齐 plink2 整数精确思路),非数量级性能。
10. **LD "灾难性精度抵消"不成立**:随机数据(等位频率 10–90%)下浮点与整数精确版 mean_r2 差 **1.7e-18**(机器精度级)。dosage∈{0,1,2} 方差天然 O(1),无大数减大数;极端近单态位点理论上可构造但实测未出现,且此类位点本被 `vx<=0` 守卫或微小方差路径处理。降级为"整数化顺带消除的理论边角"。
11. **HWE 平局浮点风险不存在**:263,577 个基因型组合网格(het 0–80,hom≤60,总样本≤140)对 `Fraction` 精确有理数,最大相对误差 **1.5e-15,零案例 >1e-9**;概率单调衰减使 `x<=obs_prob` 判断在浮点下完全稳定。原"加 1e-12 容差"建议**画蛇添足且有害**(会把应排除项错误纳入)。mid-p 建议仍有效(保守性是真问题,plink2 佐证)。
12. **HWE 堆分配成本需按 MAC 分层**(修正):3M 位点×348 样本,全 MAC 范围(mac 2..349)池化 scratch buffer 仅 **1.01x**(1111 vs 1099 ns/位点,分配不是热点,递推主导);小 MAC(mac 2..19,MAF 过滤后典型形态)**1.79x**(125 vs 70 ns/位点)。结论:池化是低成本值得做的改进(几行代码换小 MAC 场景 1.8x),但非热点,P0 之后做。

**校准后的修复优先级**(替代原 P0–P3 估级,P3 性能项全部按实测数字重标):

- **P0(正确性)**:PBS missing 记账;watterson 空窗口守卫;PBS 语义文档化(§11)。
- **P1(结构)**:元数据补齐;value() 补齐;WattersonResult 拆 total/per-bp(`WattersonResult::value()` 现除以 window_bp,而 `TajimaDResult.theta_w` 是窗口总量,两者单位差 L 倍,易误读);availability 门覆盖。
- **P2(口径)**:θ_W 拒绝 VO(Unknown);`PopulationView::validate_against(window)`(§12.2-8 遗留)。
- **P3(性能,实测校准)**:LD dosage 物化+整数算子(1.2~2x);HWE scratch 池化(小 MAC 1.8x);spectrum 族三函数合一次扫描(同窗口同群体连算 sfs/θ/tajima 有 3x 扫描浪费,方差字面重复计算 spectrum.rs:200 与 213-227)。

实验产物:`/tmp/reed-verify`(crate:reed-verify,bin:main/ld_bench/hwe_bench/hwe_bench_small + exact_hwe.py);金标准 venv `/tmp/allelvenv`(scikit-allel 1.3.13)仍在位。

### 12.4 v0.2 实施记录(2026-10-04,同日第五轮;P0/P1/P2 全部落地)

40 测试全绿(9 套件:unit 4 + core 8 + contracts 8(新增) + families 10 + unequal 3 + read_vcf 7),clippy 0 警告,全部金标准对拍不变(families/unequal 逐位一致)。版本 reed-core/reed-io → 0.2.0。

- **P0-1 PBS 审计线补全**:循环重构为"审计先行"——六个 PairCounts 对**全部**位点累计(含被 `n_a==0||n_b==0||n_c==0` 跳过的位点),`missing = chrom×(chrom−1)/2 − comps` 全额记账,与 fst.rs 行为对齐;回归:contracts.rs `pbs_audit_trail_covers_skipped_and_partial_sites`(C 全缺失位点:within_c.missing=6、between_ac/bc.missing=16;半缺失 A:3/3;并断言与独立 `fst(A,B)` 调用的审计线逐字段相等)。
- **P0-2 watterson 空窗口守卫**:`value()` 增加 `n_sites > 0` 条件(空窗口/全缺失窗口 → None,完整观测的单态窗口 → Some(0.0) 保留),对齐 §4.4 与 pi 的行为;回归:contracts.rs `watterson_reports_no_value_for_empty_or_all_missing_windows`。
- **P1 元数据**:九种窗口级结果(Sfs/Watterson/TajimaD/LdSummary/F2/F3/F4/DStat/AlleleCounts)补 `contig/start/end`;三种逐位点结果(HweSite/HetSite/LdPair)补 `contig`;回归:contracts.rs `window_level_results_carry_coordinates`。
- **P1 value() 补齐**:`PbsResult::value()`/`LdSummary::value()` 落地——至此全部窗口级结果类型提供 `value() -> Option<f64>`。
- **P1 双口径**:`WattersonResult` 拆 `total_theta()`(S/a_n,Tajima D 对照 π 的量)与 `theta_per_bp()`(allel windowed 语义);`value()` 保持 per-bp 别名;`window_bp` 字段改为方法(与 PiResult 一致);回归:total_theta 与 `TajimaDResult.theta_w` 逐位一致(0.7713498622589533)。
- **P1 门禁覆盖(逐一判定,§4.2)**:
  - `watterson_theta`:四态矩阵与 π·AccessibleBp 完全同构,`VO(Unknown)` **拒绝**、`VO(Filtered)` 有偏注记、其余可用——函数签名改为 `Result`;
  - **单态位点不变量**(推导+文档化):SFS/Tajima's D/LD r²/D-stat 的有效位点/位点对集合与单态位点存在与否无关(单态位点对任何 bin、S、π、num/den 贡献均为零),故 `VO(Complete)` 结果与全位点输入**逐位一致**(无偏),`Filtered`/`Unknown` 携带有偏注记;结果携带 `bias: Option<&'static str>`,无 Err 路径;
  - **f2/f3/f4 位点均值**:全位点语义下单态位点计入分母(和加 0、计数加 1),任何 VO 输入都改变量本身(位点密度依赖)→ 全部 VO 有偏注记;
  - **hwe/heterozygosity 位点局域豁免**、**allele_counts 普查豁免**:无跨位点分母,模块文档写明豁免依据。
- **P2 PopulationView 契约收口**:新增 `validate_against(window)`(SampleCountMismatch 错误)与非 panic 的 `get(group)`;新增 `samples_of` 边界助手,**全部 13 个估计器入口统一 `Result<T, Unavailable>`**——Err 三类:清单拒绝(仅 watterson)、群体/窗口样本数错配、组越界;深层索引 panic 类(§12.2-8)整体消除。
- **PBS 交集语义文档化**(§11 荐项之"先文档化"):pbs.rs 模块文档写明 complete intersection 与 allel 并集语义的分歧及后果;Intersection/PerPair 开关仍留 v0.3 评估。
- **Tajima 方差去重**(§12.1 降级项):variance 计算一次,守卫与结构体字段共享。

未随本轮:§12.3-P3 性能项(LD 物化+整数算子、HWE scratch 池化、spectrum 三函数合扫)、`DetectedFacts`、reader 错误路径测试、§11 L1.5 过滤器层(前置:后端抽象,见 §16.2/§11)。

## 13. 仓库工作流约定(2026-10-04)

- 本仓库(`reed-rs/`)即 git 仓库根,远端 `WWz33/reed`(私有,SSH over 443)。
- **内部分区(只提交、不推送)**:`DESIGN.md`、`docs/statistics-reference-map.md`、`bench/`——内部分析/审计/基准,留在本地提交历史。
- 分支模型:本地唯一分支 `dev`(日常开发与提交,**永不推送**);远端 `main` 只由 `scripts/publish.sh` 更新——它用临时 index 只装公共路径(`.gitignore Cargo.toml Cargo.lock README.md crates scripts`),带内部路径泄漏自检,force 更新远端 main(远端 main 是发布快照,不做增量历史)。
- 结构性防误推:本地无 `main` 分支,`git push` / `git push origin dev` 均不会把内部分区带上远端。
- 初次发布已把含内部文件的初始提交(3ea9b70)从远端引用上刷掉(远端现为 ed3970d 起始的干净子集;旧对象在 GitHub 侧待 GC,私有库无外部协作者,无泄露面)。

## 14. 主流统计族覆盖缺口全景 v2(2026-10-04,14 项目深挖后)

> 来源:对已下载的 14 个参考项目**逐一重读源码清单**(非印象盘点):allel 全部 10 个 stats 模块的公开函数、vcftools `parameters.cpp` 全部统计开关、plink2 `plink2.cc:145-184` 帮助块、vcflib 全部 148 个源文件的统计工具、snputils `fstats.py` `__all__`、GCTA `main.cpp` 命令面、pixy/SNPGenie/vcfstats 各自模块。结论:**"已实现"占比约 40~50%,四大数据拓扑空白之外又发现一批主流统计量缺口**。

### 14.1 覆盖现状总表

```
✅ 完整覆盖(金标准对拍通过)
   频率谱: AC/AF、SFS 1D(折叠/展开)
   多样性: π(双估计器)、θW、Tajima's D
   分化:   dxy(双估计器)、Fst×3(Hudson/WC84/Nei86)、PBS×2(Yi2010/PBSn1)
   HWE:    精确检验(Wigginton)、Hobs/Hexp/Fis(位点级)
   LD:     Rogers-Huff r²(未相位,成对剔除)
   f-统计: f2/f3/f4/D(朴素,无校正无刀切)

⚠️ 半覆盖
   f-统计: 缺 f4-ratio、Martin fd、Patterson 无偏校正、block jackknife
   LD:     缺 D/D′、相位模式(含三次方程解)、LD score、LD blocks(Gabriel)
   近交:   位点级 Fis 有,个体级 F 无(allel hw.py:103 / plink2 --het / vcftools --het)
   SFS:    1D 有,joint 2D 无(allel sf.py:203-330 全家)

❌ 整族空白(按数据拓扑,对应 §2 三分类 + 两个新识别维度)
   A. 全基因组累加型: PCA、GRM、KING、LD matrix      → MatrixAccumulator 层(v1.0)
   B. 路径延伸型: EHH/iHS/XP-EHH/nSL、单倍型 H1/H12、ROH、
      haplotype_diversity(moving_*)                 → HaplotypeStream 通道(v0.4,§15)
   C. 注释依赖型: πN/πS、dN/dS、MK 检验              → GFF/FASTA 注释层(远期)
   D. 基因型不确定性型(新识别): bFst/pFst/pVst/hapLrt/gl-XPEHH(vcflib GPAT 系)——
      全部基于 GL/PL 似然而非 hard-call;引入即输入层决策(§11 未决)
   E. 质控/边缘族: Mendel 错误(vcftools/plink2/allel 三家有)、单例计数
      (vcftools --singletons)、indv-freq-burden、LDhat/LDhelmet 格式输出、
      standardize_by_allele_count 分箱标准化(allel selection.py:1165,
      iHS/nSL 后处理必经,须独立于窗口的宿层)
```

### 14.2 本轮新发现的缺口清单(此前盘点遗漏,file:line 均经源码确认)

| 统计量 | 源码证据 | 地位 | 处置 |
|---|---|---|---|
| joint SFS(2D 联合频谱) | allel `sf.py:203-330`(folded/scaled/折叠缩放全家) | 群体历史推断(dadi/fastsimcoal)事实标准输入 | v0.3,MIT 直译 |
| f4-ratio(混合比例) | snputils `fstats.py:1222`(`__all__` 六导出之一) | AdmixTools 核心查询 | v0.3,BSD 直译 |
| 个体近交系数 F | allel `hw.py:103`;plink2 `--het`(plink2.cc:152);vcftools `--het` | 个体级质控标配 | v0.3,MIT 直译 |
| LD D/D′(Lewontin) | vcftools `--geno-r2/--hap-r2` 双模式;plink2 `--r-phased/--r2-phased/--r-unphased` | 重组分析标配(r² 受 MAF 强烈影响,D′ 反映无重组断裂) | v0.3,论文重写(Lewontin 1964) |
| LD Score | plink2 `--ld-score`(plink2.cc:160) | LDSC/GWAS 偏差校正输入 | 远期 |
| LD blocks(Gabriel 2002) | plink2 `--blocks`(plink2.cc:161) | 单倍型结构分析 | 远期,需相位 |
| 单倍型多样性 H1/H2 | allel `selection.py:738,768` | 单倍型族基础 | v0.4 随 HaplotypeStream |
| standardize_by_allele_count | allel `selection.py:1165`;vcflib `normalize-iHS` | iHS/nSL 后处理必经 | v0.4,宿层设计见 §11 |
| Mendel 错误 | vcftools `--mendel`;plink2 `--mendel`;allel `mendel.py` | 家系质控 | 视需求(非群体统计) |
| 单例计数 | vcftools `--singletons`(parameters.cpp:117) | 稀有变异负荷 | 视需求,小函数 |
| GPAT 系(bFst/pFst/pVst/hapLrt/gl-XPEHH/permute*) | vcflib `src/*.cpp`(bFst=贝叶斯 Fst 建模基因型不确定性;pFst=似然比;pVst=CNV Vst;hapLrt=单倍型长度指数 LRT;gl-XPEHH=GL 版 XPEHH) | 基因型不确定性路线 | §11 未决(输入层决策) |
| ROH 双家族 | vcflib/allel `roh.py:10,149`(multinomial/Poisson HMM);plink2 `--homozyg`(滑窗计数) | 近交/纯合片段 | v0.4 同 HaplotypeStream 拓扑 |
| PCA/GRM/KING | plink2 `--pca/--make-grm-bin/--make-king-table`(plink2.cc:163-169);GCTA 命令面;flashpca Lanczos | 结构三件套 | v1.0 MatrixAccumulator |

**刻意不跟进的版块**:关联分析族(GWAS/GLM/MLMA/PRS:plink2 `--glm/--score/--clump`、GCTA `--fastGWA*`)——属统计遗传学而非群体遗传学,reed 定位外;vcfstats 的 `Y~X` Lark 公式语法(`formula.py`/`macros.py` 三类装饰器注册)——**架构范本而非统计量来源**,未来统计量 API 扩展设计时参照(统计即插拔定义)。

### 14.3 落库优先级(结合许可证白名单 §8)

- **P1(纯函数、白名单源、现架构可挂)**:joint SFS(allel MIT)、f4-ratio(snputils BSD)、个体近交 F(allel MIT)、D/D′(Lewontin 1964 论文重写,vcftools/plink2 只读参照)、Patterson 校正 + block jackknife(snputils 五助手 BSD)。
- **P2(新架构层一:HaplotypeStream,v0.4)**:`ehh_decay` + `ihs` 先行(allel cython 核 MIT 可译);XP-EHH/nSL/Garud H/单倍型多样性随后;ROH HMM(allel roh.py MIT)同层。
- **P3(新架构层二:MatrixAccumulator,v1.0)**:GRM(Yang 2010 论文重写)、PCA(Halko 2011)、KING(Manichaikul 2010 论文)。
- **P4(视需求)**:基因型不确定性族(输入层决策)、质控边缘族(Mendel/单例,小函数)、密码子族(注释层,最重)。

## 15. HaplotypeStream 架构提案(衰减型统计 vs 固定窗口的权衡裁决,2026-10-04)

> 问题:EHH/iHS 族的单倍型一致性衰减跨度 500kb~3Mb,由数据涌现而非参数决定;固定 100kb 窗口装不下。裁决依据是**参考实现的源码行为**,不是直觉。

### 15.1 源码裁决:金标准们都不用窗口

- **scikit-allel iHS 无窗口**(`selection.py:324` + `opt/stats.pyx:726`):输入是**整条 contig** 的单倍型矩阵 `h:(n_variants, n_haplotypes)` + 全长 pos;前向一趟 + 反向一趟(`h[::-1]`)扫描;每到位点用**共享后缀长度 ssl 数组的增量更新**(`stats.pyx:761-796`:等位基因对逐对比较后 `ssl[u]+1` 或清零)维护"从当前位点向左的一致性";EHH 衰减由 `ssl2ihh`(`stats.pyx:434-533`)在位点索引空间梯形积分,`ehh<=min_ehh` 即截断——**衰减跨度是数据涌现的**。大 gap 用 `gap_scale=20000/max_gap=200000` 处理,`is_accessible` 数组处理参考序列空洞;**明确拒绝缺失**(`selection.py:41`,missing calls not supported)。
- **vcflib iHS 更极端**(`iHS.cpp:179-253`):每条染色体每个单倍型物化成 `std::string`,`calcEhh` 直接 `substr` 数相同子串,内存 O(单倍型×染色体长)。其文档自认窗口法是补丁(`iHS.cpp:57`):genome-wide 推荐重叠窗口+丢弃积分失败值+`sort/uniq` 去重——**这是对窗口架构的打补丁,不是设计**。
- **佐证一**:iHS 的标准化(按等位基因频率分箱,`standardize_by_allele_count` `selection.py:1165`;vcflib 配套工具 `normalize-iHS`)天然需要跨区域全量值——窗口装不下的第二证据。
- **佐证二**:allel `roh.py` 两套 HMM(multinomial `roh_mhmm`/Poisson `roh_poissonhmm`)同样是整 contig 输入。
- **并行性证据**:allel 用 ThreadPool(2) 前后向并行;各核心 SNP 的 EHH 积分互相独立,天然可按核心 SNP 分块 rayon 并行。

### 15.2 三方案对比

| 方案 | 做法 | 裁决 |
|---|---|---|
| A. 塞进现有窗口 | 窗口边界截断,或重叠窗口+去重(vcflib 式补丁) | ❌ 结构性错误:边界效应系统性低估积分;标准化需跨区全量值;vcflib 自己都要打三层补丁 |
| B. 无限大窗口 | window_size = contig_len | ❌ 内存 O(全染色体单倍型)(WindowBuffer 按 n_sites×n_samples×2 的 Vec\<i8\> 分配),违反"内存 O(窗口)"承诺 |
| C. **并行第二通道** | 固定窗口族不动;新增 HaplotypeStream(按需延伸的相位视图),选择/ROH 族挂其上 | ✅ 唯一与两个金标准实现同构的方案 |

### 15.3 设计要点(对照 allel 逐项可移植,【荐】)

1. **L1.5 层新增 `HaplotypeStream`**:与 WindowBuffer 平行的物化通道,**不做窗口切割**;保证 (a) 按 contig 流式读取;(b) 拒绝非相位输入(对齐 allel 前置检查;注意 reed 的 WindowBuffer 已有逐基因型相位位 `phased: Bits`,比 allel 的输入级检查更细,可直接复用判定);(c) 单倍型以位点为列、样本单倍型为行流式供给。
2. **核心算子照抄 allel ssl 增量模式**:O(n_pairs) 的 ssl 数组 + 前后两趟扫描,每位点 O(H²) 成对更新。成本量级:348 样本=696 单倍型=241,860 对,每位点 ~24 万次比较;全基因组 3.3M 位点 ≈ 8×10¹¹ 次操作,**必须并行**(核心 SNP 分块 rayon,各块独立)。
3. **窗口归属规则**(若需窗口级输出与窗口族 join):**核心位点归属**——每个核心 SNP 的 iHS 值归属于该 SNP 所在窗口;衰减积分可跨窗口边界,**落点不跨**。
4. **标准化独立成层**:`ihs_raw()` + `standardize_by_allele_count()` 两段式(allel 同构),标准化需整染色体/大区段全量值,宿层设计见 §11 未决项。
5. **缺失语义**:v0.4 与 allel 对齐(拒绝缺失输入,调用方先过滤);missing-aware 变体远期单独评估(allel ssl 对缺失是"assume sharing continues"的处理,`stats.pyx:775-778`,有偏,不宜默认)。

## 16. 外部提案评审记录(2026-10-04,"最终综合最优架构(无任何成本代价版)")

> 一次外部提案提出:四层 FilterPipeline(Global/Site/Genotype/Window)+ MultiScaleView(fixed/haplotype/roh/pca 四视图挂一结构)+ EstimatorRegistry(builder 注册)+ 统一 for_each_window_with_optimal_pipeline,宣称"零侵入、零拷贝、零额外开销、iHS 完美支持、所有 vcftools/bcftools 过滤器全部实现"。经两轮评审(第二轮为对第一轮批评本身的复审,含同日同缓存态实测)裁决如下。方法论延续 §12.3:**评审意见本身也要被审计——撤回、修正、加强都有发生**。

### 16.1 逐条裁决

**撤回的批评(第一轮说错,第二轮纠正)**:

1. ~~"SiteFilter::hwe 需要群体结构、存在次序依赖"~~——**撤回**。vcftools 的位点级 HWE 过滤就是全体样本口径:`filter_sites_by_HWE_pvalue` 直接 `get_genotype_counts(b11,b12,b22)`(`entry_filters.cpp:922,941`),无群体参数,cohort 级标准语义,在窗口物化之前执行完全正确。提案此项无毛病。

**修正的批评(结论维持,机制描述纠正)**:

2. **LD pruning 不能做成 per-window stage**——维持,但机制修正:plink2 `--indep-pairwise` 实际是 **bp 基滑窗 + 状态跨窗持续**(非"多遍迭代"):`LdPruneSubcontigSplitAll`(`plink2_ld.cc:2203`)仅在相邻位点间隔 > 窗口大小时切分 subcontig(条件 `:2238`)做并行,数学上与顺序扫描等价;窗口推进时 `cur_window_removed` 移除状态跨窗口携带(`:673`)。提案做成 `WindowFilter` 每窗口独立执行的真正问题:(a) 剪枝结果变成 I/O 网格(window_size/step)的函数——标准工具无此行为;(b) 不重叠窗口下跨界位点对从不被联合评估;(c) 移除状态不共享。正确位置:variant 流上的一遍独立滑窗 pass(plink2 同构)。

**实测加强的批评(从"待重测"升级为"实测翻转")**:

3. **"零额外开销"证伪 + 后端选型翻转**——见 §16.2。第一轮仅推断"需重测、可能是收窄",实测是**翻转**(软化错了方向)。

**维持的批评(第二轮复读证据后不变)**:

4. **FilterStage 谓词装不下三类过滤器**:`test(&self, item: &T) -> bool` 无法表达位点裁剪、基因级置缺失(GQ)、LD prune 返回子集——提案五个示例 stage 三个类型不成立;且提案自相矛盾:§3 说过滤在 `build_next_window` 之前(L1.5,层位正确),§2/§5 却是 `FilterPipeline<WindowBuffer>` 整窗谓词(L2,层位错误,一个 MAF 失败丢整窗)。
5. **Estimator registry 缺群体组合**:现有 `fst(&w,&pops,a,b,est)` 是(统计量,估计器,群体组合)三元组,`register_fst()` 无参无法指定算哪一对;`compute -> Option<O>` 抹掉 `Result<_, Unavailable>`(拒绝)与"算出但分母空"(None)的区分;统一返回类型抹掉各族精确计数结果(§4.3)与 bias 标记。提案全文未出现 PopulationView/SiteCompleteness/Ascertainment 任何一个词——对 reed 核心语义资产零感知。
6. **HaplotypeView 装不下 iHS**:`{buffer: Arc<WindowBuffer>, core_index, radius}` = §15.2 已否决的方案 A;radius 是构造期定值,EHH 衰减由数据涌现(allel min_ehh 截断);相位要求与缺失语义(allel `selection.py:41` 直接拒绝缺失)未处理。"iHS 完美支持"与自身类型定义矛盾。
7. **MultiScaleView 拓扑压扁**:PCA/GRM 是全基因组累加器不是 per-window view;ROH HMM 需 contig 级状态——违反 §2 拓扑正交原则。
8. **样板代码硬伤(未经代码库对表)**:示例 `for_each_window_with_optimal_pipeline` 的 registry 参数被 `EstimatorRegistry::default()` 遮蔽;`opts.filter` 声明未使用;`est.compute(&view)` 返回值丢弃(collector 从未接线,"完成所有统计"无处落地);`ReadOptions` 现为 `#[derive(Clone, Copy)]`(`reed-io/src/lib.rs:130`),加 `filter: Option<...>` 直接破坏 Copy。
9. **"所有过滤器全部实现"严重低估**:QD/FS/MQ/ReadPosRankSum 是 INFO 字段而 reed-io 从不触碰 INFO(grep 证实);GQ/DP 需扩 FORMAT 解码(§16.2 成本);多等位过滤与 v0.1 二元硬拒绝冲突。
10. **示例自我矛盾**:`VariantOnly(Complete)` 与 `.maf(0.05)` 并用——MAF 过滤丢位点后 Complete 声明即失实;过滤器必须携带清单降级规则,且 WithInvariants+位点过滤是三值枚举表达不了的中间态(单态位点还在、部分变异位点没了)。

**第二轮新发现(第一轮遗漏)**:

11. **多尺度聚合:多数族可行,PBS 是例外**。reed 计数式结果天然支持跨分辨率聚合(PairCounts、WC84 的 a/b/c、Nei 的 ht/hs、Hudson 的 num/den 均对外暴露且可加);**PBS 不行**——逐位点归一化 num/den 是函数局部变量(`pbs.rs:155,189-210`),`PbsResult` 只暴露 clip 后的 `fst_ab`、log 后的 `t_ab` 与最终 `pbs`,无法从结果重聚合。若做"多分辨率"(1kb~全基因组),先统一成分暴露(可与 §12.3-P1"补 value()/元数据"合并考虑)。
12. **样本级过滤(keep/remove/subset)是真空白**:`PopulationView` 要求每样本归组,全仓库无子集能力——提案的"Header 全局过滤器"点到了真需求,应吸收。
13. **"缓存 EHH 曲线"未定生命周期**:按核心 SNP 缓存全体曲线是 O(核心数×窗口) 内存,提案一句带过,不可行。

### 16.2 实测:字段访问模式决定后端选型(clean.bcf,同日同缓存态,2026-10-04)

clean.bcf(8GB,348 样本×3,315,604 位点)FORMAT 字段实测为 `[GT, AD, DP, GQ, PGT, PID, PL, PS]`——**GQ 在列**。工装:仓库原版 bench(`/tmp/bench-target/release/`,10-04 构建)+ 改装 bench-noo(/tmp/reed-bench-verify,增加按需解第二列,其 GT-only 模式与原版一致 111.0 vs 110.5s,证明工装无开销):

| 访问模式 | noodles | rust-htslib | 裁决 |
|---|---|---|---|
| GT only(当前 reed 负载) | **110.5s**(110.39/110.55 + 改装复跑 111.0/111.3) | 129.1s(129.10/129.15/129.26) | noodles 快 14.5%,维持 §2 选型 |
| **GT + GQ**(GQ 过滤所需) | **142.5s**(142.44/142.56/142.46) | **129.1s(零增量)** | **htslib 反超 9.2%,选型翻转** |

机制:rust-htslib 读取路径每记录 `bcf_unpack(BCF_UN_ALL)`(`rust-htslib/src/bcf/mod.rs:213-218`),全字段解码、耗时与访问无关;noodles 按列解、每列单独付费,GQ 列 +28.3%。n_extra=1,141,401,698 与 3,315,604×348 吻合,确认逐样本解码发生。

**结论与适用范围**:§2 的 noodles 选型**仅在 GT-only 访问下成立**(本机 clean.bcf 实测;VCF.gz 负载未测,历史数据 §2 表显示 hts 在 .gz 上本就更强)。任何多解 FORMAT 字段的特性(GQ/DP 过滤、GL/PL 统计)立项前必须先重测选型——Reader trait 双后端抽象(§11)由"观察项"升级为该类特性的**前置依赖**。

### 16.3 净判决与吸收清单

**总评(两轮一致,更精确)**:该提案是一份合格的**"过滤器层 + 多统计批处理"需求书**,不是架构定稿。宣称的"零成本/零侵入/iHS 完美/全过滤器免费"无一条经得起代码库核对(其唯一被实测的要件 GQ 过滤 +28.3% 且翻转后端选型)。

**吸收进路线图的部分**:

1. **L1.5 过滤器层立项**(提案最大贡献,§11 已登记):真实空白,但要 (a) 先解后端(§16.2);(b) 与 SiteCompleteness 打通携带降级规则;(c) 谓词/变换两类 stage 分开;(d) LD pruning 独立 pass;(e) 样本级 keep/remove 一并做。
2. **Builder 式统计请求 API**:`with_pi().with_fst(a,b)` 的声明式工效值得吸收,但按库形态约束修正:签名保留 `(统计量,估计器,群体组合)` 三元组与 `Result<_, Unavailable>`;静态分发(编译期可枚举群体对)优于 `dyn Estimator` 对象安全负担;请求集是**值类型**(Clone/Debug/serde 可选),错误在构建期以 `Result` 报告;**不提供"一键全统计"入口**(那是 CLI 的形态——库的用户在自己的循环里组合,默认值替他们做决定违背 §4.2"用户永远选不到错的组合"原则)。
3. **跨分辨率聚合**:以"计数式成分可加"为设计锚(§16.1-11),PBS 的成分暴露缺陷先修。

**驳回的部分**:MultiScaleView 四视图合一结构、FilterPipeline<WindowBuffer> 层位、EstimatorRegistry 统一 trait、HaplotypeView(radius 参数化)、"全部过滤器免费"的 scope、Block-Jackknife 塞进 per-window collector(刀切需跨窗口聚合层,与 v0.3 f-统计升级合并考虑)。

### 16.4 评审纪律备忘(累积三轮的教训)

1. 性能直觉不可信,必须实测:LD"20-50x"实为 1.2x(§12.3-9)、GQ"可能收窄"实为翻转(§16.2)——两次方向性错误都发生在未实测就下结论时。
2. 评审意见本身要被复审:三轮累计撤回 1 条、修正 2 条、降级 1 条;第一轮对本提案的 HWE 批评就是复审才纠正的。
3. "零成本"类宣称是可证伪承诺:拿最便宜的可运行实验去打它(本例:改装 bench 加一列解码,10 行代码)。
