# 01 复审：外部审查报告核实 + 缺陷修复规范

> 复审对象：工作区 HEAD（`reed-core` / `reed-io` / `reed-engine` / `reed-stats` / `reed-filter` / `reed` / `reed-cli` / `oracle`，约 3.8 万行）。
> 触发：一份外部分片审查报告（7 智能体并行 + 逐条读码）声称发现 `H1–H16`、`M1–M7` 及低危若干。
> 方法：**不信任报告结论，逐条读码复验**；复核行号以当前代码内容为准（报告行号有少量漂移）。
> 对 H4（`bcftools stats` 格式）与 H14（`--het` 半分型）额外用本机 **bcftools 1.24**、**plink2 v2.0.0-a.7LM**、
> **scikit-allel** 做**实机对拍**。判定只依据代码与实机结果，不依据报告转述。
> 状态：本目录不修改任何代码。
>
> 判定标记：**【成立】** = 确是真 bug；**【部分成立】** = 说对一部分（已指明对错）；
> **【不成立】** = 报告误读代码；**【新】** = 核实中新增、报告未列。

---

## 0. 总评

报告对**"静默失效"与"输出损坏"这一类缺陷的判断质量很高**：H1（`--mind` + BCF 写出损坏文件）、
H4（`--summary` 四连错）、H5/H6（LD 绕过过滤与掩码）、H7（报告模式丢过滤）、
H12（stdin + `--indep-pairwise` 读空输入）、H13（`reed index` 三问题）**全部成立**，
且 H4 的 IDD 比报告描述更严重。

报告对**统计公式的量纲与分母**（H3、H11，以及它标注"未复核"的 `f4_ratio` 分母）也判断准确。

但报告有一个**方向性错误**：H15 声称"13 个公式文件零对拍测试、全仓无 `tests/`、无 `gold/`、
无 parity 脚本、`reed_oracle` 零引用"，据此把"补 parity 门"列为第一优先级。
核实结果是**恰好相反**——仓库有 `crates/reed-stats/tests/oracle_parity.rs`（1e-12 逐位对拍）、
7 个 `tests/` 集成目录、`gold/` 下 9 个 parity `.py` 脚本、`reed_oracle` 在 5 个测试中被真实引用。
（该错误很可能源自一个按 gitignore 过滤的 `grep` 包装函数产生的假阴性；改用 `command grep` 复验即被推翻。）

因此本轮的修复方向不是"从零建数值护栏"，而是**填补具体覆盖面缺口**
（`bcftools stats` 渲染零测试、`gold` 管道只跑 SNP 而不读 PSC 第 9 列、LD 语义无 plink2 对拍）。

---

## 1. 核实范围与判据

- **成立**：代码与报告描述一致，且用户可达或语义确错。
- **部分成立**：主体属实，但报告的某个推论/后果被高估或需限定。
- **不成立**：报告误读，代码实际无此问题。
- **可达性**：优先判断"用户能否触发"。若缺陷仅存在于库 API 层而当前 CLI 不可达，降级并注明。

本轮**未**重新核实的部分（沿用报告自身的"已核实无误"清单，本轮抽取样本抽检，未发现反例）：
`unwrap` 零漏网、`unsafe` 仅 1 处、Tajima 常数、Hudson Fst 代数、窗口 bin 划分、引擎保序与 permit 释放。
其中 **bgzf 虚拟偏移与 tabix/csi 随机访问边界**是报告自认的覆盖缺口，本轮仅覆盖到 `H13`（`run_index`），
`VoMap::voff` 的非零路径、CSI 大位置边界**仍待补测**（见 06 §4）。

---

## 2. 缺陷清单（修复规范）

> 严重度：**P0** = 产出损坏文件 / 静默错误输出 / 崩溃；**P1** = 数值语义错误或健壮性缺陷；
> **P2** = 契约、文档、清理。
> 每条给 `file:line`、缺陷、后果、修法。修法中的"门"指 06-roadmap.md §1 的通用验收门。

| # | 级 | 位置 | 缺陷 | 后果 | 修法 |
|---|---|---|---|---|---|
| D7-1 | **P0** | `reed/src/filter/mod.rs:451-476`（无条件写 `BCF\x02\x02` 头）、`:380`（`sketch_out=Some`）、`:561-597`（sketch `EmitJob`）；`reed/src/filter/emit.rs` 全篇 | `--mind`（内存充裕→`SketchThenEmit`）或 `--mind --indep-pairwise`（强制 sketch）走 sketch 路径，header 已按 BCF 写出，record 却只经 `edit_vcf_line` 按 VCF 文本写出，**无任何 `sketch && output_bcf` 拒绝守卫**；`:504-509` 的 `let _ = header_level` 另证 header 级别计算被丢弃 | `reed filter … --mind -Ob/-Ou/-o x.bcf` 产出"BCF 魔数 + 头文本"后接 VCF 文本行的**损坏文件** | plan 阶段拒绝该组合（`Err("BCF output requires the streaming path; drop --mind or use -Oz")`），或让 sketch 路径接入 BCF 编码器；删 `let _ = header_level`。附：`-Ou` 经 CLI 会把 `compress_level` 置 0，故"header 被压成 level 6"经 CLI 不可达，仅隐藏 `--compress-level`+`.bcf` 推断路径有 |
| D7-2 | **P0** | `reed/src/filter/mod.rs:161-177`（`StdinSpill` 仅 `spec.mind.is_some()` 创建）、`:282-285`（pass0 无条件 `take`）、`:587`/`:635`（pass1 `open_pass1(replay_path=None)`） | `--indep-pairwise` 不要求 `--mind`，indep-only 时不建 spill；pass0 抽干 stdin 后 pass1 以 `replay_path=None` 重开 stdin 得 EOF | `cat x.vcf.gz \| reed filter - --indep-pairwise …` 报 `input is empty`（报错、非静默，但用法合法却失败） | `indep_requested` 时也创建 spill；或在 plan 里把"stdin + indep"标为需要 spill。门：stdin 溢出测试 |
| D7-3 | **P0** | `reed/src/index.rs:483-486`（`input::open` 不校验 BGZF）、`:89-100`（`VoMap::voff` 在 `blocks` 空时 `coff<<16`=0） | plain VCF / 普通 gzip 被 `input::open` 接受，`PlainChunkSource` 的 `blocks` 为空 → 每条记录虚拟偏移恒 0，全程不报错 | `reed index plain.vcf` **静默生成全 0 虚拟偏移的无效索引** | `run_index` 校验输入为 BGZF/BCF（否则 `Err`）；门：`reed index` 对 plain/gzip 输入报错、对 bgzf 索引可用 |
| D7-4 | **P0** | `reed/src/filter/pass0/matrix.rs:226-262`（`survivor_body` 向 4 个基础平面写 header **全部**样本原始 `GtCode`，掩码只写 rule 平面）、`:888`/`reed/src/filter/mod.rs:411`（`n_samples` = header 总数）、`reed/src/filter/pass0/prune.rs:57-95`（只读 4 个基础平面） | `--indep-pairwise` 的 r² 与主等位频率在**全部输入样本的原始未掩码**基因型上计算，忽视 `-s/-S/--keep/--remove/--mind` 与 `--min-dp/--min-gq/--set-gts .` | LD 剪枝结果与 plink2「先 mask 后 prune」语义不符（子集数据上系统性偏差） | 构建矩阵时只写过滤后样本并应用掩码；或让 `site_prune_stats` 读 rule 平面叠加 keep 掩码。门：与 plink2 `--indep-pairwise` 逐 ID 对拍 |
| D7-5 | **P0** | `reed/src/stats/driver.rs:169-173`（剔除→`gt_groups[i]=None`）、`reed-stats/src/matrix.rs:69-78`（None 样本打包进尾部 slice 而非丢弃）、`:308-322`（`decay_range` 遍历含尾部 slice 并 `total.add_assign`） | 被剔除样本不写丢弃、照常进矩阵平面，且进 cohort 汇总 | `--ld-decay`/`--ld-pairs` 的 cohort（`pop=all`）曲线被被剔除样本污染（命名组不受影响） | 剔除样本不写入矩阵平面；或 `decay_range`/`collect_pairs` 跳过尾部 slice。门：`--keep` 子集的 LD 曲线 = 直接对该子集跑的结果 |
| D7-6 | **P0** | `reed-cli/src/translate/stats.rs:259-275`（`keep_samples`/`site` 已解析）、`:325-341`（`ReportRequest` 无这两个字段）、`:376`（报告分支 `return Ok(())`）、`:410`（`--sites-file` 在其后解析）；`reed/src/stats/writers/reports.rs:38-60`、`run.rs`(reports) 掩码硬编码 `SampleMask::all` | `--freq/--hardy/--het/--missing/--depth/--summary` 分支解析了 `-s/-S/--keep/--remove` 与 `--geno/--maf/--min-dp/-i/-e`，但请求结构体没有这些字段、下游也完全不应用；`--sites-file/--bed-file` 永远解析不到 | 报告模式的样本/位点过滤被**静默忽略**（数字覆盖全部样本与位点） | `ReportRequest`/`ReportJob` 增加并应用样本与位点过滤；把 `--sites-file/--bed-file` 解析前移；不支持的组合显式报错。门：`--freq -s` 与"先 filter 再 freq"逐数一致 |
| D7-7 | **P0** | `reed/src/stats/writers/summary.rs`：IDD `:397-408`、`:671-676`；nIndels 渲染 `:692`（`n_indels` 全仓仅 `report.rs:160/277` 定义与 merge，无自增）；QUAL `:552-588`；样本名 `:687`；PSC `:216-233` | ① IDD：`len<0` 时按 `-len` 次把 **0** 推入桶 0（桶内存的是**长度之和**而非位点数），`len>0` 时 `let _ = len;` 全丢，渲染 `:673` 桶下标 `-1` → 整节只输出一行假 `length=-1`；② PSC `nIndels` 恒 0；③ QUAL 缺失桶 `beg==0` 返回 `Some(-1.0)` 打 `-1.0` 而非 `.`（`None=>"."` 分支为死代码）；④ PSC 样本名硬编码 `sample_{i}`；⑤ **【新】** PSC hom/het 无 SNP 行门控（注释自称 SNP only，实现无判断） | `--summary` 输出与 `bcftools stats` 全面不一致（IDD 整节损坏、nIndels 恒 0、QUAL 缺桶值错、样本名错、het 计数含 indel） | 按 `vcfstats.c` 重写：IDD 缺失/插入**分离双桶**并按位点计数、自增 `n_indels`(`var_type&VCF_INDEL && gt!=HOM_RR`)、`beg==0` 打 `.`、注入真实样本名（入参改带 `sample_names`）、hom/het 加 `VCF_SNP\|VCF_REF` 门控。门：与 `bcftools stats` 逐行 md5 对拍（含 indel 与缺失 QUAL 夹具） |
| D7-8 | **P1** | `reed-cli/src/translate/stats.rs:829`（`split_roles(spec, 3, "--f4-ratio")`）、`:832-838`（`b=a`、`c=o` 硬编码）；`reed/src/stats/render.rs:160-163`（`to_args` 用 `..` 丢 b/c）；`reed/src/stats/request.rs:486/520/595-599`；`README.md:59`、`reed-cli/src/cli/mod.rs:277` | CLI 只收 3 群体，`b`/`c` 退化；引擎侧 `F4Ratio{o,a,x,b,c}` 按 5 群体算、文档/help 声明 `A,B;C,D;X` | f4-ratio 分母退化为 `f4(a,o,a,o)`，数值错误；`A,B;C,D;X` 输入直接报错（且报错示例串还是 f3 拼法） | 按 5 角色解析 `b`/`c`；补 `to_args` 的 b/c；对齐 README 与 help。附：CLI 可构造的请求往返其实等价，"往返不等价"仅对绕过 CLI 直接构造独立 b/c 成立 |
| D7-9 | **P1** | `reed/src/stats/request.rs:586-591`（F2/F3/F4 每站点返回 `(val, 1.0)`）、`:620-648`（`row()`）；`reed-stats/src/jackknife.rs:135/207-219`（`total_num/total_den`）；`reed/src/stats/writers/fstats.rs:23,36-40` | estimate = `num.value()`（Σnum，**sum 尺度**），se/z 来自 `den≡1` 的 jackknife（**per-site 均值**尺度），无重标定 | `estimate/se = n_sites·z ≠ z`，两列量纲差 n_sites 倍；`request.rs:544` 的 `/// estimate / se.` 注释对 f2/f3/f4 为假（z 列因尺度不敏感"侥幸"正确） | 统一到同一尺度：estimate 取均值（除以 `n_sites`）或 jackknife 用 sum；修注释。仅 f2/f3/f4 中招，Dstat/F4Ratio 已一致 |
| D7-10 | **P1** | `reed-stats/src/patterson.rs:353-356`：`let d = den.den.value();` | `f4_ratio` 分母应取第二个 f4 的**分子**（`den.num.value()`），却取了其 D 统计量分母；`reed/src/stats/request.rs:595-598` 自己实现时用的是分子，两处公式矛盾 | 任何调用 `reed_stats::f4_ratio` 的路径得到错误 f4-ratio | 改为 `den.num.value()`；并入 D7-22 的累加器合并。门：与 vcftools/allel f4-ratio 已知值对拍 |
| D7-11 | **P1** | `reed/src/stats/run.rs:341-344`（`n_alt>1` 跳过、`n_multiallelic_skipped+=1`、不入累加器）；`reed-stats/src/pi.rs:135-136`、`dxy.rs:89-97`（`absent = window_bp − n_present()`，`comparable += absent·total_pairs`） | 被跳过的多等位位点文件里**存在**，却不计入 `n_present`，被 `absent` 按"全调用单态"补进分母；`pi.rs` rustdoc 自称补的是"文件不含的位置"，与实现矛盾 | variant-only 下 π/dxy 分母系统性抬高、数值偏低（多等位越多偏差越大） | 多等位位点计入 `n_present`，或补分母时扣除；修 rustdoc。注意 θ 无此机制（其偏差来自分子丢失），报告对 θ 的推论不成立 |
| D7-12 | **P1** | `reed-cli/src/translate/stats.rs:153-157`（`file_keep.and_then(|p| std::fs::read_to_string(p).ok())…`）、`:168-173`（`--remove` 同构） | `.ok()` 吞掉读取错误 → `None` → `touched=false` → `Ok(None)`，不过滤也不报错；`#` 注释行未跳过；`--force-samples` 未接线（对照 `reed/src/filter/mod.rs:743` 正常 `map_err`） | `-S/-S/--keep/--remove` 文件名写错 → **静默不过滤** | 用 `map_err` 传播错误；复用 filter 侧 `keep_list_from_file` 语义（跳 `#`、接线 `--force-samples`）。门：缺失文件必须报错 |
| D7-13 | **P1** | `reed-cli/src/translate/stats.rs:767`（`build_windowed_request` 硬编码 `skip_empty: false`；对照 `:532` 用 `args.skip_empty`）；`reed/src/stats/sink.rs:181,310-343,405`（只有 `levels.genome`/`levels.contig` 门控，无 `levels.window`） | 窗口布局路径丢弃 `--skip-empty`；sink 从不检查 `levels.window`，窗口行无条件写出 | `--skip-empty` 失效；`--levels contig`/`--levels genome` 时仍写窗口行 | 传 `args.skip_empty`；sink 增加 `levels.window` 门控。附：Bp grid 缺 Regions grid 那样的 `levels.window==true` 守卫 |
| D7-14 | **P1** | `reed-stats/src/het.rs:59-71`（仅 `s0m && s1m` 跳过；`s0a==s1a` 判纯合）；`reed-core/src/genotype/gt_code.rs:201-209`（缺失槽 alt 位为 0）；`policy.rs:36-43`（默认 `AlleleLevel`） | `./0` 判纯合、`./1` 判杂合、二者都进 `n_called`；实测**等价 plink2 `--vcf-half-call reference`**，但未声明也无开关；与 `policy.rs` 文档"missing slot as missing"矛盾；`count/mod.rs` 的 `n_het` 与 `fst h_bar` 又是第三套读法；`--het` 走 `PloidyPolicy::HaploidAsHalfCall` 使单倍型也进入同路径 | O(HOM)/N 与 E(HOM) 分子分母采用不同倍性假设，不匹配任何 plink2 模式；同二进制内 `--het`/`--hardy`/`--fst` 对半分型三套解读 | 明确并统一半分型语义：`--vcf-half-call` 语义显式化（默认选 reference 或 missing 并写文档），`het`/`hardy`/`fst` 共用一种解读。门：与 plink2 `--vcf-half-call {reference,missing}` 逐列对拍 |
| D7-15 | **P1** | `reed/src/index.rs:464-542`（`is_bcf` 只用于选 `BcfIndexJob`，无 Tbi 拒绝）；对照 `reed/src/filter/mod.rs:485-490` 有拒绝 | `reed index -tbi x.bcf` 不被拦截，`Indexer::new(Tbi)` 一律建 VCF 格式 tabix | 写出格式为 VCF 的**假 tabix** 索引 BCF | 增加 BCF+Tabix 拒绝（复用 filter 的错误文本） |
| D7-16 | **P1** | `reed/src/index.rs:329-338`（`declared_total` 缺 `text.len() < need` 守卫，`:336` `&text[..need]` 越界）；对照 `reed/src/bcf_input.rs:314-321` 有守卫 | 截断 BCF 的跨块记录前缀短于 8 字节时 `&text[..need]` **越界 panic**，早于本应给出的"truncated"错误 | 截断 BCF → 崩溃而非干净报错 | 补 `if text.len() < need { return None; }`；门：截断 BCF 夹具必须报错不崩 |
| D7-17 | **P1** | `reed-stats/src/theta.rs:117-132`（`value_total = raw/n_sites`、`value_per_bp = (raw/n_sites)/window_bp`）；`oracle/src/estimators/spectrum.rs:168-171`（`(S/a_n)/window_bp`，无 `n_sites` 归一） | 默认路径多除 `n_sites`，与 oracle 分歧；`total` 字段名两义（站点均值 vs 窗口总量） | θ_w 默认口径与冻结 oracle 不一致（variant-only 路径反而与 oracle 一致） | 裁决：按 oracle 对齐（去掉 `/n_sites`）或保留 legacy 并改名；写文档说明两义 |
| D7-18 | **P1** | `reed-cli/src/translate/filter.rs:20-45`（用"是否含 `.`"区分 step 与 r²） | `--indep-pairwise 50 5 1` 中整数 `1` 被当 step 静默覆盖，循环结束 r² 仍 `None` → 报"needs an r² threshold"；整数 r²（含 1）无法作为裸 token 传入 | 合法用法 `50 5 1` 失败；`50 5 0.2` 可用 | 按**位置**解析（第三个 token 为 r²），允许整数；渲染端已有 `.0` 兜底 |
| D7-19 | **P1** | `reed/src/filter/pass0/matrix.rs:849-859`（`close_pending` 每 contig 立即 `thread::spawn`，`budget_bytes` 只查单个矩阵，`prune_jobs` 到 `finish()` 才 join） | 无并发上限、无信号量；M 个 contig 可同时驻留 M 份矩阵 | `--memory` 预算被突破，峰值内存可远超预期（注释自称"one matrix resident"与实现不符） | 用信号量限制 in-flight 矩阵数；门：峰值 ≤ `--memory` |
| D7-20 | **P1** | `reed-stats/src/fst.rs:148-156`（成对：`n1==0`/`n_c<=0` 两级守卫）vs `:284-295`（多群体：只 `an.contains(&0)`+`n_bar<=1`/`n_c<=0`，缺 `n_i==0`） | 多群体路径缺 `n_i==0` 检查；r=2 时被 `n_c<=0` 掩盖（现有测试只覆盖 r=2），**r≥3 且某群体 `an_i==1` 时**多群体路径把该位点计入而该群体对 s² 无贡献（半分型/单倍型必然产生奇数 `an_i`） | 多群体 Fst 在含单样本群体/单倍型的位点上偏倚 | 多群体路径补 `n_i==0`（或 `an_i<2`）跳过；门：新增 r≥3 且 `an_i==1` 的 `weir_cockerham` 位点用例 |
| D7-21 | **【新】P1** | `reed-stats/src/pi.rs:135-136`（`absent` 为窗口内未覆盖 bp，`absent * self.total_pairs` 无保护） | `total_pairs ≈ 2N²`，`absent·total_pairs` 无 checked 保护；`window_bp=1e9` 时 `n_group>9.6万` 即 u64 回绕 | 百万级队列 + 整染色体窗口时**静默**产出错误的 π（比 `pi.rs:57` 的 30 亿门槛低约 4 个数量级，且进入 variant-only 分母） | `checked_mul` 或用 u128 中间量；门：大 `window_bp`×大样本溢出测试 |
| D7-22 | **P1** | `reed/src/lib.rs:52`（`#![allow(unreachable_pub)]`）、`:63-92`（未重导出 7 类型：`IndexOutcome`/`OpenedInput`/`WindowRow`/`FStatRow`/`DecayBin`/`StatAcc`/`QcReportAcc`）；`StatsRequest::regions`（`request.rs:304` 死字段，全仓无读取点）；`reed-core/src/grid.rs:83-87`（`WindowGrid::default() = Bp{1,1}`）；`reed/src/filter/mod.rs:55-70`（`compress_level` 派生 Default=0，文档写"default 6"）；`request.rs:165`（`Levels::default()` 全 false → `to_args` 产 `--levels ""`）；`reed-filter/src/spec/render.rs:119-125`（`to_args` 丢 `min_af_kind`/`max_af_kind` 与 `ld_prune`） | 公开签名中出现不可命名的类型；`regions` 设了不生效；危险 Default：1bp 窗口切全基因组、compress 0、非法 `--levels ""`；`FilterSpec` 往返有损 | 库用户无法命名返回类型/无法设置 regions；`..Default::default()` 示例行为意外；round-trip 不闭合 | 重导出 7 类型（或降为 `pub(crate)` 并改签名）；删 `regions` 或接线；修正 `Default`（或 `#[non_exhaustive]`+显式构造）；补 `to_args` 缺字段。门：`cargo doc -p reed` 公开项自洽，往返测试 |
| D7-23 | P2 | `reed-stats/src/request.rs:327-334`（variant-only 分支 `{ .. }` 忽略 `per_bp`）；`reed/src/stats/writers/windows.rs:20-27` | variant-only 下 Pi/Dxy/Theta 一律用 `value_variant_only`，忽略 `per_bp`；同请求 `per_bp:false`+`true` 时 `theta_w.total` 与 `theta_w.bp` 同值、`pi.bp` 列装 per-pair 量纲 | 库 API 层标签与数值量纲不符 | **当前 CLI 不可达**（`render_tsv` 调用点硬编码 `variant_only=false`，窗口模式只产 `per_bp:false`）。裁决：按设计文档（variant-only 只有一种取值）删双 key，或按 `per_bp` 分派并测试 |
| D7-24 | P2 | `reed-io/src/error.rs:34,64-67`（`IoError::Io(String)` 经 `From<io::Error>` 压成字符串，路径/`ErrorKind` 全丢）；`reed-filter/src/error.rs:49,52-54`（`FilterError` 未实现 `source()`，`From<IoError>` 用 `io::Error::other(to_string())` 再压） | 错误链断裂，CLI `eprintln!("reed: {err}")` 看不到文件名 | 用户定位不了是哪个文件出错 | 用 `PathBuf` 字段承载路径（对照 `reed-core/src/group.rs:97,157` 的 `cannot open {}`）；实现 `source()` |
| D7-25 | P2 | `reed/src/filter/pass0/prune.rs:26-38`（keep 集键为 `(contig,pos0)`，无记录级区分）；`single_pass.rs:323-328`、`emit.rs:143-147`（整条记录按该键丢弃） | 同 contig 同 POS 的多条记录共享一次查询：在集合中则全留、不在则全剪；contig 不在 `per_contig` 时静默全剪（sketch 路径还不计数） | 同位点多记录（罕见）被一并剪/留；未知 contig 记录静默丢失 | 用 `(contig,pos0,record_index)` 或明确"同位点同处置"语义并文档化+计数 |
| D7-26 | P2 | `reed-stats/src/het.rs:13-16` 文档 | 文档写 `E_hom = N·(1−2p(1−p))`，实现是逐站点 `N − Σ_k 2p_k(1−p_k)`；且 `:13` 的 `Σ_k p_k(1−p_k)` 漏因子 2 | 纯文档记号草率（实现经 plink2 对拍数值**正确**） | 修文档；`het.rs:24-26` 字段级注释已准确，可作模板 |
| D7-27 | P2 | `reed-stats/src/pi.rs:57`（`n_chrom * n_chrom.saturating_sub(1)/2`，同形在 `pi.rs:69`、`fst.rs:103-104`）；`reed-stats/src/fst.rs:164,313`（`h_bar` 用 `n_het`，只计双槽全调用） | `:57` 需 `n_group>30.37亿`才溢出（理论）；`h_bar` 对半分型系统性偏低（与 D7-14 同根） | 前者物理不可达；后者含半分型数据上 Fst 的 `h_bar` 偏低 | 统一用 `checked_mul`；`h_bar` 与 D7-14 的语义裁决一起修 |
| D7-28 | P2 | `reed-filter/src/funnel.rs:236-262`（`merge` 未合并 `n_ld_pruned`）、`:224-235`（`reset` 未清零） | `merge` 丢字段（不被 `mem::take` 掩盖，聚合恒 0）；`n_ld_pruned` 全仓无读取者 | 接报告时 `n_ld_pruned` 静默为 0（当前是死数据，无可见后果） | 补 `merge`；或删该字段；若接线则加测试 |
| D7-29 | P2 | `reed-stats/src/*.rs`（13 个公式文件无内联单测）、`reed/src/stats/writers/summary.rs`（零测试覆盖）、`gold/scripts/run_bcftools.sh:19`（`VIEW_ARGS=(-m2 -M2 -v snps …)`）、`gold/scripts/lib.sh:57-69`（PSC 不读第 9 列）、`reed-core/src/grid.rs:139`（`bp_windows`/`sites_windows` 生产闲置）、`crates/oracle/Cargo.toml`（`version=0.2.0` 脱 workspace lockstep）、根 `Cargo.toml:49`（`assert_cmd` 声明未用） | 文档/测试漂移：`docs/rounds/v5/05-algorithms.md:19` 的"九节逐行一致"对 `--summary` 不成立（gold 未覆盖 indel/缺失 QUAL/PSC 全列）；`lib.rs:12,16`、`fst.rs:3-5` 的 parity 声明对部分量不成立；闲置代码与死依赖 | 声称的对拍覆盖与实际不符，易误导后续修改 | 扩 gold 管道覆盖 indel/缺失 QUAL/PSC 全列；给 `summary` 渲染与 LD 语义加对拍；删闲置/死依赖；oracle 版本并入 lockstep 或注明隔离。**注：oracle-parity 集成测试与 `gold/` 脚本确实存在（见 §3-1），此处只补覆盖缺口** |

---

## 3. 原报告中不成立 / 需更正的条目

### 3-1. H15【不成立】"统计公式零对拍测试；oracle-parity 门不存在" — 被直接证据推翻

报告称"13 个公式文件无一 `#[cfg(test)]`；全仓无 `tests/` 集成目录；`reed_oracle` 全仓 `.rs` 零引用；
无 `gold/`、无 parity 脚本、无 `.py`；唯一 1e-12 断言是自比"。核实结果：

- **`crates/reed-stats/tests/oracle_parity.rs` 真实存在**（约 14 KB），首行即声明是 "M4 acceptance gate"，
  真实 `use reed_oracle::estimators::*` 对 π/dxy/Fst/PBS/SFS 逐项比对，容差 `1e-12 * y.abs().max(1.0)`
  —— 这正是 `lib.rs:12,16` 声称的 parity 门。
- `find` 全仓有 **7 个 `tests/` 集成目录**；`crates/reed-stats/tests/` 含 `oracle_parity.rs`、
  `fst_multi.rs`、`het.rs`、`ld.rs`、`patterson_fstats.rs`、`sfs_joint.rs`、`accumulator_laws.rs` 等。
- **`gold/` 存在**（40 个文件），含 `check_parity.py`、`run_all.sh`、`compare_tools.py`、`run_allel_fstats.sh`、
  `run_pixy.sh`；全仓 **9 个 `.py`**。
- `reed_oracle` 在 `oracle_parity.rs` 与 `crates/oracle/tests/*.rs`（5 处）被**真实引用**；`1e-12` 全仓 70+ 处。

**只有三条子句成立**：13 个公式文件确实无内联 `#[cfg(test)]`；`assert_cmd` 确实未使用；
`jackknife.rs:285-288` 那条 1e-12 确实是自比。**结论应改写**：数值护栏存在但覆盖不均，
真缺口是 `--summary`（bcftools stats 渲染）与 LD 语义零对拍（并入 D7-29）。

### 3-2. H1 附带子断言【部分成立/需更正】

"sketch 路径用 `opts.compress_level`（默认 6）、`-Ou` 的 header 块被压缩"——经 CLI **不可达**：
`translate/filter.rs:374-377` 在 `-Ou` 时把 `compress_level` 置 0、`-Ob<level>` 时置为 level。
仅隐藏 `--compress-level` + `.bcf` 扩展名推断那一路（`mod.rs:451-461` 硬编码 6）才有不一致。主结论仍成立。

### 3-3. H2 的"to_args 往返不等价"【需限定】

对**任何 CLI 可构造**的请求，`b≡a`/`c≡o` 恒成立，`to_args` 渲染 `a;o;x` 再解析回来仍得 `b=a,c=o`，**往返等价**。
不等价只对**绕过 CLI 直接构造**、b/c 独立的 `FStatRatio` 成立。退化与丢字段本身成立。

### 3-4. M7-2 的"取整口径"【不成立】

报告称"多群体 `p_i` 用未取整 `an_i`、`n_i` 用 `an_i/2` 取整，奇数时口径偏"。核实：
成对路径（`fst.rs:151,158-159`）口径**完全相同**，且这正是 scikit-allel `weir_cockerham_fst`
的语义（`n = an // 2`、`p = ac / an`）的**刻意复刻**（`:242-243,:279` 注释自陈）。守卫不一致（缺 `n_i==0`）
成立（见 D7-20），取整口径不成立。

### 3-5. M4 需补充限定

"`crates/reed` 重复实现 reed-stats 的 f 统计累加器"成立，但并非完全不复用：
`request.rs:9-10` 复用了 `reed_stats` 的 `f2_site/f3_site/f4_site` 与 `patterson::PattersonKind`，
`request.rs:567` 复用 `BpBlockJackknife`。被弃用的只是 `F2Acc/F3Acc/F4Acc/f4_ratio` 这层**累加器包装**。
报告标注"未复核"的 `f4_ratio` 分母问题**确实成立**（已升为 D7-10）。

### 3-6. 低危中被推翻/需更正的条目

| 报告条目 | 判定 | 更正 |
|---|---|---|
| "文档引用仓库外文件（`docs/spec/architecture.md`、`doc_v3 04 §C7`…），仓库无 docs 目录" | **不成立** | `docs/` 在仓库内且内容充实（`docs/spec/`、`docs/rounds/v1..v6/`）；105 个 `.rs` 对 `docs/...` 的引用全部指向现存文件 |
| "proptest 死 dev-dependency" | **不成立** | proptest 在 `reed-core/tests/{genotype_planes,grid}.rs`、`reed-io/tests/genotype_scan_fuzz.rs`、`reed-stats/tests/{accumulator_laws,summary_merge_laws}.rs` 共 5 处使用 |
| "oracle 整 crate 全仓零代码引用" | **不成立** | `oracle` 是 workspace 成员（`members=["crates/*"]`），被 5 个测试引用；`publish=false`、`version=0.2.0` 脱 lockstep 成立 |
| "FunnelPartial 问题被 `single_pass.rs:326` 的 `mem::take` 掩盖" | **部分成立** | 行号漂移：326 是 `n_ld_pruned += 1` 增量点，`mem::take` 在 `:157/:179/:298`。`reset` 未清零确被掩盖；`merge` 丢字段不被掩盖，但该字段全仓无读取者（死数据） |
| "planes.rs:107 assert 是高危" | **不成立（已下调信息级）** | 该 assert 带 `/// # Panics` 文档，是调用方不变量守卫 |
| "reed-cli 未继承 workspace.lints" | **不成立（报告的"打回"正确）** | `reed-cli/Cargo.toml:21-22` 有 `[lints] workspace = true` |
| "`bp_windows`/`sites_windows` 生产闲置" | **成立** | 调用点全在 `#[cfg(test)]` |
| "criterion 仅 reed-io 有 bench" | **成立** | 根 `bench/` 是独立 workspace，无 criterion |

---

## 4. 核实中新增的缺陷

1. **【新】PSC hom/het 无 SNP 行门控**（并入 D7-7⑤）：`summary.rs:216-233` 注释自称 "SNP rows only"，
   实现无任何 class 判断；bcftools 用 `var_type&VCF_SNP||VCF_REF` 包住整段。含 indel 时 reed 的 `nHets` 偏大。
2. **【新】`pi.rs:135-136` `absent * total_pairs` 溢出**（D7-21）：比 `pi.rs:57` 的溢出门槛低约 4 个数量级，
   且进入 variant-only 分母，**静默**产出错误 π，非 panic。
3. **【新】fst `h_bar` 半分型口径**（D7-27）：`n_het` 只统计双槽全调用，与 allel `count_het` 不同，含半分型时偏低。
4. **【新】gold 管道盲区**（D7-29）：`gold/scripts/run_bcftools.sh:19` 的 `VIEW_ARGS` 只跑 SNP、无 indel 无缺失 QUAL，
   `lib.sh:57-69` 读 PSC 只取 `$3-$6/$14` 不读 `$9`；`docs/rounds/v5/05-algorithms.md:19` 的"九节逐行一致"因此不可信。

---

## 5. 防回归测试清单

每个 P0/P1 缺陷**先写出复现它的失败测试，再修复**（红→绿），测试名描述缺陷本身。

| # | 测试（夹具 → 断言） |
|---|---|
| D7-1 | `--mind -Ob` / `--mind --indep-pairwise -Ou` → CLI 必须报错（不得产出）；若实现编码则 `bcftools view` 可读回 |
| D7-2 | `cat clean.vcf.gz \| reed filter - --indep-pairwise 50 5 0.2` == 文件路径结果 |
| D7-3 | `reed index plain.vcf` / `.vcf.gz`（普通 gzip）→ 报错；bgzf → 索引随机访问命中 |
| D7-4 | `--keep` 子集 + `--indep-pairwise` 的剪枝 ID == 对该子集单独跑 == plink2（mask 后） |
| D7-5 | `--pop`（含被剔除样本）+ `--ld-decay` 的 all 曲线 == 对被剔除后的样本集合跑 |
| D7-6 | `reed stats --freq -s keep.txt` == 先 `reed filter -s` 再 `--freq`；`--sites-file` 生效 |
| D7-7 | 含 indel + 缺失 QUAL 的夹具 → `--summary` 与 `bcftools stats -s -` 逐行一致（含 PSC 样本名/第 9 列） |
| D7-8 | `--f4-ratio A,B;C,D;X` 解析成功且数值 == 手算 `f4(A,O;X,C)/f4(A,O;B,C)` |
| D7-9 | 输出满足 `estimate / se == z`（容差 1e-9） |
| D7-10 | `reed_stats::f4_ratio` == vcftools/allel 已知值 |
| D7-11 | 含多等位位点的窗口，variant-only π == 把多等位排除在窗口外重算的值 |
| D7-12 | `-S nonexistent.txt` → 报错（非静默） |
| D7-13 | `--pi --window-layout … --skip-empty` 输出无空窗；`--levels contig` 无窗口行 |
| D7-14 | 半分型夹具 → `--het` 与 plink2 `--vcf-half-call {reference,missing}` 逐列对应（并声明选哪一种） |
| D7-15 | `reed index -tbi x.bcf` → 报错 |
| D7-16 | 截断 BCF → `reed index` 报"truncated"不 panic |
| D7-17 | θ_w 默认路径 == oracle `spectrum.rs` 值 |
| D7-18 | `--indep-pairwise 50 5 1` 解析出 r²=1；`50 5 1 0.2` 报错 |
| D7-19 | 100 个 contig × 小预算 → 峰值内存 ≤ `--memory`（或 in-flight 矩阵数 ≤ N） |
| D7-20 | r=3、`an=[1,4,4]` 的 WC84 位点 → 与成对路径/allel 一致（跳过该位点） |
| D7-21 | `window_bp=1e9`、大 `n_group` → 不静默回绕（报错或 u128） |

---

## 6. 值得保留并推广的做法

- **oracle 冻结差分基准**：`crates/reed-stats/tests/oracle_parity.rs` 与 `crates/oracle` 是数值正确性的基石，
  本轮核实证明它存在且有效——后续新增公式必须接入它，而不是另建。
- **gold 管道的 note 字段**：`gold/scripts/` 的写法正确，只需把覆盖扩到 indel/缺失 QUAL/PSC 全列。
- **`# Panics`/`# Errors` 合同与 `expect` 的 invariant 注释**：本轮抽检的 panic/expect 均带合同，方向正确。
- **实机对拍纪律**：本轮的 H4/H14 结论之所以可信，是因为用 `bcftools`/`plink2` 实机对拍而非读码推断；
  建议把"引用参考工具输出"写入修复 PR 的验收证据。
