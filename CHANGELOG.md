# Changelog

All notable changes to reed are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the
lockstep workspace versions, and the CLI flags plus the output file formats
are the stability promise (`docs/spec/cli.md` §0).

## [Unreleased] — doc_v23 复核修复：P0 崩溃/数据损坏与 6 项拍板落地（2026-11-03）

仓库外全库审查（`reed-rs-review-2026-11-03.md`，9 区域约 8.9 万行）经
`docs/rounds/v23/` 逐条实证后，8 条 Must-fix 中 7 条属实、1 条驳回；本轮把
A 组 6 条 P0/正确性项、两条 P1（NaN 排序、特征向量符号）与 6 项拍板全部
落地。每条先写失败测试再修复。

### Fixed
- **io/bcf**：`Cursor::typed` 与 `typed_int` 相互递归。描述符 count 半字节为
  15 时，旧实现把扩展描述符又送回 `typed`，连续 `0xF1` 字节（BGZF 压缩后
  不到 100 字节）产生数万层递归，758 字节的损坏 BCF 即 SIGABRT。现在按
  htslib `bcf_dec_size_safe` 只允许**一层**扩展，二次扩展、
  非整数/非单元扩展均为 `MalformedRecord`（V23-01）。
- **io/bcf**：CHAR 类型的 GT 让渲染与列解码读越界 panic（1 字节元素按 i32
  读 4 字节）。`fmt_gt` 对非整数类型返回错误（htslib `bcf_format_gt_v2` 同样
  报错），`float_element_column` 对非数值类型记 missing；`encode_gt_column_bcf`
  的类型守卫一并收紧（V23-02）。
- **io/vcf**：样本区编辑会越过 header 的样本数索引 keep/GT 掩码位集
  （`Bitset::get` 断言），行尾多一个空样本列 + `-s`/`--min-gq` 即 panic。
  `edit_samples` 现在以 header 样本数为上界，多余空列按扫描阶段的口径丢弃
  （V23-03）。
- **io/vcf**：`format_g6_into` 的科学计数分支用 `v.abs()` 格式化却不补负号，
  负 INFO/FORMAT 浮点在 BCF→VCF 渲染时被静默改成正数（`-1.5e-07` →
  `1.5e-07`）。现在指数分支补 `-`，与 `bcftools view -Ov` 逐字节一致
  （V23-04）。
- **filter**：表达式解析器无嵌套上限，10,000 层括号直接栈溢出 abort。加
  `MAX_NESTING`（按实测标定：debug 下 2 MiB 栈在 ~430 层崩、1 MiB 在 ~200
  层崩，取 128）与 `MAX_CHAIN`（2,000，约束左深 `||`/`&&` 树的 Drop 与求值
  递归）；超限返回错误（V23-05）。
- **core**：`round_half_even` 用 `EPSILON * max(|x|,1)` 判定平局，把相邻
  可表示值误判为平局：`round_half_even(1.4999999999999998)` 得 2（R 得 1）、
  `0.35*90` 得 32（R 得 31，DAPC 训练集大小 off-by-one）。改为精确
  `frac == 0.5`（V23-06）。
- **reduce**：`isomap`、`knn_from_distances`、`isomap_from_distances` 三处
  NaN 距离排序仍用 `partial_cmp(..).unwrap()`，缺失基因型即 panic；统一
  `total_cmp`（V23-08）。
- **reduce**：特征向量无符号约定，同一对称矩阵 fallback 得 `-0.998…`、faer
  得 `+0.998…`，切后端/升级 faer 会逐字节改 `.eigenvec`。`SymEigResult::top_k`
  与 `eig::SymEig::top_k` 统一"最大 |分量| 取正"（并列取首行、全零列不动），
  PCA/PCoA/kernel/Isomap/Laplacian/diffusion 共用（V23-09，拍板 A；接受输出
  符号变化）。
- **cli**：`stats`/`convert`/`pca`/`distance` 对 `-o -` 创建字面 `-.afreq`
  一类文件并退出 0，`distance` 还漏了空前缀检查；`popg` 的检查发生在整轮
  计算之后。四条命令与 `popg` 共用 `require_output_prefix`，在打开输入前
  拒绝 `-`/空/空白前缀（V23-21/V23-31）。

### Changed
- **stats（`--ld-pairs` 流式化）**：`LdRequest::pairs_path`（新字段）指向最终
  `PREFIX.ld` 时，sink 在 decay 作业 join（reducer 线程、contig 顺序）时逐
  contig 追加，首个 pair 到达才懒创建并写表头；全无 pair 仍不建文件，内存
  从 O(全部 pair) 降为 O(单 contig)。表头与行格式由 `writers::ld` 的
  `PAIRS_HEADER`/`render_pair_row` 单点定义，流式与内存两条路径逐字节相同；
  CLI 的 `--ld-pairs` 自动设置该路径；库调用方不设时行为不变
  （V23-14，拍板 C）。

### Documentation
- **spec/cli.md** §5：`AN=0` 时 `AF`/`MAF` 为 NaN 并按 IEEE 比较
  （`AF!=0` 为真），与 bcftools 一致，不改写为 missing；"与缺失值比较为假"
  仅指字面 `"."`（V23-50，拍板 A）。
- **oracle**：SFS/θ_W/Tajima 明示为**逐位点 complete-case**；记录与
  scikit-allel 逐位点删除语义的差异（allel 对同一输入 S 更大、bin 不同），
  并用同一 fixture 同时钉住两种口径（V23-47，拍板 A）。

### Tests
- 新增 25 项：BCF 二次扩展/非整数扩展/typed_int 边界、CHAR GT 函数级与
  record 级、尾部空列 × `-s`/GT 掩码、负科学计数与 INFO 编辑、表达式
  128/129 层与 2,000/2,001 项、`round_half_even` 近半值与 0.05…0.95×1…500
  扫描、isomap/kNN NaN GRM、`top_k` 跨后端与并列/全零规则、流式 `.ld`
  逐字节与无 pair 不建文件、五命令前缀拒绝（3 种前缀 × 5 命令 + 对照）、
  AF/MAF NaN IEEE 四条比较、oracle 缺失 golden 两项。
- 全量：`ci/test-parallel.sh` 893 passed / 0 failed（56 二进制，15 s）；
  `cargo clippy --workspace --all-targets -D warnings` 0 告警；
  `ci/check-doc-drift.sh` 97 项全 ok。真实 `subvcf.vcf.gz` 的 `--summary`
  DP/PSC 仍与 `bcftools stats -s -` 逐行相同；真实 BCF 的负浮点
  `reed filter -Ov` == `bcftools view -Ov`。

### 第二批：决策 A/B 与内存三项（2026-11-03）

拍板后继续落地的一批，明细见 `docs/rounds/v23/02-implementation.md`。

#### Fixed
- **filter（缺失值比较，V23-59）**：表达式路径的 `!=` 与 bcftools 方向相反
  ——`-i 'DP!=3'` 在 bcftools 保留 `DP=.`/`DP` 缺席的记录，reed 却丢弃。
  现按 `filter.c:4449-4455` 的 `missing_logic`：一侧缺失时仅 `!=` 为真，
  两侧都缺失时 `==`/`!=` 为真，排序算符恒假。三处入口统一：表达式
  `any_pair`、INFO field 规则的 `on_missing`（`pass` = 算符表，`fail` =
  显式丢弃）、`Number=A` 的 `.` 元素。`TAG=.` 与 tag 缺席同路径（V23-23）。
- **filter/cli（区域参数，V23-60/61；含 V23-19/20）**：`-r`/`-R` 与
  `-t`/`-T` 是两个 AND 的谓词（bcftools 交集），同族重复（`-r`+`-R`、
  `-t`+`-T`、多次 `-t`）取并集且与顺序无关。旧实现把四者塞进一个
  `include` 列表，`-r chr1 -R f` 只剩 `f`，`-r chr2 -t chr1:10` 忽略 `-t`。
  `RegionRule` 新增 `indexed`；`to_args` 把 `indexed` 渲染回 `-r`。
- **filter（驻留预算，V23-10）**：sketch pass0 旧检查逐个比较单矩阵与
  预算，多个各占大半预算的 contig 全部驻留。`ResidentBudget` 改为累计
  （`saturating_add`），拒绝的矩阵不计数。
- **distance（GRM 内存，V23-12）**：`matrices::grm` 不再复制 `n × m`
  的 `x` 缓冲（2000 样本 × 1e6 位点 = 16 GB），直接读 `st.w`，内积循环
  顺带累计 `called_both`；数值与旧算法逐值一致。
- **stats（LD 剪枝内存，V23-13）**：`LdPruner` 的 `pos_log` 在 bp 窗口
  推进后丢弃 `window_start` 之前的前缀（搜索从 `window_start + 1` 开始，
  死前缀不会再读），内存从 O(单 contig 位点数) 降为 O(窗口内位点数)。
- **reduce（`tql2` 不收敛，V23-15）**：QL 迭代 100 次仍未收敛时不再静默
  返回未变换的尾块，`sym_eig` 以带尺寸的断言终止；正常输入（2–5 次收敛）
  不受影响。
- **engine（看门狗顺序，V23-16）**：`FrameChain::wait_inner` 先取已就绪的
  cell 再判全局静默，慢存储恢复后不会把健康运行误报为 "chain stalled"。
- **distance（非有限剂量，V23-17）**：`standardize` 用 `is_finite()` 判断
  观测，Inf/-Inf 按缺失处理，不再让 `mean`/`sd` 变成 Inf/NaN 并污染整条
  标准化链。
- **reduce（核 PCA 负特征值，V23-18）**：核 PCA/GRM 核 PCA 的坐标缩放从
  `abs().sqrt()` 改为 `max(0.0).sqrt()`，数值误差产生的负特征值映射为 0。
- **stats（LD decay 缺省，V23-22）**：`decay_max_bp = None` 不再合成
  `[0, bin)` 一行；`LdDecayAcc` 的空 edges 即零 bin，pairs-only 库请求
  不写 `.ld.decay`。
- **stats（decay 分箱，V23-27）**：距离小于首条 edge 时 `bin_index`
  返回 `None` 并跳过，不再静默计入 bin 0。
- **stats（merge 不变式，V23-28）**：`StatAcc::merge`、`ContigBins::merge_bins`、
  `BpBlockJackknife::merge` 的形状断言改为无条件；release 下不再经
  `zip`/`_ => {}` 静默丢数据。
- **io（诊断与空 ALT，V23-29/V23-30）**：`MalformedRecord`/`MalformedBgzf`
  的 `reason` 改为 `String`；样本列不匹配报 `found/expected`；`count_alts("")`
  为 0。
- **cli（`--format`，V23-32）**：表类 `--type` 显式给 `--format` 是 usage
  错误；`--contig-lengths` 读取失败带文件名；`tree.rs` 的内部节点 id 上限
  文档改为 `2n−1`。
- **reduce（diffusion，V23-33）**：平稳模式判定改为相对阈值
  `1e-9 · λ_max`；`pca::backend()` 文档改为编译期常量。
- **reduce（PCA 的 GRM，V23-38）**：`pca` 的 mean-impute GRM 在 faer
  构建下走 `backend::gram_divide`（此前零调用的死代码），与标量
  `Grm::mean_impute` 逐值一致。
- **cli（样本来源校验，V23-40）**：`common::split_sample_sources` 成为
  `-s/-S/--keep/--remove` 冲突规则的唯一实现，`filter` 与
  `stats/popg/pca/distance` 共用，删除两份约 40 行的重复逻辑。

#### Fixed
- **filter（`to_args` 两字符算符，V23-24）**：`<=`/`>=` 的掩码边界不再渲染成
  `"= 10"`，provenance 命令行可重放。
- **filter（空 `sample_select`，V23-25）**：JSON 的空选择不再被当作 cohort
  变化，`-q/-Q` 保留 INFO AC/AN 短路。

#### Tests
- 护栏：`pruned_header` 保序/短 keep、IUPAC 15 个子集全表、`pretty_sequence(n≠10)`
  固定回退、WcMulti r=3 半调用组对拍 allel 1.3.13、F2/F3/F4 拆分合并恒等、
  faer 1/4 线程与 fallback 逐值对拍。
- 全量：`ci/test-parallel.sh` **931 passed / 0 failed**（56 二进制，15 s）；
  clippy `-D warnings` 0 告警；doc-drift 97 ok。

#### Changed
- **convert（`--reference` 流式，V23-11）**：FASTA 改为 `BufReader` 逐行
  解析，峰值从"整文件 + 全部 contig 的 Arc"降为"当前 contig 的短暂两份"；
  CRLF、空行、带描述的 header 行为不变。

#### Tests
- 新增 8 项：缺失值算符表（表达式/head/Number=A/双缺失）、区域组合
  （9 组 + 同族并集）、驻留预算 3 项单测 + 多 contig e2e、流式 FASTA
  （CRLF/空行/多记录）、GRM 四组逐值对拍。

## [Unreleased] — doc_v24 结构性清理（2026-11-04）

### Changed
- **io（FORMAT 双代缓存，V24-06）**：`FormatLayoutCache` 改为 current +
  previous 两代（各 256），满时轮换而不是整表清空；`intern` 返回带 epoch
  的 `LayoutId`，旧代句柄在 `get` 处显式失败，不会别名到新 layout。
- **stats（WC84 唯一公式，V24-05）**：pair 与 multi 路径合并到
  `wc84_site`（pair = r=2，multi = r=N）；pair 分量从 2× 归一为 multi
  尺度，`Fst` 比值不变；新增 `FstAcc::wc_components()`。删除 PBS 三重
  冗余守卫与 `projected_sfs` 的死函数。
- **filter（区域唯一解析，V24-04）**：`TagTable` 携带 header 的 contig 词表，
  `Region::parse_in` 先整串匹配再按最后一个合法区间切分；`filter` 的 `-r/-R`
  谓词与索引 seek 由引擎一次解析，CLI 不再重复解析；`stats`/`popg` 的
  `-r/-R` 同样 header-aware。`-r` 指向未声明 contig 现在报错（bcftools
  语义），`-t` 保持静默空选择。
- **io（死代码与 header 词法，V24-03）**：删除 `RegionChunkSource` 的
  死 `pending_block` 状态机、`OutputFormat`/`for_format`/`is_empty`/
  `eof_block`/`n_samples`/`clone_for_worker`/`PlainChunkSource::first_offset`；
  `##TAG=<...>` 解析改为引号/转义感知的 `AnglePairs`，`ID=` 不再在
  `Description` 字符串内被误匹配。

### Documentation
- **oracle 契约**：新增 `docs/rounds/v24/02-oracle-contract.md`，逐估计器
  固定位点集、半呼叫、缺失与 allel 分歧；WC84 全缺失与 hobs/hexp 半呼叫
  补命名 golden（V24-01）。
- **bench 协议**：新增 `docs/rounds/v24/03-bench-protocol.md` 与
  `bench/README.md`；`run_matrix.sh` 改为 warmup/交替顺序/中位数/JSONL/
  保留 stderr，并在 `n_rec` 不一致时失败；vs-noodles 的 BCF 分支改用
  与 VCF 分支相同的 `select("GT")`（V24-02）。

## [Unreleased] — doc_v22 真实 callset 的 I/O 兼容与深度回退（2026-11-03）

真实 callset 复核（`bsaseq/raw.bsa.vcf.gz`：4 样本、`GT:PL:AD`、无逐样本 DP）
暴露 4 类问题：深度列恒缺失、sites-only 文件读不了、单倍体 GT 无开关、尾随
空样本列整文件报错。明细与复现命令：`docs/rounds/v22/`。

### Added
- **io**：`IntSpec`（主键按标量读、可选回退键按逗号向量求和）；文本扫描器
  `GenotypeScanner::new_with_specs` 与 BCF 原生解码 `decode_format_columns`
  都按它解析。`DP` 的深度回退（记录无 DP 键 → 该样本 `AD` 之和）由调用方
  显式声明，IO 不再按字段名推断（doc_v22 D22-1）；htslib 文本口径：`.`
  跳过、向量尾截断、全缺失记 missing。
- **filter（显式深度策略）**：`GenotypeMaskRule.depth`（serde 默认 false）标记
  `--min-dp/--max-dp` 的深度语义；显式 `FMT/DP` 表达式保持字面 DP（header
  未声明 DP 时报错，与 bcftools 一致）。同一 plan 同时含两者时生成两列。
- **io**：`LineRef` 接受 8 列 sites-only 记录（`bcftools view -G` 形态）：
  `format_text()`/`samples_text()` 为空，`n_sample_cols() == 0`；`filter`、
  `index`、`stats`、`rename`、`convert` 的 9 列假设一并改为 8 列下限。
  sites-only 输出保持 8 列（VCF 与 BCF 输出）。
- **filter**：`--min-dp/--max-dp` 在 header 只声明 `FORMAT/AD` 时也可用；
  记录 FORMAT 带 DP 键时仍只用 DP（该样本 DP 为 `.` 不回退）。
- **filter**：`--ploidy require-diploid|allow-haploid`（默认
  `require-diploid`）。spec 层早已支持 `policy.ploidy`，但标志层写死默认，
  mtDNA / 雄性 X / 单倍体生物无法过滤；`stats` 原本就用
  `HaploidAsHalfCall`，现在 `filter` 也可选。
- **io**：样本区行尾多余的空列（`...\t0/1\t1/1\t`）按 htslib 语义忽略；
  非空的额外列仍是 malformed。
- **stats**：`--summary` 的 DP/PSC 与 `--depth`、`--hist dp` 使用同一深度列；
  bcftools `calc_sample_depth` 口径已对拍（只计深度 > 0，AD 向量全量求和）。
- **test infra**：`ci/test-parallel.sh`（一次 build，按测试二进制并发；
  libtest 的并行只在单二进制内，cargo 本身串行跑各二进制）；
  `gold/scripts/run_parallel.py`（gold 对拍脚本的独立 case 并发）。
  gold 脚本默认输入改为 `subvcf.vcf.gz/.bcf`，全量 `clean.*` 只显式使用。
- **io（R22-3）**：`LineRef::shape()`（`SitesOnly` / `FormatOnly` /
  `WithSamples`）与 `validate_sample_cols(expected)` 固化列形态契约；文本扫描
  器、`PCA`、`convert`、行编辑共用同一校验，不再各自数 tab。超出 header 样本数
  的空列仍忽略，非空额外列仍 malformed；9 列 FORMAT-only 记录在 header 声明
  样本时按 bcftools 语义报 malformed。
- **io/bcf（R22-2）**：`bcf/mod.rs` 增加读/写/适配模块地图；`render.rs` 只保留
  渲染，FORMAT 列解码拆到 `columns.rs`；`encode.rs` / `native_encode.rs` 改名为
  `encode_noodles.rs` / `encode_native.rs`。
- **test（R22-4）**：新增 `crates/reed-io/tests/format_columns_vs_bcf.rs`：
  500 个随机记录（FORMAT 布局、AD 宽度、`.` 与短向量、DP 有无）经 native BCF
  编码后，原生列解码、原文扫描、渲染文本扫描三条路径逐值相同；另有 2,000-case
  的 `LineRef` 列形态 proptest。

### Fixed
- **stats/summary**：`FMT/DP` 为 0 的基因型此前进入 DP 直方图与平均深度，
  bcftools 只计 `dp > 0`（`vcfstats.c:1119`）；修正后 DP/PSC 与 bcftools 一致。
- **stats/summary**：模块文档写的 "AD fallback: the maximum non-ref AD" 从未实现
  且口径错误；实际回退是 AD 之和。
- **io/bcf**：`n_fmt == 0` 的记录渲染为 FORMAT `.` + 每样本 `.`，与 htslib
  `vcf_format` / `bcftools view` 逐字节一致（原先输出空 FORMAT 列和多一个空
  样本列，htslib 回读报 `FORMAT '' is not defined in the header`）。真实
  `bcftools view -Ob` 生成的 `FORMAT=.` BCF 上，`reed filter -Ov` 的数据行
  与 `bcftools view` 逐字节相同。

### Fixed — V22-R2 第二轮评审（字段解析、索引与 stats 消费）
- **io/filter**：FORMAT 子字段索引从 `u8` 改为 `usize`（`IntLayout`、
  `Layout::gt_index/float_indices`、`gt_subfield_index`、`line_edit` 的 GT
  位置）。第 256 个及之后的 `GT`/`DP`/`GQ`/`AD` 原先被 `u8::try_from().ok()`
  静默当作缺失，掩码与统计随之 no-op；`filtered.rs` 的重复 GT 位置实现
  收敛到 `record_pipeline::gt_subfield_index`。
- **io/stats**：整数 FORMAT 解析统一为 `info::parse_info_int`（可选 `+`/`-`）
  再限制在 `0..=i32::MAX`：超出 i32、负数、非数字均记 missing，与 htslib
  对文本整数的处理一致；不再出现"饱和到 `u32::MAX` 恰好等于缺失哨兵"或
  `5000000000` 变成有效值的情况。`parse_int_sum` 与 BCF `int_sum_column`
  用 `u64` 累加并以 `clamp_int_sum` 钳到 `INT_MISSING-1`，AD 向量求和不会
  撞哨兵。`--hist` 的整数 token 与 `--summary` 的 INFO/DP 共用同一解析器，
  `+50` 在 filter/stats 两条路径同值。
- **stats**：`reports::validate_report_tags` 在 `run_stats_opened` 入口校验
  头部类型：深度消费者（`--depth`、`--summary`、`--hist dp`）需要 Integer
  `FORMAT/DP`（无 DP 时 Integer `FORMAT/AD`），`--hist gq` 需要 Integer
  `FORMAT/GQ`，`--hist dp-site` 需要 Integer `INFO/DP`；非 Integer 报 usage
  （exit 2）而不是静默按缺失。
- **cli**：`--hist dp` 的 header gate 接受 `FORMAT/DP` 或 `FORMAT/AD`
  （AD-only callset 不再被 CLI 拒绝，与 `--summary`/`--depth` 的深度口径
  一致）。
- **test**：新增 257-key 布局/扫描/GT 掩码测试、整数边界（`i32::MAX`、
  `i32::MAX+1`、`u32::MAX`、三个最大值求和、`+50`、负数）测试、文本/BCF
  边界对拍、`--hist dp` AD-only 与 Float DP/AD/GQ 的 CLI 契约测试；
  `parse_float` 的 `inf`/`nan` 语义（htslib `strtod` 兼容）由测试固化。
  Float NaN 哨兵与有符号整数列留后续轮（设计选项见 `03-review-pass2.md`）。

### Changed — 测试速度与最小夹具
- `crates/reed-io/tests/bcf_vs_noodles.rs`：本地差分改用工作区上级目录的
  `subvcf.bcf`（3,000 记录 × 348 样本，8.6 MB；`REED_SUBVCF_BCF` 可覆盖，
  缺失时跳过），并把整文件解压改为流式
  前缀解压、noodles 端改用惰性 `Record`（原 `clean.bcf` 8 GB 全量解压 +
  `RecordBuf` 使该单测约 400 s；现 4.5 s）。
- `crates/reed-graph/src/newick.rs`：`deep_tree_does_not_overflow_the_stack`
  直接构造 2000 深的 caterpillar，不再用 2000×2000 UPGMA 生成（23.5 s →
  0.01 s）。
- 新增 `ci/test-parallel.sh`：`cargo test --workspace` 逐二进制串行，它一次
  build 后并发跑全部测试二进制；本机 860 项 run 阶段 50 s → 15 s。
- gold 对拍：`gold/scripts/run_parallel.py` + `run_bcftools_af.sh` 的 40 个
  独立 case 并发（顺序 ~40 s → 5.8 s 墙钟 / 40 s CPU）；`run_pixy.sh` 的
  4 个独立 run 并发，`n_cores = nproc/JOBS`；`run_all.sh` 原本三工具已并行。
  `lib.sh` 默认输入改为 `subvcf.vcf.gz/.bcf`。

### 验收
- `cargo test --workspace`：867 个测试通过、0 失败；`cargo clippy --workspace
  --all-targets`：0 告警；`ci/test-parallel.sh`：867 passed / 0 failed（56 个
  测试二进制，run 阶段 15 s；doctest 由 `--doc` 单独跑）；
  `ci/check-doc-drift.sh`：97 项全部 ok；`cargo doc -p reed-io --no-deps`：
  0 warning。
- 真实数据 `raw.bsa.vcf.gz`（chr1:1-200000）：`--min-dp 10` 掩码 1111 个原调用
  基因型 + 174 个原缺失，与 AD 求和判据逐样本一致；`--summary` 的 DP/PSC 与
  `bcftools stats -s -` 逐行相同。
- `bcftools view -G` 的 sites-only 文件：`filter -Oz/-Ob`、`index`、`stats
  --freq/--summary` 全部通过，输出 8 列；`GT` 为单倍体的文件在
  `--ploidy allow-haploid` 下与 `bcftools view` 输出一致。
- 真实 `subvcf.vcf.gz`：`stats --summary --hist dp --depth` 的 DP/PSC/GDB
  与 `bcftools stats -s -` 逐行相同；`.dp.hist` 最小箱为 3（无 0 深度箱）。

## [Unreleased] — doc_v21 真实数据复审 + 修复（2026-10-28）

在仓库**真实**数据上复核了评审 `reed-rs-review-2026-10-27.md`（HEAD fa8131f）：
`subvcf.vcf.gz` / `clean.vcf.gz` / `clean.bcf`，本机构建重现每一条。评审的 10 条
must-fix 全部属实，另翻出一条评审判为"误报"的真 bug（M3）。11 条 must-fix +
12 条 bug 类 should-fix 全部先写失败测试再修。明细：`docs/rounds/v21/`。

### Fixed
- **io/region（M3，评审漏判）**：区域读（`-r/-R`）读到索引范围末端的整块，块跨行
  时留半条记录，引擎报 `input ends inside a record (truncated file)`。真实
  `subvcf.vcf.gz` 上 40 个区间 36 个失败（`clean.vcf.gz` 因块块行对齐 0 个，
  故评审测不出）。`ChunkBytes.truncate_to` 在索引结束虚拟偏移（记录边界）截断；
  范围终点恰在块边界时不再读该块。
- **io/region（M2）**：区域 tid 改从**索引自身名字表**解析（`.tbi` 与部分 `.csi`
  的 ref 顺序是首现顺序，非头顺序），仅在索引无名字表时退回头顺序。修正评审
  "CSI unaffected" 的说法：两种索引都会因 tid 错位失败。
- **pca/dapc/distance（M1）**：BCF 输入未过 `InputJob`，二进制被按 `\n` 切分而报
  `malformed record`。改与其余入口一致地包装。
- **pca/dapc（M11）**：`--memory` 门控未计 dense 的 `n×n` GRM（pairwise 两个），
  显式 `--solver dense` 会通过门控后 OOM；按路径分别计峰值（`2nm+n²` / `2nm+2n²`）。
- **filter（M4/M5/M7）**：`-r/-R/-t/-T` 的文件按行解析 BED（0-based 半开）/`CHROM
  POS`；裸位置 `CHROM:POS` 接受；选中零样本时记录去掉 FORMAT+样本列（8 列，对齐
  8 列头）并 warn 一次。
- **io/vcf（M6）**：空 CHROM 命中 `ContigTable` 的 last-hit 哨兵，`-W` 触发 790 GB
  分配中止；快路径加哨兵守卫。
- **rename（M8/S20/S33）**：`--clear-ids` 与 `--set-missing-var-ids` 互斥（原来
  静默覆盖）；空 `--set-id` 模板报错；空 `--add-meta`/`-s` 报错。
- **graph（M9）**：UPGMA 合并进较低槽位，对齐 R `hclust` 的 tie 处理（整数矩阵
  穷举 25/25 vs 旧规则 19/25）；更正错误注释。
- **core/stats（M10）**：`--pbs` 的组名与运行期用**同一解析器**（`--force-samples`
  丢组时不再错标或 panic）；删除第二个解析器 `names_from_tsv`。
- **cli/stats（S11–S22、S29、S30、S32–S34）**：`--hwe` 第二 token 校验、`midp`
  去掉内部码、`--maf/--mac` 越界、`--indep-pairwise` r² 范围、`-O` 级别、`-W`
  无有效目标、popg 四个静默丢弃的 flag、`--xval-rep 0`、空 `-o` 前缀、空窗 θ 打
  NA、零分母 D/f4-ratio 出 NA 行、`--snps-only=just-acgt` 真正排除非 ACGT、`--max-mac
  0`、`-S`+`--keep`、空值 no-op、`pca/distance` 补 `--verbosity`。

### Docs
- 更正偏差：`docs/spec/cli.md`（`--hwe midp`、`-W`）、`README.md`（`--hwe`、
  BED 区域、`-W`、`--max-mac 0`、`--indep-pairwise` r² 范围）；仓库外评审报告
  顶部加"Corrections"块（#2 的 CSI 说法、末端"不可复现"判定）；新增
  `docs/rounds/v21/`（README/01-review/STATUS/06-roadmap）。

### 验收
- `cargo test --workspace`：exit 0，839 个测试通过、0 失败。
- `cargo clippy --workspace --all-targets`：0 告警。
- 数值锚点仅 `reed-graph` s36 支持向量随 M9 变化（旧值即旧 tie 规则产物，测试内注明）。

## [Unreleased] — doc_v20 全库评审 must-fix 修复（2026-10-24）

11 个只读评审覆盖全部 crate（报告在仓库外 `../reed-rs-review.md`），14 个
must-fix 全部先写失败测试复现、再修复（测试名 `r1_*` … `r14_*`）。明细与
验证记录：`docs/rounds/v20/README.md`。

### Fixed

- **filter**：`-f/--apply-filters` 之前匹配整个 FILTER 列，`LowQual;LowMQ`
  被 `-f LowQual` 静默丢弃；现在按 `;` 拆分逐值匹配（bcftools 语义，R1）。
- **io**：空 FORMAT 列命中格式缓存哨兵 `(空, usize::MAX)` 导致 scanner
  越界 panic、整跑中止；现在按无 GT 布局处理（R2）。
- **index**：BCF 记录的 contig id 越界直接索引 `contig_names` panic（顺序
  路径无 catch_unwind）；现在报 job 错误（R3）。
- **cli**：`-O ''` 在 filter/rename/convert 三处 `split_at(1)` panic
  （exit 101）；现在报 usage 错误（R4）。
- **cli**：`-o` 与输入同路径时先读后截断，四个子命令（filter/rename/
  convert/index）实测 exit 0 并毁掉输入；现在拒绝并保持输入不动（R5）。
- **pca/distance**：畸形记录（列数不符、不可解析）被静默 `continue`，
  截断文件产出错误矩阵零提示；现在按 scanner 同款 `MalformedRecord` 报错，
  顺带把每位点死内存 `sites: Vec<(Vec<u8>,u64)>` 收成计数器（R6）。
- **stats**：投影 SFS 权重递推从左端锚定，`--proj ≳ 1024`（p≈0.5）时
  溢出为 inf 并静默丢弃位点质量；改为众数锚定双向递推，任意 n 不溢出，
  小 n 数值逐位不变（R7）。
- **reduce**：t-SNE 梯度用了未归一核（整体大 `sum_q≈O(n²)` 倍，输出不是
  t-SNE）；改为归一 `q_ij`（R8）。
- **distance**：单态（fixed-ALT）列把均值而非 0 写进标准化矩阵，给每个
  GRM 元素加 `4k/m` 并给 PCA 加秩一常数方向；改为中心化 0（plink2/GCTA
  约定，与 flashpca 的均值填充在 fixed-ALT 上有意分歧）（R9）。
- **pca**：stream 求解器静默忽略 `--missing pairwise`（auto 在大 N 时必然
  走到）；现在拒绝该组合并给出三个出路，模块文档中"已对 plink2 做 e2e
  校验"的不实陈述已删（R10）。
- **rename**：`--set-id` 早返回静默禁用 `--update-name`；现在互斥拒绝
  （`--set-missing-var-ids` 与 `--update-name` 本就兼容，保持可用）（R11）。
- **convert**：`--msa` 对重访 contig（乱序输入）游标倒退、位点重复发射，
  静默产出错位 FASTA；现在报 "input is not contig-sorted"，接缝处同位点
  last-wins 语义保留（R12）。
- **distance**：估计量族全程串行 O(n²·m) 且矩阵驻留 2–4 份、
  `--threads/--memory` 被忽略；现在按行 rayon 并行 + dosage 行预计算，
  任意线程数结果逐位一致（1 vs 4 线程 `to_bits()` 相等锁定）。流式分块
  与 `--threads` 全链路接入留作 backlog（R13）。
- **core**：`RegionSet::stab` 每位点线性扫描全部 BED 区间（自带 TODO），
  外显子 BED × WGS 实际跑不完；改为合并区间二分 + 成员精确过滤，
  等价性由旧线性实现做 oracle 的穷举测试锁定，200k 区间 × 2000 查询
  从分钟级降到毫秒级（R14）。

### Fixed（第二批：should-fix 中的 bug 类 16 条，同轮追加）

- **filter**：cohort 变化后（`--keep`/`--mind`/掩码），表达式
  AC/AN/AF/MAC/MAF 仍读全群体陈旧 INFO，与同 run 的 `--maf` 口径分裂；
  现按 `info_trusted` 门控从基因型重算（S58）。
- **reduce**：UMAP 学习率去掉 1/4 下限，照 umap-learn 线性衰减到 0（S28）；
  t-SNE 困惑度搜索非收敛行不再静默清零，`perplexity <1 或 ≥ n−1` 报
  类型化错误（`tsne` 变 `Result`）（S29）；三处 kNN 排序改
  `f64::total_cmp`，NaN 基因型不再 panic（S31）；`xval` 组 id 越界
  返回类型化 LDA 错误而非 panic（S32）。
- **graph/cli**：bootstrap 根节点不再恒标 100% 支持（跳过全 1 键，
  support 数量对齐 ape 的 n−2/n−3）（S36）；`--tree --bootstrap` 的
  discordance 重复矩阵改用用户的 `--mismatch-count/--error-probability`
  （修前恒默认参——Count 树配 Probability(40) 的 support）（S37）；
  同批把 `--threads` 经 scoped pool 接通 distance 估计量与 bootstrap
  （R13 收尾，1 vs 4 线程输出逐字节一致锁定）。
- **pca**：剂量矩阵分配前对照 `--memory`，超限报 usage 并指向
  `reed filter --indep-pairwise`（修前 5000×3M 直接 OOM）；`run_pca_opened`
  收敛进 `read_genotype_matrix`（S76）。
- **filter**：keep-mask 长度契约统一——record_pipeline 五处短掩码容忍
  改直接索引 + debug 断言，与 sketch 位平面路径口径一致（S68）。
- **stats**：WC84 单群 an∈{0,1}（如单个半调用）时 pair 与 multi 两条
  路径统一走 keep-c 分支（修前 pair 整位点丢弃，`--fst wc84` 与
  `--fst-global` 在该角落分歧；评审"multi 保留"的说法在 r=2 也不成立，
  一并修正）（S46）。
- **io**：FORMAT 整数越界由钳制改 missing 哨兵（对齐 INFO 路径与
  htslib）（S39）；损坏索引虚拟偏移越界报 `MalformedBgzf` 而非切片
  panic（S40）；POS 超出 BCF i32 范围与负 pos0 报 `MalformedRecord`
  而非静默钳制（S43）；符号等位 + INFO/END 的索引跨度改用 END（三处
  span 站点收敛为共享 `indexed_span_len`，区域查询不再漏跨 bin 的
  符号缺失）（S42）。
- **oracle**：PBS 模块文档改为与代码一致的 per-pair 规则（"三群交集"
  的说法不实；代码零改动，契约测试钉死）（S79）。
### Fixed（第三批：对前两批的复审修复）

- **reduce**：**撤回 R8 的梯度"修正"**——末因子应为原始核 `w` 而非归一
  `q`（对 KL 代价做中心差分验证：原始核误差 1e-10，归一核 0.58），前一版
  让梯度小了 `sum_q≈n²` 倍、几百样本时嵌入停在初始化量级；测试改为对代价
  函数做有限差分（旧测试是梯度自身的复写，恒真）。同批：困惑度搜索的
  beta 从 `1/mean(d²)` 起步（真实剂量距离下 `exp(-d²)` 下溢导致整行清零）。
- **io**：符号跨度改为 htslib 的 `get_rlen` 规则——INFO/END 与 ALT 形状
  无关（`ALT=.` 的 gVCF 块、多等位同样延伸），并补上默认单遍写入器与
  `rename -W` 两条漏掉的路径（此前只有 sketch 写入器正确）；BCF 的
  `pos0 = −1` 恢复为 htslib 的端粒哨兵（渲染 POS 0、原样传递），只有
  < −1 报错；FORMAT 越界警告不再自称 "INFO"。
- **graph**：**撤回 S36**——本机 ape 5.8.1 实测 `prop.clades`/`write.tree`
  保留根节点（Nnode 个值），前一版丢弃根并声称 ape 也丢弃是事实错误；
  支持向量恢复为每个内部节点一项，根恒 100%。
- **pca/filter**：`--memory` 上限移入逐块收集（超限即中止扫描，不再先攒
  满）；`stream + pairwise` 与 `--set-id`×`--update-name` 的拒绝移到读输入
  之前；`INFO/AC` 等带前缀的内置字段在 `--set-gts` 掩码表达式中编译期
  拒绝（此前静默求值为 Missing）。
- **filter/distance/convert**：`-f .` 只匹配整列为 `.`（对齐 bcftools 的
  filter id −1 语义）；dosage 提升缓存由 f64 改 i8（省 16 GB @ 2000×1M，
  数值逐位不变）；MSA 乱序检查覆盖窗口外帧；keep-mask 断言在 release 也生效。
- **安全**：stdin 落盘临时文件改 `create_new`（O_EXCL）+ 有界重试——
  修前共享主机上可预测路径可被预置文件/符号链接劫持并截断受害者文件
  （实测复现）（SEC-1）。

## [Unreleased] — doc_v18 随机数基础设施（2026-10-21）

`xvalDapc` 与 `aboot` 的共同前置，加上两个消费者本身。没有可分割的 RNG，
任何并行重抽样的结果都会依赖线程调度，违反 `architecture.md` §7.6 的
确定性契约。

### Added

- **`reed-core::rng::Rng`**：PCG32（O'Neill 2014），可分割。
  - `from_seed` / `from_seed_stream` / **`derive(index)`**：按**重复序号**
    派生独立子流，因此第 `i` 次重复在 1 线程与 64 线程下取到同一批随机数。
  - `below` 用**拒绝采样**而非 `% n`：模偏在置换检验里是系统性误差。
  - `shuffle`（Fisher–Yates）、`sample_n`（Floyd 算法，代价
    `O(min(k, n-k))`，`poppr::genotype_curve` 用的就是它）、`resample_n`
    （有放回）、`multinomial`（`poppr` 参数化 bootstrap 用的；负权重直接
    拒绝，否则会静默偏置）。
  - 与独立的 PCG32 参考实现逐位对拍，参考序列锁进测试。

### Changed

- `reed-reduce` 里两处私有 xorshift（`umap.rs`、`extra.rs`）迁到共享的
  `Rng`。UMAP/t-SNE/NMF/PPCA 的随机序列因此改变；它们的测试只断言结构性
  质（分离度、对称性、重构误差），没有锁死数值，所以全部通过。

### Added — `xvalDapc`（RNG 的第一个消费者）

- **`reed-core::math::r_compat`**：R 的 `round()`（银行家舍入）与
  `pretty()` 步长表。`xvalDapc` 的候选 PC 列表与训练集大小都依赖这两者，
  而 Rust 的 `f64::round` 是 half-away-from-zero，不等价。
- **`reed-reduce::xval`**：按群体分层抽样、每档 `n_rep` 次重复、以
  `sqrt(mean((success-1)^2))` 选最优（照抄 adegenet 的非标准 RMSE），
  并列取最大 `n_pca`。PCA 在全数据上算一次、不随重复重算——重算会把验证
  样本泄漏进主成分。
  - 每组训练集大小是 `round(training.set * n_g)`，**不夹取到 `n_g-1`**，
    与参照一致（该情形在 `training.set` 被压到 `(popmin-1)/popmin` 后不可达，
    但去掉夹取让算术与参照逐项相同）。
  - 随机基线是 **`1/ng`**，即 `adegenet` 标签置换基线的期望值；不是
    `Σ p_g²`（后者是按先验猜的准确率，在不等组数据上高估 17%）。
- **`reed pca --dapc --xval`**：`--xval-rep`、`--xval-training`、
  `--xval-result`、`--xval-seed`、`--xval-pca`；写 `PREFIX.dapc.xval`。

### Added — `aboot`（bootstrap 树支持率）

- **`reed-graph::bootstrap`**：按**位点**（列）有放回重采样、每次重算距离
  与树、统计每个分支在重复中出现的比例。有根树按 clade 计，无根树按
  bipartition 计（规范化到含 tip 0 的一侧，避免同一划分记两次）。
- **`reed distance --tree upgma|nj --bootstrap N`**：`--bootstrap-seed`、
  `--bootstrap-cutoff`；支持率写进 `PREFIX.nwk` 的节点标签。
- 与 `ape::prop.clades` 的确定性对拍：固定三个重采样索引，计数
  `3 3 3 2 1` 逐项一致。

### Verified

| 检查 | 结果 |
|---|---|
| PCG32 参考序列（seed=42，前 5 个 u32） | 与独立实现逐位一致 |
| bootstrap clade 计数 vs `ape::prop.clades` | 逐项一致（3/3/3/2/1） |
| 候选 n.pca 列表 vs `pretty()` | 7 个范围逐项一致 |
| `round(0.9*n)` 银行家舍入 vs R | 一致 |
| 弱信号夹具成功率 vs adegenet | 同为"PC 越多越差"，量级一致 |
| 串行 vs 4 线程分片（256 次重复） | 逐元素相同 |
| 重复 `i` 是否受其它重复影响 | 不受（烧掉 1000 次抽样后仍相同） |
| `below` 均匀性（70k 次，7 桶） | 每桶偏差 < 5% |
| `sample_n` 位置均匀性（30k 次，10 选 3） | 每位偏差 < 5% |
| `sample_n` 抽样预算（20 选 5） | 恰好 5 次，不是 20 次 |
| `derive` 流间相关（2 万条流） | 相邻流 r = -0.0005 |
| `multinomial` 权重跟随（50k 次） | 每类偏差 < 5% |

## [Unreleased] — doc_v17b 判别分析与距离树（2026-10-21）

面向自然群体的一轮：收**有监督判别**（DAPC）与**由距离矩阵建树**
（UPGMA/NJ/newick）两族。

### Added

- **`reed-reduce::lda`**：Fisher 线性判别分析，照 `MASS:::lda.default` 的
  `method="moment"` 分支（`adegenet::dapc` 调的就是它）。
  - `fit(x, groups, levels, tol)` / `Lda::project(x, dimen)` →
    `LdaProjection { coords, posterior, assign }`。
  - 单边 Jacobi SVD（工作区原本没有 SVD）；不形成 `A^T A`，精度更好。
  - 与 `MASS::lda` 对拍：`means` 2.2e-15、`posterior` 4.9e-16、`assign`
    逐样本一致；`scaling` 按 `|dot|` 归符号后 1.0e-15。
- **`crates/reed-graph`**：新 crate，只依赖 `reed-core`（L3 纪律）。
  - `Tree`（ape `phylo` 形状）、`to_newick`（`ape::write.tree` 格式，
    显式栈，深度上限 10 000）。
  - `upgma`：`hclust(method="average")` 的 Lance–Williams 递推，含并列
    扫描序与存活槽（取较新簇）规则——两者都是逐位一致的必要条件。
  - `nj`：`ape::nj` 的 `C_nj` 语义（`Q=(n−2)d−S_i−S_j`，边数 `2n−3`，
    负枝长原样保留）。
  - 与 `hclust` 对拍：合并高度 4.9e-13；与 `ape::nj` 对拍：RF=0、
    13 条枝长逐位相同。
- **`reed pca --dapc --pop FILE [--n-da N]`**：PCA → LDA，写
  `PREFIX.dapc.ind`（`#IID POP ASSIGN LD1..LDk`）、`PREFIX.dapc.eig`、
  `PREFIX.dapc.var`。`--n-da` 默认 `#{svd² > 1e-10}`，照 adegenet 的
  `lda.dim`。assign 与 `adegenet::dapc` 逐样本一致。
- **`reed distance --tree upgma|nj [--digits N]`**：由已算好的矩阵建树，
  写 `PREFIX.nwk`（不写矩阵文件）。
  - `--tree` 只接受**距离**矩阵：对角线非零的相似度/亲缘矩阵（`grm` 默认值、
    `cov`、`king`、`pihat`、`phi`、`ibs`）被拒绝并提示改用 `1-ibs`。
  - `--digits` 是 **R `format()` 语义的有效数字**（默认 10），不是小数位；
    `0` 被拒绝（会把所有枝长舍成 0）。
  - UPGMA 的枝长按 `ape::as.phylo.hclust` 的口径折半：`hclust$height` 是
    合并距离（节点高度的 2 倍），不折半会让每条枝长变成正确值的 2 倍。
    超度量输入下树的 patristic 距离与输入矩阵逐对相等，有性质测试锁死。
- 门面：`reed::run_dapc` / `run_dapc_opened` / `DapcRequest` /
  `DapcReport`；`reed::Mat` 转出（`reed-graph` 与 CLI 共用）。

### Notes

- `docs/rounds/v17b/`（README / 01-research / 02-spec / 06-roadmap /
  STATUS）。核心决定：DAPC 的 PCA 沿用 reed 口径而非 `ade4::dudi.pca`
  的 `li = U·d`；两者每轴差一个常数因子，而 LDA 对特征空间整体缩放不变，
  故 `posterior`/`assign` 一致、`ind.coord` 仅尺度不同（01 §1.5 有推导）。
- 明确不做：`find.clusters`、`snapclust`、`dist.genpop` 的五种距离、
  全部绘图层、Bruvo/SSR/GenAlEx。（`xvalDapc` 与 `aboot` 后来在 doc_v18
  完成。）

## [Unreleased] — doc_v17 个体间距离矩阵族（2026-10-20）

新增 `reed-distance` crate 与 `reed distance` 命令；PCA 的 GRM 构建迁入
该 crate，成为所有矩阵的统一底座。

### Added

- **`crates/reed-distance`**：10 种个体×个体矩阵，每种复刻一个参照工具。
  - `allele_ct_distance` / `ibs` / `one_minus_ibs`（plink `--distance`）
  - `grm` / `covariance`（plink2 `--make-rel`、GCTA）
  - `king`（plink2 `--make-king`，KING-robust）
  - `pi_hat`（plink `--genome`）
  - `phi`（vcftools `--relatedness2`，Manichaikul）
  - `ajk`（vcftools `--relatedness`，Yang）
  - `discordance`（bcftools `gtcheck`）
  - 4 种输出格式：`triangle`（含对角线）、`square`、`square0`、`bin`
- **`reed distance`**：`--type`、`--format`、`--stand`、`--mean-impute`、
  `--unbounded`/`--nudge`、`--mismatch-count`/`--error-probability`。
- **`reed::read_genotype_matrix`**：VCF→剂量矩阵的共享桥接，PCA 与
  distance 共用。
- **`docs/rounds/v17/`**：三子代理的源码调研台账（file:line + 公式）、
  规范、路线、状态。

### Changed

- `Mat` 从 `reed-reduce::linalg` 下沉到 `reed-core::linalg`，打破
  `reed-reduce` ↔ `reed-distance` 的循环依赖。
- `Grm`/`Stand`/`Standardized`/`standardize` 迁入 `reed-distance`；
  `reed-reduce::pca` 改为 re-export，v16 API 不变（PCA 数值不变：
  HapMap3 λ1 = 26.200944934）。

### Verified

| 矩阵 | 参照 | 结果 |
|---|---|---|
| `grm` | plink `--make-rel` | max abs 5e-6 |
| `king` | plink2 `--make-king` | max abs 5e-7 |
| `allele-ct` | plink `--distance` | max rel 5e-6 |
| `phi` | vcftools `--relatedness2` | 915849 对，计数零差异 |
| `discordance`（count） | bcftools `gtcheck -E 0` | 逐对相同 |

### Known issues

- **`pihat` 未对齐 plink**：plink 用 PLINK 1.07 的 Nei 偏差修正，reed 用
  理论期望值。输出保留但标记未对齐。
- `discordance` 的 HWE 列、`--genome` 的 DST/PPC/RATIO 未实现。
- `ajk` 与 `cov` 未对拍。

## [Unreleased] — doc_v16 CLI 接入：`reed pca`（2026-10-20）

把 `reed-reduce` 接进命令面，新增 `reed pca`。

### Added

- **`reed pca in.vcf.gz --pcs 20 -o out`**：样本关系矩阵（GRM）的主成分。
  - `--stand binom2|binom`（默认 binom2）、`--missing pairwise|meanimpute`
    （默认 pairwise，对齐 plink2/GCTA）、`--solver auto|dense|stream`
    （auto 按内存预算选）、`-o PREFIX` 写 plink 兼容的
    `PREFIX.eigenvec`/`.eigenval`/`.pve`。
  - `--samples/-S/--keep/--remove` 样本子集；`--threads`/`--memory` 复用
    既有预算机制。
  - 门面新增 `reed::run_pca` / `reed::run_pca_opened` / `PcaRequest` /
    `PcaReport`（`run_pca_opened` 让 CLI 只开一次输入，stdin 可流式）。
- **`crates/reed-cli/tests/e2e_pca.rs`**：6 个用例——文件与表头形状、
  binom/binom2 精确 2 倍、无缺失时两种策略一致、样本子集、5 个用法错
  （exit 2）、dense 与 stream 一致（rel < 1e-6）。

### Verified

HapMap3（957×14266，k=10）与 plink2 逐位一致：

| stand + missing | reed λ1 | plink2 |
|---|---|---|
| binom2 + meanimpute | 26.200944934 | 26.2009 |
| binom2 + pairwise | 26.280628904 | 26.2806 |
| binom + meanimpute | 52.401889869 | = binom2 × 2 |

`cargo test --workspace`：59 个测试二进制全绿。

### Not included

`--missing em`、`--use-gl`、`sd/center/none` 标准化、`--window`（Local
PCA）、`--tw`（Tracy–Widom）、`--indep-pairwise`（上游预处理）。理由见
`docs/rounds/v16/06-roadmap.md` §P2 与 `docs/spec/cli.md` §11b。

## [Unreleased] — doc_v16 降维族（2026-10-20）

新增 crate `reed-reduce`（零外部依赖，3119 行，29 个单元测试），并完成
PCA 求解器与降维族的源码级调研。**未接入 CLI**；命令面接入见
`docs/rounds/v16/06-roadmap.md` P2。

### Added

- **`crates/reed-reduce`**：零依赖的降维原型库。
  - `pca`：binom2 标准化、GRM（mean-impute / pairwise）、稠密 PCA、
    scores/loadings。HapMap3（957×14266）上 λ1/2/3 =
    26.200944934 / 23.394076211 / 6.932494547，与 plink2 exact 逐位一致。
  - `lanczos`：矩阵-free 块子空间迭代（top-k），与稠密解特征值一致。
  - `graph`：PCoA、Kernel PCA（linear/RBF/IBS）、Laplacian eigenmaps、
    Diffusion map、Isomap。
  - `grm_reduce`：上述五种方法的 GRM 驱动入口（缺失策略在 GRM 构建时
    确定一次，全部方法继承）。
  - `umap`：fuzzy simplicial set + SGD，按 umap-learn 0.5.12 逐行实现；
    `find_ab_params(1.0, 0.1)` = (1.576943618, 0.895060710) 与参考实现
    差 < 1e-7。
  - `extra`：t-SNE、NMF、LLE、PPCA-EM（缺失值原生）。
  - `eig`：稠密对称特征分解（EISPACK `tred2` + `tql2`），以及**部分谱
    求解器**（Householder 三对角化 + Sturm 二分 + 逆迭代），只提取前 k
    个特征对。N=4000 时比稠密路径快 2×（170.6 s vs 334 s），λ 逐位一致，
    残差 < 1e-13。
- **`docs/rounds/v16/`**：README、01-research（12 条判决 + 产物索引）、
  02-spec（API 与精度契约）、06-roadmap（P1–P4）、STATUS。
- **`docs/pca-algorithms.md` §8–§15**：flashpca `stand` 五选项、跨工具
  对拍、GRM 角色、求解器族谱、速度/精度实测、纯 Rust 实现路线。

### Verified

- 与 numpy `eigh` 相对误差 4.5e-16（机器精度）。
- 与 plink2 / GCTA / plink 1.9 / flashpca 在同集合同缺失策略下数值一致
  （2e-6 = 文本位数极限）。
- 与 sklearn 对拍：PCA 2.91 vs 2.91、KernelPCA 1.80 vs 1.80、
  Isomap 1.93 vs 1.95、NMF 重构误差 0.1907 vs 0.1905。

### Known limitations

- 块子空间迭代在 HapMap3 上 `max_rel_residual ≈ 1e-2`（特征值已收敛，
  Ritz 向量残差未达生产级）；稠密路径是精确参照。
- 稠密路径 O(N³) 单线程（4000×20000 需 334 s），因零依赖约束无 BLAS。
- UMAP 用精确 kNN、t-SNE 无 Barnes-Hut、Isomap 用 Floyd-Warshall。

## [Unreleased] — doc_v15 缺陷修复（2026-10-18）

对一份全库评审（6 子代理静态评审 + 逐条实证复核，40 条发现）的修复；
规范与证据见 `docs/rounds/v15/`。

### Fixed

- **filter**：Float 类型 FORMAT 字段的基因型掩码不再整体失效（D15-1）；
  掩码与 `--set-gts .` 改变 cohort 后 `-q/-Q` 一律按掩码后基因型重算，
  不再读陈旧 INFO/AC,AN（D15-2，cli.md §3.2）；`--proj` 超过群体染色体
  数、pop 文件把同一样本分配进两个群体、`--levels` 空列表均报用法错并
  点名（D15-9/10/21）；`--bed-file` 与 `--step` 解析期互斥、`-o -` 写
  stdout 而非字面文件 `-`（D15-22/39）。
- **popg**：缺少 `-o` 在打开输入前即报错（D15-13）；stdout 管道关闭由
  panic 改为普通运行错误（D15-14）；裸 popg 不再忽略
  `--variant-only/--bypass-invariant-check`，不变位点检查按窗口路由同款
  执行（D15-15）；`--f2/--ld-decay` 不再静默把 π/dxy/θ 切到 variant-only
  口径（D15-18）；窗口发射改为"严格更后位点才 emit"，跨 BGZF 块接缝的**合法同 POS 重复**（VCFv4.3："It is permitted to have multiple records with the same POS"）落入仍开窗口，窗口和与 contig 行逐项一致（D15-19）；
  尾随三态记录之后的窗口正常输出（D15-20）；`--pbs` **可重复**：单三元组写 `PREFIX_pbs.txt`（逐字节不变），多三元组每份
  `PREFIX_pbs.<A>-<B>-<C>.txt` 各持自身表头（D15-8，命名规则见 cli.md §7）。
- **convert**：`--phylip` 组合运行不再吞掉 `--012/--hapmap/--msa` 的位点
  （D15-4）；`--012` 双等位门修正，三态位点按 vcftools 语义跳过（D15-5）；
  仅 `--nexus-binary` 时 `--outgroup` 生效且未知名报错（D15-6）；MSA 长
  等位不再以 u8 截断（D15-26）；`POS=0` 按 malformed 跳过（D15-27）；
  gzip hapmap 可走 stdin（D15-29）。
- **rename**：`--rename-chrs` 撞名报错而非写出重复 `##contig` 的非法
  VCF（D15-25）；`--idspace-to` 与 `--set-missing-var-ids` 互斥（D15-30）；
  `--rm-dup-ids` 改为两阶段（串行定序 + 批次并行 deflate：worker 不再做
  被丢弃的压缩、索引键随记录并行算好，有序 sink 按字节封口批次并行压缩）：
  3M 记录 `-Oz` 1.39s → 0.23s（对照 15.4× → 2.5×），压缩流与 sidecar 对
  任意 `--threads` 逐字节相同，decompressed 输出与改前逐字节一致（D15-28）。
- **engine**：输入源线程的 panic 与 worker/reducer 一样收敛到单一错误
  出口，报 `reader panicked: …`（`EngineError::ReaderPanicked`；此前该闭包
  不在 `catch_unwind` 下，panic 会绕过引擎的错误通道，D15-40）。
- **io**：BCF 头恰落 BGZF 块边界（头文本长度 ≡ 0 mod 65280）不再整体
  错位一字节（D15-3）；`##contig` 行的显式 `IDX=` 按 htslib 规则生效
  （D15-24）。
- **stats**：Hudson Fst/PBS 在窗口内无可用分子时输出 0.0 而非 NA
  （scikit-allel `nansum` 口径，D15-11）；oracle 的 PBS 与 f2/f3/f4 改为
  与 shipped/allel 同口径（D15-12/31）。
- **性能**：`--indep-pairwise` 的保留位点查询由逐 contig 线性扫描改为
  二分（1 万 contig、30 万位点单线程 ≈11×，D15-17）；`--rm-dup-ids` 的
  去重输出由"sink 串行重压"改为两阶段并行（3M 记录 `-Oz` 1.39s →
  0.23s，D15-28）；GT 浮点渲染 `format_g6_into` 改为栈缓冲零分配
  （242.9 → 182.9 ns/次，D15-35）；`sketch_size_bytes`
  补上 mask 位平面行，`--mind`+`--indep-pairwise` 的 sketch 分支加预算
  检查（D15-16，预防性）。

### Changed

- `reed popg` 不再接受 `--tmp-dir`（该命令不落临时文件）；`reed stats`
  保留该参数并在帮助中写明"单趟报告从不落盘，此参数仅为兼容接受"
  （doc_v15 D15-23，cli.md §2 同步）。
- 内部收敛：Patterson 只保留校正式（`PattersonKind::Naive` 与相关死结构删除，
  `f2_site`/`f3_site` 去掉 kind 参数）；QC summary 只剩一个实现
  （`summary::SummaryAcc` 再导出 `report` 的活结构，`QcReportAcc::merge` 委托）；
  fst/PBS 审计改用 `PairAcc::add_partial`，不再假装跟踪 missing（D15-38）；
  `notify_abort` 持锁通知关闭丢失唤醒窗口（D15-32）；`--indep-pairwise` 的
  保留位点查询与 convert/popg 的 `--proj` 校验收口到共享层。

### Changed

- `convert` 的输出缓冲改为两级落盘（doc_v15 D15-7）：单缓冲超 4 MiB 或
  sink 缓冲总量超 8 MiB 即落临时文件、渲染时流回；覆盖 hapmap 行、
  `--012` 矩阵与 `.pos`、phylo 列与 `used_sites`、MSA 帧流（渲染改为逐样本遍历，
  并按样本 seek 只读该样本切片——不做整帧重读，延续 v13 D13-2 的禁止重复 IO）。
  800k 位点/50 样本夹具
  单线程峰值内存：hapmap 123.4 → 15.6 MB、`--012` 70.8 → 26.0 MB、MSA 151.4 → 15.8 MB，
  输出逐字节不变；临时文件在成功、失败与中止路径均清理。MSA 的
  per-contig 缓冲仍驻留内存（留待单独设计）。

### Removed

- 死代码：`JointSfsAcc`（连同 `tests/sfs_joint.rs`）、`PattersonKind::Naive`、
  `Cursor::skip`、`SiteView::info_int/info_float`、`HeaderScan::n_samples`、
  `WindowGrid::bin_of_pos0`、`keep_mode` 死参、`header_level` 弃值、stats sink 的
  死 `or(take)`、convert 的空臂（doc_v15 D15-36）。

## [0.4.0] - 2026-10-17

The doc_v14 command-surface round: `reed stats` is split into the QC command
(`stats`) and the windowed statistics command (`popg`); `--hist` and
`reed rename` are new; `reed --help` groups its commands seqkit-style. The
split is byte-identical for every existing output
(`gold/popgen/run_split_parity.sh`, 22 cases).

### Changed

- **`reed stats` is now the QC-report command** (`--freq/--hardy/--het/
  --missing/--depth/--summary`); the windowed population-genetics statistics
  moved to **`reed popg`** (alias `popgen`; doc_v14 §2). The moved flags
  (`--pi --dxy --theta --tajima-d --fst --fst-global --pbs --pbs-normed
  --f2 --f3 --f4 --dstat --f4-ratio --block-size --ld-decay --ld-pairs
  --sfs --proj --sfs-int --variant-only --bypass-invariant-check
  --contig-lengths --window --step --levels --skip-empty --bed-file`) on
  `reed stats` are parse-time errors naming the exact rewrite (doc_v14
  D14-2) — never silent aliases.
- A bare `reed stats` (no report flag) is a usage error naming the flags.
- The library facade renamed its drivers to follow the commands:
  `run_stats` → `run_popgen` (`StatsRequest` → `PopgenRequest`,
  `StatsReport` → `PopgenReport`), `run_reports` → `run_stats`
  (`ReportRequest` → `StatsReportRequest`, `ReportsOutcome` →
  `StatsOutcome`). The library has no stability promise (cli.md §0.9).
- `reed --help` lists the commands in groups (Filtering, Statistics,
  Transform, File Operations, Miscellaneous); command paths are unchanged
  (doc_v14 §1).

### Added

- Expressions (`-i/-e`) accept negative number literals (`MQRankSum<-12.5`,
  the GATK hard-filter shape; bcftools' unary minus, `filter.c:3902`). A
  `-` before a digit starts the number; a leading `+` stays an error, as in
  bcftools, which lexes it as the addition operator (`filter.c:332`); `-`
  inside identifiers stays part of the name, and there is no arithmetic, so
  `A - B` remains an error.
- `reed popg`: the windowed statistics (π, dxy, Fst, θ, Tajima's D, PBS,
  SFS, Patterson f-statistics, LD decay), with the 0.3.x flags and
  byte-identical outputs.
- `reed stats --hist LIST`: histograms from the same pass (doc_v14 §3) —
  built-ins `qual`, `dp`, `dp-site`, `gq`, `idl`, `af`, plus generic
  `INFO/<TAG>` and `FMT/<TAG>` requests; one `PREFIX.<m>.hist` TSV per
  metric with exact mean/median/p95 for integer metrics.
- `reed stats --singletons`: the per-allele singleton/private-doubleton
  table in vcftools' layout (`CHROM POS SINGLETON/DOUBLETON ALLELE INDV`,
  `PREFIX.singletons`), with a POP column under `--pop`.
- `reed stats --pop FILE` (with `--force-samples`): stratified QC
  (doc_v14 §2.2) — `--freq`/`--missing` write `PREFIX.afreq.pop` and
  `PREFIX.vmiss.pop` (one row per population, first-appearance order),
  `--singletons` gains a POP column (`.` for a carrier outside every
  group; counting still runs over the whole cohort). The stratified values
  equal a per-population plink2 `--freq/--missing` run
  (`gold/popgen/run_popqc.sh`); the cohort files are unchanged. `--pop`
  without one of `--freq/--missing/--singletons` is a usage error.
  Multiallelic ALT frequencies stay collapsed (the `.afreq` convention).
- `reed rename` (alias `reheader`): identifier edits without decoding
  genotypes (doc_v14 §4) — samples (`-s/-N/--replace-sample-name`), contigs
  (`--rename-chrs`, `--chr-style plain|ucsc|mt-is-m`), variant IDs
  (`--set-id`, `--set-missing-var-ids`, `--update-name`, `--idspace-to`,
  `--rm-dup-ids`, `--clear-ids`), header lines (`-h -f --add-meta
  --remove-meta`); VCF byte editing with `-Ov/-Oz`.

### Fixed

- The `-q/-Q` stale-INFO warning fired on every run whose header declares
  `INFO/AC` and `INFO/AN` with an unchanged cohort (every GATK header does),
  asserting "-q/-Q reads INFO/AC,INFO/AN" of commands that never passed
  `-q/-Q`. It now fires only when a `nref`-form threshold actually takes
  the INFO path (spec cli.md §3.2: 取到时).
- `reed rename` wrote a blank line at every chunk boundary (the seam record
  arrives with its terminator and got a second one), so every output over
  one chunk carried stray empty lines; with `--rm-dup-ids` the same doubled
  newline made the dedup pass fail outright (`job failed: dedup: unparsable
  record`). The suffix probe no longer restarts at `_2` for every duplicate,
  so a run of identical IDs is linear instead of quadratic.
- `reed rename` ran its record planes on inputs whose CHROM is not declared
  in `##contig`: a contig-less VCF (legal VCF) and a sites-only 8-column VCF
  both failed mid-run. Plain output now passes them through, `-W` says what
  the header is missing, and the record planes are refused up front on
  sites-only input.
- `reed rename -h FILE -W` built the index from the *input* header's contig
  order while writing the replacement header, so region queries could
  resolve to the wrong contig; the index now follows the written header.
- `reed rename -f FILE` appended the missing `##contig` lines after the
  `#CHROM` line (not a valid VCF header); they are inserted before it.
- `--rename-chrs` and `--chr-style` overwrote each other instead of
  stacking (the style was applied to the pre-rename name and discarded the
  explicit rename).
- A CRLF header's last sample kept its `\r`, so sample renames counted the
  edit but silently did not apply it.
- `--replace-sample-name` used to accept regex syntax and silently do
  nothing (or strip the wrong text); the flag is documented as literal, so
  metacharacters are now a usage error.
- `--update-name` did not check that a mapped name exists (the one
  validation-table row that was missing): a typo renamed nothing and exited
  0; the run now fails naming the key.
- `--no-command` was parsed and ignored by `reed rename` (the
  `##reed_command` line was written anyway).
- `--hist` float tags were labelled one 0.1 step high (`iqual` divided by
  10 instead of `(iqual−1)/10`), so `SC=1.0` printed as the bin `1.1 1.2`;
  a missing QUAL became a real `0.0` bin instead of `#missing`; the `af`
  bins used the wrong denominator and dropped indels; integer metrics
  passed the overflow bucket off as an exact value (`DP 600` → row `501`,
  mean 501). Bins are now labelled in the quantity's own units, capped
  observations print `>max` and stay out of the exact quantiles, and only
  occupied bins are written (a single huge value used to materialize a
  half-million-row file).
- `--hist gq`/`dp`/`dp-site` did not check that the header declares the
  field and wrote an all-zero histogram with exit 0; an all-empty
  `--hist ""` list did the same. Both are usage errors now.
- `--hist FMT/<TAG>` counted one missing observation per record when the
  FORMAT key was absent (it is one per sample) and never split
  Number=A/G/R subfields, so `FMT/AD` looked empty.
- `StatsOutcome::hist` (the library's histogram output) was always empty.
- `reed popg --ld-pairs` without `--ld-decay` exited 0 and produced no
  output; it is a usage error now.
- INFO values survive the BCF round trip: an Integer outside the 32-bit
  range becomes missing instead of wrapping (`3000000000` used to be
  written as `-1294967296`; htslib and VCFv4.5 §1.3 both forbid the value),
  a Flag carrying a value keeps it (`DB=1` lost it, and `DB=0` — bcftools'
  "clear the flag" — came back as a set flag), and a String value of `.` is
  no longer turned into an empty string (doc_v12 D12-3/4/5).
- `-q/-Q` follow the genotypes once the sample set changed: after a static
  subset, a mask or `--mind`, the input's `INFO/AC,AN` describe a different
  cohort, so the decision is recomputed from the kept genotypes like
  bcftools does (vcfview.c:341-350); the promised stale-INFO warning is
  emitted on the INFO path (doc_v12 D12-6).
- `--mind` no longer evaluates a tail `-i/-e` expression blind: the sketch
  path has no record bytes, so `--mind` with such an expression takes the
  two-pass path, and `--mind --indep-pairwise` with one is a usage error
  instead of silently dropping every record (doc_v12 D12-1).
- Expressions see the header's tag types: `DB=1`/`DB=0` test a flag's
  presence (they used to match nothing), a multi-valued INFO tag compares
  any element (`DP4=0,0,30,40` at `DP4>10` used to miss), and a lone tag
  with no comparison is a compile-time error instead of a silent
  constant (bcftools evaluates it to false; doc_v12 D12-7).
- One INFO number grammar for the filter and the writer: `XI=1e3` no longer
  filters as 1000 while writing `.` to BCF, and an unparsable numeric tag
  is missing rather than compared byte-wise (`XI=abc` used to pass
  `XI>4`; doc_v12 D12-14).
- A GT-less record keeps its input `AC/AN` under a sample subset (it used
  to be recomputed from nothing as `AC=0;AN=0`; bcftools' `bcf_get_fmt`
  null check, vcfview.c:343; doc_v12 D12-2).
- `ALT='.'` counts zero ALT alleles, so `-m 2`/`-M 1`/`N_ALT=0` match
  bcftools on invariant sites; `-M 1` means "invariant only" instead of
  being dropped as "off" (doc_v12 D12-10).
- Records whose INFO/FORMAT counts do not fit the BCF header fields are
  refused instead of packed with the overflow leaking into the next field
  (65536 INFO keys or 300 FORMAT keys produced an unreadable file;
  doc_v12 D12-9).
- `--fill-tags` validates its list: an unsupported tag (AF, F_MISSING, a
  typo) is a usage error naming the supported set, like `bcftools
  +fill-tags`; it used to be a silent no-op that the default `-x` then
  turned into a deletion (doc_v12 D12-11).
- `[[site.info]]`'s `quantifier = "any"` requires an element to satisfy the
  comparison (it was true whenever an element parsed; doc_v12 D12-8).
- An INFO `AC` vector whose length does not match the ALT count falls back
  to the genotype counts (htslib's `bcf_calc_ac` check, vcfutils.c:56-64;
  doc_v12 D12-13).

### Added

- `reed convert` (docs/rounds/v13): one parallel pass over the input, any
  combination of targets. Phylogenetic matrices with vcf2phylip 2.9
  semantics (Ortiz 2019, DOI 10.5281/zenodo.2540861) — `--phylip`,
  `--fasta`, `--nexus`, `--nexus-binary` (SNAPP), plus `--outgroup`,
  `--min-samples-locus`, `--resolve-iupac` (deterministic, seed 0),
  `--write-used-sites`; the vcftools `--012` triple; the TASSEL hapmap
  (`.hmp.txt`); a reference-driven whole-site MSA with vcf2msa semantics
  (`--msa --reference`, `--region`, `--flank`, `--collapse-suffix`); and
  `--sfs`, byte-identical to `reed popg --sfs` (one shared request builder,
  one SFS fold). Outputs are independent of `--threads`; parity is gated by
  `gold/popgen/run_convert_parity.sh` (vcf2phylip five files, vcftools
  `--012`, TASSEL read-back, vcf2msa sequences). Deliberate divergences from
  the reference tools are registered in `docs/spec/cli.md` §8b: malformed
  records are skipped and counted, `--012` writes `-1` for half calls and
  haploids (vcftools writes 1), the hapmap `alleles` column lists the
  declared REF/ALT alleles, the MSA pads indels deterministically instead of
  calling an external aligner, and an allele set outside the IUPAC table
  skips the record where vcf2phylip crashes.
- `reed-core`'s `iupac` module: the IUPAC ambiguity table (62 entries), the
  binary-NEXUS genotype table, and the inverse lookup (`ambig_expand`: `R` →
  `AG`), the single implementation shared by every convert target.
- `reed convert --vcf`: the hapmap → VCF reverse direction (v13 C4), so the
  two formats convert both ways. The input is a hapmap (`.hmp.txt`, one- or
  two-letter cells auto-detected); `chrom`/`pos`/`rs#`/`alleles` are copied
  verbatim (`ID` = `rs#`, REF/ALT in declared order), each IUPAC cell expands
  to a `GT` (`R` → `0/1`), only `N`/`?` read as missing, and `alleles=NA`
  becomes an all-missing record. `-O v` (default) or `-Oz[level]` writes
  plain text or BGZF, `-o` omitted streams to stdout; the header carries no
  `##contig` (TASSEL writes none) and its `##reed_*` provenance lines are
  controlled by `--no-version`/`--no-command`. Exclusive with every
  VCF-input target and with sample subsetting; malformed rows are skipped and
  counted. On the lossless domain (biallelic diploid SNPs, no missing) a VCF
  → hapmap → VCF round trip reproduces every record, and TASSEL 5 reads
  reed's VCF back (`gold/popgen/run_convert_parity.sh` gate 6).
- `-g/--genotype [^](hom|het|miss)` (bcftools `-g`) for `reed filter` and
  `reed popg`: a post-mask site rule over the genotype classes
  (`TailRules::gt_class` → `TailRule::GtClass`). A half-call (`0/.`) reads
  as missing (bcftools `bcf_gt_type` UNKN), so `SiteCounts` gains
  `n_half_called`, counted by `count_site`, the plane sketch and the
  single-pass pipeline. Repeated flags AND (bcftools keeps only the last
  one); exact duplicates collapse and opposite polarities (`hom` + `^hom`)
  are a compile-time error. Funnel ids are `after_mask.gt_class.(hom|het|
  miss|^hom|^het|^miss)`. Diploid parity against bcftools 1.24 PASS over
  every call shape; haploid sites stay policy-gated (`require_diploid`
  rejects, `allow_haploid` counts them as half calls, so `-g hom`/`^miss`
  diverge from bcftools there by the same convention as `--maf`/`--mac`).
  A record whose FORMAT lacks the GT field skips every `-g` rule and
  survives (`SiteCounts::has_gt`; bcftools' `bcf_get_fmt` null check does
  the same) — distinct from a GT whose values are all `./.`. Class names
  are matched exactly, case-insensitively (bcftools `strcasecmp`; no
  surrounding whitespace is trimmed). Like bcftools, the rules evaluate
  after subsetting, over the kept samples only.

### Changed

- WC84 Fst (`--fst wc84` pairwise and `--fst-global`): a site where a
  population has no called alleles now contributes only the finite
  Weir-Cockerham `c` component (h̄/2 over the typed populations) instead of
  being skipped — pixy / scikit-allel `nansum` semantics, for both the
  pairwise and the r-population forms (doc_v7 D7-20). vcftools and plink2
  drop such sites (previously reed's behavior); hierfstat recomputes them
  with the typed populations only. Parity: pairwise == pixy to 13
  significant digits on a masked Chr01 with 1–3% all-missing sites
  (0.0086203 / 0.0020595 / −0.0005171), global == scikit-allel r=3
  (0.003237); the vcftools gold gate still passes (the benchmark data has
  no such sites, and values there are bit-identical).

- The engine's abort signal is an atomic flag (doc_v6 A6-1, Y5): `is_aborted`
  no longer flips back when a `select!` waiter consumes a channel token, the
  wake-up channel is sized by `--threads` (which may exceed the core count),
  and registered byte-permit pools are notified on trigger, so
  `BytePermits::acquire` waits on a condvar instead of polling every 10 ms
  (`docs/spec/architecture.md` §7.5). Every engine entry point now builds its
  abort channel from the run's thread budget (`abort_for`).
- In-flight chunk permits live in a drop guard (A6-1, Y5): a worker failure
  or panic, an aborted send, or a dropped pending message all return the
  permit and the in-flight accounting, so the reader can no longer block
  forever on a leaked permit. `StageMetrics::in_flight_final_bytes` makes the
  residue assertable.
- `threads <= 2` now takes a real sequential path (A6-3, Y5): one thread
  reads, inflates, computes and writes, with the same seam stitcher and sink
  as the parallel path; output is byte-identical (locked by a test). The
  `is_sequential` selector documented in §7.3 finally has a caller.
- The seam stitcher is bounded (A6-2, Y5): a record reassembled across chunk
  seams may not exceed `reed_engine::pipeline::MAX_RECORD_BYTES` (1 GiB,
  shared with the BCF frame chain), so a newline-free input fails the run
  instead of growing the reducer's fragment to the file size.
- `StageMetrics` gains `reader_wait`, `reducer_wait`, `inflate_time`,
  `job_time`, `deflate_time` and `in_flight_final_bytes` (A6-4, Y5). Output
  deflate time rides `ChunkOutput` to the sink (`OrderedSink::deflate_time`),
  so the three worker stages are measured instead of inferred; `--timing`
  and `--verbosity 2` print `reader_wait_ms`, `reduce_wait_ms`, `inflate_ms`,
  `job_ms` and `deflate_ms`.

- One pass0 serves every filter-side global analysis (doc_v6 03 §7, Y1e):
  the `MatrixJob`/`MatrixSink` pair fills per-contig site-major plane
  matrices that `--mind`, `--indep-pairwise` and the sketch path's mid
  stage all read back. The three previous normalizations
  (`MissingJob`, `SketchJob`, `LdPruneJob`) and `reed_filter::sketch` are
  deleted. `--mind` counts come from the miss-plane column popcounts
  (`called(s) = 2·sites − pop(slot0_miss ∨ mask) − pop(slot1_miss ∨ mask)`);
  `--mind` + `--indep-pairwise` now combine in one run (they previously
  could not share the first open).
- The pruning runs after every site rule: without `--mind` the pass0
  evaluates the tail rules inline and the pruner sees only post-tail
  sites (plink2's load-then-prune order, 02-cli §6 step 6b); with `--mind`
  the mid stage's survivor bitmap filters the pruner. The report's
  `ld_prune` counts are post-tail (`n_seen` shrinks when tail rules drop
  sites).
- `PipelineOptions::scan_head_failures` and `RecordView::head_passed` are
  gone (their only user was the sketch pass).
- The window-layout output writers are `reed::windows` (`WindowLayoutOptions`,
  `window_kinds`) instead of the previous vendor-named module (doc_v6 05 §7.2,
  E6-13); the output bytes are unchanged.
- The unsound AVX-512 counting kernel and `CountTable` were deleted
  (doc_v6 03 §5, A6-22); the mmap is now the workspace's only audited
  `unsafe` site.
- `ContigId` is a single newtype shared by every layer (doc_v6 A6-23).
- The whole tree is `cargo fmt`-clean under the pinned 1.99 toolchain, and the
  CI clippy command (`--all-targets -- -D warnings`) passes (doc_v6 05 §1–§2).

- I/O hot paths reuse their buffers (doc_v6 03 §3, A6-5/6/8, Y2): the
  deflater keeps a resident output scratch instead of zero-filling the
  worst-case bound per block, `inflate_chunk` walks the chunk's block table
  into exact-size windows (no second header parse, no per-block 64 KiB
  memset), and the chunk tail scan uses `memrchr`. The filter gold run
  drops from 460 to 400 CPU-s.
- The streaming BGZF source keeps empty members in the chunk (A6-7, Y2):
  dropping them used to shift every later chunk's `compressed_start`, so
  virtual offsets computed from stdin/pipe input could be wrong. The
  streaming and mmap sources now produce identical chunk metadata, locked
  by a test with an empty-member fixture.
- Per-record allocations left the hot paths (A6-11/12, Y2): the error
  `Locus` is built from borrowed bytes only when an error carries it
  (`LocusSource`), INFO edits write into the output buffer directly, and
  the index writer no longer builds a contig-name `String` per record on
  the (default) CSI path. The BCF GT decode table is one `Arc` per run
  instead of one per worker (A6-10).
- The filter's survivor emit is one method (`emit_survivor`) for the text
  and native BCF paths, the stale-tag list builds once per chunk, and the
  mask rules run columnar — one contiguous pass per rule over the scanned
  FORMAT columns, integer thresholds compared in the integer domain
  (A6-15/17/20, Y3). The `--set-gts .` expression path reuses the
  pipeline's scanned FORMAT table in place (no per-sample allocation) and
  expressions compile to operator enums at parse time (A6-18). INFO reads
  go through a per-record fragment index shared by the typed head rules
  and the expressions (`reed_io::vcf::info::InfoIndex`, A6-18). The
  exact-test HWE (`--fill-tags HWE`, `--hwe`, `--hardy`) is memoized per
  worker through `HweCache` (A6-16); the emitted p-values are unchanged.
  Full-gold CPU moved 398.2 → 395.6 s; the measurements behind the
  roadmap's larger estimate are recorded in `docs/rounds/v6/STATUS.md`
  (Y3 deviation).

### Added

- Streaming per-contig pruning (Y1f): each contig's matrix feeds its own
  prune thread as soon as the stream closes it, so one matrix is resident
  at a time; a single contig over the `--memory` budget fails with
  `MemBudgetExceeded`.

- `LICENSE-MIT` and `LICENSE-APACHE` (the manifests declared
  `MIT OR Apache-2.0` with no license text until now).
- `reed --version` now prints the source state the binary was built from
  (`reed 0.3.0 (v0.3.0-12-g8709f1f 2026-10-10, rustc 1.99.0)`), injected by
  `build.rs` (E6-11).
- README sections: Installation, Citation, Input requirements, Exit codes
  (E6-12).

## [0.3.0] - 2026-10-10

The current development line (the W0–W6 / X0–X7 rounds in `docs/rounds/`).

### Added

- `reed stats`: π, dxy, Fst (`wc84`/`hudson`, pairwise and global), Watterson
  θ, Tajima's D, PBS, projected SFS, f-statistics (`--f2/--f3/--f4/--dstat/--f4-ratio`),
  LD (`--ld-decay`, `--ld-pairs`), per-site/per-sample reports
  (`--freq/--hardy/--het/--missing/--depth`) and the bcftools-stats
  compatible `--summary`, all in one pass.
- `reed filter`: bcftools-view-style hard filtering, genotype masking,
  INFO recomputation (`--fill-tags`), `-q/-Q` with the `:type` suffixes,
  LD pruning (`--indep-pairwise`) as a site rule, index writing (`-W`/`--write-index`).
- `reed index`: CSI (default) and TBI, parallel.
- Window layout output (`PREFIX_pi.txt`, `PREFIX_dxy.txt`, …) verified
  column-for-column against pixy 2.2.4; `docs/rounds/v5/05-algorithms.md`.
- BCF read and write (native, no text round-trip for the filter decision),
  BGZF chunk streaming with mmap and streaming sources.
- The deterministic chunk engine: ticket window, ordered reducer, byte-budget
  backpressure; output byte-identical across thread counts and strategies.

### Verified against

bcftools, plink2, vcftools, pixy, scikit-allel, easySFS, admixtools-style
f-statistics; see `gold/` and the `Verification` section of the README.

## [0.2.0]

The pre-migration architecture (`docs/rounds/v1/design.md`), frozen as the
differential-testing oracle (`crates/oracle`). Not a published release; the
history is recorded in `docs/rounds/v2/STATUS.md`.
