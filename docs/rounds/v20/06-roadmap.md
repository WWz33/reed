# 06 后续路线（v20 批次二已清 bug 类 16 条）

> 原则沿 v2 06：**先修正确性与会撒谎的门,再补测试金字塔,最后清 lean**。
> 来源:全库评审的 should-fix(~66 条)与 lean 清单,明细在仓库外
> `../../../reed-rs-review.md`(编号 S15–S80 对应其 Should fix 列)。
> 每条 bug 类条目沿用 v20 纪律:先写失败测试,再修。
> **状态(2026-10-24 第二批后)**:标 ~~删除线~~ 的已修,验证与修法见
> STATUS.md「第二批」表;其余待做。

---

## 1. 通用验收门（每步都过）【规】

1. `cargo test --workspace` 全绿;改动 crate `cargo clippy --all-targets`
   无新告警。
2. 数值类改动:既有 gold 对拍/oracle parity/手算锚点不变,或变动可归因于
   消除的 bug 并在测试内注明。
3. 确定性:涉并行的改动加 1 vs N 线程逐位一致测试(R13 的模式)。
4. 上游行为对齐类(UMAP/t-SNE/PPCA/plink2):引用上游源码行号写进注释,
   crosscheck example 的差值要有解释。

---

## 2. P0：会撒谎的门与静默错果

| 项 | 内容 | 回归测试 |
|---|---|---|
| CI-1 | data-gates 不构建被测二进制(跑 PATH 里的旧 `reed`),`run_all.sh` 不调 `run_reed.sh`——绿灯无证据力;先 build 再 parity | 在 CI 里故意注入错值,门必须红 |
| CI-2 | `cargo fmt --check` 与 `cargo deny check` 入 gates;先清现存漂移(含 doc_v18 WIP 文件) | CI 自身 |
| ~~SEC-1~~ | stdin spill 可预测临时路径 + `File::create`(共享 HPC 上可被符号链接劫持);改 `create_new` + 重试 | 已修(`sec1_*` ×3,修前受害者被截断为空) |
| ~~S58~~ | 表达式 AC/AN/AF/MAC/MAF 在 cohort 变化后仍读陈旧 INFO(`info_or_gt` 无视 `info_trusted`),与同 run 的 `--maf` 口径不一致 | 已修(`s58_*`;`from_line` 收参数,head/INFO 直读不受影响) |
| ~~S28~~ | UMAP 学习率下限 `.max(0.25)` 与 umap-learn(线性衰减到 0)及其自身注释矛盾 | 已修(`s28_*`;CLI 无 umap 出口,无锚点变动) |
| ~~S29~~ | t-SNE 困惑度二分 64 次不收敛时整行静默清零(n≤31 且默认 30 = 全噪声);非收敛也要写最终权重 + 拒绝 `perplexity ≥ (n−1)/3` | 已修(`s29_*` ×3;`tsne` 变 Result + 类型化 TsneError) |
| ~~S31/S32~~ | NaN 基因型 panic 三处 kNN 排序(`partial_cmp().unwrap()`,改 `total_cmp` 或入口拒绝);`xval::run` 组 id 越界在 LDA 类型化错误前 panic | 已修(`s31_*` ×3 + `s32_*`;graph.rs:326 isomap 第四点留 backlog) |
| ~~S36/S37~~ | bootstrap 根节点恒 100% 支持(`branch_keys` 跳过全 1 键);`--tree --bootstrap` 的 discordance 用默认参而非用户 flags(R13 剩余的 `--threads` 接入同批) | 已修(`s36_*` ×3 support 对齐 ape;`s37_*` ×3 e2e 含 1v4 线程逐字节一致,`--threads` 已接通) |
| ~~S76~~ | pca/dapc 全矩阵无内存预算检查(5000×3M ≈ 120 GB 直接 OOM);对照 `--memory` 报 usage 并指向 `reed filter --indep-pairwise` | 已修(`s76_*` ×2;`run_pca_opened` 顺带收敛,消掉 P2 重复块) |
| ~~S46~~ | WC84 pair 在 an=1 时整位点丢弃而 multi 保留——同数据 `--fst wc84` ≠ `--fst-global`;pair 路由走 keep-c 分支 | 已修(`s46_*` ×2;**评审有误**:multi 在 r=2 同样跳过,两路径一起对齐;oracle parity 按回加项调整并注明) |
| ~~S68~~ | keep-mask 长度约定四处不一致(sketch 尾样本被丢 / 单遍尾样本被留 / core 文档说留 / pca 报错);生产两循环改直接索引 + `debug_assert` 长度相等 | 已修(`s68_*` ×2;契约:短 keep 是 debug 期编程错误) |
| ~~S79~~ | oracle PBS 模块文档与自身代码矛盾(文档说三群交集,代码是 per-pair);冻结 oracle 的契约以文档为准会被误移植 | 已修(doc 改为 per-pair + `s79_*` 契约测试钉死;oracle 代码零改动) |

## 3. P1：验证面、io 边角、测试金字塔

| 项 | 内容 | 回归测试 |
|---|---|---|
| CLI 批 | popg `-o -` 跑完全程才拒(检查上移);`--hwe midp` 文档存在但必拒(parse 期拒绝);popg 静默 no-op(`--proj` 无 `--sfs` 等,照 convert 模式拒绝);pca/distance 无视 `--verbosity`;`--indep-pairwise` r² 越界;`--xval-rep 0`;filter/rename 的 `-O` 级别不校验(convert 已校验);convert 默认 prefix 丢目录 | 每条一个 usage-error/行为测试 |
| io 批 | ~~FORMAT 整数越界钳制而非置 missing(对齐同文件 INFO 路径)~~;~~索引 uoffset 越界 panic(`split_text` 边界)~~;多区域 `-r` 折叠成 min..max 全量扫描(区间间断 seek);~~符号等位忽略 INFO/END(rlen=1,区域查询漏跨度删除)~~;~~POS>2³¹ BCF 钳制~~;0 样本文件误导性报错;FORMAT 缓存每记录 `to_vec`(arena 化) | S39/S40/S42/S43 已修(见 STATUS 第二批;S42 三处 span 收敛为共享 `indexed_span_len`);其余各一条畸形输入/对拍测试 |
| stats 批 | `--variant-only` 空窗 θ 打 0.0 而 π/dxy 打 NA(统一 NA);`empty_group_is_rejected_at_build_time` 空断言(补真校验或改测试口径);`QcReportAcc::merge` 无定律测试(JSON 序列化对比最廉);f2/f3/f4 不在定律套件 + f3 无手造混合夹具;`pbs.rs` a 槽文档反了;jackknife 丢弃零分母块无信号 | 定律/手算夹具各一 |
| core 批 | `SampleMask::from_ids` 重复 id 双计 `n_selected`(pop 文件重复行 → 公共 API 分母错,CLI 现不受影响);`start+size` 溢出(saturating);`HweCache` 每命中克隆向量;`from_tsv`/`names_from_tsv` 双解析器(共享行迭代器 + 同文件一致性测试);`restrict` 长度不符静默保留全组 | 各一条单元测试 |
| filter 批 | `--hwe` 规则零测试(含边界/退化 p=1/midp 拒绝);`eval/head.rs` 无单元测试(QUAL==阈值、MinAfInfo 边界、AC 长度回退);全掩码位点过 `--max-maf` 但被 `-Q` 丢(口径统一或文档);`min_qual = NaN` 通过校验且全丢 | head/tail 单元测试模块 |
| input 批 | `reed index` 对 VCF.gz 每 chunk 双重解压(job 收膨胀视图);跨 contig 乱序 CSI 报 noodles 原始错/TBI 静默接受(保持 last_tid 报"needs sorted input");BCF `-r` POS 谓词 vs bcftools 重叠语义(改谓词或 README 注明);seam 路径行尾 `\n` 不一致(契约化) | 各一条;双解压以计数/断言可测 |
| 测试金字塔 | pca/dapc collection seam(`genotype.rs`)零单元测试(dosage 边角 token:半调用/单倍/相位/allele>n_alt/空 GT + keep 重映射);`reed pca`/`--dapc` 零 facade e2e(锚定几个特征值/ASSIGN 格);`--replace-sample-name`、`-N/--samples-file` 零测试;`reed-reduce/tests/` 为空(已知答案特征值 + 跨后端容差 + PCA 线程数不变) | 命名的测试本身 |
| 其余跨 crate | `count_site` 死而文档自称权威(删或三实现等价测试);GraphError 被压成 exit 2(退出码契约);`Pos0`/`Pos1` 新类型在 facade 零使用(两处 `checked_sub(1)` 换新类型);SFS-only convert 冗余逐样本解析(`has_row_targets` 早退);PCA/DAPC 请求上死的 `threads/memory` 字段(删);R2 同型:空 CHROM 哨兵 | 各一条 |

## 4. P2：lean(评审明细见外报告 Nice-to-have 节)

| 项 | 内容 |
|---|---|
| 死代码 | reed-core:`Sites` 窗格族、`AccessibilityMask`、`resize_sites`、`write_genotype`、`Pos1/Interval1` 族、`GroupTable::masks`、`SampleTable::name`、`Rng::{sample_n,multinomial}`;reed-stats:`FstMatrixAcc`(文档宣称的优化不存在)、`BlockJackknife/BlockPartial`、`IntHistogram::{quantile,weighted_sum}`、`SiSAcc.allele_count`;reed:`PcaScratch`(整型)、~~`run_pca_opened` 与 `read_genotype_matrix` 同文件重复块~~(已随 S76 收敛) |
| 逐行重复 | `isomap_from_distances` ≡ `graph::isomap`(O(n³) 两份);`-O` 三解析器已由 R4 收敛,TSV 双解析器在 P1 |
| 文档失实 | `MatrixFormat::Triangle` 对角线方向;`KingCounts.het1hom2` 标签;ld.rs TIE_EPSILON(1e-6 → 实为 2⁻⁴⁴)与"精确整数"断言(仅到 n≈3500);lib.rs "双后端结果相同"(特征向量不然);`RegionChunkSource` 模块文档;Ajk 对角线;`Rng::sample_n` 内存界;`.afreq/.hardy` 多等位坍缩未注明;`PREFIX.ld` ID 列硬编码 `.` |
| 小性能/小修 | `WcMultiAcc::push` 每位定三 Vec(构造期 scratch);~~PCA GRM 对角测试容差 0.6(换干净夹具收至 1e-9,它正掩盖过 R9)~~(已随第一批 R9 换夹具);swap 重命名误导性报错;hapmap POS 0 逐字透出;零共享位点距离 0 而非 NA(README 承诺 NA);`ContigTable` 空列哨兵(并入输入批亦可) |

## 5. 明确不做(评审裁决记录)

- t-SNE 的 Barnes–Hut/FFT 加速:n 增大前不做,先保正确(R8 已修)。
- distance 的 sites 流式分块:等 P1 的 driver `--threads` 接入后再评估,
  避免与 WIP 冲突的两次重写。
- reed-distance/reed-reduce 是否并入 workspace 锁步(lint/版本/publish):
  【决】留用户——两者读起来是刻意保留的独立库形态。
