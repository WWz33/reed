# 01 复审：全局 review 与源头复核

> 复审对象：HEAD `6fd061d`（v24 收尾）。
> 方法：9 个子代理分 3 批并行审查（范围见 README），汇总后对要害主张逐条对照源码
> 复核。复核分三类：**F**（复核确认，本轮条目输入）、**R**（子代理主张，实施对应
> 条目时先复现/确认）、**已修复**（v24 落地窗口内读到的旧态，勿再立项）。
> 标【据】的行号已在 `6fd061d` 复核；标【述】的行号以实施时 grep 为准。

---

## 0. 总评

记录管线的主链路纪律良好：错误枚举、exit code 映射、BGZF 写出、`parse_size`、
区域解析（v24 后）均单点；不变量注释密度高；差分/parity/proptest 测试齐全；
`expect("invariant: …")` 纪律覆盖了全部非测试 `unwrap`。review 的问题集中在四类：

1. **语义多实现点**：GT 子字段位置 6 份、计数语义 5 份、缺失值比较 3 份、AC/AF
   读取 2 份——V22-R2-1（u8 下标）当年一个 commit 修三处即为实证【述】；
2. **跨层契约用裸原语承载**：`&[bool]` keep-mask（S68 的修复是 mid_stage 加
   `assert!`，现象级）、`when: String` 阈值、可预测 spill 路径；
3. **不变式靠调用点手写 assert**：merge 形状一致 15 处 zip，部分仅 `debug_assert`
   （release 无操作）——V23-28 实证过的失效模式在其余位点仍开放；
4. **facade 边界击穿未回填**：reed-cli 直连三个 crate，其错误被字符串化进
   `Error::usage`，退出码契约按子命令失真。

这四类有一个共同的土壤：**职责混杂**。18 组判定（20 个 ≥670 行文件）中
9 组判"乱"、8/9 乱组与补丁热点重合；文件/函数把多个阶段装在一起，语义才有多
入口，补丁才只能打在调用点。因果链与拆分顺序见
[08-large-files.md](08-large-files.md)。

---

## 1. 逐领域评价

### 1.1 reed-io（12.4k 行）

- ✅ 分层清晰（line/site/scan 文本路径、bcf/* 原生路径、chunk/write 并行框架）；
  差分测试（bcf_vs_noodles、genotype_scan_fuzz）齐全；V24-06 双代缓存带边界与
  旧代 id 显式 panic 测试，是"有界结构替代阈值悬崖"的范本。
- ❌ F1：GT allele 域无校验【据】（§2.1）。
- ❌ F2（convert 侧）：GT 位置 6 份实现（convert + PCA `run_chunk`/`run_seam` 都取首字段）【据】（§2.1）。
- ⚠ POS=0 三种命运【述】：`fixed.rs:14` `parse_pos1` 接受 0（测试 pin）；
  convert 记 malformed；filter/index 写 start1=0 索引键。
- ⚠ NaN 既是缺失哨兵又是合法值【述】：`genotype_scan.rs:520-528`（V22-R2-7
  deferred），FORMAT 浮点 `nan` 与缺失不可区分。
- ⚠ 样本列契约只有 scanner 校验【述】：`encode_indiv` 短行静默填 ABSENT
  （`encode_native.rs:453-470`）；convert `parse_gts` 用 `begin.min(len)` 同样静默。
- ⚠ header 词法两套【述】：`bcf/dict.rs:38-60` 手写 vs noodles `HeaderFacts`
  （V24-03 统一了 `##TAG`，dict 的引号/转义仍靠肉眼同步）。

### 1.2 reed-core + reed-graph

- ✅ `GtCode` 单字节 + 显式解码 policy；PCG32 有 determinism 测试；UPGMA/NJ/Newick
  测试密度高、有 ape 对拍。
- ❌ F7：`restrict` 静默保留、`count_site` 生产死代码、幽灵计数器【据】（§2.1）。
- ⚠ `VariantClass` 词表三份镜像【述】：core `site.rs:58`（`*` → Other）、filter
  `head.rs:231-240`（`class_matches` 不含 Other）、cli `translate/filter.rs:99`
  （用户输入 `"other"` 别名到 symbolic）——`-V other` 永远选不中 core 语义的 Other。
- ⚠ `matches!` 元组模式 + `_` 兜底【述】：`expr/mod.rs:115`——core 加枚举变体时
  编译器不再提醒。
- ⚠ `n_missing_alleles(n_group)` 裸减法【述】：`count/mod.rs:37-39`，u32 下溢
  release 回绕；调用方（`sfs.rs:55`）负责传对分母。
- ⚠ bootstrap `built == 0` 时全 0 支持而非报错【述】（`bootstrap.rs:169-175`）。

### 1.3 reed-cli（9.1k 行）

- ✅ typed request、round-trip 测试、`ErrorKind`→exit code 单点（`main.rs:21-25`）、
  `parse_size` 经 value_parser 复用。
- ❌ F6：`--threads 0` 语义随命令翻转【据】（§2.1）。
- ⚠ `-O` 解析三份【述】（R5）；budget 构造 7 份 + `available_parallelism().expect`
  panic 风险【述】（R6）；`--ld-decay` 双解析 + 死分支【述】（R7）。
- ⚠ stats 的 usage/runtime 退出码随输入文件是否存在而翻转【述】（R11）。
- ⚠ popg 裸路由绕过自家契约【述】（R8）。

### 1.4 reed-filter（6.4k 行）

- ✅ spec→plan→eval 分层；record_pipeline 收敛 text/native 控制流且有
  `format_columns_vs_bcf` parity 钉住。
- ❌ F4：head 表达式 INFO 路由缺口【据】（§2.1）。
- ⚠ 缺失值比较三处实现【述】（R3，三处存在已【据】确认）。
- ⚠ AC/AF 读取 head/tail 两份、畸形防护不对称【述】（并入 V25-04）。
- ⚠ mask 归因顺序 single-pass 与 sketch 不一致【述】（R10）。
- ⚠ LD pruning 同 POS 账目 + r² 判定式与统计路径不同【述】（R9）。
- ⚠ `--max-mac 0` 哨兵静默当 off【述】（`plan/mod.rs:801-803`）。
- ⚠ `~`/`!~` 对多值 INFO 只读首元素；`Type=String` 可解析时按数值比较【述】。

### 1.5 reed-stats + reed-distance

- ✅ 累加器满足 merge 律（proptest）、Neumaier 求和、oracle parity、线程不变性
  测试；v24 后 WC84 单公式。
- ❌ F5：merge 形状不变式无类型承载【据】（§2.1）。
- ⚠ 公式孪生（Hudson ×2、Tajima 常数 ×2、jackknife SE ×4、π/dxy 骨架）【述】
  → 待拍板（07 §2）。
- ⚠ distance：dosage 解码 ×6、GRM ×3（`grm.rs` 注释自称 single implementation
  名存实亡）、`MissingPolicy` ×2【述】→ 待拍板。
- ⚠ `pi_hat` 全纯合输入 NaN 矩阵无告警【述】（`matrices.rs:519-521`）。

### 1.6 reed-reduce + reed-engine

- ✅ 引擎骨架质量高且有集成测试：票据窗口、`PermitGuard` drop-only 释放、
  `PublishGuard` 毒化、`ReorderBuffer` 窗口断言、FirstError 单出口。
- ❌ F10：Lanczos 残差判据 + 观测缺口【据】（§2.1）。
- ⚠ stream 符号约定疑似绕过 `SymEigResult::top_k` 归一【述】（R1）。
- ⚠ `lda::fit` 的 `fac_w = 1/(n-ng)` 在 n==ng 时除零【述】（`lda.rs:202`）：
  `dapc --xval` + 每组 2 样本 → 错误诊断 Degenerate、整个命令失败。
- ⚠ `reed-reduce/tests/` 为空【述】（钉子户，07 §2）。

### 1.7 crates/reed 门面（31k 行含测试）

- ❌ F3：空行跨 seam 硬失败【据】；F9：convert spill 未加固【据】（§2.1）。
- ❌ F2（PCA 侧，同 F2 类）：`PcaJob` 两份自写 "GT = 首个冒号字段"
  （`genotype.rs:171-174` 的 `run_chunk`、`:247-248` 的 `run_seam`）；
  `PcaJob::dosage`（`:106`）对越界 allele 返回 `None` → NaN，而不是
  `MalformedRecord`【据】。
- ⚠ `run_filter` ~680 行单体、两条 pass1 分支复制【述】（R12 前置确认）。
- ⚠ `run_convert_opened` 310 行装配+校验+渲染混合，`--min-samples-locus` 的
  gate 概念六处落地【述】。
- ⚠ rename 互斥校验"靠注释维持"的双份镜像【述】（doc_v21 M8 的修法是再抄一份
  提前，`rename/mod.rs:774-832` vs `603-639`）。
- ⚠ keep-mask 五处长度检查【述】→ V25-07b。

### 1.8 横切一致性

- ❌ F8：facade 击穿 + 退出码失真 + timing 键名漂移【据】（§2.1）。
- ⚠ 表格类输出无原子性（全仓 `fs::rename` 为 0）：长跑命令中断留半文件【述】（R13）。
- ⚠ provenance 三套 + argv 捕获 5 份无引号【述】→ V25-08c。
- ⚠ 库 crate 残留 `STREAM_DEBUG` env 后门【据】→ V25-10。
- ✅ 正面对照（值得推广，见 §4）：`parse_size`、exit code 映射、`MEM_FACTOR`、
  oracle 冻结、determinism 测试。

---

## 2. 缺陷清单

### 2.1 F：复核确认（v25 条目输入）

| # | 级 | 位置 | 缺陷 | 触发条件 / 后果 | 修法 |
|---|---|---|---|---|---|
| F1 | **P0** | `reed-io/src/vcf/genotype_scan.rs`（`parse_allele`，全文 0 处 `n_alt`）；`bcf/encode_native.rs:1342-1355`（`push_gt_code` 只界 i32） | GT allele 索引不对照 `n_alt` 校验【据】 | `ALT=G` 时 `GT 0/9`：filter 照常处理、`-Ob` 把 allele 9 编为 `(9+1)<<1`，htslib 读回报越界；逐 allele 管道混入幽灵计数 | scanner 对照 `site.n_alt()` 报 `MalformedRecord`；encoder 同钳制；域检查与 GT 解码同居一处（02 §1） |
| F2 | **P0** | `reed/src/convert/job.rs:954-955`（`parse_gts`）；`reed/src/genotype.rs:171-174`/`:247-248`（PCA `run_chunk`/`run_seam`）、`:106`（`dosage` 域检查） | GT 位置 **6 份**实现，convert 与 PCA 两份都取首个冒号子字段【据】；PCA 还把越界 allele 静默变 NaN | `FORMAT=DP:GT` 的合法 VCF：convert 把 DP 当 GT（`--012` 输出 -1，phylo/msa 计 malformed，presence 照记）；`reed pca` 全位点 dosage=NaN 后被 `any_alt` 静默跳过；`GT 0/9` 在 PCA 中当 missing 而非报错；V22-R2-1 当年一个 commit 修三处 | `Layout`/GT 解码提升 reed-io 公开 API，**6→1**（02 §2）；PCA 两条路径 + `dosage` 域检查同批收口（02 §1–§2） |
| F3 | **P0** | `convert/job.rs:1353-1354`（`line_is_sane` 对空行 true）+ `:759-766`（`run_seam`）+ `stats/run.rs:621-622`（`LineRef::parse` 硬报错） | 空行判"sane"喂统计 job【据】 | 空行横跨 chunk 边界（seam 为 `b"\n"`）或 CRLF 空行（`b"\r"`）→ `convert --sfs` 整体失败；是否触发随 `--memory/--threads` 变化，间歇复现 | `RecordClass { Record, Malformed, Blank }` 单点分类，五入口只选策略（02 §3） |
| F4 | **P0** | `reed-filter/src/eval/head.rs:135`（`from_line(..., true)`）；`expr/parser.rs:251`（`_ => false` 吞 `Expr::Info`）；`plan/mod.rs:562/577/881`（按该谓词路由） | head 表达式硬编码信任 INFO【据】 | 100 样本 `INFO/AC=80`，`-s s1 -e 'INFO/AC>1'`：子集后 AC=1 应保留，head 用全队列 80 丢弃；裸 `AC>1` 路由 tail 反而正确——同一量两种拼写两种队列 | 路由统一 `expr_reads_gt_fields`；`info_trusted` 一律取 plan 单字段（03 §1） |
| F5 | P1 | reed-stats 15 处 `iter_mut().zip`【据】：`het.rs:85/88/91`、`ld.rs:206/209`、`matrix.rs:241`、`fst.rs:488`、`sfs.rs:73/291/479`、`theta.rs:88/91`、`window.rs:108/137/164` | merge 形状不变式无类型承载；`sfs.rs:289/479` 仅 `debug_assert`（release 无操作），`het.rs:85`、`sfs.rs:73` 无保护 | 调用方 bug（不同 `n_samples`/LD 边界的 partial）→ release 下尾部静默丢弃、输出错数——V23-28 已实证此失效模式 | `merge_zip_eq`（release 也检查），15 处收敛；theta resize-then-zip 收编（04 §1） |
| F6 | P1 | `reed-cli/src/translate/distance.rs:102-105`；`reed-engine/src/budget.rs:62` | `--threads 0` 语义随命令翻转【据】 | distance：`num_threads(0)` = rayon 默认 = 全核；filter/stats/popg/pca：`threads <= 2` 走顺序单线程。共享集群上限线程的用户静默占满节点 | `params` 模块 `Threads` newtype 拒绝 0（05 §1–2） |
| F7 | P1 | `reed-core/src/group.rs:177-181`；`count/site_counts.rs`；`bcf/gt.rs:66/199/217` | `restrict` 长度不匹配静默保留旧分组；`count_site` 生产死代码（唯一调用者在冻结 oracle）；`n_haploid_as_missing` 位打包/解包自成环、无生产者，报告列恒零【均据】 | 9 样本 keep mask restrict 10 样本表 → 分组统计在未过滤 cohort 上算完，零报错（S68 同类）；用户按报告排查 haploid 去向被误导 | `restrict` 返回 `LengthMismatch`；`SiteCounts::from_codes` 唯一入口；计数器删或补生产者（04 §2–4） |
| F8 | P1 | `reed-cli/Cargo.toml:16-18`；`translate/distance.rs` `Error::usage`（:213/258/284【述】）；timing 渲染 4 套 | facade 边界击穿：graph/reduce 错误字符串化进 usage（数据故障 exit 2，违反 `error.rs:31-33` 自家契约）；`--timing` 键名已漂移（convert `chunks` vs 其他 `chunks_read`，busy/chain_wait 缺失）【据】 | NaN 距离 `--tree nj` → exit 2，脚本按"改命令重试"死循环；性能回归脚本要按子命令写 4 个解析器 | `reed` re-export 三 crate；`Error` 增 Runtime 载荷；`StageMetrics::render_tsv` 单点（06 §1–2） |
| F9 | P1 | `convert/job.rs:380-384`（`spill_path()` 可预测名，该文件 `create_new` 0 处） | convert spill 未加固【据】；input.rs 的 O_EXCL 加固（doc_v20 SEC-1，`input.rs:202-276` + 测试 `:416-476`）只盖 stdin | 共享主机预置 symlink `/tmp/reed-convert-{pid}-{n}.spill` → 截断/写受害者文件；spill 内容流回输出构成投毒 | `TempSpill` 共享基础设施（06 §4） |
| F10 | P1 | `reed-reduce/src/lanczos.rs:461`（残差 `.max(1e-300)`）、`:464`（`STREAM_DEBUG`）、`:521`（`max_residual` 返回）；`reed/src/pca/mod.rs`（0 处引用 `max_residual`） | Lanczos 观测缺口【据】 | 秩亏输入（重复个体、LD 剪枝后位点数 < `--pcs`）：相对残差不可达 → 固定烧满 200 轮（每轮 O(n·m²) 再正交化）；非收敛无声；env 后门绕过 verbosity | 混合尺度 `max(|λ|, λ_max·eps)`；非收敛告警；删后门（06 §5） |

### 2.2 R：子代理主张，实施对应条目时先确认

| # | 主张 | 位置【述】 | 处置 |
|---|---|---|---|
| R1 | stream 求解器绕过 `SymEigResult::top_k` 的符号归一（`backend.rs:55-64` 确有归一【据】） | `lanczos.rs:521`；`pca/mod.rs` | 独立小项：实跑 dense vs stream 对比 `.eigenvec`；分叉则修 |
| R2 | mask `<=`/`>=` 经 `to_args()` 回放成 `<`/`>`（V23-24 修一半） | `spec/render.rs:126-143` | **并入 V25-04**（03 §4）：round-trip 复现后修；`when` 携带 `Cmp` 而非靠字符串前缀 |
| R3 | 缺失值比较三处行为分叉（三处存在已【据】：`expr/mod.rs:525`、`head.rs:301`、`record_pipeline.rs:524-560`） | 同左 | 并入 V25-04：先 failing 复现 |
| R4 | 计数语义 5 份镜像行号 | `record_pipeline.rs:600-630`、`mid_stage.rs:92-135`、`reports.rs:568-581`、`count/mod.rs:56`、`site_counts.rs:124` | 并入 V25-07a：grep 刷新后迁移 |
| R5 | `-O` 解析三份 | `translate/filter.rs:393-431`、`rename.rs:46-91`、`convert.rs:55-79` | 并入 V25-06 |
| R6 | budget 构造 7 份 + `available_parallelism().expect` ×6 | `translate/{common,pca,distance,convert,rename,filter}.rs` | 并入 V25-06 |
| R7 | `--ld-decay` 双解析 + 死分支 + 0 值静默空输出 | `cli/mod.rs:686-688`、`popg.rs:615-627` | 并入 V25-06 |
| R8 | popg 裸路由绕过 `-o`/`--pop` 契约；`-o ""` 退出码与 stats 不同 | `popg.rs:25-37/235-286` | 待拍板（07 §2） |
| R9 | LD pruning：同 POS 双记录账目矛盾；pruner 乘法判定与统计路径除法 r² 边界不同 | `pass0/prune.rs:33-57`；`ld.rs:528-551` vs `:98-112` | 待拍板（07 §2） |
| R10 | mask 归因顺序 single-pass（表达式先）与 sketch（声明序）不一致 | `record_pipeline.rs:483-516` vs `mid_stage.rs:141-155` | **本轮拍板**（03 §5）：A 并入 V25-04（荐，限单文件改动）/ B 保持现状并加测试 |
| R11 | stats 的 usage/runtime 退出码随输入文件是否存在而翻转 | `translate/stats.rs:206-225` | **并入 V25-08**（06 §1）：全部 usage 检查移到 `prelude()` 前，补退出码 e2e |
| R12 | `linalg.rs` 198 行死代码（V24-03 后是否已删需 grep） | `reed-reduce/src/linalg.rs` | **并入 V25-10**（04 §5）：文件仍在、`lib.rs` 无 `mod linalg`、全仓引用 `reed_core::linalg`；删除 |
| R13 | 表格类输出无原子性（全仓 `fs::rename` 0 处） | pca/dapc/sfs/ld/distance 写出点 | 待拍板（07 §2） |
| R14 | piHat 表格 `rsplit_once('_')` 伪造 FID/IID | `translate/distance.rs:394-399` | **并入 V25-08**（06 §1.4）：两列都写样本名，列为验收项 |
| R15 | `##reed_command` 用 `env::args().join(" ")` 无引号拼接 ×5 | `translate/{filter,rename,popg,convert}.rs` | 并入 V25-08c |

### 2.3 v24 已修复（review 读到旧态，勿再立项）

| review 主张 | 实际 |
|---|---|
| FORMAT 缓存"满即清空"（CACHE_LIMIT=64） | `a0b24df` V24-06：双代缓存各 256，"never clears everything"，带边界测试 |
| WC84 pair/multi 双实现 | `248eabc` V24-05：唯一 `wc84_site` |
| 区域解析多入口 | `74090da` V24-04：header-aware 唯一解析 |
| genotype.rs 每记录 `with_capacity` 分配 | `76996f3` V24-07：现 0 处 |
| 逐值 writer / LLE 无契约 | `1047c43` V24-08：`BufWriter` 统一 + `try_lle` 上限 5000 |

### 2.4 复核方法

- 复核在 `6fd061d` 上逐条 grep/read：F1–F10 全部成立（其中 F7 的幽灵计数器，
  子代理表述"无任何 `+= 1`"不准确——`bcf/gt.rs:218` 有递增，但 `TRIG_HAPLOID_MISSING`
  位全仓仅定义/打包/解包三处、自成环无生产者，"恒零"结论成立）；
- v24 落地窗口（01:09–02:55）与子代理阅读窗口重叠，故【述】行号可能有偏移；
  实施纪律：每个条目动手前先 failing 测试/grep 重新锚定。

---

## 3. 补丁史证据（考古摘要）

236 commit：fix 63（27%）/ docs 53 / feat 46 / perf 22 / refactor 16 / test 9 / other 7。
**正确性补丁占大头且最反复；每个反复案例的根源都是结构缺口**——历轮正确性补丁
一直在为缺结构付利息。反复模式：

| 模式 | 补丁轨迹 | 根源 | 状态 |
|---|---|---|---|
| NaN 距离排序 `total_cmp` | v20 S31（3 处）→ v23 V23-08（3 处 + 2 同风格） | 无共享比较子；修复纪律是"按评审点名修"而非"按模式扫仓" | 待拍板 |
| 哨兵缓存 `usize::MAX` | v20 R2（FormatLayoutCache panic）→ v21 M6（ContigTable 790 GB abort） | last-hit + 哨兵惯用法逐结构手抄 | 待拍板 |
| `-o`/`-O`/前缀校验 | v20 R4 → v21 S15/S19 → v23 ⑨4/⑨5 → V23-21/31 收敛 | 无参数校验层 | V25-06 |
| 区域语义 | v21 M4/M5 → v23 V23-60/61/26 → V24-04 | 无 header 上下文解析 | **已修** |
| sites-only 8 列 | v22 V22-1…4/11（12 条目）→ R22-3 `LineRef::shape()` 收口 | 各命令自己数 tab | **已修（正面案例）** |
| 内存预算 | v20 S76 → v21 M11 → v23 V23-10 | 无集中资源核算 | 待拍板 |
| per-record 分配 | Y2/Y3 → v23 D 组 → V24-07 | scratch 所有权 | **已修** |
| 双实现漂移 | 组文件双解析器（#56 → M10 删）；样本源双份（S32 → V23-40）；GRM 双份 | 复制代替委托 | 部分 |
| revert 循环 | V23-34（64→1024，5 分钟内 revert）→ V24-06 根治 | 阈值悬崖（调参数不换结构必然回弹） | **已修** |

**钉子户规律**（滚 ≥3 轮未修）：正确性 P0 在 1–2 轮内必被修；挂住的全部是 CI 门、
文档诚实度、纯开销、测试金字塔——这些条目没有"拍板"机制（07 §2 列全）。

---

## 4. 补丁热点与下沉地图

> §3 按补丁轨迹看历史；本节换成两个工程视角：**热点**（同一语义被修/被抄的
> 次数排行）与**下沉地图**（每类语义该去哪个层）。二者合起来回答：下一个补丁
> 最可能落在哪，以及修哪里能让它不再落。

### 4.1 热点排行（按反复证据）

| # | 热点 | 反复证据 | 处置 |
|---|---|---|---|
| 1 | CLI 参数校验（`-o`/`-O`/`-W`/`--threads`/`--memory`） | 4 轮 15+ 补丁（v20 R4→v21 S15/S19→v23 ⑨4⑨5→V23-21/31）；`-O` 解析 3 份、budget 构造 7 份 | V25-06 |
| 2 | GT 语义（位置/域/解码） | V22-R2-1 一 commit 修三处；位置 **6 份**（convert + PCA `run_chunk`/`run_seam` 都取首字段，F2）；PCA `dosage` 是绕过 policy 的第 3 个域检查 | V25-01/02 |
| 3 | 缺失值比较真值表 | V23-59 只修 2/3；现存 3 份且已分叉（R3） | V25-04 |
| 4 | 计数语义 | 5 份镜像，core 本体 `count_site` 是生产死代码（F7） | V25-07a |
| 5 | merge 形状不变式 | V23-28 修 3 处后仍开放 15 处（F5） | V25-05 |
| 6 | 公式孪生（stats） | D7-21 打两次（dxy.rs 注释自认 "drawing the same fix"）；Hudson×2、Tajima×2、SE×4、π/dxy 骨架 | v26（site_terms） |
| 7 | NaN 排序 / 哨兵缓存 | `total_cmp` 两轮 6+ 处；哨兵 bug 两个手抄缓存各咬一次（R2 panic / M6 790GB） | v26（kNN / LastHit） |
| 8 | 行合法性 | 5 入口各自推导；B13-2 只对齐一对，空行分叉仍在（F3） | V25-03 |
| 9 | 内存预算 | 6 处检查、3 种公式、3 种超限行为（v20 S76→v21 M11→v23 V23-10） | v26（`budget::admit`） |
| 10 | facade 直连区 | `grm.rs` "single implementation" 名存实亡（×3）；isomap×2（V23-08 打两遍）；退出码/timing/provenance 漂移（F8） | V25-08 |

### 4.2 下沉地图（按目标层）

**→ reed-core（纯类型/语义）**

| 下沉物 | 从哪来 | 消掉 |
|---|---|---|
| `SiteCounts::from_codes(&[GtCode], &SampleMask)` | 5 份计数镜像 | 热点 4 |
| `KeepMask`（构造即校验长度/名称） | 5 处长度检查 + mid_stage assert | S68 类 |
| `Provenance`（自 reed-filter 提升） | 3 套实现 + argv 捕获 ×5 无引号 | 热点 10 |
| `PloidyPolicy` 封闭枚举 + dosage 视图为 GtCode 方法 | PCA 的第 3 个 GT 解析器 | 热点 2 |
| VariantClass `from_token`/`as_str` 双射 | 3 份词表镜像 + `_` 兜底 | 词表漂移 |
| （v26）`LastHit<T>` interner、`knn_indices`+比较子 | 两个哨兵缓存、6 处 total_cmp | 热点 7 |

**→ reed-io（记录解码视图）**

| 下沉物 | 从哪来 | 消掉 |
|---|---|---|
| `Layout` + `decode_gt`（带 n_alt 域检查）公开 API | 6 份 GT 位置（含 PCA `run_chunk`/`run_seam`）+ 域无校验 | 热点 2 |
| `RecordClass::classify_line` | 5 入口各自推导 | 热点 8 |
| `validate_sample_cols` 公开 | encoder/convert 静默填 ABSENT | 样本列契约 |
| （v26）数值域+哨兵类型 | text/BCF 双实现（`parse_int_sum` vs `int_sum_column`、NaN 双语义） | 数值域双修类 |

**→ reed-filter `semantics` 模块**：`compare_missing`（3 份真值表）、`info_af`
（head/tail 两份 AC/AF）、`expr_reads_gt_fields` 唯一化——热点 3 + F4。

**→ reed-stats**：`merge_zip_eq`（15 处）——热点 5；（v26）`site_terms` 层
（`hudson_terms`/`tajima_constants`/`se_from_pseudovalues`/π-dxy 泛型）——热点 6。

**→ reed-engine**：`ThreadBudget::from_cli`（7–8 份构造）、
`StageMetrics::render_tsv`（4 套漂移渲染）；（v26）`budget::admit` 唯一执行点——热点 9。

**→ reed-cli `params` 模块**（唯一"新建而非下沉"的层）：`Threads`/`MemMiB`/
`OutputPrefix`/`OutputType` newtype + value_parser——热点 1 的另一半。

**→ reed facade**：re-export 三 crate + distance 请求构造下沉 `run_*_opened`
（照 pca/dapc 模式）——热点 10 的根源。

### 4.3 不下沉清单（刻意保留的重复，防过度收敛）

- **oracle 的独立 `count_site`**：冻结的 parity 参照，必须与生产实现独立；
- **text 与 BCF 双 reader**：性能设计——下沉的是**语义**（值域/哨兵/GT 语法），
  不是合并 reader；
- **reed-stats 与 reed-distance 的 per-pair 累积**：输入类型与分母约定不同
  （等位级 vs 基因型级），sealed 设计的有意分叉，登记在 statistics-map 即可；
- **word kernel（mid_stage）**：性能豁免，但必须补与 `from_codes` 的等价性
  测试（当前没有，见 04 §2）。

### 4.4 对账

热点 = "每个消费者各持一份"；下沉 = "每类语义一个唯一入口"。v25 覆盖热点
**1/2/3/4/5/8/10**；**6/7/9** 留 v26（规模大、无正确性缺口，先让正确性落地）。

### 4.5 补丁下沉审计

§3/§4 的清单已整理为 [00-patch-sinking.md](00-patch-sinking.md) 的 P1–P6 判据
与五类动作（Sink/Delete/Keep/Verify/Split）；每个条目实施前先填 "retired patch"
与 "sink target" 两栏，完成后按 00 §6 的四问验收。补丁族与下沉物的对照见 00 §4。

### 4.6 结构因果（大文件审查）

18 组判定（20 个 ≥670 行文件；9 个子代理 × 3 批）的判定与复核见
[08-large-files.md](08-large-files.md)。结论：**9 组判乱、8/9 乱组与补丁热点
重合**；
但文件混杂本身不直接产生 bug，它产生"语义多入口"，多入口让补丁只能打在调用点。
因此拆分只作为语义下沉的载体，顺序固定为：

1. `encode_native.rs` → `bcf/gt_decode.rs`（随 V25-01/02，GT 代码不再第四次搬家）；
2. `filter/mod.rs` pass1 合一 + BCF magic 下沉（随 V25-02/03）；
3. `convert/job.rs` spill → `convert/spill.rs`（随 V25-09）；
4. `reports.rs` 计数/policy 收口（随 V25-07）；
5. `plan/mod.rs compile` 纯搬家（V25-04 之后）；
6. `request.rs` / `matrix.rs` / `popg.rs` → v26。

明确不拆：`pipeline.rs`（v9 决策）、`genotype_scan.rs`、`expr/mod.rs`、
`cli/mod.rs`、`index.rs`、`record_pipeline.rs`（正面对照）。

---

## 5. 值得保留并推广的做法

- **单点纪律的正面样本**：`parse_size`（units.rs + value_parser）、exit code 映射
  （main.rs:21-25）、`MEM_FACTOR`（budget.rs:118）——校验只写一次的参数从不反复；
- **oracle 冻结 + parity 门**：历轮重写正确性的锚（V24-01 又补了 golden）；
- **determinism / 线程不变性测试**：锁住"策略间字节一致"的承诺；
- **V24-06 双代缓存**：有界结构替代阈值悬崖的范本——v25 每个条目都应以
  "修完后该类补丁失去繁殖地"为达标线；
- **v22 R22-3 `LineRef::shape()` 收口**：一次结构投入消 12 个散点条目，是
  "修底层、消补丁"杠杆率最高的历史实证；
- `expect("invariant: …")` 纪律覆盖全部非测试 `unwrap`。
