# doc_v23 复核状态（2026-11-03）

> 来源：仓库外 `reed-rs-review-2026-11-03.md`。本文记录**复核判定**、
> **拍板结果**与**实施状态**。实施明细见
> [02-implementation.md](02-implementation.md)。
>
> 验证环境：本机 debug 构建 + `bcftools 1.24-2`/htslib 1.24-1；崩溃类用
> 最小 fixture 实跑；数值类用 `/tmp/v23probe`（path 依赖，未改仓库）；
> 行号按当前工作区（含 doc_v22 已提交修复 `8b54588`）。证据见
> [01-verification.md](01-verification.md)。

## 拍板结果（2026-11-03）

| 编号 | 拍板 | 落地要点 | 状态 |
|---|---|---|---|
| V23-09 | **按结果正确性改**（选项 A：`SymEigResult::top_k` 统一符号规范化） | 每列最大 \|分量\| 取正（并列取行序最先；全零列不动），PCA/PCoA/kernel/Isomap/Laplacian/diffusion 共用；`eig::SymEig::top_k` 同规则；落地需补 V23-48 双 backend 对拍；接受 `.eigenvec`/坐标符号变化 | **已实施** |
| V23-14 | **选项 C**：sink 按 contig 顺序直接追加到最终 `PREFIX.ld` | **首个 pair 到达时懒创建并写表头**（全无 pair 时仍不创建 `.ld`，与现状一致）；逐 contig 追加；不落临时文件，输出逐字节一致，内存 O(单 contig)；落地补逐字节对比测试 | **已实施** |
| V23-07 | **丢弃（不修）**：不接受多等位 → 该条从 v23 修复清单丢弃 | 现状保留：默认 `BiallelicGate` 已丢弃多等位；显式 `-M >2` 时 AC 仍按首 ALT（不支持路径，不实现 any-pair） | 已拍板：丢弃 |
| V23-50 | **按推荐 A**：AF/MAF 的 NaN 语义不变 | 保持与 bcftools 一致；文档写明 AN=0 时 AF/MAF 为 NaN、按 IEEE 比较，补测试钉住 | **已实施** |
| V23-21 | **按推荐 A**：前缀命令统一拒绝 `-o -` 与空前缀 | `stats`/`convert`/`pca`/`distance` 共用 `require_output_prefix`；`popg` 的检查上移到参数校验；`filter` 的 stdout 语义不变 | **已实施** |
| V23-47 | **按推荐 A**：oracle 保持 complete-case | oracle 与生产/spec 同为整站跳过；补含缺失位点的 golden；文档记录与 allel 的差异（不改为逐位点删除） | **已实施**（文档 + golden） |

> 上一行是历史记录。实施后：A 组 6 条、V23-08/09/10/11/12/13、V23-14/21/31/47/50、
> V23-59/60/61 已落地（见 [02-implementation.md](02-implementation.md)）；其余
> 条目（V23-15…20、22…30、32…46、48、49）仍为**待拍板**且未改代码。

## A. 建议 P0（崩溃/静默写错数据）

| 编号 | 来源 | 判定 | 建议修法 | 状态 |
|---|---|---|---|---|
| V23-01 | ①3-1 | **属实** | `typed()` 的 count=15 只允许一层扩展，拒绝二次扩展 | **已修**（`bcf/record.rs::extended_count`；4 项测试） |
| V23-02 | ②2 | **属实（filter 有守卫；convert/stats panic）** | `fmt_gt` 非整数类型 → malformed；`float_element_column` 非数值 → NaN；函数级测试 | **已修**（`fmt_gt` 返回 Result；encode 守卫收紧；3 项测试） |
| V23-03 | ①1 | **属实（需≥1 个额外空列）** | `edit_vcf_line` 传 `n_samples`，越界 break；`-s`/掩码两种编辑各补测试 | **已修**（样本数取自 plan 的 bitset，公开签名不变；`-s`/掩码 + 单 tab 对照） |
| V23-04 | ①2 | **属实（数据损坏，非潜伏）** | 指数分支补 `-`；INFO/FORMAT 两个入口 + 负值/边界测试 | **已修**（9 值直接测试 + INFO 编辑；真实 BCF 与 bcftools 逐字节一致） |
| V23-05 | ④3 | **属实** | parser 加 depth（如 500），超限 Err；499/500/501 测试 | **已修**（`MAX_NESTING=128` + `MAX_CHAIN=2000`；拍板的 500 经实测在 2 MiB debug 栈不安全，按 430/200 层边界取 128） |
| V23-06 | ③1 | **属实** | `frac == 0.5` 精确判平局；三组实测值 + DAPC 回归 | **已修**（精确比较 + 0.05…0.95×1…500 扫描） |

## B. 建议 P1（正确性/复现性/内存）

| 编号 | 来源 | 判定 | 建议修法 | 状态 |
|---|---|---|---|---|
| V23-07 | ④2 | **属实（与 bcftools 分歧）** | **已拍板：丢弃（不修）**——不接受多等位，该条从修复清单丢弃；默认门已丢弃，显式 `-M >2` 时 AC 按首 ALT（不支持路径，不实现 any-pair） | 已拍板：丢弃 |
| V23-08 | ⑧1 | **属实（库 API；CLI 无入口）** | `graph.rs:326`、`grm_reduce.rs:164/207` 改 `total_cmp`；NaN 测试 | **已修**（3 处 + 2 处同风格；2 项测试，另钉住 `grm_distances` 的 NaN→0 clamp） |
| V23-09 | ⑧2 | **属实（backend 间实测整体变号）** | 已拍板：按结果正确性改——`SymEigResult::top_k` 统一做最大 \|分量\| 取正（LDA 同款）；接受输出符号变化 | **已实施**（A） |
| V23-10 | ⑤1 | **属实（代码）** | sketch pass0 累计驻留字节与预算比较；多 contig 小预算测试 | **已修**（`ResidentBudget` 累计和 + 饱和；3 项单测 + 多 contig e2e） |
| V23-11 | ⑤4 | **属实（代码）** | `read_reference` 流式解析，避免整文件 + Arc 双份 | **已修**（`BufReader::split` 逐行；CRLF/空行/多记录 e2e） |
| V23-12 | ⑧6 | **属实（代码）** | `matrices::grm` 直接读 `st.w`，内积循环顺带累加 `called_both`；逐值对拍 | **已修**（单遍内积 + 计数；`reference_grm` 逐值对拍） |
| V23-13 | ⑥2 | **属实（单 contig 内；报告"全 contig"需修正）** | 窗口推进时截断 `pos_log` 前缀；大 span 测试 | **已修**（`pos_base` + 前缀 drain；2 项单测） |
| V23-14 | ⑤2 | **属实（`--ld-decay + --ld-pairs` CLI 可达）** | 已拍板：选项 C——sink 按 contig 顺序直接追加到最终 `PREFIX.ld`（首个 pair 懒创建并写表头，全无 pair 仍不建文件；无临时文件），输出逐字节不变，内存 O(单 contig) | **已实施**（C） |
| V23-15 | ⑧4 | **属实（低概率）** | tql2 不收敛 → 至少 `debug_assert`/`Result`，或提高上限 | 待拍板 |
| V23-16 | ⑧3 | **代码顺序属实，误报可达性待定** | cell 就绪检查移到 watchdog 之前 | 待拍板 |
| V23-17 | ⑧11 | **属实（代码）** | `standardize` 用 `is_finite()` 而非 `!is_nan()` | 待拍板 |
| V23-18 | ⑧10 | **属实（代码）** | kernel PCA 负特征值 `max(0.0).sqrt()` | 待拍板 |
| V23-19 | ⑨1 | **属实（实测丢 chr2）** | 同族 `-r`/`-R` 改为可重复并集，跨族 AND；组合矩阵测试 | **已实施**（随 V23-60/61 与 V23-49） |
| V23-20 | ⑨2 | **属实（实测 `-t` 被忽略）** | 跨族交集（两条 head 规则 AND）；补测试 | **已实施**（随 V23-61 与 V23-49） |
| V23-21 | ⑨3/⑨5 | **属实（实测文件名为 `-.xxx`/`.xxx`）** | 已拍板：选项 A——`stats/convert/pca/distance` 共用 `require_output_prefix`，拒绝 `-o -` 与空前缀；`popg` 检查上移；`filter` 的 stdout 语义不变 | **已实施**（A） |
| V23-22 | ⑤3 | **属实（仅库 API）** | `decay_max_bp=None` 时统一约定（不建 LD plan 或同用缺省） | 待拍板 |

## C. 建议 P2（语义/诊断/边界）

| 编号 | 来源 | 判定 | 建议修法 | 状态 |
|---|---|---|---|---|
| V23-23 | ④4 | **属实（语义不一致）** | `TAG=.`/空值统一走 `on_missing`；三形态测试 | **已实施**（随 V23-59；`.`/缺席统一走算符表，Number=A 同规则） |
| V23-24 | ④5 | **属实** | `to_args` 先匹配两字符算符；`< <= > >=` 重放测试 | 待拍板 |
| V23-25 | ④11 | **属实（JSON 边界）** | `is_some_and(!is_empty())` | 待拍板 |
| V23-26 | ④12 | **属实（HLA contig 限制）** | 先整串匹配 header contig 再切 `:` | 待拍板 |
| V23-27 | ⑥1 | **属实（库 API）** | `bin_of` 返回 Option/哨兵，首 edge 前跳过 | 待拍板 |
| V23-28 | ⑥4 | **属实（release 静默）** | 内部不变式改 `assert!` 或 `Result` | 待拍板 |
| V23-29 | ①6 | **属实（诊断）** | 错误信息带 `found/expected` | 待拍板 |
| V23-30 | ①7 | **属实（畸形输入）** | `count_alts("") == 0` | 待拍板 |
| V23-31 | ⑨4 | **代码复核（顺序）；已拍板（随 V23-21）** | `-o -` 检查上移到参数校验 | **已实施** |
| V23-32 | ⑨7/⑨8/⑨9 | **属实（低）** | `--format` 对表类报错/文档；错误带文件名；tree 文档 id 上限 2n−1 | 待拍板 |
| V23-33 | ⑧9/⑧13 | **代码复核** | diffusion 相对阈值；`pca::backend()` 文档 | 待拍板 |

## D. 性能 / 精简（不影响当前正确结果）

| 编号 | 来源 | 判定 | 建议 |
|---|---|---|---|
| V23-34 | ①3 | 属实（性能） | FORMAT 缓存 64 → LRU 或提高上限；先 profile |
| V23-35 | ②3/②4/②6/②7/②8 | 属实（性能/精简/潜伏） | 热路径 scratch、删死 API、leak 仅在未来逐记录容错时成问题 |
| V23-36 | ⑤5–⑤10 | 属实（性能/精简） | mid stage 缓冲、`auto` Vec、prune 位图、`PcaScratch`、`append_sites`、rm-dup HashSet |
| V23-37 | ⑥3/⑥6–⑥8 | 属实（性能/清理） | WcMulti scratch、PBS 冗余守卫、pair/multi 尺度统一、测试死代码 |
| V23-38 | ⑧5 | 属实（性能） | `Grm::mean_impute` 在 faer 下委托 `gram_divide`；补 backend 对拍 |
| V23-39 | ⑧7/⑧12 | 属实（性能） | 矩阵输出 BufWriter；LLE 稀疏化 |
| V23-40 | ⑨ Lean | 属实（精简） | `filter.rs` 与 `common.rs` 的样本来源冲突逻辑合并 |
| V23-41 | ⑨6/⑨10 | 基准方法 | vs-noodles 的 BCF 分支改 `select`；脚本 warmup/交替顺序/保留 stderr |

## E. 测试缺口（修代码前/后补护栏）

| 编号 | 来源 | 内容 |
|---|---|---|
| V23-42 | ①4 | `pruned_header` 保序 + 短 keep |
| V23-43 | ②5 | BCF 畸形类型/扩展计数的失败测试 |
| V23-44 | ③9/③10 | IUPAC 全表断言；`pretty_sequence(n≠10)` 行为 |
| V23-45 | ④9 | reed-filter 的 `TAG=.`、`to_args`、深度上限（AC 向量测试随 V23-07 丢弃） |
| V23-46 | ⑥5 | WcMulti r≥3 半调用组、F2/F3/F4 merge 律 |
| V23-47 | ⑦1–⑦5 | 拍板：oracle 保持 complete-case（与生产/spec 一致）；补含缺失位点 golden；文档记录与 allel 的逐位点删除差异。其余 oracle 项（WC84 全缺失实测、hobs/hexp 口径）仍待处理 |
| V23-48 | ⑧8 | faer vs fallback、1/4 线程逐值对拍 |
| V23-49 | ⑨1/⑨2 | `-r/-R/-t/-T` 组合矩阵 |

## F. 驳回 / 不建议改

| 编号 | 来源 | 判定 | 理由 |
|---|---|---|---|
| V23-50 | ④1 | **不成立（已拍板：不改语义）** | bcftools 对 AN=0 的 AF 同样是 NaN，`AF!=0` 同样保留；改 Missing 反而偏离参照实现。落地：文档写明 NaN 语义 + 测试钉住，不改代码行为 |
| V23-51 | ⑧1 的 P0 定级 | 定级修正 | 仅库 API，reed CLI 无调用方（P1 而非 P0） |
| V23-52 | ⑥2 的"全 contig" | 表述修正 | `finish_contig` 会 `clear()`；按单 contig 计算 |
| V23-53 | ⑤3 的 CLI 可达性 | 表述修正 | CLI 拒绝 pairs-only；仅库 API 的 `decay=None` 组合 |
| V23-54 | ②2 的 filter 路径 | 表述修正 | `filter` 有 GT 整数类型守卫；panic 在 convert/stats 渲染 |
| V23-55 | ③5/③11/②7/②8 | 加固/低优先 | 真实基因组不可达或仅畸形输入；不列为 bug |

## 决策点（第二批，2026-11-03 讨论后拍板）

> 讨论依据：本机实测 bcftools 1.24-2 与源码核对（`bcftools/filter.c:4449-4455`
> 的 `missing_logic`、`CMP_VECTORS`），见 [01-verification.md](01-verification.md)
> §"补充实测"。

### 决策 A：缺失值与 `!=` —— 照 bcftools

**问题比复核文档更严重**：复核只记了"`TAG=.` 绕过 `on_missing`"，实测发现
**表达式路径的 `!=` 与 bcftools 方向相反**：

```text
record: DP=. (POS 10) / DP 缺席 (POS 20) / DP=5 (POS 30)
-i 'DP!=3'      bcftools=[10 20 30]   reed=[30]      ← 方向相反
-i 'INFO/DP!=3' bcftools=[10 20 30]   reed=[30]      ← 同上
-i 'DP>3'       bcftools=[30]         reed=[30]      ← 一致
-e 'DP!=3'      bcftools=[]           reed=[10 20]   ← 反向多保留
```

机制（`filter.c:4449-4455`）：`missing_logic` 对 `==` 是 `{0,0,1}`（两侧都
缺失→真、一侧缺失→假），对 `!=` 是 `{0,1,1}`（**一侧缺失→真**、两侧都缺
失→真）。reed 把 `Missing` 在 `cmp_vectors` 里 `continue` 掉，于是 `!=`
丢掉缺失记录。

**拍板：照 bcftools**——`!=` 在一侧缺失时为真，`==`/`<`/`<=`/`>`/`>=`
在任一侧缺失时为假（`==` 仅在两**侧都**缺失时为真）。同时修正 V23-50 的
措辞：`docs/spec/cli.md` 的"与缺失值比较为假"改为按算符分列的真值表。

落地：`V23-59`（见下表）。

### 决策 B：`-r/-R` 与 `-t/-T` —— 跨族交集；同族重复取并集（有意分歧）

**先纠正讨论时的三处误读**（源码 `bcftools/vcfview.c:649-652,791-795` +
重新实测）：

1. `-r`/`-R` 与 `-t`/`-T` 各自写**同一个变量**（`regions_list` /
   `targets_list`）。同一族给两次时**后写覆盖**（`-r chr1 -r chr2` → 只有
   chr2），且 `is_file` 标志不复位——`-r chr1 -R f` 只剩 `f`，
   `-R f -r chr1` 会把 `chr1` 当文件名打开（实测报
   `Could not open file "chr1"`）。所以"并集"是我先前的误判；bcftools 的
   同族重复是**顺序相关的覆盖**，属于实现瑕疵。
2. `-T` 不是"排除"，而是 targets **文件**（包含）；排除靠条目上的 `^` 前缀。
3. 排除与包含同族并存时，bcftools 是"包含列表 + 排除列表"的集合差
   （`-t chr1 -t ^chr1:10` 两头都生效），不是"排除存在时忽略包含"。

**拍板**：

| 组合 | 语义（reed） | bcftools |
|---|---|---|
| `-r/-R` + `-t/-T` | **交集**（两条 head 规则 AND） | 一致 |
| 同族重复（`-r`+`-R`、`-t`+`-T`、多次 `-t`） | **并集**，与顺序无关 | 后写覆盖、顺序相关 |
| `-T` 的 `^` 条目 | 排除；与 `include` 取集合差 | 一致 |

同族重复处 reed **有意分歧**：复制"后写覆盖 + `is_file` 泄漏"等于把一个
顺序相关的静默丢数据行为固化成契约，与仓库既有原则（不把 htslib 的偶然
行为当规范）一致。已在 `docs/spec/cli.md` 的分歧表登记。

落地：`RegionRule` 新增 `indexed`，`-r/-R/-t/-T` 改为可重复的 `Vec`，
`compile` 生成两条 `HeadRule::Regions`；`to_args` 把 `indexed` 渲染回 `-r`。
测试 `V23-49`（9 组组合 + 同族并集）。

**实测对拍**（bcftools 1.24-2，`two.vcf.gz`：chr1:10、chr2:10）：

| 命令 | bcftools | reed |
|---|---|---|
| `-r chr2` | chr2 | chr2 |
| `-r chr1,chr2 -t chr2:10` | chr2 | chr2 |
| `-r chr2 -t chr1:10` | 空 | 空 |
| `-t ^chr1` | chr2 | chr2 |
| `-t chr1 -t ^chr1:10` | 空 | 空 |
| `-T c2only.tab` | chr2 | chr2 |
| `-r chr1,chr2 -T c2only.tab` | chr2 | chr2 |

### 决策 C：本轮范围

按讨论的排序：①缺失值 `!=`；②`-r/-R/-t/-T`；③内存三项
（V23-10/11/12）；④低成本确定性项（V23-13/15/16/17/18/22/27/28/29/30/
32/33）；⑤测试护栏（V23-42/44/46，45/48 复核覆盖）；⑥只做 V23-38/40，
其余性能/精简项挂起（先 profile 再动）。

### 新增条目

| 编号 | 来源 | 判定 | 落地 | 状态 |
|---|---|---|---|---|
| V23-59 | 讨论 A | **属实（比 V23-23 更严重：表达式路径 `!=` 方向相反）** | 表达式与 head 规则的缺失比较按 bcftools 的 `missing_logic` 真值表（含 Number=A 的 `.` 元素）；修 V23-50 的文档措辞；8 项测试 | **已实施** |
| V23-60 | 讨论 B | **属实（`-r` 存在时 `-R` 被忽略）** | 同族（`-r`+`-R`、`-t`+`-T`、多次 `-t`）可重复并集，与顺序无关（有意分歧，见决策 B） | **已实施** |
| V23-61 | 讨论 B | **属实（`-t` 被 `-r` 覆盖）** | 跨族交集（两条 head 规则 AND）；`-T` 是包含文件，`^` 才是排除 | **已实施** |

## 决策点（已拍板）

1. **V23-09 特征向量符号**：已拍板，按结果正确性改（`SymEigResult::top_k`
   与 `eig::SymEig::top_k` 统一规范化）；接受 `.eigenvec`/坐标符号变化，
   换取跨 backend 可复现；需 V23-48 双 backend 对拍。
2. **V23-14 `--ld-pairs`**：已拍板，选项 C——sink 按 contig 顺序直接追加
   最终 `.ld`；首个 pair 懒创建并写表头（全无 pair 仍不建文件），无临时
   文件、输出不变、内存 O(单 contig)。
3. **V23-07 `AC` 多等位**：已拍板，**丢弃（不修）**——不接受多等位，该条
   从 v23 修复清单丢弃；默认门已丢弃多等位，显式 `-M >2` 时 AC 按首 ALT
   属不支持路径，不实现 any-pair。
4. **V23-50 AF/NaN**：已拍板，保持 NaN 语义不变，补文档/测试；不改
   bcftools 兼容行为。
5. **V23-21 `-o -`**：已拍板，选项 A——前缀型命令统一拒绝 `-`/空前缀
   （usage），`filter` 的 stdout 语义不变；V23-31（popg 检查上移）随此
   落地。
6. **V23-47 oracle SFS**：已拍板，保持 complete-case，补缺失 golden 并
   文档记录与 allel 的差异。

其余未列入上表的条目仍待拍板。

## 实施记录（2026-11-03，已完成）

A 组 6 条（V23-01…06）+ V23-08/09 + 拍板 5 项（V23-14/21/31/47/50）+ 第二批
（V23-59/60/61 与内存三项 V23-10/11/12）+ V23-13 已落地；V23-07 按拍板丢弃。
逐条实现与测试见 [02-implementation.md](02-implementation.md)。验收：

| 检查 | 结果 |
|---|---|
| `ci/test-parallel.sh` | 905 passed / 0 failed（56 二进制，15 s；v22 基线 867） |
| clippy `-D warnings` | 0 告警 |
| `ci/check-doc-drift.sh` | 97 项全 ok |
| 真实 `subvcf.vcf.gz --summary` | DP/PSC 与 `bcftools stats -s -` 逐行相同 |
| 真实 BCF 负浮点 | `reed filter -Ov` == `bcftools view -Ov` |

三条实施中的**偏离/修正**（均有实测依据）：

1. **V23-05 的上限从 500 改为 128**：debug 构建下 `parse + Drop` 在 2 MiB
   线程栈上 ~430 层即溢栈、1 MiB 上 ~200 层；拍板文档的"如 500"不安全。
   同时新增 `MAX_CHAIN=2000` 约束左深 `||`/`&&` 链（Drop 同样递归）。
2. **V23-03 未改公开签名**：`n_samples` 直接从 plan 的
   `SampleMask`/`Bitset`（都由 header 构造）读取，三个调用点不动，避免把
   签名变更扩散到 `reed`/`reed-io` 的全部调用方。
3. **V23-14 用新增字段而非改语义**：`LdRequest::pairs_path` 为 `None` 时
   仍是内存累积（库调用方行为不变），CLI 的 `--ld-pairs` 才设路径。

第二批的两处**自查修正**：

1. **V23-59 的 `Number=A` 分支**：第一版只加了注释，`.` 元素仍 `all =
   false`，`!=` 不通过；已改为按算符表计入 `any`/`all` 并补测试。
2. **V23-60 的落地文本**：新增条目原写"同族重复 usage 错误"，与决策 B 的
   并集规则矛盾；以决策 B 为准，已改正。

其余待拍板条目（V23-15…20、22…30、32…46、48、49）未改代码；其中
V23-19/20/23 已随本批落地。
