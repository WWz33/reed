# 01 复核清单：reed-io 内部正确性、输出发布与下沉边界

> 立项 2026-11-10。审查基于当前工作树（HEAD `7e37a05` 加未提交修改），
> 不是干净提交；工作区在长审期间发生过并发修改，因此每条均已按当前源码复核。
>
> 标记：`confirmed` = 已用当前源码/独立夹具复现；`code` = 确定性代码路径可证明，
> 但未构造端到端夹具；`design` = 工程收敛建议，不是缺陷；`withdrawn` = 上轮主张经复核不成立。

## A. P0：会静默产生错误数据或覆盖好结果

| ID | 状态 | 位置 | 触发与证据 | 最小处理方向 |
|---|---|---|---|---|
| V28-01 | **confirmed** | `reed-io/src/bcf/encode_noodles.rs:195-233`、`reed-io/src/bcf/encode_native.rs:79-124` | `bcf_header_text` 只给 FILTER/INFO/FORMAT 行补 `IDX=`，没给 `##contig` 补。编码器却按 `StringMaps::contigs()` 的字典索引写记录。输入 `##contig=<ID=chr2,IDX=1>` 后接 `##contig=<ID=chr1,IDX=0>` 时，VCF→BCF 和 BCF→BCF 均交换染色体；`IDX=7` 生成 htslib 拒绝的 header。 | 统一 BCF 输出 header 与记录编码的 contig 字典：序列化 `##contig` 时写实际字典索引，或统一重编号并重映射记录。补显式 `IDX=1/0` 与 `IDX=7` 测试，用 bcftools 回读。 |
| V28-02 | **confirmed** | `reed/src/index.rs:329-342`、`:392-414`、`:606-621` | CSI 的 name/rid 表按 VCF header 行序建，记录里的 tid 却来自 BCF 字典。使用仓内 `idx.bcf`（header 行序 chr1、chr2；字典 `0→chr2`、`1→chr1`）建索引后，bcftools `-r chr1` 返回 chr1 和 chr2，`-r chr2` 为空。 | BCF 的 CSI 参考序列表必须按 **BCF 字典 tid 顺序**建（不是 header 行序，也不是重映射记录 tid）；TBI 仍只用于文本输入。用 `idx.bcf` 的 region 读断言 chr1/chr2 各两条。 |
| V28-03 | **code / confirmed CLI** | `reed-io/src/chunk/region.rs:99-106` | `read_block` 把任意 `Ok(0)` 当正常 EOF。索引覆盖的第二个数据块只剩 10 字节时，`-r chr1:1-20` 应返回 POS 10、20，实际 exit 0 只返回 10。 | `filled == 0` 才是 EOF；`0 < filled < 18` 返回 `MalformedBgzf("stream ends inside a block header")`。补短 header、短 payload、合法 EOF marker 三类测试。 |
| V28-04 | **confirmed** | `reed-io/src/bcf/render.rs:424-440`、`gt.rs:169-203`、`encode_native.rs:958-1016` | GT 元素除 vector-end 外必须非负。`fmt_gt` 对负数 `unwrap_or(0)`，`decode_pair` 把负值当 missing 放行。INT8 `[0x80,0x80]` render 成 `0/0`，native 得到 REF，BCF→BCF 原样保留负 code。 | 建立共享 GT 元素校验：先识别 vector-end，再拒绝其余 `v < 0`；render、table、encoder 共用。补 INT8/16/32 回归。 |
| V28-05 | **confirmed** | `reed-io/src/bcf/columns.rs:158-177` | 整数转 float 列只处理 missing sentinel，没有处理 vector-end。INT8 `0x81` 在文本路径是缺失，direct integer 列是 `INT_MISSING`，float 列却是 `-127.0`。 | integer-backed float 同时把 missing 和 vector-end 转 `NaN`；补 render/direct parity 测试。 |
| V28-06 | **confirmed** | `reed/src/filter/mod.rs:202-219`、`:1052-1063`；`reed/src/index.rs:533-554`；`reed/src/stats/writers/reports.rs:1125-1248` | 数据/索引/报告在 `job.finish` 的 BCF 尾部完整性检查前提交。截断 BCF 返回错误，但旧 filter 输出、旧 CSI 或旧 `.afreq` 已被部分新结果替换。 | 拆出“结束写入/flush”与“发布 commit”。filter 的 `OrderedOutput`、index 的 `Indexer::write`、reports 的 `AtomicSink` 先暂存，`job.finish` 成功后统一发布。保持逐文件原子，不承诺多文件事务。 |

## B. P1：正确性/兼容性缺口

| ID | 状态 | 位置 | 触发与证据 | 最小处理方向 |
|---|---|---|---|---|
| V28-07 | **confirmed / 例子已更正** | `reed-io/src/bcf/gt.rs:94-138`、`:209-233`；`render.rs:446-457` | `GtTable` 固定按 `vcf44=false` 构造，完整 renderer 按 header 输出首等位前缀。真差异是 `GT=|0/1` 和 `GT=|.`，不是上轮写的 `/0|1`（后者两路径等价，`genotype_scan.rs:589-590` 已覆盖）。 | 查表后按原始 int8 pair 与 `BcfDict::vcf44()` 重建首等位 phasing；slow path 同享。补 `|0/1`、`|.`、`/0|1` 三组 parity。 |
| V28-08 | **confirmed (htslib interop)** | `reed-io/src/bcf/render.rs:408-417`、`gt.rs:217-227`、`encode_native.rs:938-943` | 合法/htslib 可读的 `(NULL,count=0)` GT 被三处拒绝；bcftools 将其读成所有样本 `.`。 | 仅对 `(ty==NULL && count==0)` 按一个 haploid missing allele 处理；其他零长类型仍拒绝。 |
| V28-09 | **code (公开 API)** | `reed-io/src/vcf/genotype_scan.rs:349-368` | 3 字节 GT 快速路径直接 `from_states`，绕过 `HalfCallPolicy::TreatAsMissing`。`0/.` 在 `TreatAsMissing` 下仍保留 REF，计数为 0。当前 CLI 未启该非默认策略。 | 策略无需改写时保持快路径；否则把两个 slot 交给 `from_alleles`。补两种策略与计数测试。 |
| V28-10 | **code (底层 API)** | `reed-io/src/bcf/encode_native.rs:178-208` | `site_header`/renderer 拒绝 `n_allele==0`，`encode_bcf` 自行解包 shared header，不检查。独立 harness 对零等位 record 成功输出。正常 filter 上游会先拒绝。 | `encode_bcf` 复用 `site_header` 或等价校验，在写 shared/indiv 前拒绝。 |
| V28-11 | **confirmed** | `reed-io/src/vcf/line.rs:105-107`、`line_edit.rs:86-133` | `LineRef` 从字段剥离 CRLF，editor 所有分支补 LF；默认空 `EditPlan` 也把 `0/1\r\n` 改成 `0/1\n`。 | 保存并复制原行尾；完全无修改时整行复制。补 8/9/10 列 CRLF 测试。 |

## C. P2：测试、重复扫描与清理

| ID | 状态 | 位置 | 问题 | 最小处理方向 |
|---|---|---|---|---|
| V28-12 | **confirmed** | `reed-io/src/chunk/plain.rs:127-180`、`chunk/bgzf.rs:911-933` | header 上限测试只到约 2 MiB，生产上限 256 MiB；删掉上限检查仍通过。 | 内部读取函数接收测试上限；断言 `header exceeds`，并补跨 BGZF 块合法大 header 正例。 |
| V28-13 | **code** | `reed-io/src/vcf/info.rs:104-150`、`:192-240`；`reed-filter/src/eval/head.rs:70-83` | INFO 过滤先 `validate_info` 再 `InfoIndex::build`，同记录两遍查重，两套规则需同步。 | 仅对需要 INFO 索引的计划预建并复用；普通分类路径保持无分配。 |
| V28-14 | **code** | `reed-io/src/chunk/region.rs:242-280` | 多 region 折叠为 min..max，远端稀疏查询读取中间全部块。 | 返回排序/合并后的非连续压缩区间，分别 seek；保持 chunk id 与 contig 顺序。 |
| V28-15 | **code** | `reed/src/bcf_input.rs:297-355`、`reed/src/genotype.rs:375-384` | PCA/DAPC/distance 的 BCF 输入先完整文本化，再重扫 GT。 | 复用 `InputJob::with_bcf_records` 与 `decode_record_gt`，文本路径不变，补 BCF/VCF 矩阵等价。 |
| V28-16 | **code** | `reed-io/src/bcf/gt.rs:285-358`、`columns.rs:42-115` | `decode_record_gt` 与 `decode_format_columns` 各自遍历 FORMAT，后者每条记录分配字段表。 | 复用 `NativeScratch` 字段表，一次遍历定位 GT/主键/fallback；字典 id 在 plan 初始化时解析。 |
| V28-17 | **code** | `reed-io/src/vcf/line.rs:258-285`、`genotype_scan.rs:163-252` | 文本样本区在分类、扫描、编辑阶段重复 O(样本区字节) 遍历。 | 在 `LineRef` 缓存边界并给 scanner 提供“结构已验证”入口；直接 scanner 调用仍完整校验。 |
| V28-18 | **confirmed** | `reed-io/src/bcf/encode_noodles.rs:1-7`、`:153-159`；`bcf/record.rs:126-138` | 旧泄漏注释已失实；`new()` 永失败不了的 `Result` 与未用 `Default` 是 YAGNI；`record_ref` 第二次截断检查不可达。 | 修注释；`new` 返回 `Self`；删除 `Default` 与死分支。 |
| V28-19 | **code** | `reed-io/src/chunk/plain.rs:90-118`、`chunk/mod.rs:181-191` | 未压缩 VCF chunk 在 reader 与 worker 各复制一次。 | plain chunk 直接借用 `ChunkPayload`；reader 复用 `Arc<Vec<u8>>` 容量。 |
| V28-20 | **code** | `reed-io/src/bcf/encode_native.rs:364-394` | VCF→BCF 的每个 Float INFO 值分配临时 `Vec`。 | 复用 `Vec<f32>` 或扫描时直接写 `to_le_bytes`。 |

## D. 下沉（design，不是缺陷）

| ID | 能力 | 当前重复 | 边界 |
|---|---|---|---|
| V28-21 | BCF shared-header / FORMAT walk / GT element 基础解码 | `site_header`、renderer、`encode_bcf` 各自解 shared；`gt.rs` 与 `encode_native.rs` 各自 walk FORMAT | 只拥有 framing、type width、sentinel、allele/phasing；策略仍在上层 |
| V28-22 | Contig vocabulary（header 顺序 / BCF `IDX` / index tid） | header zip、`ContigTable::new`、rename 手扫、index tid | 三种 id 不互换；不拥有 region/window 语义 |
| V28-23 | 拥有 `BufWriter` 的 `AtomicWriter` | convert 私有实现；filter/stats 各自 flush→close→commit | 何时发布仍由调用方决定；只保证单文件原子 |
| V28-24 | `LineRef` 样本列 / FORMAT 子列 / field span | genotype、convert、reports、rename 各写 tab/colon 游标 | 只返回借用字节和范围，不决定 dosage/缺失/样本过滤 |
| V28-25 | 跨 chunk BCF 长度拼接与完整记录边界 | `bcf_input.rs` 与 `index.rs` 各一份 | `Carried`、并发链、virtual offset 仍留上层 |
| V28-26 | `Input::open`：sniff + source + sidecar | `input.rs` 文件/stdin 四分支近复制，region 打开三次 | stdin replay、region predicate、调度仍在调用方 |
| V28-27 | `VirtualOffset` / `BlockTable` | `VoMap`、`VirtualMap`、index writer 手写 shift/OR | 只表达 BGZF 地址，不知道索引策略 |

## E. 撤回或降级的上轮主张

| 上轮主张 | 裁决 | 当前证据 |
|---|---|---|
| LD decay 每条 contig 建 Rayon 池，可产生数千线程 | **withdrawn** | 当前 `reed/src/stats/sink.rs:377-430` 已是单个 `LdExecutor` pool；内部块并行使用同一 pool。对应 V27-16 已修，不能再作为本轮任务。 |
| `GT=/0|1` 在 native 与文本路径得到不同 GtCode | **例子不成立** | `/0|1` 两路径等价；真差异是 `|0/1` 与 `|.`。V28-07 已更正。 |
| `NULL,count=0` 是规范保证的正常 GT | **降级** | 当前定位为 htslib 兼容边界，不是主路径保证；仍应支持但列 P1。 |
| VCF4.4 phasing 会污染现有统计 | **降级** | 目前 `GtCode.phased()` 只进入 genotype planes，现有统计不消费相位；按 parity/API 修，不称数值错误。 |
| CRLF 规范化是数据错误 | **降级** | 语义不受影响；按“未编辑区域逐字节保留”契约修，列 P1/契约问题。 |
| 所有未覆盖模块都无 Bug | **不成立** | 本轮深扫新增 V28-01…V28-11；测试全绿不构成否定证明。 |

## F. 复现与验证环境

- 子代理用 `/tmp` 独立 harness 复现：contig IDX 交换、负 GT、vector-end、零等位、CRLF、half-call。
- 用仓库 `crates/reed-io/tests/fixtures/bcf/idx.bcf` 与 bcftools 1.24 复现 region/index tid。
- `cargo test -p reed-io`：164 unit、17 integration/property、1 doctest 通过；这些边界无现有覆盖。
- 未跑完整 workspace、6.6 GB 全量数据、真实大 BCF 压测或全部 gold 对拍。

## F. 实施期新发现（V28-28）

| ID | 状态 | 位置 | 触发与证据 | 处理方向 |
|---|---|---|---|---|
| V28-28 | **confirmed（bcftools 交叉验证）** | `reed/src/filter/mod.rs:742-755`、`reed-io/src/write/index.rs`、`reed/src/index.rs:615-646` | V28-02 只修了 `reed index`。`reed filter -Ob -W` 仍用 `##contig` 行序建 CSI 参考表，而原生编码器按 `StringMaps`（显式 `IDX=`）写记录 tid。用反序 `IDX=`（chr2=1、chr1=0）的输入：记录 tid 为 1/0，CSI 名表为 `[chr2, chr1]`，`bcftools -r chr1` 与 `-r chr2` 都返回 **0 条**（`bcftools view` 全读正常）。reed 自读走同一张错表，故既有 e2e 全绿——原先的回归用了 `-Oz`，文本输出丢弃 `IDX=`，行序恰好正确。**同一处还叠着第二个限制**：noodles 的 `Indexer` 拒绝参考 id 递减，而 htslib 只要求同一染色体的记录连续；仓内 `idx.bcf`（chr1 IDX=1 在前、chr2 IDX=0 在后）因此让 `reed index` 直接报错，`bcftools index` 却能建。 | ① `VcfHeader::contig_names_in_dict_order` 作为唯一字典序来源，filter 与 `reed index` 都按它建表（并保留稀疏 `IDX=` 的空槽）。② `Indexer` 改为每 contig 一个 noodles indexer、写盘时按 ref id 顺序合并，解除跨 contig 递增约束。 |

**反向验证**：`bcf_index_reference_table_follows_the_record_tid_space` 在旧 writer 上失败（`tid_of("chr1")` 得到 `Some(1)`，期望 `Some(0)`）；`index_accepts_descending_record_tids_like_htslib` 在旧 indexer 上直接报 `out of order reference sequence ID`。两者现均通过，并用 `bcftools -r chr1`/`-r chr2` 外部核对。
