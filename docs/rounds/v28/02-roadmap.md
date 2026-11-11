# 02 实施路线与验收门

> 原则：先修会静默产错/覆盖旧结果的 P0，再修兼容性；下沉只在同一批修复需要共用底层能力时落地，
> 不做大一统 `Record`/事务层。每条修复必须先有反向测试或独立夹具证明修复前失败。

## 批次

### Batch 1 — P0 数据正确性（V28-01…V28-05）

1. **V28-01 + V28-22 的 contig 部分**：BCF 输出 header 写 contig `IDX`，与 `StringMaps` 共用字典。
   - 验收：显式 `IDX=1/0` 的 chr2/chr1 输入，BCF 回读仍对应原染色体；`IDX=7` 可被 bcftools 读。
2. **V28-02**：BCF 的 CSI 参考序列表改用 BCF 字典顺序（记录 tid 原生空间），不再用 header 行序。
   - 验收：`idx.bcf` 的 chr1/chr2 region 读取各只返回自己的两条记录；CSI 仍能被 bcftools 正确查询。
3. **V28-03**：region reader 对短 BGZF header 报 `MalformedBgzf`。
   - 验收：截断第二块时 exit 非 0；完整输入与合法 EOF 仍正常。
4. **V28-04 + V28-21 的 GT 部分**：共享 GT 元素校验，拒绝负 code。
   - 验收：INT8/16/32 负 code 均 `MalformedRecord`，vector-end 仍截断；`0/0` 不受影响。
5. **V28-05**：integer vector-end 在 float 列转 `NaN`。
   - 验收：renderer、整数列、float 列三者对 missing/vector-end 一致。

### Batch 2 — P0 输出发布（V28-06）

- `OrderedOutput` 增加 `finish_pending()`；filter sink 暂存 atomic + indexer，`job.finish` 通过后发布。
- `IndexSink` 暂存 `Indexer + path`，BCF `chain.finish` 通过后写索引。
- `ReportSink` 暂存所有 `AtomicSink` 的 atomic owner，`job.finish` 通过后统一 flush/close/commit。
- 验收：截断 BCF + 预写 sentinel 旧输出；命令失败后旧数据/旧索引/旧报告逐字节不变、无临时残留。
- 边界：逐文件原子，不做多文件事务；stdout 不涉及发布。

### Batch 3 — P1 兼容性（V28-07…V28-11）

- phasing：`|0/1`、`|.` 与完整 renderer 对齐；`/0|1` 保持现有等价。
- NULL GT：`(NULL,0)` 按 haploid missing；其他非整数仍拒绝。
- half-call：快速路径尊重非默认策略。
- zero allele：`encode_bcf` 在写任何字节前拒绝。
- CRLF：无编辑/局部编辑都保留原 terminator。

### Batch 4 — P2 测试与清理（V28-12、V28-13…V28-20）

- header 上限用可注入上限的测试读取器真正越界。
- 删除失实注释、无失败 `Result`、未用 `Default`、不可达截断分支。
- INFO 索引复用、多 span region、BCF→PCA 原生矩阵、FORMAT 表复用、样本区单趟切分、
  plain chunk 借用、Float scratch 复用（V28-13…17、19、20）随本批落地，
  各自以“修复前失败”的回归或既有 parity 测试为门。

### Batch 4b — V28-28（实施期新发现的 P0）

- 索引参考表必须与记录 tid 同空间：`filter -Ob -W` 用字典序，`reed index` 同步；
  解除 noodles 的“参考 id 必须递增”限制（htslib 只要求同染色体连续）。
- 验收：CSI 自身参考表 + `bcftools -r` 双 contig 精确匹配；`idx.bcf` 的 tid 递减可索引。

### Batch 5 — 下沉（V28-21…V28-27，已全部落地）

按以下顺序，每次只删除已有重复，不新增配置：

1. V28-21 BCF shared/FORMAT/GT 原语（Batch 1 已依赖）。
2. V28-22 contig vocabulary（Batch 1 已依赖）。
3. V28-23 AtomicWriter（Batch 2 已依赖）。
4. V28-24 `LineRef` 访问器。
5. V28-25 BCF framing helper。
6. V28-26 Input opener。
7. V28-27 VirtualOffset/BlockTable。

## 验收门

- `cargo test -p reed-io`
- `cargo test -p reed`（至少 `e2e_bcf`、`e2e_filter`、`e2e_index`、`e2e_reports`）
- `cargo clippy -p reed-io -p reed-filter -p reed-engine -p reed -p reed-cli --all-targets -- -D warnings`
- `cargo fmt --all -- --check`
- 截断 BCF 夹具下旧输出/索引/报告 sentinel 不变（反向验证：修复前失败）
- 显式 contig IDX 的 VCF→BCF→bcftools 回读逐条 CHROM 正确
- CSI 自身的参考序列表与记录 tid 同空间（V28-28 的直接断言，不再只靠自读回环）
- `idx.bcf` 的 chr1/chr2 region 查询集合精确匹配
- 修复不改变普通双等位、无显式 IDX、LF 输入的正常输出

## 不做

- 不把 LD 线程池问题再修一遍：当前 `LdExecutor` 已是单 pool，V27-16 已落地。
- 不为多文件输出制造事务/回滚假保证。
- 不把统计策略、过滤策略、调度策略塞进 reed-io。
- 不以“测试全绿”代替边界测试覆盖。
