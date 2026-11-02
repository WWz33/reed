# doc_v22 后续路线：架构收口（2026-11-03）

> 触发：doc_v22 修复完成后的架构复盘（INFO/FORMAT 解析的集中度、文本/BCF
> 平行实现、隐式策略）。本文件记录**决策**与**待办**；完成一项就更新
> STATUS.md 的验收列。

## 决策 D22-1：DP/AD 深度回退放在哪一层

### 上游事实（本机实测）

| 实现 | 行为 |
|---|---|
| htslib | 库，无深度策略；`bcf_get_fmt`/`bcf_get_format_int32` 只返回字段本身 |
| bcftools `stats -s -` | `calc_sample_depth`（`vcfstats.c:958`）：记录 FORMAT **有 DP 键** → 该样本 DP（值为 `.` 也不回退）；**无 DP 键** → 该样本 AD 向量之和；`dp > 0` 才入直方图（`:1119`） |
| bcftools 表达式 | `bcftools view -i 'FMT/DP>=10'` 在 header 未声明 DP 时直接报错 `Error: the tag "DP" is not defined in the VCF header`；`FMT/AD` 需要 `AD[sample:subfield]` 下标，**不做**隐式求和 |
| vcftools `--minDP` | `filter_genotypes_by_depth` 只在 `have_depth_info`（FORMAT/DP 存在）时生效；`get_indv_DEPTH` 无 DP 返回 `-1`，无 AD 回退 |
| plink2 `--vcf-min-dp` | 同上，严格 FMT/DP |

结论：**"DP 缺失用 AD 之和"只是 bcftools `stats` 的深度统计口径，不是
"FMT/DP" 这个字段名的通用语义。**

### 决定

1. **显式 `FMT/DP` 表达式（`-i/-e` + `--set-gts .`）= 字面 DP**，与 bcftools
   一致：header 未声明 DP → 规则编译期报错；记录 FORMAT 无 DP 键 → 该样本
   缺失。**不**回退。
2. **深度掩码/深度报告 = 显式的"深度"概念**：
   - `filter --min-dp/--max-dp`（plink2/vcftools 语义的 reed 扩展）；
   - `stats --depth`、`--summary` 的 DP/PSC、`--hist dp`（bcftools
     `calc_sample_depth` 口径）。
   这些消费者显式构造 `IntSpec::depth()`（`DP`，回退 `AD` 求和）；`--min-dp`
   由 spec 字段 `GenotypeMaskRule.depth = true` 标记（CLI 设置，JSON spec
   可显式开启），**不再按字段名在 IO 层特判**。
3. **删除 `format_layout::spec_for_key`**：IO 只接受显式 `IntSpec`；策略在
   spec/plan（filter）或调用点（stats）决定。
4. 一个 plan 同时含 `--min-dp` 和显式 `FMT/DP` 表达式时，生成**两列**：
   DP-with-AD-fallback（掩码）与 DP-strict（表达式）；`int_fields` 的去重键
   从"名字"改为 `(名字, 回退)`。

### 兼容性

- 有 DP 的 callset：`--min-dp/--max-dp`、`--depth`、`--summary` 的输出逐位不变。
- 无 DP 的 BSA callset：`--min-dp` 仍可用（reed 扩展，文档已声明）。
- 显式 `FMT/DP` 表达式在无 DP header 上从"静默使用 AD 之和"变为报错——
  与 bcftools 对齐，属于修正而非回归。

## 待办（R22 系列）

| ID | 项 | 内容 | 验收 |
|---|---|---|---|
| **R22-1** | 显式深度 spec | 实现 D22-1：`GenotypeMaskRule.depth`、`IntSpec` 显式化（`Cow` 拥有型）、`FilterPlan.int_fields` 带 `(name, fallback)`、BCF `decode_format_columns` 收 `&[IntSpec]`、删 `spec_for_key` | 单测：掩码列有回退、表达式列无；`e2e_filter`/`e2e_bcf` 的 AD 回退测试仍绿；BSA 掩码数仍 1285 |
| **R22-2** | BCF 模块收口 | `bcf/mod.rs` 加"读/写/适配"模块地图；`render.rs` 拆出 `columns.rs`（`decode_format_columns` + 元素解码 + `int_sum_column`，共享 `indiv_fields` 描述符遍历）；`encode.rs`→`encode_noodles.rs`、`native_encode.rs`→`encode_native.rs`；`reed/src/bcf_input.rs` 在模块地图里点名 | `render.rs` 只保留渲染；全量测试/对拍绿 |
| **R22-3** | `LineRef` 列形态契约 | 显式 `Shape`（8 列 / FORMAT 无样本 / FORMAT+样本）或 `validate_sample_cols(n)`；consumer 不再各自数 tab；proptest 覆盖 8 列、9 列、尾随空列、超列、CRLF | proptest 不 panic，三形态语义固定；`v22_accepts_sites_only_haploid_and_trailing_tabs` 仍绿 |
| **R22-4** | 文本/BCF 编解码一致性 | proptest：随机记录（AD 宽度、`.`、短向量、DP 有无）→ BCF（native 编码）→ 原生列解码 == 渲染为文本后 `GenotypeScanner` 解码 | 差异 0；覆盖 DP/AD、GQ、PL 无关注 |
| **R22-5** | 头部事实 | 复核结论：`reed-io/src/header.rs:52` 已是 `TagTable` 的唯一构造点（noodles Header → TagTable），无需新抽象；在代码/文档中写清"解析→声明→编译"三段，避免误以为三份实现 | 文档；无代码改动 |
| R22-6 | 错误上下文（可选） | `MalformedRecord` 增加 FORMAT key / 样本号的结构化字段 | 暂缓，按需 |
| **R22-7** | 测试速度与并行 | ① `bcf_vs_noodles` 改用工作区上级目录的本地子集 `subvcf.bcf`（3000 记录 × 348 样本，8.6 MB；`REED_SUBVCF_BCF` 可覆盖，缺失时跳过），并改流式解压 + noodles 惰性 `Record`：单测 ~400 s → 4.5 s；② `reed-graph` 的深树测试改直接构造 caterpillar（原 2000×2000 UPGMA 占 23.5 s → 0.01 s）；③ 新增 `ci/test-parallel.sh`：build 一次后按测试二进制用 `xargs`/线程池并发（libtest 的并行只在单二进制内），整库 run 阶段 50 s → 15 s；④ 新增 `gold/scripts/run_parallel.py`，gold 对拍脚本的独立 case 可并发（`run_bcftools_af.sh` 40 case：顺序 ~40 s → 5.8 s 墙钟 / 40 s CPU）；⑤ gold 脚本默认输入改为 `subvcf.vcf.gz/.bcf`，全量 `clean.*` 只在显式 `INPUT=`/参数时使用 | `ci/test-parallel.sh`：当前 860 passed / 0 failed；`run_bcftools_af.sh` 40/40 PASS on subvcf |

## 进度（2026-11-03）

| 项 | 状态 | 备注 |
|---|---|---|
| R22-1 显式深度 spec | **完成** | `GenotypeMaskRule.depth`、`FilterPlan.int_fields: Vec<IntField>`、`IntSpec` 拥有型（`Cow`）、BCF `decode_format_columns(&[IntSpec])`、删除 `spec_for_key`；`depth_policy_tests` 4 项 + 全量回归绿 |
| R22-7 测试速度 | **完成** | 见上表；`bcf_vs_noodles` 的本地输入变更为 `/sri/home/sri2025201067/reed/subvcf.bcf`（缺失时跳过） |
| R22-2 BCF 模块收口 | **完成** | `bcf/mod.rs` 模块地图；`render.rs` 只保留渲染，列解码拆到 `columns.rs`；编码器改名为 `encode_native.rs` / `encode_noodles.rs` |
| R22-3 `LineRef` 形态契约 | **完成** | `Shape`（`SitesOnly` / `FormatOnly` / `WithSamples`）+ `validate_sample_cols`；扫描器、`PcaJob`、`ConvertJob`、`line_edit` 共用；2,000-case 列形态 proptest |
| R22-4 编解码一致性 proptest | **完成** | `crates/reed-io/tests/format_columns_vs_bcf.rs`：500 个随机 FORMAT/AD/DP 记录，native 列解码 == 原文扫描 == 渲染文本扫描 |
| R22-5 头部事实 | 完成（文档） | 复核后确认 `header.rs` 已是单点，无需改动 |
| R22-6 错误上下文 | 暂缓 | 可选，按需 |

## 不做

- 不把 `bcf/` 改成 `vcf/` 的 1:1 文件镜像：帧方式、读写方向、性能分层
  （noodles vs native）决定了不对称是合理的；用 R22-2 的模块地图和 R22-4 的
  一致性测试代替形式对称。
- 不改变 `--ploidy` 的默认值（`require-diploid`），`stats` 的
  `HaploidAsHalfCall` 保持。
- 不引入 `bcf1_t` 式统一记录模型（会放弃零拷贝快路径，与 reed 的性能目标冲突）。