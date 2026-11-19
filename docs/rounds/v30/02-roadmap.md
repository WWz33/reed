# 02 路线：修复优先级与验收门

> 本文件分两部分：**第一部分**是 sgkit 覆盖/科学符合性的 V30-01…09；**第二部分**是
> API 工程质量（V30-Q01…Q36）。
>
> 本轮（doc_v30）是**审计轮**：下面每条都**尚未实施**。实施时按本文件的验收门执行，
> 状态另记 [STATUS.md](STATUS.md)。ID 与 [01-review.md](01-review.md) 一一对应。

## 1. 优先级

### A. 必须修（正确性）

| ID | 问题 | 方向 | 验收门 |
|---|---|---|---|
| **V30-01** | PI_HAT 期望 IBS 系数错位（01 §2.3） | 把 `expected_ibs_coefficients` 改成 plink 1.9 的系数（含 `dxx1/dxx2/dyy1/dyy2` 与 `num_allelesf2/f3` 有限样本项，或至少在无限样本极限下 `e01=4pq(p²+q²)`、`e02=(p²+q²)²`、`e11=2pq`、`e12=p²+q²`）。决定是否复刻 Nei 修正：**要么完整复刻并与 plink 对拍，要么明确改成无修正的理论系数并在文档声明与 plink 不一致** | ① 新增对拍：一个已知无关对（如两个不同个体）PI_HAT 应 ≈0 而非 0.5；② `gold/` 增加 plink `--genome` 对拍脚本，`Z0/Z1/Z2/PI_HAT` 到 1e-4；③ 现有 `matrices.rs` 测试升级（当前只断言 `0≤PI_HAT≤1` 且 Z 和为 1，抓不到本缺陷） |

### B. 应该修（语义/契约）

| ID | 问题 | 方向 | 验收门 |
|---|---|---|---|
| **V30-02** | `filter -I` 文档与行为相反（01 §2.3） | 二选一并对齐文档：**(a)** 让 `-I` 真正"保留"AC/AN（`-I` 时把 AC/AN 移出默认删除清单，只保留用户显式 `-x` 的项）；或 **(b)** 若维持现行为，改 help/`filter.md` 为"`-I` 关闭重算；默认 `-x` 仍会删除 AC/AN，需 `-x ''` 保留"。**推荐 (a)**，因为 bcftools `view -I` 的语义就是保留原值 | ① `-s S -I` 输出中 AC/AN 保留且值为**原值**；② e2e 断言 `-I` 与 `-I -x ''` 等价；③ `filter.md`、`cli.md`、`--help` 三处一致 |
| **V30-03** | 多等位三路径不一致（01 §2.3） | 落实 V27-39 的"不做多等位"到 `stats`：QC 报告对多等位**拒绝或显式警告+计数**，不再静默 `CollapseAlt`；或反向决定"允许折叠"并把 `popg`/`filter` 也统一为折叠。三者必须一致 | ① 多等位夹具下三个命令行为相同（同拒绝或同折叠+同数值）；② `.afreq` 对多等位要么报错要么输出 per-ALT 列（消除单值列） |

### C. 文档收口（不改数值）

| ID | 问题 | 方向 | 验收门 |
|---|---|---|---|
| **V30-04** | f-stat `estimate=Σ per-site` 未说明 | 在 `popg.md` 的 f-统计小节显著标注"estimate/se 为**位点和尺度**（对齐 ADMIXTOOLS）；z 与比值不受影响；要位点均值须自行除以 `n_snps`" | `popg.md` 有该句；`--help` 的 `--f2` 行提示（可选） |
| **V30-05** | `--f4-ratio` 角色序 | 在 help 与 `popg.md` 写死 `A,B;C,D;X ↦ f4(A,O;B?)…` 的精确映射（A, O=第二项, C=第三项, B_ref=第四项, X=第五项），并给一个例子 | 文档映射与 `translate/popg.rs:671-681` 逐项一致 |
| **V30-06** | `f3*` 的 ΣT/ΣB 承诺不可达 | 要么接上归一化输出，要么从 `popg.md:242` 删掉该承诺 | 文档与 `--f3` 实际输出一致 |
| **V30-07** | MSA 全缺失→REF | 在 `convert.md` 的 `--msa` 小节加显著警示：全缺失基因型写 REF，非 `N` | 文档含警示；若改语义需另立 ID |
| **V30-08** | `identity_by_state` 命名易误导 | 在 `distance.md` 说明 reed 的 `ibs` 是 plink IBS（het-het=1），与 sgkit `identity_by_state`（het-het=0.5）**不是同一个量** | `distance.md` 有对照句 |

### D. 验证能力（可选）

| ID | 方向 | 验收门 |
|---|---|---|
| **V30-09** | 新增 **sgkit 数值夹具**：用 sgkit 对 `diversity`/`divergence`/`Fst(hudson)`/`Tajimas_D`/`ld_matrix` 生成小夹具，纳入 `gold/`，把 §2.1 的"公式层等价"升为"数值等价" | 新 `gold/sgkit/` 脚本 + 对拍容差表；缺口项（Nei、Garud_H）在脚本里显式 skip |

## 2. 明确不做（本轮画线）

- **不实现 GWAS/regenie/genee**（4 项）与**完整系谱模块**（8 项）：属独立产品，非本轮或近期范围。
- **不实现概率/混合倍性数据层**（`convert_probability_to_call`、`infer_*_ploidy`、`call_allele_frequencies`）：
  除非有明确的 GP/DS 输入需求，否则维持硬调用模型。
- **不改**已对拍通过的估计器数值口径（π/dxy/Fst-Hudson/WC84/PBS/Tajima/HWE/LD r²/f-stat/GRM/KING/phi/discordance）：
  它们是**有意**的 pixy/allel/plink2 口径，`--variant-only` 等属设计选择，不因"与 sgkit 不同"而改。
- **不重开** v17 的 distance 对拍结论、v27 的多等位拍板（V30-03 只把拍板落实到 `stats`）。

## 3. 实施纪律（沿用仓库既有）

1. 每条修复**先写失败测试/对拍再改**，并在 STATUS 记修前输出。
2. V30-01 是数值缺陷：必须给"无关对 PI_HAT≈0"的确定性回归，而非只断言范围。
3. 文档类（V30-04…08）改动需同步 `--help`（若涉及）与 docsite；`ci/check-doc-drift.sh` 若
   需要新守卫（例如"`popg.md` 含 f-stat 量纲说明"）一并加。
4. 验收统一出口：`cargo test --workspace` 全绿 + `clippy` 0 告警 + `ci/check-doc-drift.sh` +
   `ci/check-structure.sh`，并在 STATUS 附命令与输出。


---

## 第二部分：API 工程质量路线（V30-Q01…V30-Q36）

> 2026-11-11 并入 doc_v30。**下列均未实施。**
> ID 与 [01-review.md](01-review.md) 第二部分对应。原则：先加机器门，再收口公开面，
> 最后补文档；每一步可独立验收。批次前缀 `Q-` 用于与第一部分（V30-01…09）区分。

### 1. 批次与依赖

```
Q-A0 机器门（无依赖，先做）
 ├─ Q-A1 公开面收口 ── 依赖 Q-A2（distance 先有驱动）
 ├─ Q-A2 run_distance
 ├─ Q-A3 发布策略
 ├─ Q-A4 文档垂直切片
 └─ Q-A5 请求构造/校验
Q-B1–Q-B7 一致性（A 之后或并行）
Q-C1–Q-C9 卫生（随批）
```

#### Q-A0 机器门（先做，成本最低，保护后续所有改动）

| 动作 | 覆盖 | 验收 |
|---|---|---|
| 公共 API 快照：`cargo public-api --simplified -p reed > ci/api/reed.txt`（及下层公开 crate）入库，CI diff | V30-Q01/V30-Q02/V30-Q12 | 新增/删除公共项必须显式更新快照 |
| `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --workspace` 进门 | V30-Q35 | 当前先修 3 个 crate 的 9 处断链（私有 `MAX_DEPTH`、`covariance`、`0,1`、`grm` 歧义） |
| `cargo fmt --check`、`cargo deny check` 进门 | V30-Q32/V30-Q36 | 两条命令在 CI 绿 |
| CI 条件修正：PR→main 也跑 `gates`；main push 至少跑 test+clippy | V30-Q32 | `ci.yml` 的 `build`/`gates` 不再按 main 二选一 |
| `--no-default-features` 矩阵 | V30-Q31 | `cargo check -p reed-reduce --no-default-features` 绿 |

#### Q-A1 公开面收口（V30-Q01/V30-Q02/V30-Q03/V30-Q04/V30-Q06）

1. 门面 `lib.rs` 改为**逐项**再导出；删除裸 `pub use <crate|module> as <name>;`
   （`distance/graph/reduce/bgzf/chunk/header/budget/metrics/atomic_file`）。
2. `StatOutput` 与 `QcReportAcc` 的字段类型正式导出；`bcf_dict()` 删除或 opaque 化。
3. 新增 `pub mod prelude`（请求/报告/`Error`/`ThreadBudget`/`WindowGrid`/`Mat`）。
4. `#[doc(hidden)]` 收口：`open_regions`/`ExecStrategy` 要么正式导出并文档化，要么移入
   `#[doc(hidden)] pub mod frontend`。
5. 加守卫：`! grep -E '^pub use .* as [a-z_]+;' crates/reed/src/lib.rs`；公共签名类型可达性测试。

**依赖**：`distance/graph/reduce` 三个别名要等 Q-A2 完成后再删；`bgzf/chunk/header/budget`
可先逐项替换（CLI 只需 1–4 个符号）。

#### Q-A2 `run_distance`（V30-Q16）

- 新增 `DistanceRequest`（kind/missing/tree/bootstrap/prefix/…）、`DistanceReport`（矩阵 +
  bootstrap 支持度 + 写出路径）、`run_distance(input, &request, &budget)` 与 `_opened`。
- 把 `reed-cli/src/translate/distance.rs` 的矩阵选择、bootstrap、写盘、`GraphError` 分类
  搬进门面；前端只构造请求。
- 验收：CLI 输出与当前逐字节一致（gold 对拍）；门面上 `reed::distance` 别名可删除。

#### Q-A3 发布策略（V30-Q29/V30-Q30）

必须先由用户裁决"发布什么"：

- **选项 1（推荐）**：发布门面 `reed` + 其依赖闭包（core/io/engine/stats/filter/distance/
  reduce/graph）；全部 `version.workspace`、`edition 2024`、`[lints] workspace = true`、
  补 `repository/readme/keywords/categories/docs.rs`、crate 级 LICENSE；`reed-cli`/`oracle`
  保持 `publish=false`。
- **选项 2**：明确"仅 git 依赖"，则删除 distance/reduce/graph 的可发布暗示（补
  `publish=false`），README 写清 `[dependencies] reed = { git = ... }`。

无论哪项：`cargo package --dry-run` 对每个发布目标必须成功（当前 distance/reduce 因
path 依赖无 `version` 必失败）。

#### Q-A4 文档垂直切片（V30-Q24/V30-Q25/V30-Q26）

先做一条完整链路，再复制：

1. `reed` crate 级示例：构造 `PopgenRequest` → `run_popgen` → 消费 `windows`（可运行 doctest）。
2. `FilterSpec` 与 `PcaRequest` 各一个 `# Examples`。
3. README 增 "Library use"（`cargo add` 或 git 依赖、`run_*` 三段示例、feature flags、MSRV）。
4. `docsite/.../library.md` 由代码生成或至少与 `lib.rs` 的 `pub use` 清单对齐（当前写"六个
   入口"，实际 ≥8；称下层内部却整树导出）。
5. 清理公开 rustdoc 的 `doc_vNN` 票号（保留一节 `# Implementation notes` 链接内部文档）。

#### Q-A5 请求构造与校验（V30-Q07/V30-Q08/V30-Q09）

- 删除 `PcaRequest/DapcRequest.threads/memory`（或让 `run_*` 真读，并用 `Threads`/`MemMiB`）。
- 请求字段私有化 + `new()`/builder + `validate()`；非法组合（maf 区间、policy/alt 数、
  `Levels` 全 false）在构造期拒绝。
- `PopgenRequest::default()` 不再静默 0 行：要么返回 `Usage` 错，要么有明确的
  `PopgenRequest::empty()` 语义。
- `String` 枚举改 typed：`Stand/Missing/Solver/XvalResult`。

#### Q-B 一致性（P1）

| 批次 | 项 | 验收 |
|---|---|---|
| Q-B1 | V30-Q17 统一 `run_*` 签名 + 补 `run_filter_opened` + 合并 `open`/`open_regions` | 14 个入口签名表一致；API 快照更新 |
| Q-B2 | V30-Q18 `run_stats` 改数据返回，写出显式 `write_*` | `.afreq` 等可由返回的 `Vec<Row>` 直接消费；gold 输出逐字节不变 |
| Q-B3 | V30-Q11/V30-Q10 错误分类与 `FromStr` | `ErrorKind`/`FilterError` non_exhaustive；`Error::locus()`；核心类型 `FromStr<Err=reed::Error>` |
| Q-B4 | V30-Q05/V30-Q20 唯一 `Region` + `RegionSource`；`Length` 下沉 | 请求中不再有两个互斥字段；CLI `parse_size` 删除 |
| Q-B5 | V30-Q19 便捷构造器 | `StatKind::pi(g)`、`PopgenRequest::pi_scan(...)`；`window_kinds` 改 `&[WindowStat]` |
| Q-B6 | V30-Q33 公开树是否带 tests/docs（产品决策） | 若保留剥离，README 明说"完整开发树在 dev"；若改，`publish.sh` 纳入 docs/tests |
| Q-B7 | V30-Q31 修复 `reed-reduce` fallback | `cargo check --no-default-features` 绿；CI 矩阵 |

#### Q-C 卫生（P2）

V30-Q06 prelude（并入 Q-A1）、V30-Q12 `#[non_exhaustive]`、V30-Q13 `Mat` 封装、
V30-Q14 `FilterOptions` 清理、V30-Q15 `Populations` 统一、V30-Q21 `to_args` 统一、
V30-Q22 前端关注点移出请求、V30-Q23 表达式存 AST、V30-Q27 oracle 文档、
V30-Q28 lint 提级（`missing_docs=deny` 对公开 crate）、V30-Q34 删除
`crates/*/.gitignore` 的 `tests/` 并给 `docs/rounds/` 加否定规则（或明确 `add -f` 流程）、
V30-Q36 补 crate LICENSE / 依赖策略。

### 2. 验收门（第二部分，全部批次完成后）

1. `cargo public-api` 快照入库，且 `reed` 根公共项 ≤ 目标清单（建议 ≤80 个命名项，无模块别名）。
2. `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --workspace` 绿。
3. 每个 lib crate ≥1 个可运行 doctest；`reed` ≥4 个；README 示例在编译门内。
4. `cargo test --workspace` + `cargo test --doc --workspace` 绿；PR→main 也跑。
5. `cargo check -p reed-reduce --no-default-features` 绿。
6. 每个发布目标 `cargo package --dry-run` 绿；元数据表完整。
7. `cargo fmt --check`、`cargo deny check` 绿。
8. 公共签名中的每个类型都能从 `reed::*` 或 `reed::prelude::*` 命名。

### 3. 明确不做（第二部分）

- 不改任何统计公式/输出格式（v30 已审；本轮只动 API 形状）。
- 不在未裁决发布策略前改 `publish`/版本号。
- 不为了"更 Rust"重写算法模块；只改边界、命名、构造与文档。
- 不把内部测试（oracle/gold）改成公共 API；它们可以继续用下层 crate。


---

## 第三部分：根因修复（V30-R01…V30-R14）

> 2026-11-11 追加，来自四项根因审查（见 [01-review.md](01-review.md) §14）。
> **均未实施。** 原则：先修"规则归属"（R-A），再修"能力与资源所有权"（R-B），
> 最后用机器门把归属固定下来（R-C）；传播类（R-D）依赖 R-A。

### R-A 语义归属（先做：R01/R02/R03/R04）

| ID | 动作 | 验收门 |
|---|---|---|
| **V30-R01** | `reed-core`：`GtCode::has_called_allele()` + `AltCount` newtype；convert 的 `present` 改用它；消费方的 5 份 ALT 计数收敛为 1 | `--phylip/--012/--msa` 与 `stats --missing` 在半缺失上一致；`rg -n "fn count_alts\|fn n_alt_of" crates` 命中 1 |
| **V30-R02** | `reed-core::count::InfoCounts` 唯一解析 INFO 的 AC/AN（长度域 + checked 求和）；`stats` 明确 INFO/GT 单一来源 | 溢出夹具在 filter 与 stats 返回同一 typed 错误；`.afreq` 与 `.stats` 的 AF 同源 |
| **V30-R03** | `reed-io::bcf::record::pos1()` 唯一处理哨兵；内部统一 `reed_core::coord` | `rg -n "try_from\(header\.pos0\)\|pos0 == -1" crates` 命中 1；四命令哨兵行为一致 |
| **V30-R04** | `Region` 移到 `reed-io` 并持有 `coord::Interval1`；解析统一 `parse_in(text, contigs)`；convert 传 `reference.keys()` | `rg -n "fn contains\|fn region_allows" crates` 命中 1；`convert --msa --region` 用例通过 |

### R-B 能力与资源所有权（R05–R11）

| ID | 动作 | 验收门 |
|---|---|---|
| **V30-R05** | input 层拆 `StreamedInput`/`ReplayableInput`；策略选择接收 replay 能力 | stdin + `--strategy two_pass_reparse` 报 usage；`-r` 时输入只 open 1 次 |
| **V30-R06** | `IoError::Capability`，映射 Usage；`open_regions` 构造期校验 `Format` | plain/无索引/stdin 三条命令 exit code 与文案入 `cli_contract` |
| **V30-R07** | `InputMeta { format, framing, bcf_dict }` 唯一带出；删 `source_mut` | `bcf\.is_some\(\)` 与 `source_mut\(\)` 收敛；BCF 夹具断言 `Format::Bcf` |
| **V30-R08** | `PrunePool`（Drop 先取消再 join）；`prune_one` 接取消信号 | 故障注入后线程计数归零 |
| **V30-R09** | `ThreadBudget`/`MemBudget` 单向注入；convert spill 用预算（8 MiB 为下限） | `--memory 64 --mind` 的 strategy 反映 64 MiB；`default_half_memory` 收敛为构造处 |
| **V30-R10** | `TempSpill` 拥有 `File`/inode；注册进 `AtomicFile` 的存活表与启动扫除 | symlink 同址测试不被改写；孤儿 spill 被扫除 |
| **V30-R11** | `OrderedOutput::create` / `from_parts` 二分；`adopting` 加断言 | `create(p).adopting(a)` 后写入可见于 `p` |

### R-C 机器门（把归属固定下来，与 A0 同批）

| 门 | 覆盖 |
|---|---|
| `rg` 计数门：ALT 计数、`pos0` 转换、`region` 判定、`EditPlan {`、`source_mut`、`bcf.is_some()` 各收敛为 1/0 | R01/R03/R04/R07/R14 |
| 标签矩阵门：5 个 `--fill-tags` 标签 × 6 个枚举点全绿，且 `cli.md` 集合 == 表 | R13 |
| 架构门：`crates/reed-filter/src` 内 `reed_io::` 出现 0 次（或显式豁免清单） | R12 |
| 故障注入测试：prune 线程归零、spill 孤儿扫除、`create(p).adopting(a)` 写入可见 | R08/R10/R11 |

### R-D 传播与结构（依赖 R-A）

| ID | 动作 | 验收门 |
|---|---|---|
| **V30-R12** | 语义收回 L0/L1；`reed-filter` 去掉 `reed-io` 依赖（或改成 L1 trait） | `cargo tree -p reed-filter` 不含 `reed-io`；架构门绿 |
| **V30-R13** | `RecomputeTag` 静态表；filter/io/cli 皆由表生成 | 加一个标签只改表 + 测试；矩阵门绿 |
| **V30-R14** | `SurvivorEdit` 纯值；两个 job 只消费 | `EditPlan {` 组装点 = 1；三路径输出逐字节一致 |

### 第三部分明确不做

- 不因为"更分层"重写 IO/统计算法；只动规则归属、能力接口与资源所有权。
- 不在 R12 拍板前强拆 `reed-filter` 的依赖（先给门与豁免清单，再拆）。
- 不把 `#[cfg(test)]` 的故障注入工具变成产品代码路径。


---

## 第四部分：工具化全库审计的修复路线（V30-T01…T20）

> 2026-11-11 追加，对应 [01-review.md](01-review.md) 第四部分。**均未实施。**
> 原则：先把"门"修对（S-A/S-E），再修契约与结构（S-B/S-C），文档随后（S-D）。

### S-A API 与 semver（T01/T02/T03）

| ID | 动作 | 验收门 |
|---|---|---|
| T01 | CI 用逐 crate `cargo semver-checks -p <crate> --baseline-rev <tag>`；发布点打 tag | 日志 `Skipping` = 0；每个公开 crate 有 202 checks |
| T02 | 导出 `StatOutput`/`EngineError`（或统一 `reed::Error`）；`bcf_dict` 删除/opaque；补全 `QcReportAcc` 字段类型 | 公共签名类型可达性编译门；12 个旁路测试 → 0 |
| T03 | 逐项再导出，删除 13 条整 crate/模块别名（先补 `run_distance`）；`ci/api/*.txt` 入库 + CI diff | `grep '^pub use .* as [a-z_]+;'` = 0 |

### S-B 发布与外部契约（T04/T05/T06/T08）

| ID | 动作 | 验收门 |
|---|---|---|
| T04 | 裁决发布策略；补 repository/documentation/keywords/categories/README/LICENSE；修 reduce 的 version | 3 个 crate `cargo package --dry-run` 绿 |
| T05 | `FileIndex`/`Indexer` opaque 化；`ChunkPayload` map 私有；serde 放 feature 或 allow-list | external-types allow-list 只剩批准项；门面探针可命名所有类型 |
| T06 | CI 固定 `nightly-2026-03-20`（或升级工具） | 9 份 external-types 表入库且 exit 0 |
| T08 | 删除 `bstr`、`bitflags`(io)、`tracing`(engine)、`sha2`(filter)、`thiserror`(graph/stats)；`sha2`(reed) 移 dev-deps | `cargo machete` 0 命中 |

### S-C 结构与分层（T09/T10/T11）

| ID | 动作 | 验收门 |
|---|---|---|
| T09 | filter 的最小传输视图下沉 `reed-core`；删 `reed-io` 依赖；断 2 处模块环 | `cargo tree -p reed-filter` 不含 reed-io；模块 SCC = 0 |
| T10 | architecture.md 补齐 distance/reduce/graph/oracle 的层级与 reed-qc 去向 | 11/11 成员都有 Lx 归属 |
| T11 | `convert/job.rs` 按格式拆分；28 个 >800 行文件逐步治理 | `job.rs` < 800 行；结构守卫加文件上限 |

### S-D 文档（T12/T13/T14/T15）

| ID | 动作 | 验收门 |
|---|---|---|
| T12 | 修 19 处 rustdoc 错误；消除 lib/bin 同名；CI 加 doc 门 | `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --workspace` exit 0 |
| T13 | distance/reduce 继承 lints；公开 crate `missing_docs=deny`；每 crate ≥1 示例 | `-W missing_docs` = 0；doctest ≥1/crate |
| T14 | 文档入口表由导出生成；README 库示例；examples 移到 bench/ 或新建门面示例 | 入口数 == 14；README 示例可编译 |
| T15 | 票号移入 Implementation notes；加 `# Stability` 并统一三处口径 | 公开 rustdoc `doc_vNN` = 0；三处一致 |

### S-E 门与运行时（T07/T16/T17/T18/T19/T20）

| ID | 动作 | 验收门 |
|---|---|---|
| T07 | 修 `lanczos.rs:146` 的 `parallel`；CI 加 `cargo hack --each-feature` 与 `--no-default-features` | 两条命令 exit 0 |
| T16 | main/PR→main 也跑 test/clippy/doc；补 fmt/deny/structure；`--doc --workspace` | CI 的 main 分支有 test job |
| T17 | distance/reduce 补 lints 并修 28 条；oracle 显式 deny unsafe；expect invariant 守卫 | `cargo clippy -p reed-reduce --lib -- -D warnings` exit 0 |
| T18 | 修 fallback 后清理 inactive_code；其余分支删除或覆盖 | faer 相关 inactive = 0 |
| T19 | `prune_one`/`prune_matrices` 都接 `AbortSignal`，`PrunePool::drop` join | 故障注入后线程归零 |
| T20 | `AtomicFile` 持有 `File`；删除/弃用 `reed-distance` 路径式 writer | 按名重开为 0；`File::create` 为 0 |

### 第四部分验收门（全部完成）

1. CI 的 main 与 PR→main 都跑 test/clippy/doctest/doc/fmt/deny/structure。
2. 每个公开 crate 在 semver-checks 中有结果（无 `Skipping`）。
3. `cargo public-api` 快照入库，门面无整 crate 别名。
4. `cargo hack check --workspace --each-feature` 与 `--no-default-features` 绿。
5. `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --workspace` 绿。
6. `cargo test --doc --workspace` 每 crate ≥1。
7. 3 个可发布 crate 的 `cargo package --dry-run` 绿。
8. external-types 9 份表绿（固定工具链）。


---

## 第五部分：综合方案（P0–P4，替代逐条修复顺序）

> 2026-11-11 追加。第一~四部分是"证据与条目"，**不是工作分解**。综合诊断见
> [03-synthesis.md](03-synthesis.md)：79 条 → **5 个根因 + 1 个产品裁决**。
> 实施顺序按本部分；旧条目只作为验收清单。

| 阶段 | 内容 | 关闭的根因 | 依赖 |
|---|---|---|---|
| **P0** | 契约裁决 ADR（发布面/稳定面/外部依赖/MSRV/兼容承诺） | F（4 条） | — |
| **P1** | 语义内核：GT/ALT/AC-AN/坐标/区域/多等位/SiteFields 收进 `reed-core`；filter 去 `reed-io` | B（15）+ 部分 A | P0 |
| **P2** | 能力与状态类型：`Streamed/Replayable/SeekableInput`、`InputMeta`、`RunContext`、`PrunePool`、`TempFile`、`OutputTarget`、预算注入 | C（18） | P1 |
| **P3** | 门面收口 + 六工具门：逐项导出/prelude/可命名性、public-api 快照、semver、external-types、hack、package、CI | A（17）+ F（4） | P1/P2 |
| **P4** | 派生文档与结构：doctest/入口表/`# Stability`、`RecomputeTag`、`SurvivorEdit`、拆分 `convert/job.rs`、architecture 更新 | D（16）+ E（9） | P3 |

**验收**：`03-synthesis.md` §6 的 12 条工具门。每条旧 issue 只有在对应门命令通过后才关闭；
未通过的继续以测试用例存在。

**明确不做**：不逐条修 79 条；不新造"第 80 条"；不在 P0 前发布/打 tag。
