# STATUS：doc_v30 审计状态与复现记录

 **2026-11-11 立项。本轮是审计轮：未改任何代码。** 本文件 §1–§8 的表格和基线是 2026-11-11 的历史记录，不能单独当作当前状态；后续实施记录见 §10–§13，当前权威收口见 [06-final-audit.md](06-final-audit.md)。第二部分
（API 工程质量）的基线与实测记录见 §5。

## 0. 基线

| 项 | 值 |
|---|---|
| 工作树 | HEAD `e920663`（2026-11-10）+ 49 个文件未提交改动（含 `reed-io`、`reed-filter`、`reed/src/stats/writers/reports.rs`、`reed/src/filter/mod.rs`） |
| sgkit 侧 | `/sri/home/sri2025201067/reed/sgkit/docs/api.rst` + `sgkit/**`（当前工作树） |
| 复现二进制 | `reed 0.4.0 (v0.3.0-177-ge920663-dirty 2026-11-10)`，位于 `/tmp/reed-target/release/reed` |
| 注意 | `.cargo/config.toml` 把 `target-dir` 重定向到 `/tmp/reed-target`；仓库内 `target/release/reed` 是 **11-03 的旧二进制**（`v0.3.0-117-g02109e9-dirty`），第一轮实测误用它，已用当前源码重建并复跑 |

## 1. 条目状态

| ID | 优先级 | 状态 |
|---|---|---|
| V30-01 PI_HAT 期望 IBS 系数错位 | P0 | **未实施**（已复现） |
| V30-02 `filter -I` 文档与行为相反 | P1 | **未实施**（已复现） |
| V30-03 多等位三路径不一致 | P1 | **未实施**（代码层证明） |
| V30-04 f-stat `estimate=Σ` 未文档化 | P2 | **未实施** |
| V30-05 `--f4-ratio` 角色序 | P2 | **未实施** |
| V30-06 `f3*` 承诺不可达 | P2 | **未实施** |
| V30-07 MSA 全缺失→REF 未警示 | P2 | **未实施** |
| V30-08 `identity_by_state` 命名 | P2 | **未实施** |
| V30-09 新增 sgkit 数值夹具 | P2 | **未实施** |

## 2. 复现记录

### 2.1 V30-01 PI_HAT（`repro`）

```
$ /tmp/reed-target/release/reed distance /sri/home/sri2025201067/reed/subvcf.vcf.gz --type pihat -o pih2
reed distance: 348 samples, 3000 sites, type=pihat, format=triangle

$ awk -F'\t' 'NR>1{print $7" "$8" "$9" "$10}' pih2.genome | sort | uniq -c | sort -rn | head -1
  60346 0.0000 1.0000 0.0000 0.5000
$ awk 'NR>1' pih2.genome | wc -l
  60378
```

60378 对里 **60346 对**输出 `Z0=0 Z1=1 Z2=0 PI_HAT=0.5`；其余 32 对为近重复样本。
代码：`reed-distance/src/matrices.rs:570-592`（系数）、`:519-521`（求解器）。
plink 侧对照：`plink-ng/1.9/plink_calc.c:4863-4867`。

### 2.2 V30-02 `filter -I`（`repro`）

```
$ R filter subvcf.vcf.gz -x '' --fill-tags AC,AN -Ov -o c_withAC.vcf     # 造出 AC=331 的输入
$ R filter c_withAC.vcf -s FENGGWS193        -Ov -o a.vcf   # AC_rows=3000  first=AC=1
$ R filter c_withAC.vcf -s FENGGWS193 -I     -Ov -o b.vcf   # AC_rows=0
$ R filter c_withAC.vcf -s FENGGWS193 -I -x '' -Ov -o c.vcf # AC_rows=3000  first=AC=331
$ R filter c_withAC.vcf -I                   -Ov -o d.vcf   # AC_rows=0   samples=348
```

代码：`reed-filter/src/spec/recompute.rs:40-54`（默认陈旧清单含 AC/AN）、
`reed/src/filter/mod.rs:520-538`（`-I` 只关重算）、`reed-filter/src/plan/compile/finish.rs:81`。
文档承诺：`docsite/src/content/docs/filter.md:155-156`（"Leave INFO/AC and INFO/AN alone"，
**与实际不符**）、`reed-cli/src/cli/mod.rs:365-366`（"Do not update …"，措辞更弱但不完整）。

### 2.3 V30-03 多等位（`code`，未构造端到端夹具）

- `stats` 折叠：`reed/src/stats/run.rs:298-304`（`MultiallelicPolicy::CollapseAlt`）；
  `reports.rs:1383` 单一 `ALT_FREQS` 列。
- `popg` 跳过：`reed/src/stats/run.rs:452`。
- `filter` 拒绝：`reed/src/filter/mod.rs:543-554`。

**未做的验证**：没有用多等位 VCF 实测三命令输出的差异（只做了代码路径核对）。

## 3. "不偏移"核对

### 3.1 文档 ↔ 代码

写每个引用的 file:line 时都用当前工作树重读；`git status` 显示相关文件（`filter/mod.rs`、
`stats/run.rs`、`reports.rs`、`translate/distance.rs`）处于未提交状态，**核对的是未提交后的内容**，
不是 HEAD 旧版。

### 3.2 守卫脚本

| 命令 | 结果 |
|---|---|
| `bash ci/check-doc-drift.sh` | **exit 0**（全部 ok；本轮为纯文档新增，未触发任何既有守卫） |
| `bash ci/check-structure.sh` | **exit 0**（16/16 ok） |

### 3.3 结论的分级（诚实标注）

| 证据级 | 条目 | 说明 |
|---|---|---|
| `repro`（当前二进制实测） | V30-01、V30-02 | 命令与输出见 §2 |
| `code`（确定性代码路径） | V30-03、以及 §1 表中所有"未实现"行（全仓 `rg` 0 命中或驱动显式拒绝） | 未跑端到端夹具 |
| `formula`（公式层） | §2.1 全部"站得住"、§2.2 全部口径变体 | **未跑 sgkit 数值夹具**；对拍对象是 pixy/allel/plink2/vcftools/ADMIXTOOLS |
| 主代理亲自重读 | PI_HAT 系数与求解器、`-I` 全链路、`stats_policies`/`popg` 多等位门、`ALT_FREQS` 头、`FstKind`/`PbsKind`、10 个 `--type`、`reed-reduce` 模块表、`reed-graph` tree/bootstrap、`window_by_variant` 拒绝、convert 目标、index CSI/TBI、IBS/GRM 主体 | — |
| 子代理核、主代理抽查未逐行 | `reports.rs:843-909` 细节列、`ld.rs:452-590` 剪枝、`matrix.rs:284-357`、`summary/report.rs:165-192`、`pca.rs:58-77` 等 | 这些行号来自子代理，主代理只核了其中一部分；实施修复前应再核一次 |

### 3.4 已知不精确处

1. §1 的 56 项覆盖计数由 6 份子代理表汇总而来；主代理未逐行重读全部 56 行的实现。
   实施任何"补齐缺口"的工作前，需按 01 的 file:line 再核一次对应行。
2. `docs/rounds/v17/STATUS.md` 把 PI_HAT 差异归因于"Nei 偏差修正"；本轮的结论是
   **系数排列本身错**（§2.3），该归因应视为已修正（v17 文件本身不改，由本轮沿用）。
3. 本轮未复核 v27 对多等位的拍板原文，只引用其结论。

## 4. 产物

| 文件 | 说明 |
|---|---|
| [README.md](README.md) | 轮次总结 |
| [01-review.md](01-review.md) | §1–§5：56 项覆盖表、科学符合性、真缺陷复现、实现但 sgkit 无 API 的能力清单；§6–§13：API 工程质量 36 条 |
| [02-roadmap.md](02-roadmap.md) | V30-01…09 优先级与验收门 |
| `docs/README.md` | 索引新增 v30 一行 |


## 5. 第二部分：API 工程质量审计（V30-Q01…V30-Q36）

> 2026-11-11 并入 doc_v30。标准与证据分级同 01 §6；
> 36 条全部为"未实施"，修复计划见 [02-roadmap.md](02-roadmap.md) 第二部分。

### 5.1 基线（第二部分）

| 项 | 值 |
|---|---|
| 工作树 | 分支 `dev`，HEAD `f190708`（2026-11-11），含 60+ 未提交改动 |
| 公开树 | `origin/main`（`ci/publish.sh` 生成）259 文件；`origin/dev` 579 文件 |
| 工具链 | Rust 1.91；workspace edition 2024（distance/reduce 仍 2021） |
| 构建目录 | `.cargo/config.toml`（已被 `.gitignore:6` 忽略、从未提交）把 target-dir 指向 `/tmp/reed-target` |

### 5.2 条目状态

| 批次 | ID | 状态 |
|---|---|---|
| Q-A0 机器门 | — | **未实施** |
| Q-A1 公开面收口 | V30-Q01/V30-Q02/V30-Q03/V30-Q04/V30-Q06 | **未实施** |
| Q-A2 `run_distance` | V30-Q16 | **未实施** |
| Q-A3 发布策略 | V30-Q29/V30-Q30 | **未实施**（需用户裁决） |
| Q-A4 文档 | V30-Q24/V30-Q25/V30-Q26 | **未实施** |
| Q-A5 请求构造 | V30-Q07/V30-Q08/V30-Q09 | **未实施** |
| Q-B 一致性 | V30-Q05/V30-Q10/V30-Q11/V30-Q17/V30-Q18/V30-Q19/V30-Q20/V30-Q31/V30-Q32/V30-Q33 | **未实施** |
| Q-C 卫生 | V30-Q12/V30-Q13/V30-Q14/V30-Q15/V30-Q21/V30-Q22/V30-Q23/V30-Q27/V30-Q28/V30-Q34/V30-Q35/V30-Q36 | **未实施** |

### 5.3 实测记录（`repro`）

#### 5.3.1 公开树被剥离 tests/docs（V30-Q33）

```
$ git ls-tree -r --name-only origin/main | wc -l
259
$ git ls-tree -r --name-only origin/main | grep -c '/tests/'
0
$ git ls-tree -r --name-only origin/dev | wc -l
579
$ git ls-tree -r --name-only origin/dev | grep -c '/tests/'
65
$ git ls-tree -r --name-only origin/dev | grep -c '^docs/'
151
```

机制：`ci/publish.sh:22-43` 的 `PUB_PATHS` 不含 `docs/`/`bench/`，并显式
`':!crates/*/tests/'`；脚本内还有泄漏检查。`origin/main` 无 `CHANGELOG.md`。

#### 5.3.2 `--no-default-features` 编译失败（V30-Q31）

```
$ cargo check -p reed-reduce --no-default-features
error[E0063]: missing field `parallel` in initializer of `Operator<'a>`
  --> crates/reed-reduce/src/lanczos.rs:146:13
  146 |             Self { w }
warning: unused variable: `par`  --> crates/reed-reduce/src/backend.rs:155
warning: unused variable: `par`  --> crates/reed-reduce/src/backend.rs:206
```

#### 5.3.3 rustdoc 失败（V30-Q35）

```
$ RUSTDOCFLAGS='-D warnings' cargo doc --no-deps -p reed-distance -p reed-graph -p reed-io
error: public documentation for `to_newick` links to private item `MAX_DEPTH`
error: unresolved link to `covariance`
error: unresolved link to `0,1`
error: `grm` is both a function and a module
exit=101
```

#### 5.3.4 新测试被静默忽略（V30-Q34）

```
$ git check-ignore -v crates/reed/tests/new_contract.rs
crates/reed/.gitignore:2:tests/*  crates/reed/tests/new_contract.rs
```

7 个 crate 的 `.gitignore` 为 `tests/`；已跟踪的 65 个测试文件是历史遗留，
新增文件默认不会进 `git add`。

同一条规则也命中 `docs/`（根 `.gitignore:19`）：v30 文档在 `git status` 中不可见：

```
$ git check-ignore -v docs/rounds/v30/README.md
.gitignore:19:/docs/  docs/rounds/v30/README.md
$ git ls-files docs/rounds/v30 | wc -l
0
```

本轮已用 `git add -f` 显式索引 v30 的 4 个文件（未提交），以免记录丢失。

#### 5.3.5 死字段 / 死导出 / 缺失驱动（V30-Q01/V30-Q02/V30-Q07/V30-Q16）

```
$ rg -n "reed::reduce" crates            # 0 命中（整 crate 别名零使用）
$ rg -n "pub fn run_distance|pub struct DistanceRequest" crates   # 0 命中
$ rg -n "request\.(threads|memory)" crates/reed/src/pca crates/reed/src/dapc   # 0 命中
$ rg -n "StatOutput" crates/reed/src/lib.rs      # 0 命中（只在 request.rs:264 作为字段类型）
$ grep -c "^pub use" crates/reed/src/lib.rs      # 30
```

#### 5.3.6 其他静态计数

| 指标 | 值 | 命令 |
|---|---|---|
| `impl FromStr` | 2 | `rg -n "impl (std::str::)?FromStr" crates \| wc -l` |
| `#[doc(hidden)]` | 4 | `rg -n "doc\(hidden\)" crates/*/src \| wc -l` |
| `# Examples` | 2 | `rg -c "# Examples" crates/*/src` |
| `pub mod prelude` | 0 | `rg -n "pub mod prelude" crates` |
| `open_regions` 前端调用 | 6 | `rg -c open_regions crates/reed-cli/src` |
| 测试文件绕过门面 | 13 | `rg -c "reed_(core\|io\|engine\|stats\|filter)" crates/reed/tests/*.rs` |
| CI 调用 `check-structure.sh` | 0（仅脚本自身注释） | `rg -n check-structure .github ci` |

> 本表是 2026-11-11 审计时的静态计数，不是当前状态；`pub mod prelude` 已在 P3 补齐（见 §10）。

#### 5.3.7 发布标记

`publish = false`：oracle、reed-cli、reed-core、reed-engine、reed-filter、reed-io、
reed-stats（8 个，审计时的状态）。`reed` 已在 P3 改为发布目标；`reed-distance`/
`reed-reduce` 已 bump 到 0.5.0 并补 `version`，`reed-graph` 保持 `publish = false`。
`ci/check-package.sh` 现在验证清单与依赖版本（完整 tarball 仍待内部 crate 上架）。

### 5.4 "不偏移"核对

#### 5.4.1 主代理亲自复核（写进 01 的依据）

| 结论 | 复核方式 |
|---|---|
| 门面整树别名 / `reed::reduce` 零引用 | `sed` 读 `lib.rs:138-176`；`rg` 计数 |
| `StatOutput` 不可命名 | `rg StatOutput` 仅命中 `request.rs:264`；`lib.rs` 无导出 |
| 测试绕过门面 | 13 个文件命中，最高 `e2e_filter.rs` 79 处 |
| `PcaRequest` 死字段 | `rg request.(threads\|memory)` 零命中；读 `Default` 确认 `prefix="out"` |
| `PopgenRequest::default()` 静默 0 行 | 读 `driver.rs:204-209`：`kinds.is_empty()` → `GridMode::None`，无报错 |
| `run_index`/`run_filter`/`run_rename` 签名 | 读 `index.rs:524`、`filter/mod.rs:244`、`rename/mod.rs:374` |
| `StatsReportRequest.prefix` / `StatsOutcome.files` | 读 `reports.rs:39-110,194-207` |
| 区域三份 | 读三处定义；`rg` 两条 `Region` |
| `ErrorKind`/`FilterError` 非 non_exhaustive | 读 `error.rs:40-46`、`reed-filter/src/error.rs:8-16` |
| `Mat` 字段 pub | 读 `reed-core/src/linalg/mod.rs:7-20` |
| `to_args` 四种签名 | `rg "pub fn to_args"` 四条 |
| `window_kinds` 7 参数 | 读 `stats/request.rs:94-112` |
| `FilterOptions::strategy` 隐藏字段 | 读 `filter/mod.rs:60-72` |
| oracle 文档首行 | `head -3 crates/oracle/src/lib.rs` |
| docsite "六个入口" | `grep -n six docsite/.../library.md` |
| CI main/PR 只 build | 读 `.github/workflows/ci.yml:12-26` |
| rustdoc 失败 | 实跑，见 §2.3 |
| `--no-default-features` 失败 | 实跑，见 §2.2 |
| 公开树无 tests/docs | 实跑 `git ls-tree`，见 §2.1 |
| `.cargo/config.toml` 未入库 | `git check-ignore` 命中、`git ls-files .cargo/` 为空 |

#### 5.4.2 已撤回的上轮担忧（避免下轮重复排查）

- **`.cargo/config.toml` 不是"个人环境配置入库"**：它被 `.gitignore:6` 忽略、`git ls-files`
  为空、无提交历史。属于本地开发配置，**不是缺陷**。
- **`crates/*/.gitignore` 的 `tests/` 不影响已跟踪文件**：65 个测试文件已在索引中；
  缺陷只在于**新增**测试会被静默忽略（V30-Q34 已按此措辞记录）。

#### 5.4.3 未核实（写进 01 时已标注）

1. 未跑 `cargo package --dry-run`（V30-Q29 的"打包必失败"由 manifest 的 path 无 `version`
   推出，标 `code`）。P3 已用 `ci/check-package.sh` 固化该检查（`--list` + 依赖版本）。
2. 未跑全量 `cargo test --workspace`（耗时；本轮只做 API 静态审计）。
3. 子代理结论、主代理未逐行复核：`QcReportAcc` 的完整字段树、`FilterSpec` 的
   policy/alt 矛盾的全部组合、`deny.toml` 的 license 白名单、`run_stats` 的
   `compact_window_accumulators` 分支。
4. 未裁决：公开树是否应带 tests/docs、是否发布 crates.io——属产品决策，见 02 A3/B6。

### 5.5 产物

| 文件 | 说明 |
|---|---|
| [README.md](README.md) | 轮次总结（8 性质、36 条、6 条 P0） |
| [01-review.md](01-review.md) | 逐条发现（位置/现状/原因/建议/证据级） |
| [02-roadmap.md](02-roadmap.md) | A/B/C 批次、依赖、验收门、非目标 |
| `docs/README.md` | 索引 v30 一行（覆盖 sgkit + API 工程质量两部分） |


## 6. 第三部分：根因审查记录（V30-R01…V30-R14）

> 2026-11-11 追加。四项根因审查（语义归属 / 能力接口 / 生命周期与所有权 / 变更传播），
> 共 **14 条：4 P0（R01/R02/R03/R12）/ 10 P1**，审查时全部未实施；当前实施状态见
> §6.2 与 §10，计划见 [02-roadmap.md](02-roadmap.md) 第三部分。

### 6.1 基线（第三部分）

| 项 | 值 |
|---|---|
| 工作树 | 分支 `dev`，HEAD `f190708`（含 60+ 未提交改动） |
| 复现二进制 | `/tmp/reed-target/release/reed`，`reed 0.4.0 (v0.3.0-180-gf190708-dirty 2026-11-11)`；overflows 用 `/tmp/reed-target/debug/reed` |
| 代码检索 | `rg` + `codegraph`（索引 274 文件 / 6,187 节点 / 21,652 边，Rust 为 preview 级） |
| **未安装** | `rust-analyzer`（仅有 rustup 代理，组件缺失）、`cargo-modules`、`cargo-public-api`、`cargo-semver-checks`、`cargo-check-external-types`、`cargo-hack`、`cargo-fuzz`、`cargo-deny` |
| 共变统计 | `git log --no-merges --name-only`：触碰 `crates/` 的提交 211 个，fix/feat/refactor 162 个，其中 132 个改 ≥3 个 `.rs`；**未排除批量改名/格式化** |

### 6.2 条目状态

| 批次 | ID | 状态 |
|---|---|---|
| R-A 语义归属 | V30-R01/R02/R03/R04 | **已实施**（04 §2） |
| R-B 能力与资源 | V30-R05/R06/R07/R08/R09/R10/R11 | **已实施**（04 §3） |
| R-C 机器门 | 计数门 / 标签矩阵门 / 架构门 / 故障注入 | **部分实施**：计数门、标签矩阵门（R13 帮助/注册表断言）、架构门已入 `check-structure.sh`；故障注入未做 |
| R-D 传播与结构 | V30-R12/R13/R14 | **已实施**（04 §2/§3） |

### 6.3 主代理实测记录（`repro`）

夹具目录 `/tmp/rcfee/`（8 样本 VCF）。命令与输出：

1. **半缺失判定分叉（R01）**
   ```
   $ reed convert half.vcf --phylip --min-samples-locus 4 -o hp   # 只留 POS 20 的 0/. 位点
   $ head -1 hp.phy        -> 8 1
   $ reed convert half.vcf --012 --min-samples-locus 4 -o g012
   $ cat g012.012.pos      -> chr1	20          # `./1` 位点（POS 10）被整行剔除
   $ reed stats half.vcf --missing -o hm
   $ cat hm.vmiss          -> 两位点均 MISSING_CT=4 OBS_CT=8 F_MISS=0.5
   ```
2. **AC/AN 溢出（R02）**
   ```
   $ /tmp/reed-target/debug/reed stats ovf.vcf --summary -o ovfd
   thread '<unnamed>' panicked ... attempt to add with overflow
   reed: worker panicked on chunk 0: attempt to add with overflow
   $ /tmp/reed-target/release/reed stats ovf.vcf --summary -o ovfout   # 不报错
   $ grep '^AF' ovfout.stats   -> AF 0 0.000000 ...   # 求和回绕
   ```
3. **同一命令两个 AC（R02）**：`GT` 全 `0/0`、`INFO AC=5;AN=8`：
   ```
   $ reed stats disagree.vcf --freq --summary -o dis
   $ cat dis.afreq           -> chr1  chr1_10  A  T  0  8        # ALT_FREQS=0（GT 来源）
   $ grep '^AF' dis.stats    -> AF 0 0.610000 ...                # INFO 来源
   ```
4. **two-pass 策略在 stdin 上静默单趟（R05）**
   ```
   $ cat small.vcf | reed filter - --strategy two_pass_reparse --timing t.txt -o out
   reed filter: 3 -> 3 records (4 of 4 samples), strategy two_pass_reparse, 0.0s
   $ cat t.txt   -> strategy two_pass_reparse / phase.pass1_ms 1.256   # 无 pass0
   ```
5. **能力缺失三种错误类别（R06）**
   ```
   $ reed filter small.vcf    -r chr1:10-20 -> malformed BGZF data: ...            exit 1
   $ reed filter small.vcf.gz -r chr1:10-20 -> unsupported input: regions need an index ...  exit 1
   $ cat small.vcf | reed filter - -r chr1:10-20 -> regions need a seekable file, not stdin  exit 2
   ```
6. **静态计数（R01/R04/R07/R09/R11/R14）**
   ```
   $ rg -n "starts_with\(b\"\.\"\)" crates                      -> 1（convert/job.rs:851）
   $ rg -n "regions need a seekable file" crates                -> 3（input.rs:102、filter/mod.rs:269、common.rs:378）
   $ rg -n "\.bcf\.is_some\(\)" crates/reed/src                 -> 5
   $ rg -n 'EditPlan \{' crates/reed/src                        -> 3
   $ rg -n "default_half_memory" crates/reed/src                -> 8
   $ rg -n "reed_io::" crates/reed-filter/src                   -> 6 文件 23 处
   $ rg -n "to_pos0|to_pos1|Interval1" crates                   -> 仅 coord.rs 定义处（零调用）
   ```
7. **变更传播（R12/R13/R14）**
   ```
   $ git show --stat e182ef6   -> 14 文件 / 4 crate（GT 定位）
   $ git show --stat 61273ab   -> 12 文件（GT 索引越界）
   $ git show --stat 74090da   -> 10 文件（区域解析）
   $ git show --stat d610f88   -> 6 文件（GT-less 保留 AC/AN）
   $ git show --stat 39ba997   -> 6 文件（--fill-tags 校验）
   $ git log --no-merges --name-only ...  -> 211 提交 / 162 fix|feat|refactor / 132 改 ≥3 .rs
        crate 对：reed×reed-cli 57、reed×reed-filter 47、reed×reed-io 41
        文件对：reed-filter/src/plan/mod.rs×reed/tests/e2e_filter.rs 23
   ```

### 6.4 子代理实测、主代理未复跑

| 结论 | 状态 |
|---|---|
| R03 BCF `pos0=-1` 四命令行为（pca/stats/filter/convert） | 子代理实测，主代理只复核代码路径 |
| R04 `convert --msa --region 'chr1:100'` 失败 | 子代理实测，主代理复核 `Region::parse` 与调用点 |
| R08 故障注入（prune 线程归零） | **未编写**，结论为代码推理 |
| R10 symlink 同址 / 孤儿 spill 清扫 | **未编写**，结论为代码推理 |
| R11 `create(p).adopting(a)` 写入丢失 | **未编写**，结论为代码推理（无现有调用点触发） |

### 6.5 "不偏移"核对（第三部分）

| 结论 | 复核方式 |
|---|---|
| convert `present` 唯一且为 `starts_with(b".")` | `rg` 命中 1；读 `convert/job.rs:851` |
| summary 第二份 AC/AN 无 checked 求和 | 读 `summary.rs:294-328`；debug panic 实测 |
| `pos0` 四个转换点 | 逐个读 `render.rs`/`stats/run.rs`/`record_pipeline.rs`/`genotype.rs` |
| `region_allows` vs `Region::contains` 双份判定 | 读 `region.rs:113-125,241-245`；`contains` 仅测试调用 |
| `Region`/`MsaRegion` 三类型 | `rg` 定义 + 读 `convert/mod.rs:107-116`、`job.rs:127-131` |
| 可重读能力靠路径 | 读 `input.rs:30-41`、`filter/mod.rs:130-145,345-374`；stdin 实测 |
| 三种能力错误类别 | 实跑三条命令，exit 1/1/2 |
| `.bcf.is_some()` 5 处、`source_mut` 生产在用 | `rg` 计数 + 读调用点 |
| prune 线程只有 finish join | 读 `matrix.rs:773,810-820,1027`；`rg "impl Drop for MatrixSink"` 空 |
| `TempSpill` 只有 path、append 按名重开 | 读 `spill.rs:117-152`、`convert/spill.rs:44-58` |
| `adopting` 无断言、`create` 自建 atomic | 读 `ordered_output.rs:88-95,112-126` |
| `reed-filter` 依赖 `reed-io` 违反架构 §14.6 | 读 `Cargo.toml:12` 与 `docs/spec/architecture.md:1843`；`rg` 23 处 |
| 五个提交的文件集合 | `git show --stat` 逐个核对 |
| 共变统计 | 主代理用 `git log --no-merges --name-only` 现算（口径见 §6.1） |

### 6.6 产物（第三部分）

| 文件 | 说明 |
|---|---|
| [01-review.md](01-review.md) §14–§16 | 14 条根因（证据 → 根因 → 建议归属 → 如何验证） |
| [02-roadmap.md](02-roadmap.md) 第三部分 | R-A/R-B/R-C/R-D 批次与验收门 |
| `docs/README.md` | v30 索引行补充根因审查摘要 |


## 7. 第四部分：工具化全库审计记录（V30-T01…V30-T20）

> 2026-11-11 追加。六个子代理（API 表面 / semver 发布 / 外部类型与依赖 / 模块结构 / 文档可验证性 /
> 卫生与运行时）覆盖 11 个 workspace crate；主代理复核关键数字。共 **20 条：3 P0 / 12 P1 / 5 P2**，
> **全部未实施**，计划见 [02-roadmap.md](02-roadmap.md) 第四部分。

### 7.1 工具矩阵（实测）

| 工具 | 版本 | 命令 | 结果 |
|---|---|---|---|
| rust-analyzer | 1.99.1 | `analysis-stats .` / `diagnostics .` | 109,819 LOC；`panics: 0`；620 WeakWarning（618 `inactive_code`） |
| cargo-modules | 0.28.0 | `structure/dependencies/orphans --lib -p <crate>` | 9 crate；orphans 全 0；模块级环 2 处 |
| cargo-public-api | 0.52.0 | `cargo public-api -p <crate> -s` | 9 crate 快照；reed 1120 行 / 249 声明项 / 32 `pub use` |
| cargo-semver-checks | 0.51.0 | `check-release --workspace --baseline-rev HEAD -v` | 只检查 3/11；7 个 `Skipping` + cli 无 lib |
| cargo-check-external-types | 0.5.0 | `+nightly-2026-03-20 check-external-types ...` | 9 份表；默认 nightly（JSON v61）不可用 |
| cargo-hack | 0.6.46 | `check --workspace --each-feature / --no-default-features / --all-features` | each-feature 14 plans / 1 fail；all-features 11/11 通过 |

辅助：`cargo tree -d`（宿主 7 组/14 条目；`--target all` 10 组/20 条目）、
`cargo package`（3 个"可发布" crate 全失败）、`cargo doc -D warnings`（exit 101，5 crate/19 处）、
`cargo test --doc --workspace`（6 doctest，3 执行）、`cargo clippy --workspace --all-targets`
（exit 0 / 0 warning，仅覆盖 8 个继承 lints 的 crate）。

### 7.2 raw 输出（仓库外，临时）

| 路径 | 内容 |
|---|---|
| `/tmp/audit/public-api-reed.txt`、`pa-*.txt` | public-api 快照（9 crate） |
| `/tmp/audit/semver.txt` | semver-checks `--workspace` 输出（3 crate） |
| `/tmp/audit/modules/*.structure/.deps/.orphans.txt` | cargo-modules（9 crate） |
| `/tmp/audit/hack-each.txt`、`hack-all.txt` | cargo-hack 矩阵 |
| `/tmp/audit/doctests.txt`、`rustdoc-final.txt` | doctest 与 rustdoc 门 |
| `/tmp/audit-a3/et-*.md` | external-types 9 份表（`nightly-2026-03-20`） |
| `/tmp/audit-a5/rustdoc.txt`、`missing-docs-*.txt` | 文档门与 missing_docs 强制结果 |
| `/tmp/audit-a6/clippy2.*`、`ra-stats.*`、`ra-diag.out` | clippy / RA stats / RA diagnostics |
| `/tmp/audit-a2/cargo-package.txt`、`metadata.json` | 打包失败与元数据 |

### 7.3 主代理亲自复核

| 结论 | 复核命令/方式 |
|---|---|
| semver 只检查 3/11 | `cargo semver-checks ... -v` 输出 7 条 `Skipping` |
| 门面 30 条 `pub use`、13 条别名 + 7 writers | `grep -c '^pub use' crates/reed/src/lib.rs`；逐条读 |
| `EngineError` 在 4 个公开返回值 | 读 `fstats.rs:22`、`ld.rs:38`、`sfs.rs:110/176` |
| `QcReportAcc` 已导出 | `rg QcReportAcc crates/reed/src/lib.rs` → `:134` |
| 旁路门面测试 12 个、`e2e_filter` 79 处 | `rg -l 'reed_(core|io|engine|stats|filter)::' crates/reed/tests/*.rs` |
| `reed-graph` 有 rust-version | 读 `crates/reed-graph/Cargo.toml:6` |
| distance/graph 因 registry 无 reed-core 失败 | `cargo package` 报错原文 |
| reduce 因缺 version 失败 | `cargo package` 报错原文 |
| 重复依赖 7 组/14 条目（宿主），10 组/20 条目（all-target） | `cargo tree -d [--target all]` |
| `convert/job.rs` 2105 行；28 个 >800 行 | `wc -l` + `find/xargs` |
| rustdoc 5 crate/19 处 | 实跑 `-D warnings`，按 `--> crates/...` 计数 |
| 6 个 doctest（3 执行） | `cargo test --doc --workspace` |
| unsafe 5 处操作（1 fn + 4 block） | `rg 'unsafe\s*\{|unsafe fn' crates/*/src` |
| `AtomicFile::commit` 按名重开 | 读 `crates/reed-io/src/write/atomic.rs:289` |

### 7.4 子代理实测、主代理未逐条复跑

- T05 的 9 份 external-types 表（工具链 `nightly-2026-03-20`）与 noodles/serde 项数；
- T08 的未用依赖清单（源码 `rg` 0 命中）；
- T10 的 architecture.md 对照（8 个已定义 crate vs 11 个成员）；
- T14 的 library.md 入口表与 examples/ 自述；
- T15 的票号渲染计数（1,366）与三处稳定性口径；
- T16 的 CI 条件与 gold/popgen 未接线；
- T18 的 620 条 RA 诊断明细；
- T19 的 `prune_matrices` 调用点。

### 7.5 与前三部分的修正

见 [01-review.md](01-review.md) §19：Q01（11→20 条模块再导出）、Q02（新增 `EngineError`、旁路 12 个）、
Q12（含 1 个非错误枚举）、Q24（3→6 个 doctest）、Q29（打包根因、readme 计数）、
Q30（**删除"reed-graph 无 rust-version"这一错误**）、Q35（3/9→5/19）、Q36（重复依赖口径）、
R12（补第三方白名单违规）、R08（第二处取消缺口）、R10/R11（AtomicFile::commit 与 reed-distance writer）。

### 7.6 产物

| 文件 | 说明 |
|---|---|
| [01-review.md](01-review.md) §17–§20 | 工具矩阵、逐 crate 基线、20 条 T 发现、修正表 |
| [02-roadmap.md](02-roadmap.md) 第四部分 | S-A…S-E 批次与 8 条验收门 |
| `docs/README.md` | v30 索引行补充第四部分摘要 |


## 8. 综合诊断（03-synthesis.md）

- **79 条 → 5 个根因 + 1 个产品裁决**：A 契约面 17、B 语义 15、C 状态与资源 18、
  D 演进 16、E 文档与结构 9、F 裁决 4（合计 79）。完整归属表见
  [03-synthesis.md](03-synthesis.md) §2。
- **判定规则**：删除根因后症状按构造消失、同簇共用同一修复动作、有一条工具/编译门可验证。
- **根本方案 P0–P4**：先裁决发布面/稳定面（P0），再收回语义（P1）、能力与资源（P2）、
  门面与六工具门（P3）、派生文档与结构（P4）；顺序与依赖见 03 §5。
- **验收**：12 条工具门（03 §6）；旧条目只在对应门命令通过后关闭，未通过的保留为测试用例。
- 本部分**无新实测**：所有证据引用 §5（Q）、§6（R）、§7（T）与 03 §10 的工具输出。

## 9. 独立复核（2026-11-11 晚，见 [04-recheck.md](04-recheck.md)）

本节不改代码。复核用二进制 `reed 0.4.0 (v0.3.0-180-gf190708-dirty 2026-11-11)`，产物在 `/tmp/v30chk/`。

| 项 | 结果 | 命令/证据 |
|---|---|---|
| V30-01 现象 | 成立：60,346/60,378 对 `0/1/0/0.5` | `reed distance subvcf.vcf.gz --type pihat` |
| V30-01 修法 | **v30 的 `e02` 错**；plink 有限样本系数 + 全位点 `nn` 的 numpy 复算与 plink 1.9 `.genome` 60,378 对 4 位小数 **0 差** | `plink --genome full`；04 §2.1 |
| V30-02 | 成立；无样本变化时默认 `-x` 也删 AC/AN | 六组参数，04 §3.4 |
| V30-04 依据 | 不成立：qpDstat f4mode、qp3Pop、ADMIXTOOLS 2、allel 都是位点均值 | `qpDstat.c:1170-1213` 等 |
| V30-R02 | 成立且更宽：htslib 把超出 int32 的 INFO/AC 置缺失走 GT，reed 用 u64 读——单 ALT 无溢出时 reed AF 0.74、bcftools 0.00；`.afreq`≠`.stats` 是各自对齐参照 | `ovf.vcf`、`ovf1.vcf`、`disagree.vcf` |
| V30-R01 | convert `present` = vcf2phylip 口径；`--012` 半调用 reed 写 −1、vcftools 写 1（parity 脚本已知并剔除） | `half.vcf` |
| **N-1** | **新缺陷**：`--skip-empty` 首个非空窗口 `no_sites` 124→1，π 差 67 倍；10 位点夹具可复现 | `gap.vcf`、`subvcf` |
| N-2 | PI_HAT 的 `nn` 不含单态位点，plink 含（首对 0.6316→0.6652） | `mono.vcf` |
| N-3 | `.genome` 的 DST/PPC/RATIO 恒 0 | `distance.rs:427-432` |
| N-4 | 本机 plink（2025-08）早于上游 2026-09-07 全 ALT 权重修复；全 ALT 位点权重设 1 后差回 5e-7 | `fixalt.vcf` |
| N-5 | `--summary` 无 PSI/HWE 节 | 对 `bcftools stats -s -` |
| 滑窗 | `--indep-pairwise` 5 组参数与 plink2 逐 ID 相同；popg 滑窗坐标与 vcftools 21/21 相同 | 04 §5 |


## 10. 实施记录（2026-11-12，见 [04-implementation.md](04-implementation.md)）

本批落地 P0 契约裁决与 CI 门、P1 的 R01/R02/R03/R04/R12、P2 的
R05/R06/R07/R08/R09/R10/R11/R13/R14，以及 P3 的门面收口与机器门。
P2 §2 的 `RunContext { budget, abort }` 已统一传入 filter/popgen/stats/convert/
rename/pca/dapc/index 八个 driver；`run_budget.rs` 与 `FilterOptions` 的第二预算
删除；workspace bump 到 0.5.0，`reed` 成为发布目标。
验证：`cargo fmt`、`cargo clippy --workspace --all-targets -- -D warnings`、
`cargo test --workspace`（exit 0，全部 test-result 行 passed）、`cargo test --doc --workspace`、
`cargo tree -p reed-filter`（无 `reed-io`）、
`strace -e openat`（`filter -r --strategy single_pass` 对输入只 1 次 open）、
`rg default_half_memory crates/reed/src`（0 处）、
`rg "EditPlan \{" crates/reed/src`（0）、`rg "pub const RECOMPUTE_TAGS" crates`（1 处）、
`ci/check-api.sh`（1976 行，`prelude` 展开；13 别名 + 7 writers 移入 `#[doc(hidden)] frontend`）、
`ci/check-nameable.sh`（22 个签名类型）、`ci/check-external-types.sh`（0 error）、
`ci/check-features.sh`（cargo hack 全绿）、`ci/check-doc-drift.sh`、
`ci/check-structure.sh`（31/31）、`ci/check-docs.sh`（rustdoc/doctest/orphan/入口表）、
`ci/check-semver.sh`（0.5.0 vs v0.3.0 判定 major）、`ci/check-package.sh`（发布目标 + 依赖版本）、
`cargo check -p reed-reduce --no-default-features` 全部通过。P0–P4 全部完成。

## 11. 第二批实施记录（2026-11-13，见 [05-review-followup.md](05-review-followup.md) §4）

按 05 §4 的顺序全部落地；§5 的两个裁决按推荐项执行（phylo 跟 vcf2phylip、
`--fill-tags` 恢复 R13 前的 AC,AN,NS,MAF,HWE 顺序）。

| 序 | 项 | 改动 | 验证 |
|---|---|---|---|
| 1 | N-1 `--skip-empty` | `stats/sink.rs` 跳转后重查 `s + size > watermark` 并 break | 新测试 `skip_empty_keeps_the_full_first_nonempty_window`（1/4/56 线程、4 组 window/step）；移除修复后该测试失败（首窗口 `no_sites` 1 vs 3） |
| 2 | CI 工具链 | `gates` 装 `nightly-2026-03-20` + pinned `cargo-{public-api,hack,semver-checks,modules}`，安装步骤移到调用之前；`check-api.sh` 显式 `+nightly-2026-03-20`；`check-docs.sh` 缺工具时给出安装命令 | `bash ci/check-api.sh` exit 0（快照以同一 nightly 重生成，1976 行）；`bash ci/check-docs.sh` exit 0 |
| 3 | V30-01 PI_HAT | `expected_ibs_coefficients` 换成 plink 的有限样本公式（`plink_calc.c:4848-4867`），分母用参与位点数；N-2：`distance` 路径保留单态位点计入 `nn` | `gold/popgen/run_genome_parity.sh`（新增）：subvcf 与 subvcf+6 单态位点两份夹具，60,378 对 Z0/Z1/Z2/PI_HAT **逐位相同（最大差 0）**；单测 `expected_ibs_coefficients_match_the_plink_finite_sample_formulas`、`pi_hat_counts_monomorphic_sites_in_nn` |
| 4 | V30-02 `-I` | `finish.rs` 在 `no_update && !drop_stale_explicit` 时把 AC/AN 移出 `stale_tags` | 新 e2e `no_update_keeps_input_ac_an`：`-I` 保留 `AC=3;AN=4`、`-I` ≡ `-I -x ''`、默认 update 仍重算 `AC=1;AN=2`；移除修复后该测试失败 |
| 5 | fill 顺序回归 | `RecomputeFlags::enable` push 后按 `RECOMPUTE_TAGS` 位置排序 | 单测 `enable_keeps_registry_order_regardless_of_call_order`；e2e `subset_fill_tags_keeps_registry_order`（`-s s0 --fill-tags MAF` → `AC;AN;MAF`）；移除排序后两测失败（`MAF;AC;AN`） |
| 6 | R01 phylo 门 | `gt.phylo_present = !raw.starts_with(b".")`，只用于 phylo 的 `--min-samples-locus`；count 家族仍用 `has_called_allele` | `gold/fixtures/half.vcf` + `run_convert_parity.sh` 第 7 段：`-m 4` 与 vcf2phylip 逐字节相同（`8 1`）；`convert parity: PASS` |
| 7 | V30-03 多等位 | 每站点报告（`--freq/--hardy/--missing/--depth`）对 `n_alt > 1` 跳过并计数，`StatsOutcome::n_multiallelic_skipped` 暴露给调用者 | e2e `multiallelic_records_are_skipped_and_counted`：`.afreq` 只剩 2 条双等位行、计数 1 |
| 8 | msa 测试 | `e2e_convert.rs` 新增 `msa_region_bounds_are_usage_errors`（`chr1:0-10`、`chr1:10-2`） | 两分支均为 `ErrorKind::Usage`，错误文本固定 |
| 9 | LD idle 取消 | `driver.rs` 在无 trigger handle 时用**同一对** `abort_channel_for` 的 handle+signal | 单测 `dropping_an_unfinished_executor_sets_the_flag_the_workers_poll`：drop 后 worker 侧 `is_aborted()` 为真、后续扫描立即 `aborted`；换回旧接线该测试失败 |
| 10 | AC 单循环 | `info_counts.rs` 合成一次遍历，去掉不可达分支，`got` 语义与旧实现逐值一致 | 既有 6 个 core 单测 + 新增 `a_pathological_ac_field_is_cut_off_at_the_first_excess_element`（10,000 元素在第二个处返回） |
| 11 | Framing | `input.rs` 的区域能力检查改读 `Framing::of(format)`，不再直接 `matches!` | `rg "\.framing\(\)" crates/reed/src` 与生产路径一致；`regions_on_plain_vcf_are_a_capability_error` 仍通过 |
| 12 | 描述/文档 | `crates/reed/Cargo.toml` 描述不再列举函数名；`stats.md` 改多等位跳过说明并登记 `--summary` 的 AC 回退；`convert.md` 写明 phylo 的 presence 规则 | `ci/check-docs.sh` 入口表门通过 |

保留未做（05 未要求或属规范选择）：
- `.genome` 的 DST/PPC/RATIO 仍写 0（`04-recheck.md` N-3，未在 05 §4 排入）。
- bcftools 的 fill 标签顺序本就与 reed 不同（`NS→AN→AF→AC→HWE`）；05 §3.2 明确只要求恢复 reed 自身旧序。
- `--summary` 对 `AC=.` / 超范围 AC 回退到 GT（05 §3.5 要求登记而非跟随 htslib 怪癖）。

验证命令与结果（本机）：
```
cargo fmt --all -- --check                  exit 0
cargo clippy --locked --workspace --all-targets -- -D warnings   exit 0
cargo test --locked --workspace             全部 test-result ok
cargo test --release --locked --workspace   全部 test-result ok
cargo test --locked --doc --workspace       exit 0
bash ci/check-api.sh                        exit 0（快照 1976 行，新增 StatsOutcome::n_multiallelic_skipped）
bash ci/check-docs.sh                       exit 0
bash ci/check-structure.sh                  exit 0
bash ci/check-nameable.sh                   exit 0
bash gold/popgen/run_genome_parity.sh       PASS（两份夹具 0 差）
bash gold/popgen/run_convert_parity.sh      PASS（含新的 half-call 门）
```

## 12. 第三批（2026-11-15）：N-3 `.genome` 的 DST/PPC/RATIO

`04-recheck.md` N-3 记录 `.genome` 的 DST/PPC/RATIO 恒写 0。本轮按该条的两个选项中
更诚实的一支处理：

- **DST 按 plink 的公式实现**（`plink_calc.c:4454`，`1 − (IBS1 + 2·IBS0)/(2·nn)`），
  `reed-distance::genome_columns` 暴露；`AuxCounts` 增加 `n_valid`/`n_ibs1`/`n_ibs0`。
- **PPC/RATIO 输出 `NA`**，不再伪造数值。原因写进代码注释与 `distance.md`：plink 的
  `HOMHOM`/`HETHET` 是 `g_marker_window`（固定大小块缓冲，`plink_calc.c:756`）**最后一个
  窗口**的命中数，不是全位点统计——同一份 3,000 位点夹具上 plink 对每一对都打印
  `HOMHOM=0 HETHET=2`。这两个数取决于 plink 的内部块划分，无法用全位点计数复刻。
  写一个"看起来合理"的 PPC 会得到与 plink 不同的数字，比 NA 更糟。

验证：
- `gold/popgen/run_genome_parity.sh` 扩展为同时比对 DST：
  `subvcf` 与 `mono` 两份夹具各 60,378 对，`Z0/Z1/Z2/PI_HAT` 与 `DST` 全部
  `max_abs_diff = 0`，`dst_mismatch = 0`。
- 单测 `genome_dst_is_plinks_ibs_distance_and_ppc_ratio_are_na`（手算 DST=0.74；
  空 `nn` 时 DST 也 NA）。
- `docsite/src/content/docs/distance.md`：删掉"placeholder 0"的说法，改述 DST 已对齐、
  PPC/RATIO 为 NA 及原因；`kin.genome` 示例换成真实输出。
- `cargo clippy --workspace --all-targets -D warnings` exit 0；
  `cargo test --workspace` 69 个 test-result 全 ok、0 失败；
  `ci/check-api.sh`、`check-doc-drift`、`check-structure`、`check-nameable`、
  `check-docs`、`check-package` 全部 exit 0（API 快照零漂移）。

仍未做：`04-recheck.md` 的 N-4（对拍基线要固定含 9-07 修复的 plink 1.9，属环境要求）
与 N-5（`--summary` 缺 PSI/HWE 两节，需另开一节补实现）。

## 13. 第四批（2026-11-15）：N-5 `--summary` 节清单登记

`04-recheck.md` N-5：`--summary` 缺 `PSI`/`HWE` 两节且文档未写明。按该条给出的
"或"分支处理（补节需要 bcftools `-s -` 那样的逐样本/逐位点独立遍历，与 reed 的
单次共享扫描设计冲突，属功能扩展而非本轮修复）：

- `docsite/src/content/docs/stats.md` 新增 "`--summary` sections" 小节，逐项列出
  已实现的 9 节（SN/TSTV/SiS/AF/QUAL/IDD/ST/DP/PSC），并明确 **不输出** `PSI` 与
  `HWE`，给出替代（`--hardy` 提供逐位点 HWE 表；`--het` 提供逐样本近交系数；
  `PSI` 无等价物）。这消除"文档只说 QC summary → PREFIX.stats、用户不知道缺什么"
  的隐性缺口。

至此 `04-recheck.md` 的 N-1/N-2/N-3/N-5 全部有结论；N-4 是环境要求（对拍基线需固定
含 2026-09-07 修复的 plink 1.9），不改代码。

## 14. 当前权威状态（2026-11-19）

本节只做索引，详细结论、证据、剩余问题和发布流程见
[06-final-audit.md](06-final-audit.md)。

- IO、filter、并行/取消：本次未发现新的 P0/P1 运行时 bug。
- V30-01、V30-02、V30-03、N-1、N-2、N-3、N-5、fill-tags 顺序、R01 phylo 口径和
  LD idle 取消：已修复、已测试或已明确文档化。
- `.genome` 的 PPC/RATIO/RT/EZ/PHE、summary 的 PSI/HWE、HWE mid-p、site-count window：
  是明确的功能边界，不应写成已完成的完整兼容。
- 仍需处理的工程项：发布脚本的 clean/并发保护、`run_distance` facade 决策、opened
  filter 入口、支持 crate 稳定面、`frontend` 逃生舱和 header 重复解析优化。
