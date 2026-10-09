# 06 目录树规范（v5 执行版；承自 doc_v4/07-layout.md，内容同步）

> 用户要求：**子文件夹树混乱**，要一版可直接改的架构树。本文范围是 `reed-rs/` 仓库内；§2.2 的 `~/reed/` 为可选附录。v5 把它列为路线 X2（07-roadmap.md）。
> 本文只写方案；实现者按 §5 执行。**范围是 `reed-rs/` 仓库内的树**；§2.2 与 §4 表末四行关于 `~/reed/` 仓库外的部分是可选附录，不在用户要求内。规则：**代码在 `crates/`，验证在 `gold/`，度量在 `bench/`，规范在 `docs/`，生成物一律进 `target/`，大数据与参照仓库在仓库外**。

---

## 1. 现状（实测）

| 顶层条目 | 内容 | 体积 / 问题 |
|---|---|---|
| `crates/` | 8 个 crate | 含 `reed-core-v02-oracle`（**版本后缀**，与"包/类型不带版本后缀"的命名规则冲突） |
| `bench/` | 独立 cargo workspace：`e2e/`、`hts/`、`noo/`、`run_matrix.sh`，自带 `target/` | `hts`/`noo` 是缩写；`bench/target` 131 M |
| `docs/` | `architecture-spec.md`(1936 行)、`filter-architecture.md`(412)、`naming-conventions.md`(320)、`statistics-reference-map.md`(159)、`INDEX.md`(43) | 与 `doc_vN/` 并列，读者不知道先看哪个 |
| `doc_v2/ doc_v3/ doc_v4/` | 三轮复审与裁决 | **三个顶层目录**只为历史；每轮又各有一份 `README/STATUS` |
| `DESIGN.md` | v1 设计（49 KB） | 与 `docs/architecture-spec.md` 重复，无人知道哪个有效 |
| `gold/` | `config/`、`scripts/`、`popgen/`、`README.md`，以及 `work/` **54 G**、`out/` **25 G** 生成物 | 夹具、脚本、79 G 数据混在一个目录；`gold/out/work/` 还是嵌套的同名目录 |
| `scripts/` | `check-doc-drift.sh`、`publish.sh` | 名字是杂物抽屉；无 CI 目录 |
| `README.md` / `README.zh-CN.md` | 双语 | 没问题，保留（GitHub 惯例） |
| `target/` | 372 M | 已 gitignore |
| 根配置 | `Cargo.toml/lock`、`rust-toolchain/rustfmt/clippy/deny.toml`、`.cargo/config.toml` | 没问题 |

仓库外（`~/reed/`）：16 个参照仓库（bcftools、htslib、noodles、pixy、plink-ng、scikit-allel、vcftools、sgkit、snputils、vcflib、vcfstats、GCTA、flashpca、easySFS、SNPGenie、rust-htslib）、输入数据 `clean.vcf.gz`(7.0 G)、`clean.bcf`(8.1 G) 与索引、`clean_bg.*`(256 M)、`make_3pop.py`、`clean.3pop.txt`，以及一个 **20 G 的散落文件 `=30`**（重定向事故），全部与 `reed-rs/` 平级。

**问题归成四类**：① 文档有两个根（`docs/` + `doc_vN/` + `DESIGN.md`）；② 生成物与脚本同目录（`gold/`，79 G）；③ 顶层目录名不说明用途（`scripts/`、`hts/`、`noo/`）且 `gold/out/work` 套娃；④ 仓库外把代码、参照仓库、15 G+ 数据、垃圾文件混在一层。

---

## 2. 目标树

### 2.1 `reed-rs/`（6 个顶层目录）

```
reed-rs/
├─ Cargo.toml  Cargo.lock
├─ rust-toolchain.toml  rustfmt.toml  clippy.toml  deny.toml
├─ .cargo/config.toml                  # 机器相关，gitignored（保留）
├─ README.md  README.zh-CN.md          # 双语留在根（GitHub 惯例）
├─ crates/
│   ├─ reed-core/                      # 类型/数学/计数（不动）
│   ├─ reed-io/                        # BGZF/VCF/BCF 编解码（§2.3 C-4：noodles 编码器退出生产路径）
│   ├─ reed-engine/                    # 并行引擎（不动）
│   ├─ reed-filter/                    # 过滤语义（不动）
│   ├─ reed-stats/                     # 统计量（§2.3 C-5：individual.rs 拆开）
│   ├─ reed/                           # 门面 lib（§2.3 C-1/3/6/7：删 M3 qc、qc→summary、顶层散件归位、job.rs 拆分）
│   ├─ reed-cli/                       # 唯一的二进制（§2.3 C-2/8：删 `reed qc`、删未用依赖）
│   └─ oracle/                         # ← crates/reed-core-v02-oracle（去版本后缀，publish = false）
├─ bench/                              # 度量：只放代码与结果日志，不放构建产物
│   ├─ suites/
│   │   ├─ e2e/{results.tsv,timing/}   # ← bench/e2e（results.tsv 是要跟踪的度量台账）
│   │   ├─ vs-htslib/                  # ← bench/hts
│   │   └─ vs-noodles/                 # ← bench/noo
│   └─ run_matrix.sh
├─ gold/                               # 验证：夹具规则与对拍脚本（不含大产物）
│   ├─ README.md  config/filters.conf
│   ├─ scripts/                        # run_{reed,vcftools,plink2,bcftools}.sh、check_parity.py、lib.sh…
│   ├─ popgen/                         # run_pixy.sh、compare*.py、make_fixtures.py、make_3pop.py ← 仓库外
│   └─ fixtures/                       # 新增：小体积边界夹具，随仓库跟踪（02 §7.1 的 N=10/7 等）
├─ docs/
│   ├─ README.md                       # ← INDEX.md：先看这里，一页说清谁取代谁
│   ├─ spec/                           # 活规范（唯一"当前有效"入口，每轮先改这里）
│   │   ├─ architecture.md             # ← docs/architecture-spec.md
│   │   ├─ filter.md                   # ← docs/filter-architecture.md
│   │   ├─ cli.md                      # ← doc_v2/02 + doc_v4/02 合并后的最新版（待办）
│   │   ├─ naming.md                   # ← docs/naming-conventions.md
│   │   └─ statistics-map.md           # ← docs/statistics-reference-map.md
│   └─ rounds/                         # 历史，只增不改
│       ├─ v1/{design.md}              # ← DESIGN.md
│       ├─ v2/  v3/  v4/               # ← doc_v2/ doc_v3/ doc_v4/
└─ ci/                                 # ← scripts/（+ .github/workflows/ci.yml）
    ├─ check-doc-drift.sh
    └─ publish.sh
```

生成物去向：`target/`（已 gitignore）。`gold/work/`、`gold/out/` 与 `bench/target/` 全部移到 `target/gold/{work,out}`、`target/bench/`，脚本用 `REED_SCRATCH`（默认 `target/gold`）定位，不再在仓库里堆 79 G。

### 2.3 `crates/` 内部（用户追问"现在不乱么"）

**不乱的部分。** crate 的集合与分层与 doc_v2 03 §7 完全一致，依赖图是严格的单向分层，没有环、没有越层：

```
reed-core ← reed-io ← reed-engine
reed-core ← reed-stats
reed-core, reed-io ← reed-filter
全部 ← reed（门面 lib） ← reed-cli（唯一二进制）
oracle：仅作 reed-core / reed-stats 的 dev-dependency
```

体量（src 行数）：reed 7.5k、reed-io 7.0k、reed-filter 4.2k、reed-stats 4.0k、reed-core 4.0k、oracle 2.9k、reed-cli 1.9k、reed-engine 1.5k。

**乱的部分**是三类历史残留（两代并存、名字漂移、散件），不是分层错误：

| # | 位置 | 问题（已核实） | 处理 |
|---|---|---|---|
| C-1 | `reed/src/stats/writers/qc.rs`（146 行，自述 "M3 first consumer"） | 与 `qc_report.rs`（759 行，M6 版）是同一功能的两代；`qc.rs` 只在 `writers/mod.rs:5` 被 `pub mod` 声明，无调用者 | 删除 |
| C-2 | `reed-cli/src/cli.rs:49`、`translate/qc.rs`、`main.rs:25`：`reed qc` 子命令 | 已标 deprecated，注释还指向已合并掉的 reed-qc crate；CLI 规范只有 filter / stats / index | 删子命令与 `translate/qc.rs`；`reed stats --summary` 已覆盖 |
| C-3 | "qc" 这个词 | CLI 叫 `--summary`，库叫 `reed_stats::summary`，驱动与产物却叫 `qc_report.rs`、`<prefix>.qc.json` | `qc_report.rs` → `stats/writers/summary.rs`；`.qc.json` → `.summary.json` |
| C-4 | `reed-io/src/bcf/encode.rs`（noodles 编码器，doc_v3 B2） | B4 已用 `native_encode.rs` 取代，但 `reed/src/filter/single_pass.rs:31` 仍以 `Option<BcfEncoder>` 带在生产路径 | 生产路径只留 `NativeEncoder`；noodles 版只作 `reed-io/tests/bcf_encode.rs` 的 oracle（移进 `tests/` 或 `#[cfg(test)]`） |
| C-5 | `reed-stats/src/individual.rs` | 一个文件装两个统计量（individual F 与 joint SFS），违反"每个统计量一个文件、按方法命名"（doc_v2 03 §7） | 拆为 `het.rs`（individual F，plink `--het` 的矩估计）与 `sfs/{projected,joint}.rs` |
| C-6 | `reed/src/{filtered.rs, bcf_input.rs, input.rs}` 散在 crate 顶层 | 顶层应只有三条命令的驱动（`filter/`、`stats/`、`index.rs`）加输入层 | `input/{mod.rs, bcf.rs}`；`filtered.rs` → `stats/filtered.rs`（它是"边过滤边统计"的 stats 侧包装） |
| C-7 | `reed/src/stats/job.rs` 1,167 行 | 仓库最大文件：建请求、跑流水线、写输出三件事在一起 | 拆 `stats/{request.rs, run.rs, sink.rs}`，与 `filter/` 的 `single_pass / reduce / emit` 对称 |
| C-8 | `reed-cli/Cargo.toml` | 声明依赖 reed-core / reed-io / reed-stats，源码只 `use reed_engine`、`reed_filter`（其余经 `reed` 门面） | 删未用依赖 |
| C-9 | `reed-core-v02-oracle` | 版本后缀 | `crates/oracle`，包名 `reed-oracle`，`publish = false`（§4） |

**不改**：`reed-core/src/grid.rs`（窗口网格模型：`WindowGrid`、`BpWindow`、`RegionSet`）与 `reed-stats/src/window.rs`（按网格分箱累加：`BpBinning::from_grid`、`ContigBins`）是分层不是重复；`reed-filter/src/{sketch,funnel}.rs` 放顶层可接受；`reed-io` 的 `bgzf/ chunk/ vcf/ bcf/ write/` 分法清楚。

**顺序与门**：C-1 / C-2 / C-8 / C-9 零风险，随目录搬迁同一个 PR；C-3 / C-5 / C-6 只改名与归属，门是 `cargo test` 全绿 + 全部输出 md5 不变；C-4 / C-7 触及生产路径，各自单独 PR，门是 gold 正文 md5 与 `-Ob` 正文 sha256 不变。

### 2.2 附录（可选，不在本次要求内）：`~/reed/` 仓库外

```
~/reed/
├─ reed-rs/          # 唯一的代码仓库（真身）
├─ data/             # 输入：clean.vcf.gz、clean.bcf、*.csi（只读）
├─ refs/             # 16 个参照仓库，只读（bcftools/htslib/noodles/pixy/plink-ng/scikit-allel/vcftools/…）
├─ work/             # 中间产物与日志：clean_bg.*、clean.3pop.txt、clean_bg.log
└─ README.md         # 一页：每个目录是什么、数据从哪来、refs/ 为何在仓库外
```

`=30`（20 G）是 `cmd >=30` 之类的重定向事故，确认后删除。

---

## 3. 放置规则（以后新东西放哪）

| 新东西 | 去处 |
|---|---|
| 新 crate / 新模块 | `crates/<crate>/src/<domain>.rs`；crate 名不带版本后缀与 `vN` |
| 新的端到端测试 | `crates/reed/tests/`（等价/确定性）或 `gold/fixtures/` + `gold/scripts/`（对拍参照物） |
| 新的一轮复审 | `docs/rounds/v<N>/`，同时把结论落进 `docs/spec/` 对应文件 |
| 该轮实现状态 | `docs/rounds/v<N>/STATUS.md` |
| 基准结果 | `bench/suites/*/results.tsv`（跟踪）与 `target/bench/`（大文件，不跟踪） |
| 生成物（VCF/BCF/索引/日志） | `target/` 或 `$REED_SCRATCH`，永不入库 |
| 大数据输入 | `~/reed/data/`（仓库外） |
| 参照物源码 | `~/reed/refs/`（仓库外，只读） |
| CI 校验脚本 | `ci/`（由 `.github/workflows/ci.yml` 调用） |

---

## 4. 现状 → 目标（映射）

| 现状 | 目标 |
|---|---|
| `crates/reed-core-v02-oracle/` | `crates/oracle/`（`Cargo.toml` 里 `name = "reed-oracle"`，或只改目录名保留包名） |
| `bench/e2e/`、`bench/hts/`、`bench/noo/` | `bench/suites/{e2e,vs-htslib,vs-noodles}/` |
| `bench/target/` | 删除（改由根 `target/` 或 `CARGO_TARGET_DIR` 承载） |
| `docs/architecture-spec.md` 等 4 个 + `INDEX.md` | `docs/spec/{architecture,filter,naming,statistics-map}.md` + `docs/README.md` |
| `DESIGN.md` | `docs/rounds/v1/design.md` |
| `doc_v2/` `doc_v3/` `doc_v4/` | `docs/rounds/{v2,v3,v4}/` |
| `gold/work/`（54 G）、`gold/out/`（25 G） | `target/gold/{work,out}` |
| `gold/out/work/` | 去掉套娃，统一 `target/gold/work/` |
| `scripts/` | `ci/` |
| `~/reed/make_3pop.py` | `gold/popgen/make_3pop.py` |
| `~/reed/clean*` | `~/reed/data/`（输入）、`~/reed/work/`（产物） |
| `~/reed/{bcftools,htslib,…}` 16 个 | `~/reed/refs/` |
| `~/reed/=30` | 删除 |

命令（一次性，`git mv` 保留历史）：

```sh
cd ~/reed/reed-rs
git mv crates/reed-core-v02-oracle crates/oracle
mkdir -p bench/suites docs/spec docs/rounds/v1 ci
git mv bench/e2e bench/suites/e2e && git mv bench/hts bench/suites/vs-htslib && git mv bench/noo bench/suites/vs-noodles
rm -rf bench/target
git mv docs/architecture-spec.md docs/spec/architecture.md
git mv docs/filter-architecture.md docs/spec/filter.md
git mv docs/naming-conventions.md docs/spec/naming.md
git mv docs/statistics-reference-map.md docs/spec/statistics-map.md
git mv docs/INDEX.md docs/README.md
git mv DESIGN.md docs/rounds/v1/design.md
git mv doc_v2 docs/rounds/v2 && git mv doc_v3 docs/rounds/v3 && git mv doc_v4 docs/rounds/v4
git mv scripts ci
mkdir -p gold/fixtures
sh gold/scripts/move_scratch.sh          # 见 §5.3：把 work/ out/ 挪到 target/gold/
```

---

## 5. 要同步改的地方

### 5.1 Cargo 与构建

- 根 `Cargo.toml` 的 `members`（若按目录列举）与 `bench` 的成员路径；`bench/*/Cargo.toml` 里 `path = "../../crates/…"` 的相对层级（`suites/` 多了一层）。
- `crates/oracle` 改名后同步所有 `reed_core_v02_oracle::` 引用与依赖别名（`grep -rn reed_core_v02_oracle crates`）。
- `.gitignore`：增加 `/target/gold/`；删除 `/gold/` 这一条（夹具目录改为**跟踪**，只忽略 `work/`、`out/`）；`bench/target` 条目可删。

### 5.2 文档引用（工作量最大，建议一次 sed）

`docs/spec/*` 与各轮 README 里大量 `doc_v2/02-cli.md` 形式的交叉引用要变成 `docs/rounds/v2/02-cli.md`：

```sh
grep -rln 'doc_v[0-9]' docs | xargs sed -i -E 's#\bdoc_v([0-9])/#docs/rounds/v\1/#g'
```

随后在 `ci/check-doc-drift.sh` 加一条守卫：`! grep -rq 'doc_v[0-9]/' docs/README.md`，防止回潮。

### 5.3 脚本相对路径

`gold/scripts/{lib.sh,run_all.sh,run_*.sh,check_parity.py,compare_tools.py,profile_input.py}` 与 `gold/popgen/*.py` 里的 `out/`、`work/`、`../out` 等相对路径改为 `"$REED_SCRATCH"`（默认 `target/gold`）；新增 `gold/scripts/move_scratch.sh` 做首次搬迁；`gold/README.md` 写清 `export REED_SCRATCH=target/gold` 的用法。`ci/` 下两个脚本的 `cd "$(dirname "${BASH_SOURCE[0]}")/.."` 不变（仍在仓库根下一层）。

### 5.4 其他

- 新增 `.github/workflows/ci.yml`：`cargo test --release --workspace`、`cargo clippy`、`bash ci/check-doc-drift.sh`。
- `README.md` 的"项目结构"段同步成 §2.1 的树。
- `docs/spec/cli.md` 的**合并是待办**（把 doc_v2/02 与 doc_v4/02 回退版合成一份当前有效版），不要只靠 `git mv` 假装完成。

---

## 6. 注意事项

1. **CodeGraph 索引**：`~/reed/refs/*` 的 16 个仓库各自带 `.codegraph/`。移动后先在其中一个新路径跑 `~/miniforge3/bin/codegraph explore "<符号>"` 抽查；若索引存绝对路径则需重新 `codegraph init`。**先试一个，再全挪**。
2. **数据用 mv 不用 cp**：`clean.vcf.gz`(7.0 G)、`clean.bcf`(8.1 G) 与 `.csi` 一起移动到 `~/reed/data/`。
3. **`=30`（20 G）**：先 `ls -l ~/reed/=30` 确认，再删除；删除不可逆。
4. **顺序**：先 `git mv` 小目录（`crates/oracle`、`bench/suites`、`docs`），跑一次 `cargo test`；再动 `docs/rounds` 与 sed 引用；最后动 79 G 的 `gold/work|out`（只影响脚本，不影响 git）。
5. **不做**：不动 crate 内部模块结构（本轮只动归属与名字）；不把历史轮次的内容合进 `docs/spec/`（§5.4 的单列待办）。
