# 06 facade 与输出（V25-08 / V25-09 / V25-10）

> facade 的边界承诺（"二进制只依赖 reed 一个 crate"）在加入 distance/reduce/graph
> 时被击穿且未回填：凡从直连路径长出来的功能，每个 bug 只能修在 CLI translate 层。
> 本轮回填边界，并把"表格/报告输出"的写出惯例收敛到一层。

---

## 1. V25-08a：facade 回填与错误分类

### 现状【据】

- `reed-cli/Cargo.toml:16-18` 直连 `reed-distance` / `reed-reduce` / `reed-graph`；
- `reed::Error` 只有 4 个 sealed 载荷（Core/Io/Engine/Spec）+ `Usage(String)`
  （`error.rs:22-35`）；graph/reduce 的错误在调用点全用
  `map_err(|e| Error::usage(e.to_string()))` 字符串化
  （`translate/distance.rs:213/258/284`、`dapc/mod.rs:388-390`【述】）；
- 后果：`GraphError::NonFiniteDistance / Asymmetric / TreeTooDeep`、
  `LdaError::Degenerate` 是**由输入数据决定的运行期失败**，却报 exit 2——
  违反自家契约 `error.rs:31-33`（"usage = the operator can fix on the command
  line"）；脚本按 exit 2 改命令重试会原地死循环；
- 同类：输出文件创建失败的分类随命令翻转——rename 记 usage（`rename/mod.rs:649-652`
  【述】），filter/convert/dapc 走 runtime【述】——`mkdir x && reed rename -o x/out.vcf`
  exit 2，同输入 filter exit 1。
- R11：stats 的 usage/runtime 退出码随输入文件是否存在而翻转
  （`translate/stats.rs:206-225`）：报告模式的部分 usage 校验（`-s/-S`、
  site/threshold、`--hist` 空表）在 `prelude()` 之后，输入不存在时先以 runtime
  失败，输入存在时才报 usage；同一命令行两种退出码。

### 规范【规】

1. `reed` lib.rs re-export：`pub use reed_distance as distance; pub use reed_reduce
   as reduce; pub use reed_graph as graph;`；distance 的请求构造从
   `translate/distance.rs` 下沉进 facade（照 pca/dapc 的 `run_*_opened` 模式），
   CLI 回到只依赖 reed；
2. `Error` 增加 `Foreign(String)`（或 `Runtime`）载荷，走 Runtime / exit 1；
   调用点按 variant 分类：
   - usage（exit 2）：`TooFewTaxa`、`NotSquare` 等"命令行可修"的；
   - runtime（exit 1）：`NonFiniteDistance`、`Asymmetric`、`TreeTooDeep`、
     `LdaError::Degenerate` 等"数据决定"的；
3. "输出文件创建失败"统一分类：收敛到 §3 的共享 open helper 后**只选一次**
   （推荐 runtime——文件系统问题不是命令行拼写问题）；
4. **R14**：piHat 表格的 `rsplit_once('_')` 伪造 FID/IID
   （`translate/distance.rs:394-399`【据】）改为两列都写样本名；样本名本身含
   `_` 时不得再猜分界。
5. **R11**：report mode 的全部 usage 校验（`-s/-S`、site/threshold、`--hist`
   空表）前移到 `prelude()` 之前；`prelude()` 只负责打开输入。先判命令行，
   再判输入——同一命令行不再因输入文件是否存在而换退出码。
6. `Error::Foreign(String)` 会丢掉 variant：优先 `Box<dyn Error + Send + Sync>`
   或 `enum ForeignError { Graph(GraphError), Reduce(ReduceError), … }`，至少让
   调用方按 variant 分类；若坚持 String，在 `error.rs` 写明"不可匹配"，并在
   translate 层一次性分类，不得散落 `map_err(Error::usage(e.to_string()))`。

### 验收

- NaN 距离 `--tree nj` → exit 1；taxa 不足 → exit 2；`mkdir x && reed rename /
  filter -o x/out.vcf` 退出码一致；
- R11：输入不存在 + 非法 `-s` → exit 2（usage 优先）；输入存在 + 非法 `-s` →
  exit 2；输入不存在 + 合法参数 → exit 1；
- R14：piHat 输出前两列是样本名本身，含 `_` 的名字不被拆开；
- `reed-cli/Cargo.toml` 无三个直连依赖；
- 退出码 e2e 表：每命令 × 代表性错误（usage / runtime / io）各一行。

---

## 2. V25-08b：StageMetrics 渲染单点

### 现状【据】

同一 `StageMetrics` 四套渲染，且**已经漂移**（不是未来风险）：

| 渲染点 | 键名 | 缺陷 |
|---|---|---|
| `reed-filter/src/funnel.rs:129-158` | `reader_ms / chunks_read / worker_busy_ms / chain_wait_ms / inflate_ms …` | 最全（事实标准） |
| `translate/common.rs:15-24`（popg/stats） | `engine.reader_ms / engine.chain_wait_ms …` | 前缀风格不同 |
| `translate/convert.rs:187-197 / :467-478` | **`chunks` ≠ `chunks_read`**；缺 `worker_busy_ms` / `chain_wait_ms`；硬编码 `phase.open_ms 0.0`【述】 | 两份内联拷贝 |
| `translate/index.rs:77-82` | `total_ms / n_records`（其他命令是 `n_records_in`） | 又一套 |

解析 `--timing` 的性能回归脚本必须按子命令写四个解析器。

### 规范【规】

- `StageMetrics` 自己提供 `render_tsv(phases)`（放 reed-engine 或 reed-io），
  四处渲染点全部改为传 phase 前缀的调用；
- 键名规范以 filter 版为准（`chunks_read` / `worker_busy_ms` / `chain_wait_ms` /
  `n_records_in` …），写成一张键名表进 docsite 的输出约定页；
- convert 的两份内联拷贝删除，`open_ms 0.0` 硬编码与缺失字段随之修复。

### 验收

- `grep 'chunks_read\|"chunks"' crates/reed-cli/src/translate/` 键名唯一；
- 四命令的 `--timing` 输出键集合一致（e2e 断言）。

---

## 3. V25-08c：provenance 提升【述→规】

现状：filter 的 `Provenance`（`plan/provenance.rs:44-70`，含 sha256 / canonical
JSON）是权威实现；rename 手写 `##reed_command` 行（`rename/mod.rs:528-535`）、
hapmap 又一处（`convert/hapmap.rs:546`）；argv 捕获
`std::env::args().collect::<Vec<_>>().join(" ")` 在 CLI 复制 5 份且**无引号**
（`-o "my out"` 记录的命令行不可回放，provenance 失效）；pca / distance / dapc /
stats 的输出完全没有可复现性头。

【规】

1. `Provenance` 从 reed-filter 提升为 `reed-core::provenance`（command + version +
   可选 spec）；
2. argv 捕获与 shell 引号在 `reed-cli` 单点实现一次，五个复制点收敛；
3. 【荐】所有 prefix 型输出命令默认挂头（`--no-command/--no-version` 的过滤
   语义只有一份）。

---

## 4. V25-09：TempSpill 共享

### 现状【据】

- `input.rs:202-276`：stdin spill 完整加固——`create_new`（O_EXCL）+ 计数器
  重试 + "不跟随 symlink"论证 + 攻击测试（`:416-476` 两个用例）。这是 doc_v20
  SEC-1 的落点；
- `convert/job.rs:380-384`：`spill_path()` 生成可预测名
  `/tmp/reed-convert-{pid}-{n}.spill`，由 `RowsBlob::append/spill_now`
  （`:299-303 / :350-356`【述】）用普通 `File::create` 打开——**该文件
  `create_new` 出现 0 次**。SEC-1 修在了两份实现中的一份上。

攻击路径：共享主机上预置 `/tmp/reed-convert-{pid}-{n}.spill` symlink → reed
跟随 symlink 截断/写入受害者文件；spill 内容随后流回输出文件（`write_to`），
构成"输入 → 输出被投毒"。

### 规范【规】

1. `create_spill_file` 提升为 reed 内共享的 spill 基础设施（`TempSpill` 类型或
   自由函数），保留 O_EXCL + 重试语义；
2. convert `RowsBlob` 改用它；删除 convert 私有 `spill_path()`；
3. 攻击测试直接移植 input.rs 的两个用例（预置 symlink → 创建失败而非跟随）。

### 验收

- `grep 'temp_dir\|/tmp' crates/reed/src/convert/job.rs` 无裸路径拼接；
- 移植的攻击测试在 convert 路径通过；
- 下一个需要临时文件的功能只能用 `TempSpill`（无第三份裸实现长出来）。

---

## 5. V25-10：Lanczos 观测与收敛判据

### 现状【据】

- 残差按相对值计算：`worst = num.sqrt() / vals[c].abs().max(1e-300)`
  （`lanczos.rs:461`），文档约定 tol 作用于 `||Av-λv||/|λ|`——请求的 `k` 超过
  `W W'/m` 的秩时（重复个体、LD 剪枝后位点数 < `--pcs`），第 rank+1 个特征值
  λ≈1e-17，相对残差 ≈ λ_max·eps/λ ≫ 1e-10，**永远到不了 tol**；
- 两个后果：(a) 每次运行固定烧满 `max_cycles=200`（每轮 O(n·m²) 完整再正交化
  + 全基投影），大 n 下分钟级空转；(b) `StreamEig.max_residual`（`:521` 构造
  返回）被调用方丢弃——pca/mod.rs **0 处引用**，非收敛结果无声写出；
- `STREAM_DEBUG` env 后门（`:464-468`）绕过 verbosity/tracing，用户设了该变量
  任何命令都可能喷调试行，不在文档契约里。

### 规范【规】

1. 残差分母改混合尺度：`vals[c].abs().max(lambda_max * f64::EPSILON)`——
   近零特征值的秩亏分量不再卡死收敛；
2. `max_residual > cfg.tol` 时：stderr 一行告警 + 进 report（可 grep 的字段，
   走 §2 的 metrics 通道）；
3. 删除 `STREAM_DEBUG` 块，改 `tracing::debug!`（verbosity 已初始化
   subscriber，成本为零）；env 检查移出循环【述】；
4. 轮数变化的前后对比**只认 V24-02 协议**（`03-bench-protocol.md`：warmup、
   ≥5 次中位数、JSON 输出），不写一次性结论进 changelog。

### 测试

- 秩亏 fixture（n_samples > sites，如 1000 样本 8 位点 `--pcs 10`）：协议下
  轮数显著低于 200；
- 非收敛输入（tol 设不可达值）→ stderr 有告警（e2e 断言）；
- `grep STREAM_DEBUG` 全仓清零。

### 附（R1，独立小项）

stream 求解器是否绕过 `SymEigResult::top_k` 的符号归一（`backend.rs:55-64`
确有"最大幅值分量为正"规则【据】；lanczos 路径未见处理）：先实跑
`--solver dense` vs `--solver stream` 对比 `.eigenvec`——分叉则在
`lanczos_top_k` 返回前套用同一规则（抽 `fix_column_signs` 共用），并补
stream↔dense 一致性测试（当前 `matches_dense_within_tolerance` 只比特征值）。

---

## 实施记录（2026-11-04）

### V25-08a（提交 `b9353c7`）

- facade：`reed` 增加 `pub use reed_distance as distance; ... graph; ... reduce;`，
  `reed-cli` 删除三个直连依赖（`rayon` 保留）；`ci/check-structure.sh` 加门。
- `Error::Foreign(Box<dyn Error + Send + Sync>)` + `Error::foreign`，kind 恒为
  Runtime；graph 错误按 variant 分流（`NotSquare`/`TooFewTaxa` usage，其余
  runtime），LDA 的 `Degenerate`/`ConstantWithinGroups` runtime，其余 usage；
  rename 的输出创建失败改为 runtime（与 filter/convert 一致）。
- R11：stats 的 report-mode usage 校验（`-s/-S`、site/threshold、空 `--hist`）
  前移到 prelude 之前。
- R14：piHat 表的 FID/IID 两列都写样本名，不再按 `_` 猜分界。
- e2e：too-few-taxa=2、rename/filter 输出不可写=1、R11 三种输入/参数组合
  =2/1/2、piHat 两列逐对断言。
- **08b（提交 `c5c7829`）**：`StageMetrics::render_tsv` 成为唯一引擎时序
  渲染器（11 个规范键）；filter timing、popg/stats、convert 两份内联拷贝、
  index 全部改走；convert 的 `phase.open_ms 0.0` 占位删除、重复 `chunks` 键
  改为 `chunks_read`；`IndexOutcome` 携带 metrics。e2e 五命令
  `--verbosity 2` 断言引擎键集合一致。
- **08c（提交 `c5c7829`）**：`Provenance` 提升到 `reed-core::provenance`
  （`new` + `with_spec`，canonical JSON/sha256/`##reed_*` 名称不变）；
  `params::command_line` 单点捕获 argv 并做 shell 引号（`-o "my out"` 可回放），
  五个复制点删除。

### 08 ③（提交 `0373ed5`）

- `RowsBlob`/`SpilledMsa`/`MsaFrameReader` 纯搬到 `convert/spill.rs`，
  `writers.rs` 改导入；convert 字节一致性测试全绿。结构门新增
  "convert spill types live in convert/spill.rs"。

### V25-09（提交 `db5fff5`）

- `reed/src/spill.rs`：`create_spill_file` + `TempSpill`（O_EXCL、计数器重试、
  Drop 删除）；stdin 与 convert `RowsBlob` 全部改走；convert 私有
  `spill_path()`/`File::create` 删除。
- 三个 SEC-1 攻击测试移到共享模块（含 tag=convert 的 symlink 用例），convert
  spill 测试共用 `SPILL_TEST_LOCK`；新增 Drop 清理测试。
- 门：`TempSpill::create` 两处存在、`fn spill_path()` 清零。

### V25-10（提交 `fd591f2`）

- 残差混合尺度分母；`PcaReport::stream_max_residual` + 超 tol 的
  `tracing::warn!`；`STREAM_DEBUG` 删除，改 `tracing::debug!`。
- R1：`backend::fix_column_signs` 抽出，dense `top_k` 与 `lanczos_top_k` 共用；
  stream↔dense 测试比较特征向量（含符号）。
- R12：删除 `crates/reed-reduce/src/linalg.rs`。
- 测试：1000×8、k=10 秩亏输入在 max_cycles 前收敛；pca e2e 断言 residual
  被记录且 ≤1e-10。
