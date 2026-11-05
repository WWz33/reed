# 05 CLI 参数层（V25-06）

> 历史上 `-o`/`-O`/`-W`/`--memory`/`--threads` 的每个 bug 修复都要落在 3–7 个
> translate 文件里——修复点数量 = 声明点数量，这就是"补丁堆补丁"的机制本身。
> 正面对照就在同一仓库：绑成 clap `value_parser` 的参数（`--geno`、`--window`、
> `--verbosity`）从不反复。本轮把反复出事的参数全部类型化。

---

## 1. 现状

### 1.1 `--threads 0` 语义随命令翻转【据】

- `translate/distance.rs:102-105`：用户值直传
  `rayon::ThreadPoolBuilder::num_threads(threads)`——rayon 语义 0 = 默认 =
  **全部逻辑核**；
- filter / stats / popg / pca 走 `ThreadBudget`：`budget.rs:62`
  `is_sequential() = threads <= 2` → **顺序单线程**；
- 后果：共享集群上显式 `--threads 0` 限线程的用户，filter 单核、distance 静默
  占满节点——越界方向是"多用资源"，最危险的方向；`cli_contract.rs` 无任何
  `--threads 0` 用例兜底。

### 1.2 预算构造 7 份 + panic 风险【述】

`threads.unwrap_or_else(available_parallelism) + MemBudget::from_bytes(m.saturating_mul(1<<20))`
出现在：`translate/common.rs:417-427`、`pca.rs:38-50`、`distance.rs:82-94`、
`convert.rs:158-170` 与 `:387-399`、`rename.rs:92-104`；`filter.rs:385-389` 是
第 8 种形状（`mem_bytes`）。

- `available_parallelism().expect("invariant: parallelism")` ×6 —— panic =
  exit 101，破坏 1/2 退出码契约；`budget.rs:27` 自己就是 `map_or(1, …)` 的范本；
- `--memory 0` 被接受，运行到第一个 chunk 才报 "exceeds memory budget"
  （runtime exit 1，`pipeline.rs:412`【述】），应为 usage；
- `--memory 4194304`（×2²⁰ 饱和到 `u64::MAX`）静默变成无限预算。

### 1.3 `-O` 解析三份【述】

`translate/filter.rs:393-431`、`rename.rs:46-91`、`convert.rs:55-79` 各自实现
"拆字母 + level、level 1–9、默认 6、扩展名 .gz/.bgz 推断"。R4（`-O ''` panic）
的修复注释同时出现在三个文件——**同一个修复做了三遍**。`-W` 的 BGZF 判定另有
两条"恰好等价"的实现（`filter.rs:478-487` 与 rename 的 `compressed`）。

### 1.4 `--ld-decay` 双重解析【述】

`cli/mod.rs:686-688` 用 `parse_optional_size` 校验却只存原始 `String`；
`popg.rs:627` 再 `parse_size` 一次。`popg.rs:624-626` 的 `None => None` 分支
不可达（:615 已提前返回）；`--ld-decay 0` → r² 窗口 0bp → 空输出 exit 0。

---

## 2. 规范【规】：`params` 模块

reed-cli 新建 `params` 模块，每个反复出事的参数一个 newtype + clap
`value_parser`，校验只写一次：

| 类型 | 校验规则 | 错误（usage / exit 2） |
|---|---|---|
| `Threads(u32)` | ≥ 1 | `--threads: must be ≥ 1 (0 does not mean "all cores")` |
| `MemMiB(u32)` | ≥ 1；解析期拒绝溢出（不用 saturating_mul 吞掉） | `--memory: must be ≥ 1 MiB` |
| `OutputPrefix` | 非空、非 `-` | 沿 `require_output_prefix` 文案（V23-21 已有） |
| `OutputType { kind, level }` | kind ∈ `v/z/b/u`；level 1–9，默认 6；扩展名推断在此一处 | 沿 R4 文案，一处 |
| `WriteIndexFlag` | `-W` / `-W=tbi` | 沿现文案 |

```rust
// reed-engine：唯一构造器，7 处调用一行化
impl ThreadBudget {
    pub fn from_cli(threads: Option<Threads>, mem: Option<MemMiB>) -> ThreadBudget;
}
```

配套规则：

1. `available_parallelism()` 失败 → 1（与 `budget.rs:27` 一致），CLI 层不再
   `expect`；
2. `-O` 的解析 / 校验 / 扩展名推断收敛为 `OutputType::from_flags(out, spec)`，
   CLI 与 writer 共用——`-W` 校验的 "Mirror the writer's rule" 注释消失；
3. `--ld-decay` 的 value_parser 直接产出 `Option<u64>`（与 `--window` 一致），
   删除 popg 的二次解析与不可达分支；`0` 拒绝（usage），不静默空输出；
4. 【荐】`translate/distance.rs` 的 throwaway pool 换
   `reed_engine::budget::run_auxiliary`（`budget.rs:81-105` 已有），自建池与
   策略池两条线程路径合一；
5. 【荐】"每类参数一个 owner"写进 02-cli 的规范页：数值域在 value_parser、
   语义域在 plan compile、跨命令契约在 `require_*` 家族——新 flag 落地前先问
   owner 是谁。

---

## 3. 迁移表

| 现状点 | 迁移目标 |
|---|---|
| 7 处 budget 构造（§1.2 列全） | `ThreadBudget::from_cli` |
| `-O` 三份（filter/rename/convert） | `OutputType::from_flags` |
| `-W` 镜像校验两份 | 共用 `OutputType` 的 compressed 判定 |
| `--ld-decay` 双解析 + 死分支 | value_parser 产出 `Option<u64>` |
| `--memory` saturating 溢出 | `MemMiB` 解析期拒绝 |

## 4. 验收

usage 矩阵测试（每个子命令 × 每个非法值，断言 exit 2 与文案）：

```
--threads 0        --memory 0        --memory 4194304（溢出界）
-O ''              -O z0 / -O z10   -W 与 -o - 组合
--ld-decay 0
```

- budget 构造点 grep 7→1；`available_parallelism` 的 `.expect` 清零；
- 常驻门（`ci/check-structure.sh`）：`available_parallelism` 只允许出现在
  `reed-engine/src/budget.rs` 与 `reed-cli/src/params.rs`；`-O` 的拆字母/level
  解析只允许出现在 `OutputType::from_flags`；`--ld-decay` 的 `parse_size` 只允许
  出现在 value_parser（00 P3）；
- 既有 round-trip / cli_contract 测试全绿；
- 行为对照表（变更前后）：`--threads 0` 从"分叉"变"统一 usage"；
  `--memory 0` 从 runtime exit 1 变 usage exit 2——属"静默错误 → 显式报错"，
  不违反本轮行为边界。

### 实施记录（2026-11-04，提交 `3c2deee`）

- **newtype + value_parser**：`reed_engine::budget::{Threads, MemMiB}`（`FromStr`
  一份、`Display`、`parse_threads`/`parse_mem_mib` 供 clap）；8 个 CLI 的
  `threads`/`memory` 字段改为新类型，非法值在 clap 层 exit 2。
- **`ThreadBudget::from_cli`**：7 处预算构造（common/pca/distance/convert×2/
  rename/filter）收成一行；`index` 的裸 `available_parallelism` 也改走
  `from_cli(threads, None).threads`。CLI 内 `available_parallelism` 清零。
- **`-O`**：`reed-cli/src/params.rs::parse_output_type`（v/z/b/u + level 1-9）
  与 `parse_text_output_type`（rename / hapmap→VCF）；`split_output_type`
  删除，R4 的三份拷贝变一处。
- **`--ld-decay`**：clap value_parser `parse_ld_decay` 直接产出
  `Option<u64>`（裸 flag=100 kb，`0` 拒绝）；`build_ld_request` 的二次解析与
  不可达 `None` 分支删除。
- **`ci/check-structure.sh` 新建**：四条常驻门（CLI 无 available_parallelism、
  `-O` 解析只在 params、`--ld-decay` 不再二次 parse、GT 语法只在 gt_decode）。
- **偏差**：distance 的自建 rayon pool 未换 `run_auxiliary`（【荐】项）——它把
  pool 传进多个 estimator/bootstrap 闭包，改成 `run_auxiliary` 需要重排闭包
  边界；本轮只保证 `--threads` 域统一，池的统一留给 v26。
- 验收：usage 矩阵 17 例（7 个子命令 × 非法值）全 exit 2；988 passed / 0
  failed；clippy 0；doc-drift 97。
