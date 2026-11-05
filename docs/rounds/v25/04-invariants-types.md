# 04 不变式与契约类型（V25-05 / V25-07；附 R12）

> 形状一致、mask 长度、计数口径——这些不变式目前靠调用点手写 assert 和
> 人肉同步注释维持。本轮把它们变成复用单元：`merge_zip_eq`、`from_codes`、
> `KeepMask`。V23-28（merge 截断）与 S68（短 mask）是同类事故的两次实证。

---

## 1. V25-05：merge_zip_eq

### 现状【据】

V23-28 发现 `debug_assert!` 在 release 是空操作、`zip` 静默截断较长一侧导致丢
数据，修复落在 3 处（`StatAcc::merge`、`ContigBins::merge_bins`、
`BpBlockJackknife::merge`）。**同型位点仍开放**——reed-stats 15 处
`iter_mut().zip`：

```
het.rs:85/88/91      ld.rs:206/209       matrix.rs:241
fst.rs:488           sfs.rs:73/291/479   theta.rs:88/91
window.rs:108/137/164
```

- `sfs.rs:289/479`：仅 `debug_assert_eq`（release 无操作）；
- `het.rs:85`、`sfs.rs:73`：无任何保护；
- `theta.rs:83-93`【述】：先 `resize` 再 zip——resize 掩盖调用方 bug 而不是
  暴露它。

任何调用方 bug（不同 `n_samples` 的 partial、不同 LD 边界表的 worker）在 release
下都会**静默丢弃尾部并输出错误数字**——V23-28 已实证过这个失效模式不是理论风险。

### 规范【规】

```rust
/// 形状不一致即 panic（release 同样检查），消息带模块上下文。
/// 把"两侧形状必须一致"从调用点纪律变成复用单元。
fn merge_zip_eq<T>(ctx: &'static str, dst: &mut [T], src: &[T], add: impl Fn(&mut T, &T));
```

- 15 处全部收敛到它；`theta.rs` 的 resize-then-zip 改为显式错误或并入；
- 【荐】同批统一"merge 律"的表述：形状一致是**合并的前提**，在 trait doc 里
  写一次，各累加器不再各自注释。

### 测试

- 单测：形状不匹配 → panic，消息含 `ctx` 定位；
- 既有合并律 proptest 全绿（形状一致路径行为不变）。

### 验收

- `grep 'iter_mut().zip' crates/reed-stats/src/` 仅剩 `merge_zip_eq` 内部一处。

---

## 2. V25-07a：计数唯一入口

### 现状

- `count_site(&[i8], Option<&[bool]>)` **生产死代码**【据】：唯一调用者在冻结的
  oracle（`oracle/src/estimators/site_counts.rs` 是独立实现，非调用）；
- 计数语义（half-call 计一个 allele、不进基因型类、`-g` 读作 missing）5 份镜像
  【述】：

| # | 位置 | 说明 |
|---|---|---|
| ① | `reed-core/src/count/site_counts.rs:124` | i8 路径（死代码本体） |
| ② | `reed-core/src/count/mod.rs:56` `SiteGroupCounts::add_code` | GtCode 路径 |
| ③ | `reed/src/filter/record_pipeline.rs:600-630` | 注释自称 "`count_site` semantics" |
| ④ | `reed/src/filter/mid_stage.rs:92-135` | 位平面 word kernel |
| ⑤ | `reed/src/stats/writers/reports.rs:568-581` | 注释自称 "The same per-code tally `SiteGroupCounts::add_code` runs" |

- "identically by construction" 只靠人肉同步注释——半 call 语义再调一次要改 5 处；
- ① 的"短 mask 默认保留"（`keep.get(i).copied().unwrap_or(true)`，`site_counts.rs:127`
  【述】）与 `SampleMask` 的文档契约（"no implicit shorter-mask default"，A13 修复）
  直接矛盾——S68 的根源 API 还活着。

### 规范【规】

```rust
impl SiteCounts {
    /// codes + keep-mask → counts，唯一语义点；mask 是 SampleMask（长度强制）。
    pub fn from_codes(codes: &[GtCode], keep: &SampleMask, masked: Option<&Bitset>) -> Self;
}
```

- ①删除（oracle 的独立实现不动——它是刻意的 parity 参照）；
- ③ ⑤ 改调用；② 的 match 从同一函数派生；
- ④ word kernel 保留（性能），但**必须有与 `from_codes` 的等价性测试**（随机
  codes+mask 的 proptest——当前没有，这是 5 份镜像里唯一没有交叉验证的一对）；
- "短 mask 默认保留"随 ① 的删除消失。

### 测试

- `from_codes` × word kernel 等价性 proptest；
- 既有 e2e 计数（filter funnel / stats reports）数值不变。

### 验收

- 计数语义实现点 5→1（word kernel 作为唯一豁免，带等价性测试）；
- `grep -rn 'unwrap_or(true)' crates/ --include='*.rs' | grep -v test` 无 mask
  默认值残留；
- counter 清单三处镜像（core `PolicyCounters` / io `pack_triggers`+`apply_triggers`
  / filter `funnel::PolicyCounters`）收成 core 单一定义；新增 policy 触发器只改
  一个 crate（00 P4，见 §4 附）。

---

## 3. V25-07b：GroupTable::restrict 与 KeepMask

### 现状

- `group.rs:177-181`【据】：`group.mask.and(keep)` 失败时 `unwrap_or_else` 保留
  原 mask，函数仍返回"成功"的表：

  ```rust
  let mask = group.mask.and(keep).unwrap_or_else(|_| {
      // Table and keep mask disagree on the sample count: keep
      // the original group untouched rather than fabricate one.
      group.mask.clone()
  });
  ```

  10 样本的表 × 9 样本的 keep mask → 每组 mask 原样保留、`n_samples` 仍是 10 →
  **分组统计（π/Fst 分母）在未过滤 cohort 上算完，零报错**；
- 这与 S68（doc_v20：短 mask 让两个 strategy 数出不同样本集）同类；当时的修复
  是在 mid_stage 加 `assert!`（现象级），根源 API 未动；
- keep-mask 以裸 `&[bool]` / `Arc<[bool]>` 层间传递，"长度必须等于 header 样本数"
  的检查 + 错误文案 5 处【述】：`stats/driver.rs:43-48`、`convert/job.rs:1703-1708`、
  `pca/mod.rs:118`、`filter/mod.rs:888-912`、`filter/mid_stage.rs:210-217`。

### 规范【规】

1. `restrict` 长度不匹配返回 `CoreError::LengthMismatch`（签名改
   `Result<GroupTable, CoreError>`；`#[must_use]` 语义自然保留）；
2. `KeepMask` 新类型：`KeepMask::build(header, sample_spec) -> Result<KeepMask, _>`
   ——构造时校验长度与名称，内部持 `Arc<[bool]>`；五个消费点收口，mid_stage 的
   `assert!` 自然消失；
3. 【荐】CLI 侧 `parse_sample_keep`（`translate/common.rs:256-323`【述】）与
   filter 后端 `keep_list_from_names/file` 同批合并——V23-40 只统一了 CLI 两个
   翻译点，filter 后端那份没并进来。

### 测试

- restrict 不匹配单测：返回错误（而非旧分组），错误带两个长度值；
- S68 既有用例保持；`KeepMask` 构造错误用例（未知样本名 / force-samples 警告）。

### 验收

- `unwrap_or_else` 保底分支删除（grep `keep the original group untouched` 清零）；
- mask 长度检查点 5→1（都在 `KeepMask::build`）。

---

## 4. V25-07c：幽灵计数器【决】

### 现状【据】

`n_haploid_as_missing` 的全链路接线：

- `reed-core/src/genotype/policy.rs:67`：字段声明；
- `reed-io/src/bcf/gt.rs:198`（`pack_triggers`：计数器 > 0 → 置位）与
  `:218`（`apply_triggers`：置位 → 计数器 `+= 1`）；
- `reed-filter/src/funnel.rs:183/193`：报告列 `haploid_as_missing`。

但 `TRIG_HAPLOID_MISSING` 全仓只有三处：定义（`gt.rs:66`）、打包（`:199`，从
计数器读）、解包（`:217`，写回计数器）——**自成环、无外部生产者**。
`PloidyPolicy` 只有 `RequireDiploid / HaploidAsHalfCall / HaploidAsHomozygous`
三变体【据】，不存在会把 haploid 记为 missing 的策略。

后果：报告的 `haploid_as_missing` 列永远为 0；文档描述的 "Haploid genotypes
zeroed to missing" 策略不存在；用户按报告排查 haploid 去向被误导。

### 选项

- **A（荐）**：删除字段 + 位 + funnel 列。恒零且误导的输出没有存在理由；
  删除后 counter 清单三处镜像（core / io 打包 / funnel）同步少一项；
- B：补 `PloidyPolicy::HaploidAsMissing` 变体与 `from_alleles` 的计数——
  但无用户需求（YAGNI），且要动 policy 枚举的默认值语义与既有 e2e 的计数快照。

【决】需要你拍板 A/B（README 决策表同项）。

### 附：counter 清单的三处镜像【述】

core `PolicyCounters`（5 字段）/ io `pack_triggers`+`apply_triggers`（5 位，
`gt.rs:190-226`）/ filter `funnel::PolicyCounters`（5 字段，`funnel.rs:160-199`）。
无论 A/B，【规】core 提供 `(name, value)` 迭代或打包宏的唯一定义，io/filter
消费它——每加一个 policy 触发器从改 3 个 crate 变成改 1 个（V25-07 验收项，
00 P4）。

---

## 5. R12：死代码清理（并入 V25-10 收尾）

### 现状【据】

`crates/reed-reduce/src/linalg.rs` 仍有 198 行；`reed-reduce/src/lib.rs` 没有
`mod linalg`，全仓引用的是 `reed_core::linalg`（`backend.rs`、`grm_reduce.rs`、
`pca.rs` 等）。文件是 V24-03 清理后的残留：编译器不编译它，但它会继续吸引
"两份实现同步修"的补丁（P5）。

### 规范【规】

删除 `crates/reed-reduce/src/linalg.rs`；不迁移任何代码（`reed_core::linalg`
是唯一实现）。若 grep 发现隐藏引用，改为 `reed_core::linalg` 的 re-export。

### 验收

- `test ! -e crates/reed-reduce/src/linalg.rs`；
- `grep -rn 'mod linalg\|crate::linalg' crates/reed-reduce/src` 为 0；
- `cargo test -p reed-reduce` 全绿。

### V25-05 实施记录（2026-11-04，提交 `de05e12`）

- `reed_stats::merge_zip_eq(ctx, dst, src, add)`：`assert_eq!`（release 同样
  检查）+ ctx 定位；`merge.rs` 单测覆盖形状一致与 `catch_unwind` 的 panic 消息。
- 收口 **23** 个 `iter_mut().zip`（规范列的 15 处 + `summary/report.rs` 6 +
  `summary/histogram.rs` 1 + `window.rs` 的嵌套对），`grep 'iter_mut().zip'
  crates/reed-stats/src/` 只剩 `merge.rs` 实现内 1 处（文档措辞改为 "zip over
  iter_mut()"，避免命中验收 grep）。
- `ThetaAcc::merge` 的 resize-then-zip 删除——形状不一致现在直接 panic；
  `BpWindow` 的 `should_panic` 期望消息改为 `BpWindow bin accumulators`。
- 985 passed / 0 failed；clippy 0。

### V25-07 实施记录（2026-11-04，提交 `11b9b18`）

- **计数（a）**：`GtCode::class()`（HomRef/Het/HomAlt/HalfCall/Missing）+
  `n_ref`/`n_alt` 成为唯一分类；`SiteCounts::from_codes(_kept)` 消费它，
  `SiteGroupCounts::add_code`、reports 的 per-group tally、filter 计数段全部
  从同一 class 派生。`count_site(&[i8])` 与其短 mask 默认删除（oracle 独立
  参照保留）。kernel × from_codes 的 proptest 新增（随机基因型/keep/mask）。
- **mask（b）**：`GroupTable::restrict -> Result<_, CoreError>`（不再静默保留
  未过滤分组）；`KeepMask::try_new` 单点校验长度，stats/convert/pca/filter 四个
  驱动改用它；`pruned_header`/`SampleLayout::from_assignments`/计数核心的短
  mask 默认全部删除（短 mask 现在 panic/报错）。
- **counter（c）**：幽灵 `n_haploid_as_missing` 删除；`PolicyCounters::FIELDS`
  是触发器唯一表——io 的 bit packing / apply 与 filter 报告都从它派生，
  `funnel::PolicyCounters` 直接复用 core 类型（serde 字段名保持不变）；
  reports 的 `stats_policies` 手抄改为调用。
- **验收 grep**：`fn count_site` 仅 oracle；`keep.get(..).copied().unwrap_or(true)`
  0 处；`n_haploid_as_missing`/`TRIG_HAPLOID_MISSING` 0 处；`check-structure.sh`
  新增 3 条门（共 7 条全绿）。
- 991 passed / 0 failed；clippy 0；doc-drift 97。
