# v30 最终复核、剩余问题与发布流程

> 复核日期：2026-11-19。对象是当前 `dev` 工作树中的 v30 重构及其未提交改动。
> 本文件是当前状态的权威索引；`01-review.md`、`04-recheck.md`、
> `05-review-followup.md` 保留修复前的证据和修法讨论，不应再单独解释为当前状态。

## 1. 先给结论

### 1.1 IO、filter、并行

截至本次复核，**没有发现新的 P0/P1 运行时 bug**。重点路径已覆盖：

- VCF/BCF、BGZF、索引、scanner、native BCF encoder 和 noodles 对拍；
- filter 的样本子集、INFO/AC/AN、`--fill-tags`、多等位、half-call 和策略切换；
- `--skip-empty`、窗口完整性、线程数确定性；
- pipeline abort、LD idle cancellation、budget/permit 回收和并行输出顺序。

这不等于所有功能都完整，也不等于 API 契约已经最终收口。当前剩余项主要是
发布安全、公共 API 边界、明确的未实现功能和 IO 热路径优化。

### 1.2 v30 关键问题的当前状态

| 项目 | 当前状态 | 证据/门 |
|---|---|---|
| V30-01 PI_HAT 期望 IBS 系数 | **已修复** | `reed-distance/src/matrices.rs` 使用 PLINK 有限样本系数；`run_genome_parity.sh` 对拍 |
| V30-02 `filter -I` 删除 AC/AN | **已修复** | `reed-filter/src/plan/compile/finish.rs` 保留输入 AC/AN；显式 `-x INFO/AC` 仍可删除 |
| V30-03 多等位三路径不一致 | **已修复为拒绝/跳过并计数** | stats/report 跳过并写入 `multiallelic skipped`；e2e 断言 `.afreq/.hardy` 不含多值 ALT |
| N-1 `--skip-empty` 截断窗口 | **已修复** | 跳过空窗口后重新检查完整性；性质测试确认输出等于完整输出去掉空行 |
| N-2 PI_HAT 丢弃单态位点 | **已修复** | `nn` 统计所有 jointly-called 位点，含单态位点 |
| N-3 `.genome` 的 DST 恒为 0 | **已修复** | DST 使用 PLINK 公式并纳入对拍；PPC/RATIO 不伪造，输出 `NA` |
| N-5 `--summary` 缺 PSI/HWE 说明 | **已裁决并文档化** | stats 文档列出已实现节及明确未输出的 PSI/HWE |
| `--fill-tags` 子集顺序回归 | **已修复** | registry order 固定为 `AC,AN,NS,MAF,HWE`，并有 e2e 回归 |
| R01 phylo half-call 口径回归 | **已修复** | phylo 的 `--min-samples-locus` 恢复 vcf2phylip 口径；其他 count 路径保留共享规则 |
| LD idle 取消接线错误 | **已修复** | idle 场景使用同一对 abort handle/signal，并有 drop/abort 测试 |
| R02 AC/AN 重复解析和溢出 | **已下沉** | `reed-core::count::InfoCounts` 统一长度、溢出、`AC > AN` 和回退语义 |
| R03 BCF `pos0` 哨兵 | **已下沉** | `-1` 只表示 telomere/POS 0，其他消费者委托单一转换函数 |
| R05–R07 输入能力 | **已下沉** | `StreamedInput`、`ReplayableInput`、`SeekableInput`、`InputMeta` 和 capability error |
| R09 预算/取消所有权 | **已下沉** | `RunContext { budget, abort }` 作为 driver 的统一输入 |
| R10/R11 临时文件和有序输出所有权 | **已下沉** | `TempSpill`、`OrderedOutput` 负责生命周期和清理 |
| R12 filter 依赖 reed-io | **已修复** | 语义类型下沉到 L0；`cargo tree -p reed-filter -e normal` 不再出现 `reed-io` |
| R13/R14 标签表和发射路径 | **已下沉** | `RECOMPUTE_TAGS` 与 `SurvivorEdit` 成为单一来源 |

### 1.3 当前工程判断

v30 已经脱离了“每个 issue 加一个局部补丁”的阶段：GT/AC/AN/位置哨兵、输入能力、
预算和输出编辑等共同语义已经向底层收敛。当前不能再把它称为 IO/filter 层的
明显屎山，也没有证据表明需要为了修复一个现存 P0 而立刻全面重写。

但是，公共 API 仍有几处设计决策没有冻结，发布脚本也有发布错版本的可能。下面的
问题需要在 v30 后续工作中处理，或由产品决策明确“暂不承诺”。

## 2. 仍存在的真实工程风险

### P1-1：`publish.sh` 发布的是工作树，不是最后一个提交

`ci/publish.sh` 使用临时 Git index，对当前工作树执行 `git add`，随后从工作树生成
公开快照。因此未提交的源码、manifest、CI 和 README 都可能进入 `origin/main`。
脚本不运行完整门，也不验证这些内容已经提交到 `dev`。

这不是运行时 bug，但它会让“已审查版本”和“实际发布版本”不同。发布前必须执行：

```text
dev 完整门通过
→ 提交准备发布的改动
→ 工作树 clean
→ 推送 dev
→ 运行 publish.sh
```

建议后续给脚本增加 clean-worktree 检查，并把源提交哈希写入发布提交信息或文件中。
这属于发布安全增强，不应通过重复运行 main 的完整测试来解决。

### P1-2：`publish.sh` 使用无条件 `--force`

脚本先读取远端 `origin/main` 作为父提交，再执行：

```bash
git push --force origin "$commit:refs/heads/main"
```

在读取远端和推送之间，如果另一个发布者更新了 `main`，当前脚本会覆盖对方的提交。

建议：

1. 优先改为普通 fast-forward push；或者
2. 保留快照历史时使用带旧值的 `--force-with-lease`；并在冲突时失败，不覆盖远端新提交。

### P1-3：`dev` 的 upstream 当前指向 `origin/main`

当前本地配置是：

```text
branch.dev.remote = origin
branch.dev.merge  = refs/heads/main
```

所以不要在当前配置下直接使用裸 `git push`。明确推送 dev：

```bash
git push origin HEAD:refs/heads/dev
```

如果希望以后在 `dev` 上直接运行 `git push`，再单独执行：

```bash
git branch --set-upstream-to=origin/dev dev
```

这是本地 Git 配置问题，不应混入源码提交中。

### P1-4：v30 历史文档存在“旧结论 + 新结论”并列

`STATUS.md` 的早期表格仍把 V30-01、V30-02、N-1 等写成“未实施”，后面的追加章节
又记录了它们已经修复。`04-recheck.md` 和 `05-review-followup.md` 记录的是修复前
复核与第二批修复规范，不是当前实现状态。

本文件和 `STATUS.md` 的最后一节提供当前索引。后续新增审查应只追加“当前状态”，
不要覆盖历史证据，也不要再从旧表格推断当前代码未修复。

## 3. API 和工程设计的剩余建议

以下项目已经由静态 API/工程复核确认，但它们不是本轮发现的默认输入运行时 bug。

### P1-5：CLI 有 distance，facade 没有 `run_distance`

当前 `reed distance` 由 `reed-cli/src/translate/distance.rs` 直接编排输入、矩阵、
bootstrap、建树和输出；facade 没有 `DistanceRequest`、`DistanceReport`、
`run_distance` 或 `run_distance_opened`。这使“每个 CLI 命令都有 library driver”的规则
不完整，也使库用户不能复用 distance 主流程。

**建议：**补齐 typed request/report 和 `run_distance`。在此之前必须由你拍板：

- 推荐：distance 是正式库能力，补齐 facade；
- 备选：distance 永远只是 CLI 能力，并在 API contract 和 library 文档中明确排除。

### P1-6：公开请求类型把有限枚举写成 `String`

`PcaRequest`、`DapcRequest` 的 `stand`、`missing`、`solver`，以及 `XvalRequest.result`
只接受有限字符串集合，却在公开 Rust API 中使用 `String`。非法值只能到运行时才报错，
调用者无法从类型系统得到约束。

**建议：**定义 facade 自有的 `PcaStand`、`MissingPolicy`、`PcaSolver`、
`XvalResultMode`，实现 `Display`/`FromStr` 做 CLI 转换。不要直接把下层 crate 类型
暴露到 facade。

### P1-7：支持 crate 的稳定性承诺与公开签名仍有张力

`reed` 的公开签名包含 `reed_io`、`reed_engine`、`reed_filter`、`reed_stats` 和
`reed_core` 的类型。当前 `check-external-types` allowlist 把这些 workspace 类型全部
放行，同时 API contract 又写这些支持 crate 不纳入 `reed` 的 semver 承诺。

这在 Cargo 发布图上可以工作，但如果支持 crate 改了公开类型，facade 仍会重新编译、
甚至改变公开签名。**建议二选一：**

1. facade 定义稳定的 `Input`、`Header`、`Report`、`Index` 等自有类型；或
2. 明确宣布支持 crate 的公开类型也是 lockstep 稳定面，并为它们加独立 semver/API 门。

推荐先采用第 2 项，避免在 v30 再做一次大规模类型搬迁；等 API 使用者稳定后再逐步
收窄 facade 类型。

### P1-8：`frontend` 虽 `doc(hidden)`，仍能访问整 crate 别名

`reed::frontend::distance`、`graph`、`reduce` 仍可被外部代码调用。`#[doc(hidden)]`
只影响文档显示，不是访问控制；这与“只暴露 facade”在实际可见性上不完全一致。

**建议：**补齐 `run_distance` 后，CLI 只依赖稳定入口，再删除整 crate 别名；同时让
结构门检查缩进后的 `pub use`，避免只检查行首导致漏检。若暂时保留，文档应称其为
明确的 unstable escape hatch，而不是内部不可访问实现。

### P2-1：`run_filter_opened` 缺失，opened 入口不一致

其他主要 driver 有 `run_*_opened`，filter 只有 `run_filter`。需要自行管理输入、重放或
区域能力的库用户不能以同一模式调用 filter。

**建议：**补 `run_filter_opened`，或在稳定 API 文档中明确 filter 不提供 opened 入口。

### P2-2：`reed-distance` 和 `reed-reduce` 未显式继承 workspace lints

它们是发布图中的正常 library，但 manifest 没有 `[lints] workspace = true`。建议补齐，
并让 CI 根据 `publish != false` 自动检查所有发布 crate，而不是手工列出 crate。

### P2-3：feature gate 只写死 `reed-reduce`

当前只有 `reed-reduce` 有 feature，现状可工作；但未来其他 crate 新增 feature 时，
`ci/check-features.sh` 不会自动覆盖。建议用 `cargo metadata` 找出所有带 feature 的
正常 crate，或者把“新增 feature 必须同步修改 gate”作为 CI 的结构规则。

## 4. 明确的功能边界，不要误称为 bug 已修复

### P2-4：`.genome` 不是完整 PLINK 输出

当前保证 `Z0/Z1/Z2/PI_HAT/DST`。`PPC/RATIO` 因 PLINK 内部最后窗口计数无法从当前
全位点统计复现，输出 `NA`；`RT/EZ/PHE` 仍是占位/固定值。当前实现比原来的硬编码 `0`
诚实，但不能宣称与完整 `.genome` 逐列兼容。

建议在 distance 文档和 API contract 明确这是“部分兼容”。是否实现完整列，留给后续
产品决策，不作为本轮阻断项。

### P2-5：`--summary` 不输出 PSI/HWE

这已经在 stats 文档中明确：`--hardy` 提供逐位点 HWE，`--het` 提供逐样本近交系数，
PSI 没有等价输出。它是功能边界，不是隐藏错误；不要把 `--summary` 描述成完整
`bcftools stats -s -` 替代品。

### P2-6：`--hwe P midp` 可解析但随后拒绝

帮助语法仍接受可选 `midp`，执行时才返回“not implemented yet”。建议二选一：隐藏该
语法，或实现 mid-p；在实现前不应让帮助看起来像已支持。

### P2-7：site-count windows 尚未实现

`WindowGrid::Sites` 在 driver 和 argument rendering 阶段都明确返回 usage error。文档应
继续标记为未实现，不要列入已支持窗口类型。

## 5. IO 优化建议：重复 header 解析

当前 `PendingInput::open` 会完成 sniff、header facts 扫描；`into_opened` 随后再次调用
`build_source`，source 重新扫描和解析 header。当前只打开一次文件，结果正确，也没有
发现新的 IO bug，但大 header 会重复 inflate、复制和 noodles 解析。

建议下一个 IO 小版本做单一所有权：

1. source 接受已经解析的 `HeaderFacts`；或
2. `PendingInput` 保存首块/前缀，让 source 复用第一次扫描结果；
3. 保持现有“文件只 open 一次”和 indexed region capability 检查；
4. 用 header 很大的 VCF/BCF 做 CPU、分配和 RSS 基准后再改。

这是优化，不是当前 P0/P1 bug；不要为了消除一次解析而破坏 BGZF/region 边界。

## 6. noodles 决策建议

当前 `reed-io` 仍把 noodles 作为内部 IO 后端，`reed` facade 的 `VcfHeader` 已经将
`noodles::vcf::Header` 包在私有字段中，因此 noodles 没有直接成为 facade 的第三方公开
类型。当前没有发现一个需要立即全面移除 noodles 的 correctness bug。

**建议不在 v30 一次性全面重写。** 推荐分阶段：

1. 保持 noodles 作为 dev differential oracle，固定版本并锁定关键行为；
2. 生产路径优先保留当前已经自写且有测试的 VCF line scanner、BCF primitive、编辑器和
   encoder；
3. 如果要减少生产依赖，先替换一个边界：例如 header/VCF reader 或 CSI/TBI reader，
   每次都用 noodles 对拍、gold fixture、错误注入和性能基准；
4. 迁移完成后再从 `reed-io` 和测试依赖中删除 noodles，而不是先删依赖再补解析器。

删除 noodles 的收益是版本和维护自主权；代价是重新承担 VCF/BCF/BGZF/CSI/TBI 的格式
兼容、异常输入和安全修复。当前 v30 应把“彻底自写”作为待拍板路线，而不是默认重构目标。

## 7. `publish.sh` 到底推送什么

结论明确：**当前 `ci/publish.sh` 只推送 `origin/main`，不会推送 `dev`。**

它会：

- 用临时 index 从当前工作树收集公开路径；
- 排除 `docs/`、`docsite/`、`gold/`、各 crate 的 `tests/` 和 Pages workflow；
- 以远端 `origin/main` 为父提交生成公开快照；
- `git push` 到 `refs/heads/main`。

它不会：

- 提交本地 `dev`；
- 推送 `origin/dev`；
- 运行 dev 完整门；
- 检查工作树 clean；
- 检查当前分支是否是 `dev`；
- 发布 GitHub Release 或二进制。

当前推荐命令：

```bash
# 1. 在 dev 完成完整门并提交
git add -A
git commit -m "refactor(v30): consolidate semantic ownership and API release gates"

# 2. 显式推送 dev（当前 upstream 配置不正确）
git push origin HEAD:refs/heads/dev

# 3. 由当前 clean 工作树生成并推送公开 main
bash ci/publish.sh "publish: v0.5.0 public source"
```

随后应从 `origin/main` 生成的公开提交创建匹配 `Cargo.toml` 的版本 tag，再触发
release workflow。不能从 dev 提前打 tag，否则 release workflow 的 HEAD/tag 校验会拒绝。

这不是让 main 再重复一遍 dev 的完整门：dev 完整门是发布前置条件，main 只验证公开源代码
能够构建。公开 main 是否保留 `docsite/` 和文档部署 workflow，仍是产品选择；当前脚本
排除它们，因此 README 不应把 `docsite/` 相对路径当成 main 一定存在的内容。

## 8. 需要用户拍板的事项

| 事项 | 我的建议 | 备选 |
|---|---|---|
| distance 是否进入 facade | **进入**，补 `DistanceRequest/Report/run_distance` | 明确 CLI-only，并删掉“每个命令都有 driver”的承诺 |
| 支持 crate 的稳定面 | **0.5.x lockstep 公开承诺**，加 API/semver 门 | facade 自有类型隔离，后续做较大类型迁移 |
| `frontend` 整 crate 别名 | **补齐 facade 后删除** | 保留并明确 unstable |
| `run_filter_opened` | **补齐** | 明确 filter 不支持 opened |
| `publish.sh` | **clean 检查 + `--force-with-lease`** | 保持脚本简单，由发布者手工保证 |
| main 的 docsite | **继续源码-only，README 指向外部文档站** | 把 docsite 和 Pages workflow 纳入公开树 |
| noodles | **先保留内部后端/开发 oracle，分阶段替换** | v30 立即全面自写 |
| `.genome` 完整列 | 本轮不阻断，明确部分兼容 | 实现 PPC/RATIO/RT/EZ/PHE |
| PSI/HWE summary | 本轮不阻断，保持现有替代命令 | 扩充 summary 的独立扫描 |

## 9. 当前验证和未覆盖范围

已有验证记录包括：

- stable debug/release workspace tests：通过；
- clippy、fmt、doctest、rustdoc `-D warnings`：通过；
- API snapshot、nameability、external-types、feature、semver、package、docs、structure、
  doc-drift：通过；
- genome parity、convert parity、IO/filter/engine/LD 相关回归：通过。

本次没有从 GitHub API 读取分支保护设置，不能替你确认远端是否允许 force push；也没有把
`run_distance`、完整 `.genome` 或 noodles 全面迁移作为已完成事项。它们在上文分别列为
设计决策或后续工作。
