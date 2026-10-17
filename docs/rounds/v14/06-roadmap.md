# 06 实施路线（C0–C5）

> 原则：**先钉真值（0.3.x 输出全量固化），再纯搬运拆分（不加功能，
> 对拍全等），然后 help/渲染（零语义），最后两个新面（`--hist`、
> `rename`）各自垂直落地**。每步满足 §1 通用门。与 v13 相同：本轮的
> 绿门主要是**等价对拍**（拆分步）与**参照对拍**（hist/rename 步）。
>
> 不重开 v6–v13 的任何结构项；filter/convert 语义零触碰；QC 六 flag
> 的输出契约（拼写/布局/字节）是全轮的最高优先不变量。

---

## 1. 通用验收门（每步）

1. **双 profile 全绿**：`cargo test --workspace`（debug 与 release）；
   `cargo clippy --release --workspace --all-targets -- -D warnings`
   不新增警告。
2. **确定性**：本步新增/搬迁的每个输出，threads ∈ {1, 4} 下 sha256
   相同。
3. **既有面不回归**：v13 及以前的 gold 对拍全绿——
   `run_convert_parity.sh`、`run_easysfs.sh`、`run_pixy.sh`、
   `run_pixy_bed.sh`、`run_vcftools.sh`、`run_vcftools_fst.sh`、
   `run_bcftools_af.sh`、`run_plink2_ld.sh`（C2 起命令词 stats→popg
   的替换在脚本内完成，diff 语义不变）。
4. **退出码与消息**：usage 错 2 / 运行期错 1；重定向文本按 02 §2.4。
5. **性能**：每步跑 `bench/suites/e2e/results.tsv` 追加一行；退化门
   见 01 §8。
6. **文档同步**：落地一步即更新本目录 STATUS 与
   `ci/check-doc-drift.sh` 对应守卫（01 §6 的 1–14 按步启用）；
   `docs/spec/cli.md` 随 C5 收口。

## 2. 步骤

### C0 真值与基线（不写产品代码）

| 项 | 内容 | 门 |
|---|---|---|
| C0-1 | 0.3.x `reed stats` 输出真值固化：`gold/popgen/run_split_parity.sh` 骨架——QC 六 flag、窗口统计代表组合（π/dxy/Fst wc84+hudson/θ/Tajima/PBS/SFS/LD decay × window/levels/bed），跑 0.3.x 二进制产出期望树入库（fixture 档 + subset 档） | 期望树齐；脚本此时对 0.4.0 失败即红（未来门） |
| C0-2 | 基线墙钟：QC 全报告、QC+拟新增 hist 负载（手工 bcftools stats 对照）、popg 全统计、`bcftools annotate --rename-chrs --set-id` + `reheader -s` 三进程，clean.vcf.gz 记入 STATUS | STATUS 表四行 |
| C0-3 | bcftools 参照真值：`bcftools stats` 的 QUAL/DP/IDD 直方图节（fixture 上）与 `reheader -s`/`annotate --rename-chrs/--set-id` 的输出（含头行重排差异登记进 B14 表） | 参照文件 + 差异注记齐 |

### C1 拆分：`stats` / `popg`（纯搬运，零语义）

| 项 | 内容 | 规范 |
|---|---|---|
| C1-a | `PophCli`（搬迁表的全部 flag 原拼写）+ `Command::Poph` + removed-flag 重定向（旧 stats 上的群体/窗口/清单 flag 解析期报错 + 完整改写命令） | 02 §2 |
| C1-b | translate 拆分：`translate/popg.rs`（自 `translate/stats.rs` 提取，`build_sfs_request` 共享不动）；库门面改名 `run_stats→run_popgen`、`run_reports→run_stats`（02 §6）+ `to_args` 往返更新 | 02 §2/§6 |
| C1-c | `run_split_parity.sh` 转绿：QC 六报告与窗口统计全部逐字节等于 C0-1 期望树 | 01 §7 |
| C1-d | 既有 gold 脚本命令词替换（stats→popg）后全绿；`--indep-pairwise` 在 popg 的报错/行为文本核对 | 02 §2.3 |

### C2 分组 help（零语义）

| 项 | 内容 | 规范 |
|---|---|---|
| C2-a | `GROUPS` 常量 + help 模板 + 渲染函数（about 取 clap 元数据）；`popg` alias `popgen`、`rename` 前置占位（变体可先以 `#[arg(hide = true)]` 不渲染——实现期定：或 C2 只渲染既有五动词 + help） | 02 §1.3 |
| C2-b | 守卫：变体恰好一组测试 + `--help` 快照契约（五组标题）+ README 命令表同源 | 02 §1.4、01 §6 的 1/2/13 |

### C3 `--hist`（stats 垂直切片 → 通用标签）

| 项 | 内容 | 规范 |
|---|---|---|
| C3-a | `HistPlan` + 内建 metric 渲染先行：`qual/dp/dp-site/idl/af`（复用 `QcReportAcc` 已积累的直方图，纯渲染层）+ TSV 布局 + 精确分位行；`--hist` 单独可路由 | 02 §3.1–3.2 |
| C3-b | GQ 解码 + 通用 `INFO/TAG`/`FMT/TAG`（GenotypeScanner 字段计划扩展 + 累加器定律测试） | 02 §3.3–3.4 |
| C3-c | 对拍门：`run_hist.sh`（bcftools QUAL/DP/IDD 节一致；af 按 init_iaf）+ threads {1,4} sha256 + bench 行 | 01 §7/§8 |
| C3-d | stats 的 `--pop` 分层 QC（分层 missingness/freq 行/列）+ `--singletons`（C14-3 维持进轮；vcftools 对拍 + POP 列接 --pop） | 02 §2.2 |

### C4 `rename`（新子命令垂直切片 → 全平面）

| 项 | 内容 | 规范 |
|---|---|---|
| C4-a | `RenameCli` + 门规则（§4.3 校验表）+ `translate/rename.rs` + `Command::Rename` 分发 + 分组渲染接入（Transform 组） | 02 §4.2–4.3 |
| C4-b | VCF 字节编辑路径：样本平面（`-s/-N/--replace-sample-name`，含快路径）→ 染色体平面（`--rename-chrs/--chr-style`，头+记录双改）→ ID 平面（五个 flag 全序）→ header 平面（`-h/-f/--add-meta/--remove-meta`） | 02 §4.2–4.3 |
| C4-b' | BCF 路径：头块重写 + 染色体字典重映射 + 纯样本改名快路径（B14-7） | 02 §4.4 |
| C4-c | 对拍门：`run_rename.sh`（bcftools reheader/annotate 等价集比较）+ VCF↔BCF 路径一致 e2e + threads {1,4} sha256 + provenance 头测试 | 01 §7 |
| C4-d | 索引失效 warning（C14-6）+ `-W` 建索引 + bench 行（对标 bcftools 三进程 ≤ 1/3） | 02 §4.4、01 §8 |

### C5 收口

| 项 | 内容 | 规范 |
|---|---|---|
| C5-a | 库 API 全量改名复核（`run_popgen`/`run_stats`/`run_rename`/`HistRequest`）+ README Library use 改写进编译门 | 02 §6 |
| C5-b | 文档收口：cli.md 拆节 + B14/C14 登记、README 双语、CHANGELOG 0.4.0（逐 flag 搬迁清单）、architecture §13.14、statistics-map 补录；drift 守卫 1–14 全绿 | 02 §9 |
| C5-c | 全量确定性 + 性能终门（01 §8 表填满）+ STATUS 收口 | §1 |

## 3. 顺序与依赖

```
C0（真值/基线） → C1（拆分+重定向，纯搬运） → C2（分组 help，零语义）
              → C3（--hist） → C4（rename） → C5（收口）
```

- C1 是唯一"动既有命令面"的步骤：风险集中在 removed-flag 重定向的
  完整性（搬迁表逐 flag 核对）；等价门（C1-c）是硬门，不过不进 C2。
- C2 独立于 C3/C4，可并行或提前（但组名单里 `rename` 的出现依赖 C4-a
  的变体存在——C2 先行时 `rename` 在名单中占位、变体以最小骨架进树，
  或 C2 晚于 C4-a，实现期二选一，规范不锁）。
- C3 与 C4 正交（stats 侧 vs 新子命令），排期紧张时 C3-d 的
  `--singletons`/分层 QC 可独立后补（C14-3 已留口）。
- **最小可发布面**：C0+C1+C2+C5-b 的文档部分（拆分 + 分组 help 即
  命令面重组的主体）；C3/C4 各自独立成发布面。

## 4. 本轮不做什么

| 项 | 原因 |
|---|---|
| `pca` / `relatedness` | v15 候选；网格预留位（02 §1.2），矩阵输出模型须独立规范轮 |
| `norm` / `concat` / `isec` / `merge` / `query` | v15+ 候选；各需独立语义台账（参照源码通读） |
| `genautocomplete` | clap_complete 半天工作量，随任一步顺手加即可，不占步骤 |
| filter 的 `--extract/--write-snplist/--set-all-var-ids`（X6 遗留） | ID 模板由 `rename --set-id` 承接；`--extract` 位点清单提取留 filter 轮 |
| 撤回 stats `--sfs`（convert 双入口之外的第三入口） | 由重定向自然消亡（`--sfs` 在 stats 报错引导 popg/convert），无需删除代码 |
| `rename --order-samples` | C14-4 默认不进（merge/concat 轮前置） |
| `--hist` 参数化 bin 边界（`--hist-bins`） | C14-5 默认固定边界；参数化留 backlog |
| stats/popg 语义调整（估计器、缺失语义、窗口定义） | 本轮是命令面重组，语义冻结（D14-1 的推论） |
| 嵌套子命令（`reed popgen stats` 形态） | D14-8：破坏 cli.md §0.9 契约 |

## 5. 风险与开放问题

| 风险 | 缓解 |
|---|---|
| removed-flag 重定向漏 flag（搬迁表不完整） | 搬迁表以 `StatsCli` 结构体逐字段核对（01 §3.1）；守卫 3 + 快照测试锁 `popg --help` 与 `stats --help` 的 flag 集差集 |
| `translate/stats.rs` 提取时行为漂移（1108 行） | C1 只做"移动 + 改名"，diff 审查以移动为主（v12 D12-15 先例）；等价门逐字节兜底 |
| `QcReportAcc` 与 `HistRequest` 的耦合（hist 独立可路由但 summary 累加器是 summary 的） | 实现期二选一：hist 计划挂进 `QcReportAcc`（零重复解码）或并立计划（`--hist` 不带 summary 时跳过 summary-only 节）——02 §6 留口，规范只锁"单次扫描双消费" |
| GQ 通用解码对无 GQ 文件的性能影响 | `HistPlan` 门控：不请求 `gq`/FMT 标签不解码（`--depth` 的 DP 先例）；bench 行对照 |
| BCF 头块重写路径的正确性（字典重映射） | C4-b' 独立 e2e（VCF 路径输出转 BCF 读回比对）+ noodles 差分（v13 `bcf_vs_noodles` 先例） |
| bcftools reheader/annotate 头行重排差异导致对拍假红 | C0-3 先实测差异点登记 B14 表，等价集比较（v13 B13 同款放宽） |
| `popg` 名称用户接受度 | C14-1 拍板项；`popgen` alias 兜底，README 首例用 popg |
