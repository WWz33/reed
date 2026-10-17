# reed doc_v14：命令面重组——seqkit 式分组 help + stats/popg 拆分 + `--hist` 直方图 + `rename` 标识符编辑（2026-10-17）

> 本目录是 `docs/rounds/v13/` 之后的**第十四轮规范**，性质是**新功能轮 +
> 命令面重组**：不改 filter/convert 的任何既有语义，只做四件事——
>
> 1. **`stats` / `popg` 拆分**：现 `reed stats` 里的 QC 报告族
>    （`--freq/--hardy/--het/--missing/--depth/--summary`）留在 `stats`，
>    窗口群体统计族（`--pi/--dxy/--fst/…/--sfs/--ld-decay`）搬到新命令
>    `popg`（alias `popgen`）。拆分判据 = **报告行的单位**（位点/样本 →
>    stats；群体×窗口/全基因组 → popg）。
> 2. **seqkit 式分组 help**：命令保持扁平（单一 `Command` 枚举），顶层
>    `--help` 按组渲染（Filtering / Statistics / Transform / File
>    operations / Miscellaneous），由 `Command` 元数据驱动的自渲染器 +
>    "每个变体恰好入一组"的守卫测试实现。
> 3. **`--hist` 直方图族**（stats 新增）：内建 metric（qual/dp/dp-site/
>    gq/idl/af）+ 通用 `INFO/TAG`、`FMT/TAG`；累加器基建
>    （`IntHistogram`/`QualDist`）已在每个 QC run 里积累，本轮是渲染层 +
>    GQ/通用标签解码。
> 4. **`rename` 子命令**（新）：样本改名、染色体改名、变异 ID 模板/映射/
>    去重、header 行编辑（`-h/-f/--add-meta/--remove-meta`）——只碰
>    header + 固定列、从不解码基因型的标识符编辑器；拼写照抄 bcftools
>    `reheader`/`annotate --rename-chrs/--set-id` 与 plink2
>    `--set-*-var-ids/--update-name/--idspace-to`。
>
> 方法（与 v13 同）：参照工具源码通读（本工作区 `bcftools/`、`plink-ng/`、
> `vcftools/`、`biosoft/seqkit`）+ reed 侧可复用面盘点（§4）+ 既有命令面
> 全量核实（`crates/reed-cli/src/cli/mod.rs` 2026-10-17 实读，v13 落地后
> 的四命令形态）。
>
> **先写 doc 不改代码**（v13 同款工作方式）；实施状态见 STATUS.md。

---

## 背景：为什么这四件事是一轮

v13 之后 reed 有四个子命令（filter/stats/index/convert），命令面开始到
"需要一个分类学"的临界点：`stats` 一个命令装着两个正交家族（QC 报告与
窗口群体统计，help 里 40+ 个 flag）；convert 证明了一个动词一个模型
（矩阵折叠）的划分是可持续的。同时对照三家参照工具（bcftools 19 个子
命令、vcftools 130+ 个 flag 单命令、plink2 378 个 flag 单命令）得出两个
结论：

1. **拆分时机是现在**：QC 报告族和窗口统计族在库里早已是两个门面
   （`run_reports` / `run_stats`），CLI 不拆只是在掩盖一条存在的分界线；
   拆完后再加 `--hist`/`--geno-counts`（stats 侧）或 PCA/f-stats 新量
   （popg 侧）才不会把一个命令膨胀成 vcftools。
2. **分类学要在命令数增长前定型**：本轮只加 `rename` 一个动词，但
   分组 help 的渲染机制、参数家族（Region/Sample/Site selection）的
   展示约定、"选择参数共享、动作参数各归"的判据一次钉死，后续 v15+
   （pca/relatedness/norm/concat/isec/query）直接长在网格里。

四件事的依赖关系：拆分（2）是前提——`--hist`（3）进拆分后的 stats 才
不会被窗口 flag 干扰；分组 help（1）在拆分定型后渲染才有意义；`rename`
（4）与前三件正交，但它的"只碰标识符不碰基因型"边界要靠分组（Transform
组）来表述。

## 文件

| 文件 | 内容 |
|---|---|
| [01-review.md](01-review.md) | 命令面现状核实（§1：2026-10-17 四命令实读）；bcftools/vcftools/plink2/seqkit 四家参照台账（§2：子命令/flag/分组机制，行号或实测）；拆分判据与搬迁表（§3：每个 flag 去向）；reed 可复用面盘点（§4）；差异分档 D14/B14/C14（§5）；守卫清单 1–14（§6）；对拍真值门设计（§7）；性能基线（§8） |
| [02-architecture.md](02-architecture.md) | 修订规范：分组 help 渲染器（§1）、stats/popg 命令面与搬迁重定向（§2）、`--hist` 请求模型与输出布局（§3）、`rename` 命令面与门规则（§4）、共享参数组判据（§5）、库 API 改名（§6）、错误与退出码（§7）、验收门（§8）、文档回改（§9） |
| [06-roadmap.md](06-roadmap.md) | 实施路线 C0–C5：真值与基线 → 拆分 + 重定向（纯搬运）→ 分组 help → `--hist` → `rename` → 文档收口；每步通用门；"本轮不做"；风险 |
| [STATUS.md](STATUS.md) | 实施状态（全部未开始）+ 待拍板 C14 + 基线表（待测） |

## 结论速览

| 维度 | 结论 |
|---|---|
| 拆分 | `reed stats <INPUT> --freq/--hardy/--het/--missing/--depth/--summary/--hist/--singletons`（QC，行=位点/样本；`--geno-counts` 留 backlog，C14-3）；`reed popg <INPUT> --pop FILE --pi/--dxy/--fst/--pbs/--sfs/--ld-decay/--f2/f3/f4/… --window/--step/--levels`（群体统计，行=群体×窗口）；QC flag 零改动，popg flag 原拼写搬迁；旧 stats 上的群体 flag 解析期报错并给出改写命令（removed-flag 重定向，exit 2） |
| 分组 help | 单一 `Command` 枚举不变；`reed --help` 按 `GROUPS: &[(&str, &[&str])]` 渲染五组（Filtering / Statistics / Transform / File operations / Miscellaneous）；help 模板驱动 + "每变体恰好一组"守卫测试 + README 命令表同源生成；**无嵌套子命令**（cli.md §0.9 兼容承诺：命令路径一个字符不变） |
| `--hist` | `--hist LIST`（逗号分隔 metric 与 `INFO/TAG`/`FMT/TAG`）：内建 `qual,dp,dp-site,gq,idl,af`；输出 `PREFIX.<m>.hist` TSV（bin_lower,bin_upper,count,percent）+ 整数 metric 附精确 mean/median/p95（`IntHistogram::quantile` 精确分位，bcftools 只有固定 bin——README 卖点）；累加器定律测试沿用 |
| `rename` | `reed rename <INPUT> [-s/-N/--replace-sample-name] [--rename-chrs FILE] [--chr-style …] [--set-id PATTERN\|--set-missing-var-ids PATTERN] [--update-name FILE] [--idspace-to SEP] [--rm-dup-ids] [-h/-f/--add-meta/--remove-meta] -o OUT`；求值序 header→染色体→ID→样本；纯标识符编辑，不解码基因型，单遍流式；`reheader` 作 visible alias |
| 共享判据 | **选择参数共享，动作参数各归**：`RegionArgs`/`SampleArgs`/`SiteFilterArgs` 三组继续被 filter/stats/popg flatten（现状）；`rename`/`convert` 不收 Site 组（编辑/转换无"选择"语义）；`stats` 的 `--pop` 变可选（分层 QC） |
| 库 API | `run_stats` → 驱动 popg，改称 `run_popgen`；`run_reports` → 驱动 stats，改称 `run_stats`；库无兼容承诺（cli.md §0.9），同轮 rename + `to_args` 往返测试改造 |
| 版本 | 0.4.0（pre-1.0 破坏性 minor）；CHANGELOG + `docs/spec/cli.md` 拆节 + README 双语重写命令面 |
| 有意分歧 | B14-1…6（见 01 §5）：`--hist` 通用标签（bcftools 无）、精确分位数（bcftools 固定 bin）、`rename` 不改 POS/REF/ALT、`--chr-style` 有限枚举、`--set-id` 多等位 ALT 拼接符、`--rm-dup-ids` 后缀规则 |
| 待拍板 | C14-2…6（见 01 §5）：裸 stats 默认行为、`--geno-counts/--singletons` 是否进首发轮、`rename` 收不收 `--order-samples`、hist 浮点 bin 边界策略、index 失效自动重建（C14-1 `popg` 已拍板） |
| 许可 | 无新参照代码移植；bcftools/plink2/seqkit 只做行为对照（GPL 工具只实现行为不复制代码——v13 同款纪律） |

## 与 v13 及之前各轮的关系

- **v13（convert 轮）**：本轮不碰 convert 的九个目标与任何语义；分组
  help 把 convert 归入 Transform 组是纯展示变化。`build_sfs_request`
  共享（D13-6）在拆分后变成 **三入口**（stats 遗留兼容入口不再有 `--sfs`
  → 实际是 convert/popg 两入口 + stats 报错引导），`translate` 层的
  共享纪律原样延续。
- **v12（INFO 兼容轮）**：`rename` 不触碰 INFO 值语义；`--remove-meta`
  只删行不删字段，INFO 字段操作仍是 filter 的 `-x`。
- **v8/v9（线程/错误路径轮）**：`--hist` 走既有 chunk 流水线，累加器
  满足 §8.4 定律即自动获得确定性；`rename` 的字节编辑无 seam。
- **cli.md §0.9 兼容承诺**：受影响的是"参数语法与输出文件格式是契约"
  ——QC flag 拼写与输出布局零变化；搬迁 flag 的拼写不变、只换所属命令，
  以 removed-flag **显式报错 + 精确改写提示**处理（`--max-missing` 家族
  的既有模式，非静默别名），pre-1.0 minor 版本承载。
- **架构 spec §10（QC 层）**：`--hist` 是 §10.2"可合并分布结构"的
  CLI 落地；`--geno-counts`（若进轮）是 §10 报告族的补全。

## 决策状态

| 决策 | 状态 |
|---|---|
| stats/popg 拆分方向（判据 = 报告行单位） | 已定（02 §2） |
| 分组 help 用扁平命令 + 自渲染器，不用嵌套子命令 | 已定（02 §1） |
| `--hist` metric 集与通用标签语法 | 已定（02 §3） |
| `rename` 命令面（四平面）与边界（不动 POS/REF/ALT/GT） | 已定（02 §4） |
| C14-1 `popg` 主名 + `popgen` visible alias | **已拍板**（2026-10-17，01 §5） |
| C14-2 裸 `reed stats`（无统计 flag）报错还是默认 `--summary` | **待拍板**（默认报错，保持显式家风） |
| C14-3 `--geno-counts`/`--singletons` 是否进首发轮 | **待拍板**（默认进 `--singletons`、`--geno-counts` 留 backlog） |
| C14-4 `rename --order-samples` 是否进首发轮 | **待拍板**（默认不进，留 merge/concat 轮） |
| C14-5 `--hist` 浮点标签 bin 边界（固定 10 bin vs 参数化） | **待拍板**（默认 bcftools 同款固定边界） |
| C14-6 `rename` 后旧索引失效：警告还是默认重建 | **待拍板**（默认警告 + 提示 `reed index`） |

## 标记约定

【规】必须；【荐】推荐；【决】需要拍板。
证据通道：**【据】** bcftools/plink-ng/vcftools 源码行号（本工作区副本）
或 reed 源码行号；**【实】** gold 对拍脚本实测（对拍门随 C0 建立）。
