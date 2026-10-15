# doc_v12 实施状态

> 本页记录 `01-review.md` §1 缺陷清单与 `06-roadmap.md` R0–R3 的实施状态。
> 每条附证据（测试名 / 守卫 / 变更点）。复审基线（修复前的实测证据摘录）
> 见文末。

## 实施状态（R0–R3）

| 项 | 严重度 | 状态 | 证据 |
|---|---|---|---|
| D12-1 sketch 下尾部表达式丢行字段 | **P0** | 未修 | — |
| D12-2 无 GT 记录 AC/AN 被子集更新清零 | **P0** | 未修 | — |
| D12-3 Integer INFO 超 int32 回绕 | P1 | 未修 | — |
| D12-4 Flag 带值丢值、`DB=0` 反转 | P1 | 未修 | — |
| D12-5 String 值 `.` 变空串 | P1 | 未修 | — |
| D12-6 `-q/-Q` cohort 变化后仍用陈旧 INFO | P1 | 未修 | — |
| D12-7 表达式无类型（Flag/向量/裸标签） | P1 | 未修 | — |
| D12-8 `quantifier = "any"` 恒真 | P1 | 未修 | — |
| D12-9 打包计数溢出无检查 | P2 | 未修 | — |
| D12-10 `ALT='.'` 计 1 个 ALT | P2 | 未修 | — |
| D12-11 `--fill-tags` 不校验 | P2 | 未修 | — |
| D12-12 默认 `drop_stale` 不可见 | P2 | 未修 | — |
| D12-13 AC 长度不校验、三份实现 | P2 | 未修 | — |
| D12-14 数值解析三份实现 | P2 | 未修 | — |
| D12-15 `edit_vcf_line` 死代码 | P3（lean） | 未修 | — |

## R0–R3 进度

| 步骤 | 状态 | commit |
|---|---|---|
| R0 红测试（a–j 十组夹具） | 未开始 | — |
| R1 P0+P1 修复（D12-1…D12-8 + 守卫 1–7） | 未开始 | — |
| R2 P2 收口（D12-9…D12-14 + 守卫 8–13） | 未开始 | — |
| R3 P3 + 文档收口（D12-15 + 守卫 14 + STATUS） | 未开始 | — |

## 本轮未做（有意）

见 `06-roadmap.md` §4：sites-only VCF 支持（独立立项候选）；表达式语言扩容；
`--fill-tags AF/F_MISSING` 实现（【决】）；默认 `drop_stale` 改 opt-in（【决】）；
htslib 前缀解析全对齐（【决】）；`bcf_record_check` 全量移植；读侧 PASS 注入；
stats 的 `info.<TAG>` 直方图。

## 复审基线（2026-10-14/15，**修复前**）

复审 commit：`065175c`（doc_v10 全部落地；doc_v11 门面轮文档已暂存、与本轮
正交），代码工作区干净。reed 二进制：`cargo build --bin reed`（debug，
rustc 1.99.0）；真值：`bcftools 1.24-2-g270751fc`（htslib `1.24-1-gb499e4c7`，
`/sri/home/sri2025201067/biosoft/bcftools`）；htslib/bcftools 源码：工作区
`../../htslib`、`../../bcftools`。夹具目录 `/tmp/inforev/`（D12-1…13）与
`/tmp/parsecmp/`（D12-14 的 31 形态差分）——均临时，正文内联可重建。
codegraph 索引已同步（53 文件、1 485 节点）后用于调用面定位。

| 缺陷 | 基线证据（实测输出摘录） | 通道 |
|---|---|---|
| D12-1 | `-i 'F_MISSING<=0.5 && QD>10'`：无 `--mind` → 2 条；`--mind 0.9`（`--report` 证实 `sketch_then_emit`）→ **0 条**；`--mind 0.9 --strategy single_pass/two_pass_reparse` → 2 条。`QUAL>10` 同型。`--mind 0.9 --indep-pairwise 50 5 0.2` + 表达式 → 0 条（无表达式对照 → 1 条）。`-q 0.1 --mind 0.9`（stale 夹具）→ 0 条（无 `--mind` → 2 条） | 【实】 |
| D12-2 | nogt 夹具（`AC=4;AN=6`，FORMAT=DP 无 GT）：`bcftools view -s s0,s1` → `AC=4;AN=6`；`reed filter -s s0,s1` → **`AC=0;AN=0`**；`reed … --fill-tags AC,AN` → `AC=0;AN=0`；`reed … -I -x ''` → `AC=4;AN=6`（缺口仅在自动更新） | 【实】 |
| D12-3 | `XD=3000000000`：bcftools `-Ob` 读回 `XD=.`（+警告 "Extreme INFO/XD value … set to missing"）；reed `-Ob` 读回 **`XD=-1294967296`**（`-3000000000` → `1294967296`）。`XD=-2147483641`（htslib 保留带内）reed 原样写出 | 【实】 |
| D12-4 | `DB=1` → bcftools `-Ob` 读回 `DB=1`；reed `-Ob` 读回 **`DB`**。`DB=0` → bcftools `DB=0`；reed **`DB`**（"清除"变"设置"）。reed `-Ov` 两者均逐字节保留 | 【实】 |
| D12-5 | `SOMESTR=.` → bcftools `-Ob` 读回 `SOMESTR=.`；reed `-Ob` 读回 **`SOMESTR`**（bcftools 与 reed 读 reed 的文件结果相同——损坏在写侧）。reed `-Ov` 保留 `SOMESTR=.` | 【实】 |
| D12-6 | sub 夹具（`AC=4;AN=8`，GT s0/s1=1/1、s2/s3=0/0）`-s s2,s3 -q 0.1`：bcftools → **0 条**；reed → **2 条且输出记录写 `AC=0;AN=4`**。`-q` 全程无"陈旧 INFO"警告（`-i 'AC>1'` 有，`plan/mod.rs:422` 路径） | 【实】 |
| D12-7 | flag3 夹具：`-i 'DB=1'` bcftools 留 pos10 / reed 0 条；`DB=0` bcftools 留 pos20 / reed 0 条；`DB!=0` 同 `DB=1`。vec 夹具：`-i 'DP4>10'` bcftools 留 pos10（`DP4=0,0,30,40`）/ reed 0 条。zero 夹具：`-i 'DP'` bcftools 0 条 / reed 留 DP=5、DP=-3 两条；`-e 'DP'` bcftools 3 条 / reed 1 条。`-i 'DP>0'` 双方一致 | 【实】 |
| D12-8 | 探针 crate（`reed-filter` 库路径，`AC=1,2`）：`keep_if ">= 10"` × `any` → **true**；`"<= 0"` × `any` → **true**；`all` 两例正确。`head.rs:216` `any = true` 在可解析分支无条件置位 | 【实】【据】 |
| D12-9 | 65 536 个 INFO 片段：reed `-Ob` 写出成功；bcftools 读 **"Bad BCF record … Invalid REF/ALT type 1"**；reed 自读 **"BCF string field is not CHAR-typed"**（65 535 正常）。300 FORMAT 键：写出 `n_fmt=44`（300 & 0xFF）静默截断；bcftools 对同输入解析期报错。66 001 等位（`--max-alleles 70000`）：写出后 bcftools "Invalid FILTER type"、reed 自读 "non-integer type" | 【实】 |
| D12-10 | inv 夹具（pos10 `ALT=.`、pos20 `ALT=G`）：`-m 2` bcftools {20} / reed **{10,20}**；`-M 1` bcftools {10} / reed **{10,20}**（`max_alt_alleles=0` 被当"关闭"，规则未推入，`plan/mod.rs:352`）；`-i 'N_ALT=0'` bcftools {10} / reed **{}**；`-i 'N_ALT=1'` bcftools {20} / reed **{10,20}** | 【实】 |
| D12-11 | `--fill-tags AF` → 退出码 0、无警告、输出无 AF（输入 `AF=0.5` 被默认 `-x` 删除）；`--fill-tags BOGUS` → 退出码 0 静默；`--fill-tags ac,an`（小写）→ 静默无操作。bcftools `+fill-tags -t BOGUS` → "the tag … is not supported"，exit 1（`fill-tags.c:530-535`） | 【实】 |
| D12-12 | keep 夹具裸跑 `reed filter`：`DP=30;AF=0.5;MLEAC=1;ExcessHet=0.1;InbreedingCoeff=-0.2;FS=3` → **只剩 `FS=3.0`**（header 定义保留）；`reed filter -I` → `INFO=.`；bcftools `view`/`view -I` 均原样保留；`-x ''` → 原样保留（转义口未文档化）。`--help`/README/cli.md 无默认清单字样 | 【实】【据】 |
| D12-13 | aclen3 夹具（`AC=3` 单值 + 2 ALT + `AN=4`，GT 全 0/0）`-q 0.1`（`-M 3`）：bcftools → 警告 "Incorrect number of AC fields" + **0 条**（回退 GT）；reed → **1 条**（用 3/4=0.75）。`summary.rs:302` 的第三份实现有 `vals.len() == n_alt` 校验，`head.rs:128` 没有 | 【实】【据】 |
| D12-14 | 31 形态 × 3 路径全量差分（夹具 `/tmp/parsecmp/`，`XI` Integer / `XF` Float，真值 bcftools 1.24）。**路径 A（Integer `-Ob`）10 处不一致**（其中 `+5` 一项经 VCFv4.5 对照改判：Float 正则允许 `[-+]?`（tex:97），属规范内输入，N13 顺手对齐，不计入"不改"分歧；其余 9 项见 §2.5-B）：`+5`→reed `.`(bct `5`)、`5.0`→`.`(bct `5`)、`1e3`→`.`(bct `1`)、`0x10`→`.`(bct `0`)、`5abc`→`.`(bct `5`)、`5 `(尾空格)→`.`(bct `5`)、`3000000000`→`-1294967296`(bct `.`)、`9223372036854775807`→`-1`(bct `.`)、`-2147483649`→`2147483647`(bct `.`)、`-2147483641`→原样(bct `.`)。**路径 B（Float `-Ob`）4 处**：`0x10`→`.`(bct `16`)、` 5`/`5 `/`5abc`→`.`(bct `5`)；其余 27 形态全同（含 inf/nan/1e400）。**路径 C（`-i 'XI>4'`）reed 多保留 9 条**：`abc`(字节序 `"abc">"4"`)、`inf`(Rust parse 接受)、`1e3`(按 1000)、`3000000000`/`2147483648`/`9223372036854775807`/`99999999999999999999`(f64 解析成功)、` 5`(表达式 trim)、`5,`(Number=1 却多值：bcftools 标量读作 **0**，htslib `len != 1` 不填 `v1`，`vcf.c:4286`；reed 取首元素 5)。**内部自相矛盾**：`XI= 5` 表达式判 5(留)但 `-Ob` 写 `.`；`XI=1e3` 表达式判 1000(留)、`-Ob` 写 `.`、`-Ov` 留 `1e3` | 【实】【据】 |
| D12-15 | `line_edit.rs:86-94`：`LineRef::parse(line.line)` 重解析；`LineRef` 字段私有（`line.rs:13`）且唯一构造入口 `parse`（`line.rs:20-32`）→ 入参恒 ≥9 tab，`MalformedRecord` 分支不可达 | 【据】 |
| 读侧一致（基线） | 差分夹具 5 记录（int8/16/32 边界 `-120/-32760/32767/2147483647`、float `-0`/`1e-05`、String 含空格、Flag、`.` 缺失、`XM=1,.` 向量截断）：reed 渲染 bcftools 产 BCF 与 bcftools 自渲染 `diff` 为空 | 【实】 |
| `-q :type` 一致（基线） | aclen4 夹具五种 `:type` × 阈值 0.3：reed 与 bcftools 逐一相同 | 【实】 |
| 子集更新一致（基线，有 GT） | sub 夹具 `-s s2,s3`：双方 `AC=0;AN=4`；`-I`：双方保留 `AC=4;AN=8` | 【实】 |

## 测试与守卫基线（复审期）

| 命令 | 结果 |
|---|---|
| `cargo build --bin reed`（debug） | 成功（target-dir `/tmp/reed-target`，`.cargo/config.toml`） |
| `codegraph sync` | 53 文件、1 485 节点，索引最新 |
| 本轮全部结论 | 来自针对性 CLI 夹具 / 探针 crate 实跑 + htslib/bcftools 源码对照；未重跑全量 workspace 测试与 gold 88 MB 对拍（v10 STATUS 记录其为绿） |
