# 06 修复路线（Z0–Z4）

> 原则：**先修损坏与静默错误，再修数值语义，最后清理与补覆盖**。
> 本轮的修复对象是 01-review.md §2 的缺陷清单；不重开 v6 的 Y 系列结构项。
> 每一步都必须满足 §1 通用门。

---

## 1. 通用验收门（每一步都必须满足）【规】

1. `cargo test --release --workspace` 全绿；`cargo clippy --release --workspace --all-targets` 不新增警告。
2. **每个 P0/P1 缺陷先写出复现它的失败测试，再修复**（红 → 绿），测试名描述缺陷本身（见 01 §5）。
3. 过滤对拍：全量 `bcftools query` 的 GT 与 AC/AN/NS md5 与 gold 一致（沿用 v6 门的 md5）。
4. 确定性：threads ∈ {1, 2, 3, 7, 56} 输出 sha256 相同（filter 与 stats 都要测）。
5. 性能：不回退超过 10%，并在 `bench/suites/e2e/results.tsv` 追加一行。
6. **对拍**：凡涉及统计数值的改动，必须给出**参考工具的实际输出**作为证据
   （`bcftools stats` / `plink2` / vcftools / scikit-allel 之一，命令与输出写进 PR）。
7. 同一 PR 内删除被替代的旧实现，不留"以后再删"的并行代码。
8. 本目录相应条目标记为已完成并注明 commit；`docs/rounds/v7/STATUS.md` 记录状态。

---

## 2. 步骤

### Z0 立即修复（P0：损坏输出 / 静默错误 / 崩溃）

| 项 | 内容 | 门 |
|---|---|---|
| D7-1 | plan 阶段拒绝 `sketch && output_bcf`；删 `let _ = header_level` | `--mind -Ob`/`--mind --indep-pairwise -Ou` 报错；`bcftools view` 可读回（若实现编码） |
| D7-2 | `indep_requested` 时创建 stdin spill | stdin 管道 + `--indep-pairwise` == 文件路径 |
| D7-3 | `run_index` 校验 BGZF/BCF 输入 | plain/gzip 报错；bgzf 索引随机访问命中 |
| D7-4 | 矩阵只写过滤后样本、应用掩码后再写平面（或 prune 读 rule 平面） | `--keep` + `--indep-pairwise` 剪枝 ID == 子集单独跑 == plink2 |
| D7-5 | 剔除样本不写矩阵平面 / `decay_range` 跳过尾部 slice | `--keep` 的 LD all 曲线 == 子集单独跑 |
| D7-6 | `ReportRequest`/`ReportJob` 增加并应用样本与位点过滤；`--sites-file`/`--bed-file` 解析前移 | `--freq -s` == 先 filter 再 freq |
| D7-7 | 按 `vcfstats.c` 重写 IDD 双桶、`n_indels`、QUAL `.`、真实样本名、PSC SNP gate | 含 indel/缺失 QUAL 夹具 → 与 `bcftools stats -s -` 逐行一致 |

> Z0 当天完成。D7-4 若与 Y1（全局分析矩阵化）同 PR，则并入其中，不单独实现。

### Z1 数值语义（P1：统计量错误）

| 项 | 内容 | 门 |
|---|---|---|
| D7-8 | `--f4-ratio` 按 5 角色解析 b/c；补 `to_args`；对齐 README/help | `A,B;C,D;X` 成功且 == 手算 |
| D7-9 | f2/f3/f4 的 estimate 与 se/z 统一尺度；修 `request.rs:544` 注释 | `estimate/se == z`（1e-9） |
| D7-10 | `f4_ratio` 分母改 `den.num.value()` | == vcftools/allel 已知值 |
| D7-11 | 多等位位点计入 `n_present`；修 `pi.rs` rustdoc | 含多等位窗口的 variant-only π == 排除后重算 |
| D7-14 | 半分型语义显式化（`--vcf-half-call` 默认并声明）；`het`/`hardy`/`fst` 统一 | 与 plink2 对应模式逐列一致 |
| D7-17 | θ_w 口径裁决（按 oracle 对齐或改名 legacy） | 默认路径 == oracle `spectrum.rs` |
| D7-20 | 多群体 Fst 补 `n_i==0`/`an_i<2` 跳过 | r=3、`an=[1,4,4]` 位点与成对/allel 一致 |
| D7-21 | `pi.rs` 溢出用 `checked_mul`/u128 | 大 window×大样本不静默回绕 |

### Z2 健壮性与契约（P1）

| 项 | 内容 | 门 |
|---|---|---|
| D7-12 | stats 的 keep/remove 文件读取错误传播；跳 `#`；接线 `--force-samples` | 缺失文件报错 |
| D7-13 | 传 `args.skip_empty`；sink 增 `levels.window` 门控；补 Bp grid 守卫 | `--skip-empty` 生效；`--levels contig` 无窗口行 |
| D7-15 | `reed index` 拒绝 BCF+Tabix | `-tbi x.bcf` 报错 |
| D7-16 | `declared_total` 补 `text.len()<need` 守卫 | 截断 BCF 报错不 panic |
| D7-18 | `--indep-pairwise` 按位置解析 r²，允许整数 | `50 5 1` 得 r²=1 |
| D7-19 | prune 线程信号量限制 in-flight 矩阵数 | 峰值 ≤ `--memory` |
| D7-22 | 重导出公开类型；`regions` 删或接线；修正危险 `Default`；补 `FilterSpec::to_args` | `cargo doc -p reed` 自洽；往返测试 |

### Z3 清理（P2）

- D7-23（variant-only 双 key）、D7-24（IO 错误路径与 `source()`）、D7-25（同位点多记录）、
  D7-26（het 文档）、D7-27（`pi.rs:57` checked + `h_bar` 口径）、D7-28（`FunnelPartial` 死字段）、
  D7-29（删闲置/死依赖；oracle 版本并入 lockstep 或注明）。
- 门：§1；无新增警告；删除项在 `docs/README.md` 一致性守卫中同步。

### Z4 测试与文档（贯穿，与 Z0–Z3 同步）

- 按 01 §5 逐项落地回归测试。
- **扩 gold 管道**（D7-29）：`run_bcftools.sh` 增加含 indel 与缺失 QUAL 的夹具；
  `lib.sh` 读 PSC 全列（含第 9 列）；新增 LD 语义与 plink2 的对拍。
- 修正 `docs/rounds/v5/05-algorithms.md:19` 等与实际覆盖不符的措辞；把"引用参考工具实际输出"写进验收纪律。

---

## 3. 顺序与并行

```
Z0（当天）──► Z1（数值语义，单独 PR）──► Z2（健壮性，可拆多个 PR）
                                     └─► Z3（清理）
Z4 与 Z0 起同步，每个 PR 补对应回归与 gold 覆盖行。
```

- Z0 与 Z1 均触及统计数值者（D7-7 与 D7-9/D7-11）应在同一里程碑内先定语义再改代码，避免二次返工。
- D7-4/D7-1 与 Y1/Y4 同域：若 Y1/Y4 先落地，则对应条目转为**验收项**（对拍/往返）而非重实现。
- D7-14 的半分型裁决是 D7-20/D7-27 的前置：先定"半分型算什么"，再改 `het`/`hardy`/`fst`/`h_bar`。

---

## 4. 本轮未覆盖、需补测的边界（报告自认的缺口）

- **BGZF 虚拟偏移的非零路径**：本轮只核实了 `run_index` 对 plain 输入返回 0（D7-3），
  bgzf 输入下 `VoMap::voff` 的 `coff<<16|text_off` 是否正确**未测**——需构造多块文件逐记录核对虚拟偏移。
- **tabix/CSI 随机访问边界**：位置 ≥ 2^29（TBI 上限）、≥ 2^16 个 contig、CSI 大位置区间**未测**。
- 建议在 Z4 补：`reed index` + `bcftools view -r` 随机访问与 `tabix -l` 的逐区间命中测试。
