# 01 复核清单：convert 内部正确性、契约空白与下沉边界

> 立项 2026-11-11。审查基于当前工作树（HEAD `e920663` 加未提交改动；`crates/reed/src/convert/*`
> 已含 v28 的 `LineRef` / `AtomicWriter` 迁移），不是干净提交。
>
> 标记：`confirmed` = 已用当前源码/独立夹具复现；`code` = 确定性代码路径可证明，
> 但未构造端到端夹具；`design` = 工程收敛建议，不是缺陷；`withdrawn` = 上轮主张经复核不成立。
>
> 审查通道：7 个只读子代理分别覆盖架构/并发、基因型正确性、HapMap/MSA、CLI 契约、
> 测试边界、性能资源、补丁热点；主审另做源码复核与最小复现。子代理结论凡与本文件不一致的，
> 以本文件的更正为准（见 §E）。

## A. 已实测确认

| ID | 状态 | 位置 | 触发与证据 | 最小处理方向 |
|---|---|---|---|---|
| V29-01 | **confirmed** | `crates/reed/src/convert/job.rs:699-701`、`:780-789` | `--nexus-binary` 用 `gt.raw_start` 作切片起点，但 `raw_start` 保存的是**样本列起点**，`raw_len` 才是 GT 子字段长度。独立夹具 `FORMAT=DP:GT`、样本 `30:0/0`、`40:1/1` 输出 `? ?`；同一记录 `--012` 为 `0`、`2`。 | 让 `format_layout` 返回 GT 子字段的真实 span（起点相对样本列），或单独保存 `raw_offset`；二进制 NEXUS 按该 span 查 `GEN_BIN`。把 `v25_02_gt_position_follows_the_format_order` 扩到 PHYLIP/012/binary NEXUS 三目标。 |
| V29-02 | **confirmed** | `crates/reed/src/convert/job.rs:586-591`、`:508-520`、`:552-573`；`crates/reed-io/src/vcf/genotype_scan.rs:193-207` | `sfs_accepts_gt` 检查**全部** `scratch.gts`，不看 keep mask；底层统计 scanner 也先解码全部样本。独立夹具保留 `s0=0/1`、排除 `s1=0/9`，`--sfs --proj 2 -s s0` 整条位点被丢，全零谱；先物理裁掉 `s1` 的期望是 `0 1 0`。 | keep mask（以及只属于 SFS 群体的样本）必须先于 GT 校验；只改 `sfs_accepts_gt` 不够，stats scanner 侧也要按参与样本解码。回归：SFS 的 `-s` 等价于预先裁剪文件，覆盖被排除样本的越界、多倍体、长前导零 GT。 |
| V29-04 | **confirmed（行为）** | `crates/reed/src/convert/hapmap.rs:369-395`、`:418-427` | 双字母格逐字符处理，遇到 `N` 直接 `continue`；只剩一个已调用等位时 `set.len()==1` 写成 `i/i`。实测 `AN → 0/0`、`NG → 1/1`、`NN → ./.`。 | 双字母分支独立保存两个 allele slot：1 缺失 1 已调用写 `i/.` 或 `./i`；两缺失写 `./.`。D29-02 已拍板 A（保留半缺失；bcftools 读入 half-call 原样保留）。补 `AA/AG/AN/NA/NG/GN/NN` 表驱动，禁止静默纯合化。 |
| V29-05 | **confirmed（行为）** | `crates/reed/src/convert/job.rs:1780-1810`；`crates/reed-io/src/chunk/plain.rs:95-107` | 只有 header、无 record 的 VCF 不会产生 chunk，`merged.sfs` 恒为 `None`，而 SFS 写出要求 `merged.sfs` 存在；实测 `--sfs --proj 4` exit 0 且**不生成任何 SFS 文件**。有 record 但全被跳过时（v25 测试）会写零谱。 | 请求存在即初始化零 spectrum/零 report，与 v27-05“请求决定输出集合”一致；D29-05 已拍板 A。header-only 的 plain/BGZF/BCF 与“全记录被过滤”统一断言。 |

## B. 代码路径确认（拍板后修）

| ID | 状态 | 位置 | 触发与证据 | 最小处理方向 |
|---|---|---|---|---|
| V29-03 | **code / CLI 场景可确认** | `crates/reed/src/convert/writers.rs:224`、`:209`、`:217`、`:315-316`；`crates/reed-cli/src/translate/convert.rs:182-190`；`crates/reed/src/convert/job.rs:1450-1485`；`crates/reed/src/convert/hapmap.rs:487-608` | CLI 只用 `refuse_output_eq_input(input, args.out)` 比较 `-o` **原文**与输入。输入 `foo.012` + `--012 -o foo` 时最终产物恰为 `foo.012`；MSA 重复 `--region` 或 flank 收敛后写同一路径；直接调库时 `ConvertTarget::Vcf{output: Some(input)}` 无任何别名检查。结果是静默覆盖输入/参考/同一 run 的其他产物。 | 在任何 writer 打开前计算完整输出集合（含 `.log`、派生 sidecar、每个 MSA region），检查重复路径与 input/reference/population 文件别名；CLI 与库入口共用实现。D29-04 已拍板 A（冲突 exit 2）；保持逐文件原子，不做多文件事务。 |
| V29-06 | **code** | `crates/reed/src/stats/run.rs:381-399`；`crates/reed/src/convert/job.rs:1810` | `StatsPartial::merge_sfs` 只合并 spectrum，不合并 `n_records`；convert 最后用 `partial.records()` 覆盖计数。干净单 chunk 实测 `records: 1`；多 chunk 时只保留第一个产生 SFS partial 的 chunk。上轮“恒为 0”表述已更正。 | 在 `merge_sfs` 累加 `n_records`，并明确它统计“进入统计扫描的记录数”还是“全部输入记录数”；补多 chunk 回归。 |
| V29-07 | **code** | `crates/reed/src/convert/writers.rs:315-316`；`crates/reed/src/convert/job.rs:1384-1430` | MSA 输出路径为 `PREFIX_<contig>[_start-end].fasta`，contig 只做 `from_utf8_lossy`，未限制 `/`、未查重。重复 region 会覆盖同一文件；含 `/` 的 contig 会落到带前缀的子目录（`out_../victim.fasta` 这类例子本身不会直接逃到父目录，不能把风险夸大）。 | contig 作为文件名片段做严格校验或可逆编码；`build_msa` 规划阶段就检测最终路径撞名并 exit 2。 |
| V29-08 | **code** | `crates/reed/src/convert/job.rs:1338-1370`；`crates/reed/src/convert/writers.rs:283-345` | `read_reference` 同名 FASTA record 直接 `map.insert` 覆盖；MSA 渲染从不比较 VCF REF 与参考切片，也不检查 `POS + REF.len()` 的跨度。越界 POS 会在 `pos0 >= end0` 处跳过，**未发现上轮所称的切片 panic**。 | 重复 contig 名拒绝；REF 与 FASTA 大小写不敏感比较，不一致时报错；D29-01 B 已定义重叠规则：只有相同 POS 且相同 REF span 允许 last-wins，其余重叠直接拒绝。 |
| V29-09 | **code** | `crates/reed/src/convert/job.rs:744-757`；调用点 `:619`、`:682`、`:986`、`:1031`、`:1127` | 每个 target 各自调用 `bump_malformed`，同一条 `0/9` 在 `--phylip --012 --hapmap` 下计入 3 次，可出现 `malformed > total_records`，与 `ConvertReport` 文档“skipped malformed records”冲突。 | 拍板：改成记录级计数一次 + 可选 per-target 明细，或把字段改名并明确目标是“拒绝次数”。 |
| V29-10 | **code** | `crates/reed/src/convert/mod.rs:64-66`、`:163-282`；`crates/reed/src/convert/hapmap.rs:216-241`；`crates/reed/src/convert/job.rs:1450-1495` | 库入口未统一校验：`level` 注释 1–9 但未验证，`BlockWriter::new(...).expect("CLI clamp")` 把非法库请求提升为 panic 风险；重复 SFS/MSA target 静默覆盖配置；`to_args()` 不渲染 populations、keep mask、seed 等，不能作为完整重放表示。 | 引入 `ConvertRequest::validate` 与 `ConvertPlan::compile`；压缩级别在驱动入口返回 Usage；明确重复 target 规则；`to_args` 更名/文档说明为 target-only，或补全字段。 |

## C. 上轮主张的撤回或降级

| ID | 上轮主张 | 裁决 | 当前证据 |
|---|---|---|---|
| W29-01 | MSA 长 REF 游标只前进 1，属于 P0；应改成按 `ref_len` 前进 | **withdrawn（作为 Bug；已转为契约升级 D29-01 B）** | `e2e_convert.rs:1225-1227` 明确断言 `300 + 999 = 1299`；`docs/rounds/v15/01-review.md` D15-26 把 1299 写成验收标准，`02-architecture.md` 记录实现取“REF block + 每剩余位点一列”。上轮“本应在第 6 位的变异落到第 7 位”也不成立：实测 10 bp 例子输出 `AAAAAGAAAA`，G 在输出第 6 位；长 REF 的真实现象是下游多出 `REF_len − 1` 列。**2026-11-11 拍板 D29-01 B：默认改为参考坐标对齐（示例 1000），D15-26 被显式取代；这不是按 Bug 静默修，而是产品契约升级。** |
| W29-02 | folded 一维 SFS 把 `k=0/n` 合并到 masked bin 0 是错误 | **withdrawn** | `docs/rounds/v5/05-algorithms.md` §3.1 要求非分离质量保留在数组里；`e2e_sfs.rs:136-143` 明确“k = 0/n pair 是单次贡献”。`sfs.rs:269-284` 的实现与文档一致，不是 Bug。 |
| W29-03 | 多倍体导致 SFS 静默空谱、012 写 `-1`，是 P1 数据错误 | **降级为策略空白** | `e2e_convert.rs:493-560` 的 v25 用例把 polyploid 归入“合法空谱”；`--012` 的 `n != 2 → -1` 是文档登记的 haploid/half-call 兼容选择，多倍体未单独定义。需要的是显式诊断/文档，不是直接改算法。 |
| — | MSA REF/FASTA 越界会 panic | **withdrawn** | `writers.rs` 的 `pos0 >= end0` 提前 `break`；`region.end1` 也会 `min(len)`。越界表现为跳过或输出与参考不一致，不是 panic。 |
| — | SFS `n_records` 恒为 0 | **更正为 V29-06** | 单 chunk 正常输入实测 `records: 1`；缺陷只在多 chunk 合并时丢后续计数。 |
| — | HapMap `alleles=NA` 仍写 `0/0` | **不成立** | `hapmap.rs:295-317`：`alleles.is_empty()` 时写 `REF=.`、全部 `./.`；v13 STATUS C4-b 同此。真正遗留的是 `REF=.` 的 VCF 合法性问题。 |
| — | MSA contig `../victim` 可直接覆盖任意文件 | **降级** | 输出名带 `out_` 前缀，`out_../victim.fasta` 不会跳到父目录；只有能构造出 `out_x/../../victim.fasta` 这类路径时才存在越界写入条件。核心问题仍是未限制 `/` 与未查重。 |

## D. 下沉（design，不是缺陷）

| ID | 能力 | 当前重复 | 边界 |
|---|---|---|---|
| V29-11 | 无损 GT 词法与真实子字段 span | convert `job.rs:762-810`、PCA `genotype.rs`、严格 scanner `genotype_scan.rs:311-484` 三套 GT 文本解析；Git 历史 `61273ab`、`e182ef6`、`10f7242` 曾同步修改 | 只拥有 `/`、`\|`、缺失、整数/ALT 域、真实 span、原始 phasing 前缀；不决定“跳过/NaN/fatal”，不替代各消费者的倍性、多等位、样本过滤策略。先迁移 convert，再 PCA，最后 scanner。V29-01 应直接用这个 span 修。 |
| V29-12 | 请求校验与输出计划 | CLI `translate/convert.rs:178-353` 做了完整 gate，库入口 `job.rs:1450-1495` 部分缺失；driver 把 targets 展开成多组布尔/Option，重复 target 行为不一 | `ConvertRequest::validate()`（不依赖 header）+ 内部 `ConvertPlan::compile(facts)` + `planned_outputs()`；输出命名规则留在 convert，安全的“声明/别名检查”原语可放 `reed-io`。不是 `ConversionBackend` trait，不重设计公共 builder。 |
| V29-13 | 已筛选记录流与空 SFS 报告 | convert `job.rs:494-545` 先折叠、再复制 sanitized body 给统计；空谱在 drain 里补丁式构造 | `StatsJob::run_lines(lines)` + `sfs_report_or_empty(plan)`；convert 保留 malformed 兼容策略与 `DrawKey` 行序，stats 保留严格 fatal 策略。 |
| V29-14 | HapMap 文本 chunk / gzip sniff / spill 注册 | `hapmap.rs:128-190` 私有分块器；`hapmap.rs:490-540` 两路 gzip 探测；`job.rs:424-471` 两份 blob 清单 | 复用 `PlainChunkSource::from_reader` + 统一 gzip sniff；spill 的“内存转临时文件”收敛为一个函数、注册收敛为一份清单。不要发明通用 `TextFormat` trait。 |

## E. 不采纳的子代理结论

| 结论 | 不采纳理由 |
|---|---|
| MSA 长 REF 是 P0，按 `ref_len` 前进，输出应为 1000 | **作为“Bug”不采纳**：与 D15-26 的 1299 验收冲突；但 **D29-01 B 已将其作为显式契约升级采纳**，默认改为参考坐标对齐并重开 D15-26。 |
| folded SFS 端点重复计数是 P1 | 与 v5/05 §3.1 及现有测试锚点冲突；见 W29-02。 |
| 多倍体 SFS/012 是已证实数据错误 | 属策略空白，现有测试把 polyploid 归入合法空谱；见 W29-03。 |
| `POS=0` 按普通越界报错 | VCF 4.3 对端粒允许 0 或 N+1；本轮不引入该规则。 |
| MSA REF 越界会 panic | 当前控制流会跳过，见 §C。 |
| SFS `n_records` 恒为 0 | 单 chunk 实测非 0；真缺陷是多 chunk 合并，见 V29-06。 |
| HapMap contig `../` 可直接逃逸 | 前缀使该例子不成立；保留“未限制 `/` + 未查重”的问题。 |

## F. 复现与验证环境

- 独立复现目录：`/tmp/reed-convert-review-YmOWiybHCj`（V29-01/02/04/05）。
- 构建目录：`/tmp/reed-review5-target`（避开共享 target 的 `proc-macro2` probe `Device or resource busy`）。
- 现有测试（当前工作树最终一次运行，独立 target）：`e2e_convert` 25 通过、convert lib 7 通过、
  `cli_contract` convert 4 通过、`--bin reed convert` 4 通过，另有 7 个选定契约测试通过。
- `gold/popgen/run_convert_parity.sh` 通过，但它调用 `/tmp/reed-target/release/reed`（预构建），
  **不是当前源码的等价构建验证**。
- 未运行完整 workspace、8 GB 全量、真实磁盘写满/断流、多 chunk SFS 记录数夹具；这些在
  STATUS 中列为未验证。
- 工作区在本轮审查期间有其他并行修改；V29 行号对应 2026-11-11 的当前工作树，
  提交前需按提交时点复核。
- D29-02 的 bcftools 依据（2026-11-11）：仓库 bcftools 的 `view -Ov`（header 标记
  `bcftools_viewVersion=1.24-30-g347c2732+htslib-1.25-7-g02469fbe`）读入
  `GT=0/.`、`./1` 时原样输出，不纯合化、不改写 `./.`；`query -f '[%GT]\n'`
  同样逐字返回。故 HapMap→VCF 写 `0/.` / `./1` 与 bcftools 的读入行为一致。
  该二进制不支持 `--vcf-half-call`（`view` 报 unrecognized option），因此不能拿它
  证明“改成 missing”的默认策略。
