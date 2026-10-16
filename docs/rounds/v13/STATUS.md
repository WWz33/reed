# doc_v13 实施状态

> 本页记录 `06-roadmap.md` C0–C4 的实施状态与证据。本轮是**新功能轮**
> （`reed convert` 子命令），无缺陷修复清单；复审对象是"实现 vs 参照工具
> 行为"的对拍，真值门 = vcf2phylip 2.9 / easySFS / vcftools / TASSEL /
> vcf2msa。每步完成即更新。

## 实施状态（C0–C4）

| 项 | 内容 | 状态 | 证据 |
|---|---|---|---|
| C0-1 | phylo.vcf fixture + 六跑法期望输出 | **完成** | `gold/fixtures/phylo.vcf`（make_fixtures.py `phylo()`）；期望输出固化为 `crates/reed/tests/e2e_convert.rs` 的 `REF_*` 常量（vcf2phylip 2.9 现跑） |
| C0-1b | convert012/convert_msa fixtures + 读回脚本 | **完成** | `gold/fixtures/convert012.vcf`、`convert_msa.{vcf,fa}`；hapmap 真值用 TASSEL 读回（脚本门 3） |
| C0-2b | 全量性能基线 | **完成**（vcf2phylip 全量档；012/SFS 全量档低优先未测） | `clean.vcf.gz`：27 m 27 s vs 10.4 s（§基线） |
| C0-3b | run_convert_parity.sh（六门） | **完成** | `gold/popgen/run_convert_parity.sh`：`convert parity: PASS` |
| C0-2 | vcf2phylip / stats --sfs 基线墙钟 | 未开始 | — |
| C0-3 | run_vcf2phylip.sh 门脚本 | 未开始 | — |
| C1-a | reed-core iupac 模块（AMBIG/GEN_BIN） | **完成** | `crates/reed-core/src/iupac.rs`（62 + 10 项，排序与完整性单测） |
| C1-b | ConvertCli + Command::Convert 分发 | **完成** | `crates/reed-cli/src/cli/mod.rs`（ConvertCli）、`translate/convert.rs`、`main.rs` 分发 |
| C1-c | convert 驱动 + 门面 | **完成** | `crates/reed/src/convert/{mod,job,writers}.rs`；`run_convert`/`run_convert_opened` 根导出 |
| C1-d | phylip 对拍全等 e2e | **完成** | `e2e_convert::phylo_outputs_match_vcf2phylip_reference`（五文件逐字节） |
| C2-a | fasta/nexus 渲染 + 三格式同模型 | **完成** | `writers::write_phylo`（一份列模型三渲染） |
| C2-b | nexus-binary 独立模型 | **完成** | `GEN_BIN` 查表 + 双等位门；`REF_BIN_NEXUS` 逐字节 |
| C2-c | resolve-iupac / used_sites / outgroup | **完成** | splitmix64 派生（chunk,line,sample,seed）；`min_samples_locus_gates_and_clamps`、`resolve_iupac_is_deterministic`（threads 1==4） |
| C2-d | 六跑法对拍转绿 + threads {1,4} | **完成** | `run_convert_parity.sh` 全绿；混合全目标（sfs+5 格式+012+hapmap+resolve）threads 1==4 |
| C2b-a | --012 对拍 vcftools | **完成** | 干净夹具逐字节（三文件）；B13-8 行断言 `-1`（`count_matrix_matches_vcftools_and_flags_half_calls`） |
| C2b-b | --hapmap + TASSEL 读回 | **完成** | TASSEL 5 导入 reed 输出再导出成功、17 站点保持（脚本门 3）；字母表收敛 B13-10（`hapmap_layout_and_alphabet`） |
| C2b-c | --msa 对拍 vcf2msa | **完成** | 三样本序列逐字符相等（`msa_matches_vcf2msa_on_its_domain` + 脚本门 4，区域映射 `B-(E+1)`） |
| C2b-d | 三格式确定性 + 混合跑法 | **完成** | 混合全目标 threads 1==4（e2e + 手工验证） |
| C3-a | build_sfs_request 共享 + ConvertTarget::Sfs | **完成** | `translate::build_sfs_request`（stats/convert 共用）；`StatsJob::for_sfs` 组合（方案 A，GridMode::None 不造窗口） |
| C3-b | convert ≡ stats 逐字节 + 往返测试 | **完成** | `convert_sfs_equals_stats_sfs`（1D）+ `diff -r` 同前缀树（1D/2D、folded/unfolded 全等）；`convert_and_stats_sfs_spellings_agree` |
| C3-c | 混合确定性 + 性能门 | **完成** | 混合全目标 threads 1==4；全量 158×（门 ≥8×）、五文件逐字节（见 §基线） |
| C3-d | 文档收口 | **完成** | cli.md §8b + README 双语 + CHANGELOG + 守卫 1–15（`check-doc-drift.sh` 全绿） |
| C4-a | hapmap 输入探测与装载 | **完成** | `convert/hapmap.rs`：`read_header`（`##` 跳过 + `rs#` 列行 + 首行宽度嗅探，嗅探行回放给 source）；`HapmapSource`（1 MiB plain chunk） |
| C4-b | 逐行翻译（IUPAC → GT） | **完成** | `HapmapVcfJob::translate` + `cell_gt`（`iupac::ambig_expand` 逆表）；`alleles=NA` → 全 `./.`；malformed 跳过+计数；one/two-letter 两态 |
| C4-c | `--vcf` 目标 + `-O` + CLI 门规则 | **完成** | `ConvertTarget::Vcf` + `run_hapmap_to_vcf`（`OrderedOutput` 流式）；`-O v\|z[level]`（`b`/`u` usage 错）；`-o` 扩展名推断 BGZF；规则 6 互斥；无 `##contig` |
| C4-d | 往返保真 + TASSEL 可读性 + 守卫 16–18 | **完成** | `hapmap_vcf_round_trip_is_faithful_and_deterministic`（记录逐行还原 + threads 1==8 + `-Oz` 帧）；脚本门 6（含 TASSEL `-importGuess` 读 reed 的 VCF）；守卫 16–18 全绿 |

## 落地文件清单

| 文件 | 内容 |
|---|---|
| `crates/reed-core/src/iupac.rs`（新，235 行） | AMBIG（62）/GEN_BIN（10）常量 + 查表 + `ambig_expand` 逆表 + 单测 |
| `crates/reed/src/convert/{mod,job,writers,hapmap}.rs`（新，2795 行） | 请求类型、chunk job/partial/sink、六个 VCF→X 目标、`hapmap.rs` 的 X→VCF 反向目标 |
| `crates/reed-cli/src/translate/convert.rs`（新，479 行） | 门规则（含 `--vcf` 门）、请求构建、`.log`、timing |
| `crates/reed/tests/e2e_convert.rs`（新，812 行） | 13 个对拍/契约测试（真值 = 参照工具实测输出） |
| `gold/popgen/run_convert_parity.sh`（新） | 六门对拍脚本（含 hapmap↔VCF 往返 + TASSEL 读 reed 的 VCF） |
| `gold/popgen/make_fixtures.py`（改） | `phylo()`/`convert012()`/`convert_msa()` 夹具生成 |
| `crates/reed-cli/tests/cli_contract.rs`（改） | `convert --vcf` 三条契约测试（记录/头/门规则/stdout 流式/非 hapmap 输入） |
| `crates/reed/src/stats/{driver,request,run,sink}.rs`（改） | `resolve_groups`/`build_sfs_plan`/`sfs_report` 提取共享；`GridMode::None` + `StatsJob::for_sfs`（组合复用，stats 行为零变化） |

## 差异登记落地状态（B13 → cli.md 差异表）

| 项 | 状态 |
|---|---|
| B13-1 `--phylip` opt-in | 未落地（随 C1-b 进 help，随 C3-d 进 cli.md） |
| B13-2 畸形行跳过+计数 | 未落地 |
| B13-3 resolve-iupac 固定 seed | 未落地 |
| B13-4 outgroup 报错 | 未落地 |
| B13-5 hapmap alleles 声明集 | **已落地**（C2b-b；`hapmap_layout_and_alphabet`） |
| B13-6 msa indel 确定性填充 | **已落地**（C2b-c；无外部对齐器） |
| B13-7 msa 样本名不折叠 | **已落地**（`--collapse-suffix` 显式开启；重名报 usage 错） |
| B13-8 012 半调用/haploid = -1 | **已落地**（C2b-a；脚本与 e2e 双断言） |
| B13-9 AMBIG 未命中跳过（参照崩溃） | **已落地**（C1-c；与索引越界同口径；`ambig_miss_and_index_overflow_are_skipped_and_counted`） |
| B13-10 hapmap 格收敛 TASSEL 字母表 | **已落地**（C2b-b） |
| B13-11 hmp→VCF 直抄 `chrom`/`pos`/`rs#`/`alleles` | **已落地**（C4-b；`ID` = `rs#` 原样，REF/ALT 按声明顺序） |
| B13-12 hmp→VCF 无 `##contig`、`##reed_*` 可关 | **已落地**（C4-c；`--no-version` 去 `##reed_*`，`--no-command` 只去 argv 行） |

## 待拍板（C13，默认已定可推翻）

| 项 | 默认 |
|---|---|
| C13-1 GT 掩码族不进 convert | 已定（02 §2.5） |
| C13-2 `--min-samples-locus` 独立参数 | 已定（02 §3.7） |
| C13-3 SFS 目标 `--pop` 可缺省 | 已定（02 §1.2） |
| C13-4 msa 样本名不折叠 | 已定并落地（默认不折叠；折叠重名报错） |
| C13-5 plink 不进首发轮 | 已定（02 §3.11） |
| C13-6 相位/似然格式本轮不做 | 已定（06 §4） |

## 测试与守卫（落地后）

| 命令 | 结果 |
|---|---|
| `cargo test --workspace`（debug / release） | 全绿（497 + 新增 convert/CLI 测试；两 profile 均通过） |
| `cargo clippy --workspace --all-targets -- -D warnings` | 无告警 |
| `bash ci/check-doc-drift.sh` | 全绿（87 项，含 v13 守卫 1–18） |
| `bash gold/popgen/run_convert_parity.sh`（release 二进制） | `convert parity: PASS`（六门全绿） |
| `bash gold/popgen/run_easysfs.sh` | `SFS parity: PASS (rtol 1e-6)` |

## 基线（2026-10-16 实测，`clean.vcf.gz` = 8 GB / 348 样本 / 3 315 604 记录，28 核）

| 负载 | 参照（单线程） | reed convert | 实测 |
|---|---|---|---|
| phylip+fasta+nexus+bin+used_sites 一次跑 | vcf2phylip 2.9 **27 m 27 s** | **10.4 s**（user 3 m 33 s，峰值 RSS 8.2 GiB） | **≈158×**；五文件 `cmp` 逐字节相同 |
| 线程确定性（全量） | — | threads 28 vs 4 | `.phy`/`.bin.nexus` `cmp` 相同 |
| `--012` 全量 / `--sfs --proj` 全量 | 未测（低优先） | — | fixture 档门已绿 |
| 峰值内存拆解 | — | — | mmap 输入页 ~8 GB（可回收）+ 列矩阵 ~2.3 GB + worker 缓冲 ~0.2 GB；`--memory` 只管拷贝式源 |

`gold/popgen/run_convert_parity.sh DATA.vcf.gz` 是这条全量门的脚本形态
（跑参照 + reed 各一遍并 diff 五文件；本表数字即该负载的实跑）。

## 本轮未做（有意）

见 `06-roadmap.md` §4：stats --sfs 保留；GT 掩码族；`--geno` 换算；
`--seed`；>2D 联合 SFS；sites-only VCF；interleaved PHYLIP；convert 的
JSON report；convert 区域参数；`--vcf` 的样本子集（hapmap 列头即样本集）、
`--vcf` 的 BCF 输出（`-Ob/-Ou` usage 错）、`--vcf` 的 `##contig` 补全
（TASSEL 同款不写；下游要长度走 `--contig-lengths`）。
