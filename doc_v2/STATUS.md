# 实施状态（2026-10-07）

> 对照 06-roadmap.md 的 R0–R8 与 01-review.md 的缺陷清单。每个已完成条目注明 commit。

## 已完成

| 项 | 内容 | commit | 验证 |
|---|---|---|---|
| R0 D1 | 真正的票据：reader 取令牌、reducer 按序交付后归还；重排环容量 = 票据数 | 28bc26d | 原 `output_is_invariant_across_thread_counts` 转绿，重复 20 次通过；新增 chunk 0 延迟 500 ms 回归测试 |
| R0 D2 | BGZF 空块（ISIZE=0，含 EOF 标记）一律跳过，只在底层 EOF 结束；缺 EOF 标记警告一次 | 28bc26d | 拼接文件 fixture + 全部块偏移穷举插入 EOF 测试 |
| R0 D11 | 块按 65,280 字节切、记录可跨块；`deflate_block` 超限报错；索引 chunk 终点 = 下一条记录起始（跨块跟踪） | 28bc26d | 3,000 样本、单行 >100 KB 输出经 noodles 读回逐字节一致 |
| R0 D12 | 扩展名推断输出格式（.gz/.bgz 压缩）；`-o x.bcf` 显式报错 | 28bc26d | CLI 测试 |
| R0 D4 | tid 用 u32、未知 contig 报错；`--write-index[=csi\|tbi]` 默认 CSI（depth 6，覆盖 2^32 坐标）；TBI 拒绝 ≥ 2^29 位置 | 28bc26d | 70,000 contig + 位置 6e8 的 CSI 查询测试；TBI 拒绝测试；未声明 contig 报错测试 |
| R0 D14 | GT 等位下标 checked 算术，70000 不再回绕 | 28bc26d | 单元测试 |
| R0 D19 | M4/M5b 验收行补登 results.tsv；此后每步实测入日志 | 28bc26d | — |
| D7 | contig 长度不编造；`--contig-lengths`（.fai）补充；末窗截断缺长度即报错 | a7fdc6a | 缺长度报错 + lengths 文件修复后与声明长度运行逐字节一致 |
| D8 | 位点升序（worker chunk 内 + reducer 跨 chunk）与超长位置校验；同 POS 多记录合法 | a7fdc6a | 乱序/超长 fixture 报错测试 |
| D9 | 群体文件 = pixy 两列无表头（GroupTable::from_tsv）；未知样本报错，`--force-samples` 跳过并警告；clean.3pop.txt 与 make_3pop.py 移除表头 | a7fdc6a | 全量输出与改动前逐字节一致 |
| R1 | 删除 legacy（reed-core bits/window/population/result/estimators、reed-io record//sink/Assembler、reed-filter exec/emitter/二进制/legacy spec+plan、legacy 测试）；HeaderFacts 携带 noodles header；`reed::input::open` 单一入口，stdin（`-`）可用；驱动不再调用 record_header/filter_header | 08de8a1 | 全量对拍 GT/AC/AN/NS md5 = gold；262 测试全绿 |
| R3 §4.1 | 真正的票据（同 R0 D1） | 28bc26d | 同上 |
| R3 §4.2 | 按字节内存预算 = 背压：reader 取 footprint（压缩 × 8）许可、reducer 交付后归还；仅单 chunk 超全预算才报错；驱动预算默认半物理内存（修复 D5） | f708ea0 | `memory_budget_bounds_in_flight_bytes` 在途峰值 ≤ 预算测试；tiny 预算仍报 MemBudgetExceeded |
| R3 §4.3 | run_sequential 删除，单一代码路径；watcher 线程删除，内外共用一个 Abort 通道 | f708ea0 | 全量 stats 输出 56/7/2 线程 md5 一致；filter 正文 56/7 线程一致 |

### R7（部分：π/dxy 已对齐，其余待续）

| 项 | 内容 | 验证 |
|---|---|---|
| 唯一窗口网格（04 §4） | 从位置 1 起 `[1,w]`、`[1+s,w+s]`…；末窗截断到 contig 长度；空窗输出 NA | 手算 fixture 测试（窗口 1/11/21）；pixy 与 vcftools 同为该网格 |
| π / dxy（04 §3） | Σdiffs/Σcomps（π）、Σcross/Σn_A n_B（dxy），缺失数据按 pixy 逐位点计数；整数分量 diffs/comparables/missing 全部输出 | **pixy 2.2.4 bypass 模式逐位一致**（全量 9,501 窗 × 3 群体，最大相对误差 2e-14；子集同） |
| `--variant-only` | 窗口内文件未出现的位点视为单态且全部调用：π 分母加 `(L−|S|)·C(2N,2)`，dxy 加 `(L−|S|)·2N_A·2N_B` | **vcftools `--window-pi` 逐行一致**（子集 60/60 行，按 vcftools 打印精度；全量待 pixy 完成后并入） |
| 不变位点检查（04 §3） | 记录级 `ALT="."` 计数；无不变位点且未声明 `--variant-only`/`--bypass-invariant-check` 时报错（与 pixy 的检查同口径） | 回归测试；all-sites fixture 上默认模式与 pixy 默认模式逐列一致 |
| pixy 列布局（02-cli §5.4） | `PREFIX_pi.txt` / `PREFIX_dxy.txt` / `PREFIX.log`，列名与 pixy 相同，缺失 `NA` | 手算 fixture 测试 + gold/popgen/compare.py |
| CLI 增量（R5 的 stats 部分） | `--pop/--populations`、`--pi`、`--dxy`、`--window/--size`（100k/1M 单位）、`--step`、`--variant-only`、`--bypass-invariant-check`、`-o PREFIX` | 上述对拍 |

对拍工具：`gold/popgen/`（run_pixy.sh、run_vcftools.sh、compare.py、README）。本机 pixy 2.2.4 与
scikit-allel 1.3.13 从仓库上级目录源码安装（`pip install -e pixy scikit-allel` + dask）。

R7 未完成部分（仍属 R7，按 04 §3 逐个落地）：Watterson θ（逐位点 n 的调和数）、Tajima's D
（窗口平均 n、方差一次）、Fst `wc84`/`hudson`（Bhatia 比值之和）、PBS（三对独立 Hudson）、
SFS 投影、多等位位点的逐等位推广（当前多等位记录跳过，本数据集为 0 条）、
`--sfs/--freq/--hardy/--het/--missing/--depth/--summary` 的 plink2/vcftools/bcftools 布局输出。

## 偏差记录（06-roadmap.md §4.6）

1. **§4.4 类型化 EngineError 暂缓**：45 处 `Job(String)` 全部经 Display 渲染，
   当前无消费者检查 source 链，改动为纯机械 churn。待 R2 统一解码时一并落地。
2. **spec/compat.rs（--conf）保留**：删除清单把它排在 R5（新参数语言）一起，
   提前删除会中途切断 gold 管线。R5 落地 `run_reed.sh` 翻译后删除。
3. **reed-io/vcf/gt.rs（tally_gt）保留**：唯一消费者是 stats 驱动，
   R2 的 GenotypeScanner 落地时删除。
4. **R0 性能数字的环境口径**：输入在 /sri（NFS，共享机器），wall 随页缓存
   状态在 13–42 s 波动；user CPU 时间恒定（filter 8m31s，A/B 对 R0 无回归）。
   results.tsv 中注明口径。

## 全量数据验收（2026-10-07，56 线程）

- stats：18,992 窗 × 3 群体，2.6–2.96 s（门 ≤ 3.0 s）
- filter：3,315,604 → 1,768,759，13–17 s 热缓存（门 ≤ 15.6 s）
- GT md5 `e9efcc63…`、AC/AN/NS md5 `09459e65…` = gold；`bgzip -t` 通过
- CSI 区域查询 Chr05:1.0–1.1 Mb = 119 条（与复审一致）
- 确定性：stats 输出 56/7/2 线程 md5 一致；filter 正文 56/7 线程一致
  （文件级 sha 仅差 `##reed_command` 中的 argv，含 `--threads`，与 bcftools 行为一致）

## 未完成（按 06-roadmap.md §3 顺序）

- **R2** 统一记录解码（GtCode/GenotypeScanner/FormatLayoutCache）：修复 D10、D13
  （当前数据为二倍体，D13 无实际影响）；`##reed_ascertainment` 与清单门仍在。
  多等位逐等位计数（04 §3）也依赖这一层。
- **R4** 输出层模块化（BlockWriter/OrderedOutput 移入 reed-io/src/write）：
  行为已由 R0 的按字节切块 + 跨块索引跟踪实现，模块搬家未做。
- **R5** CLI 重写：reed-cli crate、`--pi/--dxy/--fst/--theta/--tajima-d/--pbs/--sfs/
  --freq/--hardy/--het/--missing/--depth/--summary`、`--pop/--window/--step/
  --variant-only`、bcftools 表达式子集、删除 --spec/--conf/TOML/REED_THREADS。
- **R6** 输入/输出补全：未压缩 VCF、普通 gzip、BCF 输入；`-r/-R` 区域；BCF 输出。
- **R7** 统计对齐：按 04 §3 唯一定义重写估计器（pixy 缺失数据处理、多等位推广、
  窗口网格从 1 起、`--variant-only` 分母补全、不变位点检查）。注意：本机有
  vcftools/plink2/bcftools，无 pixy/easySFS/scikit-allel，04 §5 对拍表中只有
  vcftools/plink2 的行可本地执行。
- **R8** 更多参考实现（LD、f 统计、联合 SFS）。
