# 02 GT 语义域（V25-01 / V25-02 / V25-03）

> 一个语义一个入口：allele 的合法域、GT 子字段的位置、一行的合法性，
> 三者都在 reed-io 定义一次；scanner / filter / encoder / convert / pca 只消费。

---

## 1. V25-01：GT allele 域校验

### 现状【据】

- `genotype_scan.rs` 的 `parse_allele`（:426 起）只校验：token 为 `.`、空/长度 > 5、
  非数字、u16 溢出——**全文 0 处 `n_alt`**；
- encode 侧 `push_gt_code`（`encode_native.rs:1342-1355`）只界到 i32
  （"a = 1073741822; … doc_v10 D10-1"），不对照记录的 ALT 数；
- `gt_codes_into`（`encode_native.rs:1300`）只查语法（分隔符/相位），不查范围；
- PCA 侧 `PcaJob::dosage`（`crates/reed/src/genotype.rs:81-110`）在 `:106`
  把 `a/b > n_alt` 返回 `None`，调用点 `unwrap_or(f64::NAN)`——越界 allele
  静默当缺失，不是 `MalformedRecord`；`run_chunk`/`run_seam` 还用首字段
  启发式找 GT（§2）。

### 触发 / 后果

输入 `chr1 1 . A G . PASS . GT 0/9`（ALT 只有 1 个）：

- filter 路径照常处理（多等位按 CollapseAlt 当作 ALT）；
- `-Ob` 输出：allele 9 被编为 `(9+1)<<1 | bit`——htslib 读回该记录报 allele 越界，
  下游工具拒绝或误读；
- `tracking_alleles` 的逐 allele 统计（π per-allele、SFS）混入不存在的 allele 计数。

scan / encode / pca 没有显式 `n_alt` 校验；convert 的 `parse_gts` 也没有，
但它的下游行为按目标分叉：`phylo_column` / `resolve_msa_value` 用
`allele_ranges.get(index)` 间接拒绝越界（计 malformed），`--012` 的
`cell_012` 对 `n==2` 的越界返回 `None`、对 `n<2` 直接返回 `-1`（缺失）。
同一越界 allele 三种结果——正是要把域检查收进 `decode_gt` 的理由。

PCA 侧是第四条路径：同一记录在 `reed pca` 中不报错，越界 allele 被
`PcaJob::dosage` 折成 NaN；若同一列的 GT 位置也错（`DP:GT`），整个位点会被
`any_alt=false` 静默跳过。`filter`/`convert` 报错、`pca` 静默，口径分叉。

### 规范【规】

1. **scanner**：`GenotypeScanner::scan` 已持有 site 上下文；`parse_allele` 增加上界
   （`n_alt: u32`），合法域是 `0..=n_alt`（0 = REF；`ALT=.` 时 `n_alt = 0`），
   `token > n_alt` 返回 `IoError::MalformedRecord`
   （消息含 locus 与 `GT allele index N exceeds ALT count M`）。文本路径与 BCF
   快表路径同规则：BCF 侧的存储值是 `(allele+1)<<1|phased`，`decode_pair` 先
   解出 allele 再对照记录的 `n_allele`。
2. **encoder**：`gt_codes_into` / `gt_shape` 的调用点（`encode_indiv`）持有记录的
   ALT 数，越界同样经 `bad()` 报 `MalformedRecord`。
3. **同居**：域检查与 GT 语法解码（§2 的公开 API）放在一处，不散落第三份；
4. **PCA**：`PcaJob::run_chunk`/`run_seam` 的 GT 位置改走 §2 的公开 API；
   `dosage` 的越界分支改成返回错误（由 job 转 `MalformedRecord`）；多等位
   skip 保留（`allow_multiallelic=false` 是策略，不是域校验）。

### 测试

- 文本路径：`0/9`、`1|2`（ALT=G 时 2 越界）、`0/1`（合法对照）三 fixture；
- BCF 快表路径同三 fixture（`bcf_vs_noodles` 差分套件加 case）；
- PCA 路径：同三 fixture 走 `run_chunk` 与 `run_seam` 各一遍；`0/9` 报错，
  `0/1` 数值不变；
- oracle golden：沿 V24-01 的 `02-oracle-contract.md` 追加"越界即报错"契约
  （oracle 是纯计数库，不读 GT 文本；契约约束的是喂给它的 GT 域）；
- e2e：`reed filter`（文本路径与 `-Ob` 编码路径各一条）与 `reed pca` 对 fixture
  报错，退出码 1；`reed convert` 的 row targets 按 B13-9 政策把越界 GT 计
  malformed 并跳过（exit 0，日志计数），不新增第二套错误政策。

### 验收

- 四条路径（scan / encode / convert / pca）对越界 allele 都 `MalformedRecord`；
- `grep -c n_alt crates/reed-io/src/vcf/genotype_scan.rs` ≥ 1（检查点存在）；
- `grep '> n_alt' crates/reed/src/genotype.rs` 只出现在 `decode_gt` 的调用点或
  错误分支，不再出现 `unwrap_or(f64::NAN)` 吞掉域错误；
- 既有全量测试（945 基线）不回归。

---

## 2. V25-02：GT 位置唯一入口

### 现状

同一语义——"GT 在 FORMAT 的第几列"——**6 份**实现（PCA 的 body/seam 各一份副本）：

| # | 位置 | 实现 | 状态 |
|---|---|---|---|
| 1 | `reed-io/src/vcf/format_layout.rs`：`Layout`（V24-06 后经 `LayoutId{epoch,index}` 寻址，`get(id).gt_index`） | 缓存内布局 | 【据】存在——权威候选 |
| 2 | `reed/src/filter/record_pipeline.rs:735` `gt_subfield_index` | 自写查找 | 【述】 |
| 3 | `reed-io/src/bcf/encode_native.rs:446-447` | 自行 split FORMAT | 【述】 |
| 4 | `reed/src/convert/job.rs:954-955` `parse_gts` | `memchr(b':')` 取首个子字段当 GT | 【据】——**错的** |
| 5 | `reed/src/genotype.rs:171-174`（`PcaJob::run_chunk`） | `position(|&c| c == b':')` 取首个子字段 | 【据】——**错的** |
| 6 | `reed/src/genotype.rs:247-248`（`PcaJob::run_seam`） | 同上（seam 副本） | 【据】——**错的** |

### 触发 / 后果【据】

`FORMAT=DP:GT` 的合法 VCF（GT 不在首位）经 `reed convert`：

- `parse_gts` 把 DP 值当 GT：`30:0/1` → raw=`30` → `n=1`、presence=true；
  `--012` 的 `cell_012` 对 `n=1` 返回 `-1`（半调用/缺失），**静默输出缺失**；
  phylo/msa 路径随后 `allele_ranges.get(30)` 越界 → 计 malformed/跳过。
  同一格式错误在不同目标上两种结果，presence 却照记；
- `reed pca` 的 `run_chunk`/`run_seam` 同样把 DP 当 GT：`dosage("30", n_alt)`
  返回 None → NaN；全样本 NaN → `any_alt=false` → 位点**静默跳过**；全文件
  `DP:GT` 时 PCA 在零位点上运行（或报空矩阵），而不是报格式错误；
- 补丁史实证：V22-R2-1（u8 下标 bug）当年一个 commit 修了四处中的三处，
  commit message 自述 "filtered.rs's duplicate GT lookup reuses…"。

### 规范【规】

1. `Layout`（含 `gt_index`）与 GT 语法解码提升为 reed-io 公开 API（签名示例，
   实现者可调整）：

   ```rust
   /// GT 子字段的位置；FORMAT 无 GT 时为 None。
   pub fn gt_subfield<'a>(sample_col: &'a [u8], layout: &Layout) -> Option<&'a [u8]>;
   /// 一次 GT 解码 + §1 的域检查。
   pub fn decode_gt(sub: &[u8], n_alt: u16) -> Result<GtTokens, IoError>;
   ```

2. 六个消费点全部改走它：scanner（已是）、filter `record_pipeline`、encoder
   `encode_native`、convert `parse_gts`、PCA `run_chunk`、PCA `run_seam`；
   `parse_gts` 删除 `memchr(b':')` 启发式，PCA 两处删除
   `position(|&c| c == b':')`，全部改用 `layout.gt_index` 取子字段。
3. 【荐】`validate_sample_cols` 同批提升为公开 API，encoder/convert 的样本列走查
   不再各自静默填 ABSENT / `begin.min(len)`（01 §1.1）。

### 测试

- `DP:GT`、`GT:DP`、`AD:DP:GQ:GT`（GT 在末位）三 fixture 走 convert / filter /
  `-Ob`，三命令输出的 GT 一致；
- 同三 fixture 走 `reed pca`：`run_chunk`（body）与 `run_seam`（人为切在记录
  中间）两条路径的 dosage 矩阵一致，`DP:GT` 与 `GT:DP` 相同；
- 既有 V22-R2-1 回归用例保持绿。

### 验收

- GT 位置实现点 **6→1**：grep 自写查找（`gt_subfield_index`、encode 侧 split、
  convert 侧 `memchr(b':')`、PCA 两处 `position(|&c| c == b':')`）清零；
- `DP:GT` fixture 的 convert 输出与 `GT:DP` 等价（同 GT、同 presence）；
- `reed pca` 对 `DP:GT` 的 `.eigenvec` 与 `GT:DP` 逐字节一致，不再静默丢位点。

---

## 3. V25-03：RecordClass 行合法性单点

### 现状

"这行是不是合法记录"没有单一实现，五个入口各自从 `LineRef::parse +
validate_sample_cols + parse_pos1` 重新推导，策略还互相矛盾：

| 入口 | 判定 | Blank 策略 | Malformed 策略 |
|---|---|---|---|
| convert row（`job.rs:733-741`【述】） | 自写 | 跳过、不计 | 计 malformed 继续 |
| convert `line_is_sane`（`:1352-1360`【据】） | parse + 列数 + POS | **返回 true**（:1353-1354） | false |
| stats `push_line`/`run_seam`（`stats/run.rs:621-622/681-683`【据】） | `LineRef::parse` | — | **硬报错** |
| rename（`rename/mod.rs:1185`【述】） | 无（字节透传） | 透传 | 透传 |
| filter（`record_pipeline.rs:169-171`【述】） | `LineRef::parse` | — | 报错 |

### 触发 / 后果【据】

- 空行恰好横跨 chunk 边界：engine 的 seam 交出 `b"\n"` → `line_is_sane` 判 sane
  （:1353-1354）→ `run_seam` 喂给 stats（`:759-766`）→ `LineRef::parse(b"\n")`
  失败 → `EngineError::job("record line has fewer than 8 columns")` →
  **`reed convert --sfs … 整体失败**；同一文件去掉 `--sfs` 正常完成；
- CRLF 文件的空行（body 中 `b"\r"`）在任何 chunk 位置触发同一失败；
- 是否触发取决于 `--memory/--threads` 决定的 chunk 边界——**间歇性复现**，
  最难排查的一类。

### 规范【规】

1. reed-io 提供一次分类：

   ```rust
   pub enum RecordClass<'a> {
       Record(LineRef<'a>),
       Blank,                    // "" / "\n" / "\r" / "\r\n"
       Malformed(MalformedLine), // parse / 列数 / POS 任一失败
   }
   pub enum MalformedLine { FewerThanEightColumns, SampleColumns(SampleColumnsError), Position }
   pub fn classify_line(line: &[u8], n_samples: usize) -> RecordClass<'_>;
   /// chunk body（以 \n 结尾的完整行序列）里真正的行：丢掉末尾终止符后的空片，
   /// 保留正文中的真空行。
   pub fn body_lines(body: &[u8]) -> impl Iterator<Item = &[u8]>;
   ```

   `MalformedLine` 携带失败项而不是只给一个 unit，消费点不必二次 parse 就能
   保留原有诊断（convert 的一次性 warning、stats 的致命错误、filter 的错误）。
   POS 0 归入 `Position`（VCF 1-based；convert 行路径原本也按 malformed 处理）。

2. `line_is_sane` / `sanitize_body` 变成 `matches!(classify_line(..), Record(_))`；
   五个入口只声明对 Blank / Malformed 的策略，不重写判别；
3. **策略表成文（本表即契约，2026-11-04 按实现校正 rename 的 Blank 格）**：

   | 入口 | Blank | Malformed |
   |---|---|---|
   | convert row | 跳过（不计 malformed、不占 draw slot） | 计 malformed 继续（现状） |
   | convert SFS | 跳过 | 跳过（与 row 的统计口径对齐；纯 SFS 与 mixed 一致） |
   | stats（popgen `StatsJob` 与 `stats` 的 `ReportJob`） | 跳过 | **硬报错**（严格口径不变） |
   | rename | **跳过**（非记录；body 与 seam 一致） | 透传（字节编辑器不丢非空行） |
   | filter | 报错 | 报错 |
   | pca / distance（`PcaJob`） | 跳过（本轮补齐：`run_chunk`/`run_seam` 接 `classify_line`） | 报错 |

   **convert SFS 的第二道过滤（2026-11-05 复核补）**：结构合法的记录若 GT 超出
   convert row 的接受范围（越界 allele、非数字 token、polyploid、五位数以上的
   token），也不进 SFS——严格 stats 扫描器会把它们判为致命错误，而 convert 的
   政策是 exit 0（B13-9）。全站都被跳过时写出全零谱并 exit 0（`SfsPartial` 是
   惰性初始化的，不能因此 panic）。

   stats 对 Malformed 的严格性保留——本轮只修 Blank 的分叉，不放宽任何口径。
   rename 的 Blank 从"透传"改为"跳过"：body 原本就丢 LF 空行，而 seam 会多输出
   一条空行、CRLF 空行则被当作字节原样写出——三处不一致正是本节要消灭的分叉；
   透传需要给空行编造一个索引 key（`--rm-dup-ids` 的 key 流还会错位），而空行
   本就不是记录。非空 malformed 行仍逐字节透传，`--rm-dup-ids` 的 raw 流与 key
   流保持一一对应。

### 测试

- e2e：按压缩块构造 1 MiB chunk 目标边界上的空行 fixture（`seam_blank_vcf`：
  逐块 deflate、累计压缩字节到目标、下一块正好是空行）+ `--sfs`；LF 与 CRLF
  两种 seam 空行，谱与无空行 fixture 逐字节一致；
- CRLF 空行落在 body + `--sfs`；
- filter 对 LF / CRLF 空行都报 `blank record line`；
- rename：body 的 LF/CRLF 空行与 seam 空行都消失，记录数不变；
- stats：`run_seam(b"\n"/b"\r\n")` 跳过、body 的 CRLF 空行跳过、malformed 仍致命；
- 既有 B13-2 sanitize 用例保持绿。

### 验收

- 空行 / CRLF 空行在任何 chunk 位置都不再使 `convert --sfs` 失败；
- `line_is_sane` 的 blanket-true 分支消失：grep `b == b'\n' || b == b'\r'`
  只出现在 `classify_line`；
- 五入口对 Blank/Malformed 的行为与上表逐格一致（单测钉住）。

### 实施记录（2026-11-04）

- **跨 seam 是否复现由输入字节决定，不是 `--memory/--threads`**（现状段原文的
  修正）：chunk 边界只由压缩块与 1 MiB 目标决定；测试因此直接构造边界
  （`e2e_convert::seam_blank_vcf`），而不是靠预算参数碰运气。
- 首版 e2e 的两个新测试在实现前按预期失败：`convert --sfs` 报
  `Error(Engine(Job("record line has fewer than 8 columns")))`（seam 与 CRLF
  body 各一）。
- 实现后发现一个额外分叉：`convert` 的 SFS 目标在 mixed 运行会 sanitize，纯
  SFS 却把 body 原样交给 stats——同一条 malformed 行在两种组合下"跳过 vs 硬
  报错"。按表统一为 SFS 永远 sanitize（干净 body 零拷贝返回，字节不变）。
- 提交：`0cef99f`（语义，969 passed）、`9b61a37`（`filter/mod` pass1 合一 +
  `bcf_file_header` 下沉，970 passed）。
