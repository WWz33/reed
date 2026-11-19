# 04 复核：v30 结论的独立核对 + 本轮新发现

> 2026-11-11 晚。对象：本目录 01–03 与 STATUS 的结论。
> 方法：用当前二进制 `reed 0.4.0 (v0.3.0-180-gf190708-dirty 2026-11-11)`（`/tmp/reed-target/release/reed`，
> 晚于所有 `crates/**/*.rs` 的修改时间）实跑，并与参照工具的**二进制输出**和**源码**逐项比对：
> plink 1.9（v1.9.0-b.7.11）、plink2、bcftools、vcftools、scikit-allel、ADMIXTOOLS（`qpDstat.c`、`qp3Pop.c`）。
> 数据：`~/reed/subvcf.vcf.gz`（348 样本、3000 位点、单 contig Chr05）+ 本文给出的小夹具。临时产物在 `/tmp/v30chk/`。
> **本文不改代码。** 结论分三类：v30 对且修法对；v30 对但**修法错**；v30 漏掉的。

---

## 1. 结论先行

| 类别 | 条目 |
|---|---|
| **修法错（照做会破坏参照对拍）** | V30-01 的 `e02` 公式；V30-04 的"Σ 尺度对齐 ADMIXTOOLS"依据；V30-R02 的"`.afreq` 与 `.stats` 同源"验收门 |
| **定性需改** | V30-R01 不是"语义分叉"而是 convert 与 vcf2phylip/vcftools 的两个口径同时存在；V30-R12 的规则在 v1 就被 v2 推翻过；V30-Q29/Q16/Q24/T01 的前提与 v5 08 §6 的产品裁决相反 |
| **新发现（v30 漏掉）** | **N-1** `--skip-empty` 截断窗口（**数值错误**，P0）；**N-2** PI_HAT 计入单态位点的口径与 plink 不同；**N-3** `.genome` 的 DST/PPC/RATIO 恒写 0；**N-4** 本机 plink 1.9 早于上游修复，`ibs` 对拍基线需固定版本；**N-5** `--summary` 缺 PSI/HWE 两节且未写明 |
| **复核成立** | V30-01（现象）、V30-02、V30-03、V30-05、V30-07、V30-R02（溢出）、V30-R09、V30-Q07、V30-Q31/T07、V30-Q34、T11 |

---

## 2. 修法错的三条

### 2.1 V30-01：现象成立，但给出的"正确系数"里 `e02` 是错的

**现象复核（成立）**：`reed distance subvcf.vcf.gz --type pihat` 的 60,378 对里 **60,346 对**输出
`Z0=0 Z1=1 Z2=0 PI_HAT=0.5`；plink 1.9 `--genome` 对同一输入的 PI_HAT 分布是
min 0.4676 / 中位 0.5899 / max 0.7287（这批样本本来就高度相关，"无关对 PI_HAT≈0"的验收门在这份数据上不成立）。

**修法错误**：01 §2.3 与 02 的 V30-01 都写 `e02 = (p²+q²)²`。按 `plink-ng/1.9/plink_calc.c:4865` 取无限样本极限
（`dxx1,dxx2,dyy1,dyy2,(x−3)/x,(y−3)/y,f2,f3 → 1`）：

```
e02 = q⁴ + p⁴ + 4p²q²        （不是 (p²+q²)² = p⁴ + q⁴ + 2p²q²）
```

独立复算（numpy，全部 60,378 对，与 plink 打印值比较）：

| 系数方案 | 与 plink 1.9 `.genome` 的最大差（4 位小数） |
|---|---|
| plink 有限样本系数（`plink_calc.c:4848-4867`） | **0**（全部 60,378 对逐位相同） |
| 无限样本极限，`e02 = p⁴+q⁴+4p²q²` | 8.1e-4（每一对都差在第 4 位） |
| v30 给的 `e02 = (p²+q²)²` | 首对 PI_HAT 0.6408 vs plink 0.6316（差 9.2e-3） |

所以 v30 路线里的"或至少在无限样本极限下……"两条都不能对上 plink 的输出；**唯一能逐位对拍的是完整复刻有限样本系数**，
它只多 10 行算术，不需要"复刻 genome.cpp 的推导"（v17 STATUS §4 的顾虑不成立）。

**还要同时满足的口径**（复算时逐项确认，缺一个就对不上）：

1. 系数只累加 `num_alleles > 3 && 0 < p < 1` 的位点，并除以**这些位点的个数**（`ibd_prect`），不是除以 `n_sites`；
2. `num_alleles` 取该位点的已调用等位数（`2 × (n − 缺失)`），每位点不同；
3. 每对的 `nn` 是**两人都有调用的全部位点数，包括单态位点**（`plink_calc.c:4381`）——见 N-2；
4. 求解与裁剪顺序已与 `plink_calc.c:4383-4418` 同构（`clip_ibd`），不用改。

**代码示意**（替换 `crates/reed-distance/src/matrices.rs:568-594` 的 `expected_ibs_coefficients`）：

```rust
/// plink 1.9 `--genome` 的期望 IBS 系数（plink_calc.c:4848-4867，含 PLINK 1.07 的有限样本修正）。
/// 只有 `num_alleles > 3` 且 `0 < p < 1` 的位点参与，均值的分母是参与位点数（plink 的 `ibd_prect`）。
fn expected_ibs_coefficients(st: &Standardized) -> [f64; 5] {
    let n = st.n_samples();
    let (mut e, mut used) = ([0.0f64; 5], 0u64);
    for (s, &p) in st.freq.iter().enumerate() {
        let called = (0..n).filter(|&i| st.is_called(i, s)).count();
        let na = 2.0 * called as f64;
        let q = 1.0 - p;
        if !(na > 3.0 && p > 0.0 && q > 0.0) {
            continue;
        }
        let f2 = na * na / ((na - 1.0) * (na - 2.0));
        let f3 = f2 * na / (na - 3.0);
        let (x, y) = (p * na, q * na);
        let (x1, y1) = ((x - 1.0) / x, (y - 1.0) / y);
        let (x2, y2) = (x1 * (x - 2.0) / x, y1 * (y - 2.0) / y);
        e[0] += 2.0 * p * p * q * q * x1 * y1 * f3;
        e[1] += 4.0 * p * q * f3 * (p * p * x2 + q * q * y2);
        e[2] += f3 * (q.powi(4) * y2 * (y - 3.0) / y + p.powi(4) * x2 * (x - 3.0) / x + 4.0 * p * p * q * q * x1 * y1);
        e[3] += 2.0 * p * q * f2 * (p * x1 + q * y1);
        e[4] += f2 * (p.powi(3) * x2 + q.powi(3) * y2 + p * p * q * x1 + p * q * q * y1);
        used += 1;
    }
    if used == 0 { return [f64::NAN; 5]; }               // plink 此时无法给出 Z；输出 NA 而不是 1.0 哨兵
    e.map(|v| v / used as f64)
}
```

（`called` 的逐位点计数可以在 `standardize` 里顺手存下，避免 O(n·m) 重算；`p` 用 REF 还是 ALT 频率不影响结果，五个式子对 `p↔q` 对称。）

**验收门（替换 02 的 V30-01 门）**：

1. `gold/` 新增 `run_genome_parity.sh`：plink 1.9 `--genome full` vs `reed distance --type pihat`，
   Z0/Z1/Z2/PI_HAT 按 4 位小数**逐对相等**（本轮 numpy 复算已证明可达 0 差）；
2. 夹具必须含单态位点（见 N-2）和缺失；
3. 删除"无关对 PI_HAT≈0"这条门：它依赖数据，不是实现的性质。

### 2.2 V30-04：estimate = Σ 是事实，但"对齐 ADMIXTOOLS 的 sum 尺度"不成立

**事实（成立）**：`crates/reed/src/stats/request.rs:700-733` 对 f2/f3/f4 输出 `Σnum`。
本轮实跑 `--f2 P1,P2 --f3 'P3;P1,P2'`：estimate = −3.069 / −1.316（3000 个位点的和）。

**依据不成立**：代码注释（`request.rs:701-703`）、v7 D7-9、v30 §2.2-4 都说"ADMIXTOOLS 报的是位点和"。源码不是这样：

| 工具 | 估计量 | 证据 |
|---|---|---|
| qpDstat `f4mode: YES` | `Σ ztop / Σ zbot`，`zbot += swt`（每位点 +1） | `AdmixTools/src/qpDstat.c:1170-1174`、`:1210-1213`：**位点均值** |
| qp3Pop | `Σ btop / Σ bbot` 的比值 | `AdmixTools/src/qp3Pop.c:828-829` |
| ADMIXTOOLS 2 (R) | 每块 `rowMeans`（块内位点均值） | `admixtools/src/cpp_fstats.cpp:295,323-388` |
| scikit-allel `average_patterson_f3(normed=False)` | `np.nanmean(T)` | `scikit-allel/allel/stats/admixture.py:353` |

reed 的 f2/f3/f4 对拍脚本（`gold/popgen/compare_fstats.py:94-102`）比的是 `np.nansum(...)`，
即"自己定义 Σ，再拿 allel 的逐位点值求和来对"——这证明了代数一致，没有证明与任何工具的**输出尺度**一致。
ADMIXTOOLS 对拍（`run_admixtools_jackknife.sh`）只覆盖 D/比值类，Σ 与均值在比值里相消，所以从没暴露。

**后果**：f2/f3/f4 的 estimate 和 se 随位点数线性增长；同一对群体换一个 SNP 面板，数值不可比；
用户拿 0.0x 量级的文献值去比会差两三个数量级。z 不受影响。

**建议（改输出，不只是改文档）**：estimate 与 se 都改为**位点均值尺度**（`Σnum / n_snps`，jackknife 用
`jackknife_shared`，den 已经是逐位点 1.0，所以只需把 f2/f3/f4 分支并入比值分支）：

```rust
// crates/reed/src/stats/request.rs FStatAcc::row —— f2/f3/f4 与 D/f4-ratio 走同一个比值 jackknife。
// den 对加性统计量逐位点为 1.0，所以 num/den = 位点均值，与 qpDstat f4mode、ADMIXTOOLS 2、allel 同尺度。
let den = self.den.value();
if den == 0.0 { return FStatRow::na(self.test.clone(), self.n_sites); }
let (_, se, n_blocks) = self.jk.jackknife_shared(grid).unwrap_or((self.num.value() / den, f64::NAN, 0));
(self.num.value() / den, se, n_blocks)
```

门：`compare_fstats.py` 改比 `np.nanmean`（f2/f4）与 `average_patterson_f3(normed=False)`；
qpDstat `f4mode: YES` 加入 ADMIXTOOLS 对拍；`estimate/se == z` 保持。这是输出格式变化，按 `cli.md` §0 的稳定性承诺要在 CHANGELOG 记一笔。
如果决定保留 Σ，则至少把 `request.rs:701-703` 与 v7 D7-9 的"对齐 ADMIXTOOLS"改成"reed 自定义 Σ 尺度"，并在 `popg.md` 写明。

### 2.3 V30-R02：溢出是真问题；"`.afreq` 与 `.stats` 同源"不是修法

**成立的部分，且比 v30 写的更宽**：问题不只是求和溢出，而是 reed 用 `u64` 读 INFO/AC，htslib 用 int32。
htslib 解析 VCF 时把超出 int32 的 INFO 整数**置为缺失**（`htslib/vcf.c:4037-4049`，打印 `Extreme INFO/AC value … set to missing`），
`bcf_calc_ac`（`htslib/vcfutils.c:39-75`）于是读不到 AC，改走 GT。reed 的 `summary.rs:299-328` 直接 `parse::<u64>()`：

| 输入（2 样本 `0/1 1/1` 或 `0/1 0/2`） | `bcftools stats -s -` 的 AF 段 | `reed stats --summary` |
|---|---|---|
| `A>G`，`AC=18446744073709551615;AN=100`（单 ALT，**无求和溢出**） | 0.000000 一行（走 GT） | **0.740000**（u64 解析成功后按 AC/AN） |
| `A>G,T`，`AC=18446744073709551615,1;AN=100` | 0.000000 两个 ALT 同一行 | 0.000000 + 0.990000 两行（求和在 release 回绕） |
| 同上但不写 AC/AN | 0.000000 | 0.000000（一致） |
| `AC=5;AN=8`，GT 全 0/0 | 0.610000 | 0.610000（一致） |

debug 构建在第二行 panic。修法：按 htslib 的规则把超出 `i32` 的 AC/AN 当作缺失（不是 `checked_add`），然后回退到 GT 计数；
filter 的 `info_af_indexed`（`reed-filter/src/semantics.rs:195-271`）当前对同一输入**拒绝**，两条路径的判定要统一到 htslib 的"置缺失"。

**不成立的部分**：v30 把"同一命令 `.afreq` ALT_FREQS=0 而 `.stats` AF=0.61"列为缺陷，并要求两者同源。实测：

```
INFO AC=5;AN=8，GT 全 0/0
reed .afreq   ALT_FREQS 0          ← plink2 --freq 永远按 GT 算
reed .stats   AF 0.610000
bcftools stats AF 0.610000         ← bcftools stats 优先读 INFO/AC,AN（vcfstats.c:673 传 BCF_UN_INFO）
```

两个输出各自与自己的参照工具逐字一致。让它们"同源"会让其中一个偏离参照。验收门应改为：
上表四个夹具在 `reed stats --summary` 与 `bcftools stats -s -` 的 AF 段逐行相同；filter 与 stats 对超范围 AC 的判定相同；
**`.afreq` 与 `.stats` 的 AF 来源差异写进文档**，不合并。

---

## 3. 定性需改的条目

### 3.1 V30-R01：convert 的"有数据"判定是 vcf2phylip 的口径，`--012` 的半调用编码是另一个问题

`gt.present = !raw.starts_with(b".")`（`crates/reed/src/convert/job.rs:851`）逐字对应
`vcf2phylip/vcf2phylip.py:102-108` 的 `num_genotypes`（`record[i].startswith(".")` 计为缺失）。
所以 `--phylip --min-samples-locus` 上 `./1` 与 `0/.` 判定不同，是**参照工具本身**的不对称，不是 reed 内部分叉。

真正的偏差在 `--012`：

```
half.vcf：位点 10 = ./1 ×4 + ./. ×4；位点 20 = 0/. ×4 + ./. ×4
reed convert --012      全部样本两位点都写 -1
vcftools --012          前 4 个样本两位点都写 1（variant_file_format_convert.cpp:461-471：只有两个等位都缺才写 -1）
```

`gold/popgen/run_convert_parity.sh:20-21,87-88` 已经知道这一点：它**删掉半调用那一行再对拍**，并断言"reed −1 vs vcftools 1"。
所以这是一个有意的、但用户文档（`convert.md:26-27`）没写的偏离。处理：二选一——对齐 vcftools，或在 `convert.md` 写明"半调用写 −1（vcftools 写 1）"。
R01 建议的 `has_called_allele()` 统一会让 `--phylip` 偏离 vcf2phylip，**不要这样做**。

ALT 计数"5 份"成立（`summary.rs:86-92`、`reports.rs:634-639`、`convert/job.rs:824-829`、`reed-filter/src/eval/head.rs:249-254`、
`genotype.rs:263-269`），但五份对 `.` 的处理一致，只有 `genotype.rs` 的空 ALT 语义不同；这是维护债，不是数值错误。
四份都在能依赖 `reed-io` 的位置，直接调 `reed_io::vcf::site::count_alts` 即可，不需要新 `AltCount` 类型。

### 3.2 V30-R12：`architecture.md:1843` 是 v1 规则，v2 已经改过分层

`docs/spec/architecture.md:168-196`（v1，2026-10-03）把 `reed-filter` 放 L3 且"只依赖 reed-core"。
`docs/spec/filter.md:50-62`（v2，2026-10-05，"按建成现状修订"）把**传输**定为 L0 `reed-io`，`reed-filter` 的 L2 计划层建在它上面。
`reed-filter` 实际只用 9 个 `reed-io` 项，全是 VCF 文本视图与 INFO 解析：

```
InfoIndex ×7、parse_info_number ×4、count_alts ×3、LineRef ×3、parse_qual ×2、
InfoValueOut、info_fragments、IntSpec、IoError
```

这是"规格之间的矛盾"，不是"代码违规"。要么按 v2 改 `architecture.md` §3.2/§14.6（承认 filter 依赖 io 的文本视图），
要么做 R12 的下沉——后者要把 INFO 片段索引搬进 `reed-core`，成本远大于收益。**推荐前者**，并把 R12 从 P0 降为文档项。

### 3.3 V30-Q29/Q16/Q24/T01/T04：前提与 v5 的产品裁决相反

`docs/rounds/v5/08-api-surface.md:37,53` 与 `docs/spec/cli.md:26` 已定："`reed` 与各 `reed-*` crate 都是内部结构，`publish = false`，
Rust 公开项不承诺兼容；稳定性只在 CLI 与输出格式"。v30 第二、四部分把"不能 `cargo add reed`""semver 只查 3/11""零用户 doctest"
记为 P0，前提是"reed 是外部 API 库"。这是要重新裁决的产品问题（03 的根因 F），在裁决前不应作为 P0 缺陷排期。

但是有两处矛盾必须现在修，因为它们是**对外文档自己的说法**：

1. `docsite/src/content/docs/library.md:5-12` 写"The `reed` crate is the public API""The stability promise covers the crate root"，
   与 `cli.md:26` 冲突（该页 2026-10-30 随 docsite 迁移加入，没有对应裁决）；
2. `reed-distance`/`reed-reduce`/`reed-graph` 没写 `publish = false`，与 `cli.md:26`"各 `reed-*` crate 都 `publish = false`"冲突。

二者都是改一行的事：裁决前让 `library.md` 与 `cli.md` 一致，三个 crate 补 `publish = false`。

### 3.4 其余 V30 条目的复核

| 条目 | 结论 | 本轮证据 |
|---|---|---|
| V30-02 `-I` | **成立**，且更宽：默认 `-x` 在**无任何样本变化**时也删 AC/AN/AF/MLEAC/MLEAF/ExcessHet/InbreedingCoeff/DP（`reed-filter/src/spec/recompute.rs:40-54`）。实测 `reed filter withAC.vcf -Ov` 与 `--maf 0.05` 都输出 0 条带 AC 的记录。这是 gold 管线的约定（`cli.md:141`、`--help` 都写了），文档自洽；**只有 `docsite/filter.md:155` 的 `-I` 一句是错的**。修法 (a) 让 `-I` 保留 AC/AN 与 bcftools `view -I` 一致，推荐 | `/tmp/v30chk` 六组参数 |
| V30-03 多等位 | **成立**，但 `--summary` 不受影响：它对多等位与 `bcftools stats` 的 SN/AF/SiS 逐行一致。`.afreq` 对 `A G,T` 写单值 `0.5`，plink2 写 `0.25,0.25`；`.hardy` 只有一行，plink2 每个等位一行（共 3 行） | multi.vcf |
| V30-05 `--f4-ratio` 角色序 | **成立**：`translate/popg.rs:671-683` 是 `f4(A,O;X,C)/f4(A,O;B,C)`，O = 第二个名字；`popg.md:39,245` 只写 `A,B;C,D;X` | 源码 |
| V30-06 `f3*` | **成立**：`popg.md:242` 承诺 ΣT/ΣB，CLI 无此输出（`--help` 无 normalize 选项） | `popg --help` |
| V30-R09 `--memory` | **成立**：`filter/mod.rs:357-361` 用 `default_half_memory()`。实测 12 万记录、草图上界 27.5 MiB 时 `--memory 16` 仍选 `sketch_then_emit`；`--memory 1` 在更早的引擎许可处直接报错 | huge.vcf |
| V30-Q07 死字段 | **成立**：`pca/mod.rs:258-260` 的 `threads/memory` 在 `pca/`、`dapc/` 内零读取 | grep |
| V30-Q31/T07 | **成立**：`cargo check -p reed-reduce --no-default-features` → E0063 `lanczos.rs:146` | 本轮实跑 |
| V30-Q34 | **成立**：`git check-ignore` 命中 `crates/reed/.gitignore:2`、`crates/reed-core/.gitignore:1`；本目录新文件也被根 `.gitignore:19 /docs/` 隐藏 | 本轮实跑 |
| V30-R06 | **成立**：plain VCF 加 `-r` 报 `malformed BGZF data`，exit 1 | multi.vcf |
| T11 | **成立**：`convert/job.rs` 2105 行 | wc |

---

## 4. v30 漏掉的缺陷

### N-1（P0，数值错误）`--skip-empty` 把还在累积的窗口提前输出，π/θ/dxy 被截断

**复现**（10 位点，`chr1:1000-5000` 与 `chr1:205000-209000` 两簇，`--window 10k`）：

```
                       window          avg_pi        no_sites
不加 --skip-empty      200001-210000   0.00028571    5
加   --skip-empty      200001-210000   0.00005714    1      ← 只算了第一个位点
```

真实数据：`subvcf.vcf.gz --pi --variant-only --window 100k --skip-empty`，第一个非空窗口 `16000001-16100000`
输出 `no_sites=1`、π 4.995e-06；不加 `--skip-empty` 是 124 个位点、π 3.330e-04（**差 67 倍**）。
与 vcftools `--window-pi 100000 --window-pi-step 50000` 对比，21 个窗口里 20 个一致（最大相对差 < 1e-5），
唯一不一致的就是这个被截断的窗口（vcftools 29 位点，reed 1）。`--theta` 同样截断。1/4/56 线程、BGZF 与纯文本输入结果相同（确定性地错）。
contig/genome 层不受影响（它们 fold 全部 bin）。

**根因**（`crates/reed/src/stats/sink.rs:756-767`，v27 V27-30 的提速引入，`e920663`）：

```rust
while s.saturating_add(size) <= watermark {      // ① 用旧的 s 判断"窗口已完整"
    if self.skip_empty {
        let Some(next_bin) = bins.next_bin_at_or_after(s / width) else { break };
        let earliest = next_bin.saturating_mul(width).saturating_sub(size - 1);
        if earliest > s { s = earliest.div_ceil(step) * step; }   // ② 把 s 跳到后面
    }
    let end = s + size;                           // ③ 新窗口的 end 可能 > watermark，却照样输出并 drop 了 bin
```

① 的完整性检查针对跳之前的 `s`；② 之后新窗口可能还没读完，③ 仍然输出并在循环后 `drop_bins_below` 丢掉它的 bin。
第一个位点到达时 watermark 前进，空窗口被跳过，跳到含该位点的窗口，窗口里只有 1 个位点就被输出了。
现有测试 `skip_empty_drops_window_rows_only`（`crates/reed/tests/e2e_stats.rs:262`）只断言**行坐标**，不断言计数，所以没抓到。

**修法**：跳之后重新检查完整性。

```rust
while s.saturating_add(size) <= watermark {
    if self.skip_empty {
        let Some(next_bin) = bins.next_bin_at_or_after(s / width) else { break };
        let earliest = next_bin.saturating_mul(width).saturating_sub(size - 1);
        if earliest > s {
            s = earliest.div_ceil(step) * step;
            if s.saturating_add(size) > watermark { break; }   // 跳到的窗口还没读完：留给下一次
        }
    }
    // …其余不变
}
```

`close_contig`（`:827-835`）的同构跳转不需要这个检查（contig 已结束，所有 bin 都在），但要加同一组回归测试。

**门**：① 上面的 10 位点夹具：两种模式下非空窗口的 `no_sites/count_diffs/count_comparisons` 逐列相同；
② e2e：`--skip-empty` 的输出 == 不加 `--skip-empty` 的输出删掉 `no_sites=0` 的行（任意数据、任意步长）；
③ vcftools `--window-pi` 对拍加一组 `--skip-empty`。

### N-2（P1）PI_HAT 的 `nn` 不含单态位点；plink 含

`crates/reed/src/genotype.rs:290-301` 在收集剂量时**丢掉没有任何 ALT 调用的位点**（`pca` 与 `distance` 共用）。
plink 1.9 `--genome` 的每对位点数 `nn` 包含单态位点（`plink_calc.c:4381`），所以加入单态位点会改变 plink 的 PI_HAT：

```
subvcf + 300 个全 0/0 位点       plink 1.9 首对 PI_HAT 0.6316 → 0.6652（IBS2 2215 → 2515）
复算：nn 含全部位点               与 plink 4 位小数 0 差
复算：nn 只含多态位点             最大差 0.0489
```

即使 V30-01 的系数修好，只要输入含单态位点，PI_HAT 仍对不上。修法：PI_HAT 的 `nn/ibs` 计数用**过滤前**的全部双等位位点
（系数仍只用多态位点）。其他距离类型不需要改：加入全 REF 位点时，`king` 与 `ibs` 在 reed 和 plink/plink2 上都不变（全 ALT 位点的 `ibs` 见 N-4）。

注意：plink2 `--make-rel`、`--pca` 对全 REF 位点**会变**（GRM 对角线比 0.909，特征值 6.76 → 6.15），reed 不变。
这是 reed 有意跟随 plink2 "单态列标准化为 0"的结果（`reed-distance/src/standard.rs:101-108`），但 reed 是靠**丢弃**达到的，
与 plink2 实际行为不同；plink 上游 2026-09-23 已把 `--pca` 遇到单态位点改为报错（`plink_calc.c:7256-7268`）。
建议：`distance --type grm` 与 `pca` 文档写明"单态位点被丢弃，因此与 plink2 在含单态位点的输入上不同；先过 `--mac 1`"。

### N-3（P2）`.genome` 的 DST/PPC/RATIO 恒为 0，RT/PHE 与 plink 不同

`crates/reed-cli/src/translate/distance.rs:427-432` 硬写 `OT … -9 0 0 0`。plink 对同一对给 `UN … -1 0.865430 0.8413 NA`
（`plink_calc.c:4450-4463`：DST = IBS 距离，PPC = 正态近似检验，RATIO = 两类计数之比）。下游按 PI_HAT 以外的列筛关系的脚本会静默得到 0。
修法：照 plink 计算 DST/PPC/RATIO（需要的 IBS0/IBS1/IBS2/HETHET 计数在 `pi_hat` 里已有），或删掉这三列并在文档说明。

### N-4（P2，对拍基线）本机 plink 1.9 早于上游的"全 ALT 位点"修复；reed 与上游现版一致

reed 的缺失权重 `hwe_weights`（`crates/reed-distance/src/matrices.rs:228-252`）对 p=1 的位点给 `p(1−p)(p²−p+1) = 0`。
上游 plink 1.9 在 2026-09-07 把 p=1 的权重显式置 0（`plink-ng/1.9/plink_calc.c:7778-7781` "bugfix (7 Sep 2026)"），与 reed 相同。
本机的 `~/biosoft/plink` 是 v1.9.0-b.7.11（2025-08-19），**早于**该修复。

实测：subvcf 上 reed 与本机 plink 的 IBS 最大差 5e-7（打印精度）；加 300 个全 1/1 位点后最大差 1.5e-3；全 0/0 位点不影响。
把 300 个全 ALT 位点的缺失权重设为 1 后，按 reed 的公式复算与本机 plink 的差回到 5e-7（全部 60,378 对），
所以差异完全来自旧版 plink 的这个权重，reed 的实现是对的。

处理：不改 reed。`gold/` 的 distance 对拍要固定一个**含该修复**的 plink 1.9 构建（或在夹具里避开全 ALT 位点），
并在脚本里打印 plink 版本；否则换一台机器、换一个 plink 版本就会出现假的 IBS 回归。

### N-5（P2）`--summary` 缺 PSI 与 HWE 两节，用户文档没写

对照 `bcftools stats -s -`：reed 的 SN/TSTV/SiS/AF/QUAL/IDD/ST/DP/PSC 一致（多等位夹具上还差 1 处 ST 计数和 1 处 PSC 的 nTransitions），
但**没有 PSI 与 HWE 节**。`stats.md:26` 只写"QC summary → PREFIX.stats"，`cli.md` 也未登记缺节。
修法：补这两节，或在 `stats.md` 列出已实现的节。多等位上的 ST/PSC 差异随 V30-03 的裁决一起处理。

---

## 5. 关于"滑窗"的复核（回答用户 2026-10-10 的问题）

用户问 `--indep-pairwise` 是不是 plink 的语法、reed 是不是没做滑窗、是不是所有统计量都没做滑窗。实测结论：**都做了滑窗**。

**LD 剪枝**：`reed filter --indep-pairwise` 的窗口推进照 plink2（`reed-stats/src/ld.rs:454` `evaluate_window(window_start + step)`），
输出与 plink2 `.prune.in` **逐 ID 相同**：

| 参数 | reed 保留 | plink2 保留 | 只在 reed | 只在 plink2 |
|---|---|---|---|---|
| `50 5 0.2` | 1578 | 1578 | 0 | 0 |
| `50 10 0.2` | 1579 | 1579 | 0 | 0 |
| `50 25 0.2` | 1583 | 1583 | 0 | 0 |
| `100kb 1 0.2` | 1476 | 1476 | 0 | 0 |
| `500kb 1 0.2` | 1419 | 1419 | 0 | 0 |

`100kb 5 0.2` 两者都拒绝（plink2：`window-increment must be 1 when window size is in kb`；reed：`a kb window requires step 1`）。
参数写法 `W[kb] [S] R2` 是 plink 的，按用户 2026-10-10 的裁决保留短拼写；与 plink 不同的是**一条命令直接输出剪枝后的 VCF**，没有 `.prune.in/.prune.out` 两段式。

**窗口统计**：`popg --window 10k --step 5k` 的窗口起点是 `1, 5001, 10001…`（重叠滑动）；不给 `--step` 时等于窗口大小（不重叠）。
`--step` 的默认值在 v6 曾经必须显式给，现已修复（`--help` 显示 `[default: 0]` 即"等于窗口大小"）。
与 vcftools `--window-pi 100000 --window-pi-step 50000` 对拍：21 个滑动窗口坐标全部相同，20 个数值一致；不一致的那个是 N-1。

---

## 6. 建议的处理顺序

| 顺序 | 项 | 理由 |
|---|---|---|
| 1 | N-1 `--skip-empty` | 默认数据路径上的错误数值，修法 1 行 + 回归测试 |
| 2 | V30-01 + N-2 PI_HAT | 用 §2.1 的有限样本系数 + 全位点 `nn`；plink 对拍可达逐位相同 |
| 3 | V30-02 `-I`、V30-R02 溢出 | 行为修正，各有明确参照（bcftools） |
| 4 | V30-04 f 统计尺度 | 先裁决：改为位点均值（推荐）还是保留 Σ 并改掉"对齐 ADMIXTOOLS"的说法 |
| 5 | N-3、N-5、V30-03、V30-05/06/07、R01 的 `--012`；N-4 的 plink 版本固定 | 输出列/文档/对拍基线 |
| 6 | library.md 与 `cli.md:26` 对齐；三个 crate 补 `publish = false` | 两行改动，消除对外文档自相矛盾 |
| 7 | 03 的 P0–P4 | 等产品裁决（是否做外部 Rust 库）后再排 |

## 7. 本文未覆盖

- 03 的 C/D/E 根因（能力类型、资源 owner、CI、rustdoc 门）只抽查了 R06/R09/Q31/Q34，其余沿用 01 的证据，未复跑。
- 没有跑 `cargo test --workspace`，没有跑 gold 全量对拍；N-1 的修法只做了代码推理和夹具复现，未实现。
- sgkit 覆盖表（01 §1）未逐行复核。
