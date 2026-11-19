# 05 实施后复核：外部 review 12 条 + 第二轮修复规范

> 2026-11-13。对象：`04-implementation.md` 记录的 v30 第一批实施（P0–P4、R01–R14），
> 与一份外部 review（12 条，作者自行读源码复核，结论"10 条属实、1 条行号错、1 条机制要修正、1 条未核"）。
> 方法：从当前工作树重编译 `reed 0.5.0 (v0.3.0-180-gf190708-dirty)`（`/tmp/reed-target/release/reed`，2026-11-13 02:28），
> 逐条重跑、读源码、与 bcftools / vcf2phylip / vcftools / plink 1.9 输出比对。产物在 `/tmp/v30chk/`。
> **本文只写规范，不改代码。** 实施者按 §4 执行，结果写进 `STATUS.md`。

---

## 1. 结论

1. **外部 review 的判断基本成立**：我复核的 10 条（1、2、4、5、6、7、8、9、11、12）结论都对；第 3 条与它的"旁证"见 §2.3；第 10 条我也未核。
2. **review 漏了一件更重的事**：`04-recheck.md` 里的 **N-1（`--skip-empty` 截断窗口，数值错误）在当前二进制上仍然复现**，
   且 `04-implementation.md` 没有认领它。它应排在 CI 之前。
3. **v30 第一批引入了两处回归**（review 第 6 条只看到其中一半）：
   - `--fill-tags` 在样本子集时把 MAF 写在 AC/AN **前面**（旧版与 bcftools 都是 `AC;AN;MAF`）；
   - R01 把 `--phylip/--fasta/--nexus` 的 `--min-samples-locus` 判定从 vcf2phylip 口径改走，`./1` 半调用位点现在被保留，vcf2phylip 丢弃。
4. **四条审计头号问题（V30-01/02/03 与 N-1）本轮都没碰**——这不是"改坏了"，是"没排进来"；但它们都是现行数值/行为错误，
   必须在第二批里先做。

---

## 2. 外部 review 逐条复核

| # | review 结论 | 我的复核 | 证据（当前二进制/源码） |
|---|---|---|---|
| 1 PI_HAT 系数 | 属实 | **属实**；60,346/60,378 对仍为 `0/1/0/0.5` | `reed distance subvcf.vcf.gz --type pihat` |
| 2 `filter -I` 删 AC/AN | 属实 | **属实**；`-I`、`-s S -I`、无参数三种都输出 0 条带 AC 的记录 | `withAC.vcf` 三组参数 |
| 3 CI 缺 `cargo modules` | 属实 | **属实**，且范围更准确的说法见 §2.3 | `ci.yml` 的 `taiki-e/install-action` 只装 `cargo-public-api,cargo-hack,cargo-semver-checks`；`ci/check-docs.sh:21` 调 `cargo modules orphans` |
| 4 check-api 缺 nightly | 结论对，机制修正 | **review 的修正正确** | `cargo-public-api-0.52.0/src/main.rs:590-603`：stable 时打 Warning 并切 `nightly`，不安装；`public-api-0.52.2/src/lib.rs:70` 下限 `nightly-2025-11-22` |
| 5 多等位三路径 | 属实 | **属实**：`stats/run.rs:303-309` `CollapseAlt`；`:474-477` popg 跳过；`.afreq` 对 `A G,T` 仍写单值 `0.5`（plink2 `0.25,0.25`） | `multi.vcf` |
| 6 fill 顺序 | 属实 | **属实，且是本轮引入的回归**，见 §3.2 | `plan/mod.rs:308-314` `enable` 是 `push` |
| 7 msa 测试缺失 | 属实 | **属实**：`convert/job/msa.rs:240-251` 的 0 边界与 start>end 两条 usage 分支在 `crates/*/tests` 中零命中；行为本身正确（实跑 `chr1:0-10`、`chr1:20-10` 均 exit 2） | `msa.vcf` 五组 `--region` |
| 8 LD idle 不能取消 | 属实 | **属实**，机制见 §3.4 | `stats/driver.rs:376-389` |
| 9 AC 双解析 | 属实，行号错 | **review 的行号修正正确**：双循环在 `reed-core/src/count/info_counts.rs:78-103`；`:112` 起才是 `#[cfg(test)]`；第二循环的 `return Absent`（`:95-97`）不可达 | 源码 |
| 10 header 双扫 | 未核 | 未核 | — |
| 11 Framing 死 API | 属实 | **属实**：`reed/src/input.rs` 的 `InputMeta::framing()` 只在 `input.rs:651,660` 的测试里读；`reed-engine/src/pipeline.rs:538,803` 读的是 `ChunkJob::framing()`，不是同一个方法 | grep |
| 12 描述过期 | 属实 | **属实**：`crates/reed/Cargo.toml:3` 列 6 个 `run_*`，`lib.rs:112,121` 还有 `run_pca`/`run_dapc`（及 `_opened`） | 源码 |

### 2.1 第 4 条的准确表述（采纳 review 的修正，补两点）

- `ci.yml` 的 `gates` job 只装 `stable`（`components: clippy, rustfmt`）。`cargo public-api` 检测到 stable 后切到 `nightly`
  这个**名字**；runner 上没有就失败，有就看它是否 ≥ `nightly-2025-11-22`。所以门的结果取决于 runner 镜像，不取决于代码。
- 快照 `ci/api/reed.txt` 是用本机 `nightly`（`rustc 1.101.0-nightly 2026-10-24`）生成的；rustdoc JSON 随 nightly 变，
  **即使 runner 有 nightly，版本不同也可能出现非代码引起的快照漂移**。固定版本才是可复现的门。

### 2.2 第 9 条的处理

不可达分支不是 bug，但它让读者以为第二循环还会遇到非整数。改成一次遍历（见 §4 第 7 项的代码），
既去掉不可达分支，也保留"先计数后分配"的 V26-02 目的：

### 2.3 第 3 条与"本地绿、CI 红"

本机装了 `cargo-modules 0.28.0`、`cargo-public-api 0.52.0`、`cargo-hack 0.6.46`、`cargo-semver-checks 0.51.0`，
所以 `04-implementation.md:329` 的 `check-docs.sh` exit 0 是真的；CI 上 `gates` job 没装 `cargo-modules`，
`check-docs.sh` 在 `set -euo pipefail` 下第 3 段必然失败。同理 `check-api.sh` 依赖 runner 的 nightly（§2.1）。
**结论：这两道门目前是"本地可过、CI 不可复现"**。

---

## 3. 复核中新发现或需要补充的事实

### 3.1 N-1 仍然存在（P0，数值错误）

```
gap.vcf（10 位点：1k–5k 与 205k–209k），--window 10k
不加 --skip-empty   200001-210000   avg_pi 0.00028571   no_sites 5
加   --skip-empty   200001-210000   avg_pi 0.00005714   no_sites 1     ← 当前二进制仍然如此
```

根因与修法见 `04-recheck.md` §4 N-1（`crates/reed/src/stats/sink.rs:756-767` 跳转后未重查完整性）。
`04-implementation.md` 的 R01–R14 与 P0–P4 都不覆盖它。

### 3.2 `--fill-tags` 顺序：本轮回归，范围比 review 说的窄也更具体

| 场景 | 旧二进制（`02109e9`，R13 前） | 当前二进制 | bcftools |
|---|---|---|---|
| 子集 `-s s1,s2` + `--fill-tags MAF` | `AC=3;AN=4;MAF=0.25` | **`MAF=0.25;AC=3;AN=4`** | `AC=3;AN=4;MAF=0.25` |
| 不子集，`--fill-tags AC,AN,NS,MAF,HWE`（gold 规则） | `AC;AN;NS;MAF;HWE` | `AC;AN;NS;MAF;HWE` | `NS;AN;MAF;AC;HWE` |

- **回归**：`filter/mod.rs:555-556` 在 compile 之后 `enable("AC")`/`enable("AN")`，`RecomputeFlags::enable` 是 `push`
  （`reed-filter/src/plan/mod.rs:308-314`），所以 AC/AN 落在用户标签之后。R13 之前是固定字段顺序 AC,AN,NS,MAF,HWE。
- **与 bcftools 的顺序本来就不同**：bcftools 按 `NS→AN→AF/MAF→AC→…→HWE` 写（`plugins/fill-tags.c:815-880`），
  reed 一直按 AC,AN,NS,MAF,HWE。gold 对拍没抓到，因为 `gold/scripts/check_parity.py:34` 只比 `CHROM/POS/REF/ALT`，不比 INFO。
  这是规范层面的选择，本文不要求改成 bcftools 顺序；**只要求恢复 R13 之前的固定顺序**并锁住。

### 3.3 R01 改变了 phylo 矩阵的 `--min-samples-locus`，偏离 vcf2phylip

`04-implementation.md` 认为改成 `has_called_allele` 后"半缺失与 stats 一致"。代价是 phylo 输出不再是 vcf2phylip：

```
half.vcf：位点 10 = ./1 ×4 + ./. ×4；位点 20 = 0/. ×4 + ./. ×4；-m / --min-samples-locus 4
vcf2phylip 2.9   8 1   （./1 以 "." 开头，按 num_genotypes 计为缺失 → 位点 10 只有 0 个样本 → 丢弃）
reed（当前）     8 2   （两位点都保留）
```

`gold/popgen/run_convert_parity.sh` 仍 PASS，因为 `gold/fixtures/phylo.vcf` 只有 1 个半调用、且不在门限边界上。
convert 的用户文档写"vcf2phylip semantics"（`convert/mod.rs:2`、`job/mod.rs:6`），现在这个承诺对半调用不成立。
`--012` 不受影响：两位点仍全部写 −1（与改动前相同，与 vcftools 的 `1` 不同，parity 脚本已登记）。

**需要裁决**（见 §5）：phylo 的位点门跟 vcf2phylip，还是跟 stats。

### 3.4 第 8 条的机制

`stats/driver.rs:376-389`：

```rust
abort.trigger_handle().cloned()
     .unwrap_or_else(|| reed_engine::abort::abort_channel_for(pipeline_threads).0),
abort.clone(),
```

`RunContext::idle` 的信号没有 trigger handle（`reed-engine/src/abort.rs:92-110`，`handle: None`）。
此时 LD executor 拿到一个**新建通道的 handle**，但 `signal` 仍是 idle 信号的克隆；`LdExecutor::drop`（`sink.rs:575-588`）
调用的 `self.abort.trigger()` 置的是新通道的标志，而 LD 工作循环轮询的是 `signal.is_aborted()`（`sink.rs:540,612`）——
**两者不是同一个标志**，所以 idle 上下文里出错退出时，LD 块循环收不到取消，只能等 rayon 池 drop 时 join 跑完。
CLI 不受影响：七个前端都用 `abort_for(&budget)`（`reed-cli/src/translate/*.rs`）。受影响的是库调用者与 `prelude` doctest 的写法。

### 3.5 AC/AN 读取器（R02）与 bcftools 的剩余差异

R02 修掉了 debug panic 与多 ALT 求和回绕（`ovf.vcf` 现在与 `bcftools stats -s -` 逐列相同）。剩下两处：

| 输入（单 ALT，GT `0/1 1/1`） | `bcftools stats -s -` AF | reed `--summary` AF |
|---|---|---|
| `AC=18446744073709551615;AN=100` | 0.000000（htslib 把超 int32 的 INFO 整数置缺失，`vcf.c:4037-4049`） | **0.740000**（u64 解析成功，AC>AN → GT 回退） |
| `AC=.;AN=100` | 0.000000 | **0.740000** |
| `AC=3000000000;AN=4000000000` | 0.740000 | 0.740000 |
| `AC=3;AN=4` | 0.740000 | 0.740000 |

bcftools 在 AC 为缺失值时仍走 INFO 分支（`vcfutils.c:63-89` 把缺失值当作计数相加，得到 AF 桶 0），这是 htslib 的怪癖，
**不建议跟随**。本文只要求：`.stats` 的这两种输入写明"reed 回退到 GT"，并在 `stats.md` 的 `--summary` 段登记。

---

## 4. 第二批修复规范（按顺序）

| 序 | 项 | 改什么 | 门 |
|---|---|---|---|
| 1 | **N-1 `--skip-empty`** | `sink.rs` 跳转后重查完整性（代码见下） | §4.1 |
| 2 | **CI 工具链** | `gates` 装 `cargo-modules`；API 快照固定 nightly | §4.2 |
| 3 | **V30-01 PI_HAT** | `04-recheck.md` §2.1 的有限样本系数 + N-2 的全位点 `nn`，与 plink 对拍一起做 | §4.3 |
| 4 | **V30-02 `-I`** | `-I` 时 AC/AN 不进默认删除清单 | §4.4 |
| 5 | **fill 顺序回归** | 按表顺序输出，不按启用顺序 | §4.5 |
| 6 | **R01 phylo 门** | 按 §5 的裁决 | §4.6 |
| 7 | V30-03、msa 测试、LD idle 取消、AC 单循环、Framing、描述 | 各自一行到十几行 | §4.7 |

### 4.1 N-1

```rust
// crates/reed/src/stats/sink.rs, emit_complete_windows
while s.saturating_add(size) <= watermark {
    if self.skip_empty {
        let Some(next_bin) = bins.next_bin_at_or_after(s / width) else { break };
        let earliest = next_bin.saturating_mul(width).saturating_sub(size - 1);
        if earliest > s {
            s = earliest.div_ceil(step) * step;
            // The window we jumped to may still be filling: the loop guard
            // above checked the old `s`. Leave it for the next call.
            if s.saturating_add(size) > watermark {
                break;
            }
        }
    }
    // …unchanged…
}
```

门：
1. 新测试（`crates/reed/tests/e2e_stats.rs`）：10 位点两簇夹具，`--skip-empty` 与不加时，非空窗口的
   `no_sites/count_diffs/count_comparisons/count_missing` 逐列相等；1/4/56 线程相同。
2. 性质测试：对任意 `window/step`，`--skip-empty` 的输出 == 不加时的输出删掉 `no_sites == 0` 的行。
3. 现有 `skip_empty_drops_window_rows_only`（`e2e_stats.rs:262`）补断言 `no_sites`，不只断言坐标。

### 4.2 CI 工具链

```yaml
# .github/workflows/ci.yml, gates job
      - uses: actions-rust-lang/setup-rust-toolchain@v1
        with:
          toolchain: stable,nightly-2026-03-20     # 后者供 cargo public-api；与 external-types job 同一个
          components: clippy, rustfmt
      - name: Install cargo-public-api, cargo-hack, cargo-semver-checks, cargo-modules
        uses: taiki-e/install-action@v2
        with:
          tool: cargo-public-api@0.52.0,cargo-hack@0.6.46,cargo-semver-checks@0.51.0,cargo-modules@0.28.0
```

```bash
# ci/check-api.sh：显式指定工具链，不依赖 "nightly" 这个浮动名字
cargo +nightly-2026-03-20 public-api -p reed -s >"$tmp"
```

- `nightly-2026-03-20` 满足 `MINIMUM_NIGHTLY_RUST_VERSION = nightly-2025-11-22`，且与 `check-external-types.sh` 已固定的版本相同，
  一个仓库只固定一个 nightly。
- **用同一 nightly 重新生成 `ci/api/reed.txt`**（当前快照来自 `nightly 2026-10-24`），否则第一次跑就会出现工具链漂移。
- `check-docs.sh` 第 3 段不要求 CI 外的人装 `cargo-modules`：保留严格模式，但在开头检查工具并给出安装命令
  （`command -v cargo-modules || { echo "cargo install cargo-modules@0.28.0"; exit 1; }`），避免本地报一个难懂的错误。
- 门：在一个**干净的** runner（或 `docker run rust:1.99`）上从零跑 `gates`，全绿；`04-implementation.md` 的"本地 exit 0"改为 CI 链接。

### 4.3 V30-01 PI_HAT

按 `04-recheck.md` §2.1（`expected_ibs_coefficients` 的完整有限样本版本）与 N-2（`nn`/IBS 计数用过滤前的全部双等位位点）
**一起**做。只换系数而不改 `nn`，含单态位点的输入仍与 plink 差到 0.049。

门：`gold/popgen/run_genome_parity.sh`（新增）——plink 1.9 `--genome full` vs `reed distance --type pihat`，
Z0/Z1/Z2/PI_HAT 四位小数逐对相同；夹具含缺失与单态位点；脚本打印 plink 版本。
本轮 numpy 复算已证明 subvcf 与 subvcf+300 单态位点两份输入都能做到 0 差。

### 4.4 V30-02 `-I`

```rust
// crates/reed-filter/src/plan/compile/finish.rs（或 filter/mod.rs 合成 stale 清单处）
// bcftools `view -I`: keep the input's AC/AN. Only an explicit `-x AC`/`-x AN`
// still removes them.
let mut stale_tags = spec.recompute.stale_bytes();
if spec.recompute.no_update && !spec.recompute.drop_stale_explicit {
    stale_tags.retain(|t| t.as_slice() != b"AC" && t.as_slice() != b"AN");
}
```

门：`-s S -I` 输出带原值 AC/AN（`withAC.vcf` 上首条 `AC=331`）；`-I` 与 `-I -x ''` 对 AC/AN 等价；
`docsite/filter.md:155`、`cli.md:140`、`--help` 三处措辞一致（"keep the input's INFO/AC and INFO/AN"）。

### 4.5 `--fill-tags` 顺序

```rust
// crates/reed-filter/src/plan/mod.rs —— enable 后按注册表顺序重排，输出与启用顺序无关
pub fn enable(&mut self, name: &str) {
    if let Some(tag) = reed_core::count::RecomputeTag::find(name)
        && !self.wants(tag.name)
    {
        self.tags.push(tag);
        let order = |t: &&'static reed_core::count::RecomputeTag| {
            reed_core::count::RECOMPUTE_TAGS.iter().position(|r| r.name == t.name)
        };
        self.tags.sort_by_key(order);
    }
}
```

（compile 时从用户列表构建 `tags` 的地方同样按表顺序排序；或直接让 `RecomputeFlags` 存一个 `u8` 位集、`fill_info_set` 遍历 `RECOMPUTE_TAGS`。）

门：
1. 单测：`enable("AC")` 在 `MAF` 之后调用，`fill_info_set` 输出仍是 `AC;AN;MAF`。
2. e2e：`-s s1,s2 --fill-tags MAF` 输出 `AC=3;AN=4;MAF=0.25`（R13 前的字节）。
3. `gold/scripts/check_parity.py` 增加 INFO 比较（按标签集合比，不比顺序，因为 bcftools 顺序本就不同），
   并另加一个 reed 自身的"INFO 字节快照"门，让顺序变化必须显式更新快照。

### 4.6 R01 phylo 门

按 §5 的裁决二选一：

- **(a) 恢复 vcf2phylip 口径（推荐）**：phylo 目标的 `present` 用 `!raw.starts_with(b".")`（只用于 `--min-samples-locus` 计数），
  `--012/--hapmap/--msa` 继续用 `has_called_allele`。在 `convert.md` 写明"phylo 矩阵的样本计数沿用 vcf2phylip：以 `.` 开头的 GT 记为缺失"。
- **(b) 保留当前统一口径**：`convert.md`、`convert/mod.rs:2`、`job/mod.rs:6` 删掉"vcf2phylip semantics"对位点门的承诺，
  写明与 vcf2phylip 在 `./x` 上的差异。

门（两种都要）：`gold/fixtures/` 加 `half.vcf`（上面 8 样本 2 位点），`run_convert_parity.sh` 新增一段：
(a) 与 vcf2phylip `-m 4` 逐字节相同；(b) 断言 reed `8 2`、vcf2phylip `8 1`，并打印差异原因。

### 4.7 其余

| 项 | 规范 | 门 |
|---|---|---|
| V30-03 多等位 | 按 V27-39 统一为"拒绝/跳过并计数"：`stats` 的 `.afreq/.hardy` 对 `n_alt > 1` 跳过并在 `.log` 记 `multiallelic skipped`；`--summary` 保持与 bcftools 一致（它本就对齐） | `multi.vcf`：`.afreq` 只剩双等位行；`.log` 计数 1 |
| msa 测试 | `e2e_convert.rs` 加两个用例：`--region chr1:0-10`、`chr1:20-10` → exit 2 且错误文本固定 | 两个测试 |
| LD idle 取消 | idle 信号时 LD executor 自建一对 `(handle, signal)` 并用**同一对**，或 `RunContext::idle` 直接用 `abort_channel_for(budget.threads)` 建真通道（推荐，删掉"无 handle 信号"这个状态） | 库级测试：idle 上下文 + 注入错误后 LD 工作在 1 s 内返回 |
| AC 单循环 | `info_counts.rs:78-103` 合成一次遍历：逐段 `parse`（失败即 `Absent`）、计数超 `n_alt` 即 `AcLength`、`checked_add` 溢出即 `AcOverflow`、`per_allele` 用 `SmallVec<[u64; 4]>` 避免大 AC 字段的分配 | 现有 5 个 core 单测 + 新增"超长 AC 不分配"测试 |
| Framing | `InputMeta::framing()` 删除，或给它一个生产调用者（例如 `open_regions` 的 BGZF 检查改读它） | `rg "\.framing\(\)" crates/reed/src` 非测试命中 ≥1 或方法删除 |
| 描述 | `crates/reed/Cargo.toml:3` 改为不列函数名（"Facade library for reed: request types, run_* drivers and reports"），由 `check-docs.sh` 的入口表门负责列举 | `check-package.sh` |
| AC/AN 文档 | `stats.md` 的 `--summary` 段：INFO/AC 为缺失值或超范围时 reed 回退到 GT，bcftools 记入 AF=0 桶 | 文档 |

---

## 5. 需要你裁决

| 项 | 推荐 | 备选 |
|---|---|---|
| phylo 的 `--min-samples-locus` 口径（§3.3） | **跟 vcf2phylip**（convert 的定位就是 vcf2phylip 的替代品，用户拿两者结果对比时不应出现位点数差异） | 跟 `stats`，删掉"vcf2phylip semantics"承诺 |
| `--fill-tags` 输出顺序 | **恢复 R13 前的 AC,AN,NS,MAF,HWE**（gold 输出不变） | 改成 bcftools 顺序 NS,AN,MAF,AC,HWE（gold 输出字节会变） |

---

## 6. 本文未覆盖

- review 第 10 条（header 双扫）未核。
- 没有跑 `cargo test --workspace`；只跑了 `gold/popgen/run_convert_parity.sh`（PASS）与上述手工夹具。
- `04-implementation.md` 的 P3/P4 机器门除 §2 提到的两道外，未在干净环境复跑。
