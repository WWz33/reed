# reed

Rust 群体遗传学统计库。设计原则:每个 `(统计量, 估计器)` 是**封死的完整定义**——精确公式、分母约定、缺失数据处理、文献引用、金标准对比测试,四者不可拆分。设计决策全文见 [`DESIGN.md`](DESIGN.md)。

## 架构(三层)

```
L1 传输   noodles(src/io)          惰性记录流,noodles 单后端
L2 物化   reed-core::WindowBuffer  窗口缓冲:定长窗口、contig 长度为一等公民、
                                    单态记录可入缓冲、逐基因型相位位
L3 统计   reed-core::estimators    纯函数 over 不可变窗口;结果永远是精确整数计数 + 派生的 f64
```

## 已实现统计族

| 族 | 内容 |
|---|---|
| 频率 | allele counts、SFS(折叠/展开) |
| 多样性 | π(PerPairComparable / AccessibleBp)、Watterson θ、Tajima's D |
| 分化 | dxy(两估计器)、Fst(Hudson / Weir-Cockerham 84 / Nei 86)、PBS(Yi2010 / PBSn1) |
| HWE-近交 | HWE 精确检验(Wigginton 2005)、Hobs/Hexp/Fis |
| LD | Rogers–Huff r²(逐对 + 窗口汇总) |
| f-统计 | f2 / f3 / f4 / D-stat |

未实现(需新机械,见 DESIGN §12):PCA/GRM/relatedness(矩阵积累层)、EHH/iHS/Garud H(相位扫描层)、πN/πS(注释层)。

## 快速开始

```bash
export PATH=$HOME/.cargo/bin:$PATH   # rustup stable ≥ 1.91
cd reed-rs
cargo test --workspace               # 27 tests
# 扫描真实 VCF/BCF(每 window_bp 一个窗口,打印 π/Fst)
cargo run --release -p reed-io --example scan -- <file.vcf|file.vcf.gz|file.bcf> [window_bp]
```

```rust
use reed_core::estimators::fst::{fst, FstEstimator};
use reed_core::{Ascertainment, PopulationView, SiteCompleteness};

// 一次声明数据事实;估计器按清单推导可用性
let manifest = SiteCompleteness::VariantOnly(Ascertainment::Complete);
// for window in reader { ... }
// let r = fst(&window, &pops, 0, 1, FstEstimator::WeirCockerham84)?;
// r.value()          -> Option<f64>
// r.wc               -> 精确的 a/b/c 组件(可审计)
```

## 输入约束(v0.1)

- 仅 VCF / VCF.gz / BCF;二倍体;双等位基因(多等位显式报错,不静默丢弃)。
- `contig_len` 必需(缺失即报错;Beagle 输出会丢 `##contig`)。
- 窗口网格 `[start, start+size)`、步长 `step`,不跨 contig,尾部按 contig 长度截断。

## 验证

金标准对拍 scikit-allel 1.3.13(venv 见 `/tmp/allelvenv` 可重建)、HWE 用 Levene-Haldane 独立枚举交叉验证;全部常量固化为测试断言,容差 1e-12。许可证政策:仅 MIT/Apache/BSD 源可翻译(vcftools LGPL、plink2/GCTA/flashpca/SNPGenie GPL 只读)。
