# reed

[English](README.md) | [简体中文](README.zh-CN.md)

用于从 VCF 和 BCF 文件计算窗口化群体遗传学统计量的 Rust 库。

## 支持的统计量

| 类别 | 统计量 | 支持的方法与估计器 |
|---|---|---|
| 多样性 | $\pi$, Watterson's $\theta$, Tajima's $D$ | 逐样本对可比位点（PerPairComparable）、可及碱基对（AccessibleBp） |
| 群体分化 | $F_{ST}$, $d_{xy}$, PBS | Hudson (2013), Weir & Cockerham (1984), Nei (1986); Yi (2010), PBSn1 |
| 等位基因频率 | 等位基因计数, SFS | 折叠与展开位点频率谱 |
| 连锁不平衡与 HWE | Rogers-Huff $r^2$, HWE 精确检验, $F_{IS}$ | Wigginton (2005) 精确检验 |
| 渗入与 $f$ 统计量 | $f_2$, $f_3$, $f_4$, Patterson's $D$ | 样本量校正估计器 |
