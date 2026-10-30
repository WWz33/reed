---
title: index
---

## 用法

```
reed index [OPTIONS] <FILE>
```

给 VCF/BCF 建 CSI 或 TBI 索引，默认写 CSI。`filter -r/-R` 的区域查询、以及任何按区域随机访问的下游工具都依赖输入旁边的索引。

## 参数

```
<FILE>                      输入文件，必须是 BGZF 压缩的 VCF 或 BCF。stdin 被拒绝：
                            索引要记录文件的字节末偏移并回读
-c, --csi                   写 CSI 索引（默认；与 -t 互斥）
-t, --tbi                   写 TBI 索引
-o, --output FILE           索引写到别处（默认写在输入旁，加 .csi/.tbi 后缀）
-f, --force                 覆盖已存在的索引；不加且索引已存在时报错并提示
    --threads INT           工作线程数，默认全部可用核
    --verbosity INT         0 静默、1 摘要（默认）、2 阶段耗时、3 调试
```

## 行为细节

- 输入必须是文件。普通 gzip 或纯文本 VCF 报错：对它们建索引会得到全零虚拟偏移。
- BCF 不能配 TBI；BCF 只能用 CSI。
- TBI 的坐标上限是 2²⁹−1；更大的位点只能用 CSI。
- 要求排序输入：同一 contig 内位置倒退时报 `positions are not ascending (indexing needs sorted input)`。
- CSI 的 bin 结构与分块深度采用标准参数，产物可被任何读取 CSI 的程序使用。

## 例子

```bash
reed index raw.vcf.gz                  # raw.vcf.gz.csi
reed index -t raw.vcf.gz               # raw.vcf.gz.tbi
reed index -f raw.vcf.gz               # 覆盖旧索引
reed index -o /fast/disk/raw.csi raw.vcf.gz
```
