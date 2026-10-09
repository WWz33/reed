#!/usr/bin/env python3
"""gold/scripts/compare_tools.py — 三方 (bcftools/vcftools/plink2) 过滤结果一致性对拍。

比较对象是三个工具的 **common 层** 输出（同一套规范、同一份样本名单）：
  out/bcftools.common.vcf.gz
  out/vcftools.vcf.gz
  out/plink2.vcf.gz

输出：
  out/concordance.tsv   机器可读的成对/三方位点集合统计
  out/concordance.md    人读报告
  控制台摘要

位点身份用 (CHROM, POS, REF, ALT) 规范化后比较（本数据全为已归一化双等位 SNP，
故该键稳定）。基因型一致性用 `bcftools gtcheck` 在同名样本间计算，可选。
"""
from __future__ import annotations
import argparse
import json
import os
import subprocess
import sys
from itertools import combinations

TOOLS = ["bcftools", "vcftools", "plink2"]
OUTNAME = {
    "bcftools": "bcftools.common.vcf.gz",
    "vcftools": "vcftools.vcf.gz",
    "plink2": "plink2.vcf.gz",
}


def which_bcftools() -> str:
    return os.environ.get("BCFTOOLS", "bcftools")


def site_keys(vcf: str) -> set[str]:
    """流式读取 (CHROM,POS,REF,ALT) 集合。"""
    proc = subprocess.run(
        [which_bcftools(), "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", vcf],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools query 失败 on {vcf}: {proc.stderr[:300]}")
    keys = set()
    for line in proc.stdout.splitlines():
        if line:
            c, p, r, a = line.split("\t")
            keys.add(f"{c}:{p}:{r}:{a}")
    return keys


def sample_list(vcf: str) -> list[str]:
    p = subprocess.run([which_bcftools(), "query", "-l", vcf], capture_output=True, text=True)
    return [s for s in p.stdout.splitlines() if s] if p.returncode == 0 else []


def gtcheck_discordance(ref: str, qry: str) -> tuple[float | None, float | None]:
    """同名样本间的基因型不一致 (bcftools gtcheck -u GT,GT)。

    返回 (mean, max)；0 表示完全一致。失败返回 (None, None)。
    """
    import tempfile
    samps = sample_list(ref)
    if not samps:
        return None, None
    with tempfile.NamedTemporaryFile("w", suffix=".pairs", delete=False) as fh:
        for s in samps:
            fh.write(f"{s}\t{s}\n")
        pairs = fh.name
    try:
        proc = subprocess.run(
            [which_bcftools(), "gtcheck", "-u", "GT,GT", "-g", ref, "-P", pairs, qry],
            capture_output=True, text=True,
        )
    finally:
        os.unlink(pairs)
    if proc.returncode != 0:
        return None, None
    vals = []
    for line in proc.stdout.splitlines():
        if line.startswith("DCv2"):
            f = line.split("\t")
            try:
                vals.append(float(f[3]))
            except (IndexError, ValueError):
                pass
    if not vals:
        return None, None
    return sum(vals) / len(vals), max(vals)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--no-genotype", action="store_true", help="跳过基因型一致性检查")
    args = ap.parse_args()
    od = args.outdir

    sets: dict[str, set[str]] = {}
    missing = []
    for t in TOOLS:
        f = os.path.join(od, OUTNAME[t])
        if not os.path.exists(f):
            missing.append(f)
            continue
        sets[t] = site_keys(f)
    if missing:
        print("缺少输出文件，跳过对应工具:\n  " + "\n  ".join(missing), file=sys.stderr)
    present = [t for t in TOOLS if t in sets]
    if len(present) < 2:
        print("少于两个工具产出，无法对拍。", file=sys.stderr)
        return 2

    lines = []          # tsv rows
    md = []             # markdown

    # 输入路径：优先环境变量，否则从任一 stats.json 取
    inp = os.environ.get("INPUT", "")
    if not inp:
        for t in present:
            p = os.path.join(od, f"{t}.stats.json")
            if os.path.exists(p):
                try:
                    inp = json.load(open(p)).get("input", "")
                except (ValueError, OSError):
                    pass
                if inp:
                    break

    md.append("# 三方过滤一致性报告\n")
    md.append(f"- 输入: `{inp}`")
    md.append(f"- 比较工具: {', '.join(present)}\n")
    md.append("## 位点数\n")
    md.append("| 工具 | 位点数 |")
    md.append("|---|---|")
    for t in present:
        md.append(f"| {t} | {len(sets[t]):,} |")
    md.append("")

    lines.append("kind\ta\tb\tn_a\tn_b\tn_ab\tn_a_only\tn_b_only\tjaccard")

    # 成对
    md.append("## 成对一致性 (位点集合)\n")
    md.append("| A | B | A∩B | A\\B | B\\A | Jaccard |")
    md.append("|---|---|---|---|---|---|")
    pair_stats = {}
    for a, b in combinations(present, 2):
        ia = sets[a] & sets[b]
        na = sets[a] - sets[b]
        nb = sets[b] - sets[a]
        jac = len(ia) / len(sets[a] | sets[b]) if (sets[a] | sets[b]) else 1.0
        pair_stats[(a, b)] = (len(ia), len(na), len(nb), jac)
        md.append(f"| {a} | {b} | {len(ia):,} | {len(na):,} | {len(nb):,} | {jac:.6f} |")
        lines.append(f"pair\t{a}\t{b}\t{len(sets[a])}\t{len(sets[b])}\t{len(ia)}\t{len(na)}\t{len(nb)}\t{jac:.6f}")
    md.append("")

    # 三方
    if len(present) == 3:
        inter = sets[present[0]] & sets[present[1]] & sets[present[2]]
        union = sets[present[0]] | sets[present[1]] | sets[present[2]]
        md.append("## 三方\n")
        md.append(f"- 三方交集 (全部保留): **{len(inter):,}**")
        md.append(f"- 三方并集: {len(union):,}")
        md.append(f"- 交集/并集 (Jaccard): {len(inter)/len(union):.6f}")
        for t in present:
            others = set().union(*[sets[x] for x in present if x != t])
            md.append(f"- {t} 独有 (其余两工具都没有): {len(sets[t] - others):,}")
        lines.append(f"triple\t{present[0]}\t{present[1]}\t{present[2]}\t{len(sets[present[0]])}\t{len(sets[present[1]])}\t{len(sets[present[2]])}\t{len(inter)}\t{len(union)}")
        md.append("")

    # 基因型一致性
    if not args.no_genotype and len(present) >= 2:
        md.append("## 基因型一致性 (同名样本, bcftools gtcheck -u GT,GT)\n")
        md.append("| A | B | 平均不一致 | 最大不一致 |")
        md.append("|---|---|---|---|")
        for a, b in combinations(present, 2):
            mean, mx = gtcheck_discordance(os.path.join(od, OUTNAME[a]), os.path.join(od, OUTNAME[b]))
            if mean is None:
                md.append(f"| {a} | {b} | n/a | n/a |")
            else:
                md.append(f"| {a} | {b} | {mean:.4g} | {mx:.4g} |")
        md.append("")

    with open(os.path.join(od, "concordance.tsv"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(od, "concordance.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")

    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
