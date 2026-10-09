#!/usr/bin/env python3
"""gold/scripts/profile_input.py — 过滤前画像：密度 + 各字段分布 + 过滤漏斗。

在设定过滤参数之前先看清数据，避免参数落在「全保留/全删除」的无效区间。
输出 out/profile.{md,json}。

用法:
    python3 profile_input.py --input work/subset.Chr01-03.vcf.gz --outdir out
"""
from __future__ import annotations
import argparse
import json
import os
import subprocess
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
GOLD_DIR = os.path.dirname(HERE)


def bcftools() -> str:
    return os.environ.get("BCFTOOLS", "bcftools")


def sh(args: list[str]) -> str:
    p = subprocess.run(args, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{' '.join(args)}\n{p.stderr[:500]}")
    return p.stdout


def pct(vals: list[float], qs=(0, 1, 5, 25, 50, 75, 95, 99, 100)) -> dict:
    if not vals:
        return {}
    vals = sorted(vals)
    n = len(vals)

    def at(q):
        i = min(n - 1, max(0, int(round(q / 100 * (n - 1)))))
        return vals[i]
    return {f"p{q}": at(q) for q in qs}


def hist(vals, edges):
    """返回 [(lo,hi,count), ...]。"""
    out = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        c = sum(1 for v in vals if lo <= v < hi)
        out.append((lo, hi, c))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", default=os.path.join(GOLD_DIR, "out"))
    ap.add_argument("--gt-sample-records", type=int, default=2000,
                    help="用于统计逐基因型 DP/GQ 分布的记录数（抽样）")
    ap.add_argument("--min-qual", type=float, default=30)
    ap.add_argument("--gt-dp-min", type=int, default=5)
    ap.add_argument("--gt-dp-max", type=int, default=100)
    ap.add_argument("--gt-gq-min", type=int, default=20)
    ap.add_argument("--site-miss-max", type=float, default=0.10)
    ap.add_argument("--maf-min", type=float, default=0.05)
    ap.add_argument("--hwe-min", type=float, default=1e-6)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ---- 样本数 / contig 长度 ----
    header = sh([bcftools(), "view", "-h", args.input])
    nsamples = 0
    contigs = {}
    for line in header.splitlines():
        if line.startswith("#CHROM"):
            nsamples = len(line.rstrip("\n").split("\t")) - 9
        elif line.startswith("##contig="):
            try:
                seg = line.split("ID=", 1)[1]
                cid = seg.split(",")[0].split(">")[0]
                length = int(seg.split("length=")[1].split(",")[0].split(">")[0])
                contigs[cid] = length
            except (IndexError, ValueError):
                pass
    if nsamples == 0:
        print("无法确定样本数", file=sys.stderr)
        return 2

    # ---- 每条染色体位点数 ----
    chrom_counts = Counter(
        c for c in sh([bcftools(), "query", "-f", "%CHROM\n", args.input]).splitlines() if c
    )
    total_sites = sum(chrom_counts.values())

    # ---- 逐基因型 DP/GQ 抽样分布 ----
    def gt_sample(tag):
        out = sh([bcftools(), "query", "-f", f"[%{tag}\t]\n", args.input])
        vals = []
        for i, line in enumerate(out.splitlines()):
            if i >= args.gt_sample_records:
                break
            for tok in line.split("\t"):
                if tok and tok != ".":
                    try:
                        vals.append(int(tok))
                    except ValueError:
                        pass
        return vals

    dp_vals = gt_sample("DP")
    gq_vals = gt_sample("GQ")

    # ---- 主 pass: type/qual/pass -> strip stale INFO -> mask GT -> fill-tags ----
    stale = ("INFO/AC,INFO/AN,INFO/AF,INFO/MLEAC,INFO/MLEAF,"
             "INFO/ExcessHet,INFO/InbreedingCoeff,INFO/DP")
    mask = (f"FMT/DP<{args.gt_dp_min} | FMT/DP>{args.gt_dp_max} | FMT/GQ<{args.gt_gq_min}")
    p1 = subprocess.Popen(
        [bcftools(), "view", "-m2", "-M2", "-v", "snps", "-f", "PASS",
         "-i", f"QUAL>={args.min_qual}", args.input, "-Ou"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen([bcftools(), "annotate", "-x", stale, "-Ou"],
                          stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p3 = subprocess.Popen([bcftools(), "filter", "-e", mask, "-S", ".", "-Ou"],
                          stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2.stdout.close()
    p4 = subprocess.Popen([bcftools(), "+fill-tags", "-Ou", "--", "-t", "AC,AN,NS,MAF,HWE"],
                          stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p3.stdout.close()
    p5 = subprocess.Popen(
        [bcftools(), "query", "-f", "%QUAL\t%INFO/AN\t%INFO/MAF\t%INFO/HWE\n"],
        stdin=p4.stdout, stdout=subprocess.PIPE, text=True)
    p4.stdout.close()

    quals, site_miss, mafs, hwes = [], [], [], []
    n_qual = 0
    for line in p5.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        try:
            q = float(f[0]); an = float(f[1])
        except ValueError:
            continue
        n_qual += 1
        quals.append(q)
        site_miss.append(1.0 - an / (2 * nsamples))
        try:
            mafs.append(float(f[2]))
        except ValueError:
            mafs.append(0.0)
        try:
            hwes.append(float(f[3]))
        except ValueError:
            hwes.append(0.0)
    p5.wait()

    # ---- 过滤漏斗 ----
    n_total = total_sites
    n_pass_qual = n_qual  # 已经过 type+pass+qual
    ok_site = [i for i, m in enumerate(site_miss) if m <= args.site_miss_max]
    ok_maf = [i for i, m in enumerate(mafs) if m >= args.maf_min]
    ok_hwe = [i for i, h in enumerate(hwes) if h >= args.hwe_min]
    keep = set(ok_site) & set(ok_maf) & set(ok_hwe)

    report = {
        "input": args.input,
        "samples": nsamples,
        "total_sites": total_sites,
        "after_type_pass_qual": n_pass_qual,
        "after_site_miss": len(ok_site),
        "after_maf": len(ok_maf),
        "after_hwe": len(ok_hwe),
        "after_all_site_filters": len(keep),
        "density": {c: {"len_bp": contigs.get(c), "n_sites": chrom_counts[c],
                        "sites_per_kb": (chrom_counts[c] / contigs[c] * 1000) if contigs.get(c) else None,
                        "mean_gap_bp": (contigs[c] / chrom_counts[c]) if chrom_counts[c] else None}
                    for c in sorted(chrom_counts)},
        "qual": pct(quals),
        "site_missing": pct(site_miss),
        "maf": pct(mafs),
        "hwe": pct(hwes),
        "gt_dp": {"sampled": len(dp_vals), "pct": pct([float(v) for v in dp_vals]),
                  "below_min_frac": (sum(1 for v in dp_vals if v < args.gt_dp_min) / len(dp_vals)) if dp_vals else None,
                  "above_max_frac": (sum(1 for v in dp_vals if v > args.gt_dp_max) / len(dp_vals)) if dp_vals else None},
        "gt_gq": {"sampled": len(gq_vals), "pct": pct([float(v) for v in gq_vals]),
                  "below_min_frac": (sum(1 for v in gq_vals if v < args.gt_gq_min) / len(gq_vals)) if gq_vals else None},
        "params": {"min_qual": args.min_qual, "gt_dp_min": args.gt_dp_min, "gt_dp_max": args.gt_dp_max,
                   "gt_gq_min": args.gt_gq_min, "site_miss_max": args.site_miss_max,
                   "maf_min": args.maf_min, "hwe_min": args.hwe_min},
    }

    # ---- markdown ----
    md = []
    md.append("# 输入画像 (profile)\n")
    md.append(f"- 输入: `{args.input}`")
    md.append(f"- 样本数: {nsamples}")
    md.append(f"- 位点总数: {total_sites:,}\n")
    md.append("## 密度 (每条染色体)\n")
    md.append("| Chr | 长度(bp) | 位点 | sites/kb | 平均间距(bp) |")
    md.append("|---|---|---|---|---|")
    for c, d in report["density"].items():
        md.append(f"| {c} | {d['len_bp']:,} | {d['n_sites']:,} | "
                  f"{d['sites_per_kb']:.3f} | {d['mean_gap_bp']:.1f} |" if d["len_bp"] else f"| {c} | ? | {d['n_sites']:,} | ? | ? |")
    md.append("")
    md.append("## 分布 (分位数)\n")
    for k in ("qual", "site_missing", "maf", "hwe"):
        md.append(f"- **{k}**: " + ", ".join(f"{a}={b:.4g}" for a, b in report[k].items()))
    md.append(f"- **gt_dp** (抽样 {report['gt_dp']['sampled']}): " +
              ", ".join(f"{a}={b:.4g}" for a, b in report["gt_dp"]["pct"].items()) +
              f" | DP<{args.gt_dp_min} 占 {report['gt_dp']['below_min_frac']:.4f}, "
              f"DP>{args.gt_dp_max} 占 {report['gt_dp']['above_max_frac']:.4f}")
    md.append(f"- **gt_gq** (抽样 {report['gt_gq']['sampled']}): " +
              ", ".join(f"{a}={b:.4g}" for a, b in report["gt_gq"]["pct"].items()) +
              f" | GQ<{args.gt_gq_min} 占 {report['gt_gq']['below_min_frac']:.4f}")
    md.append("")
    md.append("## 过滤漏斗\n")
    md.append(f"- 位点总数: {n_total:,}")
    md.append(f"- 过 type+PASS+QUAL≥{args.min_qual}: {n_pass_qual:,} ({n_pass_qual/n_total:.4%})")
    md.append(f"- 过 位点缺失≤{args.site_miss_max}: {len(ok_site):,} ({len(ok_site)/n_total:.4%})")
    md.append(f"- 过 MAF≥{args.maf_min}: {len(ok_maf):,} ({len(ok_maf)/n_total:.4%})")
    md.append(f"- 过 HWE≥{args.hwe_min}: {len(ok_hwe):,} ({len(ok_hwe)/n_total:.4%})")
    md.append(f"- **全部位点层累积保留**: {len(keep):,} ({len(keep)/n_total:.4%})")
    md.append("")

    text = "\n".join(md)
    print(text)
    with open(os.path.join(args.outdir, "profile.md"), "w") as fh:
        fh.write(text + "\n")
    with open(os.path.join(args.outdir, "profile.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
