#!/usr/bin/env python3
"""gold/scripts/check_parity.py — 对拍工具：把任意「过滤后 VCF」与黄金标准 VCF 逐位点比对。

用法:
    python3 check_parity.py <candidate.vcf[.gz|.bcf]> [--gold PATH] [--json] [--report FILE]

判定:
    PASS  位点集合 (CHROM,POS,REF,ALT) 完全一致 且 同名样本基因型零不一致
    FAIL  否则 (退出码 1)；参数/读取错误退出码 2

这是给未来 reed 自带过滤器层用的验收门：把 reed 过滤输出喂进来，PASS 才算与
bcftools/vcftools/plink2 三方金标准一致。
"""
from __future__ import annotations
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GOLD_DIR = os.path.dirname(HERE)


def bcftools() -> str:
    return os.environ.get("BCFTOOLS", "bcftools")


def run(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(args, capture_output=True, text=True)


def site_keys(vcf: str) -> set[str]:
    p = run([bcftools(), "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", vcf])
    if p.returncode != 0:
        raise RuntimeError(f"bcftools query 失败 on {vcf}: {p.stderr[:300]}")
    out = set()
    for line in p.stdout.splitlines():
        if line:
            c, pos, r, a = line.split("\t")
            out.add(f"{c}:{pos}:{r}:{a}")
    return out


def samples(vcf: str) -> list[str]:
    p = run([bcftools(), "query", "-l", vcf])
    if p.returncode != 0:
        raise RuntimeError(f"bcftools query -l 失败 on {vcf}: {p.stderr[:300]}")
    return [s for s in p.stdout.splitlines() if s]


def gtcheck(ref: str, qry: str) -> tuple[float | None, float | None, str]:
    """同名样本 GT 不一致 (bcftools gtcheck -u GT,GT)。返回 (mean, max, err)。"""
    import tempfile
    samps = samples(ref)
    if not samps:
        return None, None, "no samples in gold"
    with tempfile.NamedTemporaryFile("w", suffix=".pairs", delete=False) as fh:
        for s in samps:
            fh.write(f"{s}\t{s}\n")
        pairs = fh.name
    try:
        p = run([bcftools(), "gtcheck", "-u", "GT,GT", "-g", ref, "-P", pairs, qry])
    finally:
        os.unlink(pairs)
    if p.returncode != 0:
        return None, None, p.stderr.strip()[:200]
    vals = []
    for line in p.stdout.splitlines():
        if line.startswith("DCv2"):
            f = line.split("\t")
            try:
                vals.append(float(f[3]))
            except (IndexError, ValueError):
                pass
    if not vals:
        return None, None, "no DCv2 lines"
    return sum(vals) / len(vals), max(vals), ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("candidate", help="待验收的过滤后 VCF")
    ap.add_argument("--gold", default=os.path.join(GOLD_DIR, "out", "gold.vcf.gz"))
    ap.add_argument("--json", action="store_true", help="额外输出机器可读 JSON")
    ap.add_argument("--report", default=None, help="把报告写入文件")
    args = ap.parse_args()

    if not os.path.exists(args.gold):
        print(f"黄金标准不存在: {args.gold}\n先运行 gold/scripts/run_all.sh", file=sys.stderr)
        return 2
    if not os.path.exists(args.candidate):
        print(f"候选文件不存在: {args.candidate}", file=sys.stderr)
        return 2

    try:
        gold_sites = site_keys(args.gold)
        cand_sites = site_keys(args.candidate)
        gold_samples = samples(args.gold)
        cand_samples = samples(args.candidate)
    except RuntimeError as e:
        print(str(e), file=sys.stderr)
        return 2

    common = gold_sites & cand_sites
    gold_only = gold_sites - cand_sites
    cand_only = cand_sites - gold_sites
    union = gold_sites | cand_sites
    jac = len(common) / len(union) if union else 1.0

    sample_ok = gold_samples == cand_samples
    missing_s = [s for s in gold_samples if s not in set(cand_samples)]
    extra_s = [s for s in cand_samples if s not in set(gold_samples)]

    disc = disc_max = None
    disc_err = ""
    if common and sample_ok:
        disc, disc_max, disc_err = gtcheck(args.gold, args.candidate)

    site_ok = (not gold_only) and (not cand_only)
    gt_ok = (disc is not None and disc == 0.0 and disc_max == 0.0)
    passed = site_ok and sample_ok and gt_ok

    rep = []
    rep.append("=== gold parity report ===")
    rep.append(f"gold      : {args.gold}")
    rep.append(f"candidate : {args.candidate}")
    rep.append(f"gold sites: {len(gold_sites):,}   candidate sites: {len(cand_sites):,}")
    rep.append(f"common    : {len(common):,}   gold-only: {len(gold_only):,}   cand-only: {len(cand_only):,}   jaccard: {jac:.6f}")
    rep.append(f"samples   : {'MATCH' if sample_ok else 'DIFFER'}  (gold={len(gold_samples)}, cand={len(cand_samples)}"
               + (f", missing={len(missing_s)}, extra={len(extra_s)}" if not sample_ok else "") + ")")
    if common and sample_ok:
        rep.append(f"genotype  : {'n/a ('+disc_err+')' if disc is None else f'mean={disc:.4g} max={disc_max:.4g}'}")
    rep.append(f"VERDICT   : {'PASS' if passed else 'FAIL'}")
    text = "\n".join(rep)
    print(text)
    if args.report:
        with open(args.report, "w") as fh:
            fh.write(text + "\n")
    if args.json:
        print(json.dumps({
            "gold": args.gold, "candidate": args.candidate,
            "gold_sites": len(gold_sites), "candidate_sites": len(cand_sites),
            "common": len(common), "gold_only": len(gold_only), "cand_only": len(cand_only),
            "jaccard": jac, "samples_match": sample_ok,
            "missing_samples": missing_s, "extra_samples": extra_s,
            "genotype_discordance": disc, "verdict": "PASS" if passed else "FAIL",
        }, indent=2))

    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
