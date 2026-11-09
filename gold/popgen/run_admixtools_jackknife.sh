#!/usr/bin/env bash
# V27-35: ADMIXTOOLS observed-site block grid and weighted jackknife parity.
#
# reed retains its natural Σnum/Σden estimate; qpDstat reports a bias-corrected
# estimate. The diagnostic reference build adds only this print after the
# unchanged getdsc weightjack call: REED_JK natural=<ymean> se=<sig>.
# Without such a qpDstat this exits 2 (unverified), never a fake PASS.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REED="${REED:-/tmp/reed-target/release/reed}"
DIAG="${ADMIXTOOLS_DIAG_DIR:-/tmp/reed-v27-admix-instrumented}"
OUT="${REED_SCRATCH:-$ROOT/target/gold}/out/admixtools_jackknife"
[[ -x "$REED" ]] || { echo "reed not executable: $REED" >&2; exit 1; }
[[ -x "$DIAG/qpDstat" ]] || {
  echo "no full-precision ADMIXTOOLS qpDstat in $DIAG" >&2
  echo "set ADMIXTOOLS_DIAG_DIR; SE parity is NOT verified" >&2
  exit 2
}
mkdir -p "$OUT"
python3 -I - "$OUT" <<'PY'
from pathlib import Path
import sys
r = Path(sys.argv[1])
samples = [f"{p}{i}" for p in "ABCD" for i in range(4)]
lines = ["##fileformat=VCFv4.3", "##contig=<ID=chr1,length=60000000>",
         "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">",
         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
with (r / "input.snp").open("w") as snp, (r / "input.geno").open("w") as geno:
    for block in range(6):
        for i in range(8 if block == 0 else 2):
            pos = 1 + block * 10_000_000 + i * 100
            gt = [(block * 17 + i * 11 + j * 7) % 19 for j in range(len(samples))]
            gt = ["9" if x == 0 else str(x % 3) for x in gt]
            snp.write(f"rs{block * 8 + i}\t1\t0.0\t{pos}\tA\tG\n")
            geno.write("".join(gt) + "\n")
            decode = {"0": "0/0", "1": "0/1", "2": "1/1", "9": "./."}
            lines.append(f"chr1\t{pos}\trs{block * 8 + i}\tA\tG\t50\tPASS\t.\tGT\t" +
                         "\t".join(decode[x] for x in gt))
(r / "input.vcf").write_text("\n".join(lines) + "\n")
with (r / "pops.txt").open("w") as out, (r / "input.ind").open("w") as ind:
    for sample in samples:
        out.write(f"{sample}\t{sample[0]}\n")
        ind.write(f"{sample}\tU\t{sample[0]}\n")
(r / "d.pop").write_text("A B C D\n")
(r / "d.par").write_text(
    f"indivname: {r}/input.ind\nsnpname: {r}/input.snp\n"
    f"genotypename: {r}/input.geno\npopfilename: {r}/d.pop\n"
    "blgsize: 0.05\ngfromp: YES\nfewblocksok: YES\n")
PY

"$DIAG/qpDstat" -p "$OUT/d.par" > "$OUT/d.log" 2> "$OUT/d.err"
"$REED" popg "$OUT/input.vcf" --pop "$OUT/pops.txt" --dstat A,B,C,D \
  --block-size 5M -o "$OUT/reed" --verbosity 0
cat "$OUT/reed.fstats"
grep 'REED_JK' "$OUT/d.err"
grep 'number of blocks' "$OUT/d.log"
python3 -I - "$OUT" <<'PY'
from pathlib import Path
import sys
r = Path(sys.argv[1])
row = (r / "reed.fstats").read_text().splitlines()[1].split("\t")
est, se, blocks, sites = float(row[2]), float(row[3]), int(row[5]), int(row[6])
line = next(x for x in (r / "d.err").read_text().splitlines()
            if x.startswith("REED_JK natural="))
ref = dict((k, float(v)) for k, v in (item.split("=") for item in line.split()[1:]))
assert abs(est - ref["natural"]) < 2e-12, (est, ref["natural"])
assert abs(se - ref["se"]) < 2e-12, (se, ref["se"])
assert blocks == 6 and sites == 18
print("PASS ADMIXTOOLS: 6 observed-site blocks, D natural estimate and weighted SE")
PY
