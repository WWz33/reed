#!/usr/bin/env bash
# gold/popgen/run_convert_parity.sh — `reed convert` parity against the
# reference tools (docs/rounds/v13/01-review.md §7).
#
# Usage:
#   bash gold/popgen/run_convert_parity.sh                 # fixture gates
#   bash gold/popgen/run_convert_parity.sh DATA.vcf.gz     # + full-data phylo
#
# The optional DATA argument adds the full-data phylogenetic gate: one
# vcf2phylip run (PHYLIP/FASTA/NEXUS/bin/used_sites, no outgroup) against one
# `reed convert` run over the same input, every file byte-equal. It needs
# both tools' runtime and is not run by default.
#
# Overridable: REED, REED_SCRATCH, VCF2PHYLIP, VCFTOOLS, TASSEL, VCF2MSA,
# PYTHON, BGZIP, TABIX, HTSLIB_DIR.
#
# Gates:
#   1. vcf2phylip 2.9 on gold/fixtures/phylo.vcf: PHYLIP, FASTA, NEXUS,
#      binary NEXUS and used_sites byte-equal (--outgroup s3, default gate).
#   2. vcftools 0.1.18 `--012`: the full-diploid rows byte-equal; the B13-8
#      rows (haploid/half call) asserted as reed `-1` vs vcftools `1`.
#   3. TASSEL 5 reads reed's `.hmp.txt` back (import + re-export succeeds,
#      taxa and site counts preserved).
#   4. vcf2msa on gold/fixtures/convert_msa.vcf: sequences equal with the
#      region convention mapped (vcf2msa `B-(E+1)` == reed `B-E`; vcf2msa's
#      end bound is exclusive, reed's is the natural inclusive one).
#   5. `reed convert --sfs` == `reed stats --sfs` byte-for-byte on the same
#      prefix, input and grouping (v13 §4.1).
#   6. hapmap ↔ VCF round trip (v13 §3.12): VCF → hapmap → VCF reproduces
#      the records on the lossless domain (biallelic diploid SNPs, no
#      missing), and TASSEL 5 reads reed's VCF back.
set -euo pipefail

DATA="${1:-}"
GOLD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$GOLD_DIR/.." && pwd)"
SCRATCH="${REED_SCRATCH:-$REPO_ROOT/target/gold}/convert"
REED="${REED:-/tmp/reed-target/release/reed}"
VCF2PHYLIP="${VCF2PHYLIP:-$HOME/reed/vcf2phylip/vcf2phylip.py}"
VCFTOOLS="${VCFTOOLS:-$HOME/biosoft/vcftools}"
TASSEL="${TASSEL:-$HOME/TASSEL5/run_pipeline.pl}"
VCF2MSA="${VCF2MSA:-$HOME/reed/vcf2msa.py/vcf2msa.py}"
PYTHON="${PYTHON:-python3}"
HTSLIB_DIR="${HTSLIB_DIR:-$HOME/biosoft/htslib-1.23}"
BGZIP="${BGZIP:-$(command -v bgzip || echo "$HTSLIB_DIR/bgzip")}"
TABIX="${TABIX:-$(command -v tabix || echo "$HTSLIB_DIR/tabix")}"

[[ -x "$REED" ]] || { echo "reed binary not executable: $REED" >&2; exit 1; }
[[ -f "$VCF2PHYLIP" ]] || { echo "vcf2phylip not found: $VCF2PHYLIP" >&2; exit 1; }
[[ -x "$VCFTOOLS" ]] || { echo "vcftools not found: $VCFTOOLS" >&2; exit 1; }

mkdir -p "$SCRATCH"
fail=0
note() { echo "== $* =="; }
check() { # check LABEL FILE_A FILE_B
    if diff -q "$2" "$3" >/dev/null 2>&1; then
        echo "  SAME  $1"
    else
        echo "  DIFF  $1" >&2
        diff "$2" "$3" | head -10 >&2 || true
        fail=1
    fi
}

# ---------------------------------------------------------------- 1. phylo --
note "vcf2phylip 2.9 vs reed convert (PHYLIP/FASTA/NEXUS/bin/used_sites)"
PHYLO_VCF="$GOLD_DIR/fixtures/phylo.vcf"
REF="$SCRATCH/vcf2phylip"
GOT="$SCRATCH/reed"
rm -rf "$REF" "$GOT" && mkdir -p "$REF" "$GOT"
$PYTHON "$VCF2PHYLIP" -i "$PHYLO_VCF" -f -n -b -w -o s3 --output-folder "$REF" >"$SCRATCH/vcf2phylip.log" 2>&1
"$REED" convert "$PHYLO_VCF" --phylip --fasta --nexus --nexus-binary \
    --outgroup s3 --write-used-sites -o "$GOT/out" --verbosity 0
for f in phy fasta nexus; do
    check "$f" "$GOT/out.$f" "$REF/phylo.min4.$f"
done
check "bin.nexus" "$GOT/out.bin.nexus" "$REF/phylo.min4.bin.nexus"
check "used_sites.tsv" "$GOT/out.used_sites.tsv" "$REF/phylo.min4.used_sites.tsv"
# The counters mirror vcf2phylip's stdout report.
grep -q "counters: total_records=20 excess_missing=0 mnps=2 accepted=17 biallelic_for_bin=14 malformed=1" "$GOT/out.log" \
    && echo "  SAME  counters" || { echo "  DIFF  counters" >&2; fail=1; }

# ------------------------------------------------------------------- 2. 012 --
note "vcftools --012 vs reed convert --012"
C012_VCF="$GOLD_DIR/fixtures/convert012.vcf"
CLEAN="$SCRATCH/convert012_clean.vcf"
# Drop the B13-8 rows (pos 5: haploid and half calls) for the diff gate.
grep -v -P '^chr1\t5\t' "$C012_VCF" >"$CLEAN"
"$VCFTOOLS" --vcf "$CLEAN" --012 --out "$SCRATCH/v012" >/dev/null 2>&1
"$REED" convert "$CLEAN" --012 -o "$SCRATCH/r012" --verbosity 0
check ".012" "$SCRATCH/r012.012" "$SCRATCH/v012.012"
check ".012.pos" "$SCRATCH/r012.012.pos" "$SCRATCH/v012.012.pos"
check ".012.indv" "$SCRATCH/r012.012.indv" "$SCRATCH/v012.012.indv"
# B13-8: the excluded row prints reed -1 where vcftools prints 1.
"$REED" convert "$C012_VCF" --012 -o "$SCRATCH/r012_full" --verbosity 0
if [[ "$(grep -P '^0\t' "$SCRATCH/r012_full.012" | cut -f5)" == "-1" ]]; then
    echo "  SAME  B13-8 (haploid/half call -> -1)"
else
    echo "  DIFF  B13-8 row" >&2; fail=1
fi

# ---------------------------------------------------------------- 3. hapmap --
if [[ -x "$TASSEL" ]] && command -v java >/dev/null 2>&1; then
    note "TASSEL 5 read-back of reed's hapmap"
    "$REED" convert "$PHYLO_VCF" --hapmap -o "$SCRATCH/hmp" --verbosity 0
    mkdir -p "$SCRATCH/hmp_back"
    if timeout 300 "$TASSEL" -Xmx2g -importGuess "$SCRATCH/hmp.hmp.txt" \
        -export "$SCRATCH/hmp_back/back" -exportType Hapmap >"$SCRATCH/tassel.log" 2>&1 \
        && [[ -s "$SCRATCH/hmp_back/back.hmp.txt" ]]; then
        want=$(($(grep -c "" "$SCRATCH/hmp.hmp.txt") - 1))
        got=$(($(grep -c "" "$SCRATCH/hmp_back/back.hmp.txt") - 1))
        if [[ "$want" == "$got" ]]; then
            echo "  SAME  TASSEL read-back ($got sites)"
        else
            echo "  DIFF  TASSEL read-back: $want sites in, $got out" >&2; fail=1
        fi
    else
        echo "  DIFF  TASSEL could not read reed's hapmap" >&2; fail=1
    fi
else
    note "TASSEL not available; hapmap read-back gate skipped"
fi

# ------------------------------------------------------------------- 4. msa --
if [[ -f "$VCF2MSA" ]] && [[ -x "$BGZIP" ]] && [[ -x "$TABIX" ]]; then
    if $PYTHON -c "import vcf, pysam" >/dev/null 2>&1; then
        note "vcf2msa vs reed convert --msa (region B-E(reed) == B-(E+1)(vcf2msa))"
        MSA_VCF="$GOLD_DIR/fixtures/convert_msa.vcf"
        MSA_REF="$GOLD_DIR/fixtures/convert_msa.fa"
        "$BGZIP" -c "$MSA_VCF" >"$SCRATCH/convert_msa.vcf.gz"
        "$TABIX" -f -p vcf "$SCRATCH/convert_msa.vcf.gz"
        rm -f "$SCRATCH"/contig_chr1_1-11.fasta
        (cd "$SCRATCH" && $PYTHON "$VCF2MSA" -f "$MSA_REF" -v convert_msa.vcf.gz -r chr1:1-11 >/dev/null 2>&1)
        "$REED" convert "$MSA_VCF" --msa --reference "$MSA_REF" --region chr1:1-10 \
            -o "$SCRATCH/msa" --verbosity 0
        $PYTHON - "$SCRATCH/contig_chr1_1-11.fasta" "$SCRATCH/msa_chr1_1-10.fasta" <<'PY'
import sys

def seqs(path):
    out = {}
    for block in open(path).read().strip().split(">"):
        if not block.strip():
            continue
        head, seq = block.strip().split("\n")[:2]
        out[head.split("_")[-1]] = seq
    return out

want, got = seqs(sys.argv[1]), seqs(sys.argv[2])
if want != got:
    print(f"  DIFF  msa: {want} vs {got}", file=sys.stderr)
    sys.exit(1)
print(f"  SAME  msa ({len(got)} samples)")
PY
    else
        note "vcf2msa python deps missing (pyvcf3/pysam); msa gate skipped"
    fi
else
    note "vcf2msa or htslib tools not available; msa gate skipped"
fi

# ------------------------------------------------------------------- 5. sfs --
note "convert --sfs vs stats --sfs (byte-for-byte)"
SFS_VCF="$GOLD_DIR/fixtures/sfs.vcf"
SFS_POPS="$GOLD_DIR/fixtures/sfs.pops"
rm -rf "$SCRATCH/sfs_a" "$SCRATCH/sfs_bPosition" "$SCRATCH/sfs_b"
mkdir -p "$SCRATCH/sfs_a" "$SCRATCH/sfs_b"
"$REED" convert "$SFS_VCF" --pop "$SFS_POPS" --sfs --proj 6 -o "$SCRATCH/sfs_a/out" --verbosity 0
"$REED" stats "$SFS_VCF" --pop "$SFS_POPS" --sfs --proj 6 -o "$SCRATCH/sfs_b/out" --verbosity 0
if diff -r "$SCRATCH/sfs_a/out_sfs" "$SCRATCH/sfs_b/out_sfs" >/dev/null 2>&1; then
    echo "  SAME  _sfs trees"
else
    echo "  DIFF  _sfs trees" >&2; fail=1
fi
# The existing easySFS gate also exercises the convert entry.
if [[ -f "$GOLD_DIR/popgen/run_easysfs.sh" ]]; then
    bash "$GOLD_DIR/popgen/run_easysfs.sh" "$SFS_VCF" "$SFS_POPS" 6 "$SCRATCH/easysfs" >/dev/null 2>&1 \
        && echo "  SAME  easySFS gate (stats entry)" \
        || { echo "  DIFF  easySFS gate" >&2; fail=1; }
fi

# ------------------------------------------------- 6. hapmap ↔ VCF round trip --
note "hapmap ↔ VCF round trip and TASSEL read-back of reed's VCF"
RT="$SCRATCH/rt"
rm -rf "$RT" && mkdir -p "$RT"
cat >"$RT/in.vcf" <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s0	s1	s2	s3
chr1	10	chr1:10	A	G	.	PASS	.	GT	0/0	0/1	1/1	0/1
chr1	20	chr1:20	C	T	.	PASS	.	GT	0/1	0/0	1/1	1/1
chr1	30	chr1:30	A	C	.	PASS	.	GT	1/1	1/1	0/0	0/1
chr1	40	chr1:40	G	T	.	PASS	.	GT	0/1	1/1	0/1	0/0
VCF
"$REED" convert "$RT/in.vcf" --hapmap -o "$RT/hmp" --verbosity 0
"$REED" convert "$RT/hmp.hmp.txt" --vcf -o "$RT/back.vcf" --verbosity 0
grep -v '^#' "$RT/in.vcf" >"$RT/in.rows"
grep -v '^#' "$RT/back.vcf" >"$RT/back.rows"
check "round-trip records" "$RT/back.rows" "$RT/in.rows"
if head -c2 "$RT/back.vcf" | od -An -tx1 | grep -q "1f 8b"; then
    echo "  DIFF  default -O must be plain text" >&2; fail=1
fi
"$REED" convert "$RT/hmp.hmp.txt" --vcf -Oz -o "$RT/back.vcf.gz" --verbosity 0
if [[ "$(head -c2 "$RT/back.vcf.gz" | od -An -tx1 | tr -d ' \n')" == "1f8b" ]]; then
    echo "  SAME  -Oz is BGZF"
else
    echo "  DIFF  -Oz is not BGZF" >&2; fail=1
fi
if [[ -x "$TASSEL" ]] && command -v java >/dev/null 2>&1; then
    mkdir -p "$RT/tassel_back"
    if timeout 300 "$TASSEL" -Xmx2g -importGuess "$RT/back.vcf" \
        -export "$RT/tassel_back/t" -exportType Hapmap >"$RT/tassel_vcf.log" 2>&1 \
        && [[ -s "$RT/tassel_back/t.hmp.txt" ]]; then
        got=$(($(grep -c "" "$RT/tassel_back/t.hmp.txt") - 1))
        want=$(grep -c "" "$RT/in.rows")
        if [[ "$want" == "$got" ]]; then
            echo "  SAME  TASSEL read-back of reed's VCF ($got sites)"
        else
            echo "  DIFF  TASSEL read reed's VCF: $want sites in, $got out" >&2; fail=1
        fi
    else
        echo "  DIFF  TASSEL could not read reed's VCF" >&2; fail=1
    fi
else
    note "TASSEL not available; VCF read-back gate skipped"
fi

# ------------------------------------------------------------ 7. full data --
if [[ -n "$DATA" ]]; then
    [[ -f "$DATA" ]] || { echo "input not found: $DATA" >&2; exit 1; }
    note "full-data phylo parity: $DATA"
    FREF="$SCRATCH/full_ref"
    FGOT="$SCRATCH/full_got"
    rm -rf "$FREF" "$FGOT" && mkdir -p "$FREF" "$FGOT"
    # NOTE: vcf2phylip's -o is --outgroup; the prefix comes from
    # --output-prefix (``out`` keeps the two sides' names aligned).
    ( time $PYTHON "$VCF2PHYLIP" -i "$DATA" -f -n -b -w \
        --output-prefix out --output-folder "$FREF" ) \
        >"$SCRATCH/full_ref.log" 2>&1
    ( time "$REED" convert "$DATA" --phylip --fasta --nexus --nexus-binary \
        --write-used-sites -o "$FGOT/out" --verbosity 0 ) >"$SCRATCH/full_got.log" 2>&1
    for f in phy fasta nexus; do
        check "full/$f" "$FGOT/out.$f" "$FREF/out.min4.$f"
    done
    check "full/bin.nexus" "$FGOT/out.bin.nexus" "$FREF/out.min4.bin.nexus"
    check "full/used_sites.tsv" "$FGOT/out.used_sites.tsv" "$FREF/out.min4.used_sites.tsv"
    echo "  vcf2phylip: $(grep -E '^real' "$SCRATCH/full_ref.log" | head -1)"
    echo "  reed:       $(grep -E '^real' "$SCRATCH/full_got.log" | head -1)"
fi

if [[ "$fail" -eq 0 ]]; then
    echo "convert parity: PASS"
else
    echo "convert parity: FAIL" >&2
    exit 1
fi