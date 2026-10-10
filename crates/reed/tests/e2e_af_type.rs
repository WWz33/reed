//! `-q/-Q` `:type` suffix acceptance (docs/rounds/v5/02-cli.md §3.2;
//! bcftools `vcfview.c` semantics): the five allele-frequency definitions,
//! with INFO/AC+AN preferred over the genotype tallies. The expected values
//! come from bcftools on the same fixture; the parity run lives in
//! `gold/popgen/run_bcftools_af.sh`.

use reed::spec::{AfKind, FilterSpec};
use reed::{FilterOptions, run_filter};

/// 4 samples; sites: pos 1 has AC=2 (INFO stale, AN=8), pos 2 AC=6,
/// pos 3 has no INFO/AC (genotype-only), pos 4 is monomorphic.
const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
chr1\t1\t.\tA\tG\t50\tPASS\tAC=3;AN=8\tGT\t0/0\t0/0\t0/0\t0/1
chr1\t2\t.\tA\tG\t50\tPASS\tAC=6;AN=8\tGT\t1/1\t0/1\t1/1\t0/1
chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1
chr1\t4\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0
";

fn keep(tag: &str, kind: AfKind, min_af: f64) -> Vec<u64> {
    let dir = std::env::temp_dir().join(format!("reed-aftype-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("in.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let out = dir.join("out.vcf");
    let mut spec = FilterSpec::default();
    spec.site_after_mask.min_af = Some(min_af);
    spec.site_after_mask.min_af_kind = Some(kind);
    let opts = FilterOptions {
        threads: 1,
        output: Some(out.clone()),
        ..FilterOptions::default()
    };
    run_filter(Some(&vcf), &spec, &opts).expect("invariant: filter runs");
    let text = std::fs::read_to_string(&out).expect("invariant: output");
    let positions: Vec<u64> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(1).expect("POS").parse().expect("pos"))
        .collect();
    let _ = std::fs::remove_dir_all(&dir);
    positions
}

#[test]
fn af_types_use_info_counts_and_the_bcftools_definitions() {
    // The expected sets are bcftools' own on this fixture (`bcftools view
    // -q X:type`), which prefers INFO/AC+AN over the genotype tallies.
    for (tag, kind) in [("nref", AfKind::Nref), ("alt1", AfKind::Alt1)] {
        assert_eq!(keep(tag, kind, 0.3), vec![1, 2, 3], "{tag} 0.3");
        assert_eq!(keep(tag, kind, 0.5), vec![2, 3], "{tag} 0.5");
        assert_eq!(keep(tag, kind, 0.6), vec![2], "{tag} 0.6");
    }
    // minor = min over REF and ALT counts (INFO first, genotype fallback).
    assert_eq!(keep("minor", AfKind::Minor, 0.3), vec![1, 3]);
    assert_eq!(keep("minor", AfKind::Minor, 0.5), vec![3]);
    assert_eq!(keep("minor", AfKind::Minor, 0.6), [] as [u64; 0]);
    // major = max over both counts.
    assert_eq!(keep("major", AfKind::Major, 0.3), vec![1, 2, 3, 4]);
    assert_eq!(keep("major", AfKind::Major, 0.5), vec![1, 2, 3, 4]);
    assert_eq!(keep("major", AfKind::Major, 0.6), vec![1, 2, 4]);
    // nonmajor = 1 − major.
    assert_eq!(keep("nonmajor", AfKind::Nonmajor, 0.3), vec![1, 3]);
    assert_eq!(keep("nonmajor", AfKind::Nonmajor, 0.5), vec![3]);
    assert_eq!(keep("nonmajor", AfKind::Nonmajor, 0.6), [] as [u64; 0]);
}
