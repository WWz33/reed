//! R7 acceptance: the per-site/per-sample reports (02-cli.md §5.4) with
//! hand-computed values. The plink2 parity runs live in `gold/popgen/`.

use std::io::Write;

use reed::reports::{StatsReportRequest, run_stats};
use reed_engine::budget::{MemBudget, ThreadBudget};

fn bgzf(text: &str) -> Vec<u8> {
    let mut deflater =
        reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6 is valid");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(64) {
        deflater
            .deflate_block(block, &mut out)
            .expect("invariant: block fits");
    }
    out.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
    out
}

/// 4 samples; REF/ALT from the file. Sites:
///   pos 1: s0=0/1 s1=0/0 s2=0/0 s3=0/1  → AC=2 AN=8 hom_ref=2 het=2 hom_alt=0
///   pos 5: ALT="." (invariant), all 0/0 → AC=0 AN=8 hom_ref=4
///   pos 25: s0=0/1 s1=1/1 s2=1/1 s3=0/0 → AC=5 AN=8 hom_ref=1 het=1 hom_alt=2
const FIXTURE: &str = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=30>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:10\t0/0:20\t0/0:30\t0/1:40
chr1\t5\t.\tA\t.\t50\tPASS\t.\tGT:DP\t0/0:5\t0/0:.\t0/0:15\t0/0:25
chr1\t25\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:.\t1/1:50\t1/1:60\t0/0:70
";

fn fixture() -> (std::path::PathBuf, std::path::PathBuf) {
    let dir = std::env::temp_dir().join(format!("reed-reports-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("in.vcf.gz");
    let mut file = std::fs::File::create(&vcf).expect("invariant: fixture file");
    file.write_all(&bgzf(FIXTURE)).expect("invariant: write");
    drop(file);
    (dir, vcf)
}

fn lines(path: &std::path::Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .expect("invariant: report file")
        .lines()
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

#[test]
fn reports_match_hand_computed_values() {
    let (dir, vcf) = fixture();
    let prefix = dir.join("out");
    let request = StatsReportRequest {
        freq: true,
        hardy: true,
        het: true,
        missing: true,
        depth: true,
        summary: false,
        hist: Vec::new(),
        singletons: false,
        populations: None,
        force_samples: false,
        targets: Vec::new(),
        targets_exclude: Vec::new(),
        prefix: prefix.clone(),
    };
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::default_half_memory(),
    };
    let outcome = run_stats(Some(&vcf), &request, &budget).expect("invariant: reports run");
    assert_eq!(outcome.n_records, 3);
    assert_eq!(
        outcome.files.len(),
        7,
        "afreq, hardy, vmiss, ldepth, het, smiss, idepth"
    );

    // afreq: ALT_FREQS = AC/AN, OBS_CT = AN.
    let afreq = lines(std::path::Path::new(&format!("{}.afreq", prefix.display())));
    assert_eq!(afreq[0][0], "#CHROM");
    assert_eq!(&afreq[1][1..], ["chr1_1", "A", "G", "0.25", "8"]);
    assert_eq!(&afreq[2][1..], ["chr1_5", "A", ".", "0", "8"]);
    assert_eq!(&afreq[3][1..], ["chr1_25", "A", "G", "0.625", "8"]);

    // hardy: A1 = REF, HOM_A1_CT = hom-ref, TWO_AX_CT = hom-alt,
    // O(HET) = het/called, E(HET) = 2 p_ref (1 - p_ref).
    let hardy = lines(std::path::Path::new(&format!("{}.hardy", prefix.display())));
    assert_eq!(hardy[0][0], "#CHROM");
    assert_eq!(&hardy[1][2..9], ["A", "G", "2", "2", "0", "0.5", "0.375"]);
    let p1: f64 = hardy[1][9].parse().expect("HWE p");
    // n=4 called, mac=2, het=2 → the observed configuration is the most
    // likely one, p = 1.
    assert!((p1 - 1.0).abs() < 1e-12, "{p1}");
    assert_eq!(&hardy[3][2..8], ["A", "G", "1", "1", "2", "0.25"]);
    assert!((hardy[3][8].parse::<f64>().expect("E(HET)") - 0.46875).abs() < 1e-6);

    // vmiss: MISSING_CT = fully-missing genotypes, OBS_CT = n_samples.
    let vmiss = lines(std::path::Path::new(&format!("{}.vmiss", prefix.display())));
    assert_eq!(vmiss[0][0], "#CHROM");
    assert_eq!(&vmiss[1][2..], ["0", "4", "0"]);

    // het: plink's moment F per sample over the three sites.
    // Σ 2p(1−p) over called sites = 0.375 + 0 + 0.46875 = 0.84375 for every
    // sample (all sites called here).
    let het = lines(std::path::Path::new(&format!("{}.het", prefix.display())));
    assert_eq!(&het[0][..], ["#IID", "O(HOM)", "E(HOM)", "OBS_CT", "F"]);
    // s0: one hom call out of three → F = (1 − 2.15625)/0.84375.
    assert_eq!(&het[1][..3], ["s0", "1", "2.15625"]);
    // Files print %g with 6 significant digits, so compare at that precision.
    assert!((het[1][4].parse::<f64>().expect("F") - (1.0 - 2.15625) / 0.84375).abs() < 1e-4);
    // s1: three hom calls → F = (3 − 2.15625)/0.84375.
    assert_eq!(&het[2][..3], ["s1", "3", "2.15625"]);
    assert!((het[2][4].parse::<f64>().expect("F") - (3.0 - 2.15625) / 0.84375).abs() < 1e-4);

    // smiss: no missing genotypes; OBS_CT is the record count (plink2).
    let smiss = lines(std::path::Path::new(&format!("{}.smiss", prefix.display())));
    assert_eq!(&smiss[0][..], ["#IID", "MISSING_CT", "OBS_CT", "F_MISS"]);
    assert_eq!(&smiss[1][..], ["s0", "0", "3", "0"]);

    // idepth (vcftools): N_SITES = sites where the GT is called and DP is
    // present; MEAN_DEPTH over those.
    // s0: sites 1 (DP 10) and 5 (DP 5); site 25 has DP "." → 2 sites, mean 7.5.
    // s1: sites 1 (20) and 25 (50); site 5 has DP "." → mean 35.
    // s2: 3 sites (30 + 15 + 60 = 105) → 35. s3: 3 sites (135) → 45.
    let idepth = lines(std::path::Path::new(&format!(
        "{}.idepth",
        prefix.display()
    )));
    assert_eq!(&idepth[0][..], ["INDV", "N_SITES", "MEAN_DEPTH"]);
    assert_eq!(&idepth[1][..], ["s0", "2", "7.5"]);
    assert_eq!(&idepth[2][..], ["s1", "2", "35"]);
    assert_eq!(&idepth[3][..], ["s2", "3", "35"]);
    assert_eq!(&idepth[4][..], ["s3", "3", "45"]);

    // ldepth.mean (vcftools): per-site mean and sample variance (n−1).
    // Site 1: 10,20,30,40 → mean 25, var 500/3.
    // Site 5: s1 excluded → 5,15,25 → mean 15, var 100.
    // Site 25: s0 excluded → 50,60,70 → mean 60, var 100.
    let ldepth = lines(std::path::Path::new(&format!(
        "{}.ldepth.mean",
        prefix.display()
    )));
    assert_eq!(&ldepth[0][..], ["CHROM", "POS", "MEAN_DEPTH", "VAR_DEPTH"]);
    assert_eq!(&ldepth[1][..3], ["chr1", "1", "25"]);
    assert!((ldepth[1][3].parse::<f64>().expect("var") - 500.0 / 3.0).abs() < 1e-3);
    assert_eq!(&ldepth[2][..], ["chr1", "5", "15", "100"]);
    assert_eq!(&ldepth[3][..], ["chr1", "25", "60", "100"]);

    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v7 D7-7 fixture: SNPs, insertions, deletions, a missing QUAL and a
/// missing genotype, over two samples (matches `bcftools stats -s -`).
const D77: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\t.\tA\tG\t50\tPASS\tAC=1;AN=4\tGT\t0/1\t0/0
chr1\t20\t.\tA\tAT\t60\tPASS\tAC=1;AN=4\tGT\t0/1\t0/0
chr1\t30\t.\tAT\tA\t.\tPASS\tAC=1;AN=4\tGT\t0/1\t0/0
chr1\t40\t.\tATGC\tA\t40\tPASS\tAC=2;AN=4\tGT\t1/1\t0/0
chr1\t50\t.\tA\tAGG\t50\tPASS\tAC=1;AN=4\tGT\t.\t0/1
chr1\t60\t.\tACGT\tA\t.\tPASS\tAC=1;AN=4\tGT\t0/1\t./.
chr1\t70\t.\tA\tT\t30\tPASS\tAC=2;AN=4\tGT\t1/1\t1/1
";

/// doc_v7 D7-7: IDD, QUAL and PSC must match `bcftools stats -s -`. The
/// expected lines are byte-for-byte `bcftools` 1.24 output on this fixture
/// (verified); IDD lengths are signed, the QUAL missing bin is `.`, and PSC
/// carries the real sample names with nIndels and SNP-gated genotype
/// classes.
#[test]
fn summary_idd_qual_psc_match_bcftools() {
    let dir = std::env::temp_dir().join(format!("reed-reports-d77-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("in.vcf.gz");
    let mut file = std::fs::File::create(&vcf).expect("invariant: fixture file");
    file.write_all(&bgzf(D77)).expect("invariant: write");
    drop(file);
    let prefix = dir.join("out");
    let request = StatsReportRequest {
        summary: true,
        prefix: prefix.clone(),
        ..StatsReportRequest::default()
    };
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let outcome = run_stats(Some(&vcf), &request, &budget).expect("invariant: run");
    let qc = outcome.qc.expect("invariant: summary requested");
    let names = vec!["S1".to_string(), "S2".to_string()];
    let text = reed::summary::render_bcftools_stats(&qc, &names);
    for line in [
        "QUAL\t0\t.\t0\t0\t0\t2",
        "QUAL\t0\t30.0\t1\t0\t1\t0",
        "QUAL\t0\t40.0\t0\t0\t0\t1",
        "QUAL\t0\t50.0\t1\t1\t0\t1",
        "QUAL\t0\t60.0\t0\t0\t0\t1",
        "IDD\t0\t-3\t2\t0\t.",
        "IDD\t0\t-1\t1\t0\t.",
        "IDD\t0\t1\t1\t0\t.",
        "IDD\t0\t2\t1\t0\t.",
        "PSC\t0\tS1\t0\t1\t1\t1\t1\t4\t0.0\t5\t0\t0\t1",
        "PSC\t0\tS2\t4\t1\t0\t0\t1\t1\t0.0\t1\t0\t0\t1",
    ] {
        assert!(text.contains(line), "missing {line:?}\n{text}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}
