//! Projected-SFS acceptance on the shared fixture (docs/rounds/v5/05
//! §3.1): 4 diploid samples, 11 sites with ALT counts 1,1,2,2,3,3,4,4,5,6,7
//! of 8 chromosomes, projected to 6. The expectations are the exact
//! hypergeometric fractions (1/28 granularity), computed by hand; the
//! fastsimcoal2/dadi file layouts are checked against easySFS's.

use std::path::Path;

use reed::budget::{MemBudget, ThreadBudget};
use reed::{run_stats, Levels, SfsRequest, StatsRequest};

const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1
chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0
chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/1
chr1\t4\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0
chr1\t5\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/0
chr1\t6\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t0/1\t0/0\t0/0
chr1\t7\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t0/1\t0/1\t0/0
chr1\t8\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t1/1\t0/0\t0/0
chr1\t9\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t1/1\t1/0\t0/0
chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t1/1\t1/1\t0/0
chr1\t11\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t1/1\t1/1\t1/0
";

/// ξ_k for k = 0..=6 as exact fractions /28 (the hypergeometric projection
/// of ALT counts 1,1,2,2,3,3,4,4,5,6,7 from 8 to 6 chromosomes).
const EXPECTED: [f64; 7] = [
    16.0 / 28.0,
    72.0 / 28.0,
    72.0 / 28.0,
    62.0 / 28.0,
    42.0 / 28.0,
    36.0 / 28.0,
    8.0 / 28.0,
];

fn fixture(dir: &Path) -> std::path::PathBuf {
    let vcf = dir.join("sfs.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    vcf
}

fn request(folded: bool) -> StatsRequest {
    StatsRequest {
        kinds: Vec::new(),
        gt_groups: vec![Some(0); 4],
        window_size: 1000,
        step: 1000,
        levels: Levels { window: false, contig: false, genome: true },
        skip_empty: true,
        contig_lengths: std::collections::HashMap::new(),
        pixy: None,
        regions: Vec::new(),
        targets: Vec::new(),
        targets_exclude: Vec::new(),
        keep_samples: None,
        filter: None,
        sfs: Some(SfsRequest { projections: vec![6], folded }),
        fstats: None,
        ld: None,
    }
}

#[test]
fn projected_sfs_matches_the_hand_computed_spectrum() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = fixture(&dir);
    let budget = ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) };

    let report = run_stats(Some(&vcf), request(false), &budget).expect("invariant: run");
    let sfs = report.sfs.expect("invariant: SFS requested");
    assert_eq!(sfs.projections, vec![6]);
    assert!(!sfs.folded);
    assert_eq!(sfs.n_sites, vec![11]);
    assert_eq!(sfs.n_skipped_m_lt_n, vec![0]);
    assert_eq!(sfs.spectra.len(), 1);
    for (got, want) in sfs.spectra[0].iter().zip(&EXPECTED) {
        assert!((got - want).abs() <= 1e-12 * want.max(1.0), "got {got} want {want}");
    }

    // The dadi/fastsimcoal2 files: header, values and mask.
    let prefix = dir.join("out");
    let files = reed::sfs::write_projected_sfs(
        &prefix,
        &["pop1".to_string()],
        &sfs,
        false,
    )
    .expect("invariant: write");
    assert_eq!(files.len(), 2);
    let dadi = std::fs::read_to_string(format!("{}_sfs/dadi/pop1-6.sfs", prefix.display()))
        .expect("invariant: dadi file");
    let lines: Vec<&str> = dadi.lines().collect();
    assert_eq!(lines[0], "7 unfolded");
    assert_eq!(lines[2], "1 0 0 0 0 0 1");
    assert_eq!(lines[1].split(' ').count(), 7);
    let fsc2 =
        std::fs::read_to_string(format!("{}_sfs/fastsimcoal2/pop1_DAFpop0.obs", prefix.display()))
            .expect("invariant: fsc2 file");
    assert!(fsc2.starts_with("1 observation\nd0_0\td0_1\td0_2\td0_3\td0_4\td0_5\td0_6\n"));

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn folded_spectrum_pads_the_upper_half_and_masks_it() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-fold-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = fixture(&dir);
    let budget = ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) };

    let report = run_stats(Some(&vcf), request(true), &budget).expect("invariant: run");
    let sfs = report.sfs.expect("invariant: SFS requested");
    assert!(sfs.folded);
    // Folded: ξ'_k = ξ_k + ξ_{n−k} for 0 < k < n/2; the middle bin (n even)
    // and the k = 0/n pair are single contributions.
    let want = [
        EXPECTED[0] + EXPECTED[6],
        EXPECTED[1] + EXPECTED[5],
        EXPECTED[2] + EXPECTED[4],
        EXPECTED[3],
    ];
    assert_eq!(sfs.spectra[0].len(), 4);
    for (got, w) in sfs.spectra[0].iter().zip(&want) {
        assert!((got - w).abs() <= 1e-12 * w.max(1.0), "got {got} want {w}");
    }

    let prefix = dir.join("out");
    let files =
        reed::sfs::write_projected_sfs(&prefix, &["pop1".to_string()], &sfs, false)
            .expect("invariant: write");
    assert!(files.iter().any(|p| p.to_string_lossy().ends_with("pop1_MAFpop0.obs")));
    let dadi = std::fs::read_to_string(format!("{}_sfs/dadi/pop1-6.sfs", prefix.display()))
        .expect("invariant: dadi file");
    let lines: Vec<&str> = dadi.lines().collect();
    assert_eq!(lines[0], "7 folded");
    // easySFS keeps the n+1 length for folded spectra, zeroing the top.
    assert_eq!(lines[1].split(' ').count(), 7);
    assert!(lines[1].ends_with("0 0 0"), "{}", lines[1]);
    assert_eq!(lines[2], "1 0 0 0 1 1 1");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn spectra_are_identical_across_thread_counts() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-det-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = fixture(&dir);
    let mut digests = Vec::new();
    for threads in [1_usize, 3, 7] {
        let budget = ThreadBudget { threads, mem: MemBudget::from_bytes(u64::MAX) };
        let report = run_stats(Some(&vcf), request(false), &budget).expect("invariant: run");
        let sfs = report.sfs.expect("invariant: SFS requested");
        digests.push(sfs.spectra[0].clone());
    }
    assert!(digests.windows(2).all(|w| w[0] == w[1]), "{digests:?}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn missing_chromosomes_are_skipped_and_counted() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-skip-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("half.vcf");
    // One site, one sample called of two: m = 2 < n = 4 → skipped.
    std::fs::write(
        &vcf,
        VCF.replace(
            "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1",
            "chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t./.\t./.\t0/1",
        )
        .replace("\nchr1\t", "\nchr1\t"), // no-op; keeps the literal obvious
    )
    .expect("invariant: fixture write");
    let budget = ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) };
    let report = run_stats(Some(&vcf), request(false), &budget).expect("invariant: run");
    let sfs = report.sfs.expect("invariant: SFS requested");
    assert_eq!(sfs.n_skipped_m_lt_n, vec![1]);
    // The other ten sites (m = 8) still project.
    assert_eq!(sfs.n_sites, vec![10]);

    let _ = std::fs::remove_dir_all(&dir);
}
