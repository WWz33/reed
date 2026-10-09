//! LD acceptance (docs/rounds/v5/05-algorithms.md §3.4): the shared
//! fixture's per-pair r² and the plink2 `--indep-pairwise` prune sets.
//! `gold/popgen/run_plink2_ld.sh` cross-checks both against plink2; the
//! values here are plink2's (printed at six significant digits).

use reed::budget::{MemBudget, ThreadBudget};
use reed::{run_stats, InvariantMode, LdPruneRequest, LdRequest, LdWindow, Levels, PixyOptions, StatsRequest};

const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\ts4\ts5
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t1/1\t0/1\t0/0
chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t1/1\t0/1\t0/0
chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t0/1\t0/0\t0/1\t1/1
chr1\t4\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t0/1\t0/0\t0/1\t1/1
chr1\t5\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/1\t0/0\t0/1
chr1\t6\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t0/1\t0/0\t0/1\t0/0\t0/1
";

/// plink2 `--r2-unphased` on the fixture: (pos_a, pos_b) → r².
const PLINK2_R2: [(u64, u64, f64); 15] = [
    (1, 2, 1.0),
    (1, 3, 0.584_775),
    (1, 4, 0.584_775),
    (1, 5, 0.058_823_5),
    (1, 6, 0.086_505_2),
    (2, 3, 0.584_775),
    (2, 4, 0.584_775),
    (2, 5, 0.058_823_5),
    (2, 6, 0.086_505_2),
    (3, 4, 1.0),
    (3, 5, 0.058_823_5),
    (3, 6, 0.169_55),
    (4, 5, 0.058_823_5),
    (4, 6, 0.169_55),
    (5, 6, 0.058_823_5),
];

/// plink2 `--indep-pairwise 10 1 0.5` on the fixture.
const PLINK2_PRUNE_IN: [u64; 3] = [1, 5, 6];
const PLINK2_PRUNE_OUT: [u64; 3] = [2, 3, 4];

fn request(ld: LdRequest) -> StatsRequest {
    StatsRequest {
        kinds: Vec::new(),
        gt_groups: vec![Some(0); 6],
        window_size: 1000,
        step: 1000,
        levels: Levels { window: false, contig: false, genome: true },
        skip_empty: true,
        contig_lengths: std::collections::HashMap::new(),
        pixy: Some(PixyOptions { invariant: InvariantMode::VariantOnly }),
        regions: Vec::new(),
        targets: Vec::new(),
        targets_exclude: Vec::new(),
        keep_samples: None,
        filter: None,
        sfs: None,
        fstats: None,
        ld: Some(ld),
    }
}

#[test]
fn pairwise_r2_and_prune_match_plink2() {
    let dir = std::env::temp_dir().join(format!("reed-ld-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let budget = ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) };

    let report = run_stats(
        Some(&vcf),
        request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
            prune: Some(LdPruneRequest {
                window: LdWindow::Variants(10),
                step: 1,
                r2: 0.5,
            }),
            pairs: true,
        }),
        &budget,
    )
    .expect("invariant: run");
    let ld = report.ld.expect("invariant: LD requested");

    assert_eq!(ld.pairs.len(), PLINK2_R2.len(), "pairs: {:?}", ld.pairs);
    for (contig, a, b, r2) in &ld.pairs {
        assert_eq!(contig, "chr1");
        let want = PLINK2_R2
            .iter()
            .find(|(wa, wb, _)| wa == a && wb == b)
            .unwrap_or_else(|| panic!("unexpected pair ({a}, {b})"));
        assert!(
            (r2 - want.2).abs() <= 2e-6,
            "pair ({a}, {b}): reed {r2} vs plink2 {}",
            want.2
        );
    }

    // Decay: 15 pairs inside [0, 1000) bp (the fixture's sites are ≤ 6 bp
    // apart), so the first bin holds all of them.
    assert_eq!(ld.decay[0].2, 15);
    assert_eq!(ld.decay[0].0, 0);
    assert_eq!(ld.decay[0].1, 1000);
    assert!(ld.decay[0].3.is_some());

    let mut in_pos: Vec<u64> = ld.prune_in.iter().map(|(_, p)| *p).collect();
    let mut out_pos: Vec<u64> = ld.prune_out.iter().map(|(_, p)| *p).collect();
    in_pos.sort_unstable();
    out_pos.sort_unstable();
    assert_eq!(in_pos, PLINK2_PRUNE_IN, "prune.in");
    assert_eq!(out_pos, PLINK2_PRUNE_OUT, "prune.out");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn ld_is_identical_across_thread_counts() {
    let dir = std::env::temp_dir().join(format!("reed-ld-det-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");

    let mut values = Vec::new();
    for threads in [1_usize, 3, 7] {
        let budget = ThreadBudget { threads, mem: MemBudget::from_bytes(u64::MAX) };
        let report = run_stats(
            Some(&vcf),
            request(LdRequest {
                decay_max_bp: Some(100_000),
                decay_bin_bp: 1000,
                prune: None,
                pairs: true,
            }),
            &budget,
        )
        .expect("invariant: run");
        let ld = report.ld.expect("invariant: LD requested");
        values.push(ld.pairs.clone());
    }
    assert!(values.windows(2).all(|w| w[0] == w[1]), "{values:?}");

    let _ = std::fs::remove_dir_all(&dir);
}
