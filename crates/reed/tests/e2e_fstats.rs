//! Patterson f-statistics acceptance (docs/rounds/v5/05-algorithms.md §3.3).
//! The expected values are scikit-allel `admixture.py`'s on the shared
//! 3-population fixture; the parity script cross-checks them on real data
//! (`gold/popgen/run_allel_fstats.sh`).

use reed::budget::{MemBudget, ThreadBudget};
use reed::{
    FStatTest, FStatsRequest, InvariantMode, Levels, StatsRequest, WindowLayoutOptions, run_stats,
};

const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\ts4\ts5
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/1\t1/1\t1/1
chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/1\t1/1\t0/1
chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t0/1\t0/1\t0/0\t0/1
chr1\t4\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t0/1\t0/0\t0/0\t0/1\t0/0
chr1\t5\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\t0/1\t0/1
";

/// scikit-allel `patterson_f2/f3/patterson_d` on the fixture (sums).
const ALLEL_F2_AB: f64 = 0.375;
const ALLEL_F3_C_AB: f64 = 0.291_666_666_666_666_7;
const ALLEL_F4_AB_CC: f64 = 0.0;
const ALLEL_D_AB_AC: f64 = 0.509_803_921_568_627_4;

fn request(tests: Vec<FStatTest>, block_bp: u64) -> StatsRequest {
    StatsRequest {
        kinds: Vec::new(),
        gt_groups: vec![Some(0), Some(0), Some(1), Some(1), Some(2), Some(2)],
        window_size: 1000,
        step: 1000,
        levels: Levels {
            window: false,
            contig: false,
            genome: true,
        },
        skip_empty: true,
        contig_lengths: std::collections::HashMap::new(),
        layout: Some(WindowLayoutOptions {
            invariant: InvariantMode::VariantOnly,
        }),
        regions: Vec::new(),
        targets: Vec::new(),
        targets_exclude: Vec::new(),
        keep_samples: None,
        filter: None,
        sfs: None,
        fstats: Some(FStatsRequest { tests, block_bp }),
        ld: None,
    }
}

#[test]
fn fstats_match_scikit_allel() {
    let dir = std::env::temp_dir().join(format!("reed-fstats-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("fstats.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };

    let tests = vec![
        FStatTest::F2 { a: 0, b: 1 },
        FStatTest::F3 { c: 2, a: 0, b: 1 },
        FStatTest::F4 {
            a: 0,
            b: 1,
            c: 2,
            d: 2,
        },
        FStatTest::Dstat {
            a: 0,
            b: 1,
            c: 0,
            d: 2,
        },
    ];
    let report = run_stats(Some(&vcf), request(tests, 1), &budget).expect("invariant: run");
    let fstats = report.fstats.expect("invariant: f-stats requested");
    assert_eq!(fstats.rows.len(), 4);

    let want = [ALLEL_F2_AB, ALLEL_F3_C_AB, ALLEL_F4_AB_CC, ALLEL_D_AB_AC];
    for (row, w) in fstats.rows.iter().zip(&want) {
        assert!(
            (row.estimate - w).abs() <= 1e-12 * w.abs().max(1.0),
            "{:?}: reed {} vs allel {w}",
            row.test,
            row.estimate
        );
        assert_eq!(row.n_sites, 5);
        // 1 bp blocks: one block per site (5 sites at distinct positions).
        assert_eq!(row.n_blocks, 5);
        assert!(row.se.is_finite() && row.se >= 0.0, "se {}", row.se);
    }

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn fstats_are_identical_across_thread_counts() {
    let dir = std::env::temp_dir().join(format!("reed-fstats-det-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("fstats.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");

    let mut values = Vec::new();
    for threads in [1_usize, 3, 7] {
        let budget = ThreadBudget {
            threads,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        let report = run_stats(
            Some(&vcf),
            request(
                vec![
                    FStatTest::F2 { a: 0, b: 1 },
                    FStatTest::Dstat {
                        a: 0,
                        b: 1,
                        c: 0,
                        d: 2,
                    },
                ],
                1,
            ),
            &budget,
        )
        .expect("invariant: run");
        let rows = report.fstats.expect("invariant: f-stats").rows;
        values.push(
            rows.iter()
                .map(|r| (r.estimate, r.se, r.n_blocks))
                .collect::<Vec<_>>(),
        );
    }
    assert!(values.windows(2).all(|w| w[0] == w[1]), "{values:?}");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn f4_ratio_divides_the_two_f4_estimates() {
    let dir = std::env::temp_dir().join(format!("reed-f4ratio-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("fstats.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };

    // f4-ratio = Σf4(A, O; X, C) / Σf4(A, O; B, C) with A = popA, O = popB,
    // X = popC, B = popA, C = popC. The two f4 sums come from independent
    // F4 tests over the same data, so the ratio is checkable by construction.
    let ratio = run_stats(
        Some(&vcf),
        request(
            vec![FStatTest::F4Ratio {
                o: 1,
                a: 0,
                x: 2,
                b: 0,
                c: 2,
            }],
            1,
        ),
        &budget,
    )
    .expect("invariant: run");
    let rows = ratio.fstats.expect("invariant: f-stats").rows;
    assert_eq!(rows.len(), 1);

    let parts = run_stats(
        Some(&vcf),
        request(
            vec![
                FStatTest::F4 {
                    a: 0,
                    b: 1,
                    c: 2,
                    d: 2,
                },
                FStatTest::F4 {
                    a: 0,
                    b: 1,
                    c: 0,
                    d: 2,
                },
            ],
            1,
        ),
        &budget,
    )
    .expect("invariant: run");
    let sums = parts.fstats.expect("invariant: f-stats").rows;
    let want = sums[0].estimate / sums[1].estimate;
    assert!(
        (rows[0].estimate - want).abs() <= 1e-12 * want.abs().max(1.0),
        "ratio {} vs {want}",
        rows[0].estimate
    );

    let _ = std::fs::remove_dir_all(&dir);
}
