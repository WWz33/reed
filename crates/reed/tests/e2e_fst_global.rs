//! Global (r-population) WC84 Fst acceptance on the shared fixture
//! (docs/rounds/v5/05 §3.2). The expected values are vcftools' per-site
//! WEIR_AND_COCKERHAM_FST column over three populations (6 samples,
//! 2 per population), which the parity script cross-checks on real data;
//! here they are pinned as literals.

use reed::budget::{MemBudget, ThreadBudget};
use reed::{InvariantMode, Levels, StatKind, StatsRequest, WindowLayoutOptions, run_stats};

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

/// vcftools `--weir-fst-pop` ×3 on the fixture, per site (pos → Fst).
const VCFTOOLS_PER_SITE: [(u64, f64); 5] = [
    (1, 0.75),
    (2, 0.076_923_1),
    (3, -0.090_909_1),
    (4, 0.384_615),
    (5, -0.090_909_1),
];

fn request(window: u64, levels: Levels) -> StatsRequest {
    StatsRequest {
        kinds: vec![StatKind::FstGlobal],
        gt_groups: vec![Some(0), Some(0), Some(1), Some(1), Some(2), Some(2)],
        window_size: window,
        step: window,
        levels,
        skip_empty: true,
        contig_lengths: std::collections::HashMap::new(),
        // Structured window rows are filled in pixy-layout mode.
        layout: Some(WindowLayoutOptions {
            invariant: InvariantMode::VariantOnly,
        }),
        regions: Vec::new(),
        targets: Vec::new(),
        targets_exclude: Vec::new(),
        keep_samples: None,
        filter: None,
        sfs: None,
        fstats: None,
        ld: None,
    }
}

#[test]
fn global_fst_matches_vcftools_per_site_and_genome_wide() {
    let dir = std::env::temp_dir().join(format!("reed-fstg-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("fstats.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };

    // Per-site windows (size 1): one value per site.
    let report = run_stats(
        Some(&vcf),
        request(
            1,
            Levels {
                window: true,
                contig: false,
                genome: false,
            },
        ),
        &budget,
    )
    .expect("invariant: run");
    assert_eq!(report.windows.len(), 5);
    for (row, (pos, want)) in report.windows.iter().zip(&VCFTOOLS_PER_SITE) {
        assert_eq!(row.start1, *pos);
        let value = row.accs[0]
            .value(&StatKind::FstGlobal, 1)
            .expect("invariant: window has sites");
        assert!(
            (value - want).abs() <= 1e-6 * want.abs().max(1.0),
            "pos {pos}: reed {value} vs vcftools {want}"
        );
    }

    // Genome-wide: the Σa/Σ(a+b+c) over all sites (vcftools' weighted Fst).
    let report = run_stats(
        Some(&vcf),
        request(
            100_000,
            Levels {
                window: false,
                contig: false,
                genome: true,
            },
        ),
        &budget,
    )
    .expect("invariant: run");
    let value = report
        .windows
        .first()
        .and_then(|row| row.accs[0].value(&StatKind::FstGlobal, 100_000))
        .expect("invariant: genome row");
    assert!(
        (value - 0.25).abs() <= 1e-12,
        "genome Fst {value} vs vcftools 0.25"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn hudson_has_no_global_form() {
    // The CLI rejects `--fst hudson --fst-global`; the accumulator itself is
    // WC84-only, so a two-group request still builds (r = 2 WC84 is the
    // pairwise estimator).
    let dir = std::env::temp_dir().join(format!("reed-fstg-hud-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("fstats.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };

    let mut two = request(
        1,
        Levels {
            window: true,
            contig: false,
            genome: false,
        },
    );
    two.gt_groups = vec![Some(0), Some(0), Some(1), Some(1), None, None];
    let report = run_stats(Some(&vcf), two, &budget).expect("invariant: run");
    // Two populations: the global value equals the pairwise WC84 value.
    let global = report.windows[0].accs[0]
        .value(&StatKind::FstGlobal, 1)
        .expect("invariant: window has sites");
    let pairwise = run_stats(
        Some(&vcf),
        StatsRequest {
            kinds: vec![StatKind::Fst {
                a: 0,
                b: 1,
                kind: reed::FstKind::WeirCockerham84,
            }],
            gt_groups: vec![Some(0), Some(0), Some(1), Some(1), None, None],
            window_size: 1,
            step: 1,
            levels: Levels {
                window: true,
                contig: false,
                genome: false,
            },
            ..request(
                1,
                Levels {
                    window: true,
                    contig: false,
                    genome: false,
                },
            )
        },
        &budget,
    )
    .expect("invariant: run");
    let pairwise_value = pairwise.windows[0].accs[0]
        .value(
            &StatKind::Fst {
                a: 0,
                b: 1,
                kind: reed::FstKind::WeirCockerham84,
            },
            1,
        )
        .expect("invariant: window has sites");
    assert!(
        (global - pairwise_value).abs() <= 1e-12,
        "two-population global {global} != pairwise {pairwise_value}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn global_fst_is_identical_across_thread_counts() {
    let dir = std::env::temp_dir().join(format!("reed-fstg-det-{}", std::process::id()));
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
                1,
                Levels {
                    window: true,
                    contig: false,
                    genome: false,
                },
            ),
            &budget,
        )
        .expect("invariant: run");
        values.push(
            report
                .windows
                .iter()
                .map(|row| row.accs[0].value(&StatKind::FstGlobal, 1).expect("sites"))
                .collect::<Vec<_>>(),
        );
    }
    assert!(values.windows(2).all(|w| w[0] == w[1]), "{values:?}");

    let _ = std::fs::remove_dir_all(&dir);
}
