//! LD acceptance (docs/rounds/v5/05-algorithms.md §3.4): the shared
//! fixture's per-pair r² and the plink2 `--indep-pairwise` prune sets.
//! `gold/popgen/run_plink2_ld.sh` cross-checks both against plink2; the
//! values here are plink2's (printed at six significant digits).

use std::fmt::Write as _;

use reed::budget::{MemBudget, ThreadBudget};
use reed::{
    InvariantMode, LdRequest, Levels, Populations, StatsRequest, WindowLayoutOptions, run_stats,
};

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

/// Two populations of three samples each (s0..s2 -> A, s3..s5 -> B).
fn pops_two_file(dir: &std::path::Path) -> std::path::PathBuf {
    let path = dir.join("pops.tsv");
    std::fs::write(&path, "s0\tA\ns1\tA\ns2\tA\ns3\tB\ns4\tB\ns5\tB\n").expect("invariant: pops");
    path
}

fn request(ld: LdRequest) -> StatsRequest {
    StatsRequest {
        kinds: Vec::new(),
        populations: Populations::All,
        grid: reed::WindowGrid::Bp {
            size: 1000,
            step: 1000,
        },
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
        site: None,
        force_samples: false,
        sfs: None,
        fstats: None,
        ld: Some(ld),
    }
}

#[test]
fn pairwise_r2_matches_plink2() {
    let dir = std::env::temp_dir().join(format!("reed-ld-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };

    let report = run_stats(
        Some(&vcf),
        request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
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

    let _ = std::fs::remove_dir_all(&dir);
}

/// The `--indep-pairwise` site rule lives in `reed filter`
/// (docs/rounds/v5/05-algorithms.md §3.4.1); this is the unit-level check of
/// its end-to-end contract on the fixture (plink2's `prune.in` = sites
/// 1, 5, 6 for `10 1 0.5`).
#[test]
fn filter_indep_pairwise_keeps_the_plink2_set() {
    let dir = std::env::temp_dir().join(format!("reed-ldf-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let out = dir.join("pruned.vcf");
    let mut spec = reed::spec::FilterSpec::default();
    spec.ld_prune = Some(reed::spec::LdPruneSpec {
        window: 10,
        kb: false,
        step: 1,
        r2: 0.5,
    });
    let opts = reed::FilterOptions {
        threads: 2,
        output: Some(out.clone()),
        ..reed::FilterOptions::default()
    };
    let report = reed::run_filter(Some(&vcf), &spec, &opts).expect("invariant: filter runs");
    assert_eq!(report.ld_prune, Some((6, 3)), "6 survivors, 3 pruned");
    let text = std::fs::read_to_string(&out).expect("invariant: output");
    let positions: Vec<u64> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(1).expect("POS").parse().expect("pos"))
        .collect();
    assert_eq!(positions, vec![1, 5, 6], "plink2 prune.in");
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
        let budget = ThreadBudget {
            threads,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        let report = run_stats(
            Some(&vcf),
            request(LdRequest {
                decay_max_bp: Some(100_000),
                decay_bin_bp: 1000,
                pairs: true,
            }),
            &budget,
        )
        .expect("invariant: run");
        let ld = report.ld.expect("invariant: LD requested");
        values.push((ld.pairs.clone(), ld.decay.clone(), ld.group_decay.clone()));
    }
    assert!(values.windows(2).all(|w| w[0] == w[1]), "{values:?}");

    let _ = std::fs::remove_dir_all(&dir);
}

/// With `--pop`, the report carries one decay curve per population next to
/// the cohort's: each population counts a pair only when the pair is
/// defined over that population's samples, so its counts never exceed the
/// cohort's (docs/rounds/v6/02-models.md §2.2's group-major matrix).
#[test]
fn ld_decay_carries_per_population_curves() {
    let dir = std::env::temp_dir().join(format!("reed-ld-pop-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };

    let mut req = request(LdRequest {
        decay_max_bp: Some(100_000),
        decay_bin_bp: 1000,
        pairs: false,
    });
    // Two populations of three samples each.
    req.populations = Populations::File(pops_two_file(dir.as_path()));
    let report = run_stats(Some(&vcf), req, &budget).expect("invariant: run");
    let ld = report.ld.expect("invariant: LD requested");

    assert_eq!(ld.group_decay.len(), 2, "one curve per population");
    for (bin, cohort) in ld.decay.iter().enumerate() {
        for pop in &ld.group_decay {
            let n = pop[bin].2;
            assert!(
                n <= cohort.2,
                "bin {bin}: population count {n} exceeds the cohort's {}",
                cohort.2
            );
        }
    }
    // The cohort pools every sample, so its bin holds all 15 pairs.
    assert_eq!(ld.decay[0].2, 15);

    let _ = std::fs::remove_dir_all(&dir);
}

/// The pruning runs after every site rule (02-cli §6 step 6b): the pruner
/// sees only post-tail sites, so a tail-failing site neither pairs nor
/// occupies a window slot (plink2 loads its filters before
/// `--indep-pairwise`). Fixture: sites 1 and 3 are perfect complements
/// (r² = 1, tied major frequency 0.5 → the later one prunes); site 2 is
/// 60% missing (fails `--geno 0.5`) with an orthogonal called pattern, so
/// no pair with it reaches the threshold.
#[test]
fn pruner_sees_only_post_tail_sites() {
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\
         s0\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\ts9\n",
    );
    // Site 1: alt on s0..s4 (MAF 0.5).
    let site1 = "0/1\t0/1\t0/1\t0/1\t0/1\t0/0\t0/0\t0/0\t0/0\t0/0";
    // Site 2: called only on s0..s3, alternating (fails max_missing 0.5).
    let site2 = "1/1\t0/0\t1/1\t0/0\t./.\t./.\t./.\t./.\t./.\t./.";
    // Site 3: the complement of site 1 (r² = 1, same MAF).
    let site3 = "0/0\t0/0\t0/0\t0/0\t0/0\t0/1\t0/1\t0/1\t0/1\t0/1";
    for (pos, gts) in [(1, site1), (2, site2), (3, site3)] {
        let _ = writeln!(text, "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t{gts}");
    }
    let dir = std::env::temp_dir().join(format!("reed-ld-tail-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, text).expect("invariant: fixture write");
    let spec_with = |max_missing: bool| reed::spec::FilterSpec {
        site_after_mask: reed::spec::TailRules {
            max_missing: max_missing.then(|| reed::spec::MissingRate::parse("0.5").expect("valid")),
            ..reed::spec::TailRules::default()
        },
        ld_prune: Some(reed::spec::LdPruneSpec {
            window: 3,
            kb: false,
            step: 1,
            r2: 0.5,
        }),
        ..reed::spec::FilterSpec::default()
    };
    let out = |name: &str| dir.join(name);
    let run = |spec: &reed::spec::FilterSpec, out: &std::path::Path| {
        reed::run_filter(
            Some(&vcf),
            spec,
            &reed::FilterOptions {
                threads: 2,
                output: Some(out.to_path_buf()),
                ..reed::FilterOptions::default()
            },
        )
        .expect("invariant: filter runs")
    };
    let positions = |out: &std::path::Path| -> Vec<u64> {
        std::fs::read_to_string(out)
            .expect("invariant: output")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).expect("POS").parse().expect("pos"))
            .collect()
    };

    // Without the tail rule: all three sites reach the pruner; the (1,3)
    // pair prunes site 3 (frequency tie → the later one).
    let report = run(&spec_with(false), &out("plain.vcf"));
    assert_eq!(report.ld_prune, Some((3, 1)), "3 seen, site 3 pruned");
    assert_eq!(positions(&out("plain.vcf")), vec![1, 2]);

    // With `--geno 0.5`: site 2 fails the tail and never reaches the
    // pruner (the pre-Y1e pruner counted it: n_seen would stay 3).
    let report = run(&spec_with(true), &out("tail.vcf"));
    assert_eq!(report.ld_prune, Some((2, 1)), "only post-tail sites seen");
    assert_eq!(positions(&out("tail.vcf")), vec![1]);
    let _ = std::fs::remove_dir_all(&dir);
}

/// `--mind` and `--indep-pairwise` share one pass0 (docs/rounds/v5/05
/// §3.4.1): the sample rule resolves first, then the pruner works over the
/// tail survivors, then pass1 emits both effects.
#[test]
fn mind_and_indep_pairwise_combine() {
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\
         s0\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\ts9\n",
    );
    // s0 is missing in 9 of 10 records (missing 0.9 > 0.5 → dropped).
    // Sites 1-5 share one pattern, 6-10 an orthogonal one: within each
    // group every pair has r² = 1 and tied frequencies, so the later
    // member prunes; across groups r² = 0.
    for i in 1..=10_u64 {
        let s0 = if i <= 9 { "./." } else { "0/0" };
        let (a, c) = if i <= 5 {
            ("0/1", "0/0")
        } else {
            ("0/0", "0/1")
        };
        let gts = [s0, a, a, "0/0", "0/0", c, c, "0/0", "0/0", "0/0"];
        let _ = writeln!(
            text,
            "chr1\t{i}\t.\tA\tG\t50\tPASS\t.\tGT\t{}",
            gts.join("\t")
        );
    }
    let dir = std::env::temp_dir().join(format!("reed-ld-mind-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, text).expect("invariant: fixture write");
    let spec = reed::spec::FilterSpec {
        mind: Some(reed::spec::MissingRate::parse("0.5").expect("valid")),
        ld_prune: Some(reed::spec::LdPruneSpec {
            window: 10,
            kb: false,
            step: 1,
            r2: 0.5,
        }),
        ..reed::spec::FilterSpec::default()
    };
    let out = dir.join("pruned.vcf");
    let report = reed::run_filter(
        Some(&vcf),
        &spec,
        &reed::FilterOptions {
            threads: 2,
            output: Some(out.clone()),
            ..reed::FilterOptions::default()
        },
    )
    .expect("invariant: filter runs");
    assert_eq!(report.n_samples_out, 9, "s0 dropped by --mind");
    assert_eq!(report.ld_prune, Some((10, 8)), "one survivor per group");
    let body: Vec<u64> = std::fs::read_to_string(&out)
        .expect("invariant: output")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(1).expect("POS").parse().expect("pos"))
        .collect();
    assert_eq!(body, vec![1, 6], "the earliest site of each group");
    let _ = std::fs::remove_dir_all(&dir);
}
