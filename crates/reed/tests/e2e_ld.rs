//! LD acceptance (docs/rounds/v5/05-algorithms.md §3.4): the shared
//! fixture's per-pair r² and the plink2 `--indep-pairwise` prune sets.
//! `gold/popgen/run_plink2_ld.sh` cross-checks both against plink2; the
//! values here are plink2's (printed at six significant digits).

use std::fmt::Write as _;

use reed::budget::{MemBudget, ThreadBudget};
use reed::{
    InvariantMode, LdRequest, Levels, PopgenRequest, Populations, WindowLayoutOptions, run_popgen,
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

fn request(ld: LdRequest) -> PopgenRequest {
    PopgenRequest {
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
        window_output: None,
        compact_window_accumulators: false,
        skip_empty: true,
        contig_lengths: std::collections::HashMap::new(),
        layout: Some(WindowLayoutOptions {
            invariant: InvariantMode::VariantOnly,
        }),
        regions: Vec::new(),
        select: reed::spec::region::SiteSelector::default(),
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

    let report = run_popgen(
        Some(&vcf),
        &request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
            pairs: true,
            pairs_path: None,
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

/// doc_v7 D7-4: `--indep-pairwise` must run on the post-filter data — the
/// static sample subset and the genotype mask both drop samples from every
/// r²/major-frequency computation, so pruning a `-s` subset (or a set the
/// mask empties) must equal pruning a VCF that holds only those samples.
const LD_SEL_VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30\t0/1:30\t0/0:5\t1/1:5
chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30\t0/1:30\t0/0:5\t1/1:5
chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30\t0/0:30\t0/1:5\t0/0:5
";

/// The same three sites with only s0 and s1.
const LD_SEL_SUBSET_VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1
chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1
chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0
";

fn prune_positions(
    vcf: &std::path::Path,
    spec: &reed::spec::FilterSpec,
    out: &std::path::Path,
) -> Vec<u64> {
    let opts = reed::FilterOptions {
        threads: 2,
        output: Some(out.to_path_buf()),
        ..reed::FilterOptions::default()
    };
    reed::run_filter(Some(vcf), spec, &opts).expect("invariant: filter runs");
    std::fs::read_to_string(out)
        .expect("invariant: output")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(1).expect("POS").parse().expect("pos"))
        .collect()
}

fn prune_spec() -> reed::spec::FilterSpec {
    let mut spec = reed::spec::FilterSpec::default();
    spec.ld_prune = Some(reed::spec::LdPruneSpec {
        window: 10,
        kb: false,
        step: 1,
        r2: 0.5,
    });
    spec
}

#[test]
fn indep_pairwise_respects_the_static_sample_subset() {
    let dir = std::env::temp_dir().join(format!("reed-ld-subset-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let full = dir.join("full.vcf");
    std::fs::write(&full, LD_SEL_VCF).expect("invariant: write");
    let subset_vcf = dir.join("subset.vcf");
    std::fs::write(&subset_vcf, LD_SEL_SUBSET_VCF).expect("invariant: write");

    // Whole cohort: sites 1 and 2 are identical over all four samples.
    let all = prune_positions(&full, &prune_spec(), &dir.join("all.out.vcf"));

    // `-s s0,s1`: over the two kept samples sites 1 and 2 carry zero
    // variance, so nothing prunes.
    let mut spec = prune_spec();
    spec.sample_select = Some(reed::spec::SampleSelect {
        keep: Some(reed::spec::SampleList::Names(vec![
            "s0".into(),
            "s1".into(),
        ])),
        remove: None,
    });
    let kept_subset = prune_positions(&full, &spec, &dir.join("subset.out.vcf"));

    // The same result as pruning a VCF that holds only s0 and s1.
    let reference = prune_positions(&subset_vcf, &prune_spec(), &dir.join("ref.out.vcf"));

    assert_eq!(
        kept_subset, reference,
        "-s subset must prune like the subset VCF"
    );
    assert_ne!(all, reference, "the full cohort must prune differently");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn indep_pairwise_respects_the_genotype_mask() {
    let dir = std::env::temp_dir().join(format!("reed-ld-mask-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let full = dir.join("full.vcf");
    std::fs::write(&full, LD_SEL_VCF).expect("invariant: write");
    let subset_vcf = dir.join("subset.vcf");
    std::fs::write(&subset_vcf, LD_SEL_SUBSET_VCF).expect("invariant: write");

    // `--min-dp 20` masks s2 and s3 (DP=5), leaving s0/s1 — the subset.
    let mut spec = prune_spec();
    spec.genotype = reed::spec::GenotypeRules {
        mask_if: vec![reed::spec::GenotypeMaskRule {
            field: "DP".into(),
            when: "< 20".into(),
            on_missing: reed::spec::OnMissing::Keep,
            depth: true,
        }],
        ..reed::spec::GenotypeRules::default()
    };
    let masked = prune_positions(&full, &spec, &dir.join("masked.out.vcf"));
    let reference = prune_positions(&subset_vcf, &prune_spec(), &dir.join("ref.out.vcf"));
    assert_eq!(masked, reference, "the mask must drop s2/s3 from the prune");
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
        let report = run_popgen(
            Some(&vcf),
            &request(LdRequest {
                decay_max_bp: Some(100_000),
                decay_bin_bp: 1000,
                pairs: true,
                pairs_path: None,
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
        pairs_path: None,
    });
    // Two populations of three samples each.
    req.populations = Populations::File(pops_two_file(dir.as_path()));
    let report = run_popgen(Some(&vcf), &req, &budget).expect("invariant: run");
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

/// Truncates every record to the first `keep` sample columns.
fn drop_columns(text: &str, keep: usize) -> String {
    let mut out = String::new();
    for l in text.lines() {
        if l.starts_with("##") {
            out.push_str(l);
        } else {
            let mut f: Vec<&str> = l.split('\t').collect();
            f.truncate(9 + keep);
            out.push_str(&f.join("\t"));
        }
        out.push('\n');
    }
    out
}

/// doc_v7 D7-5: the LD cohort must not include samples dropped by
/// `-s/-S/--keep/--remove`. With `Populations::All` every sample is
/// ungrouped, so it lands in the trailing slice; a dropped sample must get
/// no plane bit, or it pollutes the cohort sums. Pruning the six-sample
/// fixture down to s0..s3 must equal the four-sample VCF.
#[test]
fn ld_cohort_excludes_removed_samples() {
    let dir = std::env::temp_dir().join(format!("reed-ld-keep-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let full = dir.join("full.vcf");
    std::fs::write(&full, VCF).expect("invariant: write");
    let four = dir.join("four.vcf");
    std::fs::write(&four, drop_columns(VCF, 4)).expect("invariant: write");

    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let ld = || LdRequest {
        decay_max_bp: Some(100_000),
        decay_bin_bp: 1000,
        pairs: true,
        pairs_path: None,
    };

    let mut kept = request(ld());
    kept.keep_samples = Some(vec![true, true, true, true, false, false]);
    let a = run_popgen(Some(&full), &kept, &budget)
        .expect("invariant: run")
        .ld
        .expect("invariant: LD requested");

    let b = run_popgen(Some(&four), &request(ld()), &budget)
        .expect("invariant: run")
        .ld
        .expect("invariant: LD requested");

    assert_eq!(
        a.decay, b.decay,
        "removed samples must not enter the cohort decay"
    );
    assert_eq!(a.pairs, b.pairs, "removed samples must not enter the pairs");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v7 D7-19: the streaming prune caps concurrent per-contig threads.
/// With `threads` below the contig count the oldest job is joined early, so
/// the multi-contig result must still be thread-count independent.
#[test]
fn prune_multi_contig_is_thread_count_independent() {
    let dir = std::env::temp_dir().join(format!("reed-ld-multi-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=c1,length=1000>\n\
         ##contig=<ID=c2,length=1000>\n\
         ##contig=<ID=c3,length=1000>\n\
         ##contig=<ID=c4,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\ts4\ts5\n",
    );
    for c in ["c1", "c2", "c3", "c4"] {
        let _ = write!(
            text,
            "{c}\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t1/1\t0/1\t0/0\n\
             {c}\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t1/1\t0/1\t0/0\n\
             {c}\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t0/1\t0/0\t0/1\t1/1\n"
        );
    }
    let vcf = dir.join("multi.vcf");
    std::fs::write(&vcf, text).expect("invariant: write");
    let run = |threads: usize, out: &std::path::Path| -> Vec<(String, u64)> {
        let opts = reed::FilterOptions {
            threads,
            output: Some(out.to_path_buf()),
            ..reed::FilterOptions::default()
        };
        let mut spec = prune_spec();
        spec.ld_prune = Some(reed::spec::LdPruneSpec {
            window: 10,
            kb: false,
            step: 1,
            r2: 0.5,
        });
        reed::run_filter(Some(&vcf), &spec, &opts).expect("invariant: filter runs");
        std::fs::read_to_string(out)
            .expect("invariant: output")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let mut it = l.split('\t');
                (
                    it.next().expect("CHROM").to_string(),
                    it.next().expect("POS").parse().expect("pos"),
                )
            })
            .collect()
    };
    let reference = run(1, &dir.join("t1.vcf"));
    assert_ne!(reference.len(), 0, "the one-thread prune kept nothing");
    for threads in [2_usize, 4, 8] {
        let out = dir.join(format!("t{threads}.vcf"));
        assert_eq!(
            run(threads, &out),
            reference,
            "threads={threads} must match the one-thread prune"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v8 D8-2: the sink's in-flight decay jobs are capped at
/// `budget.auxiliary()`, so a many-contig run joins the oldest job mid-run
/// instead of keeping every contig's matrix resident. The cap must not
/// change the output: joins stay in contig order, so the merged decay and
/// pair tables are identical across thread counts (the single-thread
/// budget has auxiliary() == 1, exercising the join path on every contig).
#[test]
fn ld_sink_caps_in_flight_decay_jobs() {
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    );
    for contig in 0..10 {
        let _ = writeln!(text, "##contig=<ID=chr{contig},length=10000>");
    }
    text.push_str(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\ts4\ts5\n",
    );
    for contig in 0..10 {
        for pos in 1..=3 {
            let _ = writeln!(
                text,
                "chr{contig}\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t\
                 0/1\t0/1\t0/0\t1/1\t0/1\t0/0"
            );
        }
    }
    let dir = std::env::temp_dir().join(format!("reed-ld-cap-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, text).expect("invariant: fixture write");

    let mut values = Vec::new();
    for threads in [1_usize, 8] {
        let budget = ThreadBudget {
            threads,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        let report = run_popgen(
            Some(&vcf),
            &request(LdRequest {
                decay_max_bp: Some(100_000),
                decay_bin_bp: 1000,
                pairs: true,
                pairs_path: None,
            }),
            &budget,
        )
        .expect("invariant: run");
        let ld = report.ld.expect("invariant: LD requested");
        assert_eq!(ld.pairs.len(), 10 * 3, "every contig's pairs survive");
        values.push((ld.pairs.clone(), ld.decay.clone(), ld.group_decay.clone()));
    }
    assert_eq!(values[0], values[1], "the cap must not change the output");

    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v23 V23-14 (decision C): `--ld-pairs` streams each contig's rows to
/// the final `PREFIX.ld` as the decay jobs join, instead of accumulating
/// every pair in memory. The streamed file must be byte-identical to what
/// `write_ld` renders from an accumulated report, and the report's
/// `pairs` must be empty (the rows already left).
#[test]
fn streamed_pairs_file_is_byte_identical_to_the_accumulated_report() {
    let dir = std::env::temp_dir().join(format!("reed-ld-stream-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, VCF).expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };

    // Accumulated (the pre-v23 library path): rows come back in the report.
    let accumulated = run_popgen(
        Some(&vcf),
        &request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
            pairs: true,
            pairs_path: None,
        }),
        &budget,
    )
    .expect("invariant: run")
    .ld
    .expect("invariant: LD requested");
    assert!(!accumulated.pairs.is_empty(), "control run must carry rows");

    // Streamed: the sink writes the file during the run.
    let prefix = dir.join("out");
    let streamed = run_popgen(
        Some(&vcf),
        &request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
            pairs: true,
            pairs_path: Some(prefix.with_extension("ld")),
        }),
        &budget,
    )
    .expect("invariant: run")
    .ld
    .expect("invariant: LD requested");
    assert!(
        streamed.pairs.is_empty(),
        "streamed rows must not also be held in memory"
    );
    assert_eq!(streamed.decay, accumulated.decay, "decay is unaffected");

    // Byte identity: render the accumulated report through the writer, and
    // compare with what the sink streamed.
    let mut expected = String::from("#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n");
    for (contig, a, b, r2) in &accumulated.pairs {
        reed::ld::render_pair_row(&mut expected, contig, *a, *b, *r2);
    }
    let got = std::fs::read_to_string(prefix.with_extension("ld")).expect("invariant: .ld written");
    assert_eq!(
        got, expected,
        "streamed pairs differ from the report render"
    );

    // The decay writer must not recreate or truncate the streamed file.
    let decay_prefix = dir.join("decay_only");
    let mut report = accumulated.clone();
    report.pairs.clear();
    let files = reed::ld::write_ld(&decay_prefix, &report, &[], false).expect("invariant: write");
    assert!(
        files.iter().all(|p| !p.to_string_lossy().ends_with(".ld")),
        "-.ld must not be written when the rows were streamed: {files:?}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// A run with no pair at all (a single site) publishes the header-only
/// pair table: the request determines the output, not whether a pair arrived.
#[test]
fn streamed_pairs_file_is_header_only_without_pairs() {
    let dir = std::env::temp_dir().join(format!("reed-ld-nopair-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let text = "##fileformat=VCFv4.2\n\
        ##contig=<ID=chr1,length=100000>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
        chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n";
    let vcf = dir.join("one.vcf");
    std::fs::write(&vcf, text).expect("invariant: write");
    let prefix = dir.join("out");
    let path = prefix.with_extension("ld");
    let budget = ThreadBudget {
        threads: 1,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    run_popgen(
        Some(&vcf),
        &request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
            pairs: true,
            pairs_path: Some(path.clone()),
        }),
        &budget,
    )
    .expect("invariant: run");
    assert_eq!(
        std::fs::read_to_string(&path).expect("invariant: pair table"),
        String::from_utf8_lossy(reed::ld::PAIRS_HEADER),
        "no pair means a valid header-only .ld, never a stale previous table"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// The streaming writer's header is opened at run start; a pairless leading
/// contig contributes no rows, and a later contig with pairs appends them in
/// contig order.
#[test]
fn streamed_pairs_skip_a_pairless_leading_contig() {
    let dir = std::env::temp_dir().join(format!("reed-ld-lead-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let text = "##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0
chr2\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0
chr2\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1
";
    let vcf = dir.join("ld.vcf");
    std::fs::write(&vcf, text).expect("invariant: write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    // Accumulated reference.
    let accumulated = run_popgen(
        Some(&vcf),
        &request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
            pairs: true,
            pairs_path: None,
        }),
        &budget,
    )
    .expect("invariant: run")
    .ld
    .expect("invariant: LD requested");
    assert_eq!(
        accumulated.pairs.len(),
        1,
        "only chr2's one pair exists: {:?}",
        accumulated.pairs
    );

    let path = dir.join("out.ld");
    let streamed = run_popgen(
        Some(&vcf),
        &request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
            pairs: true,
            pairs_path: Some(path.clone()),
        }),
        &budget,
    )
    .expect("invariant: run")
    .ld
    .expect("invariant: LD requested");
    assert_eq!(streamed.pairs.len(), 0, "streamed rows stay out of memory");
    let mut expected = String::from("#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n");
    for (contig, a, b, r2) in &accumulated.pairs {
        reed::ld::render_pair_row(&mut expected, contig, *a, *b, *r2);
    }
    assert_eq!(
        std::fs::read_to_string(&path).expect("the later contig creates the file"),
        expected
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-02: a run that fails *after* the sink finished used to publish
/// its (truncated) streaming `.ld` over the previous run's good file — the
/// sink committed the `AtomicFile`, and only then did `job.finish` detect a
/// BCF record truncated at the end of the stream.
///
/// The fixture is a legal BGZF stream whose last record is incomplete, so
/// every check before `job.finish` passes. After the fix the previous `.ld`
/// must be byte-for-byte unchanged and no temp file may survive.
#[test]
fn failed_run_leaves_the_previous_ld_untouched() {
    let dir = std::env::temp_dir().join(format!("reed-ld-v2702-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");

    // Enough records to span several BGZF members, so the truncation lands
    // inside the last member's payload and is only detectable once the
    // stream ends (the record sits in the framing carry).
    let mut text = String::from(
        "##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
",
    );
    for pos in 1..=800 {
        let gts = if pos % 2 == 0 { "0/1\t0/0" } else { "0/0\t0/1" };
        let _ = writeln!(text, "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t{gts}");
    }
    let vcf = dir.join("in.vcf");
    std::fs::write(&vcf, text).expect("invariant: write vcf");

    // A BCF through the facade's own filter writer, then recompressed with
    // the last few payload bytes dropped: still one valid BGZF member, with
    // the final record incomplete.
    let good = dir.join("good.bcf");
    reed::run_filter(
        Some(&vcf),
        &reed::spec::FilterSpec::default(),
        &reed::FilterOptions {
            threads: 1,
            output: Some(good.clone()),
            write_index: reed::WriteIndex::None,
            compress_level: 6,
            strategy: None,
            command_line: String::new(),
            timing: None,
            report_json: None,
            force_samples: false,
            regions: Vec::new(),
            tmp_dir: None,
            output_compressed: Some(true),
            output_bcf: Some(6),
            mem_bytes: None,
            no_version: true,
            no_command: true,
        },
    )
    .expect("invariant: build the BCF fixture");

    let mut raw = Vec::new();
    {
        let bytes = std::fs::read(&good).expect("invariant: read bcf");
        let mut inflater = reed::bgzf::inflate::Inflater::new();
        let mut pos = 0;
        while pos + 18 <= bytes.len() {
            // `inflate_block` wants exactly one member: slice by BSIZE.
            let bsize = u16::from_le_bytes([bytes[pos + 16], bytes[pos + 17]]) as usize + 1;
            if pos + bsize > bytes.len() {
                break;
            }
            let block = &bytes[pos..pos + bsize];
            if reed::bgzf::block::is_eof_block(block) {
                break;
            }
            inflater
                .inflate_block(block, &mut raw)
                .expect("invariant: inflate one member");
            pos += bsize;
        }
    }
    assert!(raw.len() > 32, "the BCF payload must be non-trivial");

    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut truncated = Vec::new();
    deflater
        .deflate_block(&raw[..raw.len() - 8], &mut truncated)
        .expect("invariant: deflate");
    truncated.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);
    let cut = dir.join("cut.bcf");
    std::fs::write(&cut, &truncated).expect("invariant: write cut bcf");

    let prefix = dir.join("out");
    let pairs_path = prefix.with_extension("ld");
    let budget = ThreadBudget {
        threads: 1,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let ld = |pairs_path: Option<std::path::PathBuf>| LdRequest {
        decay_max_bp: Some(100_000),
        decay_bin_bp: 1000,
        pairs: true,
        pairs_path,
    };

    let accumulated = run_popgen(Some(&good), &request(ld(None)), &budget)
        .expect("invariant: clean run")
        .ld
        .expect("invariant: ld requested");
    assert!(
        !accumulated.pairs.is_empty(),
        "the fixture must yield one pair"
    );
    let mut expected = String::from_utf8_lossy(reed::ld::PAIRS_HEADER).into_owned();
    for (contig, a, b, r2) in &accumulated.pairs {
        reed::ld::render_pair_row(&mut expected, contig, *a, *b, *r2);
    }
    std::fs::write(&pairs_path, &expected).expect("invariant: seed the old .ld");

    let err = run_popgen(Some(&cut), &request(ld(Some(pairs_path.clone()))), &budget)
        .expect_err("invariant: the truncated BCF must fail");
    assert!(
        err.to_string().contains("inside a record"),
        "unexpected error: {err}"
    );

    let after = std::fs::read_to_string(&pairs_path).expect("invariant: .ld readable");
    assert_eq!(
        after, expected,
        "a failed run must leave the previous .ld byte-for-byte unchanged"
    );

    let leftovers: Vec<String> = std::fs::read_dir(&dir)
        .expect("invariant: read dir")
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.contains("reed-tmp"))
        .collect();
    assert!(leftovers.is_empty(), "no temp may survive: {leftovers:?}");

    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-07: LD's genome-wide matrix analysis has no window-clipping
/// dependency, so it accepts a contig declaration without a length.
#[test]
fn genome_wide_ld_does_not_require_contig_lengths() {
    let dir = std::env::temp_dir().join(format!("reed-ld-nolen-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("ld.vcf");
    std::fs::write(
        &vcf,
        VCF.replace("##contig=<ID=chr1,length=100000>", "##contig=<ID=chr1>"),
    )
    .expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let report = run_popgen(
        Some(&vcf),
        &request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
            pairs: true,
            pairs_path: None,
        }),
        &budget,
    )
    .expect("invariant: no-length LD run");
    assert_ne!(report.ld.expect("invariant: LD report").pairs.len(), 0);
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-05/V27-20: a successful zero-pair run publishes a header-only
/// pair table and the corresponding 0/NA decay bins. Reusing a prefix must
/// replace, not preserve, the previous run's pair rows.
#[test]
fn an_empty_ld_run_publishes_header_only_pairs_and_zero_decay() {
    let dir = std::env::temp_dir().join(format!("reed-ld-empty-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("one.vcf");
    let text = "##fileformat=VCFv4.2\n\
        ##contig=<ID=chr1,length=100000>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
        chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n";
    std::fs::write(&vcf, text).expect("invariant: write");
    let path = dir.join("out.ld");
    std::fs::write(&path, "stale pair rows\n").expect("invariant: stale file");
    let budget = ThreadBudget {
        threads: 1,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let report = run_popgen(
        Some(&vcf),
        &request(LdRequest {
            decay_max_bp: Some(2_000),
            decay_bin_bp: 1000,
            pairs: true,
            pairs_path: Some(path.clone()),
        }),
        &budget,
    )
    .expect("invariant: empty LD run");
    assert_eq!(
        std::fs::read_to_string(&path).expect("invariant: pair table"),
        String::from_utf8_lossy(reed::ld::PAIRS_HEADER)
    );
    assert_eq!(report.pairs_path, Some(path));
    let ld = report.ld.expect("invariant: LD report");
    assert_eq!(ld.decay.len(), 2);
    assert!(
        ld.decay
            .iter()
            .all(|&(_, _, n, mean)| n == 0 && mean.is_none())
    );
    let files = reed::ld::write_ld(&dir.join("out"), &ld, &["all".into()], false)
        .expect("invariant: decay writer");
    assert!(files.iter().any(|p| p.ends_with("out.ld.decay")));
    assert!(
        files.iter().all(|p| !p.ends_with("out.ld")),
        "streamed pairs are not rewritten"
    );
    let decay = std::fs::read_to_string(dir.join("out.ld.decay")).expect("invariant: decay");
    assert!(decay.contains("all\t0\t1000\t0\tNA"), "{decay}");
    assert!(decay.contains("all\t1000\t2000\t0\tNA"), "{decay}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-16/V27-17: LD work shares the run's thread budget and is
/// cancelled when the library call fails. This is a library-level check:
/// a CLI's immediate process exit would hide a leaked worker pool.
#[test]
fn ld_error_cancels_workers_and_leaves_no_auxiliary_threads() {
    let dir = std::env::temp_dir().join(format!("reed-ld-abort-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("bad.vcf");
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=30000>\n\
         ##contig=<ID=chr2,length=100>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n",
    );
    for pos in 1..=20_000 {
        writeln!(
            text,
            "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1"
        )
        .expect("invariant: write to String");
    }
    text.push_str("chr2\t999999\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n");
    std::fs::write(&vcf, text).expect("invariant: fixture");
    let threads_before = std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find_map(|line| line.strip_prefix("Threads:")?.trim().parse::<usize>().ok())
        })
        .unwrap_or(1);
    let err = run_popgen(
        Some(&vcf),
        &request(LdRequest {
            decay_max_bp: Some(100_000),
            decay_bin_bp: 1000,
            pairs: true,
            pairs_path: None,
        }),
        &ThreadBudget {
            threads: 8,
            mem: MemBudget::from_bytes(u64::MAX),
        },
    )
    .expect_err("the out-of-range record must fail");
    assert!(err.to_string().contains("exceeds contig"), "{err}");
    // Dropping the sink triggers abort and drops the one LD pool. Allow a
    // small settling margin for rayon workers to observe the signal; the
    // assertion is about ownership, not a wall-clock performance promise.
    std::thread::sleep(std::time::Duration::from_millis(100));
    let threads_after = std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find_map(|line| line.strip_prefix("Threads:")?.trim().parse::<usize>().ok())
        })
        .unwrap_or(1);
    assert!(
        threads_after <= threads_before + 1,
        "LD workers survived the failed run: before={threads_before}, after={threads_after}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
