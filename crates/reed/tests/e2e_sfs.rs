//! Projected-SFS acceptance on the shared fixture (docs/rounds/v5/05
//! §3.1): 4 diploid samples, 11 sites with ALT counts 1,1,2,2,3,3,4,4,5,6,7
//! of 8 chromosomes, projected to 6. The expectations are the exact
//! hypergeometric fractions (1/28 granularity), computed by hand; the
//! fastsimcoal2/dadi file layouts are checked against easySFS's.

use std::path::Path;

use reed::budget::{MemBudget, ThreadBudget};
use reed::{Levels, PopgenRequest, Populations, SfsRequest, run_popgen};

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

fn request(folded: bool) -> PopgenRequest {
    PopgenRequest {
        kinds: Vec::new(),
        populations: Populations::All,
        force_samples: false,
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
        layout: None,
        regions: Vec::new(),
        select: reed::spec::region::SiteSelector::default(),
        keep_samples: None,
        site: None,
        sfs: Some(SfsRequest {
            projections: vec![6],
            folded,
        }),
        fstats: None,
        ld: None,
    }
}

#[test]
fn projected_sfs_matches_the_hand_computed_spectrum() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = fixture(&dir);
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };

    let report = run_popgen(Some(&vcf), &request(false), &budget).expect("invariant: run");
    let sfs = report.sfs.expect("invariant: SFS requested");
    assert_eq!(sfs.projections, vec![6]);
    assert!(!sfs.folded);
    assert_eq!(sfs.n_sites, vec![11]);
    assert_eq!(sfs.n_skipped_m_lt_n, vec![0]);
    assert_eq!(sfs.spectra.len(), 1);
    for (got, want) in sfs.spectra[0].iter().zip(&EXPECTED) {
        assert!(
            (got - want).abs() <= 1e-12 * want.max(1.0),
            "got {got} want {want}"
        );
    }

    // The dadi/fastsimcoal2 files: header, values and mask.
    let prefix = dir.join("out");
    let files = reed::sfs::write_projected_sfs(&prefix, &["pop1".to_string()], &sfs, false)
        .expect("invariant: write");
    assert_eq!(files.len(), 2);
    let dadi = std::fs::read_to_string(format!("{}_sfs/dadi/pop1-6.sfs", prefix.display()))
        .expect("invariant: dadi file");
    let lines: Vec<&str> = dadi.lines().collect();
    assert_eq!(lines[0], "7 unfolded");
    assert_eq!(lines[2], "1 0 0 0 0 0 1");
    assert_eq!(lines[1].split(' ').count(), 7);
    let fsc2 = std::fs::read_to_string(format!(
        "{}_sfs/fastsimcoal2/pop1_DAFpop0.obs",
        prefix.display()
    ))
    .expect("invariant: fsc2 file");
    assert!(fsc2.starts_with("1 observation\nd0_0\td0_1\td0_2\td0_3\td0_4\td0_5\td0_6\n"));

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn folded_spectrum_pads_the_upper_half_and_masks_it() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-fold-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = fixture(&dir);
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };

    let report = run_popgen(Some(&vcf), &request(true), &budget).expect("invariant: run");
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
    let files = reed::sfs::write_projected_sfs(&prefix, &["pop1".to_string()], &sfs, false)
        .expect("invariant: write");
    assert!(
        files
            .iter()
            .any(|p| p.to_string_lossy().ends_with("pop1_MAFpop0.obs"))
    );
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
        let budget = ThreadBudget {
            threads,
            mem: MemBudget::from_bytes(u64::MAX),
        };
        let report = run_popgen(Some(&vcf), &request(false), &budget).expect("invariant: run");
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
        ),
    )
    .expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let report = run_popgen(Some(&vcf), &request(false), &budget).expect("invariant: run");
    let sfs = report.sfs.expect("invariant: SFS requested");
    assert_eq!(sfs.n_skipped_m_lt_n, vec![1]);
    // The other ten sites (m = 8) still project.
    assert_eq!(sfs.n_sites, vec![10]);

    let _ = std::fs::remove_dir_all(&dir);
}

/// The two-population joint (2D) spectrum on the same fixture split
/// `s0,s1 | s2,s3`, projected `4,4`. Expectations verified against
/// easySFS's `dadi_twoD_sfs_combinations` on an identical VCF (folded and
/// unfolded, dadi and fastsimcoal2 layouts byte-equal modulo the dadi
/// population-label header, which reed omits like the 1D writer).
#[test]
fn joint_spectrum_matches_easysfs() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-joint-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = fixture(&dir);
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let mut request = request(true);
    request.populations = Populations::File({
        let pops = dir.join("pops.tsv");
        std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: pops");
        pops
    });
    request.sfs = Some(SfsRequest {
        projections: vec![4, 4],
        folded: true,
    });
    let report = run_popgen(Some(&vcf), &request, &budget).expect("invariant: run");
    let sfs = report.sfs.expect("invariant: SFS requested");
    let joint = sfs.joint.clone().expect("invariant: joint spectra");
    assert_eq!(joint.pairs, vec![(0, 1)]);
    assert_eq!(joint.projections, vec![(4, 4)]);
    assert!(joint.folded);
    // easySFS's folded dadi matrix for the same 11-site fixture (row-major).
    let want = [
        0.0, 2.0, 2.0, 1.0, 0.5, //
        1.0, 1.0, 0.0, 0.5, 0.0, //
        0.0, 1.0, 0.0, 0.0, 0.0, //
        1.0, 0.5, 0.0, 0.0, 0.0, //
        0.5, 0.0, 0.0, 0.0, 0.0,
    ];
    assert_eq!(joint.spectra[0].len(), 25);
    for (got, w) in joint.spectra[0].iter().zip(&want) {
        assert!((got - w).abs() <= 1e-12, "got {got} want {w}");
    }

    // Files: dadi `<A>-<B>.sfs` and the transposed fastsimcoal2 layout.
    let prefix = dir.join("out");
    let files = reed::sfs::write_projected_sfs(
        &prefix,
        &["popA".to_string(), "popB".to_string()],
        &sfs,
        false,
    )
    .expect("invariant: write");
    assert!(files.iter().any(|p| p.ends_with("dadi/popA-popB.sfs")));
    let dadi = std::fs::read_to_string(format!("{}_sfs/dadi/popA-popB.sfs", prefix.display()))
        .expect("invariant: dadi file");
    let lines: Vec<&str> = dadi.lines().collect();
    assert_eq!(lines[0], "5 5 folded");
    assert_eq!(
        lines[2],
        "1 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 1 1 1 0 1 1 1 1"
    );
    let fsc = std::fs::read_to_string(format!(
        "{}_sfs/fastsimcoal2/out_jointMAFpop1_0.obs",
        prefix.display()
    ))
    .expect("invariant: fsc2 read");
    assert_eq!(
        fsc,
        "1 observation\n\td0_0\td0_1\td0_2\td0_3\td0_4\nd1_0\t0 1 0 1 0.5\nd1_1\t2 1 1 0.5 0\nd1_2\t2 0 0 0 0\nd1_3\t1 0.5 0 0 0\nd1_4\t0.5 0 0 0 0\n"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-08: SFS file names carry the population name, which comes
/// from the operator's populations file and used to be appended verbatim.
/// A group named `../../../../ESCAPED-<pid>` therefore wrote outside the
/// `-o` directory tree — in this layout all the way to the temp dir's
/// parent. The name is pid-unique so the assertion cannot be satisfied (or
/// broken) by a leftover from another run.
#[test]
fn sfs_group_names_cannot_escape_the_output_directory() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-v2708-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");

    let tag = format!("ESCAPED-{}", std::process::id());
    let name = format!("../../../../{tag}");
    // The escaped names a sanitizing writer must never produce. The 1-D
    // spectrum projects to 6 (the fixture's projection).
    let outside: Vec<std::path::PathBuf> = [
        std::env::temp_dir()
            .parent()
            .expect("invariant: temp dir has a parent")
            .join(format!("{tag}-6.sfs")),
        std::env::temp_dir()
            .parent()
            .expect("invariant: temp dir has a parent")
            .join(format!("{tag}_DAFpop0.obs")),
    ]
    .into_iter()
    .collect();
    for path in &outside {
        let _ = std::fs::remove_file(path);
    }

    let vcf = fixture(&dir);
    let budget = ThreadBudget {
        threads: 1,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let pops = dir.join("esc.tsv");
    std::fs::write(
        &pops,
        format!("s0\t{name}\ns1\t{name}\ns2\t{name}\ns3\t{name}\n"),
    )
    .expect("invariant: write pops");
    let mut req = request(false);
    req.populations = Populations::File(pops);
    let report = run_popgen(Some(&vcf), &req, &budget).expect("invariant: run");
    let sfs = report.sfs.expect("invariant: SFS requested");

    let prefix = dir.join("nested/out");
    std::fs::create_dir_all(dir.join("nested")).expect("invariant: nested");
    let files = reed::sfs::write_projected_sfs(&prefix, std::slice::from_ref(&name), &sfs, false)
        .expect("invariant: write");
    assert_eq!(files.len(), 2);

    let root = dir.join("nested");
    for path in &files {
        assert!(
            path.starts_with(&root),
            "{} escaped the -o tree ({})",
            path.display(),
            root.display()
        );
        assert!(
            !path.components().any(|c| c.as_os_str() == ".."),
            "{} still carries a parent component",
            path.display()
        );
    }
    for path in &outside {
        assert!(
            !path.exists(),
            "{} was written outside the tree",
            path.display()
        );
    }

    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-09: the 1-D `<name>-<n>.sfs` and the pairwise `<A>-<B>.sfs`
/// share `_sfs/dadi`. Group `A` projected to 4 and a second group named `4`
/// used to address the same path, so one spectrum silently overwrote the
/// other. The writer now refuses instead of picking a winner.
#[test]
fn sfs_one_dimensional_and_joint_names_cannot_collide() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-v2709-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = fixture(&dir);
    let budget = ThreadBudget {
        threads: 1,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let pops = dir.join("c.tsv");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\t4\ns3\t4\n").expect("invariant: write pops");
    let mut req = request(false);
    req.populations = Populations::File(pops);
    req.sfs = Some(SfsRequest {
        folded: false,
        projections: vec![4, 4],
    });
    let report = run_popgen(Some(&vcf), &req, &budget).expect("invariant: run");
    let sfs = report.sfs.expect("invariant: SFS requested");

    let prefix = dir.join("out");
    let err =
        reed::sfs::write_projected_sfs(&prefix, &["A".to_string(), "4".to_string()], &sfs, false)
            .expect_err("invariant: A-4 collides with the joint A/4 spectrum");
    let text = err.to_string();
    assert!(
        text.contains("A-4.sfs") && text.contains("would write"),
        "unexpected error: {text}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-07: a genome-wide SFS does not have trailing windows to
/// clip, so a valid header without contig lengths is sufficient. The
/// projection must still be validated against the population size.
#[test]
fn genome_wide_sfs_does_not_require_contig_lengths() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-nolen-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("sfs.vcf");
    std::fs::write(
        &vcf,
        VCF.replace("##contig=<ID=chr1,length=100000>", "##contig=<ID=chr1>"),
    )
    .expect("invariant: fixture write");
    let budget = ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(u64::MAX),
    };
    let report = run_popgen(Some(&vcf), &request(false), &budget).expect("invariant: run");
    let sfs = report.sfs.expect("invariant: SFS requested");
    assert_eq!(sfs.n_sites, vec![11]);
    for (got, want) in sfs.spectra[0].iter().zip(&EXPECTED) {
        assert!((got - want).abs() <= 1e-12 * want.max(1.0));
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-05: an empty successful run still has a request-shaped zero
/// spectrum, not `None`. Its writer therefore replaces the previous files.
#[test]
fn an_empty_sfs_run_returns_and_writes_a_zero_spectrum() {
    let dir = std::env::temp_dir().join(format!("reed-sfs-empty-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = fixture(&dir);
    let prefix = dir.join("out");
    let first = run_popgen(
        Some(&vcf),
        &request(false),
        &ThreadBudget {
            threads: 1,
            mem: MemBudget::from_bytes(u64::MAX),
        },
    )
    .expect("invariant: first run")
    .sfs
    .expect("invariant: first SFS");
    reed::sfs::write_projected_sfs(&prefix, &["all".to_string()], &first, false)
        .expect("invariant: first write");
    let dadi = dir.join("out_sfs/dadi/all-6.sfs");
    assert!(
        std::fs::read_to_string(&dadi)
            .expect("first dadi")
            .contains("0.571428")
    );

    let mut empty = request(false);
    empty.select.sites = vec![("chr2".to_string(), Some(1), Some(1))];
    let report = run_popgen(
        Some(&vcf),
        &empty,
        &ThreadBudget {
            threads: 2,
            mem: MemBudget::from_bytes(u64::MAX),
        },
    )
    .expect("invariant: empty run");
    let sfs = report
        .sfs
        .expect("invariant: an empty run still has an SFS");
    assert_eq!(sfs.n_sites, vec![0]);
    assert_eq!(sfs.spectra, vec![vec![0.0; 7]]);
    reed::sfs::write_projected_sfs(&prefix, &["all".to_string()], &sfs, false)
        .expect("invariant: zero write");
    let text = std::fs::read_to_string(&dadi).expect("invariant: replacement");
    assert_eq!(text, "7 unfolded\n0 0 0 0 0 0 0\n1 0 0 0 0 0 1\n");
    let _ = std::fs::remove_dir_all(&dir);
}
