//! `run_stats` end-to-end: hand-computed window values, thread-count
//! determinism (§7.6), the contig-revisit refusal, `--skip-empty`
//! (§8.5) and the window/contig/genome levels (§8.7).

use std::collections::HashMap;
use std::io::Write;

use reed::stats::{assign_groups, default_kinds, run_stats, Levels, StatsRequest};
use reed_engine::budget::{MemBudget, ThreadBudget};
use reed_stats::fst::FstKind;
use reed_stats::request::StatKind;

/// Compresses `text` into a BGZF stream cut at `block_data`-byte blocks —
/// boundaries fall mid-line, forcing seams.
fn bgzf(text: &str, block_data: usize) -> Vec<u8> {
    let mut deflater =
        reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6 is valid");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(block_data) {
        deflater
            .deflate_block(block, &mut out)
            .expect("invariant: block fits");
    }
    out.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
    out
}

fn write_fixture(name: &str, text: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-stats-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join(name);
    let mut file = std::fs::File::create(&path).expect("invariant: fixture file");
    file.write_all(&bgzf(text, 64)).expect("invariant: write");
    drop(file);
    path
}

fn fixture_text() -> String {
    // 3 diploid samples, one per group; windows of 100 bp.
    // chr1:50 — A=0/1, B=1/1, C=0/0; chr1:80 — A=0/0, B=0/1, C=./.
    // chr1:90 — multiallelic (skipped); chr1:250 — A=1/1, B=0/0, C=0/1.
    "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
chr1\t50\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\t0/0
chr1\t80\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t./.
chr1\t90\t.\tA\tG,T\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0
chr1\t250\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t0/0\t0/1
"
    .to_string()
}

fn request(kinds: Vec<StatKind>, levels: Levels, skip_empty: bool) -> StatsRequest {
    request_with(kinds, levels, skip_empty, HashMap::new())
}

fn request_with(
    kinds: Vec<StatKind>,
    levels: Levels,
    skip_empty: bool,
    contig_lengths: HashMap<Vec<u8>, u64>,
) -> StatsRequest {
    StatsRequest {
        kinds,
        gt_groups: vec![Some(0), Some(1), Some(2)],
        window_size: 100,
        step: 100,
        levels,
        skip_empty,
        contig_lengths,
        pixy: None,
    }
}

/// Splits the TSV into rows of cells.
fn rows(tsv: &str) -> Vec<Vec<&str>> {
    tsv.lines()
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect()
}

/// Asserts one cell is `.` or within 1e-12 of the expected f64.
fn assert_cell(row: &[&str], column: usize, expected: Option<f64>, what: &str) {
    match (expected, row[column]) {
        (None, ".") => {}
        (Some(x), cell) => {
            let got: f64 = cell.parse().unwrap_or_else(|_| panic!("{what}: {cell}"));
            assert!(
                (got - x).abs() <= 1e-12 * x.abs().max(1.0),
                "{what}: got {got}, expected {x}"
            );
        }
        _ => panic!("{what}: expected {expected:?}, got {:?}", row[column]),
    }
}

const KINDS_TEST: [StatKind; 8] = [
    StatKind::Pi { group: 0, per_bp: false },
    StatKind::Pi { group: 0, per_bp: true },
    StatKind::Dxy { a: 0, b: 1, per_bp: false },
    StatKind::Fst { a: 0, b: 1, kind: FstKind::Hudson },
    StatKind::Sfs { group: 0, folded: false },
    StatKind::ThetaW { group: 0, per_bp: false },
    StatKind::ThetaW { group: 0, per_bp: true },
    StatKind::TajimaD { group: 0, min_sites: 0 },
];

#[test]
fn windows_match_hand_computation_and_thread_counts_agree() {
    let path = write_fixture("hand.vcf.gz", &fixture_text());

    let mut expected = String::from(
        "contig\tstart\tend\tpi.pair.0\tpi.bp.0\tdxy.pair.0.1\tfst.hudson.0.1\t\
         sfs.unfolded.0\ttheta_w.total.0\ttheta_w.bp.0\ttajima_d.0\n",
    );
    // Window [1,100]: π = 1/2; π_bp = (1 + 0)/100; dxy = 4/8; Hudson = 0/1;
    // SFS alt=1 once; θ (pixy: per-site n) = (1/a(2) + 0/a(2))/2 sites = 0.5,
    // per bp 0.005; n = 2 gives Tajima variance 0 → '.'.
    expected.push_str("chr1\t1\t100\t0.5\t0.01\t0.5\t0\t1\t0.5\t0.005\t.\n");
    // Window [101,200]: no sites — undefined statistics stay '.', the SFS
    // bins are all zero.
    expected.push_str("chr1\t101\t200\t.\t.\t.\t.\t0\t.\t.\t.\n");
    // Window [201,300]: the single 1/1 site is monomorphic-ALT for A
    // (k = 1 distinct allele → 0 mutations), so π = 0 and θ = 0;
    // dxy = 4/4; Hudson num/den = 1/1.
    expected.push_str("chr1\t201\t300\t0\t0\t1\t1\t0\t0\t0\t.\n");

    for threads in [1_usize, 2, 4, 7] {
        let budget = ThreadBudget { threads, mem: MemBudget::from_bytes(u64::MAX) };
        let report = run_stats(
            Some(&path.as_path()),
            request(KINDS_TEST.to_vec(), Levels::default(), false),
            &budget,
        )
        .expect("invariant: fixture run succeeds");
        assert_eq!(report.tsv, expected, "threads = {threads}");
        assert_eq!(report.n_windows, 3);
        assert_eq!(report.n_records, 3);
        assert_eq!(report.n_multiallelic_skipped, 1);
    }
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

#[test]
fn skip_empty_drops_window_rows_only() {
    let path = write_fixture("skip.vcf.gz", &fixture_text());
    let budget = ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) };
    let report = run_stats(
        Some(&path.as_path()),
        request(KINDS_TEST.to_vec(), Levels::default(), true),
        &budget,
    )
    .expect("invariant: fixture run succeeds");
    let got = rows(&report.tsv);
    assert_eq!(got.len(), 2, "the empty middle window is skipped");
    assert_eq!(got[0][..3], ["chr1", "1", "100"]);
    assert_eq!(got[1][..3], ["chr1", "201", "300"]);
    assert_eq!(report.n_windows, 2);
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

#[test]
fn contig_and_genome_levels_aggregate_the_bins() {
    let path = write_fixture("levels.vcf.gz", &fixture_text());
    let budget = ThreadBudget { threads: 3, mem: MemBudget::from_bytes(u64::MAX) };
    let levels = Levels::parse("window,contig,genome").expect("invariant: levels");
    let report = run_stats(Some(&path.as_path()), request(KINDS_TEST.to_vec(), levels, false), &budget)
        .expect("invariant: fixture run succeeds");

    let got = rows(&report.tsv);
    assert_eq!(got.len(), 5, "3 window rows + contig row + genome row");
    // Window order is preserved, then the contig row, then the genome row.
    assert_eq!(got[0][0], "chr1");
    assert_eq!(got[2][..3], ["chr1", "201", "300"]);
    assert_eq!(got[3][..3], ["chr1", "1", "300"], "contig span");
    assert_eq!(got[4][..3], ["genome", "0", "300"], "genome row: summed span");

    // Contig totals by hand: π diffs 1 over comps 1+1+1 = 3 (each site
    // has n = 2 chromosomes → one comparable pair);
    // dxy 8/12; Hudson num 1 over den 2; SFS [1]; θ (pixy) = raw 1 over the
    // 3 sites with calls = 1/3; per-bp divides by the 300 bp span.
    let contig = &got[3];
    assert_cell(contig, 3, Some(1.0 / 3.0), "contig π");
    assert_cell(contig, 4, Some(1.0 / 300.0), "contig π_bp");
    assert_cell(contig, 5, Some(2.0 / 3.0), "contig dxy");
    assert_cell(contig, 6, Some(0.5), "contig Hudson");
    assert_eq!(contig[7], "1", "contig SFS");
    assert_cell(contig, 8, Some(1.0 / 3.0), "contig θ");
    assert_cell(contig, 9, Some(1.0 / 900.0), "contig θ_bp");
    // The genome row's values equal the contig row's (one contig in the
    // fixture); the label/span columns differ by construction.
    for (a, b) in got[3][3..].iter().zip(&got[4][3..]) {
        assert_eq!(a, b, "genome values mirror the single contig");
    }
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

#[test]
fn levels_parse_rejects_unknown_tokens() {
    assert!(Levels::parse("window,genome").is_ok());
    assert!(Levels::parse("bogus").is_err());
    assert_eq!(
        Levels::parse("window").expect("invariant: valid"),
        Levels { window: true, contig: false, genome: false }
    );
}

#[test]
fn contig_revisit_is_refused() {
    let text = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
chr1\t50\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0
chr2\t50\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t1/1\t0/0
chr1\t60\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t1/1
";
    let path = write_fixture("revisit.vcf.gz", text);
    let budget = ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) };
    let err = run_stats(
        Some(&path.as_path()),
        request(KINDS_TEST.to_vec(), Levels::default(), false),
        &budget,
    )
    .expect_err("invariant: revisit is refused");
    assert!(err.to_string().contains("revisited"), "{err}");
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

#[test]
fn assign_groups_maps_header_order() {
    let dir = std::env::temp_dir().join(format!("reed-groups-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join("groups.tsv");
    // Pixy populations-file format: two columns, no header row.
    std::fs::write(&path, "s1\tpop2\ns0\tpop1\ns2\tpop1\n")
        .expect("invariant: write");
    let (gt_groups, names) =
        assign_groups(&path, &["s0".into(), "s1".into(), "s2".into(), "s9".into()], false)
            .expect("invariant: parse");
    // First appearance order: pop2 (s1) then pop1 (s0, s2); s9 unlisted.
    assert_eq!(names, ["pop2", "pop1"]);
    assert_eq!(gt_groups, [Some(1), Some(0), Some(1), None]);
    let _ = std::fs::remove_dir_all(&dir);
}

/// Regression (01-review D9): a file row naming a sample absent from the
/// header is an error (pixy behavior), not a silent skip; --force-samples
/// downgrades it to a skip. A header row is DATA now: `sample` is not a
/// real sample, so the header form must fail too.
#[test]
fn assign_groups_rejects_unknown_samples() {
    let dir = std::env::temp_dir().join(format!("reed-groups-unknown-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join("groups.tsv");
    std::fs::write(&path, "s0\tpop1\nghost\tpop2\n").expect("invariant: write");
    let err = assign_groups(&path, &["s0".into(), "s1".into()], false).unwrap_err();
    assert!(err.to_string().contains("ghost"), "{err}");
    let (gt_groups, names) =
        assign_groups(&path, &["s0".into(), "s1".into()], true).expect("invariant: forced");
    assert_eq!(names, ["pop1"]);
    assert_eq!(gt_groups, [Some(0), None]);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn default_kinds_covers_every_family() {
    let kinds = default_kinds(3);
    let keys: Vec<String> = kinds.iter().map(|k| k.key().0).collect();
    assert!(keys.iter().any(|k| k.starts_with("pi.")));
    assert!(keys.iter().any(|k| k.starts_with("dxy.")));
    assert!(keys.iter().any(|k| k.starts_with("fst.hudson.")));
    assert!(keys.iter().any(|k| k.starts_with("fst.weir_cockerham_84.")));
    assert!(keys.iter().any(|k| k.starts_with("fst.nei_86.")));
    assert!(keys.iter().any(|k| k.starts_with("pbs.yi_2010.")));
    assert!(keys.iter().any(|k| k.starts_with("sfs.unfolded.")));
    assert!(keys.iter().any(|k| k.starts_with("theta_w.")));
    assert!(keys.iter().any(|k| k.starts_with("tajima_d.")));
    // Keys are unique.
    let mut sorted = keys.clone();
    sorted.sort();
    sorted.dedup();
    assert_eq!(sorted.len(), keys.len());
}

/// Regression (01-review D7): a contig without a declared length must fail
/// the run (the trailing-window clip and per-bp denominators would be
/// fabricated from the last position otherwise), and --contig-lengths must
/// repair it with identical window values.
#[test]
fn missing_contig_length_is_an_error_and_lengths_file_repairs_it() {
    // Same records as fixture_text() but the header omits ##contig lengths.
    let text = "\
##fileformat=VCFv4.3
##contig=<ID=chr1>
##contig=<ID=chr2>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
chr1\t50\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\t0/0
chr1\t250\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t0/0\t0/1
";
    let path = write_fixture("d7-no-length", &text);
    let gt_groups = vec![Some(0), Some(1), Some(2)];

    let err = run_stats(
        Some(&path.as_path()),
        request(vec![StatKind::Pi { group: 0, per_bp: true }], Levels::default(), false),
        &ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) },
    )
    .expect_err("invariant: missing length is refused");
    assert!(err.to_string().contains("no ##contig length"), "{err}");

    // The .fai-style supplement repairs it: same values as a declared run.
    let mut lengths = HashMap::new();
    lengths.insert(b"chr1".to_vec(), 1_000_u64);
    let with_file = run_stats(
        Some(&path.as_path()),
        request_with(vec![StatKind::Pi { group: 0, per_bp: true }], Levels::default(), false, lengths),
        &ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) },
    )
    .expect("invariant: lengths supplement works");

    let declared = write_fixture("d7-declared", &fixture_text());
    let reference = run_stats(
        Some(&declared.as_path()),
        request(vec![StatKind::Pi { group: 0, per_bp: true }], Levels::default(), false),
        &ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) },
    )
    .expect("invariant: declared run");
    assert_eq!(
        with_file.tsv, reference.tsv,
        "supplemented lengths must match a header-declared run"
    );
    let _ = gt_groups;
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
    let _ = std::fs::remove_dir_all(declared.parent().expect("invariant: temp dir path"));
}

/// Regression (01-review D8): non-ascending positions and positions beyond
/// the declared contig length are refused instead of silently dropping
/// sites into already-emitted bins.
#[test]
fn position_contract_violations_are_refused() {
    let unsorted = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
chr1\t80\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0
chr1\t50\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0
";
    let path = write_fixture("d8-unsorted", &unsorted);
    let err = run_stats(
        Some(&path.as_path()),
        request(vec![StatKind::Pi { group: 0, per_bp: false }], Levels::default(), false),
        &ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) },
    )
    .expect_err("invariant: unsorted input is refused");
    assert!(err.to_string().contains("ascending"), "{err}");
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));

    let beyond = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
chr1\t2000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0
";
    let path = write_fixture("d8-beyond", &beyond);
    let err = run_stats(
        Some(&path.as_path()),
        request(vec![StatKind::Pi { group: 0, per_bp: false }], Levels::default(), false),
        &ThreadBudget { threads: 2, mem: MemBudget::from_bytes(u64::MAX) },
    )
    .expect_err("invariant: position beyond the contig length is refused");
    assert!(err.to_string().contains("exceeds contig"), "{err}");
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}
