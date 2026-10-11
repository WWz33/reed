//! R7 acceptance: the pixy-layout π/dxy path (04-algorithms.md §3–§4).
//! Hand-computed expectations; the software parity runs live in
//! `gold/popgen/`.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::Write;

use reed::{
    InvariantMode, Levels, Populations, StatsRequest, WindowLayoutOptions, run_stats, window_kinds,
};
use reed_engine::budget::{MemBudget, ThreadBudget};
use reed_stats::request::StatKind;

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

fn write_fixture(name: &str, text: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-pixy-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join("in.vcf.gz");
    let mut file = std::fs::File::create(&path).expect("invariant: fixture file");
    file.write_all(&bgzf(text)).expect("invariant: write");
    drop(file);
    path
}

/// 4 samples, 2 populations (pop1 = s0,s1; pop2 = s2,s3).
/// Sites: pos 1 (pop1 0/1,0/0 | pop2 0/0,0/1), pos 5 (all 0/0), pos 25 (pop1 0/1,1/1 | pop2 1/1,0/0).
/// Window 10 over contig length 30 → windows [1,10], [11,20], [21,30].
const FIXTURE: &str = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=30>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/1
chr1\t5\t.\tA\t.\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0
chr1\t25\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\t1/1\t0/0
";

fn pops_two_file() -> std::path::PathBuf {
    static PATH: std::sync::OnceLock<std::path::PathBuf> = std::sync::OnceLock::new();
    PATH.get_or_init(|| {
        let dir = std::env::temp_dir().join(format!("reed-pixy-pops2-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let path = dir.join("pops.tsv");
        std::fs::write(&path, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: pops");
        path
    })
    .clone()
}

fn request(invariant: InvariantMode) -> StatsRequest {
    StatsRequest {
        kinds: window_kinds(2, true, true, false, false, None, None),
        populations: Populations::File(pops_two_file()),
        grid: reed::WindowGrid::Bp { size: 10, step: 10 },
        levels: Levels {
            window: true,
            contig: false,
            genome: false,
        },
        skip_empty: false,
        contig_lengths: HashMap::new(),
        layout: Some(WindowLayoutOptions { invariant }),
        ..StatsRequest::default()
    }
}

fn budget() -> ThreadBudget {
    ThreadBudget {
        threads: 2,
        mem: MemBudget::default_half_memory(),
    }
}

#[test]
fn variant_only_pi_dxy_hand_computed() {
    let path = write_fixture("vonly", FIXTURE);
    let report = run_stats(Some(&path), request(InvariantMode::VariantOnly), &budget())
        .expect("invariant: variant-only run");

    // Window [1,10]: present sites 1 and 5.
    // pop1: site1 n=4 k=1 → diffs 3, comps 6; site5 → 0, 6; |S|=2, L=10
    //   → comps + 8·C(4,2) = 12 + 48 = 60 → π = 3/60 = 0.05
    // pop2: site1 diffs 3, site5 0 → same = 0.05
    // dxy pop1×pop2: site1 a=(3 ref,1 alt) b=(3 ref,1 alt) → cross = 3·1+1·3 = 6,
    //   comps = 4·4 = 16; site5 cross 0, comps 16; |S|=2, absent 8
    //   → 6 / (32 + 8·16) = 6/160 = 0.0375
    let w0 = &report.windows[0];
    let pi = |g: usize| {
        let kind = StatKind::Pi {
            group: g,
            per_bp: false,
        };
        w0.accs[g].value_variant_only(&kind, w0.bp)
    };
    assert_eq!(w0.start1, 1);
    assert_eq!(w0.end0, 10);
    assert!((pi(0).expect("pi") - 0.05).abs() < 1e-15, "{:?}", pi(0));
    assert!((pi(1).expect("pi") - 0.05).abs() < 1e-15, "{:?}", pi(1));
    let dxy_kind = StatKind::Dxy {
        a: 0,
        b: 1,
        per_bp: false,
    };
    let dxy = w0.accs[2]
        .value_variant_only(&dxy_kind, w0.bp)
        .expect("dxy");
    assert!((dxy - 0.0375).abs() < 1e-15, "{dxy}");
    // Window [11,20] has no site: 8→ absent 10 → π = 0/60 = 0 for present
    // sites absent... the accumulator is empty, so the completed value is 0
    // and the writer prints NA (checked in the writer test below).
    assert_eq!(report.windows[1].accs[0].n_present(), Some(0));
    // Grid from 1: window starts are 1, 11, 21 regardless of site positions.
    let starts: Vec<u64> = report.windows.iter().map(|w| w.start1).collect();
    assert_eq!(starts, [1, 11, 21]);
    assert_eq!(report.n_invariant_sites, 1, "site 5 is site-wide invariant");
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

#[test]
fn all_sites_default_matches_the_all_sites_definition() {
    let path = write_fixture("allsites", FIXTURE);
    let report = run_stats(Some(&path), request(InvariantMode::Require), &budget())
        .expect("invariant: the fixture has an invariant site");
    // Window [1,10] with invariant site 5 counted: comps 12 (unchanged),
    // diffs 3 → π = 0.25 for pop1 (no completion: the file is all-sites).
    let w0 = &report.windows[0];
    let pi = w0.accs[0].value(
        &StatKind::Pi {
            group: 0,
            per_bp: false,
        },
        w0.bp,
    );
    assert!((pi.expect("pi") - 0.25).abs() < 1e-15, "{pi:?}");
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

#[test]
fn variant_only_input_without_declaration_is_refused() {
    let mut text = String::from(FIXTURE);
    text = text.replace("chr1\t5\t.\tA\t.\t", "chr1\t5\t.\tA\tG\t"); // no invariant site left
    let path = write_fixture("nodecl", &text);
    let err = run_stats(Some(&path), request(InvariantMode::Require), &budget())
        .expect_err("invariant: missing invariant sites must fail");
    assert!(err.to_string().contains("no invariant sites"), "{err}");
    // Declaring it works.
    run_stats(Some(&path), request(InvariantMode::VariantOnly), &budget())
        .expect("invariant: declared variant-only run");
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

#[test]
fn writer_emits_pixy_columns_with_na_for_empty_windows() {
    let path = write_fixture("writer", FIXTURE);
    let report = run_stats(Some(&path), request(InvariantMode::VariantOnly), &budget())
        .expect("invariant: run");
    let prefix = path.parent().expect("invariant: temp dir path").join("out");
    let kinds = window_kinds(2, true, true, false, false, None, None);
    reed::windows::write(
        &prefix,
        &report,
        &["pop1".to_string(), "pop2".to_string()],
        &kinds,
        WindowLayoutOptions {
            invariant: InvariantMode::VariantOnly,
        },
        "reed stats test",
        "test",
    )
    .expect("invariant: write");
    let pi = std::fs::read_to_string(format!("{}_pi.txt", prefix.display()))
        .expect("invariant: pi file");
    let mut lines = pi.lines();
    assert_eq!(
        lines.next().expect("header"),
        "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing"
    );
    let first: Vec<&str> = lines.next().expect("row").split('\t').collect();
    assert_eq!(&first[..4], ["pop1", "chr1", "1", "10"]);
    assert!((first[4].parse::<f64>().expect("avg_pi") - 0.05).abs() < 1e-15);
    assert_eq!(&first[5..], ["2", "3", "12", "0"]);
    // Rows run window-major: pop2 of window [1,10] comes next.
    let second: Vec<&str> = lines.next().expect("row").split('\t').collect();
    assert_eq!(&second[..4], ["pop2", "chr1", "1", "10"]);
    // Window [11,20] is empty → NA with zero counts (pixy layout).
    let empty: Vec<&str> = lines.next().expect("row").split('\t').collect();
    assert_eq!(&empty[..4], ["pop1", "chr1", "11", "20"]);
    assert_eq!(&empty[4..], ["NA", "0", "0", "0", "0"]);
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

/// Hand-computed Watterson θ and Tajima's D on the pixy definitions
/// (04-algorithms.md §3), using only the first site's window.
///
/// Site 1 for pop1 (s0 = 0/1, s1 = 0/0): n = 4 observed alleles, k = 2 →
/// one mutation; a(4) = 1 + 1/2 + 1/3 = 11/6 → raw θ = 6/11 = 0.545454…;
/// window [1,10] also contains the invariant site 5 (n = 4, no mutation),
/// so no_sites = 2 and avg θ = raw/2 = 3/11.
/// raw π = Σ per-site π = (3·1)/6 + 0 = 0.5; S = 1; mean n = round(8/2) = 4
/// → e1 = 0.0055092…, stdev = 0.0742270 → D = (0.5 − 0.5454545)/0.0742270.
#[test]
fn watterson_theta_and_tajima_d_hand_computed() {
    let path = write_fixture("theta", FIXTURE);
    let mut request = request(InvariantMode::Bypass);
    request.kinds = window_kinds(2, false, false, true, true, None, None);
    let report = run_stats(Some(&path), request, &budget()).expect("invariant: theta run");

    let w0 = &report.windows[0];
    let theta_kind = StatKind::ThetaW {
        group: 0,
        per_bp: false,
    };
    let theta = w0.accs[0].value(&theta_kind, w0.bp).expect("theta");
    assert!((theta - 3.0 / 11.0).abs() < 1e-15, "{theta}");
    let (raw, no_sites, no_var, _weighted) = w0.accs[0].theta_columns().expect("columns");
    assert!((raw - 6.0 / 11.0).abs() < 1e-15, "{raw}");
    assert_eq!((no_sites, no_var), (2, 1));

    let tajima_kind = StatKind::TajimaD {
        group: 0,
        min_sites: 0,
    };
    // kinds = [theta pop1, theta pop2, tajima pop1, tajima pop2].
    let d = w0.accs[2].value(&tajima_kind, w0.bp).expect("tajima");
    let expected = (0.5 - 6.0 / 11.0) / 0.074_226_961_902_520_39;
    assert!((d - expected).abs() < 1e-12, "{d} vs {expected}");
    let (raw_pi, raw_theta, stdev, sites) = w0.accs[2].tajima_columns().expect("columns");
    assert!((raw_pi.expect("raw_pi") - 0.5).abs() < 1e-15);
    assert!((raw_theta.expect("raw_theta") - 6.0 / 11.0).abs() < 1e-15);
    assert!((stdev.expect("stdev") - 0.074_226_961_902_520_39).abs() < 1e-15);
    assert_eq!(sites, 2);
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

/// Hand-computed Fst on the pixy definitions (04-algorithms.md §3).
///
/// Pop1 = s0, s1; pop2 = s2, s3 (2 diploids each).
/// Site 1: ac1 = ac2 = [3 ref, 1 alt] → p = (0.75, 0.25); h_between =
///   1 − (0.75² + 0.25²) = 0.375; h1 = h2 = (4/3)·0.375 = 0.5 →
///   Hudson num = 0.375 − 0.5 = −0.125, den = 0.375.
/// Site 5: all ref → num = 0, den = 0.
/// Site 25: ac1 = ac2 = [1 ref, 3 alt] → same as site 1 → num −0.125, den 0.375.
/// Window [1,10] = sites 1 + 5 → Hudson = −0.125/0.375 = −1/3.
/// WC84 at site 1 (n̄ = n_C = 2, p̄ = 0.75, s² = 0, h̄ = 0.5): a = −0.125,
///   b = 0, c = 0.5 → −0.125/0.375 = −1/3 (site 5 adds zeros).
#[test]
fn fst_hudson_and_wc84_hand_computed() {
    let path = write_fixture("fst", FIXTURE);
    for (kind, name) in [
        (reed::FstKind::Hudson, "hudson"),
        (reed::FstKind::WeirCockerham84, "wc84"),
    ] {
        let mut request = request(InvariantMode::Bypass);
        request.kinds = window_kinds(2, false, false, false, false, Some(kind), None);
        let report = run_stats(Some(&path), request, &budget()).expect("invariant: fst run");
        let w0 = &report.windows[0];
        let fst_kind = StatKind::Fst { a: 0, b: 1, kind };
        let value = w0.accs[0].value(&fst_kind, w0.bp).expect("fst value");
        assert!(
            (value + 1.0 / 3.0).abs() < 1e-12,
            "{name}: got {value}, expected -1/3"
        );
        let (_, no_snps) = w0.accs[0].fst_columns(&fst_kind).expect("columns");
        assert_eq!(no_snps, 2, "{name}: window [1,10] holds 2 sites");
    }
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}

/// PBS regression on an explicit 3-population fixture. Expected values come
/// from scikit-allel (`allel.pbs(..., window_size=3, normed=...)` and
/// `allel.hudson_fst` per window), i.e. the reference named by
/// 04-algorithms.md §5; `gold/popgen/compare_pbs.py` re-runs that check.
const PBS_FIXTURE: &str = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tp1s0\tp1s1\tp1s2\tp1s3\tp2s0\tp2s1\tp2s2\tp2s3\tp3s0\tp3s1\tp3s2\tp3s3
chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\t1/1\t0/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0
chr1\t20\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/1\t1/1\t0/0\t0/0\t0/1\t0/1\t0/0\t0/0
chr1\t30\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/0\t0/0\t0/0\t0/0\t0/0\t1/1\t0/1\t0/0\t0/0
chr1\t110\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t0/0\t1/1\t0/0\t0/0\t0/1\t0/0\t0/1\t1/1\t0/0\t0/0
chr1\t120\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1\t0/1\t0/0
chr1\t130\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0\t0/1\t0/0\t1/1\t0/0\t0/0\t0/0\t0/0\t0/0
";

#[test]
fn pbs_matches_scikit_allel_reference() {
    let path = write_fixture("pbs", PBS_FIXTURE);
    // Groups: pop1 = samples 0..4, pop2 = 4..8, pop3 = 8..12.
    let pops = {
        let dir = std::env::temp_dir().join(format!("reed-pixy-pbs-pops-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let path = dir.join("pops.tsv");
        let mut text = String::new();
        for i in 0..12 {
            writeln!(text, "p{}s{}\tpop{}", i / 4 + 1, i % 4, i / 4 + 1).expect("invariant: pops");
        }
        std::fs::write(&path, text).expect("invariant: pops");
        path
    };
    for (kind, expected, expected_fst) in [
        (
            reed::PbsKind::Yi2010,
            [0.002_309_472_928_147_316_f64, -0.051_839_214_737_420_786],
            [
                (
                    0.142_857_142_857_142_9,
                    0.004_608_294_930_875_617,
                    0.142_857_142_857_142_9,
                ),
                (
                    0.058_823_529_411_764_73,
                    -0.124_423_963_133_640_51,
                    0.151_515_151_515_151_55,
                ),
            ],
        ),
        (
            reed::PbsKind::PbsN1,
            [0.001_997_019_026_245_494, -0.046_598_561_709_660_744],
            [
                (
                    0.142_857_142_857_142_9,
                    0.004_608_294_930_875_617,
                    0.142_857_142_857_142_9,
                ),
                (
                    0.058_823_529_411_764_73,
                    -0.124_423_963_133_640_51,
                    0.151_515_151_515_151_55,
                ),
            ],
        ),
    ] {
        let mut request = request(InvariantMode::Bypass);
        request.populations = Populations::File(pops.clone());
        request.grid = reed::WindowGrid::Bp {
            size: 100,
            step: 100,
        };
        request.kinds = window_kinds(3, false, false, false, false, None, Some(kind));
        let report = run_stats(Some(&path), request, &budget()).expect("invariant: pbs run");
        assert_eq!(report.windows.len(), 2);
        let pbs_kind = StatKind::Pbs {
            a: 0,
            b: 1,
            c: 2,
            kind,
        };
        for (row, (want, want_fst)) in report
            .windows
            .iter()
            .zip(expected.iter().zip(&expected_fst))
        {
            let (value, no_snps, (f_ab, f_ac, f_bc)) =
                row.accs[0].pbs_columns(&pbs_kind).expect("pbs columns");
            assert!(
                (value.expect("pbs") - want).abs() < 1e-12,
                "{kind:?}: {} vs {want}",
                value.expect("pbs")
            );
            assert_eq!(no_snps, 3);
            let got = (
                f_ab.expect("f_ab"),
                f_ac.expect("f_ac"),
                f_bc.expect("f_bc"),
            );
            for (g, w) in [got.0, got.1, got.2]
                .iter()
                .zip([want_fst.0, want_fst.1, want_fst.2])
            {
                assert!((g - w).abs() < 1e-12, "pair Fst {g} vs {w}");
            }
        }
    }
    let _ = std::fs::remove_dir_all(path.parent().expect("invariant: temp dir path"));
}
