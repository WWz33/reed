//! R7 acceptance: the pixy-layout π/dxy path (04-algorithms.md §3–§4).
//! Hand-computed expectations; the software parity runs live in
//! `gold/popgen/`.

use std::collections::HashMap;
use std::io::Write;

use reed::stats::{pixy_kinds, run_stats, InvariantMode, Levels, PixyOptions, StatsRequest};
use reed_engine::budget::{MemBudget, ThreadBudget};
use reed_stats::request::StatKind;

fn bgzf(text: &str) -> Vec<u8> {
    let mut deflater =
        reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6 is valid");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(64) {
        deflater.deflate_block(block, &mut out).expect("invariant: block fits");
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

fn request(invariant: InvariantMode) -> StatsRequest {
    StatsRequest {
        kinds: pixy_kinds(2, true, true),
        gt_groups: vec![Some(0), Some(0), Some(1), Some(1)],
        window_size: 10,
        step: 10,
        levels: Levels { window: true, contig: false, genome: false },
        skip_empty: false,
        contig_lengths: HashMap::new(),
        pixy: Some(PixyOptions { invariant }),
    }
}

fn budget() -> ThreadBudget {
    ThreadBudget { threads: 2, mem: MemBudget::from_half_memory_for_tests() }
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
        let kind = StatKind::Pi { group: g, per_bp: false };
        w0.accs[g].value_variant_only(&kind, w0.bp)
    };
    assert_eq!(w0.start1, 1);
    assert_eq!(w0.end0, 10);
    assert!((pi(0).expect("pi") - 0.05).abs() < 1e-15, "{:?}", pi(0));
    assert!((pi(1).expect("pi") - 0.05).abs() < 1e-15, "{:?}", pi(1));
    let dxy_kind = StatKind::Dxy { a: 0, b: 1, per_bp: false };
    let dxy = w0.accs[2].value_variant_only(&dxy_kind, w0.bp).expect("dxy");
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
    let pi = w0.accs[0].value(&StatKind::Pi { group: 0, per_bp: false }, w0.bp);
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
    let kinds = pixy_kinds(2, true, true);
    reed::pixy::write(
        &prefix,
        &report,
        &["pop1".to_string(), "pop2".to_string()],
        &kinds,
        PixyOptions { invariant: InvariantMode::VariantOnly },
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
