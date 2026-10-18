//! `--hist` and `--singletons` end-to-end gates (doc_v14 02 §3): the
//! built-in metrics, the generic `INFO/`/`FMT/` tags, the exact quantile
//! lines, the accumulator laws under merging, the vcftools singleton layout
//! and thread determinism.

use std::path::{Path, PathBuf};

use reed::histogram::{HistKind, HistMetric, render_hist};
use reed::reports::{StatsReportRequest, run_stats};
use reed_engine::budget::ThreadBudget;

/// A VCF with GT/DP/GQ and an INFO/DP, so every histogram path has data.
fn fixture_text() -> String {
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=100000>\n\
         ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"site depth\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
         ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"quality\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n",
    );
    // Deterministic pseudo-random values: no RNG dependency, and the
    // expected quantiles below are hand-checked from this table.
    let rows = [
        (10, 5, 10, 0, 20),
        (20, 15, 20, 30, 40),
        (30, 25, 30, 60, 80),
        (40, 35, 40, 90, 99),
        (50, 45, 50, 50, 50),
        (60, 55, 60, 10, 5),
    ];
    for (i, &(dp0, dp1, info_dp, gq0, gq1)) in rows.iter().enumerate() {
        use std::fmt::Write as _;
        let _ = writeln!(
            text,
            "chr1\t{}\t.\tA\tG\t50\tPASS\tDP={info_dp}\tGT:DP:GQ\t0/1:{dp0}:{gq0}\t0/0:{dp1}:{gq1}",
            i + 1
        );
    }
    text
}

fn dir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("reed-hist-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("invariant: temp dir");
    d
}

fn write_vcf(path: &Path) {
    std::fs::write(path, fixture_text()).expect("invariant: fixture");
}

fn budget(threads: usize) -> ThreadBudget {
    ThreadBudget {
        threads,
        mem: reed_engine::budget::MemBudget::default_half_memory(),
    }
}

fn run_hist(vcf: &Path, prefix: &Path, metrics: Vec<HistMetric>, threads: usize) {
    let request = StatsReportRequest {
        hist: metrics,
        prefix: prefix.to_path_buf(),
        ..StatsReportRequest::default()
    };
    run_stats(Some(vcf), &request, &budget(threads)).expect("invariant: hist run");
}

fn read(path: &Path) -> String {
    std::fs::read_to_string(path).expect("invariant: read")
}

/// doc_v14 02 §3.1–3.2: the built-in metrics write one TSV per metric, the
/// integer ones carry exact quantiles, and `--hist` alone is a valid request.
#[test]
fn builtin_metrics_write_hist_files() {
    let d = dir("builtin");
    let vcf = d.join("in.vcf");
    write_vcf(&vcf);
    let prefix = d.join("qc");
    run_hist(
        &vcf,
        &prefix,
        vec![HistMetric::Qual, HistMetric::Dp, HistMetric::DpSite, HistMetric::Gq, HistMetric::Idl, HistMetric::Af],
        2,
    );
    let dp = read(&PathBuf::from(format!("{}.dp.hist", prefix.display())));
    assert!(dp.starts_with("bin_lower\tbin_upper\tcount\tpercent\n"), "{dp}");
    // 12 genotypes with DP 5,15,25,35,45,55,10,20,30,40,50,60 → median 30.
    assert!(dp.contains("#median\t30"), "{dp}");
    assert!(dp.contains("#mean\t32.500000"), "{dp}");
    // QUAL is labelled in quality units, not bcftools' internal iqual index
    // (`iqual = 1 + QUAL·10`, so the six QUAL=50 records land in [50.0, 50.1)).
    let qual = read(&PathBuf::from(format!("{}.qual.hist", prefix.display())));
    assert!(qual.contains("50.000000\t50.100000\t6\t"), "{qual}");
    // AF bins use bcftools' printed frequency ((i−1)/(m_af−1)): every site
    // here has one ALT copy, so the singletons fold into the first row.
    let af = read(&PathBuf::from(format!("{}.af.hist", prefix.display())));
    assert!(af.contains("0.000000\t0.010000\t6\t"), "{af}");
    // GQ is an integer metric with exact quantiles (median of the twelve
    // values 0,5,10,20,30,40,50,50,60,80,90,99 under the ceiling rule).
    let gq = read(&PathBuf::from(format!("{}.gq.hist", prefix.display())));
    assert!(gq.contains("#median\t40"), "{gq}");
    let _ = std::fs::remove_dir_all(&d);
}

/// doc_v14 02 §3.3 (guard 6): a generic `INFO/DP` request resolves against
/// the header and renders exact quantiles.
#[test]
fn generic_info_tag_histogram() {
    let d = dir("infotag");
    let vcf = d.join("in.vcf");
    write_vcf(&vcf);
    let prefix = d.join("qc");
    run_hist(
        &vcf,
        &prefix,
        vec![HistMetric::Tag {
            format: false,
            tag: "DP".into(),
            kind: HistKind::Int,
        }],
        2,
    );
    let text = read(&PathBuf::from(format!("{}.INFO_DP.hist", prefix.display())));
    assert!(text.contains("10\t10\t1"), "{text}");
    // The six INFO/DP values are 10,20,30,40,50,60: mean 35, median 30.
    assert!(text.contains("#mean\t35.000000"), "{text}");
    assert!(text.contains("#median\t30"), "{text}");
    let _ = std::fs::remove_dir_all(&d);
}

/// doc_v14 02 §3.4 (guard 7): the merge laws hold — threads 1 and 4 produce
/// byte-identical files.
#[test]
fn hist_is_deterministic_across_threads() {
    let d = dir("threads");
    let vcf = d.join("in.vcf");
    write_vcf(&vcf);
    let p1 = d.join("t1");
    let p4 = d.join("t4");
    let metrics = vec![
        HistMetric::Qual,
        HistMetric::Dp,
        HistMetric::Gq,
        HistMetric::Tag {
            format: false,
            tag: "DP".into(),
            kind: HistKind::Int,
        },
    ];
    run_hist(&vcf, &p1, metrics.clone(), 1);
    run_hist(&vcf, &p4, metrics, 4);
    for suffix in [".qual.hist", ".dp.hist", ".gq.hist", ".INFO_DP.hist"] {
        let a = read(&PathBuf::from(format!("{}{suffix}", p1.display())));
        let b = read(&PathBuf::from(format!("{}{suffix}", p4.display())));
        assert_eq!(a, b, "{suffix} differs across thread counts");
    }
    let _ = std::fs::remove_dir_all(&d);
}

/// `--singletons` (doc_v14 02 §2.2, C14-3): the vcftools layout
/// (`CHROM POS SINGLETON/DOUBLETON ALLELE INDV`), one row per singleton
/// allele; every site in this fixture has one ALT copy, so each is an `S`.
#[test]
fn singletons_match_the_vcftools_layout() {
    use reed::reports::{StatsReportRequest, run_stats};
    let d = dir("singletons");
    let vcf = d.join("in.vcf");
    write_vcf(&vcf);
    let prefix = d.join("qc");
    let request = StatsReportRequest {
        singletons: true,
        prefix: prefix.clone(),
        ..StatsReportRequest::default()
    };
    run_stats(Some(&vcf), &request, &budget(2)).expect("invariant: singletons run");
    let text = read(&PathBuf::from(format!("{}.singletons", prefix.display())));
    assert!(
        text.starts_with("CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\n"),
        "{text}"
    );
    assert!(text.contains("chr1\t1\tS\tG\ts0\n"), "{text}");
    assert_eq!(text.lines().count(), 7, "header + 6 singleton rows: {text}");
    let _ = std::fs::remove_dir_all(&d);
}

/// `--pop` (doc_v14 02 §2.2): the stratified afreq/vmiss rows equal a
/// per-population run's values, the singletons gain the POP column, and the
/// cohort files keep their 0.3.x bytes. The pops file puts s0 in popA and
/// s1 in popB, so each group sees one sample.
#[test]
fn pop_stratifies_freq_missing_and_singletons() {
    use reed::reports::{StatsReportRequest, run_stats};
    let d = dir("pop");
    let vcf = d.join("in.vcf");
    write_vcf(&vcf);
    let pops = d.join("pops.tsv");
    std::fs::write(&pops, "s0\tpopA\ns1\tpopB\n").expect("invariant: pops");
    let run = |prefix: &Path, populations: Option<reed::Populations>| {
        let request = StatsReportRequest {
            freq: true,
            missing: true,
            singletons: true,
            populations,
            prefix: prefix.to_path_buf(),
            ..StatsReportRequest::default()
        };
        run_stats(Some(&vcf), &request, &budget(2)).expect("invariant: pop run");
    };
    run(&d.join("cohort"), None);
    run(&d.join("strat"), Some(reed::Populations::File(pops)));
    // The cohort afreq/vmiss are byte-identical with and without --pop (the
    // singletons file gains the POP column by design).
    for suffix in [".afreq", ".vmiss"] {
        let a = read(&PathBuf::from(format!("{}{suffix}", d.join("cohort").display())));
        let b = read(&PathBuf::from(format!("{}{suffix}", d.join("strat").display())));
        assert_eq!(a, b, "{suffix}: --pop must not touch the cohort files");
    }
    // Stratified afreq: every site is 0/1 in s0 (popA) and 0/0 in s1
    // (popB), so ALT_FREQS is 0.5/0 with OBS_CT 2 per group.
    let afreq = read(&PathBuf::from(format!("{}.afreq.pop", d.join("strat").display())));
    assert!(
        afreq.starts_with("#CHROM\tID\tPOP\tREF\tALT\tALT_FREQS\tOBS_CT\n"),
        "{afreq}"
    );
    assert!(afreq.contains("chr1\tchr1_1\tpopA\tA\tG\t0.5\t2\n"), "{afreq}");
    assert!(afreq.contains("chr1\tchr1_1\tpopB\tA\tG\t0\t2\n"), "{afreq}");
    assert_eq!(afreq.lines().count(), 13, "header + 6 sites × 2 groups: {afreq}");
    // Stratified vmiss: no missing genotypes in the fixture.
    let vmiss = read(&PathBuf::from(format!("{}.vmiss.pop", d.join("strat").display())));
    assert!(
        vmiss.starts_with("#CHROM\tID\tPOP\tMISSING_CT\tOBS_CT\tF_MISS\n"),
        "{vmiss}"
    );
    assert!(vmiss.contains("chr1_1\tpopA\t0\t1\t0\n"), "{vmiss}");
    // Singletons: the POP column names the carrier's group.
    let singletons = read(&PathBuf::from(format!("{}.singletons", d.join("strat").display())));
    assert!(
        singletons.starts_with("CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\tPOP\n"),
        "{singletons}"
    );
    assert!(singletons.contains("chr1\t1\tS\tG\ts0\tpopA\n"), "{singletons}");
    let _ = std::fs::remove_dir_all(&d);
}

/// `--pop` rows naming samples absent from the header: an error without
/// `--force-samples`, a skip with it (the popg rule, doc_v14 02 §2.2).
#[test]
fn pop_unknown_samples_follow_the_force_rule() {
    use reed::reports::{StatsReportRequest, run_stats};
    let d = dir("popforce");
    let vcf = d.join("in.vcf");
    write_vcf(&vcf);
    let pops = d.join("pops.tsv");
    std::fs::write(&pops, "s0\tpopA\nghost\tpopB\n").expect("invariant: pops");
    let request = StatsReportRequest {
        freq: true,
        populations: Some(reed::Populations::File(pops.clone())),
        prefix: d.join("qc"),
        ..StatsReportRequest::default()
    };
    let err = run_stats(Some(&vcf), &request, &budget(2))
        .expect_err("an unknown sample without --force-samples");
    assert!(err.to_string().contains("ghost"), "{err}");
    let request = StatsReportRequest {
        freq: true,
        force_samples: true,
        populations: Some(reed::Populations::File(pops)),
        prefix: d.join("qc"),
        ..StatsReportRequest::default()
    };
    run_stats(Some(&vcf), &request, &budget(2)).expect("invariant: force run");
    let afreq = read(&PathBuf::from(format!("{}.afreq.pop", d.join("qc").display())));
    assert!(
        afreq.contains("popA"),
        "the known group still reports: {afreq}"
    );
    assert!(
        !afreq.contains("popB"),
        "the ghost's group has no samples and no rows: {afreq}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

/// The stratified outputs are deterministic across thread counts (the
/// ordered sink merges the rendered rows in chunk order).
#[test]
fn pop_reports_are_deterministic_across_threads() {
    use reed::reports::{StatsReportRequest, run_stats};
    let d = dir("popthreads");
    let vcf = d.join("in.vcf");
    write_vcf(&vcf);
    let pops = d.join("pops.tsv");
    std::fs::write(&pops, "s0\tpopA\ns1\tpopB\n").expect("invariant: pops");
    for threads in [1, 4] {
        let request = StatsReportRequest {
            freq: true,
            missing: true,
            singletons: true,
            populations: Some(reed::Populations::File(pops.clone())),
            prefix: d.join(format!("t{threads}")),
            ..StatsReportRequest::default()
        };
        run_stats(Some(&vcf), &request, &budget(threads)).expect("invariant: pop run");
    }
    for suffix in [".afreq.pop", ".vmiss.pop", ".singletons"] {
        let a = read(&PathBuf::from(format!("{}{suffix}", d.join("t1").display())));
        let b = read(&PathBuf::from(format!("{}{suffix}", d.join("t4").display())));
        assert_eq!(a, b, "{suffix} differs across thread counts");
    }
    let _ = std::fs::remove_dir_all(&d);
}

/// The renderer's own contract: the TSV shape and the `#missing` line.
#[test]
fn render_hist_shape() {
    let out = reed::histogram::HistOut {
        metric: HistMetric::Dp,
        int_bins: vec![(1, 3), (2, 1)],
        overflow: None,
        stats: Some(reed::histogram::HistStats {
            mean: 1.25,
            median: 1,
            p95: 2,
        }),
        float_bins: Vec::new(),
        n_missing: 4,
    };
    let text = render_hist(&out);
    assert!(text.starts_with("bin_lower\tbin_upper\tcount\tpercent\n"));
    assert!(text.contains("1\t1\t3\t75.000000"), "{text}");
    assert!(text.contains("#mean\t1.250000"), "{text}");
    assert!(text.contains("#median\t1"), "{text}");
    assert!(text.contains("#p95\t2"), "{text}");
    assert!(text.contains("#missing\t4"), "{text}");
}
