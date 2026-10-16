//! The facade's public-surface contract (doc_v11): what a crate that depends
//! on `reed` **alone** can reach, and how the API behaves at its defaults.
//!
//! This file is the compile-time gate for the v11 round. It imports only
//! `reed` (plus the test-only `reed_engine`/`reed_stats`/`reed_io` for
//! fixtures), so a name that disappears from the facade fails the build here
//! instead of only failing downstream users.

use std::io::Write;

use reed::{
    FilterReport, Populations, ReportRequest, ReportsOutcome, StatsRequest,
    assign_groups_in_memory, default_kinds, run_reports, run_stats,
};

// ---- A11-6: the four drivers and their reports live at the root ----
fn _drivers_and_reports_are_at_the_root(
    _: fn(
        Option<&std::path::Path>,
        &reed::spec::FilterSpec,
        &reed::FilterOptions,
    ) -> Result<FilterReport, reed::Error>,
    _: fn(
        Option<&std::path::Path>,
        &StatsRequest,
        &reed::budget::ThreadBudget,
    ) -> Result<reed::StatsReport, reed::Error>,
    _: fn(
        Option<&std::path::Path>,
        &ReportRequest,
        &reed::budget::ThreadBudget,
    ) -> Result<ReportsOutcome, reed::Error>,
    _: fn(
        &std::path::Path,
        reed::WriteIndex,
        Option<&std::path::Path>,
        usize,
    ) -> Result<reed::IndexOutcome, reed::Error>,
) {
}

// ---- A11-5: `Error` is sealed; classification goes through `kind()` ----
fn _error_is_classified_not_matched(e: &reed::Error) -> reed::ErrorKind {
    // No `match e { reed::Error::Io(..) => .. }` is possible: the variants
    // are private. `Display` + `kind()` + `source()` are the surface.
    let _ = e.to_string();
    let _ = std::error::Error::source(e);
    e.kind()
}

/// A BGZF fixture writer using only the facade's re-exported helpers.
fn bgzf(text: &str) -> Vec<u8> {
    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level 6");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(64) {
        deflater
            .deflate_block(block, &mut out)
            .expect("invariant: block fits");
    }
    out.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);
    out
}

const TEXT: &str = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
chr1\t50\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1
chr1\t250\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1
";

fn fixture(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-api-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join(name);
    let mut file = std::fs::File::create(&path).expect("invariant: fixture file");
    file.write_all(&bgzf(TEXT)).expect("invariant: write");
    path
}

fn budget() -> reed::budget::ThreadBudget {
    reed::budget::ThreadBudget {
        threads: 2,
        mem: reed::budget::MemBudget::from_bytes(u64::MAX),
    }
}

/// A11-1: a request built from `StatsRequest::default()` still emits its
/// window rows — the default `Levels` is the window level, not all-false.
#[test]
fn default_request_emits_window_rows() {
    let vcf = fixture("default.vcf.gz");
    let request = StatsRequest {
        kinds: vec![reed::StatKind::Pi {
            group: 0,
            per_bp: false,
        }],
        populations: Populations::All,
        grid: reed::WindowGrid::Bp {
            size: 100,
            step: 100,
        },
        ..StatsRequest::default()
    };
    assert_eq!(
        request.levels,
        reed::Levels {
            window: true,
            contig: false,
            genome: false
        },
        "the derived Default must be the window level"
    );
    let report = run_stats(Some(&vcf), &request, &budget()).expect("invariant: run");
    assert!(
        !report.windows.is_empty(),
        "a default request must emit window rows, got {}",
        report.windows.len()
    );
}

/// A11-2: the README "Library use" example, verbatim in spirit — the request
/// is borrowed by `run_stats` and still usable afterwards.
#[test]
fn readme_library_example_compiles_and_runs() {
    let vcf = fixture("readme.vcf.gz");
    let request = StatsRequest {
        kinds: vec![reed::StatKind::Pi {
            group: 0,
            per_bp: false,
        }],
        populations: Populations::All,
        grid: reed::WindowGrid::Bp {
            size: 100,
            step: 100,
        },
        ..StatsRequest::default()
    };
    let report = run_stats(Some(vcf.as_path()), &request, &budget()).expect("invariant: run");
    // The borrow survives the call: this line is the A11-2 regression.
    let tsv = reed::windows::render_tsv(&report.windows, &request.kinds, false);
    assert!(tsv.starts_with("contig\tstart\tend"));
    assert!(tsv.lines().count() > 1, "expected data rows:\n{tsv}");
}

/// A11-3: `to_args` reports "no flag spelling" as a usage error instead of
/// panicking, so it composes with the crate's own kind builders.
#[test]
fn to_args_errors_instead_of_panicking_on_unrenderable_kinds() {
    let names = vec!["A".to_string(), "B".to_string()];
    let request = StatsRequest {
        kinds: default_kinds(2),
        populations: Populations::All,
        grid: reed::WindowGrid::Bp {
            size: 100,
            step: 100,
        },
        ..StatsRequest::default()
    };
    let err = request
        .to_args(&names)
        .expect_err("default_kinds carries per-bp and SFS kinds: no flag spelling");
    assert_eq!(err.kind(), reed::ErrorKind::Usage, "{err}");

    // The windowed form still renders.
    let windowed = StatsRequest {
        kinds: reed::window_kinds(2, true, false, false, false, None, None),
        populations: Populations::All,
        grid: reed::WindowGrid::Bp {
            size: 100,
            step: 100,
        },
        ..StatsRequest::default()
    };
    let args = windowed.to_args(&names).expect("invariant: windowed form");
    assert!(args.contains(&"--pi".to_string()), "{args:?}");
}

/// A11-4: every facade constructor's error folds into `reed::Error`.
#[test]
fn facade_errors_converge_on_reed_error() {
    let levels = reed::Levels::parse("nope").expect_err("invariant: unknown level");
    assert_eq!(levels.kind(), reed::ErrorKind::Usage);
    let af = reed::spec::AfKind::parse("nope").expect_err("invariant: unknown af type");
    let as_reed: reed::Error = af.into();
    assert_eq!(as_reed.kind(), reed::ErrorKind::Usage);
    let gt = reed::spec::GtClassRule::parse("nope").expect_err("invariant: unknown class");
    let as_reed: reed::Error = gt.into();
    assert_eq!(as_reed.kind(), reed::ErrorKind::Usage);
}

/// A11-7: a summary-only request is not "empty".
#[test]
fn summary_only_report_request_is_not_empty() {
    let request = ReportRequest {
        summary: true,
        ..ReportRequest::default()
    };
    assert!(!request.is_empty());
    assert!(ReportRequest::default().is_empty());
}

/// A11-11: the in-memory assignment agrees with the file form byte for byte.
#[test]
fn in_memory_assignment_matches_the_file_form() {
    let vcf = fixture("assign.vcf.gz");
    let pairs = vec![
        ("s0".to_string(), "pop1".to_string()),
        ("s1".to_string(), "pop2".to_string()),
    ];
    let (masks, names) = assign_groups_in_memory(&pairs, &["s0".into(), "s1".into()], false)
        .expect("invariant: assignment");
    assert_eq!(names, ["pop1", "pop2"]);
    assert_eq!(masks, [Some(0), Some(1)]);

    let from_memory = StatsRequest {
        kinds: vec![reed::StatKind::Pi {
            group: 0,
            per_bp: false,
        }],
        populations: Populations::Assignment(pairs),
        grid: reed::WindowGrid::Bp {
            size: 100,
            step: 100,
        },
        ..StatsRequest::default()
    };
    let dir = vcf.parent().expect("invariant: dir");
    let pops = dir.join("pops.tsv");
    std::fs::write(&pops, "s0\tpop1\ns1\tpop2\n").expect("invariant: pops");
    let from_file = StatsRequest {
        populations: Populations::File(pops),
        ..from_memory.clone()
    };
    let a = run_stats(Some(&vcf), &from_memory, &budget()).expect("invariant: memory run");
    let b = run_stats(Some(&vcf), &from_file, &budget()).expect("invariant: file run");
    let kinds = from_memory.kinds.clone();
    assert_eq!(
        reed::windows::render_tsv(&a.windows, &kinds, false),
        reed::windows::render_tsv(&b.windows, &kinds, false),
        "in-memory and file assignments must produce the same table"
    );
}

/// A11-6: `run_reports` is reachable from the root and returns its report.
#[test]
fn run_reports_is_reachable_from_the_root() {
    let vcf = fixture("reports.vcf.gz");
    let prefix = vcf.parent().expect("invariant: dir").join("out");
    let request = ReportRequest {
        freq: true,
        prefix,
        ..ReportRequest::default()
    };
    let outcome = run_reports(Some(&vcf), &request, &budget()).expect("invariant: reports");
    assert_eq!(outcome.n_records, 2);
    assert_eq!(outcome.files.len(), 1);
}
