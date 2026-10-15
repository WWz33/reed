//! doc_v12 R0: the INFO-field defect fixtures (01-review §1 D12-1…D12-14).
//!
//! Every test here reproduces one defect from the review with the smallest
//! fixture that showed it; each was verified red before its fix landed.
//! The expected values come from `bcftools 1.24` on the same fixture
//! **except** where 01-review §2.5-B marks a divergence as "not aligned"
//! (permanent, reed-stricter) — those assert reed's own documented value
//! and never bcftools parity (see 06-roadmap §1 gate 3 and guard 15).

use reed::{FilterOptions, run_filter};
use reed_filter::spec::{FilterSpec, SiteRules};

/// Compresses `text` into a BGZF stream cut at `block_data`-byte blocks.
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

/// Inflates a BGZF fixture back to text.
fn inflate_to_string(path: &std::path::Path) -> String {
    let bytes = std::fs::read(path).expect("invariant: read");
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut out = Vec::new();
    let mut offset = 0;
    while offset < bytes.len() {
        let block = &bytes[offset..];
        let Some(size) = reed_io::bgzf::block::parse_block_size(block) else {
            break;
        };
        if reed_io::bgzf::block::is_eof_block(&block[..size]) {
            break;
        }
        inflater
            .inflate_block(&block[..size], &mut out)
            .expect("invariant: inflates");
        offset += size;
    }
    String::from_utf8(out).expect("invariant: utf8")
}

/// Per-test scratch directory.
fn dir(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-info-r0-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    dir
}

fn options(strategy: Option<reed_filter::plan::ExecStrategy>, out: &std::path::Path) -> FilterOptions {
    FilterOptions {
        threads: 4,
        output: Some(out.to_path_buf()),
        write_index: reed::WriteIndex::None,
        compress_level: 6,
        strategy,
        command_line: "test".to_string(),
        timing: None,
        report_json: None,
        force_samples: false,
        regions: Vec::new(),
        tmp_dir: None,
        output_compressed: None,
        output_bcf: None,
        mem_bytes: None,
        no_version: false,
        no_command: false,
    }
}

/// Writes `text` as a BGZF fixture and returns its path.
fn fixture(dir: &std::path::Path, name: &str, text: &str) -> std::path::PathBuf {
    let path = dir.join(format!("{name}.vcf.gz"));
    std::fs::write(&path, bgzf(text, 256)).expect("invariant: write");
    path
}

/// Runs the filter and returns the surviving record lines (no header).
fn run_keep_lines(
    path: &std::path::Path,
    spec: &FilterSpec,
    strategy: Option<reed_filter::plan::ExecStrategy>,
    out_name: &str,
) -> Vec<String> {
    let out = path.parent().expect("invariant: parent").join(out_name);
    run_filter(Some(path), spec, &options(strategy, &out)).expect("invariant: filter runs");
    inflate_to_string(&out)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// Positions (column 2) of the surviving records.
fn kept_positions(lines: &[String]) -> Vec<u64> {
    lines
        .iter()
        .map(|l| l.split('\t').nth(1).expect("POS").parse().expect("pos"))
        .collect()
}

// ---------------------------------------------------------------------------
// R0-a — D12-1: sketch must not evaluate tail expressions without the record
// ---------------------------------------------------------------------------

const MIND_HEADER: &str = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=QD,Number=1,Type=Float,Description=\"qd\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
";

/// Four records, QD=200/200/5/5; the tail expression reads INFO/QD, so only
/// the first two survive. `--mind` must not change that (the sketch path has
/// no record bytes, so it must not be chosen while a tail expression exists).
#[test]
fn mind_with_tail_expression_keeps_records() {
    let d = dir("mind-tail-expr");
    let mut text = String::from(MIND_HEADER);
    for (i, qd) in [200, 200, 5, 5].iter().enumerate() {
        use std::fmt::Write as _;
        writeln!(
            text,
            "Chr01\t{}\t.\tA\tG\t50\tPASS\tQD={qd}\tGT\t0/1\t0/0\t0/0\t0/0",
            i + 1
        )
        .expect("invariant: string write");
    }
    let path = fixture(&d, "in", &text);

    // Reference: the single-pass answer (the plan's own strategy choice).
    let baseline = run_keep_lines(&path, &mind_tail_spec(), None, "single.vcf.gz");
    assert_eq!(
        kept_positions(&baseline),
        vec![1, 2],
        "single-pass baseline keeps the two QD=200 records"
    );

    // The defect: --mind re-picks the strategy after compile and used to
    // select sketch, whose mid stage has no record bytes -> every record
    // fails `QD>10` and the output is empty.
    let mut spec = mind_tail_spec();
    spec.mind = Some(reed_filter::spec::MissingRate::parse("0.9").expect("valid"));
    let got = run_keep_lines(&path, &spec, None, "mind.vcf.gz");
    assert_eq!(
        kept_positions(&got),
        vec![1, 2],
        "D12-1: --mind with a tail expression must keep the same records as single pass"
    );
    let _ = std::fs::remove_dir_all(&d);
}

/// The `-q` variant of D12-1: the strategy must not change the verdict.
/// INFO says AF=0.5 but the genotypes say 0, and `--mind` is in play, so
/// the mid stage (which has no INFO) and the single pass must agree on
/// which sites survive. Red before the fix: the single pass read the INFO
/// AF while the sketch fell back to the genotypes.
#[test]
fn mind_with_min_af_matches_single_pass() {
    let d = dir("mind-min-af");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
Chr01\t1\t.\tA\tG\t50\tPASS\tAC=4;AN=8\tGT\t0/0\t0/0\t0/0\t0/0
Chr01\t2\t.\tA\tG\t50\tPASS\tAC=4;AN=8\tGT\t0/0\t0/0\t0/0\t0/0
";
    let path = fixture(&d, "in", text);
    let mut spec = FilterSpec::default();
    spec.site_after_mask.min_af = Some(0.3);
    spec.mind = Some(reed_filter::spec::MissingRate::parse("0.9").expect("valid"));

    let auto = kept_positions(&run_keep_lines(&path, &spec, None, "auto.vcf.gz"));
    let single = kept_positions(&run_keep_lines(
        &path,
        &spec,
        Some(reed_filter::plan::ExecStrategy::SinglePass),
        "single.vcf.gz",
    ));
    let two_pass = kept_positions(&run_keep_lines(
        &path,
        &spec,
        Some(reed_filter::plan::ExecStrategy::TwoPassReparse),
        "two.vcf.gz",
    ));
    assert_eq!(
        auto, single,
        "D12-1: -q must not depend on the strategy (auto vs single pass)"
    );
    assert_eq!(
        two_pass, single,
        "D12-1: -q must not depend on the strategy (two pass vs single pass)"
    );
    let _ = std::fs::remove_dir_all(&d);
}

/// A `--mind` tail-expression spec (`F_MISSING<=0.5 && QD>10`).
fn mind_tail_spec() -> FilterSpec {
    let mut spec = FilterSpec::default();
    spec.site_after_mask.include = vec!["F_MISSING<=0.5 && QD>10".to_string()];
    spec
}

// ---------------------------------------------------------------------------
// R0-b — D12-2: GT-less records keep their AC/AN under a sample subset
// ---------------------------------------------------------------------------

/// FORMAT carries DP but no GT; the input declares AC=4;AN=6 for 3 samples.
/// Subsetting to two samples must not recompute AC/AN from nothing
/// (bcftools `view -s s0,s1` keeps `AC=4;AN=6`, vcfview.c:343).
#[test]
fn gt_less_records_keep_their_ac_an_under_subset() {
    let d = dir("gtless-subset");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
Chr01\t10\t.\tA\tG\t50\tPASS\tAC=4;AN=6\tDP\t9\t9\t9
";
    let path = fixture(&d, "in", text);
    let mut spec = FilterSpec::default();
    spec.sample_select = Some(reed_filter::spec::SampleSelect {
        keep: Some(reed_filter::spec::SampleList::Names(vec![
            "s0".into(),
            "s1".into(),
        ])),
        remove: None,
    });
    let lines = run_keep_lines(&path, &spec, None, "out.vcf.gz");
    assert_eq!(lines.len(), 1, "the record survives");
    let info = lines[0].split('\t').nth(7).expect("INFO");
    assert!(
        info.contains("AC=4") && info.contains("AN=6"),
        "D12-2: a GT-less record keeps its input AC/AN (bcftools view -s parity), got {info:?}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

// ---------------------------------------------------------------------------
// R0-c — D12-3/4/5: BCF value fidelity (integers, flags, strings)
// ---------------------------------------------------------------------------

/// Reads a BCF file back through bcftools (the htslib truth gate).
fn bcftools_view(path: &std::path::Path) -> Vec<String> {
    let out = std::process::Command::new("bcftools")
        .arg("view")
        .arg("-H")
        .arg(path)
        .output()
        .expect("invariant: bcftools runs");
    assert!(
        out.status.success(),
        "bcftools rejected {}: {}",
        path.display(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout)
        .expect("invariant: utf8")
        .lines()
        .map(str::to_string)
        .collect()
}

/// Runs the filter writing BCF and reads it back with bcftools.
fn run_bcf_roundtrip(path: &std::path::Path, out_name: &str) -> Vec<String> {
    let out = path.parent().expect("invariant: parent").join(out_name);
    let mut opts = options(None, &out);
    opts.output_bcf = Some(0); // -Ou
    run_filter(Some(path), &FilterSpec::default(), &opts).expect("invariant: filter runs");
    bcftools_view(&out)
}

/// Integer INFO beyond the BCF range must become missing (htslib
/// `vcf.c:4004-4013`; VCFv4.5 §1.3 forbids `[-2^31, -2^31+7]` and values
/// outside the 32-bit signed range in VCF **and** BCF).
#[test]
fn out_of_range_integers_become_missing_in_bcf() {
    let d = dir("int-range");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=XD,Number=1,Type=Integer,Description=\"x\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0
Chr01\t1\t.\tA\tG\t50\tPASS\tXD=3000000000\tGT\t0/1
Chr01\t2\t.\tA\tG\t50\tPASS\tXD=-3000000000\tGT\t0/1
Chr01\t3\t.\tA\tG\t50\tPASS\tXD=-2147483641\tGT\t0/1
Chr01\t4\t.\tA\tG\t50\tPASS\tXD=2147483647\tGT\t0/1
Chr01\t5\t.\tA\tG\t50\tPASS\tXD=-2147483640\tGT\t0/1
";
    let path = fixture(&d, "in", text);
    let records = run_bcf_roundtrip(&path, "out.bcf");
    let info: Vec<&str> = records
        .iter()
        .map(|l| l.split('\t').nth(7).expect("INFO"))
        .collect();
    assert_eq!(info[0], "XD=.", "3000000000 wraps to -1294967296 today");
    assert_eq!(info[1], "XD=.", "-3000000000 wraps to 1294967296 today");
    assert_eq!(
        info[2], "XD=.",
        "-2147483641 is in htslib's reserved band (BCF_MIN_BT_INT32)"
    );
    assert_eq!(info[3], "XD=2147483647", "the int32 max is in range");
    assert_eq!(info[4], "XD=-2147483640", "BCF_MIN_BT_INT32 is in range");
    let _ = std::fs::remove_dir_all(&d);
}

/// A Flag carrying a value is stored as a string (htslib vcf.c:3949-3951);
/// `DB=0` therefore survives as `DB=0` instead of flipping to `DB`.
#[test]
fn flag_with_a_value_round_trips_through_bcf() {
    let d = dir("flag-value");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"db\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0
Chr01\t1\t.\tA\tG\t50\tPASS\tDB=1\tGT\t0/1
Chr01\t2\t.\tA\tG\t50\tPASS\tDB=0\tGT\t0/1
Chr01\t3\t.\tA\tG\t50\tPASS\tDB\tGT\t0/1
";
    let path = fixture(&d, "in", text);
    let records = run_bcf_roundtrip(&path, "out.bcf");
    let info: Vec<&str> = records
        .iter()
        .map(|l| l.split('\t').nth(7).expect("INFO"))
        .collect();
    assert_eq!(info[0], "DB=1", "a Flag with a value keeps it (htslib)");
    assert_eq!(info[1], "DB=0", "DB=0 must not flip to DB");
    assert_eq!(info[2], "DB", "a bare flag stays a flag");
    let _ = std::fs::remove_dir_all(&d);
}

/// A String value of `.` is a literal dot, not a missing marker (the whole
/// INFO column is `.` when empty; a value `.` round-trips as `TAG=.`).
#[test]
fn string_dot_value_round_trips_through_bcf() {
    let d = dir("string-dot");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=SS,Number=1,Type=String,Description=\"s\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0
Chr01\t1\t.\tA\tG\t50\tPASS\tSS=.\tGT\t0/1
";
    let path = fixture(&d, "in", text);
    let records = run_bcf_roundtrip(&path, "out.bcf");
    let info = records[0].split('\t').nth(7).expect("INFO");
    assert_eq!(info, "SS=.", "a String value of `.` is literal (htslib)");
    let _ = std::fs::remove_dir_all(&d);
}

// ---------------------------------------------------------------------------
// R0-d — D12-6: -q/-Q fall back to genotypes once the cohort changed
// ---------------------------------------------------------------------------

/// INFO says AC=4;AN=8 (AF 0.5) but the kept samples are all hom-ref after
/// subsetting to s2,s3 (AF 0). bcftools recomputes from the subset GT and
/// drops both sites at `-q 0.1` (vcfview.c:341-350).
#[test]
fn min_af_after_subset_uses_the_genotypes() {
    let d = dir("min-af-subset");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
Chr01\t10\t.\tA\tG\t50\tPASS\tAC=4;AN=8\tGT\t1/1\t1/1\t0/0\t0/0
Chr01\t20\t.\tA\tG\t50\tPASS\tAC=4;AN=8\tGT\t1/1\t1/1\t0/0\t0/0
";
    let path = fixture(&d, "in", text);
    let mut spec = FilterSpec::default();
    spec.sample_select = Some(reed_filter::spec::SampleSelect {
        keep: Some(reed_filter::spec::SampleList::Names(vec![
            "s2".into(),
            "s3".into(),
        ])),
        remove: None,
    });
    spec.site_after_mask.min_af = Some(0.1);
    let lines = run_keep_lines(&path, &spec, None, "out.vcf.gz");
    assert!(
        lines.is_empty(),
        "D12-6: after subsetting the GT AF is 0, so -q 0.1 drops both sites; got {}",
        lines.len()
    );
    let _ = std::fs::remove_dir_all(&d);
}

// ---------------------------------------------------------------------------
// R0-e — D12-7: expression layer sees the header type
// ---------------------------------------------------------------------------

const EXPR_HEADER: &str = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"db\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">
##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"vec\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
Chr01\t1\t.\tA\tG\t50\tPASS\tDB;DP=10;DP4=0,0,30,40\tGT\t0/1\t0/0
Chr01\t2\t.\tA\tG\t50\tPASS\tDP=20;DP4=1,2,3,4\tGT\t0/1\t0/0
Chr01\t3\t.\tA\tG\t50\tPASS\tDB;DP=0;DP4=0,0,0,0\tGT\t0/1\t0/0
";

/// `DB=1` keeps the DB records (bcftools `filters_set_info_flag`; docs
/// "1 (or 0) to test the presence (or absence) of a flag").
#[test]
fn expression_flag_presence_is_one_zero() {
    let d = dir("expr-flag");
    let path = fixture(&d, "in", EXPR_HEADER);
    let mut spec = FilterSpec::default();
    spec.site.include = vec!["DB=1".to_string()];
    let kept = kept_positions(&run_keep_lines(&path, &spec, None, "pos.vcf.gz"));
    assert_eq!(kept, vec![1, 3], "DB=1 keeps the records carrying DB");

    spec.site.include = vec!["DB=0".to_string()];
    let kept = kept_positions(&run_keep_lines(&path, &spec, None, "neg.vcf.gz"));
    assert_eq!(kept, vec![2], "DB=0 keeps the records lacking DB");
    let _ = std::fs::remove_dir_all(&d);
}

/// A multi-valued INFO compares any-element (bcftools `CMP_VECTORS`).
#[test]
fn expression_vector_info_compares_any_element() {
    let d = dir("expr-vector");
    let path = fixture(&d, "in", EXPR_HEADER);
    let mut spec = FilterSpec::default();
    spec.site.include = vec!["DP4>10".to_string()];
    let kept = kept_positions(&run_keep_lines(&path, &spec, None, "out.vcf.gz"));
    assert_eq!(
        kept,
        vec![1],
        "DP4=0,0,30,40 has an element above 10 (position 2's max is 4)"
    );
    let _ = std::fs::remove_dir_all(&d);
}

/// A bare tag is always false in bcftools (`filter_test` never sets
/// `pass_site` for a lone value token); reed rejects it at compile time
/// (01-review §2.5-C1 default: an error beats a silent constant).
#[test]
fn expression_bare_tag_is_rejected_at_compile_time() {
    let d = dir("expr-bare");
    let path = fixture(&d, "in", EXPR_HEADER);
    let mut spec = FilterSpec::default();
    spec.site.include = vec!["DP".to_string()];
    let out = d.join("out.vcf.gz");
    let err = run_filter(Some(&path), &spec, &options(None, &out))
        .expect_err("a bare tag must not silently filter");
    let msg = err.to_string();
    assert!(
        msg.contains("bare") || msg.contains("always false"),
        "the error must explain the bcftools semantics, got {msg:?}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

// ---------------------------------------------------------------------------
// R0-f — D12-8: `quantifier = "any"` requires a satisfying element
// ---------------------------------------------------------------------------

/// `any` means "at least one ALT value satisfies the comparison"
/// (predicate.rs:168-176; architecture.md:1210).
#[test]
fn quantifier_any_requires_a_satisfying_element() {
    use reed_core::tag_table::{TagNumber, TagTable, TagType};
    use reed_filter::spec::{InfoPredicate, OnMissing, Quantifier};

    let mut tags = TagTable::new();
    tags.add_info("AC", TagType::Integer, TagNumber::A);
    let evaluate = |keep_if: &str, q: Quantifier, ac: &str| {
        let mut spec = FilterSpec::default();
        spec.policy.multiallelic = reed_filter::spec::MultiallelicPolicy::Allow;
        spec.site = SiteRules {
            info: vec![InfoPredicate {
                tag: "AC".into(),
                keep_if: keep_if.into(),
                on_missing: OnMissing::Pass,
                quantifier: Some(q),
            }],
            ..SiteRules::default()
        };
        let plan = reed_filter::plan::compile(&spec, &tags, "test", false).expect("compiles");
        let line = format!("Chr01\t10\t.\tA\tG,T\t50\tPASS\tAC={ac}\tGT\t0/1\t1/2\n");
        let lr = reed_io::vcf::line::LineRef::parse(line.as_bytes()).expect("parses");
        let mut hits = vec![0_u64; plan.head.len()];
        let mut scratch = reed_filter::eval::head::HeadScratch::new();
        reed_filter::eval::head::eval_head(&lr, 10, &plan.head, &mut hits, &mut scratch, &tags)
    };

    assert!(
        !evaluate(">= 10", Quantifier::Any, "1,2"),
        "no element of AC=1,2 satisfies >= 10"
    );
    assert!(
        !evaluate("<= 0", Quantifier::Any, "1,2"),
        "no element of AC=1,2 satisfies <= 0"
    );
    assert!(
        evaluate(">= 10", Quantifier::Any, "1,20"),
        "the second element satisfies >= 10"
    );
    assert!(
        !evaluate(">= 10", Quantifier::All, "1,20"),
        "the first element fails >= 10"
    );
    assert!(
        evaluate(">= 10", Quantifier::All, "11,20"),
        "both elements satisfy >= 10"
    );
}

// ---------------------------------------------------------------------------
// R0-g — D12-9: packed-count overflow errors instead of corrupting
// ---------------------------------------------------------------------------

/// 65 536 distinct INFO keys (arbitrary keys are permitted; §1.6.1) do not
/// fit the 16-bit n_info field; htslib refuses them at parse time
/// (vcf.c:3912) and reed must refuse them at write time instead of packing
/// the overflow into n_allele.
#[test]
fn info_count_overflow_errors_instead_of_corrupting() {
    let d = dir("info-overflow");
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=Chr01,length=100000>\n\
         ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
         ##INFO=<ID=K0,Number=1,Type=Integer,Description=\"k\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n",
    );
    let info: Vec<String> = (0..65_536).map(|i| format!("K0={i}")).collect();
    {
        use std::fmt::Write as _;
        writeln!(text, "Chr01\t1\t.\tA\tG\t50\tPASS\t{}\tGT\t0/1", info.join(";"))
            .expect("invariant: string write");
    }
    let path = fixture(&d, "in", &text);
    let out = d.join("out.bcf");
    let mut opts = options(None, &out);
    opts.output_bcf = Some(0); // -Ou
    let err = run_filter(Some(&path), &FilterSpec::default(), &opts)
        .expect_err("65 536 INFO fields cannot be packed into n_info");
    assert!(
        err.to_string().contains("INFO") || err.to_string().contains("too many"),
        "the error must name the limit, got {err:?}"
    );
    let _ = std::fs::remove_dir_all(&d);
}

// ---------------------------------------------------------------------------
// R0-h — D12-10: ALT='.' counts zero ALTs
// ---------------------------------------------------------------------------

/// `ALT=.` (no variant) is zero ALT alleles: `-m 2` drops it, `-M 1` keeps
/// only it, and `N_ALT=0` matches it (bcftools `filters_set_nalt`,
/// vcfview.c:418-434).
#[test]
fn invariant_sites_count_zero_alts() {
    let d = dir("alt-dot");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
Chr01\t10\t.\tA\t.\t50\tPASS\t.\tGT\t0/0\t0/0
Chr01\t20\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0
";
    let path = fixture(&d, "in", text);

    let mut spec = FilterSpec::default();
    spec.site.min_alt_alleles = 1; // -m 2
    let kept = kept_positions(&run_keep_lines(&path, &spec, None, "m2.vcf.gz"));
    assert_eq!(kept, vec![20], "-m 2 drops the invariant site");

    // `-M 1` is `max_alt_alleles = Some(0)`; the `Option` distinguishes it
    // from "no -M flag" (None).
    let mut spec = FilterSpec::default();
    spec.site.max_alt_alleles = Some(0);
    let kept = kept_positions(&run_keep_lines(&path, &spec, None, "M1.vcf.gz"));
    assert_eq!(kept, vec![10], "-M 1 keeps only the invariant site");

    let mut spec = FilterSpec::default();
    spec.site.include = vec!["N_ALT=0".to_string()];
    let kept = kept_positions(&run_keep_lines(&path, &spec, None, "n0.vcf.gz"));
    assert_eq!(kept, vec![10], "N_ALT=0 matches the invariant site");
    let _ = std::fs::remove_dir_all(&d);
}

// ---------------------------------------------------------------------------
// R0-i — D12-13: AC length mismatch falls back to the genotypes
// ---------------------------------------------------------------------------

/// AC declares one value for two ALTs: htslib warns and falls back to the
/// genotype counts (vcfutils.c:56-64), so `-q 0.1` drops the site (GT is
/// all hom-ref) instead of using the malformed 3/4 = 0.75.
#[test]
fn ac_length_mismatch_falls_back_to_genotypes() {
    let d = dir("ac-length");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=1000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
Chr01\t5\t.\tA\tG,T\t50\tPASS\tAC=3;AN=4\tGT\t0/0\t0/0
";
    let path = fixture(&d, "in", text);
    let mut spec = FilterSpec::default();
    spec.policy.multiallelic = reed_filter::spec::MultiallelicPolicy::Allow;
    spec.site_after_mask.min_af = Some(0.1);
    let lines = run_keep_lines(&path, &spec, None, "out.vcf.gz");
    assert!(
        lines.is_empty(),
        "D12-13: the malformed AC falls back to GT (AF=0) so -q 0.1 drops the site"
    );
    let _ = std::fs::remove_dir_all(&d);
}

// ---------------------------------------------------------------------------
// R0-j — D12-14: one value, one answer (internal consistency)
// ---------------------------------------------------------------------------

/// The same record must yield the same number to the filter, to `-Ob` and
/// to `-Ov`. `XI=1e3` used to filter as 1000, write `.` to BCF and keep
/// `1e3` in text. The assertions are reed's own spec (02 §13.3 reed column),
/// **not** bcftools parity: the strict-vs-prefix divergence is permanent
/// (01-review §2.5-B3).
#[test]
fn info_value_is_the_same_for_filtering_and_writing() {
    let d = dir("value-consistency");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=XI,Number=1,Type=Integer,Description=\"i\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0
Chr01\t1\tid1\tA\tG\t50\tPASS\tXI=5\tGT\t0/1
Chr01\t2\tid2\tA\tG\t50\tPASS\tXI=1e3\tGT\t0/1
Chr01\t3\tid3\tA\tG\t50\tPASS\tXI=5abc\tGT\t0/1
Chr01\t4\tid4\tA\tG\t50\tPASS\tXI= 5\tGT\t0/1
";
    let path = fixture(&d, "in", text);

    // Which records does `XI>4` keep?
    let mut spec = FilterSpec::default();
    spec.site.include = vec!["XI>4".to_string()];
    let kept = kept_positions(&run_keep_lines(&path, &spec, None, "filtered.vcf.gz"));

    // What does `-Ob` write for them? (`XI=1e3`/`5abc`/` 5` are malformed
    // per the strict rule, so they must read back as missing — the record
    // that passed the filter must therefore be the one that writes a
    // value satisfying the comparison.)
    let out = d.join("all.bcf");
    let mut opts = options(None, &out);
    opts.output_bcf = Some(0); // -Ou
    run_filter(Some(&path), &FilterSpec::default(), &opts).expect("invariant: runs");
    let records = bcftools_view(&out);
    let by_id = |id: &str| -> String {
        records
            .iter()
            .find(|l| l.split('\t').nth(2) == Some(id))
            .and_then(|l| l.split('\t').nth(7))
            .expect("record")
            .to_string()
    };

    for pos in &kept {
        let id = format!("id{pos}");
        let info = by_id(&id);
        assert!(
            info == "XI=." || info.starts_with("XI=5"),
            "D12-14: {id} passed `XI>4`, so its BCF value must satisfy the same comparison; got {info:?}"
        );
    }
    // The strict rule: `1e3`, `5abc` and ` 5` are not valid integers.
    assert!(!kept.contains(&2), "XI=1e3 is not a valid Integer");
    assert!(!kept.contains(&3), "XI=5abc is not a valid Integer");
    assert!(!kept.contains(&4), "XI= 5 has a leading space (no trim)");
    let _ = std::fs::remove_dir_all(&d);
}

// ---------------------------------------------------------------------------
// R2-c/d — D12-11/12: --fill-tags validation and the visible default -x
// ---------------------------------------------------------------------------

/// An unsupported `--fill-tags` tag is a spec error (bcftools `+fill-tags`
/// refuses it, plugins/fill-tags.c:530-535); it used to be a silent no-op
/// that the default `-x` then turned into a deletion of the input's tag.
#[test]
fn unsupported_fill_tags_are_rejected() {
    let d = dir("fill-tags");
    let path = fixture(&d, "in", EXPR_HEADER);
    for tag in ["AF", "BOGUS", "F_MISSING"] {
        let mut spec = FilterSpec::default();
        spec.recompute.tags = vec![tag.to_string()];
        let out = d.join("out.vcf.gz");
        let err = run_filter(Some(&path), &spec, &options(None, &out))
            .expect_err("an unsupported fill-tags entry must not be a silent no-op");
        let msg = err.to_string();
        assert!(
            msg.contains(tag) && msg.contains("supported"),
            "the error names the tag and the supported set, got {msg:?}"
        );
    }
    // Case-insensitive like bcftools' strcasecmp.
    let mut spec = FilterSpec::default();
    spec.recompute.tags = vec!["ac".to_string()];
    let out = d.join("lower.vcf.gz");
    let report = run_filter(Some(&path), &spec, &options(None, &out)).expect("lowercase AC works");
    assert_eq!(report.n_output_records, 3, "the lowercase spelling runs");
    let _ = std::fs::remove_dir_all(&d);
}

/// A GT-less record under the automatic subset update must not lose its
/// AC/AN even when the tag is also in the default `-x` list (the automatic
/// entries are exempt; an explicit `-x` is not).
#[test]
fn explicit_remove_still_deletes_on_gt_less_records() {
    let d = dir("gtless-remove");
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
Chr01\t10\t.\tA\tG\t50\tPASS\tAC=4;AN=6\tDP\t9\t9\t9
";
    let path = fixture(&d, "in", text);
    let mut spec = FilterSpec::default();
    spec.sample_select = Some(reed_filter::spec::SampleSelect {
        keep: Some(reed_filter::spec::SampleList::Names(vec![
            "s0".into(),
            "s1".into(),
        ])),
        remove: None,
    });
    // The user asked for the removal (the CLI sets `drop_stale_explicit`
    // when `-x` is given): it applies regardless of GT.
    spec.recompute.drop_stale = vec!["AC".into(), "AN".into()];
    spec.recompute.drop_stale_explicit = true;
    let lines = run_keep_lines(&path, &spec, None, "out.vcf.gz");
    let info = lines[0].split('\t').nth(7).expect("INFO");
    assert_eq!(info, ".", "an explicit -x removes the tags on GT-less records too");
    let _ = std::fs::remove_dir_all(&d);
}
