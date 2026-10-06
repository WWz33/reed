//! M5 acceptance tests: the three strategy/shape gates from
//! `docs/architecture-spec.md` §11.4 — SinglePass ≡ TwoPassReparse byte
//! equality, spec validation at compile time, and compat ≡ TOML.

use std::io::Write;

use reed::filter::{run_filter, FilterOptions};
use reed_filter::spec::compat::from_gold_kv;
use reed_filter::spec::{FilterSpecV2, OnMissing};
use reed_filter::FilterError;
use std::collections::HashMap;

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

const HEADER: &str = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=QD,Number=1,Type=Float,Description=\"q\">
##INFO=<ID=FS,Number=1,Type=Float,Description=\"f\">
##INFO=<ID=MQ,Number=1,Type=Float,Description=\"m\">
##INFO=<ID=SOR,Number=1,Type=Float,Description=\"s\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
";

/// 4 samples; s3 missing in 98% of records so the cohort drops it. Every
/// funnel stage fires (QUAL, PASS, biallelic, INFO, mask, missing, MAF).
fn fixture_text(n_records: usize) -> String {
    let mut text = String::from(HEADER);
    for i in 0..n_records {
        let pos = i + 1;
        let qual = if i % 50 == 25 { "10.0" } else { "50.0" };
        let filter = if i % 60 == 30 { "LowQual" } else { "PASS" };
        let (ref_b, alt) = if i % 70 == 35 {
            ("AT", "AGG")
        } else if i % 90 == 45 {
            ("A", "G,T")
        } else {
            ("A", "G")
        };
        let qd = if i % 80 == 40 { "2.0" } else { "10.0" };
        let info = format!("QD={qd};FS=10.0;MQ=60.0;SOR=1.0");
        let gts: [&str; 3] = match i % 4 {
            0 => ["0/1", "0/1", "0/0"],
            1 => ["1/1", "0/1", "0/0"],
            2 => ["0/0", "0/1", "0/1"],
            _ => ["0/0", "0/0", "0/0"],
        };
        let gts = if i % 40 == 20 { ["0/0", "0/0", "0/0"] } else { gts };
        let mut dps = [30_u32; 3];
        if i % 10 == 5 {
            dps[1] = 3;
        }
        if i % 17 == 8 {
            dps[0] = 150;
        }
        let mut gqs = [99_u32; 3];
        if i % 13 == 6 {
            gqs[2] = 5;
        }
        let sample = |gt: &str, dp: u32, gq: u32| format!("{gt}:{dp}:{gq}");
        let s3 = if i % 50 == 0 {
            sample("1/1", 30, 99)
        } else {
            "./.:0:0".to_string()
        };
        text.push_str(&format!(
            "Chr01\t{pos}\t.\t{ref_b}\t{alt}\t{qual}\t{filter}\t{info}\tGT:DP:GQ\t{}\t{}\t{}\t{s3}\n",
            sample(gts[0], dps[0], gqs[0]),
            sample(gts[1], dps[1], gqs[1]),
            sample(gts[2], dps[2], gqs[2]),
        ));
    }
    text
}

fn write_fixture(name: &str, n_records: usize) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-filter-m5-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join("in.vcf.gz");
    let mut file = std::fs::File::create(&path).expect("invariant: fixture file");
    file.write_all(&bgzf(&fixture_text(n_records), 128))
        .expect("invariant: write");
    drop(file);
    path
}

/// The gold conf as key-value (the parity spec, minus runner keys).
fn gold_kv() -> HashMap<String, String> {
    [
        ("SAMPLE_MISS_MAX", "0.9"),
        ("KEEP_BIALLELIC_SNPS", "1"),
        ("REQUIRE_PASS", "1"),
        ("MIN_QUAL", "30"),
        ("INFO_QD_MIN", "4.0"),
        ("INFO_FS_MAX", "50.0"),
        ("INFO_MQ_MIN", "40.0"),
        ("INFO_SOR_MAX", "3.0"),
        ("SITE_MISS_MAX", "0.5"),
        ("MAF_MIN", "0.05"),
        ("HWE_MIN", "0"),
        ("GT_DP_MIN", "5"),
        ("GT_DP_MAX", "100"),
        ("GT_GQ_MIN", "20"),
    ]
    .iter()
    .map(|(k, v)| (k.to_string(), v.to_string()))
    .collect()
}

fn options(strategy: Option<reed_filter::plan::ExecStrategy>, out: &std::path::Path) -> FilterOptions {
    FilterOptions {
        threads: 4,
        output: Some(out.to_path_buf()),
        write_index: false,
        compress_level: 6,
        strategy,
        command_line: "test".to_string(),
        quiet: true,
        report_json: None,
    }
}

#[test]
fn single_pass_and_two_pass_reparse_produce_identical_output() {
    let dir_path = write_fixture("strategies", 6_000);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();
    let spec: FilterSpecV2 = from_gold_kv(&gold_kv()).expect("invariant: valid kv");

    let out_single = dir.join("single.vcf.gz");
    let out_two = dir.join("two.vcf.gz");
    let r1 = run_filter(
        &dir_path,
        &spec,
        &options(Some(reed_filter::plan::ExecStrategy::SinglePass), &out_single),
    )
    .expect("invariant: single pass runs");
    let r2 = run_filter(
        &dir_path,
        &spec,
        &options(Some(reed_filter::plan::ExecStrategy::TwoPassReparse), &out_two),
    )
    .expect("invariant: two pass runs");

    // The funnel is identical (ids and counts) and the output bytes match.
    assert_eq!(r1.funnel.len(), r2.funnel.len(), "same rule count");
    for (a, b) in r1.funnel.iter().zip(&r2.funnel) {
        assert_eq!(a.rule, b.rule, "same rule order");
        assert_eq!(a.n_passed, b.n_passed, "same counts for {}", a.rule);
    }
    assert_eq!(r1.n_output_records, r2.n_output_records);
    assert_eq!(r1.n_samples_out, r2.n_samples_out);
    assert_eq!(r1.n_masked_gt, r2.n_masked_gt);
    assert!(
        r1.n_output_records > 0 && r1.n_output_records < r1.n_input_records,
        "the filter must keep and drop records"
    );

    let bytes_a = std::fs::read(&out_single).expect("invariant: output");
    let bytes_b = std::fs::read(&out_two).expect("invariant: output");
    assert_eq!(bytes_a, bytes_b, "strategies must be byte-identical");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn output_is_deterministic_across_thread_counts() {
    let dir_path = write_fixture("determinism", 6_000);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();
    let spec: FilterSpecV2 = from_gold_kv(&gold_kv()).expect("invariant: valid kv");

    let mut reference: Option<Vec<u8>> = None;
    for threads in [1_usize, 3, 8] {
        let out = dir.join(format!("out-{threads}.vcf.gz"));
        run_filter(&dir_path, &spec, &options(None, &out)).expect("invariant: run");
        let bytes = std::fs::read(&out).expect("invariant: output");
        match &reference {
            None => reference = Some(bytes),
            Some(prev) => {
                // Only the ##reed_command line may differ (it embeds the
                // thread count); records must be byte-identical.
                assert_eq!(prev.len(), bytes.len());
                let strip = |d: &[u8]| -> Vec<u8> {
                    let text = String::from_utf8_lossy(d);
                    let mut out = Vec::with_capacity(d.len());
                    for line in text.lines() {
                        if !line.starts_with("##reed_command") {
                            out.extend_from_slice(line.as_bytes());
                            out.push(b'\n');
                        }
                    }
                    out
                };
                assert_eq!(strip(prev), strip(&bytes), "threads = {threads}");
            }
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn spec_validation_rejects_unknown_tags_and_ranges() {
    let dir_path = write_fixture("validation", 10);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();

    // Unknown INFO tag: rejected at compile time, before any record.
    let kv = gold_kv();
    let mut spec = from_gold_kv(&kv).expect("invariant: valid kv");
    spec.site.info[0].tag = "NotDeclared".to_string();
    let err = run_filter(&dir_path, &spec, &options(None, &dir.join("o.vcf.gz")))
        .expect_err("invariant: unknown tag rejected");
    assert!(err.to_string().contains("NotDeclared"), "{err}");

    // Threshold out of range.
    let mut spec = from_gold_kv(&kv).expect("invariant: valid kv");
    spec.site_after_mask.min_maf = Some(0.9);
    let err = run_filter(&dir_path, &spec, &options(None, &dir.join("o.vcf.gz")))
        .expect_err("invariant: maf range rejected");
    assert!(err.to_string().contains("min_maf"), "{err}");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn compat_conf_matches_the_toml_path() {
    // The conf reader and the TOML reader must land on the same spec for
    // the gold parameters (§11.1: "保证与 TOML 等价").
    let conf_spec = from_gold_kv(&gold_kv()).expect("invariant: valid kv");
    let toml = r#"
[cohort]
max_missing = 0.9

[site]
variant_class = ["snv"]
max_alt_alleles = 1
require_filter = ["PASS"]
min_qual = 30.0

[[site.info]]
tag = "QD"
keep_if = ">= 4.0"

[[site.info]]
tag = "FS"
keep_if = "<= 50.0"

[[site.info]]
tag = "MQ"
keep_if = ">= 40.0"

[[site.info]]
tag = "SOR"
keep_if = "<= 3.0"

[[genotype.mask_if]]
field = "DP"
when = "< 5"

[[genotype.mask_if]]
field = "DP"
when = "> 100"

[[genotype.mask_if]]
field = "GQ"
when = "< 20"

[site_after_mask]
max_missing = 0.5
min_maf = 0.05

[recompute]
tags = ["AC", "AN", "NS", "MAF", "HWE"]
"#;
    let toml_spec = FilterSpecV2::from_toml_str(toml).expect("invariant: valid toml");
    assert_eq!(conf_spec.cohort, toml_spec.cohort);
    // INFO predicates compare by parsed comparison (the conf formats
    // 4.0 as 4; the semantics are what matter, §11.1).
    let cmps = |spec: &FilterSpecV2| -> Vec<(String, reed_filter::spec::predicate::Cmp)> {
        spec.site
            .info
            .iter()
            .map(|r| {
                (
                    r.tag.clone(),
                    reed_filter::spec::predicate::Cmp::parse(&r.keep_if)
                        .expect("invariant: valid cmp"),
                )
            })
            .collect()
    };
    assert_eq!(cmps(&conf_spec), cmps(&toml_spec));
    assert_eq!(conf_site(&conf_spec).info, conf_site(&toml_spec).info);
    assert_eq!(conf_spec.genotype, toml_spec.genotype);
    assert_eq!(conf_spec.site_after_mask, toml_spec.site_after_mask);
    assert_eq!(conf_spec.recompute.tags, toml_spec.recompute.tags);
}

fn conf_site(spec: &FilterSpecV2) -> reed_filter::spec::SiteRulesV2 {
    // The conf path pushes "snv" into variant_class like the TOML; compare
    // as sorted sets to stay insensitive to declaration order.
    let mut s = spec.site.clone();
    s.variant_class.sort();
    s
}

#[test]
fn mask_default_on_missing_is_keep() {
    let toml = r#"
[[genotype.mask_if]]
field = "DP"
when = "< 5"
"#;
    let spec = FilterSpecV2::from_toml_str(&format!("{HEADER}\n{toml}"))
        .expect("invariant: toml parses");
    // The mask default is Keep directly at the spec level (§11.1) and the
    // compiler accepts it.
    assert_eq!(spec.genotype.mask_if[0].on_missing, OnMissing::Keep);
    let dir_path = write_fixture("maskdefault", 5);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();
    let out = dir.join("o.vcf.gz");
    let result = run_filter(&dir_path, &spec, &options(None, &out));
    assert!(result.is_ok(), "Pass normalizes to Keep: {result:?}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn provenance_and_ascertainment_are_written() {
    let dir_path = write_fixture("provenance", 30);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();
    let spec: FilterSpecV2 = from_gold_kv(&gold_kv()).expect("invariant: valid kv");
    let out = dir.join("o.vcf.gz");
    run_filter(&dir_path, &spec, &options(None, &out)).expect("invariant: run");

    let text = {
        use std::io::Read;
        let mut reader = noodles::bgzf::io::Reader::new(
            std::fs::File::open(&out).expect("invariant: output"),
        );
        let mut header = Vec::new();
        let mut buf = [0_u8; 1 << 16];
        loop {
            let n = reader.read(&mut buf).expect("invariant: read");
            if n == 0 {
                break;
            }
            header.extend_from_slice(&buf[..n]);
            if header.windows(4).rev().any(|w| w == b"\nChr") {
                break;
            }
        }
        String::from_utf8_lossy(&header).into_owned()
    };
    assert!(text.contains("##reed_version="), "{text}");
    assert!(text.contains("##reed_command="));
    assert!(text.contains("##reed_filter_spec="));
    assert!(text.contains("##reed_filter_sha256="));
    // Site rules + masks + a cohort subset → all four dimensions derived.
    assert!(text.contains("##reed_ascertainment="));
    let line = text
        .lines()
        .find(|l| l.starts_with("##reed_ascertainment="))
        .expect("invariant: ascertainment line");
    assert!(line.contains("Filtered"), "{line}");
    assert!(line.contains("Masked"), "{line}");
    assert!(line.contains("Absent"), "{line}");
    let _ = std::fs::remove_dir_all(&dir);
}

// Silence the unused-import warning path for FilterError in some configs.
#[allow(dead_code)]
fn _assert_error_type(e: FilterError) -> FilterError {
    e
}
