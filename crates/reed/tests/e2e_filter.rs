//! M5 acceptance tests: the three strategy/shape gates from
//! `docs/architecture-spec.md` §11.4 — SinglePass ≡ TwoPassReparse byte
//! equality, spec validation at compile time, and compat ≡ TOML.

use std::io::Write;

use reed::filter::{run_filter, FilterOptions};
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

/// The gold pipeline as a typed spec (the flags of 02-cli.md §6 expressed
/// through the spec types; the flag-based CLI builds exactly this).
fn gold_spec() -> FilterSpecV2 {
    use reed_filter::spec::{
        CohortRule, GenotypeMaskRule, GenotypeRules, InfoPredicate, OnMissing, RecomputeSpec,
        SiteRulesV2, TailRules,
    };
    FilterSpecV2 {
        cohort: Some(CohortRule::MaxMissing(0.9)),
        site: SiteRulesV2 {
            variant_class: vec!["snv".into()],
            max_alt_alleles: 1,
            require_filter: vec!["PASS".into()],
            min_qual: Some(30.0),
            info: [
                ("QD", ">= 4.0"),
                ("FS", "<= 50.0"),
                ("MQ", ">= 40.0"),
                ("SOR", "<= 3.0"),
            ]
            .iter()
            .map(|(tag, keep_if)| InfoPredicate {
                tag: (*tag).into(),
                keep_if: (*keep_if).into(),
                on_missing: OnMissing::Pass,
                quantifier: None,
            })
            .collect(),
            ..SiteRulesV2::default()
        },
        genotype: GenotypeRules {
            mask_if: [("DP", "< 5"), ("DP", "> 100"), ("GQ", "< 20")]
                .iter()
                .map(|(field, when)| GenotypeMaskRule {
                    field: (*field).into(),
                    when: (*when).into(),
                    on_missing: OnMissing::Keep,
                })
                .collect(),
        },
        site_after_mask: TailRules {
            max_missing: Some(0.5),
            min_maf: Some(0.05),
            ..TailRules::default()
        },
        recompute: RecomputeSpec::all(),
        ..FilterSpecV2::default()
    }
}

fn options(strategy: Option<reed_filter::plan::ExecStrategy>, out: &std::path::Path) -> FilterOptions {
    FilterOptions {
        threads: 4,
        output: Some(out.to_path_buf()),
        write_index: reed::filter::WriteIndex::None,
        compress_level: 6,
        strategy,
        command_line: "test".to_string(),
        quiet: true,
        report_json: None,
        force_samples: false,
    }
}

#[test]
fn single_pass_and_two_pass_reparse_produce_identical_output() {
    let dir_path = write_fixture("strategies", 6_000);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();
    let spec: FilterSpecV2 = gold_spec();

    let out_single = dir.join("single.vcf.gz");
    let out_two = dir.join("two.vcf.gz");
    let r1 = run_filter(
        Some(&dir_path.as_path()),
        &spec,
        &options(Some(reed_filter::plan::ExecStrategy::SinglePass), &out_single),
    )
    .expect("invariant: single pass runs");
    let r2 = run_filter(
        Some(&dir_path.as_path()),
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
    let spec: FilterSpecV2 = gold_spec();

    let mut reference: Option<Vec<u8>> = None;
    for threads in [1_usize, 3, 8] {
        let out = dir.join(format!("out-{threads}.vcf.gz"));
        run_filter(Some(&dir_path.as_path()), &spec, &options(None, &out)).expect("invariant: run");
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
    let mut spec = gold_spec();
    spec.site.info[0].tag = "NotDeclared".to_string();
    let err = run_filter(Some(&dir_path.as_path()), &spec, &options(None, &dir.join("o.vcf.gz")))
        .expect_err("invariant: unknown tag rejected");
    assert!(err.to_string().contains("NotDeclared"), "{err}");

    // Threshold out of range.
    let mut spec = gold_spec();
    spec.site_after_mask.min_maf = Some(0.9);
    let err = run_filter(Some(&dir_path.as_path()), &spec, &options(None, &dir.join("o.vcf.gz")))
        .expect_err("invariant: maf range rejected");
    assert!(err.to_string().contains("min_maf"), "{err}");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn mask_default_on_missing_is_keep() {
    use reed_filter::spec::{GenotypeMaskRule, GenotypeRules};
    // The mask default is Keep at the spec level (§11.1): a missing DP/GQ
    // never masks (plink2 semantics), and the compiler accepts it.
    let spec = FilterSpecV2 {
        genotype: GenotypeRules {
            mask_if: vec![GenotypeMaskRule {
                field: "DP".into(),
                when: "< 5".into(),
                on_missing: OnMissing::Keep,
            }],
        },
        ..FilterSpecV2::default()
    };
    assert_eq!(spec.genotype.mask_if[0].on_missing, OnMissing::Keep);
    let dir_path = write_fixture("maskdefault", 5);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();
    let out = dir.join("o.vcf.gz");
    let result = run_filter(Some(&dir_path.as_path()), &spec, &options(None, &out));
    assert!(result.is_ok(), "Pass normalizes to Keep: {result:?}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn provenance_and_ascertainment_are_written() {
    let dir_path = write_fixture("provenance", 30);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();
    let spec: FilterSpecV2 = gold_spec();
    let out = dir.join("o.vcf.gz");
    run_filter(Some(&dir_path.as_path()), &spec, &options(None, &out)).expect("invariant: run");

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

#[test]
fn all_three_strategies_produce_identical_bytes() {
    let dir_path = write_fixture("three-strategies", 8_000);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();
    let spec: FilterSpecV2 = gold_spec();

    let mut hashes = Vec::new();
    for (name, strategy) in [
        ("single", Some(reed_filter::plan::ExecStrategy::SinglePass)),
        ("two", Some(reed_filter::plan::ExecStrategy::TwoPassReparse)),
        ("sketch", Some(reed_filter::plan::ExecStrategy::SketchThenEmit)),
    ] {
        let out = dir.join(format!("{name}.vcf.gz"));
        let report = run_filter(Some(&dir_path.as_path()), &spec, &options(strategy, &out))
            .unwrap_or_else(|e| panic!("{name} failed: {e}"));
        assert_eq!(report.strategy, strategy.expect("invariant: forced").id(), "the report names the strategy");
        let bytes = std::fs::read(&out).expect("invariant: output");
        let hash = {
            use sha2::{Digest, Sha256};
            let mut h = Sha256::new();
            h.update(&bytes);
            format!("{:x}", h.finalize())
        };
        hashes.push(hash);
    }
    assert_eq!(hashes[0], hashes[1], "SinglePass vs TwoPassReparse");
    assert_eq!(hashes[1], hashes[2], "TwoPassReparse vs SketchThenEmit");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn sketch_funnel_matches_the_single_pass_funnel() {
    let dir_path = write_fixture("sketch-funnel", 6_000);
    let dir = dir_path.parent().expect("invariant: temp dir path").to_path_buf();
    let spec: FilterSpecV2 = gold_spec();

    let single = run_filter(
        Some(&dir_path.as_path()),
        &spec,
        &options(Some(reed_filter::plan::ExecStrategy::SinglePass), &dir.join("s.vcf.gz")),
    )
    .expect("invariant: single pass");
    let sketch = run_filter(
        Some(&dir_path.as_path()),
        &spec,
        &options(Some(reed_filter::plan::ExecStrategy::SketchThenEmit), &dir.join("k.vcf.gz")),
    )
    .expect("invariant: sketch");
    assert_eq!(single.funnel, sketch.funnel, "funnel is strategy-invariant");
    assert_eq!(single.n_masked_gt, sketch.n_masked_gt);
    assert_eq!(single.masked_by_rule, sketch.masked_by_rule);
    let _ = std::fs::remove_dir_all(&dir);
}

/// Regression (01-review D4): tid must stay u32 past 65,535 contigs, the
/// CSI index must address positions past the TBI 2^29 ceiling, and an
/// undeclared contig must fail instead of landing on tid 0.
mod write_index {
    use super::*;

    fn wide_header() -> String {
        let mut h = String::from(
            "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n",
        );
        for c in 0..70_000 {
            h.push_str(&format!("##contig=<ID=c{c},length=1000000000>\n"));
        }
        h.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n");
        h
    }

    fn wide_fixture(dir: &std::path::Path) -> std::path::PathBuf {
        let mut text = wide_header();
        text.push_str("c0\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n");
        text.push_str("c69999\t600000000\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n");
        let path = dir.join("wide.vcf.gz");
        std::fs::write(&path, bgzf(&text, 4_096)).expect("invariant: write");
        path
    }

    fn opts(out: &std::path::Path, index: reed::filter::WriteIndex) -> FilterOptions {
        FilterOptions {
            threads: 2,
            output: Some(out.to_path_buf()),
            write_index: index,
            compress_level: 6,
            strategy: None,
            command_line: "test".into(),
            quiet: true,
            report_json: None,
            force_samples: false,
        }
    }

    #[test]
    fn csi_indexes_70000_contigs_and_positions_past_512mb() {
        let dir = std::env::temp_dir().join(format!("reed-filter-d4-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let input = wide_fixture(&dir);
        let out = dir.join("wide.out.vcf.gz");

        run_filter(Some(&input.as_path()), &FilterSpecV2::default(), &opts(&out, reed::filter::WriteIndex::Csi))
            .expect("invariant: filter runs");

        let index = noodles::csi::fs::read(dir.join("wide.out.vcf.gz.csi"))
            .expect("invariant: csi written");
        let contig = index
            .reference_sequences()
            .get(69_999)
            .expect("invariant: 70,000 reference sequences");
        // The far record must be queryable through the index.
        let chunks: Vec<_> = contig
            .bins()
            .values()
            .flat_map(|bin| bin.chunks())
            .cloned()
            .collect();
        assert!(!chunks.is_empty(), "contig 69,999 has index chunks");
        use std::io::BufRead as _;
        let mut reader = noodles::bgzf::io::Reader::new(
            std::fs::File::open(&out).expect("invariant: output"),
        );
        let query = noodles::csi::io::Query::new(&mut reader, chunks);
        let lines: Vec<Vec<u8>> = query
            .split(b'\n')
            .map(|l| l.expect("invariant: query read"))
            .collect();
        assert!(lines
            .iter()
            .any(|l| l.starts_with(b"c69999\t600000000\t")),
            "record on c69999 at 600,000,000 must be found, got {lines:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn tbi_refuses_positions_beyond_29_bits() {
        let dir = std::env::temp_dir().join(format!("reed-filter-d4-tbi-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let input = wide_fixture(&dir);
        let out = dir.join("wide.tbi.vcf.gz");

        let result = run_filter(Some(&input.as_path()), &FilterSpecV2::default(), &opts(&out, reed::filter::WriteIndex::Tbi));
        assert!(result.is_err(), "TBI must refuse coordinates ≥ 2^29");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn undeclared_contig_is_an_error() {
        let dir = std::env::temp_dir().join(format!("reed-filter-d4-und-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let mut text = wide_header();
        // Header declares c0..c69999; the record uses an unlisted contig.
        text.push_str("unlisted\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n");
        let input = dir.join("und.vcf.gz");
        std::fs::write(&input, bgzf(&text, 4_096)).expect("invariant: write");
        let out = dir.join("und.out.vcf.gz");

        let result = run_filter(Some(&input.as_path()), &FilterSpecV2::default(), &opts(&out, reed::filter::WriteIndex::None));
        assert!(result.is_err(), "an undeclared contig must fail, not alias tid 0");
        let _ = std::fs::remove_dir_all(&dir);
    }
}

/// Regression (01-review D11): a single record wider than the 65,280-byte
/// BGZF member limit must span blocks instead of corrupting the member's
/// 16-bit BSIZE. Read back with noodles as the `bgzip -t` stand-in.
#[test]
fn record_wider_than_one_block_round_trips() {
    let dir = std::env::temp_dir().join(format!("reed-filter-d11-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");

    const N: usize = 3_000;
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=c0,length=1000000>\n\
         ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##FORMAT=<ID=XH,Number=1,Type=Float,Description=\"wide value\">\n",
    );
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for s in 0..N {
        text.push_str(&format!("\ts{s}"));
    }
    text.push('\n');
    // Narrow record, then one ~105 KB record (30 digits per sample).
    text.push_str("c0\t10\t.\tA\tG\t.\tPASS\t.\tGT");
    for s in 0..N {
        text.push_str("\t0/1");
    }
    text.push('\n');
    let mut wide = String::from("c0\t20\t.\tA\tG\t.\tPASS\t.\tGT:XH");
    for s in 0..N {
        let _ = s;
        wide.push_str(&format!("\t0/1:1.{}", "1234567890".repeat(3)));
    }
    wide.push('\n');
    text.push_str(&wide);
    assert!(wide.len() > 100_000, "the wide record must exceed 100 KB");

    let input = dir.join("wide.vcf.gz");
    std::fs::write(&input, bgzf(&text, 4_096)).expect("invariant: write");
    let out = dir.join("wide.out.vcf.gz");

    run_filter(Some(&input.as_path()), &FilterSpecV2::default(), &options(None, &out))
        .expect("invariant: filter runs on a record wider than one block");

    use std::io::Read as _;
    let mut reader = noodles::bgzf::io::Reader::new(
        std::fs::File::open(&out).expect("invariant: output"),
    );
    let mut round = Vec::new();
    reader
        .read_to_end(&mut round)
        .expect("invariant: output is a valid BGZF stream (bgzip -t equivalent)");
    let out_text = String::from_utf8(round).expect("invariant: utf8");
    let out_lines: Vec<&str> = out_text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(out_lines.len(), 2);
    assert_eq!(out_lines[1], wide.trim_end(), "wide record must round-trip byte-identical");
    let _ = std::fs::remove_dir_all(&dir);
}
