//! M5 acceptance tests: the three strategy/shape gates from
//! `docs/spec/architecture.md` §11.4 — SinglePass ≡ TwoPassReparse byte
//! equality, spec validation at compile time, and compat ≡ TOML.

use std::io::Write;

use reed::{FilterOptions, run_filter};
use reed_filter::FilterError;
use reed_filter::spec::{FilterSpec, OnMissing};
use std::fmt::Write as _;

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
        let gts = if i % 40 == 20 {
            ["0/0", "0/0", "0/0"]
        } else {
            gts
        };
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
        let _ = writeln!(
            text,
            "Chr01\t{pos}\t.\t{ref_b}\t{alt}\t{qual}\t{filter}\t{info}\tGT:DP:GQ\t{}\t{}\t{}\t{s3}",
            sample(gts[0], dps[0], gqs[0]),
            sample(gts[1], dps[1], gqs[1]),
            sample(gts[2], dps[2], gqs[2])
        );
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
fn gold_spec() -> FilterSpec {
    use reed_filter::spec::{
        GenotypeMaskRule, GenotypeRules, InfoPredicate, MissingRate, OnMissing, RecomputeSpec,
        SiteRules, TailRules,
    };
    FilterSpec {
        // s3 is missing in 98% of the records; s1's `--min-dp` mask leaves it
        // at exactly 0.1 (the boundary is kept).
        mind: Some(MissingRate::parse("0.1").expect("valid")),
        site: SiteRules {
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
            ..SiteRules::default()
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
            ..GenotypeRules::default()
        },
        site_after_mask: TailRules {
            max_missing: Some(MissingRate::parse("0.5").expect("valid")),
            min_maf: Some(0.05),
            ..TailRules::default()
        },
        recompute: RecomputeSpec::all(),
        ..FilterSpec::default()
    }
}

fn options(
    strategy: Option<reed_filter::plan::ExecStrategy>,
    out: &std::path::Path,
) -> FilterOptions {
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

#[test]
fn single_pass_and_two_pass_reparse_produce_identical_output() {
    let dir_path = write_fixture("strategies", 6_000);
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    let spec: FilterSpec = gold_spec();

    let out_single = dir.join("single.vcf.gz");
    let out_two = dir.join("two.vcf.gz");
    let r1 = run_filter(
        Some(dir_path.as_path()),
        &spec,
        &options(
            Some(reed_filter::plan::ExecStrategy::SinglePass),
            &out_single,
        ),
    )
    .expect("invariant: single pass runs");
    let r2 = run_filter(
        Some(dir_path.as_path()),
        &spec,
        &options(
            Some(reed_filter::plan::ExecStrategy::TwoPassReparse),
            &out_two,
        ),
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
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    let spec: FilterSpec = gold_spec();

    let mut reference: Option<Vec<u8>> = None;
    for threads in [1_usize, 3, 8] {
        let out = dir.join(format!("out-{threads}.vcf.gz"));
        run_filter(Some(dir_path.as_path()), &spec, &options(None, &out)).expect("invariant: run");
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
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();

    // Unknown INFO tag: rejected at compile time, before any record.
    let mut spec = gold_spec();
    spec.site.info[0].tag = "NotDeclared".to_string();
    let err = run_filter(
        Some(dir_path.as_path()),
        &spec,
        &options(None, &dir.join("o.vcf.gz")),
    )
    .expect_err("invariant: unknown tag rejected");
    assert!(err.to_string().contains("NotDeclared"), "{err}");

    // Threshold out of range.
    let mut spec = gold_spec();
    spec.site_after_mask.min_maf = Some(0.9);
    let err = run_filter(
        Some(dir_path.as_path()),
        &spec,
        &options(None, &dir.join("o.vcf.gz")),
    )
    .expect_err("invariant: maf range rejected");
    assert!(err.to_string().contains("min_maf"), "{err}");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn mask_default_on_missing_is_keep() {
    use reed_filter::spec::{GenotypeMaskRule, GenotypeRules};
    // The mask default is Keep at the spec level (§11.1): a missing DP/GQ
    // never masks (plink2 semantics), and the compiler accepts it.
    let spec = FilterSpec {
        genotype: GenotypeRules {
            mask_if: vec![GenotypeMaskRule {
                field: "DP".into(),
                when: "< 5".into(),
                on_missing: OnMissing::Keep,
            }],
            ..GenotypeRules::default()
        },
        ..FilterSpec::default()
    };
    assert_eq!(spec.genotype.mask_if[0].on_missing, OnMissing::Keep);
    let dir_path = write_fixture("maskdefault", 5);
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    let out = dir.join("o.vcf.gz");
    let result = run_filter(Some(dir_path.as_path()), &spec, &options(None, &out));
    assert!(result.is_ok(), "Pass normalizes to Keep: {result:?}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn provenance_and_ascertainment_are_written() {
    let dir_path = write_fixture("provenance", 30);
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    let spec: FilterSpec = gold_spec();
    let out = dir.join("o.vcf.gz");
    run_filter(Some(dir_path.as_path()), &spec, &options(None, &out)).expect("invariant: run");

    let text = {
        use std::io::Read;
        let mut reader =
            noodles::bgzf::io::Reader::new(std::fs::File::open(&out).expect("invariant: output"));
        let mut header = Vec::new();
        let mut buf = vec![0_u8; 1 << 16];
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
    // The manifest-v2 ascertainment line was withdrawn (doc_v2 04 §1,
    // doc_v3 C1): provenance carries the spec hash and command instead.
    assert!(!text.contains("##reed_ascertainment"), "{text}");
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
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    let spec: FilterSpec = gold_spec();

    let mut hashes = Vec::new();
    for (name, strategy) in [
        ("single", Some(reed_filter::plan::ExecStrategy::SinglePass)),
        ("two", Some(reed_filter::plan::ExecStrategy::TwoPassReparse)),
        (
            "sketch",
            Some(reed_filter::plan::ExecStrategy::SketchThenEmit),
        ),
    ] {
        let out = dir.join(format!("{name}.vcf.gz"));
        let report = run_filter(Some(dir_path.as_path()), &spec, &options(strategy, &out))
            .unwrap_or_else(|e| panic!("{name} failed: {e}"));
        assert_eq!(
            report.strategy,
            strategy.expect("invariant: forced").id(),
            "the report names the strategy"
        );
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
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    let spec: FilterSpec = gold_spec();

    let single = run_filter(
        Some(dir_path.as_path()),
        &spec,
        &options(
            Some(reed_filter::plan::ExecStrategy::SinglePass),
            &dir.join("s.vcf.gz"),
        ),
    )
    .expect("invariant: single pass");
    let sketch = run_filter(
        Some(dir_path.as_path()),
        &spec,
        &options(
            Some(reed_filter::plan::ExecStrategy::SketchThenEmit),
            &dir.join("k.vcf.gz"),
        ),
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
            let _ = writeln!(h, "##contig=<ID=c{c},length=1000000000>");
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

    fn opts(out: &std::path::Path, index: reed::WriteIndex) -> FilterOptions {
        FilterOptions {
            threads: 2,
            output: Some(out.to_path_buf()),
            write_index: index,
            compress_level: 6,
            strategy: None,
            command_line: "test".into(),
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

    #[test]
    fn csi_indexes_70000_contigs_and_positions_past_512mb() {
        let dir = std::env::temp_dir().join(format!("reed-filter-d4-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let input = wide_fixture(&dir);
        let out = dir.join("wide.out.vcf.gz");

        run_filter(
            Some(input.as_path()),
            &FilterSpec::default(),
            &opts(&out, reed::WriteIndex::Csi),
        )
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
            .copied()
            .collect();
        assert!(!chunks.is_empty(), "contig 69,999 has index chunks");
        use std::io::BufRead as _;
        let mut reader =
            noodles::bgzf::io::Reader::new(std::fs::File::open(&out).expect("invariant: output"));
        let query = noodles::csi::io::Query::new(&mut reader, chunks);
        let lines: Vec<Vec<u8>> = query
            .split(b'\n')
            .map(|l| l.expect("invariant: query read"))
            .collect();
        assert!(
            lines.iter().any(|l| l.starts_with(b"c69999\t600000000\t")),
            "record on c69999 at 600,000,000 must be found, got {lines:?}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn tbi_refuses_positions_beyond_29_bits() {
        let dir = std::env::temp_dir().join(format!("reed-filter-d4-tbi-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let input = wide_fixture(&dir);
        let out = dir.join("wide.tbi.vcf.gz");

        let result = run_filter(
            Some(input.as_path()),
            &FilterSpec::default(),
            &opts(&out, reed::WriteIndex::Tbi),
        );
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

        let result = run_filter(
            Some(input.as_path()),
            &FilterSpec::default(),
            &opts(&out, reed::WriteIndex::None),
        );
        assert!(
            result.is_err(),
            "an undeclared contig must fail, not alias tid 0"
        );
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
        let _ = write!(text, "\ts{s}");
    }
    text.push('\n');
    // Narrow record, then one ~105 KB record (30 digits per sample).
    text.push_str("c0\t10\t.\tA\tG\t.\tPASS\t.\tGT");
    for _s in 0..N {
        text.push_str("\t0/1");
    }
    text.push('\n');
    let mut wide = String::from("c0\t20\t.\tA\tG\t.\tPASS\t.\tGT:XH");
    for s in 0..N {
        let _ = s;
        let _ = write!(wide, "\t0/1:1.{}", "1234567890".repeat(3));
    }
    wide.push('\n');
    text.push_str(&wide);
    assert!(wide.len() > 100_000, "the wide record must exceed 100 KB");

    let input = dir.join("wide.vcf.gz");
    std::fs::write(&input, bgzf(&text, 4_096)).expect("invariant: write");
    let out = dir.join("wide.out.vcf.gz");

    run_filter(
        Some(input.as_path()),
        &FilterSpec::default(),
        &options(None, &out),
    )
    .expect("invariant: filter runs on a record wider than one block");

    use std::io::Read as _;
    let mut reader =
        noodles::bgzf::io::Reader::new(std::fs::File::open(&out).expect("invariant: output"));
    let mut round = Vec::new();
    reader
        .read_to_end(&mut round)
        .expect("invariant: output is a valid BGZF stream (bgzip -t equivalent)");
    let out_text = String::from_utf8(round).expect("invariant: utf8");
    let out_lines: Vec<&str> = out_text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(out_lines.len(), 2);
    assert_eq!(
        out_lines[1],
        wide.trim_end(),
        "wide record must round-trip byte-identical"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// R6 gate (03-architecture.md §6): the same data as plain VCF and as BGZF
/// VCF must produce identical filter output; stdin reads both through the
/// same entry.
#[test]
fn plain_and_bgzf_inputs_agree() {
    let dir = std::env::temp_dir().join(format!("reed-formats-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let text = fixture_text(400);

    let plain = dir.join("in.vcf");
    std::fs::write(&plain, &text).expect("invariant: write plain");
    let bgzf_path = dir.join("in.vcf.gz");
    std::fs::write(&bgzf_path, bgzf(&text, 512)).expect("invariant: write bgzf");

    let spec = FilterSpec::default();
    let out_plain = dir.join("plain.out.vcf.gz");
    let out_bgzf = dir.join("bgzf.out.vcf.gz");
    let rp =
        run_filter(Some(&plain), &spec, &options(None, &out_plain)).expect("invariant: plain run");
    let rb = run_filter(Some(&bgzf_path), &spec, &options(None, &out_bgzf))
        .expect("invariant: bgzf run");
    assert_eq!(rp.n_output_records, rb.n_output_records);
    assert_eq!(
        std::fs::read(&out_plain).expect("invariant: read"),
        std::fs::read(&out_bgzf).expect("invariant: read"),
        "plain and BGZF inputs must produce byte-identical output"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// R2 / D13: the declared decode policies execute. A called haploid is a
/// policy rejection under the default `require_diploid` (it used to be
/// silently truncated to a half call), and a polyploid genotype errors
/// instead of dropping its third allele.
#[test]
fn default_policies_reject_haploid_and_polyploid_calls() {
    let text = format!(
        "{HEADER}Chr01\t1\t.\tA\tG\t50\tPASS\tQD=10;FS=1;MQ=60;SOR=1\tGT:DP:GQ\t1:30:99\t0/1:30:99\t0/0:30:99\t0/0:30:99\n"
    );
    for (name, gt) in [("haploid", "1"), ("polyploid", "0/1/1")] {
        let dir = std::env::temp_dir().join(format!("reed-policy-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let path = dir.join("in.vcf.gz");
        let body = text.replace("\t1:30:99", &format!("\t{gt}:30:99"));
        std::fs::write(&path, bgzf(&body, 128)).expect("invariant: write");
        let out = dir.join("o.vcf.gz");
        let err = run_filter(Some(&path), &FilterSpec::default(), &options(None, &out))
            .expect_err("the default policies must reject this input");
        let message = err.to_string();
        assert!(
            message.contains("haploid") || message.contains("polyploid"),
            "{name}: unexpected error {message}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}

/// R2 / D13: under `allow_haploid` the transform is executed and counted —
/// the report's `policy_counters` become non-zero, and a multi-allelic
/// allele index collapses onto the single ALT bit (counted per allele).
#[test]
fn relaxed_policy_counts_half_calls_and_collapsed_alleles() {
    use reed_filter::spec::{PloidyPolicy, PolicyRules};
    let body = format!(
        "{HEADER}\
Chr01\t1\t.\tA\tG\t50\tPASS\tQD=10;FS=1;MQ=60;SOR=1\tGT:DP:GQ\t1:30:99\t0/. :30:99\t0/0:30:99\t0/0:30:99\n\
Chr01\t2\t.\tA\tG\t50\tPASS\tQD=10;FS=1;MQ=60;SOR=1\tGT:DP:GQ\t1/2:30:99\t0/0:30:99\t0/0:30:99\t0/0:30:99\n"
    )
    .replace("0/. :", "0/.:");
    let dir = std::env::temp_dir().join(format!("reed-policy-relaxed-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join("in.vcf.gz");
    std::fs::write(&path, bgzf(&body, 128)).expect("invariant: write");
    let out = dir.join("o.vcf.gz");
    let spec = FilterSpec {
        policy: PolicyRules {
            ploidy: PloidyPolicy::AllowHaploid,
            ..PolicyRules::default()
        },
        ..FilterSpec::default()
    };
    let report = run_filter(Some(&path), &spec, &options(None, &out)).expect("invariant: run");
    assert_eq!(report.n_output_records, 2);
    assert_eq!(report.policy_counters.haploid_as_half_call, 1);
    assert_eq!(report.policy_counters.multiallelic_collapsed, 1);
    // The policy governs decoding, not the record text: an untouched
    // genotype is written back byte for byte.
    let text = {
        use std::io::Read;
        let mut reader =
            noodles::bgzf::io::Reader::new(std::fs::File::open(&out).expect("invariant: output"));
        let mut s = String::new();
        reader.read_to_string(&mut s).expect("invariant: read");
        s
    };
    assert!(text.contains("\t1:30:99"), "haploid call preserved: {text}");
    assert!(
        text.contains("\t1/2:30:99"),
        "multi-allelic text preserved: {text}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// R5: the bounds added with the flag set — `-V/--exclude-types`,
/// `--max-maf`, `--max-mac` — filter the same sites as their mirrors.
#[test]
fn exclude_types_and_allele_bounds_are_mirrors() {
    use reed_filter::spec::SiteRules;
    let dir_path = write_fixture("bounds", 400);
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();

    let run = |spec: &FilterSpec, name: &str| {
        let out = dir.join(format!("{name}.vcf.gz"));
        let report = run_filter(Some(dir_path.as_path()), spec, &options(None, &out))
            .expect("invariant: run");
        let text = {
            use std::io::Read;
            let mut reader = noodles::bgzf::io::Reader::new(
                std::fs::File::open(&out).expect("invariant: output"),
            );
            let mut s = String::new();
            reader.read_to_string(&mut s).expect("invariant: read");
            s
        };
        (report.n_output_records, text)
    };

    let snv_only = FilterSpec {
        site: SiteRules {
            variant_class: vec!["snv".into()],
            ..SiteRules::default()
        },
        ..FilterSpec::default()
    };
    let exclude_indels = FilterSpec {
        site: SiteRules {
            exclude_variant_class: vec!["indel".into(), "mnv".into()],
            ..SiteRules::default()
        },
        ..FilterSpec::default()
    };
    let (n_snv, _) = run(&snv_only, "snv");
    let (n_excl, _) = run(&exclude_indels, "excl");
    assert_eq!(
        n_snv, n_excl,
        "-V indels,mnps must keep exactly the SNVs (fixture has no symbolic/ref)"
    );

    // MAF band: [0.1, 0.5] keeps exactly what max_maf leaves.
    let min_maf = FilterSpec {
        site_after_mask: reed_filter::spec::TailRules {
            min_maf: Some(0.1),
            ..reed_filter::spec::TailRules::default()
        },
        ..FilterSpec::default()
    };
    let with_max = FilterSpec {
        site_after_mask: reed_filter::spec::TailRules {
            min_maf: Some(0.1),
            max_maf: Some(0.5),
            ..reed_filter::spec::TailRules::default()
        },
        ..FilterSpec::default()
    };
    let (n_min, _) = run(&min_maf, "min");
    let (n_max, _) = run(&with_max, "max");
    assert_eq!(
        n_min, n_max,
        "MAF <= 0.5 is implied by the fixture's MAF distribution"
    );
    assert!(n_min < 400, "the MAF floor must drop something: {n_min}");

    // MAC upper bound mirrors the MAF one at 2N = 8 alleles.
    let max_mac = FilterSpec {
        site_after_mask: reed_filter::spec::TailRules {
            min_maf: Some(0.1),
            max_mac: 4,
            ..reed_filter::spec::TailRules::default()
        },
        ..FilterSpec::default()
    };
    let (n_mac, _) = run(&max_mac, "mac");
    assert!(
        n_mac <= n_max,
        "MAC <= 4 is a stricter band than MAF <= 0.5"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// R5: `--set-gts .` turns a `FMT/` expression into a per-sample mask
/// (bcftools `filter -S .`): `-i` masks where the expression is false
/// (missing values included), `-e` where it is true, and only GT changes.
#[test]
fn set_gts_masks_from_format_expressions() {
    use reed_filter::spec::{GenotypeRules, TailRules};
    let dir_path = write_fixture("setgts", 200);
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    // The fixture masks DP<5 / DP>100 / GQ<20 in the gold spec; here the
    // FMT expression does the masking.
    let run = |include: bool, name: &str| {
        let mut after = TailRules::default();
        if include {
            after.include.push("FMT/DP>10".to_string());
        } else {
            after.exclude.push("FMT/DP>10".to_string());
        }
        let spec = FilterSpec {
            genotype: GenotypeRules {
                set_gts: Some(".".to_string()),
                ..GenotypeRules::default()
            },
            site_after_mask: after,
            // The fixture carries multiallelic records; allow them so the
            // record count below stays about masking, not the default gate.
            policy: reed_filter::spec::PolicyRules {
                multiallelic: reed_filter::spec::MultiallelicPolicy::Allow,
                ..reed_filter::spec::PolicyRules::default()
            },
            ..FilterSpec::default()
        };
        let out = dir.join(format!("{name}.vcf.gz"));
        let report = run_filter(Some(dir_path.as_path()), &spec, &options(None, &out))
            .expect("invariant: run");
        let text = {
            use std::io::Read;
            let mut reader = noodles::bgzf::io::Reader::new(
                std::fs::File::open(&out).expect("invariant: output"),
            );
            let mut s = String::new();
            reader.read_to_string(&mut s).expect("invariant: read");
            s
        };
        (report.n_masked_gt, text)
    };

    let (masked_include, text_include) = run(true, "include");
    let (masked_exclude, text_exclude) = run(false, "exclude");
    assert!(masked_include > 0 && masked_exclude > 0);
    assert_ne!(
        masked_include, masked_exclude,
        "-i and -e mask different samples"
    );
    // `-i 'FMT/DP>10'` keeps samples above the threshold: the fixture's low-DP
    // samples (masked by the expression) become `./.`.
    assert!(
        text_include.contains("./.:3"),
        "fixture DP=3 must be masked under -i"
    );
    assert!(
        text_exclude.contains("./.:30") || text_exclude.contains("./.:150"),
        "DP=30/150 samples are masked under -e"
    );
    // No site is dropped: masking never removes a record.
    assert_eq!(
        text_include.lines().filter(|l| !l.starts_with('#')).count(),
        200
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// R5: `FMT/` expressions without `--set-gts .` are refused at compile time,
/// and subscripts are still reported as unsupported.
#[test]
fn format_expressions_need_set_gts() {
    use reed_filter::spec::TailRules;
    let dir_path = write_fixture("setgtsfail", 5);
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    let out = dir.join("o.vcf.gz");
    let spec = FilterSpec {
        site_after_mask: TailRules {
            include: vec!["FMT/DP>10".to_string()],
            ..TailRules::default()
        },
        ..FilterSpec::default()
    };
    let err = run_filter(Some(dir_path.as_path()), &spec, &options(None, &out))
        .expect_err("FMT/ without --set-gts must fail");
    assert!(err.to_string().contains("--set-gts"), "{err}");

    let spec = FilterSpec {
        genotype: reed_filter::spec::GenotypeRules {
            set_gts: Some(".".to_string()),
            ..reed_filter::spec::GenotypeRules::default()
        },
        site_after_mask: reed_filter::spec::TailRules {
            include: vec!["FMT/DP[0]>10".to_string()],
            ..reed_filter::spec::TailRules::default()
        },
        ..FilterSpec::default()
    };
    let err = run_filter(Some(dir_path.as_path()), &spec, &options(None, &out))
        .expect_err("subscripts are not supported yet");
    assert!(err.to_string().contains("subscripts"), "{err}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v4 02 回退版 §3.2 gate: `--mind` counts post-mask **missing alleles**
/// over the records that survive the head rules, and the denominator is
/// `2 × that record count`. Here s0 is missing in 90 records that fail the
/// head rule and once in the 10 that survive (0.1 — the head failures must
/// not drag it down); s1 is masked by `--min-dp` in half the survivors (0.5,
/// exactly the threshold → kept); s3 is missing in 9 of the 10 survivors
/// (0.9 → dropped).
#[test]
fn mind_counts_post_mask_over_head_survivors() {
    use reed_filter::plan::ExecStrategy;
    use reed_filter::spec::{GenotypeMaskRule, GenotypeRules, MissingRate, SiteRules};
    let mut text = String::from(HEADER);
    for i in 0..100_u32 {
        let pass = i >= 90;
        let filter = if pass { "PASS" } else { "LowQual" };
        let s0 = if !pass || i == 90 { "./." } else { "0/0" };
        let s1_dp = if pass && i < 95 { 3 } else { 30 };
        let s3 = if pass && i > 90 { "./." } else { "0/0" };
        let _ = writeln!(
            text,
            "Chr01\t{}\t.\tA\tG\t50\t{filter}\tQD=10;FS=1;MQ=60;SOR=1\tGT:DP:GQ\t\
             {s0}:30:99\t0/0:{s1_dp}:99\t0/0:30:99\t{s3}:30:99",
            i + 1
        );
    }
    let dir = std::env::temp_dir().join(format!("reed-sample-rate-c0-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join("in.vcf.gz");
    std::fs::write(&path, bgzf(&text, 256)).expect("invariant: write");
    let spec_with = |threshold: &str| FilterSpec {
        mind: Some(MissingRate::parse(threshold).expect("valid")),
        site: SiteRules {
            require_filter: vec!["PASS".into()],
            ..SiteRules::default()
        },
        genotype: GenotypeRules {
            mask_if: vec![GenotypeMaskRule {
                field: "DP".into(),
                when: "< 5".into(),
                on_missing: OnMissing::Keep,
            }],
            ..GenotypeRules::default()
        },
        ..FilterSpec::default()
    };
    for (name, strategy) in [
        ("single", ExecStrategy::SinglePass),
        ("two_pass", ExecStrategy::TwoPassReparse),
        ("sketch", ExecStrategy::SketchThenEmit),
    ] {
        let out = dir.join(format!("{name}.vcf.gz"));
        let report = run_filter(
            Some(&path),
            &spec_with("0.5"),
            &options(Some(strategy), &out),
        )
        .unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(report.n_samples_out, 3, "{name}: s3 dropped, s0/s1/s2 kept");
        // The rule reports under its stable id (doc_v4 02 回退版 §5.3), not
        // as a record funnel stage.
        let rule = report.sample_rule.as_ref().expect("sample rule reported");
        assert_eq!(rule.rule, "sample.mind", "{name}");
        assert_eq!(rule.threshold, "0.5", "{name}");
        assert_eq!(rule.n_dropped, 1, "{name}");
        let mut text = String::new();
        std::io::Read::read_to_string(
            &mut noodles::bgzf::io::Reader::new(
                std::fs::File::open(&out).expect("invariant: open"),
            ),
            &mut text,
        )
        .expect("invariant: read");
        let columns: Vec<&str> = text
            .lines()
            .find(|l| l.starts_with("#CHROM"))
            .expect("invariant: header")
            .split('\t')
            .collect();
        // s1's 10/20 missing alleles is exactly the threshold → kept; s3's
        // 18/20 → dropped; s0's 2/20 → kept (the 90 head failures are not
        // counted).
        assert!(!columns.contains(&"s3"), "{name}: s3 dropped: {columns:?}");
        assert!(
            columns.contains(&"s1"),
            "{name}: s1 kept at the boundary: {columns:?}"
        );
        assert!(columns.contains(&"s0"), "{name}: s0 kept: {columns:?}");

        // One step below s1's exact rate drops it too (missing > 0.5).
        let out = dir.join(format!("{name}-strict.vcf.gz"));
        let report = run_filter(
            Some(&path),
            &spec_with("0.499999999"),
            &options(Some(strategy), &out),
        )
        .unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(
            report.n_samples_out, 2,
            "{name}: s1 drops just below its rate"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// 04 §C0b: `-s ^A,B` excludes, `-S FILE` + the dynamic sample rule combine
/// (the static subset first, then `--mind`).
#[test]
fn static_subset_and_mind_combine() {
    use reed_filter::spec::{MissingRate, SampleList, SampleSelect};
    let dir_path = write_fixture("subsetmind", 1_000);
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    // Keep s0,s3 statically; drop s3 via the call-rate rule (98% missing).
    let spec = FilterSpec {
        sample_select: Some(SampleSelect {
            keep: Some(SampleList::Names(vec!["s0".into(), "s3".into()])),
            remove: None,
        }),
        mind: Some(MissingRate::parse("0.5").expect("valid")),
        ..FilterSpec::default()
    };
    let out = dir.join("o.vcf.gz");
    let report =
        run_filter(Some(dir_path.as_path()), &spec, &options(None, &out)).expect("invariant: run");
    assert_eq!(report.n_samples_in, 4);
    assert_eq!(
        report.n_samples_out, 1,
        "s3 dropped by --mind after the static subset"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// §11.6 / doc_v2 STATUS: the filter's default is plink2-style biallelic
/// input — a record with several ALTs fails the biallelic gate unless the
/// spec allows multiallelics (`multiallelic = "allow"` or
/// `max_alt_alleles > 1`).
#[test]
fn default_biallelic_gate_drops_multiallelic_records() {
    let mut text = String::from(HEADER);
    text.push_str(concat!(
        "Chr01\t1\t.\tA\tG\t50\tPASS\tQD=10;FS=1;MQ=60;SOR=1\tGT:DP:GQ\t",
        "0/1:30:99\t0/1:30:99\t0/1:30:99\t0/1:30:99\n"
    ));
    text.push_str(concat!(
        "Chr01\t2\t.\tA\tG,T\t50\tPASS\tQD=10;FS=1;MQ=60;SOR=1\tGT:DP:GQ\t",
        "1/2:30:99\t0/1:30:99\t0/1:30:99\t0/1:30:99\n"
    ));
    let dir = std::env::temp_dir().join(format!("reed-gate-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join("in.vcf.gz");
    std::fs::write(&path, bgzf(&text, 128)).expect("invariant: write");

    let out = dir.join("gated.vcf.gz");
    let gated = run_filter(
        Some(path.as_path()),
        &FilterSpec::default(),
        &options(None, &out),
    )
    .expect("invariant: run");
    assert_eq!(
        gated.n_output_records, 1,
        "the multiallelic record is gated out"
    );
    assert!(
        gated
            .funnel
            .iter()
            .any(|s| s.rule == "site.biallelic_gate" && s.n_passed == 1),
        "the gate is attributed in the funnel: {:?}",
        gated.funnel
    );

    let allowed = FilterSpec {
        policy: reed_filter::spec::PolicyRules {
            multiallelic: reed_filter::spec::MultiallelicPolicy::Allow,
            ..reed_filter::spec::PolicyRules::default()
        },
        ..FilterSpec::default()
    };
    let out = dir.join("allowed.vcf.gz");
    let report =
        run_filter(Some(path.as_path()), &allowed, &options(None, &out)).expect("invariant: run");
    assert_eq!(
        report.n_output_records, 2,
        "allow lets both records through"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// A7 / R7: with multiallelics allowed, the AC/MAF/HWE update needs per-ALT
/// counts; until R7 lands the driver refuses the combination instead of
/// writing one value where the header declares `Number=A`.
#[test]
fn ac_update_with_allowed_multiallelics_is_rejected() {
    use reed_filter::spec::{MultiallelicPolicy, PolicyRules, SampleList, SampleSelect};
    let dir_path = write_fixture("ac-multiallelic", 200);
    let dir = dir_path
        .parent()
        .expect("invariant: temp dir path")
        .to_path_buf();
    // Dropping a sample turns on the AC/AN update (bcftools view semantics).
    let spec = FilterSpec {
        sample_select: Some(SampleSelect {
            keep: Some(SampleList::Names(vec![
                "s0".into(),
                "s1".into(),
                "s2".into(),
            ])),
            remove: None,
        }),
        policy: PolicyRules {
            multiallelic: MultiallelicPolicy::Allow,
            ..PolicyRules::default()
        },
        ..FilterSpec::default()
    };
    let out = dir.join("err.vcf.gz");
    let err = run_filter(Some(dir_path.as_path()), &spec, &options(None, &out))
        .expect_err("multiallelic input + AC update must fail until R7");
    assert!(
        err.to_string().contains("R7"),
        "an actionable message: {err}"
    );

    // With the default gate the multiallelic records are dropped first, so
    // the same subset update runs and the recomputed AC stays valid.
    let gated = FilterSpec {
        policy: PolicyRules::default(),
        ..spec
    };
    let out = dir.join("ok.vcf.gz");
    let report = run_filter(Some(dir_path.as_path()), &gated, &options(None, &out))
        .expect("invariant: gated run");
    assert!(report.n_output_records > 0);
    assert_eq!(report.n_samples_out, 3);
    let _ = std::fs::remove_dir_all(&dir);
}

/// Y1e semantic lock: `--mind`'s denominator is the head-surviving record
/// count (doc_v4 02 回退版 §3.2), so a tail rule must not shift it. s1 is
/// missing only in records the tail rule keeps (8/20 = 0.4, kept at the
/// 0.5 boundary over the head survivors); over the post-tail records its
/// rate would be 0.8 and it would drop — the pre-Y1e prune matrices
/// collected post-tail and would flip it.
#[test]
fn mind_counts_ignore_tail_rules() {
    use reed_filter::spec::{MissingRate, TailRules};
    let mut text = String::from(HEADER);
    for i in 0..10_u32 {
        // Records 1-5: s2 and s3 missing → missing-allele fraction 4/8
        // fails `--geno 0.4`. Record 6: only s2 missing (tail passes).
        // Records 7-10: only s1 missing (s2's only other missing record).
        let row = |g: [&str; 4]| {
            format!(
                "Chr01\t{}\t.\tA\tG\t50\tPASS\t.\tGT\t{}\n",
                i + 1,
                g.join("\t")
            )
        };
        let gts = match i {
            0..=4 => ["0/1", "0/0", "./.", "./."],
            5 => ["0/0", "0/0", "./.", "0/0"],
            _ => ["0/0", "./.", "0/0", "0/0"],
        };
        text.push_str(&row(gts));
    }
    let dir = std::env::temp_dir().join(format!("reed-mind-tail-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join("in.vcf.gz");
    std::fs::write(&path, bgzf(&text, 256)).expect("invariant: write");
    let spec_with = |geno: bool| FilterSpec {
        mind: Some(MissingRate::parse("0.5").expect("valid")),
        site_after_mask: TailRules {
            max_missing: geno.then(|| MissingRate::parse("0.4").expect("valid")),
            ..TailRules::default()
        },
        ..FilterSpec::default()
    };
    for (name, geno) in [("plain", false), ("tail", true)] {
        let out = dir.join(format!("{name}.vcf.gz"));
        let report = run_filter(Some(&path), &spec_with(geno), &options(None, &out))
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let rule = report.sample_rule.as_ref().expect("sample rule reported");
        assert_eq!(rule.n_dropped, 1, "{name}: only s2 crosses 0.5");
        assert_eq!(report.n_samples_out, 3, "{name}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}
