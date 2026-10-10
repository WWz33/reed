//! BCF input gates (doc_v3 02-bcf.md B1/B3a): the same file in VCF.gz and
//! BCF must produce identical statistics, report files and filter output.

use std::path::{Path, PathBuf};

use reed::{FilterOptions, run_filter};
use reed::{Levels, StatsRequest, WindowLayoutOptions, run_stats};
use reed_engine::budget::{MemBudget, ThreadBudget};
use reed_filter::spec::FilterSpec;

/// The fixtures live with the decoder (reed-io); the drivers read them from
/// here so one binary fixture serves both test suites.
fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../reed-io/tests/fixtures/bcf")
        .join(name)
}

/// A per-test scratch directory (the tests run in parallel).
fn dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-bcf-e2e-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    dir
}

/// Inflates a whole BGZF buffer into bytes.
fn inflate_bgzf_bytes(bytes: &[u8]) -> Vec<u8> {
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut out = Vec::new();
    let mut rest = bytes;
    while let Some(total) = reed_io::bgzf::block::parse_block_size(rest) {
        inflater
            .inflate_block(&rest[..total], &mut out)
            .expect("invariant: inflates");
        rest = &rest[total..];
    }
    out
}

/// Re-blocks a byte stream as BGZF with `block`-byte payloads.
fn reblock(bytes: &[u8], block: usize) -> Vec<u8> {
    let mut deflater = reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut out = Vec::new();
    for piece in bytes.chunks(block) {
        deflater
            .deflate_block(piece, &mut out)
            .expect("invariant: block fits");
    }
    out.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
    out
}

/// Inflates a whole BGZF buffer (the filter writes BGZF text).
fn inflate_bgzf(bytes: &[u8]) -> String {
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut text = Vec::new();
    let mut rest = bytes;
    while let Some(total) = reed_io::bgzf::block::parse_block_size(rest) {
        inflater
            .inflate_block(&rest[..total], &mut text)
            .expect("invariant: output inflates");
        rest = &rest[total..];
    }
    String::from_utf8(text).expect("invariant: VCF text")
}

fn budget() -> ThreadBudget {
    ThreadBudget {
        threads: 4,
        mem: MemBudget::from_bytes(u64::MAX),
    }
}

fn request() -> StatsRequest {
    StatsRequest {
        kinds: reed::window_kinds(2, true, true, true, false, None, None),
        gt_groups: vec![Some(0), Some(0), Some(1)],
        window_size: 1000,
        step: 1000,
        levels: Levels {
            window: true,
            contig: false,
            genome: false,
        },
        skip_empty: true,
        contig_lengths: std::collections::HashMap::new(),
        layout: Some(WindowLayoutOptions {
            invariant: reed::InvariantMode::VariantOnly,
        }),
        ..StatsRequest::default()
    }
}

/// `stats` over the BCF fixture equals the same run over the VCF.gz fixture.
#[test]
fn stats_output_is_identical_across_formats() {
    let bcf = run_stats(Some(&fixture("fix.bcf")), request(), &budget()).expect("invariant: bcf");
    let vcf =
        run_stats(Some(&fixture("fix.vcf.gz")), request(), &budget()).expect("invariant: vcf.gz");
    assert_eq!(bcf.tsv, vcf.tsv);
    assert_eq!(bcf.n_records, vcf.n_records);
}

/// Reports read FORMAT values: the BCF renderer feeds the same scanner.
#[test]
fn reports_are_identical_across_formats() {
    let dir = dir("reports");
    let bcf_prefix = dir.join("b");
    let vcf_prefix = dir.join("v");
    let request = |prefix: &Path| reed::reports::ReportRequest {
        freq: true,
        hardy: true,
        het: true,
        missing: true,
        depth: true,
        summary: false,
        targets: Vec::new(),
        targets_exclude: Vec::new(),
        prefix: prefix.to_path_buf(),
    };
    reed::reports::run_reports(Some(&fixture("fix.bcf")), &request(&bcf_prefix), &budget())
        .expect("invariant: bcf reports");
    reed::reports::run_reports(
        Some(&fixture("fix.vcf.gz")),
        &request(&vcf_prefix),
        &budget(),
    )
    .expect("invariant: vcf reports");
    for suffix in [
        ".afreq",
        ".hardy",
        ".het",
        ".smiss",
        ".vmiss",
        ".idepth",
        ".ldepth.mean",
    ] {
        let read = |prefix: &Path| {
            let mut p = prefix.as_os_str().to_os_string();
            p.push(suffix);
            std::fs::read(PathBuf::from(p)).expect("invariant: report file")
        };
        assert_eq!(read(&bcf_prefix), read(&vcf_prefix), "report {suffix}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// `-Ob` output is BCF that reed's own decoder renders back to the exact
/// records the `-Oz` path wrote (the round trip exercises encode + decode
/// together, since tests cannot shell out to bcftools).
#[test]
fn bcf_output_round_trips_to_the_same_records() {
    let dir = dir("bcf-out");
    let spec = FilterSpec {
        site: reed_filter::spec::SiteRules {
            min_qual: Some(50.0),
            ..reed_filter::spec::SiteRules::default()
        },
        // The fixture carries a haploid call (chr2:500).
        policy: reed_filter::spec::PolicyRules {
            ploidy: reed_filter::spec::PloidyPolicy::AllowHaploid,
            ..reed_filter::spec::PolicyRules::default()
        },
        ..FilterSpec::default()
    };
    let run = |input: &Path, out: &Path, bcf: Option<u8>| {
        let opts = FilterOptions {
            threads: 2,
            output: Some(out.to_path_buf()),
            write_index: reed::WriteIndex::None,
            compress_level: 1,
            strategy: None,
            command_line: "test".into(),
            timing: None,
            report_json: None,
            force_samples: false,
            regions: Vec::new(),
            tmp_dir: None,
            output_compressed: if bcf.is_none() { Some(true) } else { None },
            output_bcf: bcf,
            mem_bytes: None,
            no_version: false,
            no_command: false,
        };
        run_filter(Some(input), &spec, &opts)
    };
    let text_out = dir.join("text.vcf.gz");
    let bcf_out = dir.join("data.bcf");
    run(&fixture("fix.vcf.gz"), &text_out, None).expect("invariant: text output");
    run(&fixture("fix.vcf.gz"), &bcf_out, Some(6)).expect("invariant: bcf output");

    // Render the BCF back through the input path and compare the records.
    let rendered = render_via_reed(&bcf_out);
    let expected = inflate_bgzf(&std::fs::read(&text_out).expect("invariant: text output"));
    let records = |text: &str| {
        text.lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(records(&rendered), records(&expected));
    let _ = std::fs::remove_dir_all(&dir);
}

/// Reads a BCF file back with reed's own decoder, rendering every record.
fn render_via_reed(path: &Path) -> String {
    let bytes = std::fs::read(path).expect("invariant: BCF output");
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut data = Vec::new();
    let mut rest = &bytes[..];
    while let Some(total) = reed_io::bgzf::block::parse_block_size(rest) {
        inflater
            .inflate_block(&rest[..total], &mut data)
            .expect("invariant: BCF output inflates");
        rest = &rest[total..];
    }
    assert_eq!(&data[..3], b"BCF", "output is BCF");
    let l_text = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
    let dict = reed_io::bcf::BcfDict::from_header(&data[9..9 + l_text]).expect("invariant: header");
    let mut off = 9 + l_text;
    let mut out = Vec::new();
    while let Some(len) = reed_io::bcf::record_len(&data[off..]) {
        let rec = reed_io::bcf::record_ref(&data[off..off + len]).expect("invariant: record");
        reed_io::bcf::render_record(&rec, &dict, &mut out).expect("invariant: renders");
        off += len;
    }
    String::from_utf8(out).expect("invariant: VCF text")
}

/// BCF output bytes do not depend on the thread count (doc_v3 02-bcf.md §4).
#[test]
fn bcf_output_is_thread_count_independent() {
    let dir = dir("bcf-threads");
    // The fixture carries a haploid call (chr2:500).
    let spec = FilterSpec {
        policy: reed_filter::spec::PolicyRules {
            ploidy: reed_filter::spec::PloidyPolicy::AllowHaploid,
            ..reed_filter::spec::PolicyRules::default()
        },
        ..FilterSpec::default()
    };
    let run = |threads: usize, out: &Path| {
        let opts = FilterOptions {
            threads,
            output: Some(out.to_path_buf()),
            write_index: reed::WriteIndex::None,
            compress_level: 1,
            strategy: None,
            command_line: "test".into(),
            timing: None,
            report_json: None,
            force_samples: false,
            regions: Vec::new(),
            tmp_dir: None,
            output_compressed: None,
            output_bcf: Some(1),
            mem_bytes: None,
            no_version: false,
            no_command: false,
        };
        run_filter(Some(&fixture("fix.bcf")), &spec, &opts)
    };
    let mut reference: Option<Vec<u8>> = None;
    for threads in [1_usize, 2, 3, 7, 56] {
        let out = dir.join(format!("t{threads}.bcf"));
        run(threads, &out).expect("invariant: run");
        let bytes = std::fs::read(&out).expect("invariant: output");
        match &reference {
            None => reference = Some(bytes),
            Some(expected) => assert_eq!(&bytes, expected, "threads={threads}"),
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// `-q/-Q` (bcftools `--min-af/--max-af`) on the fixture: the expected
/// positions come from `bcftools view` on the same file (recorded in the
/// comment), with the INFO AC/AN values deciding (bcftools semantics).
#[test]
fn af_bounds_match_bcftools_on_the_fixture() {
    let dir = dir("bcf-af");
    let run = |min_af: Option<f64>, max_af: Option<f64>| -> Vec<String> {
        let spec = FilterSpec {
            site_after_mask: reed_filter::spec::TailRules {
                min_af,
                max_af,
                ..reed_filter::spec::TailRules::default()
            },
            // The fixture carries a haploid call; bcftools accepts it. It
            // also carries one multiallelic record (chr1:30), which the
            // filter's default biallelic gate drops; allow it so the
            // assertions below stay bcftools-comparable.
            policy: reed_filter::spec::PolicyRules {
                ploidy: reed_filter::spec::PloidyPolicy::AllowHaploid,
                multiallelic: reed_filter::spec::MultiallelicPolicy::Allow,
            },
            ..FilterSpec::default()
        };
        let out = dir.join("o.vcf.gz");
        let opts = FilterOptions {
            threads: 2,
            output: Some(out.clone()),
            write_index: reed::WriteIndex::None,
            compress_level: 1,
            strategy: None,
            command_line: "test".into(),
            timing: None,
            report_json: None,
            force_samples: false,
            regions: Vec::new(),
            tmp_dir: None,
            output_compressed: Some(true),
            output_bcf: None,
            mem_bytes: None,
            no_version: false,
            no_command: false,
        };
        run_filter(Some(&fixture("fix.vcf.gz")), &spec, &opts).expect("invariant: run");
        inflate_bgzf(&std::fs::read(&out).expect("invariant: output"))
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let mut cols = l.split('\t');
                format!(
                    "{}:{}",
                    cols.next().expect("chrom"),
                    cols.next().expect("pos")
                )
            })
            .collect()
    };
    // bcftools view -q 0.2 fix.vcf.gz -> chr1:10 chr1:20 chr1:30
    assert_eq!(run(Some(0.2), None), ["chr1:10", "chr1:20", "chr1:30"]);
    // bcftools view -Q 0.2 fix.vcf.gz -> chr2:500 (AC=1/AN=6 = 0.167)
    assert_eq!(run(None, Some(0.2)), ["chr2:500"]);
    // bcftools view -q 0.3 -Q 0.4 fix.vcf.gz -> chr1:10 (nref 0.333)
    assert_eq!(run(Some(0.3), Some(0.4)), ["chr1:10"]);
    let _ = std::fs::remove_dir_all(&dir);
}

/// The length-prefix chain must frame records correctly wherever the BGZF
/// block boundaries fall (doc_v3 02-bcf.md §2.2): re-blocking the same BCF
/// at many small payload sizes must not change a single output record.
#[test]
fn chain_framing_is_boundary_independent() {
    let dir = dir("bcf-splits");
    let bytes = inflate_bgzf_bytes(&std::fs::read(fixture("fix.bcf")).expect("invariant: fixture"));
    let spec = FilterSpec {
        policy: reed_filter::spec::PolicyRules {
            ploidy: reed_filter::spec::PloidyPolicy::AllowHaploid,
            ..reed_filter::spec::PolicyRules::default()
        },
        ..FilterSpec::default()
    };
    let run = |input: &Path, out: &Path| {
        let opts = FilterOptions {
            threads: 2,
            output: Some(out.to_path_buf()),
            write_index: reed::WriteIndex::None,
            compress_level: 1,
            strategy: None,
            command_line: "test".into(),
            timing: None,
            report_json: None,
            force_samples: false,
            regions: Vec::new(),
            tmp_dir: None,
            output_compressed: Some(true),
            output_bcf: None,
            mem_bytes: None,
            no_version: false,
            no_command: false,
        };
        run_filter(Some(input), &spec, &opts).expect("invariant: run");
        let text = inflate_bgzf(&std::fs::read(out).expect("invariant: output"));
        text.lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n")
    };
    let reference = run(&fixture("fix.bcf"), &dir.join("ref.vcf.gz"));
    assert!(!reference.is_empty(), "the reference must have records");
    for block in [13_usize, 17, 24, 32, 41, 64, 97, 128, 200, 512] {
        let path = dir.join(format!("in{block}.bcf"));
        std::fs::write(&path, reblock(&bytes, block)).expect("invariant: write");
        let got = run(&path, &dir.join(format!("out{block}.vcf.gz")));
        assert_eq!(got, reference, "block size {block}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// Filtering BCF input writes the same VCF text as filtering VCF.gz input.
#[test]
fn filter_output_is_identical_across_formats() {
    let dir = dir("filter");
    let spec = FilterSpec {
        site: reed_filter::spec::SiteRules {
            min_qual: Some(50.0),
            ..reed_filter::spec::SiteRules::default()
        },
        genotype: reed_filter::spec::GenotypeRules {
            mask_if: vec![reed_filter::spec::GenotypeMaskRule {
                field: "DP".into(),
                when: "< 8".into(),
                on_missing: reed_filter::spec::OnMissing::Keep,
            }],
            ..reed_filter::spec::GenotypeRules::default()
        },
        // The fixture carries a haploid call (chr2:500).
        policy: reed_filter::spec::PolicyRules {
            ploidy: reed_filter::spec::PloidyPolicy::AllowHaploid,
            ..reed_filter::spec::PolicyRules::default()
        },
        ..FilterSpec::default()
    };
    let run = |input: &Path, out: &Path| {
        let opts = FilterOptions {
            threads: 4,
            output: Some(out.to_path_buf()),
            write_index: reed::WriteIndex::None,
            compress_level: 1,
            strategy: None,
            command_line: "test".into(),
            timing: None,
            report_json: None,
            force_samples: false,
            regions: Vec::new(),
            tmp_dir: None,
            output_compressed: Some(true),
            output_bcf: None,
            mem_bytes: None,
            no_version: false,
            no_command: false,
        };
        run_filter(Some(input), &spec, &opts)
    };
    let bcf_out = dir.join("b.vcf.gz");
    let vcf_out = dir.join("v.vcf.gz");
    let a = run(&fixture("fix.bcf"), &bcf_out).expect("invariant: bcf filter");
    let b = run(&fixture("fix.vcf.gz"), &vcf_out).expect("invariant: vcf filter");
    assert_eq!(a.n_output_records, b.n_output_records);
    // The fixtures carry different `##bcftools_viewCommand` provenance, so
    // compare the records and the rest of the header minus that one line.
    let body = |path: &Path| {
        let text = inflate_bgzf(&std::fs::read(path).expect("invariant: output"));
        text.lines()
            .filter(|l| !l.starts_with("##bcftools_viewCommand"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(body(&bcf_out), body(&vcf_out));
    let _ = std::fs::remove_dir_all(&dir);
}
