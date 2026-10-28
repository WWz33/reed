//! BCF input gates (doc_v3 02-bcf.md B1/B3a): the same file in VCF.gz and
//! BCF must produce identical statistics, report files and filter output.

use std::path::{Path, PathBuf};

use reed::{FilterOptions, run_filter};
use reed::{Levels, Populations, PopgenRequest, WindowLayoutOptions, run_popgen};
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

/// s1,s2 -> A; s3 -> B; s4 belongs to no group.
fn pops_file() -> PathBuf {
    static PATH: std::sync::OnceLock<std::path::PathBuf> = std::sync::OnceLock::new();
    PATH.get_or_init(|| {
        let dir = std::env::temp_dir().join(format!("reed-bcf-pops-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let path = dir.join("pops.tsv");
        std::fs::write(&path, "s1\tA\ns2\tA\ns3\tB\n").expect("invariant: pops");
        path
    })
    .clone()
}

fn request() -> PopgenRequest {
    PopgenRequest {
        kinds: reed::window_kinds(2, true, true, true, false, None, None),
        populations: Populations::File(pops_file()),
        grid: reed::WindowGrid::Bp {
            size: 1000,
            step: 1000,
        },
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
        ..PopgenRequest::default()
    }
}

/// `stats` over the BCF fixture equals the same run over the VCF.gz fixture.
#[test]
fn stats_output_is_identical_across_formats() {
    let bcf = run_popgen(Some(&fixture("fix.bcf")), &request(), &budget()).expect("invariant: bcf");
    let vcf =
        run_popgen(Some(&fixture("fix.vcf.gz")), &request(), &budget()).expect("invariant: vcf.gz");
    // The request uses `InvariantMode::VariantOnly`, so the flat table is
    // rendered with the variant-only denominators (as the old sink did).
    let kinds = reed::window_kinds(2, true, true, true, false, None, None);
    assert_eq!(
        reed::windows::render_tsv(&bcf.windows, &kinds, true),
        reed::windows::render_tsv(&vcf.windows, &kinds, true)
    );
    assert_eq!(bcf.n_records, vcf.n_records);
}

/// Reports read FORMAT values: the BCF renderer feeds the same scanner.
#[test]
fn reports_are_identical_across_formats() {
    let dir = dir("reports");
    let bcf_prefix = dir.join("b");
    let vcf_prefix = dir.join("v");
    let request = |prefix: &Path| reed::reports::StatsReportRequest {
        freq: true,
        hardy: true,
        het: true,
        missing: true,
        depth: true,
        summary: false,
        hist: Vec::new(),
        singletons: false,
        populations: None,
        force_samples: false,
        targets: Vec::new(),
        targets_exclude: Vec::new(),
        prefix: prefix.to_path_buf(),
    };
    reed::reports::run_stats(Some(&fixture("fix.bcf")), &request(&bcf_prefix), &budget())
        .expect("invariant: bcf reports");
    reed::reports::run_stats(
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

/// doc_v10 D10-7: an input header without `##FILTER=<ID=PASS,...>` still
/// produces BCF — the output header declares PASS (htslib injects it on
/// both read and write, vcf.c:1485-1490/1724-1726), so a `FILTER=PASS`
/// record no longer fails the encoder's dictionary lookup.
#[test]
fn pass_is_declared_even_when_the_input_lacks_it() {
    let dir = dir("pass-injection");
    // A header without the PASS line, with a PASS record in the data.
    let header = "##fileformat=VCFv4.2\n\
                  ##contig=<ID=chr1,length=1000>\n\
                  ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
                  ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";
    let records = "chr1\t100\t.\tA\tG\t50\tPASS\tDP=5\tGT\t0/1\n\
                   chr1\t200\t.\tA\tG\t50\t.\tDP=7\tGT\t1/1\n";
    let text = format!("{header}{records}");
    let mut deflater =
        reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6 is valid");
    let mut bytes = Vec::new();
    for block in text.as_bytes().chunks(64) {
        deflater
            .deflate_block(block, &mut bytes)
            .expect("invariant: block fits");
    }
    bytes.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
    let input = dir.join("no-pass.vcf.gz");
    std::fs::write(&input, &bytes).expect("invariant: write fixture");

    let out = dir.join("out.bcf");
    let opts = FilterOptions {
        threads: 2,
        output: Some(out.clone()),
        write_index: reed::WriteIndex::None,
        compress_level: 6,
        command_line: "test".into(),
        output_bcf: Some(6),
        no_version: true,
        no_command: true,
        ..FilterOptions::default()
    };
    run_filter(Some(&input), &FilterSpec::default(), &opts)
        .expect("a PASS record must encode even when the input lacks the declaration");

    // The written BCF header declares PASS, and its dictionary resolves the
    // record's FILTER index to PASS (index 0, as htslib assumes).
    let data = {
        let raw = std::fs::read(&out).expect("invariant: output");
        let mut inflater = reed_io::bgzf::inflate::Inflater::new();
        let mut data = Vec::new();
        let mut rest = &raw[..];
        while let Some(total) = reed_io::bgzf::block::parse_block_size(rest) {
            inflater
                .inflate_block(&rest[..total], &mut data)
                .expect("invariant: inflates");
            rest = &rest[total..];
        }
        data
    };
    let l_text = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
    let header_text = String::from_utf8_lossy(&data[9..9 + l_text]).into_owned();
    assert!(
        header_text.contains("##FILTER=<ID=PASS,"),
        "output header must declare PASS: {header_text}"
    );
    let dict = reed_io::bcf::BcfDict::from_header(&data[9..9 + l_text]).expect("invariant: header");
    assert_eq!(
        dict.tag(0),
        Some(b"PASS".as_slice()),
        "PASS sits at index 0"
    );
    // Both records render (the PASS one resolves its FILTER index).
    let rendered = render_via_reed(&out);
    assert!(
        rendered.contains("chr1\t100\t.\tA\tG\t50\tPASS\t"),
        "{rendered}"
    );
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

/// A6-14 (doc_v6 03 §3): VCF 4.4 prefixed phasing and the `-0.0` INFO float
/// survive every path — text and BCF input, `-Oz` and `-Ob` output — and the
/// `-Ob` record read back through reed's own decoder equals the `-Oz` body.
#[test]
fn vcf44_phasing_and_negative_zero_survive_the_round_trip() {
    let dir = dir("vcf44");
    let input = dir.join("phasing.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.4\n\
         ##FILTER=<ID=PASS,Description=\"All\">\n\
         ##contig=<ID=Chr01,length=1000000>\n\
         ##INFO=<ID=XF,Number=1,Type=Float,Description=\"float\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"dp\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n\
         Chr01\t100\trs1\tA\tG\t50\tPASS\tXF=-0\tGT:DP\t|0/1:15\t/0|1:12\t0|1:9\t0|.:4\n\
         Chr01\t200\trs2\tA\tG\t60\tPASS\tXF=0\tGT:DP\t|1/0:15\t0|.:12\t./.:9\t|.:4\n",
    )
    .expect("invariant: fixture write");
    // The fixture's own text is the expected body (bcftools agrees on these
    // lines; the unit tests pin the same rules on codes).
    let expected: Vec<String> = std::fs::read_to_string(&input)
        .expect("invariant: fixture")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect();
    let spec = FilterSpec {
        policy: reed_filter::spec::PolicyRules {
            // `|1/0` and `|.` exercise the half-call and haploid paths.
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
            output_compressed: bcf.is_none().then_some(true),
            output_bcf: bcf,
            mem_bytes: None,
            no_version: false,
            no_command: false,
        };
        run_filter(Some(input), &spec, &opts)
    };
    let text_out = dir.join("text.vcf.gz");
    let bcf_out = dir.join("out.bcf");
    run(&input, &text_out, None).expect("invariant: text output");
    run(&input, &bcf_out, Some(6)).expect("invariant: bcf output");
    let text_body: Vec<String> = inflate_bgzf(&std::fs::read(&text_out).expect("invariant: read"))
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect();
    assert_eq!(
        text_body, expected,
        "-Oz body keeps the 4.4 prefixes and -0"
    );
    let rendered: Vec<String> = render_via_reed(&bcf_out)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect();
    assert_eq!(rendered, expected, "-Ob read-back equals -Oz body");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v7 D7-1: the sketch path (`--mind`, or `--mind --indep-pairwise`)
/// renders VCF text and has no BCF encoder. Combining it with BCF output must
/// be refused up front — never written as a BCF magic followed by text lines.
#[test]
fn sketch_path_rejects_bcf_output() {
    let dir = dir("sketch-bcf");
    let spec = FilterSpec {
        mind: Some(reed_filter::spec::MissingRate::from_f64(0.1).expect("invariant: rate")),
        // The fixture carries a haploid call (chr2:500).
        policy: reed_filter::spec::PolicyRules {
            ploidy: reed_filter::spec::PloidyPolicy::AllowHaploid,
            ..reed_filter::spec::PolicyRules::default()
        },
        ..FilterSpec::default()
    };
    let opts = FilterOptions {
        threads: 4,
        output: Some(dir.join("x.bcf")),
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
        output_bcf: Some(6),
        mem_bytes: None,
        no_version: false,
        no_command: false,
    };
    let err = run_filter(Some(&fixture("fix.vcf.gz")), &spec, &opts)
        .expect_err("sketch + BCF output must be refused");
    let msg = format!("{err}");
    assert!(
        msg.contains("BCF output requires the streaming path"),
        "unexpected error: {msg}"
    );
    assert!(
        !dir.join("x.bcf").exists(),
        "no corrupt file may be written"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 D15-3: a BCF whose header text ends exactly at a BGZF block
/// boundary (the terminating NUL starts the next block) must decode like
/// its one-byte-shorter twin. Both fixtures were hand-built from the same
/// records; only the header padding differs (`crates/reed-io/tests/
/// fixtures/bcf/boundary.bcf` vs `control.bcf`).
#[test]
fn bcf_header_at_a_block_boundary_decodes() {
    let dir = dir("boundary");
    let spec = FilterSpec::default();
    let run = |name: &str, out: &Path| {
        let opts = FilterOptions {
            threads: 1,
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
            output_bcf: None,
            mem_bytes: None,
            no_version: true,
            no_command: true,
        };
        run_filter(Some(&fixture(name)), &spec, &opts)
    };
    let boundary_out = dir.join("boundary.vcf.gz");
    run("boundary.bcf", &boundary_out).expect("boundary header must decode (doc_v15 D15-3)");
    let control_out = dir.join("control.vcf.gz");
    run("control.bcf", &control_out).expect("control decodes");
    let records = |path: &Path| -> Vec<String> {
        inflate_bgzf(&std::fs::read(path).expect("invariant: read"))
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(str::to_string)
            .collect()
    };
    let boundary = records(&boundary_out);
    let control = records(&control_out);
    assert_eq!(boundary.len(), 4, "four records expected: {boundary:?}");
    assert_eq!(boundary, control, "boundary and control records must match");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 D15-24: an explicit `IDX=` on `##contig` lines wins over the
/// header line order (htslib's rule). `idx.bcf` swaps the two contigs'
/// IDX values; `idx_control.bcf` carries the same records with consistent
/// IDX. The swapped file must resolve rid 0 to chr2 and rid 1 to chr1.
#[test]
fn bcf_contig_idx_wins_over_line_order() {
    let dir = dir("idx");
    let spec = FilterSpec::default();
    let run = |name: &str, out: &Path| {
        let opts = FilterOptions {
            threads: 1,
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
            output_bcf: None,
            mem_bytes: None,
            no_version: true,
            no_command: true,
        };
        run_filter(Some(&fixture(name)), &spec, &opts)
    };
    let swapped = dir.join("swapped.vcf.gz");
    run("idx.bcf", &swapped).expect("idx.bcf decodes");
    let control = dir.join("idx-control.vcf.gz");
    run("idx_control.bcf", &control).expect("idx_control.bcf decodes");
    let chroms = |path: &Path| -> Vec<String> {
        inflate_bgzf(&std::fs::read(path).expect("invariant: read"))
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').next().unwrap_or_default().to_string())
            .collect()
    };
    assert_eq!(
        chroms(&swapped),
        ["chr2", "chr2", "chr1", "chr1"],
        "explicit contig IDX must win (doc_v15 D15-24)"
    );
    assert_eq!(chroms(&control), ["chr1", "chr1", "chr2", "chr2"]);
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v21 M1: `pca`, `dapc` and `distance` read BCF as raw text (they did
/// not wrap their job in `InputJob`, so BCF bytes were split on '\\n'), while
/// every other command accepted the same file. The matrix reader must give
/// the same samples and dosages for both formats.
#[test]
fn m1_genotype_matrix_is_format_independent() {
    let vcf = fixture("fix.vcf.gz");
    let bcf = fixture("fix.bcf");
    let budget = budget();
    let (from_vcf, ids_vcf) =
        reed::read_genotype_matrix(reed::open_regions(Some(&vcf), &[]).expect("invariant: vcf opens"), None, &budget, reed::MatrixUse::Plain)
            .expect("VCF input");
    let (from_bcf, ids_bcf) =
        reed::read_genotype_matrix(reed::open_regions(Some(&bcf), &[]).expect("invariant: bcf opens"), None, &budget, reed::MatrixUse::Plain)
            .expect("BCF input must parse like VCF (doc_v21 M1)");
    assert_eq!(ids_vcf, ids_bcf);
    assert_eq!(from_vcf.rows, from_bcf.rows);
    assert_eq!(from_vcf.cols, from_bcf.cols);
    let bits = |m: &reed::Mat| -> Vec<u64> { m.data.iter().map(|v| v.to_bits()).collect() };
    assert_eq!(bits(&from_vcf), bits(&from_bcf), "dosages must match bit for bit");
}
