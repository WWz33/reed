//! Differential test: the hand-written sample scanner must agree with a
//! noodles-based reference decode, byte for byte, on the edge-case fixture
//! (`docs/architecture-spec.md` §6.4, decision D1's correctness gate).
//!
//! Semantic scope: the reference (legacy `record::decode`) treats haploid
//! calls as half calls and allele indices above 1 as ALT — the scanner runs
//! under `HaploidAsHalfCall` + `CollapseAlt` to match. Polyploid genotypes
//! are out of scope (the legacy path silently takes the first two alleles;
//! the scanner rejects them — the A6 fix goes the other way).

use std::io::BufReader;

/// One record's expected planes and int columns (as noodles reports them).
type SampleExpectations = (Vec<i8>, Vec<Option<u32>>, Vec<Option<u32>>);
use std::path::Path;

use reed_io::Locus;
use reed_io::vcf::fixed::parse_pos1;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::sample_scan::{INT_MISSING, SampleScanner};
use reed_io::record;
use reed_core::genotype::{
    MultiallelicPolicy, PloidyPolicy, PlaneBlock, PolicyCounters,
};

const FIXTURE: &str = "tests/fixtures/scan_edge.vcf";

fn fixture_path() -> std::path::PathBuf {
    // CARGO_MANIFEST_DIR keeps this independent of the test runner's cwd.
    Path::new(env!("CARGO_MANIFEST_DIR")).join(FIXTURE)
}

/// Reference decode via noodles RecordBufs (the legacy transport path).
fn reference(
    path: &Path,
) -> Vec<SampleExpectations> {
    let file = std::fs::File::open(path).expect("invariant: fixture exists");
    let mut reader =
        noodles::vcf::io::Reader::new(BufReader::new(file));
    let header = reader.read_header().expect("invariant: header parses");
    let n_samples = header.sample_names().len();
    let mut out = Vec::new();
    for buf in reader.record_bufs(&header) {
        let buf = buf.expect("invariant: record parses");
        let view = record::decode(&buf, n_samples).expect("invariant: reference decode");
        out.push((view.gt, view.dp.clone(), view.gq.clone()));
    }
    out
}

/// Scanner decode of the same records.
fn scanned(path: &Path) -> Vec<(PlaneBlock, Vec<Vec<u32>>)> {
    let bytes = std::fs::read(path).expect("invariant: fixture exists");
    let mut scanner = SampleScanner::new(&[b"DP", b"GQ"]);
    let policies = reed_core::genotype::GenotypePolicies {
        ploidy: PloidyPolicy::HaploidAsHalfCall,
        multiallelic: MultiallelicPolicy::CollapseAlt,
        ..reed_core::genotype::GenotypePolicies::default()
    };
    let mut out = Vec::new();
    for line_bytes in bytes.split(|&b| b == b'\n') {
        if line_bytes.is_empty() || line_bytes.starts_with(b"#") {
            continue;
        }
        let line = LineRef::parse(line_bytes).expect("invariant: fixture line has 9 tabs");
        let locus = Locus { contig_name: "chr1".into(), pos1: 1 };
        let pos1 = parse_pos1(line.pos_text(), &locus).expect("invariant: POS");
        let mut planes = PlaneBlock::new(1, 3, false);
        let mut int_out = vec![Vec::new(); scanner.n_int_columns()];
        let mut counters = PolicyCounters::default();
        scanner
            .scan_record(&line, 0, pos1, &mut planes, &mut int_out, &policies, &mut counters)
            .expect("invariant: fixture line scans under the relaxed policies");
        out.push((planes, int_out));
    }
    out
}

#[test]
fn scanner_matches_noodles_on_every_fixture_record() {
    let expected = reference(&fixture_path());
    let actual = scanned(&fixture_path());
    assert_eq!(expected.len(), actual.len(), "both sides must see the same records");

    for (record_idx, ((gt, dp, gq), (planes, int_cols))) in expected.iter().zip(&actual).enumerate() {
        let n_samples = gt.len() / 2;
        assert_eq!(int_cols.len(), 2);
        assert_eq!(int_cols[0].len(), n_samples, "DP column per sample");
        assert_eq!(int_cols[1].len(), n_samples, "GQ column per sample");

        for sample in 0..n_samples {
            let (a, b) = (gt[2 * sample], gt[2 * sample + 1]);
            // Slot 0: -1 → missing; >= 1 → ALT; 0 → called REF.
            assert_eq!(planes.plane(0, 1)[sample / 64] >> (sample % 64) & 1 == 1, a == -1,
                "record {record_idx} sample {sample} slot0 missing");
            assert_eq!(planes.plane(0, 0)[sample / 64] >> (sample % 64) & 1 == 1, a >= 1,
                "record {record_idx} sample {sample} slot0 alt");
            assert_eq!(planes.plane(0, 3)[sample / 64] >> (sample % 64) & 1 == 1, b == -1,
                "record {record_idx} sample {sample} slot1 missing");
            assert_eq!(planes.plane(0, 2)[sample / 64] >> (sample % 64) & 1 == 1, b >= 1,
                "record {record_idx} sample {sample} slot1 alt");

            let expect_int = |v: &Option<u32>| v.unwrap_or(INT_MISSING);
            assert_eq!(int_cols[0][sample], expect_int(&dp[sample]),
                "record {record_idx} sample {sample} DP");
            assert_eq!(int_cols[1][sample], expect_int(&gq[sample]),
                "record {record_idx} sample {sample} GQ");
        }
    }
}
