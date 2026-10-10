//! Differential test: the unified [`GenotypeScanner`] must agree with a
//! noodles-based reference decode on the edge-case fixture
//! (`docs/rounds/v2/03-architecture.md` §3.4: the scanner is the one decoder, and
//! this test pins it).
//!
//! Semantic scope: the reference treats haploid calls as half calls and
//! allele indices above 1 as ALT, so the scanner runs under
//! `HaploidAsHalfCall` + `CollapseAlt` (the stats/report policies). Polyploid
//! genotypes are out of scope (the scanner rejects them; the reference has no
//! representation for them either).

use std::io::BufReader;
use std::path::{Path, PathBuf};

use reed_core::genotype::{
    GenotypePolicies, GtCode, MultiallelicPolicy, Ploidy, PloidyPolicy, PolicyCounters,
};
use reed_io::vcf::genotype_scan::{GenotypeScanner, INT_MISSING};
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;

const FIXTURE: &str = "tests/fixtures/scan_edge.vcf";

fn fixture_path() -> PathBuf {
    // CARGO_MANIFEST_DIR keeps this independent of the test runner's cwd.
    Path::new(env!("CARGO_MANIFEST_DIR")).join(FIXTURE)
}

fn policies() -> GenotypePolicies {
    GenotypePolicies {
        ploidy: PloidyPolicy::HaploidAsHalfCall,
        multiallelic: MultiallelicPolicy::CollapseAlt,
        ..GenotypePolicies::default()
    }
}

/// One record's raw slots per sample (2 slots, index order) with the
/// genotype's arity, plus the DP/GQ columns noodles reports. The arity is
/// what tells a haploid call from a half call.
type RawRecord = (
    Vec<([Option<u16>; 2], u8)>,
    Vec<Option<u32>>,
    Vec<Option<u32>>,
);

/// Reference decode via noodles RecordBufs.
fn reference(path: &Path) -> Vec<RawRecord> {
    use noodles::vcf::variant::record::samples::Series as _;
    use noodles::vcf::variant::record::samples::series::Value;

    let file = std::fs::File::open(path).expect("invariant: fixture exists");
    let mut reader = noodles::vcf::io::Reader::new(BufReader::new(file));
    let header = reader.read_header().expect("invariant: header parses");
    let n_samples = header.sample_names().len();
    let mut out = Vec::new();
    for buf in reader.record_bufs(&header) {
        let buf = buf.expect("invariant: record parses");
        let mut slots = vec![([None::<u16>; 2], 2_u8); n_samples];
        let mut dp: Vec<Option<u32>> = vec![None; n_samples];
        let mut gq: Vec<Option<u32>> = vec![None; n_samples];
        for series in buf.samples().series() {
            let name = series.name(&header).expect("invariant: series name");
            match name {
                "GT" => {
                    for (i, value) in series.iter(&header).enumerate() {
                        let Ok(Some(Value::Genotype(genotype))) = value else {
                            // noodles reports a lone `.` as no genotype value
                            // at all; the VCF text is the haploid missing call.
                            slots[i] = ([None, None], 1);
                            continue;
                        };
                        slots[i].1 =
                            u8::try_from(genotype.iter().count().min(2)).expect("invariant: ≤ 2");
                        for (slot, pair) in genotype.iter().enumerate().take(2) {
                            let (allele, _) = pair.expect("invariant: allele pair");
                            if let Some(index) = allele {
                                slots[i].0[slot] = Some(
                                    u16::try_from(index).expect("invariant: allele index fits u16"),
                                );
                            }
                        }
                    }
                }
                "DP" | "GQ" => {
                    let column = if name == "DP" { &mut dp } else { &mut gq };
                    for (i, value) in series.iter(&header).enumerate() {
                        let Ok(Some(Value::Integer(v))) = value else {
                            continue;
                        };
                        if i < column.len() && v >= 0 {
                            column[i] = Some(v as u32);
                        }
                    }
                }
                _ => {}
            }
        }
        out.push((slots, dp, gq));
    }
    out
}

/// Scanner decode of the same records: codes, DP/GQ columns and raw alleles.
type ScannedRecord = (Vec<GtCode>, Vec<Vec<u32>>, Vec<Vec<f32>>, Vec<Option<u16>>);

fn scanned(path: &Path) -> Vec<ScannedRecord> {
    let bytes = std::fs::read(path).expect("invariant: fixture exists");
    let n_samples = 3;
    let mut scanner =
        GenotypeScanner::new(n_samples, &[b"DP", b"GQ"], &[], policies()).tracking_alleles();
    let mut out = Vec::new();
    let mut counters = PolicyCounters::default();
    for line_bytes in bytes.split(|&b| b == b'\n') {
        if line_bytes.is_empty() || line_bytes.starts_with(b"#") {
            continue;
        }
        let line = LineRef::parse(line_bytes).expect("invariant: fixture line has 9 tabs");
        let site = SiteView::parse(&line).expect("invariant: site parses");
        let scan = scanner
            .scan(&line, &site, &mut counters)
            .expect("invariant: fixture line scans under the relaxed policies");
        out.push((
            scan.codes.to_vec(),
            scan.int_columns.to_vec(),
            scan.float_columns.to_vec(),
            scan.alleles.expect("invariant: tracking is on").to_vec(),
        ));
    }
    out
}

#[test]
fn scanner_matches_noodles_on_every_fixture_record() {
    let expected = reference(&fixture_path());
    let actual = scanned(&fixture_path());
    assert_eq!(
        expected.len(),
        actual.len(),
        "both sides must see the same records"
    );

    for (record_idx, ((slots, dp, gq), (codes, int_cols, _floats, alleles))) in
        expected.iter().zip(&actual).enumerate()
    {
        let n_samples = slots.len();
        assert_eq!(
            codes.len(),
            n_samples,
            "record {record_idx}: one code per sample"
        );
        assert_eq!(int_cols.len(), 2);
        assert_eq!(
            int_cols[0].len(),
            n_samples,
            "record {record_idx}: DP per sample"
        );
        assert_eq!(
            int_cols[1].len(),
            n_samples,
            "record {record_idx}: GQ per sample"
        );
        assert_eq!(
            alleles.len(),
            2 * n_samples,
            "record {record_idx}: raw alleles per sample"
        );

        for sample in 0..n_samples {
            let (raw, arity) = slots[sample];
            let ploidy = if arity == 1 {
                Ploidy::Haploid
            } else {
                Ploidy::Diploid
            };
            let mut scratch = PolicyCounters::default();
            let expected_code = GtCode::from_alleles(raw, ploidy, false, &policies(), &mut scratch)
                .unwrap_or_else(|_| {
                    panic!("record {record_idx} sample {sample}: reference decode")
                });
            assert_eq!(
                codes[sample].state(),
                expected_code.state(),
                "record {record_idx} sample {sample}: slot states"
            );
            assert_eq!(
                codes[sample].is_haploid(),
                expected_code.is_haploid(),
                "record {record_idx} sample {sample}: haploid flag"
            );
            assert_eq!(
                &alleles[2 * sample..2 * sample + 2],
                &raw[..],
                "record {record_idx} sample {sample}: raw allele indices"
            );
            let expect_int = |v: &Option<u32>| v.unwrap_or(INT_MISSING);
            assert_eq!(
                int_cols[0][sample],
                expect_int(&dp[sample]),
                "record {record_idx} sample {sample} DP"
            );
            assert_eq!(
                int_cols[1][sample],
                expect_int(&gq[sample]),
                "record {record_idx} sample {sample} GQ"
            );
        }
    }
}

#[test]
fn scanner_policies_are_counted_once_per_trigger() {
    let actual = scanned(&fixture_path());
    // Record 200 holds the haploid run `1`, `0`, `.`: two called haploids
    // (half calls) and one haploid missing (no trigger). Record 300 holds
    // `10/11` (two collapsed alleles); record 500 holds `1/2` (two).
    let mut expected_scratch = PolicyCounters::default();
    for (record_idx, (slots, _, _)) in reference(&fixture_path()).iter().enumerate() {
        for (raw, arity) in slots {
            let ploidy = if *arity == 1 {
                Ploidy::Haploid
            } else {
                Ploidy::Diploid
            };
            let _ = GtCode::from_alleles(*raw, ploidy, false, &policies(), &mut expected_scratch)
                .unwrap_or_else(|_| panic!("record {record_idx}: reference decode"));
        }
    }
    assert!(
        expected_scratch.n_haploid_as_half_call > 0,
        "fixture must exercise haploid calls"
    );
    assert!(
        expected_scratch.n_multiallelic_collapsed > 0,
        "fixture must exercise collapse"
    );

    // The scanner's counters must match the same triggers exactly.
    let mut scanner = GenotypeScanner::new(3, &[b"DP"], &[], policies());
    let mut counters = PolicyCounters::default();
    let bytes = std::fs::read(fixture_path()).expect("invariant: fixture exists");
    for line_bytes in bytes.split(|&b| b == b'\n') {
        if line_bytes.is_empty() || line_bytes.starts_with(b"#") {
            continue;
        }
        let line = LineRef::parse(line_bytes).expect("invariant: 9 tabs");
        let site = SiteView::parse(&line).expect("invariant: site parses");
        scanner
            .scan(&line, &site, &mut counters)
            .expect("invariant: relaxed policies");
    }
    assert_eq!(counters, expected_scratch);
    let _ = actual;
}
