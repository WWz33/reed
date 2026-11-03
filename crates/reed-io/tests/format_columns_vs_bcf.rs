//! Text/BCF parity for the filter's typed FORMAT columns (doc_v22 R22-4).
//!
//! Property: a random record text is encoded by the native BCF encoder; the
//! requested integer columns are decoded from the framed record; the same
//! record is rendered back to text and decoded by [`GenotypeScanner`]. All
//! three value vectors must agree bit for bit. The generator varies the FORMAT
//! layout (DP present or not, so the depth fallback to `AD` is exercised), the
//! `AD` vector width (including empty and all-missing vectors), and the
//! scalar values.
//!
//! Scope: integer FORMAT columns. GT decoding has its own differential gate
//! (`tests/bcf_native_input.rs`, `tests/genotype_scan_vs_noodles.rs`).

use proptest::prelude::*;

use reed_core::genotype::{GenotypePolicies, MultiallelicPolicy, PloidyPolicy, PolicyCounters};
use reed_io::bcf::{
    BcfDict, NativeEncoder, bcf_header_text, decode_format_columns, record_len, record_ref,
    render_record, site_header,
};
use reed_io::vcf::format_layout::IntSpec;
use reed_io::vcf::genotype_scan::GenotypeScanner;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::line_edit::EditPlan;
use reed_io::vcf::site::SiteView;

const N_SAMPLES: usize = 3;

const HEADER: &str = "\
##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=1000000>\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\n";

/// Every FORMAT layout the generator draws from; both `DP`/`AD` orders appear
/// so column placement is not accidentally tied to one order.
const FORMATS: &[&str] = &[
    "GT",
    "GT:DP",
    "GT:AD",
    "GT:GQ",
    "GT:DP:AD",
    "GT:AD:DP",
    "GT:DP:GQ",
    "GT:AD:GQ",
    "GT:GQ:AD",
    "GT:DP:AD:GQ",
    "GT:AD:DP:GQ",
    "GT:GQ:DP:AD",
];

/// The columns the filter/stats ask for: the depth fallback (which sums
/// `AD` when the record has no DP key), a strict `DP`, and `GQ`.
fn specs() -> Vec<IntSpec> {
    vec![
        IntSpec::depth(),
        IntSpec::scalar(b"DP"),
        IntSpec::scalar(b"GQ"),
    ]
}

fn policies() -> GenotypePolicies {
    GenotypePolicies {
        ploidy: PloidyPolicy::HaploidAsHalfCall,
        multiallelic: MultiallelicPolicy::CollapseAlt,
        ..GenotypePolicies::default()
    }
}

/// One sample's generated values, rendered as VCF subfields.
#[derive(Clone, Debug)]
struct Sample {
    gt: String,
    dp: Option<u32>,
    ad: Option<Vec<Option<u32>>>,
    gq: Option<u32>,
}

fn arb_sample() -> impl Strategy<Value = Sample> {
    (
        prop_oneof![
            Just("0/0".to_string()),
            Just("0/1".to_string()),
            Just("1/1".to_string()),
            Just("./.".to_string()),
            Just("0|1".to_string()),
            Just("0".to_string()),
            Just("1".to_string()),
            Just(".".to_string()),
        ],
        prop::option::of(0_u32..=40_000),
        prop::option::of(prop::collection::vec(prop::option::of(0_u32..=500), 0..=4)),
        prop::option::of(0_u32..=99),
    )
        .prop_map(|(gt, dp, ad, gq)| Sample { gt, dp, ad, gq })
}

fn scalar(value: Option<u32>) -> String {
    value.map_or_else(|| ".".to_string(), |v| v.to_string())
}

fn ad_text(values: Option<&[Option<u32>]>) -> String {
    match values {
        None => ".".to_string(),
        Some(values) => values
            .iter()
            .map(|value| scalar(*value))
            .collect::<Vec<_>>()
            .join(","),
    }
}

/// One sample column under `format`.
fn sample_text(format: &str, sample: &Sample) -> String {
    format
        .split(':')
        .map(|key| match key {
            "GT" => sample.gt.clone(),
            "DP" => scalar(sample.dp),
            "AD" => ad_text(sample.ad.as_deref()),
            "GQ" => scalar(sample.gq),
            other => unreachable!("undefined FORMAT key {other}"),
        })
        .collect::<Vec<_>>()
        .join(":")
}

/// Encodes `record`, decodes the native columns, and renders the record back
/// to text.
fn encode_and_decode(record: &[u8], specs: &[IntSpec]) -> (Vec<Vec<u32>>, Vec<u8>) {
    let header: noodles::vcf::Header = HEADER.parse().expect("invariant: header parses");
    let mut encoder = NativeEncoder::new(&header).expect("invariant: native encoder");
    let line = LineRef::parse(record).expect("invariant: record parses");
    let site = SiteView::parse(&line).expect("invariant: site parses");
    let mut framed = Vec::new();
    encoder
        .encode_text(&line, &site, N_SAMPLES, &EditPlan::default(), &mut framed)
        .expect("invariant: native encode");

    let len = record_len(&framed).expect("invariant: framed record");
    let rec = record_ref(&framed[..len]).expect("invariant: record ref");
    let rec_header = site_header(&rec).expect("invariant: site header");
    let header_text = bcf_header_text(&header).expect("invariant: BCF header text");
    let dict = BcfDict::from_header(&header_text).expect("invariant: dictionaries");

    let mut int_columns = vec![Vec::new(); specs.len()];
    let mut float_columns: Vec<Vec<f32>> = Vec::new();
    decode_format_columns(
        &rec,
        &rec_header,
        &dict,
        specs,
        &[],
        &mut int_columns,
        &mut float_columns,
    )
    .expect("invariant: native columns decode");

    let mut rendered = Vec::new();
    render_record(&rec, &dict, &mut rendered).expect("invariant: render");
    (int_columns, rendered)
}

/// Decodes the same columns through the text scanner.
fn scan_columns(record: &[u8], specs: &[IntSpec]) -> Vec<Vec<u32>> {
    let line = LineRef::parse(record).expect("invariant: record parses");
    let site = SiteView::parse(&line).expect("invariant: site parses");
    let mut scanner = GenotypeScanner::new_with_specs(N_SAMPLES, specs, &[], policies());
    let mut counters = PolicyCounters::default();
    scanner
        .scan(&line, &site, &mut counters)
        .expect("invariant: text scan")
        .int_columns
        .to_vec()
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 500,
        ..ProptestConfig::default()
    })]

    /// Native BCF column decode == text scan of the original line == text
    /// scan of the rendered line, for every generated FORMAT shape.
    #[test]
    fn native_columns_match_both_text_paths(
        format in prop::sample::select(FORMATS.to_vec()),
        samples in prop::collection::vec(arb_sample(), N_SAMPLES..=N_SAMPLES),
        pos in 1_u32..=1_000_000,
    ) {
        let record = format!(
            "chr1\t{pos}\t.\tA\tG\t.\tPASS\t.\t{format}\t{}\n",
            samples
                .iter()
                .map(|sample| sample_text(format, sample))
                .collect::<Vec<_>>()
                .join("\t")
        );
        let specs = specs();
        let (native, rendered) = encode_and_decode(record.as_bytes(), &specs);
        let from_original = scan_columns(record.as_bytes(), &specs);
        let from_rendered = scan_columns(&rendered, &specs);

        prop_assert_eq!(&native, &from_original, "record: {}", record.trim_end());
        prop_assert_eq!(&native, &from_rendered, "record: {}", record.trim_end());
        prop_assert_eq!(
            native[0].len(),
            N_SAMPLES,
            "every column has one value per sample"
        );
    }
}

/// V22-R2-2: the text and native paths agree at the `i32` boundaries that
/// used to saturate or collide with the missing sentinel.
#[test]
fn boundary_ints_match_text_and_native() {
    let record = "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP:AD:GQ\
                  \t0/1:2147483647:2147483647,2147483647,2147483647:2147483647\
                  \t0/1:2147483648:2147483647,2147483648,-5:5000000000\
                  \t0/1:+50:+5,+0:0007\n";
    let specs = specs();
    let (native, rendered) = encode_and_decode(record.as_bytes(), &specs);
    let from_original = scan_columns(record.as_bytes(), &specs);
    let from_rendered = scan_columns(&rendered, &specs);
    assert_eq!(native, from_original, "text scan of the source");
    assert_eq!(native, from_rendered, "text scan of the rendered BCF");
    // depth (DP present, so no AD fallback) / strict DP / GQ.
    assert_eq!(
        native[0],
        vec![2_147_483_647, reed_io::vcf::genotype_scan::INT_MISSING, 50]
    );
    assert_eq!(native[1], native[0]);
    assert_eq!(
        native[2],
        vec![2_147_483_647, reed_io::vcf::genotype_scan::INT_MISSING, 7]
    );
}
