//! Native BCF filter inputs (doc_v4 O-4): the site line the head rules and
//! expressions evaluate, and the typed FORMAT columns the mask consumes, must
//! agree with the text path bit for bit. The fixture is the same one
//! `bcf_render.rs` locks to htslib (haploid calls, multiallelics, flags,
//! strings, missing values).

use std::io::Read;

use reed_core::genotype::{GenotypePolicies, MultiallelicPolicy, PloidyPolicy, PolicyCounters};
use reed_io::bcf::{
    BcfDict, GtQuery, GtTable, decode_format_columns, decode_record_gt, format_key, record_len,
    record_ref, render_record, render_site_line, site_header,
};
use reed_io::vcf::genotype_scan::GenotypeScanner;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;

fn fixture(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/bcf")
        .join(name)
}

fn inflate_all(path: &std::path::Path) -> Vec<u8> {
    let mut reader = std::io::BufReader::with_capacity(
        1 << 16,
        std::fs::File::open(path).expect("invariant: fixture opens"),
    );
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut data = Vec::new();
    let mut header = [0_u8; 18];
    while reader.read_exact(&mut header).is_ok() {
        let total = reed_io::bgzf::block::parse_block_size(&header).expect("invariant: BGZF");
        let mut block = header.to_vec();
        block.resize(total, 0);
        reader
            .read_exact(&mut block[18..])
            .expect("invariant: full block");
        inflater
            .inflate_block(&block, &mut data)
            .expect("invariant: inflates");
    }
    data
}

/// The policies both sides run under (the fixture has haploid calls, so the
/// strict default would reject them).
fn policies() -> GenotypePolicies {
    GenotypePolicies {
        ploidy: PloidyPolicy::HaploidAsHalfCall,
        multiallelic: MultiallelicPolicy::CollapseAlt,
        ..GenotypePolicies::default()
    }
}

#[test]
fn native_inputs_match_the_text_path() {
    let data = inflate_all(&fixture("fix.bcf"));
    assert_eq!(&data[..3], b"BCF");
    let l_text = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
    let dict = BcfDict::from_header(&data[9..9 + l_text]).expect("invariant: header");
    let table = GtTable::new(policies());
    let gt_key = format_key(&dict, b"GT");

    let mut off = 9 + l_text;
    let mut checked = 0_usize;
    while let Some(len) = record_len(&data[off..]) {
        let bytes = &data[off..off + len];
        let rec = record_ref(bytes).expect("invariant: record");
        let header = site_header(&rec).expect("invariant: header fields");

        // ---- text side: render the whole record and scan it ----
        let mut full = Vec::new();
        render_record(&rec, &dict, &mut full).expect("invariant: renders");
        let full = String::from_utf8(full).expect("invariant: ASCII");
        let line = full.trim_end_matches('\n');
        let lr = LineRef::parse(line.as_bytes()).expect("invariant: parses");
        let site = SiteView::parse(&lr).expect("invariant: site parses");
        let mut scanner = GenotypeScanner::new(header.n_sample, &[b"DP"], &[], policies());
        let mut counters = PolicyCounters::default();
        let text = scanner
            .scan(&lr, &site, &mut counters)
            .expect("invariant: scans");

        // ---- native side: the site line plus the typed columns ----
        let mut site_line = Vec::new();
        render_site_line(&rec, &dict, &mut site_line).expect("invariant: site line");
        let nlr = LineRef::parse(&site_line).expect("invariant: native line parses");
        let nsite = SiteView::parse(&nlr).expect("invariant: native site parses");
        assert_eq!(nsite.contig_name, site.contig_name, "record {checked}");
        assert_eq!(nsite.pos1, site.pos1, "record {checked}");
        assert_eq!(nsite.reference, site.reference, "record {checked}");
        assert_eq!(nsite.alts, site.alts, "record {checked}");
        assert_eq!(nsite.qual, site.qual, "record {checked}");
        assert_eq!(nsite.filter, site.filter, "record {checked}");
        assert_eq!(
            nsite.info, site.info,
            "record {checked}: INFO must be identical"
        );
        assert_eq!(nlr.id(), lr.id(), "record {checked}");
        assert_eq!(
            nlr.n_sample_cols(),
            lr.n_sample_cols(),
            "record {checked}: N_SAMPLES must stay visible"
        );

        let mut int_columns = vec![Vec::new()];
        let mut float_columns: Vec<Vec<f32>> = Vec::new();
        decode_format_columns(
            &rec,
            &header,
            &dict,
            &[b"DP"],
            &[],
            &mut int_columns,
            &mut float_columns,
        )
        .expect("invariant: columns decode");
        assert_eq!(
            int_columns[0], text.int_columns[0],
            "record {checked}: the DP column must equal the scanned text"
        );

        let mut codes = Vec::new();
        let mut native_counters = PolicyCounters::default();
        decode_record_gt(
            &rec,
            &header,
            &GtQuery {
                gt_key,
                n_samples: header.n_sample,
                contig_name: nsite.contig_name,
                pos1: nsite.pos1,
            },
            &table,
            &mut codes,
            &mut native_counters,
        )
        .expect("invariant: GT decodes");
        assert_eq!(codes, text.codes, "record {checked}: GT codes must match");
        checked += 1;
        off += len;
    }
    assert!(
        checked >= 4,
        "the fixture must exercise several records, got {checked}"
    );
}
