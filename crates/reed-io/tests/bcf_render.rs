//! BCF → VCF rendering parity (doc_v3 02-bcf.md B3a): the fixture was
//! produced by `bcftools view -Ob` and `fix.expected.vcf` is
//! `bcftools view -H` of that BCF, so this locks the renderer to htslib's
//! text output (missing values, phased GT, haploid calls, multiallelics,
//! flags, strings and `%g` floats).

use std::io::Read;

fn fixture(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/bcf")
        .join(name)
}

/// Reads a whole BGZF file into one buffer.
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

/// Every record rendered, joined (the BCF fixture is a single chunk).
fn render_bcf(path: &std::path::Path) -> String {
    let data = inflate_all(path);
    assert_eq!(&data[..3], b"BCF", "fixture is BCF");
    let l_text = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
    let dict = reed_io::bcf::BcfDict::from_header(&data[9..9 + l_text]).expect("invariant: header");
    let mut off = 9 + l_text;
    let mut out = Vec::new();
    while let Some(len) = reed_io::bcf::record_len(&data[off..]) {
        let rec = reed_io::bcf::record_ref(&data[off..off + len]).expect("invariant: record");
        reed_io::bcf::render_record(&rec, &dict, &mut out).expect("invariant: renders");
        off += len;
    }
    String::from_utf8(out).expect("invariant: ASCII")
}

#[test]
fn rendering_matches_bcftools() {
    let got = render_bcf(&fixture("fix.bcf"));
    let want = std::fs::read_to_string(fixture("fix.expected.vcf")).expect("invariant: expected");
    assert_eq!(got, want);
}

#[test]
fn dictionaries_come_from_the_header() {
    let data = inflate_all(&fixture("fix.bcf"));
    let l_text = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
    let dict = reed_io::bcf::BcfDict::from_header(&data[9..9 + l_text]).expect("invariant: header");
    // PASS is 0, then LowQual, then INFO/FORMAT tags in header order.
    assert_eq!(dict.tag(0), Some(b"PASS".as_slice()));
    assert_eq!(dict.tag(1), Some(b"LowQual".as_slice()));
    assert_eq!(dict.contig(0), Some(b"chr1".as_slice()));
    assert_eq!(dict.contig(1), Some(b"chr2".as_slice()));
}
