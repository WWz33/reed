//! Text → BCF encoding gates (doc_v3 02-bcf.md B2/B4): the native encoder
//! and noodles' encoder must produce the same record bytes, and the result
//! must render back to the record text the fixture started from.

use std::io::Read;

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

/// The fixture's header (through `#CHROM`) and its record lines.
fn fixture_parts() -> (String, Vec<Vec<u8>>) {
    let data = inflate_all(&fixture("fix.vcf.gz"));
    let text = String::from_utf8(data).expect("invariant: UTF-8 VCF");
    let chrom = text.find("#CHROM").expect("invariant: header");
    let end = text[chrom..]
        .find('\n')
        .map_or(text.len(), |i| chrom + i + 1);
    let header = text[..end].to_string();
    let records = text[end..]
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let mut v = l.as_bytes().to_vec();
            v.push(b'\n');
            v
        })
        .collect();
    (header, records)
}

/// Both encoders produce identical bytes for every fixture record.
#[test]
fn native_and_noodles_encoders_agree() {
    let (header_text, records) = fixture_parts();
    let header: noodles::vcf::Header = header_text.parse().expect("invariant: header parses");
    let mut native = reed_io::bcf::NativeEncoder::new(&header).expect("invariant: native encoder");
    let mut noodles = reed_io::bcf::BcfEncoder::new().expect("invariant: noodles encoder");
    noodles.write_header(&header).expect("invariant: header");
    for (i, record) in records.iter().enumerate() {
        let mut expected = Vec::new();
        encode_text(&mut native, &header, record, &mut expected);
        let mut got = Vec::new();
        noodles
            .encode_record(&header, record, &mut got)
            .expect("invariant: noodles encode");
        assert_eq!(
            got,
            expected,
            "record {i}: {}",
            String::from_utf8_lossy(record)
        );
    }
}

/// The native encoder's text entry with a default (no-op) edit.
fn encode_text(
    encoder: &mut reed_io::bcf::NativeEncoder,
    header: &noodles::vcf::Header,
    record: &[u8],
    out: &mut Vec<u8>,
) {
    let lr = reed_io::vcf::line::LineRef::parse(record).expect("invariant: full line");
    let site = reed_io::vcf::site::SiteView::parse(&lr).expect("invariant: site parses");
    encoder
        .encode_text(
            &lr,
            &site,
            header.sample_names().len(),
            &reed_io::vcf::line_edit::EditPlan::default(),
            out,
        )
        .expect("invariant: native encode");
}

/// Encoding the fixture's records and decoding them again reproduces the
/// record text (the fixture's expected rendering).
#[test]
fn encode_round_trips_through_the_renderer() {
    let (header_text, records) = fixture_parts();
    let header: noodles::vcf::Header = header_text.parse().expect("invariant: header parses");
    let mut encoder = reed_io::bcf::NativeEncoder::new(&header).expect("invariant: encoder");
    // A synthetic BCF stream: header text + records, rendered back.
    let mut stream = Vec::new();
    for record in &records {
        encode_text(&mut encoder, &header, record, &mut stream);
    }
    // The encoder numbers tags with noodles' StringMaps; the file header a
    // real run writes carries explicit IDX for every tag, which pins htslib
    // to that same numbering. Decode with that header text.
    let out_header = reed_io::bcf::bcf_header_text(&header).expect("invariant: BCF header text");
    let dict = reed_io::bcf::BcfDict::from_header(&out_header).expect("invariant: dictionaries");
    let mut off = 0;
    let mut rendered = Vec::new();
    while let Some(len) = reed_io::bcf::record_len(&stream[off..]) {
        let rec = reed_io::bcf::record_ref(&stream[off..off + len]).expect("invariant: record");
        reed_io::bcf::render_record(&rec, &dict, &mut rendered).expect("invariant: render");
        off += len;
    }
    let expected: Vec<u8> = records.concat();
    assert_eq!(
        String::from_utf8_lossy(&rendered),
        String::from_utf8_lossy(&expected)
    );
}
