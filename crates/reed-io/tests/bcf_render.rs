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

/// Appends a typed integer (descriptor + payload), the BCF key encoding.
fn push_typed_int(out: &mut Vec<u8>, value: i64) {
    out.push(0x11); // count 1, INT8
    out.push(value as i8 as u8);
}

/// One framed BCF record with a `(NULL, count=1)` descriptor — the INFO
/// field when `null_info`, else the FORMAT (GT) field (doc_v10 D10-2).
///
/// Layout (hts-specs §6.3.1): chrom=0, pos=0, rlen=1, QUAL missing, one
/// allele (REF=A), no ALT, ID empty, FILTER empty, and exactly one INFO
/// (DP) or one FORMAT (GT) field whose value descriptor declares NULL type
/// with count 1 — which htslib rejects (`bcf_record_check`,
/// vcf.c:2231-2236) and reed used to feed into the empty-payload renderers.
fn synthetic_null_count_record(null_info: bool) -> Vec<u8> {
    let mut shared = Vec::new();
    shared.extend_from_slice(&0_i32.to_le_bytes()); // chrom
    shared.extend_from_slice(&0_i32.to_le_bytes()); // pos
    shared.extend_from_slice(&1_i32.to_le_bytes()); // rlen
    shared.extend_from_slice(&0x7F80_0001_u32.to_le_bytes()); // QUAL missing
    let n_info: u32 = u32::from(null_info);
    shared.extend_from_slice(&(1_u32 << 16 | n_info).to_le_bytes()); // n_allele=1
    let n_fmt: u32 = u32::from(!null_info);
    let n_sample: u32 = u32::from(!null_info);
    shared.extend_from_slice(&(n_fmt << 24 | n_sample).to_le_bytes());
    shared.push(0x07); // ID: CHAR, count 0 (empty)
    shared.push(0x17); // REF: CHAR, count 1
    shared.push(b'A');
    shared.push(0x01); // FILTER: INT8, count 0
    if null_info {
        push_typed_int(&mut shared, 1); // INFO key DP
        shared.push(0x10); // (NULL, count 1): the malformed descriptor
    }
    let mut indiv = Vec::new();
    if !null_info {
        push_typed_int(&mut indiv, 2); // FORMAT key GT
        indiv.push(0x10); // (NULL, count 1)
    }
    let mut out = Vec::new();
    out.extend_from_slice(&(shared.len() as u32).to_le_bytes());
    out.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
    out.extend_from_slice(&shared);
    out.extend_from_slice(&indiv);
    out
}

/// doc_v10 D10-2: a `(NULL, count > 0)` descriptor is malformed; the
/// renderer must report it, not panic slicing an empty payload. `(NULL, 0)`
/// (an INFO flag) stays legal.
#[test]
fn null_descriptor_with_nonzero_count_errors_without_panicking() {
    let dict = reed_io::bcf::BcfDict::from_header(
        b"##fileformat=VCFv4.2\n\
          ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
          ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
          ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
          ##contig=<ID=chr1,length=1000>\n",
    )
    .expect("invariant: header");
    for null_info in [true, false] {
        let bytes = synthetic_null_count_record(null_info);
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let rec = reed_io::bcf::record_ref(&bytes).expect("invariant: framed");
            let mut out = Vec::new();
            reed_io::bcf::render_record(&rec, &dict, &mut out)
        }));
        match outcome {
            Ok(Ok(())) => panic!("NULL count=1 must be rejected (null_info={null_info})"),
            Ok(Err(_)) => {}
            Err(_) => panic!("renderer panicked on a NULL descriptor (null_info={null_info})"),
        }
    }
}
