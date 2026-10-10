//! Differential test (doc_v3 02-bcf.md §4): our BCF decoder + renderer must
//! produce exactly what noodles' independent BCF reader renders — for the
//! whole fixture and, when the local data file exists, for the first
//! 100,000 records of `clean.bcf`.

use std::io::Read;

use noodles::bcf;
use noodles::vcf;
use noodles::vcf::variant::io::Write as _;

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

/// Renders `limit` records of `path` with our decoder.
fn render_ours(path: &std::path::Path, limit: usize) -> String {
    let data = inflate_all(path);
    assert_eq!(&data[..3], b"BCF");
    let l_text = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
    let dict = reed_io::bcf::BcfDict::from_header(&data[9..9 + l_text]).expect("invariant: header");
    let mut off = 9 + l_text;
    let mut out = Vec::new();
    let mut n = 0;
    while n < limit {
        let Some(len) = reed_io::bcf::record_len(&data[off..]) else {
            break;
        };
        let rec = reed_io::bcf::record_ref(&data[off..off + len]).expect("invariant: record");
        reed_io::bcf::render_record(&rec, &dict, &mut out).expect("invariant: renders");
        off += len;
        n += 1;
    }
    String::from_utf8(out).expect("invariant: VCF text")
}

/// Renders `limit` records of `path` with noodles' reader.
fn render_noodles(path: &std::path::Path, limit: usize) -> String {
    let file = std::fs::File::open(path).expect("invariant: fixture opens");
    let mut reader = bcf::io::Reader::new(std::io::BufReader::new(file));
    let header = reader.read_header().expect("invariant: header");
    let mut out = Vec::new();
    {
        let mut writer = vcf::io::Writer::new(&mut out);
        writer.write_header(&header).expect("invariant: header");
        let mut buffer = vcf::variant::RecordBuf::default();
        for _ in 0..limit {
            match reader.read_record_buf(&header, &mut buffer) {
                Ok(0) => break,
                Ok(_) => {}
                Err(err) => unreachable!("noodles read: {err}"),
            }
            writer
                .write_variant_record(&header, &buffer)
                .expect("invariant: noodles writes");
        }
    }
    let text = String::from_utf8(out).expect("invariant: VCF text");
    // Drop noodles' header: compare records only.
    let mut out = String::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        out.push_str(line);
        out.push('\n');
    }
    out
}

/// Splits off the leading numeric run of `s` (`None` when it starts with a
/// literal).
fn take_numeric(s: &str) -> Option<(&str, &str)> {
    let len: usize = s
        .chars()
        .take_while(|c| c.is_ascii_digit() || matches!(c, '.' | '-' | '+' | 'e' | 'E'))
        .map(char::len_utf8)
        .sum();
    (len > 0).then(|| s.split_at(len))
}

/// Compares two rendered records token by token.
///
/// htslib renders floats outside `[1e-4, 1e6)` in `%g` notation
/// (`1.84742e+06`) while noodles prints them in fixed notation
/// (`1847420`); our renderer follows htslib (verified byte-for-byte against
/// `bcftools view -H` over the whole file), so numeric tokens only have to
/// agree *as values* here.
fn tokens_agree(ours: &str, theirs: &str) -> bool {
    if ours == theirs {
        return true;
    }
    let mut a = ours;
    let mut b = theirs;
    loop {
        if a.is_empty() || b.is_empty() {
            return a == b;
        }
        // Numeric runs compare by value (a lone `e` fails to parse and falls
        // through to the literal step).
        if let (Some((da, ra)), Some((db, rb))) = (take_numeric(a), take_numeric(b))
            && let (Ok(va), Ok(vb)) = (da.parse::<f64>(), db.parse::<f64>())
        {
            // BCF floats are f32; comparing the widened values catches a
            // wrong byte while tolerating the text form.
            // BCF floats are f32; the text forms may carry more digits.
            #[allow(clippy::cast_possible_truncation)]
            let same = (va as f32).to_bits() == (vb as f32).to_bits();
            if !same {
                return false;
            }
            a = ra;
            b = rb;
            continue;
        }
        let (ca, ra) = a.split_at(1);
        let (cb, rb) = b.split_at(1);
        if ca != cb {
            return false;
        }
        a = ra;
        b = rb;
    }
}

fn compare(path: &std::path::Path, limit: usize) {
    let ours = render_ours(path, limit);
    let theirs = render_noodles(path, limit);
    let ours_lines: Vec<&str> = ours.lines().collect();
    let theirs_lines: Vec<&str> = theirs.lines().collect();
    assert_eq!(ours_lines.len(), theirs_lines.len(), "record count differs");
    for (i, (a, b)) in ours_lines.iter().zip(&theirs_lines).enumerate() {
        assert!(
            tokens_agree(a, b),
            "record {i} differs\n  reed:    {a}\n  noodles: {b}"
        );
    }
}

#[test]
fn decoding_matches_noodles_on_the_fixture() {
    compare(&fixture("fix.bcf"), usize::MAX);
}

/// The local full-data BCF, when present (skipped otherwise).
#[test]
fn decoding_matches_noodles_on_the_full_file_prefix() {
    let path = std::path::Path::new("/sri/home/sri2025201067/reed/clean.bcf");
    if !path.exists() {
        return;
    }
    compare(path, 100_000);
}

#[test]
fn token_walker_accepts_equivalent_float_forms() {
    assert!(tokens_agree("A\t1.84742e+06\tB", "A\t1847420\tB"));
    assert!(tokens_agree("AF=0.539;DP=58", "AF=0.539;DP=58"));
    assert!(!tokens_agree("A\t1.84742e+06", "A\t1847421"));
    assert!(!tokens_agree("AC=375", "AC=376"));
    // A full record with the QUAL form differing only.
    let ours = "Chr01\t3654463\tChr01_3654463\tT\tA\t1.84742e+06\tPASS\tAC=375\tGT\t0/1:1.5e+03";
    let theirs = "Chr01\t3654463\tChr01_3654463\tT\tA\t1847420\tPASS\tAC=375\tGT\t0/1:1500";
    assert!(tokens_agree(ours, theirs));
}
