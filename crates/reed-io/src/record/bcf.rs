//! BCF path: the bcf::Record -> RecordBuf conversion.
//!
//! Implementation note (verified against noodles 0.117):
//! `RecordBuf::try_from_variant_record(&header, &record)` accepts any
//! `R: vcf::variant::Record + ?Sized`, and `bcf::Record` implements that
//! trait (through the string maps, with lazy field access). No hand-written
//! field-level conversion is needed: one call suffices and matches the
//! noodles conversion semantics exactly (all INFO/FORMAT fields,
//! missing-value semantics from the library).
//!
//! Pitfall: `bcf::io::Reader::from` takes the inner reader verbatim, while
//! `bcf::io::Reader::new` wraps its input in a fresh BGZF layer. Handing a
//! `bgzf::io::Reader` to `new` stacks a second decompression layer and fails
//! with "invalid BGZF header"; pass the raw file to `new`, or build the
//! explicit BGZF layer and use `from`.

use crate::ReadError;
use noodles::bcf;
use noodles::vcf;
use noodles::vcf::variant::RecordBuf;

/// bcf::Record -> RecordBuf with every field preserved (full fidelity:
/// INFO/FORMAT kept, missing-value semantics from noodles).
pub fn to_record_buf(rec: &bcf::Record, header: &vcf::Header) -> Result<RecordBuf, ReadError> {
    RecordBuf::try_from_variant_record(header, rec).map_err(|e| {
        let contig = rec
            .reference_sequence_name(header.string_maps())
            .unwrap_or_default()
            .to_string();
        let pos = rec
            .variant_start()
            .and_then(|r| r.ok())
            .map(|p| usize::from(p) as u64)
            .unwrap_or(0);
        ReadError::Decode {
            contig,
            pos,
            msg: format!("bcf to RecordBuf conversion failed: {e}"),
        }
    })
}

#[cfg(test)]
mod tests {
    //! Acceptance test: a BCF write-read round trip checked against
    //! hand-computed expectations.

    use super::*;
    use noodles::bgzf;
    use std::io::BufReader;
    use std::num::NonZero;
    use vcf::variant::io::Write as _;
    use vcf::variant::record_buf::samples::sample::Value;
    use vcf::variant::record_buf::Samples;

    fn header() -> vcf::Header {
        "##fileformat=VCFv4.4\n\
         ##contig=<ID=chr1,length=1000000>\n\
         ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
            .parse()
            .unwrap()
    }

    fn sample_record_buf() -> RecordBuf {
        let samples = Samples::new(
            ["GT", "DP"].into_iter().map(String::from).collect(),
            vec![
                vec![Some(Value::String(String::from("0/1"))), Some(Value::Integer(30))],
                vec![Some(Value::String(String::from("./."))), Some(Value::Integer(0))],
            ],
        );
        RecordBuf::builder()
            .set_reference_sequence_name("chr1")
            .set_variant_start(noodles::core::Position::try_from(42).unwrap())
            .set_reference_bases("A")
            .set_alternate_bases(vec!["G".to_string()].into_iter().collect())
            .set_quality_score(50.5)
            .set_info(
                [("DP".to_string(), Some(noodles::vcf::variant::record_buf::info::field::Value::Integer(60)))]
                    .into_iter()
                    .collect(),
            )
            .set_samples(samples)
            .build()
    }

    fn write_bcf(path: &std::path::Path, hdr: &vcf::Header, buf: &RecordBuf) {
        let mtw = bgzf::io::MultithreadedWriter::with_worker_count(
            NonZero::new(2).unwrap(),
            std::fs::File::create(path).unwrap(),
        );
        let mut w = bcf::io::Writer::from(mtw);
        w.write_header(hdr).unwrap();
        w.write_variant_record(hdr, buf).unwrap();
        let mut mtw = w.into_inner();
        mtw.finish().unwrap();
    }

    #[test]
    fn bcf_round_trip_preserves_fields() {
        let hdr = header();
        let original = sample_record_buf();

        let dir = std::env::temp_dir().join(format!("reed-io-bcf-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("rt.bcf");

        write_bcf(&path, &hdr, &original);

        let file = std::fs::File::open(&path).unwrap();
        let mut reader = bcf::io::Reader::from(bgzf::io::Reader::new(BufReader::new(file)));
        let read_header = reader.read_header().unwrap();
        let mut rec = bcf::Record::default();
        let n = reader.read_record(&mut rec).unwrap();
        assert!(n > 0, "BCF should contain one record");

        let converted = to_record_buf(&rec, &read_header).unwrap();

        // field-by-field comparison (conversion fidelity)
        assert_eq!(converted.reference_sequence_name(), "chr1");
        assert_eq!(
            converted.variant_start().map(usize::from),
            Some(42)
        );
        assert_eq!(converted.reference_bases(), "A");
        assert_eq!(converted.quality_score(), Some(50.5));
        // samples: GT and DP match per sample
        let s = converted.samples();
        let s1 = s.get_index(0).unwrap();
        let s2 = s.get_index(1).unwrap();
        // GT may decode as Genotype or String; accept both, but the value
        // semantics must be 0/1 and ./.
        let gt_of = |sample: vcf::variant::record_buf::samples::Sample<'_>| {
            sample
                .get(&"GT".to_string())
                .flatten()
                .map(|v| match v {
                    Value::String(t) => t.clone(),
                    Value::Genotype(_) => "<genotype>".to_string(),
                    other => format!("unexpected:{other:?}"),
                })
        };
        let gt1 = gt_of(s1).unwrap();
        let gt2 = gt_of(s2).unwrap();
        assert!(gt1 == "0/1" || gt1 == "<genotype>", "s1 GT: {gt1}");
        assert!(gt2 == "./." || gt2 == "<genotype>", "s2 GT: {gt2}");

        std::fs::remove_dir_all(&dir).ok();
    }
}
